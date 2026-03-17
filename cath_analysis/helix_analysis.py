"""
Análise avançada de composição de hélices alpha usando DSSP.

Inclui:
- Propensão observada vs teórica (Chou-Fasman)
- Preferência por posição (N-terminal / Middle / C-terminal)
- Distribuição de comprimentos
- Padrão heptad de hidrofobicidade (calculado sobre TODOS os resíduos,
  não apenas helicais – reflecte o padrão coil-to-helix)
"""

import logging
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import PDB
from Bio.PDB import DSSP
from tqdm.auto import tqdm

from .config import (
    AA_PROPERTIES,
    ANALYSIS_PATH,
    STRUCTURES_CLEAN_PATH,
    HELIX_PROPENSITY,
    HYDROPHOBIC_AA,
    STANDARD_AMINO_ACIDS,
)

logger = logging.getLogger(__name__)


# ── Utilitários ───────────────────────────────────────────────────────────────

def _check_dssp() -> bool:
    """Verifica se mkdssp está disponível no PATH.

    Returns:
        True se encontrado, False caso contrário.
    """
    if shutil.which("mkdssp"):
        logger.info("DSSP encontrado: %s", shutil.which("mkdssp"))
        return True

    print(
        "\n[ERRO] mkdssp não encontrado no PATH.\n"
        "  macOS:  brew install brewsci/bio/dssp\n"
        "  Linux:  apt install dssp\n"
        "  Manual: https://swift.cmbi.umcn.nl/gv/dssp/\n"
    )
    return False


def _extract_helix_segments(
    ss_sequence: list[str],
    aa_sequence: list[str],
    min_length: int = 4,
) -> tuple[Counter, dict, list[int]]:
    """Extrai resíduos helicais e posições de uma sequência de estrutura secundária.

    Args:
        ss_sequence: Lista de códigos DSSP por resíduo.
        aa_sequence: Lista de aminoácidos (1-letra) correspondentes.
        min_length: Comprimento mínimo para considerar uma hélice.

    Returns:
        Tupla (helix_residues, helix_positions, helix_lengths).
    """
    helix_residues: Counter = Counter()
    helix_positions: dict = defaultdict(Counter)
    helix_lengths: list[int] = []

    current: list[str] = []

    def _flush(segment: list[str]) -> None:
        if len(segment) < min_length:
            return
        helix_lengths.append(len(segment))
        for j, res in enumerate(segment):
            helix_residues[res] += 1
            if j == 0:
                helix_positions["N-terminal"][res] += 1
            elif j == len(segment) - 1:
                helix_positions["C-terminal"][res] += 1
            else:
                helix_positions["Middle"][res] += 1

    for ss, aa in zip(ss_sequence, aa_sequence):
        if ss in ("H", "G", "I"):
            current.append(aa)
        else:
            _flush(current)
            current = []
    _flush(current)  # última hélice

    return helix_residues, dict(helix_positions), helix_lengths


# ── Análise DSSP ──────────────────────────────────────────────────────────────

def analyze_single_structure_dssp(pdb_file: Path) -> dict:
    """Analisa uma estrutura com DSSP.

    Args:
        pdb_file: Arquivo PDB limpo.

    Returns:
        Dict com status e resultados, ou dict com status='error'.
    """
    pdb_code = pdb_file.stem
    parser = PDB.PDBParser(QUIET=True)

    try:
        structure = parser.get_structure(pdb_code, str(pdb_file))
        model = structure[0]
        dssp = DSSP(model, str(pdb_file), dssp="mkdssp")

        ss_seq = [dssp[k][2] for k in dssp.property_keys]
        aa_seq = [dssp[k][1] for k in dssp.property_keys]
        all_residues = Counter(aa_seq)

        helix_residues, helix_positions, helix_lengths = _extract_helix_segments(
            ss_seq, aa_seq
        )

        return {
            "pdb_code": pdb_code,
            "status": "success",
            "helix_residues": dict(helix_residues),
            "all_residues": dict(all_residues),
            "helix_positions": helix_positions,
            "helix_lengths": helix_lengths,
        }

    except Exception as exc:  # noqa: BLE001
        logger.debug("DSSP falhou para %s: %s", pdb_code, exc)
        return {"pdb_code": pdb_code, "status": "error", "message": str(exc)}


def analyze_secondary_structure_with_dssp(
    clean_dir: Path,
    analysis_path: Path,
) -> dict | None:
    """Agrega análise DSSP de todas as estruturas.

    Args:
        clean_dir: Diretório com PDBs limpos.
        analysis_path: Diretório de saída (não utilizado aqui, reservado).

    Returns:
        Dict agregado ou None se DSSP não estiver disponível.
    """
    if not _check_dssp():
        return None

    pdb_files = list(clean_dir.glob("*.pdb"))
    global_helix: Counter = Counter()
    global_all: Counter = Counter()
    global_positions: dict = defaultdict(lambda: defaultdict(int))
    global_lengths: list[int] = []
    success = failed = 0

    print(f"\nProcessando {len(pdb_files):,} estruturas com DSSP...\n")

    with tqdm(total=len(pdb_files), desc="DSSP Analysis") as pbar:
        for pdb_file in pdb_files:
            result = analyze_single_structure_dssp(pdb_file)

            if result["status"] == "success":
                success += 1
                for aa, n in result["helix_residues"].items():
                    global_helix[aa] += n
                for aa, n in result["all_residues"].items():
                    global_all[aa] += n
                for pos_type, pos_data in result["helix_positions"].items():
                    for aa, n in pos_data.items():
                        global_positions[pos_type][aa] += n
                global_lengths.extend(result["helix_lengths"])
            else:
                failed += 1

            pbar.update(1)

    print(f"\nDSSP concluído: {success} sucessos, {failed} falhas")

    return {
        "helix_residues": global_helix,
        "all_residues": global_all,
        # Converte defaultdict para dict simples antes de retornar
        "helix_positions": {k: dict(v) for k, v in global_positions.items()},
        "helix_lengths": global_lengths,
    }


# ── DataFrames analíticos ─────────────────────────────────────────────────────

def calculate_helix_propensities(
    helix_residues: Counter,
    all_residues: Counter,
) -> pd.DataFrame:
    """Calcula propensão observada e teórica (Chou-Fasman) para cada aminoácido.

    Args:
        helix_residues: Contagem de aminoácidos em hélices.
        all_residues: Contagem global de aminoácidos.

    Returns:
        DataFrame com colunas AA, Count_Helix, Count_Total, Freq_Helix,
        Freq_Total, Propensity_Observed, Propensity_Theoretical.
    """
    total_helix = sum(helix_residues.values())
    total_all = sum(all_residues.values())
    rows = []

    for aa in sorted(STANDARD_AMINO_ACIDS):
        h_count = helix_residues.get(aa, 0)
        a_count = all_residues.get(aa, 0)
        if a_count == 0:
            continue

        freq_h = h_count / total_helix * 100
        freq_a = a_count / total_all * 100
        prop_obs = freq_h / freq_a if freq_a > 0 else 0.0

        rows.append({
            "AA": aa,
            "Count_Helix": h_count,
            "Count_Total": a_count,
            "Freq_Helix": freq_h,
            "Freq_Total": freq_a,
            "Propensity_Observed": prop_obs,
            "Propensity_Theoretical": HELIX_PROPENSITY.get(aa, 1.0),
        })

    return pd.DataFrame(rows)


def analyze_helix_positions(helix_positions: dict) -> pd.DataFrame:
    """Tabela de frequências por posição na hélice.

    Args:
        helix_positions: Dict {posição: {aa: contagem}}.

    Returns:
        DataFrame com colunas Position, AA, Count, Frequency.
    """
    rows = []
    for pos in ("N-terminal", "Middle", "C-terminal"):
        pos_counts = helix_positions.get(pos, {})
        total = sum(pos_counts.values())
        for aa in sorted(STANDARD_AMINO_ACIDS):
            count = pos_counts.get(aa, 0)
            freq = count / total * 100 if total > 0 else 0.0
            rows.append({"Position": pos, "AA": aa, "Count": count, "Frequency": freq})
    return pd.DataFrame(rows)


def analyze_hydrophobic_patterns(clean_dir: Path) -> dict:
    """Calcula fração de resíduos hidrofóbicos em cada posição do heptad.

    Nota: este cálculo é feito sobre TODOS os resíduos da cadeia A (não
    apenas helicais), capturando o contexto geral de periodicidade.

    Args:
        clean_dir: Diretório com PDBs limpos.

    Returns:
        Dict {posição 0-6: {'hydrophobic': int, 'total': int}}.
    """
    print("\nAnalisando padrões de hidrofobicidade (heptad)...")
    parser = PDB.PDBParser(QUIET=True)
    heptad: dict = {i: {"hydrophobic": 0, "total": 0} for i in range(7)}

    pdb_files = list(clean_dir.glob("*.pdb"))
    with tqdm(total=len(pdb_files), desc="Hydrophobic patterns") as pbar:
        for pdb_file in pdb_files:
            try:
                structure = parser.get_structure(pdb_file.stem, str(pdb_file))
                for model in structure:
                    for chain in model:
                        residues = [r for r in chain if r.id[0] == " "]
                        for i, residue in enumerate(residues):
                            pos = i % 7
                            resname = residue.resname.strip()
                            heptad[pos]["total"] += 1
                            if resname in HYDROPHOBIC_AA:
                                heptad[pos]["hydrophobic"] += 1
            except Exception as exc:  # noqa: BLE001
                logger.debug("Erro no padrão heptad para %s: %s", pdb_file.stem, exc)
            pbar.update(1)

    return heptad


def analyze_aa_properties_distribution(all_residues: Counter) -> pd.DataFrame:
    """Frequência relativa por grupo físico-químico.

    Args:
        all_residues: Contagem global de aminoácidos.

    Returns:
        DataFrame com colunas Property, Count, Frequency.
    """
    total = sum(all_residues.values())
    rows = [
        {
            "Property": prop,
            "Count": (n := sum(all_residues.get(aa, 0) for aa in aa_list)),
            "Frequency": n / total * 100 if total > 0 else 0.0,
        }
        for prop, aa_list in AA_PROPERTIES.items()
    ]
    return pd.DataFrame(rows)


# ── Persistência ──────────────────────────────────────────────────────────────

def save_comprehensive_report(
    analysis_path: Path,
    propensity_df: pd.DataFrame,
    position_df: pd.DataFrame,
    properties_df: pd.DataFrame,
    helix_lengths: list[int],
) -> None:
    """Gera relatório TXT e CSVs da análise avançada de hélices.

    Args:
        analysis_path: Diretório de saída.
        propensity_df: DataFrame de propensões.
        position_df: DataFrame de posições.
        properties_df: DataFrame de propriedades.
        helix_lengths: Lista de comprimentos de hélices.
    """
    report_file = analysis_path / "helix_analysis_comprehensive_report.txt"

    with report_file.open("w") as f:
        f.write("=" * 80 + "\n")
        f.write("ANÁLISE ABRANGENTE DE COMPOSIÇÃO DE HÉLICES ALPHA\n")
        f.write("Estruturas CATH Mainly-Alpha\n")
        f.write("=" * 80 + "\n\n")

        f.write("### PROPENSÃO DE AMINOÁCIDOS PARA HÉLICES ###\n\n")
        f.write(f"{'AA':<5} {'Helix':<10} {'Total':<10} {'Prop.Obs':<12} {'Prop.Teor':<12}\n")
        f.write("-" * 80 + "\n")
        for _, row in propensity_df.iterrows():
            f.write(
                f"{row['AA']:<5} {row['Count_Helix']:<10} {row['Count_Total']:<10} "
                f"{row['Propensity_Observed']:<12.2f} {row['Propensity_Theoretical']:<12.2f}\n"
            )

        f.write("\n### TOP 5 FORMADORES DE HÉLICES (Propensão Observada) ###\n")
        top5 = propensity_df.nlargest(5, "Propensity_Observed")
        for rank, (_, row) in enumerate(top5.iterrows(), start=1):
            f.write(f"  {rank}. {row['AA']}: {row['Propensity_Observed']:.2f}\n")

        f.write("\n### DISTRIBUIÇÃO POR PROPRIEDADES ###\n\n")
        for _, row in properties_df.iterrows():
            f.write(f"{row['Property']:<25}: {row['Frequency']:>6.2f}% ({row['Count']:,})\n")

        if helix_lengths:
            arr = np.array(helix_lengths)
            f.write("\n### ESTATÍSTICAS DE COMPRIMENTO DE HÉLICES ###\n\n")
            f.write(f"  Número total de hélices: {len(arr):,}\n")
            f.write(f"  Comprimento médio:        {arr.mean():.1f} resíduos\n")
            f.write(f"  Mediana:                  {np.median(arr):.1f} resíduos\n")
            f.write(f"  Desvio padrão:            {arr.std():.1f}\n")
            f.write(f"  Intervalo:                {arr.min()} – {arr.max()} resíduos\n")
            f.write("\n  Percentis:\n")
            for p in (25, 50, 75, 90, 95):
                f.write(f"    {p}%: {np.percentile(arr, p):.1f} resíduos\n")

        f.write("\n" + "=" * 80 + "\n")

    print(f"\nRelatório abrangente salvo em: {report_file}")

    propensity_df.to_csv(analysis_path / "helix_propensities.csv", index=False)
    position_df.to_csv(analysis_path / "helix_positions.csv", index=False)
    properties_df.to_csv(analysis_path / "aa_properties_distribution.csv", index=False)
    print("CSVs salvos: helix_propensities, helix_positions, aa_properties_distribution")


# ── Ponto de entrada ──────────────────────────────────────────────────────────

def run_helix_analysis() -> None:
    """Executa o pipeline completo de análise avançada de hélices."""
    # Import local para evitar dependência circular de plotting
    from .plotting import (  # noqa: PLC0415
        plot_helix_length_distribution,
        plot_helix_positions,
        plot_helix_propensities,
        plot_hydrophobic_heptad,
    )

    clean_dir = STRUCTURES_CLEAN_PATH
    analysis_path = ANALYSIS_PATH
    analysis_path.mkdir(exist_ok=True)

    print("=" * 80)
    print("ANÁLISE AVANÇADA DE HÉLICES ALPHA")
    print("=" * 80)

    dssp_results = analyze_secondary_structure_with_dssp(clean_dir, analysis_path)
    if dssp_results is None:
        return

    propensity_df = calculate_helix_propensities(
        dssp_results["helix_residues"], dssp_results["all_residues"]
    )
    position_df = analyze_helix_positions(dssp_results["helix_positions"])
    properties_df = analyze_aa_properties_distribution(dssp_results["all_residues"])
    heptad_data = analyze_hydrophobic_patterns(clean_dir)

    print("\nGerando visualizações...")
    plot_helix_propensities(propensity_df, analysis_path)
    plot_helix_positions(position_df, analysis_path)
    plot_helix_length_distribution(dssp_results["helix_lengths"], analysis_path)
    plot_hydrophobic_heptad(heptad_data, analysis_path)

    save_comprehensive_report(
        analysis_path,
        propensity_df,
        position_df,
        properties_df,
        dssp_results["helix_lengths"],
    )

    print("\n" + "=" * 80)
    print("ANÁLISE CONCLUÍDA")
    print("=" * 80)
