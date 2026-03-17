"""
Classificação e análise de tipos de hélices (H, G, I) usando DSSP.

Distingue entre:
  H – Alpha helix (i→i+4)
  G – 3-10 helix  (i→i+3)
  I – Pi helix    (i→i+5)

e compara composição de aminoácidos entre eles.
"""

import logging
import shutil
import warnings
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

warnings.filterwarnings("ignore", category=UserWarning, message="parse error")

import numpy as np
import pandas as pd
from Bio import PDB
from Bio.PDB import DSSP
from tqdm.auto import tqdm

from .config import (
    ANALYSIS_PATH,
    HELIX_CHARACTERISTICS,
    HELIX_TYPES,
    PROCESS_WORKERS,
    STANDARD_AMINO_ACIDS,
    STRUCTURES_CLEAN_PATH,
    glob_pdb,
)

logger = logging.getLogger(__name__)

_MIN_HELIX_LEN = 3  # comprimento mínimo de resíduos para registrar uma hélice


# ── Classificação ─────────────────────────────────────────────────────────────

def classify_helix_types(clean_dir: Path) -> dict:
    """Classifica hélices por tipo (H/G/I) em todas as estruturas.

    Args:
        clean_dir: Diretório com PDBs limpos.

    Returns:
        Dict com chaves:
          - helix_counts       : Counter {tipo: n_hélices}
          - helix_residues     : dict {tipo: Counter {aa: n}}
          - helix_lengths      : dict {tipo: [comprimentos]}
          - helix_positions    : dict {tipo: {posição: Counter {aa: n}}}
          - transitions        : dict {tipo_origem: Counter {tipo_destino: n}}
          - structures_analyzed: int
    """
    if not shutil.which("mkdssp"):
        raise RuntimeError(
            "mkdssp não encontrado. Instale com: brew install brewsci/bio/dssp"
        )

    pdb_files = glob_pdb(clean_dir)

    helix_counts: Counter = Counter()
    helix_residues: dict = defaultdict(Counter)
    helix_lengths: dict = defaultdict(list)
    helix_positions: dict = defaultdict(lambda: defaultdict(Counter))
    transitions: dict = defaultdict(Counter)

    successful = failed = 0

    print(f"\nClassificando tipos de hélices em {len(pdb_files):,} estruturas (workers: {PROCESS_WORKERS})...")

    with ThreadPoolExecutor(max_workers=PROCESS_WORKERS) as executor:
        futures = {executor.submit(_classify_single_pdb, f): f for f in pdb_files}
        with tqdm(total=len(pdb_files), desc="Classifying helix types") as pbar:
            for future in as_completed(futures):
                result = future.result()
                if result is not None:
                    ss_seq, aa_seq = result
                    _process_helix_sequence(
                        ss_seq, aa_seq,
                        helix_counts, helix_residues, helix_lengths,
                        helix_positions, transitions,
                    )
                    successful += 1
                else:
                    failed += 1
                pbar.update(1)

    print(f"Classificação concluída: {successful} sucessos, {failed} falhas")

    return {
        "helix_counts": dict(helix_counts),
        "helix_residues": dict(helix_residues),
        "helix_lengths": dict(helix_lengths),
        "helix_positions": {
            ht: {pos: dict(cnt) for pos, cnt in pos_data.items()}
            for ht, pos_data in helix_positions.items()
        },
        "transitions": dict(transitions),
        "structures_analyzed": successful,
    }


def _classify_single_pdb(pdb_file: Path) -> tuple[list, list] | None:
    """Roda DSSP em um único PDB e retorna (ss_seq, aa_seq) ou None em caso de erro."""
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_file.stem, str(pdb_file))
        model = structure[0]
        dssp = DSSP(model, str(pdb_file), dssp="mkdssp", file_type="PDB")
        ss_seq = [dssp[k][2] for k in dssp.property_keys]
        aa_seq = [dssp[k][1] for k in dssp.property_keys]
        return ss_seq, aa_seq
    except Exception as exc:  # noqa: BLE001
        logger.debug("Falha em %s: %s", pdb_file.stem, exc)
        return None


def _process_helix_sequence(
    ss_seq: list[str],
    aa_seq: list[str],
    helix_counts: Counter,
    helix_residues: dict,
    helix_lengths: dict,
    helix_positions: dict,
    transitions: dict,
) -> None:
    """Itera sobre a sequência SS e registra segmentos helicais por tipo.

    Mantido separado para facilitar testes unitários.
    """
    current: dict = {"type": None, "residues": []}

    def _flush(segment: dict) -> None:
        h_type = segment["type"]
        residues = segment["residues"]
        if h_type is None or len(residues) < _MIN_HELIX_LEN:
            return
        helix_counts[h_type] += 1
        helix_lengths[h_type].append(len(residues))
        for j, res in enumerate(residues):
            helix_residues[h_type][res] += 1
            if j == 0:
                helix_positions[h_type]["N-terminal"][res] += 1
            elif j == len(residues) - 1:
                helix_positions[h_type]["C-terminal"][res] += 1
            else:
                helix_positions[h_type]["Middle"][res] += 1

    for ss, aa in zip(ss_seq, aa_seq):
        if ss in ("H", "G", "I"):
            if ss == current["type"]:
                current["residues"].append(aa)
            else:
                # Mudança de tipo de hélice – registra transição
                if current["type"] is not None:
                    transitions[current["type"]][ss] += 1
                _flush(current)
                current = {"type": ss, "residues": [aa]}
        else:
            _flush(current)
            current = {"type": None, "residues": []}

    _flush(current)


# ── DataFrames ────────────────────────────────────────────────────────────────

def analyze_helix_type_composition(helix_residues: dict) -> pd.DataFrame:
    """Frequência de aminoácidos por tipo de hélice.

    Args:
        helix_residues: Dict {tipo: Counter {aa: n}}.

    Returns:
        DataFrame com colunas Helix_Type, Helix_Name, AA, Count, Frequency.
    """
    rows = []
    for h_type, residue_counts in helix_residues.items():
        total = sum(residue_counts.values())
        h_name = HELIX_TYPES.get(h_type, h_type)
        for aa in sorted(STANDARD_AMINO_ACIDS):
            count = residue_counts.get(aa, 0)
            rows.append({
                "Helix_Type": h_type,
                "Helix_Name": h_name,
                "AA": aa,
                "Count": count,
                "Frequency": count / total * 100 if total > 0 else 0.0,
            })
    return pd.DataFrame(rows)


def compare_helix_types(helix_residues: dict) -> pd.DataFrame:
    """Top-10 aminoácidos por tipo de hélice (ranking por frequência).

    Args:
        helix_residues: Dict {tipo: Counter {aa: n}}.

    Returns:
        DataFrame com colunas Helix_Type, Helix_Name, Rank, AA, Count, Frequency.
    """
    rows = []
    for h_type, residue_counts in helix_residues.items():
        total = sum(residue_counts.values())
        h_name = HELIX_TYPES.get(h_type, h_type)
        for rank, (aa, count) in enumerate(residue_counts.most_common(10), start=1):
            rows.append({
                "Helix_Type": h_type,
                "Helix_Name": h_name,
                "Rank": rank,
                "AA": aa,
                "Count": count,
                "Frequency": count / total * 100 if total > 0 else 0.0,
            })
    return pd.DataFrame(rows)


def statistical_comparison(composition_df: pd.DataFrame) -> pd.DataFrame:
    """Diferença de frequência entre Alpha (H) e 3-10 (G) para cada aminoácido.

    Args:
        composition_df: Retorno de analyze_helix_type_composition.

    Returns:
        DataFrame ordenado pela magnitude da diferença (decrescente).
        Vazio se não houver ambos H e G.
    """
    types_present = composition_df["Helix_Type"].unique()
    if "H" not in types_present or "G" not in types_present:
        return pd.DataFrame()

    rows = []
    for aa in sorted(STANDARD_AMINO_ACIDS):
        aa_data = composition_df[composition_df["AA"] == aa]
        alpha_rows = aa_data[aa_data["Helix_Type"] == "H"]["Frequency"].values
        g310_rows = aa_data[aa_data["Helix_Type"] == "G"]["Frequency"].values

        if alpha_rows.size == 0 or g310_rows.size == 0:
            continue

        diff = float(alpha_rows[0] - g310_rows[0])
        rows.append({
            "AA": aa,
            "Alpha_Freq": float(alpha_rows[0]),
            "3-10_Freq": float(g310_rows[0]),
            "Difference": diff,
            "Enriched_In": "Alpha" if diff > 0 else "3-10",
        })

    df = pd.DataFrame(rows)
    return df.sort_values("Difference", key=abs, ascending=False).reset_index(drop=True)


# ── Relatório ─────────────────────────────────────────────────────────────────

def save_helix_type_report(
    analysis_path: Path,
    results: dict,
    composition_df: pd.DataFrame,
    comparison_df: pd.DataFrame,
    stat_df: pd.DataFrame,
) -> None:
    """Gera relatório TXT e CSVs da análise por tipo de hélice.

    Args:
        analysis_path: Diretório de saída.
        results: Retorno de classify_helix_types.
        composition_df: Retorno de analyze_helix_type_composition.
        comparison_df: Retorno de compare_helix_types.
        stat_df: Retorno de statistical_comparison.
    """
    report_file = analysis_path / "helix_types_comprehensive_report.txt"
    total_helices = sum(results["helix_counts"].values())

    with report_file.open("w") as f:
        f.write("=" * 80 + "\n")
        f.write("ANÁLISE DETALHADA DE TIPOS DE HÉLICES\n")
        f.write("Estruturas Evo-Hex\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Estruturas analisadas: {results['structures_analyzed']}\n\n")

        f.write("### DISTRIBUIÇÃO DE TIPOS DE HÉLICES ###\n\n")
        f.write(f"{'Tipo':<6} {'Nome':<20} {'Contagem':<12} {'%':<10}\n")
        f.write("-" * 60 + "\n")
        for h_type in ("H", "G", "I"):
            if h_type not in results["helix_counts"]:
                continue
            count = results["helix_counts"][h_type]
            pct = count / total_helices * 100 if total_helices else 0.0
            f.write(f"{h_type:<6} {HELIX_TYPES[h_type]:<20} {count:<12,} {pct:>8.1f}%\n")

        f.write("\n### CARACTERÍSTICAS ESTRUTURAIS ###\n\n")
        for h_type in ("H", "G", "I"):
            if h_type not in results["helix_counts"]:
                continue
            char = HELIX_CHARACTERISTICS[h_type]
            f.write(f"{char['name']}:\n")
            f.write(f"  Resíduos/volta:    {char['residues_per_turn']}\n")
            f.write(f"  Aumento/resíduo:   {char['rise_per_residue']} Å\n")
            f.write(f"  Ligação H:         {char['hydrogen_bond']}\n")
            f.write(f"  Comprimento típico:{char['typical_length']}\n\n")

        f.write("### ESTATÍSTICAS DE COMPRIMENTO ###\n\n")
        for h_type in ("H", "G", "I"):
            if h_type not in results["helix_lengths"]:
                continue
            lengths = np.array(results["helix_lengths"][h_type])
            f.write(f"{HELIX_TYPES[h_type]}:\n")
            f.write(f"  Média:   {lengths.mean():.1f} resíduos\n")
            f.write(f"  Mediana: {np.median(lengths):.1f} resíduos\n")
            f.write(f"  DP:      {lengths.std():.1f}\n")
            f.write(f"  Mín-Máx: {lengths.min()} – {lengths.max()} resíduos\n\n")

        f.write("### TOP 5 AMINOÁCIDOS POR TIPO ###\n\n")
        for h_type in ("H", "G", "I"):
            if h_type not in results["helix_residues"]:
                continue
            residues = results["helix_residues"][h_type]
            total = sum(residues.values())
            f.write(f"{HELIX_TYPES[h_type]}:\n")
            for rank, (aa, count) in enumerate(residues.most_common(5), start=1):
                freq = count / total * 100 if total > 0 else 0.0
                f.write(f"  {rank}. {aa}: {freq:.2f}% ({count:,})\n")
            f.write("\n")

        if not stat_df.empty:
            f.write("### DIFERENÇAS Alpha vs 3-10 (Top 10) ###\n\n")
            f.write(f"{'AA':<5} {'Alpha (%)':<12} {'3-10 (%)':<12} {'Diff':<12} {'Enriq.':<10}\n")
            f.write("-" * 60 + "\n")
            for _, row in stat_df.head(10).iterrows():
                f.write(
                    f"{row['AA']:<5} {row['Alpha_Freq']:<12.2f} "
                    f"{row['3-10_Freq']:<12.2f} {row['Difference']:<12.2f} "
                    f"{row['Enriched_In']:<10}\n"
                )

        f.write("\n" + "=" * 80 + "\n")

    print(f"\nRelatório salvo em: {report_file}")

    composition_df.to_csv(analysis_path / "helix_type_composition.csv", index=False)
    comparison_df.to_csv(analysis_path / "helix_type_top_residues.csv", index=False)
    if not stat_df.empty:
        stat_df.to_csv(
            analysis_path / "helix_type_statistical_comparison.csv", index=False
        )
    print("CSVs salvos: helix_type_composition, helix_type_top_residues, helix_type_statistical_comparison")


# ── Ponto de entrada ──────────────────────────────────────────────────────────

def run_helix_type_analysis() -> None:
    """Executa o pipeline completo de classificação de tipos de hélices."""
    from .plotting import (  # noqa: PLC0415
        plot_helix_length_comparison,
        plot_helix_type_composition,
        plot_helix_type_distribution,
        plot_top_amino_acids_by_type,
        plot_statistical_differences,
    )

    clean_dir = STRUCTURES_CLEAN_PATH
    analysis_path = ANALYSIS_PATH
    analysis_path.mkdir(exist_ok=True)

    print("=" * 80)
    print("ANÁLISE DE TIPOS DE HÉLICES E COMPOSIÇÃO")
    print("=" * 80)

    results = classify_helix_types(clean_dir)
    composition_df = analyze_helix_type_composition(results["helix_residues"])
    comparison_df = compare_helix_types(results["helix_residues"])
    stat_df = statistical_comparison(composition_df)

    print("\nGerando visualizações...")
    plot_helix_type_distribution(results["helix_counts"], analysis_path)
    plot_helix_type_composition(composition_df, analysis_path)
    plot_top_amino_acids_by_type(comparison_df, analysis_path)
    plot_helix_length_comparison(results["helix_lengths"], analysis_path)
    if not stat_df.empty:
        plot_statistical_differences(stat_df, analysis_path)

    save_helix_type_report(analysis_path, results, composition_df, comparison_df, stat_df)

    print("\n" + "=" * 80)
    print("RESUMO")
    print("=" * 80)
    for h_type in ("H", "G", "I"):
        if h_type not in results["helix_counts"]:
            continue
        count = results["helix_counts"][h_type]
        residues = results["helix_residues"][h_type]
        top_aa = residues.most_common(1)[0] if residues else ("N/A", 0)
        lengths = results["helix_lengths"].get(h_type, [])

        print(f"\n{HELIX_TYPES[h_type]}:")
        print(f"  Total: {count:,} hélices")
        print(f"  AA mais frequente: {top_aa[0]} ({top_aa[1]:,})")
        if lengths:
            print(f"  Comprimento médio: {np.mean(lengths):.1f} resíduos")

    print("\n" + "=" * 80)
    print("CONCLUÍDO")
    print("=" * 80)
