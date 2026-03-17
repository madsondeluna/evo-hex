"""
Análise de frequência global de aminoácidos.

Processa as estruturas limpas sequencialmente (parsing PDB é CPU-bound
e o GIL limita o ganho de threads para este caso).
"""

import logging
from collections import Counter
from pathlib import Path

from Bio import PDB
from tqdm.auto import tqdm

from .config import STANDARD_AMINO_ACIDS, glob_pdb

logger = logging.getLogger(__name__)


def analyze_amino_acid_frequency(
    clean_dir: Path,
) -> tuple[Counter, dict[str, dict]]:
    """Conta aminoácidos em todas as estruturas limpas.

    Args:
        clean_dir: Diretório com PDBs limpos.

    Returns:
        Tupla (global_counter, per_structure).
        per_structure mapeia código PDB → dict {aa: contagem}.

    Raises:
        FileNotFoundError: Se clean_dir não existir ou não tiver PDBs.
    """
    pdb_files = glob_pdb(clean_dir)
    if not pdb_files:
        raise FileNotFoundError(f"Nenhuma estrutura limpa encontrada em {clean_dir}")

    print(f"\n{'='*60}")
    print("ANÁLISE DE FREQUÊNCIA DE AMINOÁCIDOS")
    print(f"{'='*60}\n")
    print(f"Analisando {len(pdb_files):,} estruturas...\n")

    global_counter: Counter = Counter()
    per_structure: dict[str, dict] = {}
    parser = PDB.PDBParser(QUIET=True)
    errors = 0

    with tqdm(total=len(pdb_files), desc="Analyzing") as pbar:
        for pdb_file in pdb_files:
            try:
                structure = parser.get_structure(pdb_file.stem, str(pdb_file))
                local_counter: Counter = Counter()

                for model in structure:
                    for chain in model:
                        for residue in chain:
                            if residue.id[0] == " ":
                                resname = residue.resname.strip()
                                if resname in STANDARD_AMINO_ACIDS:
                                    global_counter[resname] += 1
                                    local_counter[resname] += 1

                per_structure[pdb_file.stem] = dict(local_counter)

            except Exception as exc:  # noqa: BLE001
                logger.warning("Erro ao analisar %s: %s", pdb_file.stem, exc)
                errors += 1

            pbar.update(1)

    if errors:
        logger.warning("%d estruturas com erros foram ignoradas.", errors)

    return global_counter, per_structure


def save_frequency_results(
    analysis_path: Path,
    global_counter: Counter,
    per_structure: dict[str, dict],
) -> None:
    """Salva frequências globais e por estrutura em disco.

    Args:
        analysis_path: Diretório de saída das análises.
        global_counter: Contagem global de aminoácidos.
        per_structure: Contagem por código PDB.
    """
    analysis_path.mkdir(parents=True, exist_ok=True)
    total_residues = sum(global_counter.values())

    # Frequências globais (texto)
    lines = [
        "=" * 60,
        "FREQUÊNCIA GLOBAL DE AMINOÁCIDOS",
        "=" * 60 + "\n",
        f"Total de resíduos:   {total_residues:,}",
        f"Total de estruturas: {len(per_structure):,}\n",
        f"{'AA':<5} {'Contagem':<12} {'Frequência (%)':>15}",
        "-" * 60,
    ]
    for aa, count in global_counter.most_common():
        freq = count / total_residues * 100
        lines.append(f"{aa:<5} {count:<12,} {freq:>14.2f}%")

    freq_file = analysis_path / "amino_acid_frequencies_global.txt"
    freq_file.write_text("\n".join(lines) + "\n")
    print(f"\nFrequências globais: {freq_file}")

    # Por estrutura (CSV)
    all_aa = sorted(STANDARD_AMINO_ACIDS)
    header = "PDB_Code," + ",".join(all_aa) + ",Total"
    rows = [header]
    for pdb_code in sorted(per_structure):
        counts = per_structure[pdb_code]
        row = [pdb_code] + [str(counts.get(aa, 0)) for aa in all_aa]
        row.append(str(sum(counts.values())))
        rows.append(",".join(row))

    csv_file = analysis_path / "amino_acid_frequencies_per_structure.csv"
    csv_file.write_text("\n".join(rows) + "\n")
    print(f"Frequências por estrutura: {csv_file}")


def derive_frequencies_from_unified(analysis_path: Path) -> tuple[Counter, dict]:
    """Deriva global_counter dos resultados do passo DSSP unificado.

    Le o CSV de propensidades e o relatorio unificado para reconstruir um
    Counter aproximado de aminoacidos sem precisar re-parsear todos os PDBs.

    Args:
        analysis_path: Diretorio com os arquivos de analise (helix_propensities.csv,
                       unified_report.txt).

    Returns:
        Tupla (global_counter, per_structure).
        per_structure e sempre vazio — use amino_acid_frequencies_per_structure.csv
        se precisar de dados por estrutura.
    """
    import re
    import pandas as pd

    propensity_csv = analysis_path / "helix_propensities.csv"
    if not propensity_csv.exists():
        return Counter(), {}

    # Tenta ler total de residuos do relatorio unificado
    total_residues: int | None = None
    report_txt = analysis_path / "unified_report.txt"
    if report_txt.exists():
        for line in report_txt.read_text().splitlines():
            m = re.match(r"Total de residuos:\s+([\d,]+)", line)
            if m:
                total_residues = int(m.group(1).replace(",", ""))
                break

    df = pd.read_csv(propensity_csv)
    global_counter: Counter = Counter()

    from .config import ONE_TO_THREE
    for _, row in df.iterrows():
        aa3 = ONE_TO_THREE.get(row["AA"])
        if aa3:
            if total_residues:
                count = int(round(row["Freq_Total"] / 100.0 * total_residues))
            else:
                count = int(round(row["Freq_Total"] * 1000))
            if count > 0:
                global_counter[aa3] = count

    return global_counter, {}


def print_frequency_summary(global_counter: Counter, per_structure: dict) -> None:
    """Exibe resumo das análises de frequência no console.

    Args:
        global_counter: Contagem global de aminoácidos.
        per_structure: Contagem por código PDB.
    """
    total_residues = sum(global_counter.values())
    print(f"\n{'='*60}")
    print("RESUMO DA ANÁLISE DE FREQUÊNCIA")
    print(f"{'='*60}")
    print(f"\nTotal de estruturas analisadas: {len(per_structure):,}")
    print(f"Total de resíduos:              {total_residues:,}")
    print("\nTop 5 aminoácidos mais frequentes:")
    for aa, count in global_counter.most_common(5):
        freq = count / total_residues * 100
        print(f"  {aa}: {count:,} ({freq:.2f}%)")
    print(f"{'='*60}\n")
