"""
Coleta de dados evolutivos para hélices alpha em um único passo DSSP.

Todas as funções de coleta percorrem o diretório de estruturas uma única vez,
acumulando os dados necessários para todos os gráficos evolutivos.
"""

import math
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from .config import HYDROPHOBIC_AA, STANDARD_AMINO_ACIDS, glob_pdb


# Mapeamento 3→1 letra (fallback manual para garantir robustez)
_THREE_TO_ONE: dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

_EISENBERG_SCALE: dict[str, float] = {
    "ALA":  0.62, "ARG": -2.53, "ASN": -0.78, "ASP": -0.90,
    "CYS":  0.29, "GLN": -0.85, "GLU": -0.74, "GLY":  0.48,
    "HIS": -0.40, "ILE":  1.38, "LEU":  1.06, "LYS": -1.50,
    "MET":  0.64, "PHE":  1.19, "PRO":  0.12, "SER": -0.18,
    "THR": -0.05, "TRP":  0.81, "TYR":  0.26, "VAL":  1.08,
}

_HELIX_PROPENSITY_CF: dict[str, float] = {
    "ALA": 1.42, "GLU": 1.51, "LEU": 1.21, "MET": 1.45,
    "GLN": 1.11, "LYS": 1.16, "ARG": 0.98, "HIS": 1.00,
    "VAL": 1.06, "ILE": 1.08, "TYR": 0.69, "PHE": 1.13,
    "TRP": 1.08, "THR": 0.83, "SER": 0.77, "CYS": 0.70,
    "ASP": 1.01, "ASN": 0.67, "GLY": 0.57, "PRO": 0.57,
}

_CODON_DEGENERACY: dict[str, int] = {
    "ALA": 4, "ARG": 6, "ASN": 2, "ASP": 2, "CYS": 2,
    "GLN": 2, "GLU": 2, "GLY": 4, "HIS": 2, "ILE": 3,
    "LEU": 6, "LYS": 2, "MET": 1, "PHE": 2, "PRO": 4,
    "SER": 6, "THR": 4, "TRP": 1, "TYR": 2, "VAL": 4,
}

_PROTEOME_FREQ: dict[str, float] = {
    "ALA": 6.97, "ARG": 5.53, "ASN": 4.06, "ASP": 5.25, "CYS": 2.27,
    "GLN": 3.93, "GLU": 6.75, "GLY": 6.87, "HIS": 2.29, "ILE": 5.49,
    "LEU": 9.68, "LYS": 5.19, "MET": 2.32, "PHE": 3.87, "PRO": 5.02,
    "SER": 7.14, "THR": 5.57, "TRP": 1.33, "TYR": 3.21, "VAL": 6.47,
}


def _res3_to_1(res3: str) -> str | None:
    """Converte código de 3 letras para 1 letra; retorna None se não encontrado."""
    try:
        from Bio.Data.IUPACData import protein_letters_3to1
        result = protein_letters_3to1.get(res3.upper(), None)
        if result:
            return result
    except Exception:
        pass
    return _THREE_TO_ONE.get(res3.upper(), None)


def collect_evolutionary_data(clean_dir: Path) -> dict:
    """Coleta dados para gráficos evolutivos em um único passo DSSP.

    Percorre todos os arquivos .pdb em clean_dir, roda DSSP em cada um e
    acumula os dados necessários para todos os gráficos evolutivos.

    Args:
        clean_dir: Diretório com arquivos .pdb limpos.

    Returns:
        Dict com chaves:
        - helix_sequences: list[list[str]] – sequências 1-letra por hélice (tipo H, ≥4 res)
        - ncap_residues: Counter – AAs nas posições 0,1,2 de hélices H
        - ccap_residues: Counter – AAs nas posições -1,-2,-3 de hélices H
        - ncap_position: dict{0: Counter, 1: Counter, 2: Counter}
        - ccap_position: dict{-1: Counter, -2: Counter, -3: Counter}
        - helix_content_per_structure: list[float] – fração de resíduos em hélices por estrutura
        - per_helix_data: list[dict{'length': int, 'aa_counts': Counter}]
        - heptad_aa_distribution: dict{0..6: Counter}
        - transitions: Counter de (from_type, to_type)
        - helix_lengths_by_type: dict{'H': list, 'G': list, 'I': list}
    """
    from Bio.PDB import PDBParser, DSSP

    pdb_files = sorted(glob_pdb(clean_dir))
    if not pdb_files:
        raise FileNotFoundError(f"Nenhum .pdb encontrado em {clean_dir}")

    try:
        from tqdm import tqdm
        iterator = tqdm(pdb_files, desc="Coletando dados evolutivos", unit="pdb")
    except ImportError:
        iterator = pdb_files

    helix_sequences: list[list[str]] = []
    ncap_residues: Counter = Counter()
    ccap_residues: Counter = Counter()
    ncap_position: dict[int, Counter] = {0: Counter(), 1: Counter(), 2: Counter()}
    ccap_position: dict[int, Counter] = {-1: Counter(), -2: Counter(), -3: Counter()}
    helix_content_per_structure: list[float] = []
    per_helix_data: list[dict] = []
    heptad_aa_distribution: dict[int, Counter] = {i: Counter() for i in range(7)}
    transitions: Counter = Counter()
    helix_lengths_by_type: dict[str, list] = {"H": [], "G": [], "I": []}

    parser = PDBParser(QUIET=True)
    helix_dssp_codes = {"H", "G", "I"}

    for pdb_path in iterator:
        try:
            structure = parser.get_structure(pdb_path.stem, str(pdb_path))
            model = structure[0]
            dssp = DSSP(model, str(pdb_path), dssp="mkdssp")
        except Exception:
            continue

        # Recolhe sequência de DSSP na ordem dos resíduos
        dssp_records = list(dssp)
        total_residues = len(dssp_records)
        if total_residues == 0:
            continue

        helix_residue_count = sum(1 for r in dssp_records if r[2] in helix_dssp_codes)
        helix_content_per_structure.append(helix_residue_count / total_residues)

        # Coleta heptad para todos os resíduos
        for global_pos, record in enumerate(dssp_records):
            res3 = record[1].resname if hasattr(record[1], "resname") else None
            # record layout: (dssp_index, res_id, ss, acc, phi, psi, ...)
            # Para o heptad, usa 3-letra do resíduo
            # record[1] é o resíduo Bio.PDB; mais seguro usar record diretamente
            # Tentamos pegar 3-letra pelo resíduo
            try:
                res_obj = record[1]  # pode ser ResidueKey ou similar
                # Em versões modernas do Biopython, dssp retorna:
                # (dssp_idx, residue, ss, acc, phi, psi, ...)
                # onde residue é um objeto Residue
                if hasattr(res_obj, "resname"):
                    res3 = res_obj.resname.strip()
                else:
                    res3 = None
            except Exception:
                res3 = None

            aa1 = _res3_to_1(res3) if res3 else None
            if aa1:
                heptad_pos = global_pos % 7
                heptad_aa_distribution[heptad_pos][aa1] += 1

        # Segmenta em hélices contíguas
        # Percorre a lista e extrai segmentos contíguos do mesmo tipo
        i = 0
        n = len(dssp_records)
        prev_type: str | None = None

        while i < n:
            record = dssp_records[i]
            ss = record[2]

            if ss not in helix_dssp_codes:
                if prev_type is not None:
                    prev_type = None
                i += 1
                continue

            # Início de um segmento de hélice
            helix_type = ss
            segment_start = i
            segment = []

            while i < n and dssp_records[i][2] == helix_type:
                r = dssp_records[i]
                try:
                    res_obj = r[1]
                    res3 = res_obj.resname.strip() if hasattr(res_obj, "resname") else None
                except Exception:
                    res3 = None
                aa1 = _res3_to_1(res3) if res3 else None
                if aa1:
                    segment.append(aa1)
                i += 1

            seg_len = len(segment)

            # Transição
            if prev_type is not None and prev_type != helix_type:
                transitions[(prev_type, helix_type)] += 1

            prev_type = helix_type

            # Comprimento por tipo
            if helix_type in helix_lengths_by_type:
                helix_lengths_by_type[helix_type].append(seg_len)

            # Dados por hélice (todos os tipos)
            if seg_len > 0:
                per_helix_data.append({
                    "length": seg_len,
                    "aa_counts": Counter(segment),
                })

            # Apenas hélices H (alpha) com ≥4 resíduos para sequências/ncap/ccap
            if helix_type == "H" and seg_len >= 4:
                helix_sequences.append(segment)

                # N-cap: posições 0, 1, 2
                for pos in range(min(3, seg_len)):
                    aa = segment[pos]
                    ncap_residues[aa] += 1
                    ncap_position[pos][aa] += 1

                # C-cap: posições -1, -2, -3
                for offset, pos_key in enumerate([-1, -2, -3]):
                    idx = seg_len - 1 - offset
                    if idx >= 0:
                        aa = segment[idx]
                        ccap_residues[aa] += 1
                        ccap_position[pos_key][aa] += 1

    return {
        "helix_sequences": helix_sequences,
        "ncap_residues": ncap_residues,
        "ccap_residues": ccap_residues,
        "ncap_position": ncap_position,
        "ccap_position": ccap_position,
        "helix_content_per_structure": helix_content_per_structure,
        "per_helix_data": per_helix_data,
        "heptad_aa_distribution": heptad_aa_distribution,
        "transitions": transitions,
        "helix_lengths_by_type": helix_lengths_by_type,
    }


def compute_hydrophobic_moments(helix_sequences: list) -> list[float]:
    """Calcula o momento hidrofóbico de Eisenberg para cada hélice.

    Para cada resíduo na posição i, o ângulo é i * 100° (convertido para radianos).
    μH = (1/N) * sqrt( (Σ Hi*sin(θi))² + (Σ Hi*cos(θi))² )

    Args:
        helix_sequences: Lista de listas de AAs em código 1-letra.

    Returns:
        Lista de valores μH para hélices com ≥4 resíduos.
    """
    moments: list[float] = []
    angle_per_residue = math.radians(100.0)

    for seq in helix_sequences:
        if len(seq) < 4:
            continue

        sin_sum = 0.0
        cos_sum = 0.0
        count = 0

        for i, aa1 in enumerate(seq):
            # Converte 1-letra para 3-letras para consultar escala
            aa3 = next((k for k, v in _THREE_TO_ONE.items() if v == aa1), None)
            if aa3 is None or aa3 not in _EISENBERG_SCALE:
                continue
            h = _EISENBERG_SCALE[aa3]
            theta = i * angle_per_residue
            sin_sum += h * math.sin(theta)
            cos_sum += h * math.cos(theta)
            count += 1

        if count >= 4:
            mu_h = (1.0 / count) * math.sqrt(sin_sum ** 2 + cos_sum ** 2)
            moments.append(mu_h)

    return moments


def compute_aa_cooccurrence(helix_sequences: list) -> pd.DataFrame:
    """Matriz 20×20 de co-ocorrência de AAs em hélices.

    Para cada hélice, conta todos os pares (ordem irrelevante) de AAs que
    co-ocorrem. O resultado é normalizado pelo total de pares.

    Args:
        helix_sequences: Lista de listas de AAs em código 1-letra.

    Returns:
        DataFrame 20×20 com índice e colunas = AAs ordenados alfabeticamente.
    """
    all_aa_1 = sorted(set(_THREE_TO_ONE.values()))
    aa_to_idx = {aa: i for i, aa in enumerate(all_aa_1)}
    n = len(all_aa_1)
    matrix = np.zeros((n, n), dtype=float)

    total_pairs = 0
    for seq in helix_sequences:
        # Conta AAs presentes na hélice
        present = [aa for aa in seq if aa in aa_to_idx]
        for i_idx in range(len(present)):
            for j_idx in range(i_idx + 1, len(present)):
                a = aa_to_idx[present[i_idx]]
                b = aa_to_idx[present[j_idx]]
                matrix[a, b] += 1
                matrix[b, a] += 1
                total_pairs += 1

    if total_pairs > 0:
        matrix /= total_pairs

    return pd.DataFrame(matrix, index=all_aa_1, columns=all_aa_1)


def compute_helix_length_composition(per_helix_data: list) -> pd.DataFrame:
    """Composição por grupo de propriedade físico-química × comprimento da hélice.

    Bins: curta (4-9), média (10-19), longa (≥20).

    Args:
        per_helix_data: Lista de dicts {'length': int, 'aa_counts': Counter}.

    Returns:
        DataFrame com colunas: Bin, Group, Frequency.
    """
    groups = {
        "Hidrofóbico": {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP"},
        "Polar":        {"SER", "THR", "CYS", "TYR", "ASN", "GLN"},
        "Carregado+":   {"LYS", "ARG", "HIS"},
        "Carregado-":   {"ASP", "GLU"},
        "Especial":     {"GLY", "PRO"},
    }

    # 1-letter sets
    def _to_1(aa_set3: set) -> set:
        return {_THREE_TO_ONE[aa] for aa in aa_set3 if aa in _THREE_TO_ONE}

    groups_1 = {g: _to_1(aas) for g, aas in groups.items()}

    bins = {
        "Curta (4-9)":    [],
        "Média (10-19)":  [],
        "Longa (≥20)":    [],
    }

    for item in per_helix_data:
        length = item["length"]
        counts = item["aa_counts"]
        total = sum(counts.values())
        if total == 0:
            continue
        if 4 <= length <= 9:
            bin_key = "Curta (4-9)"
        elif 10 <= length <= 19:
            bin_key = "Média (10-19)"
        else:
            bin_key = "Longa (≥20)"
        bins[bin_key].append((counts, total))

    rows = []
    for bin_name, items in bins.items():
        if not items:
            for group in groups_1:
                rows.append({"Bin": bin_name, "Group": group, "Frequency": 0.0})
            continue
        for group_name, aa_set in groups_1.items():
            freqs = [
                sum(counts.get(aa, 0) for aa in aa_set) / total * 100
                for counts, total in items
            ]
            rows.append({
                "Bin": bin_name,
                "Group": group_name,
                "Frequency": float(np.mean(freqs)),
            })

    return pd.DataFrame(rows)


def compute_shannon_entropy_heptad(heptad_aa_distribution: dict) -> pd.DataFrame:
    """Entropia de Shannon por posição no heptad repeat.

    Args:
        heptad_aa_distribution: dict {0..6: Counter de AAs em 1-letra}.

    Returns:
        DataFrame com colunas: Position (a-g), Entropy, Top_AA, Top_Freq.
    """
    position_labels = list("abcdefg")
    rows = []

    for pos_idx in range(7):
        counter = heptad_aa_distribution.get(pos_idx, Counter())
        total = sum(counter.values())
        if total == 0:
            rows.append({
                "Position": position_labels[pos_idx],
                "Entropy": 0.0,
                "Top_AA": "?",
                "Top_Freq": 0.0,
            })
            continue

        entropy = 0.0
        for count in counter.values():
            p = count / total
            if p > 0:
                entropy -= p * math.log2(p)

        top_aa, top_count = counter.most_common(1)[0]
        top_freq = top_count / total * 100

        rows.append({
            "Position": position_labels[pos_idx],
            "Entropy": entropy,
            "Top_AA": top_aa,
            "Top_Freq": top_freq,
        })

    return pd.DataFrame(rows)


def compute_codon_degeneracy_vs_propensity(
    helix_residues: Counter,
    all_residues: Counter,
) -> pd.DataFrame:
    """Degenerescência de códons vs propensão para hélice (teórica e observada).

    Args:
        helix_residues: Counter de AAs (3-letras) em hélices alpha.
        all_residues: Counter de AAs (3-letras) em todas as estruturas.

    Returns:
        DataFrame: AA, Codon_Degeneracy, Propensity_Theoretical,
                   Propensity_Observed, Freq_Helix, Freq_Total.
    """
    total_helix = sum(helix_residues.values()) or 1
    total_all = sum(all_residues.values()) or 1

    rows = []
    for aa3 in sorted(STANDARD_AMINO_ACIDS):
        freq_helix = helix_residues.get(aa3, 0) / total_helix
        freq_total = all_residues.get(aa3, 0) / total_all

        obs_propensity = freq_helix / freq_total if freq_total > 0 else 0.0

        rows.append({
            "AA": aa3,
            "Codon_Degeneracy": _CODON_DEGENERACY.get(aa3, 0),
            "Propensity_Theoretical": _HELIX_PROPENSITY_CF.get(aa3, 1.0),
            "Propensity_Observed": obs_propensity,
            "Freq_Helix": freq_helix * 100,
            "Freq_Total": freq_total * 100,
        })

    return pd.DataFrame(rows)


def compute_proteome_comparison(global_counter: Counter) -> pd.DataFrame:
    """Enriquecimento de AAs em relação ao proteoma humano de referência.

    Args:
        global_counter: Counter de AAs (3-letras) observados.

    Returns:
        DataFrame: AA, Freq_Observed(%), Freq_Proteome(%), Enrichment.
    """
    total = sum(global_counter.values()) or 1
    rows = []

    for aa3 in sorted(STANDARD_AMINO_ACIDS):
        freq_obs = global_counter.get(aa3, 0) / total * 100
        freq_prot = _PROTEOME_FREQ.get(aa3, 0.0)
        enrichment = freq_obs / freq_prot if freq_prot > 0 else 0.0

        rows.append({
            "AA": aa3,
            "Freq_Observed": freq_obs,
            "Freq_Proteome": freq_prot,
            "Enrichment": enrichment,
        })

    return pd.DataFrame(rows)
