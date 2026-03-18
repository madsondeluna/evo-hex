"""
Funcoes de analise evolutiva para helices alpha.

Todas as funcoes compute_* recebem dados ja coletados pelo passo DSSP
unificado (unified_dssp.py) e produzem DataFrames para os graficos.
"""

import math
from collections import Counter

import numpy as np
import pandas as pd

from .config import (
    CODON_DEGENERACY,
    EISENBERG_SCALE,
    HELIX_PROPENSITY,
    HYDROPHOBIC_AA,
    ONE_TO_THREE,
    PROTEOME_FREQ,
    STANDARD_AMINO_ACIDS,
    THREE_TO_ONE,
)


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
            aa3 = ONE_TO_THREE.get(aa1)
            if aa3 is None or aa3 not in EISENBERG_SCALE:
                continue
            h = EISENBERG_SCALE[aa3]
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
    all_aa_1 = sorted(THREE_TO_ONE.values())
    aa_to_idx = {aa: i for i, aa in enumerate(all_aa_1)}
    n = len(all_aa_1)
    matrix = np.zeros((n, n), dtype=float)

    total_pairs = 0
    for seq in helix_sequences:
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
        "Hydrophobic": {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP"},
        "Polar":       {"SER", "THR", "CYS", "TYR", "ASN", "GLN"},
        "Charged+":    {"LYS", "ARG", "HIS"},
        "Charged-":    {"ASP", "GLU"},
        "Special":     {"GLY", "PRO"},
    }

    groups_1 = {
        g: {THREE_TO_ONE[aa] for aa in aas if aa in THREE_TO_ONE}
        for g, aas in groups.items()
    }

    bins = {
        "Short (4–9)":    [],
        "Medium (10–19)": [],
        "Long (≥20)":     [],
    }

    for item in per_helix_data:
        length = item["length"]
        counts = item["aa_counts"]
        total = sum(counts.values())
        if total == 0:
            continue
        if 4 <= length <= 9:
            bin_key = "Short (4–9)"
        elif 10 <= length <= 19:
            bin_key = "Medium (10–19)"
        else:
            bin_key = "Long (≥20)"
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
            "Codon_Degeneracy": CODON_DEGENERACY.get(aa3, 0),
            "Propensity_Theoretical": HELIX_PROPENSITY.get(aa3, 1.0),
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
        freq_prot = PROTEOME_FREQ.get(aa3, 0.0)
        enrichment = freq_obs / freq_prot if freq_prot > 0 else 0.0

        rows.append({
            "AA": aa3,
            "Freq_Observed": freq_obs,
            "Freq_Proteome": freq_prot,
            "Enrichment": enrichment,
        })

    return pd.DataFrame(rows)
