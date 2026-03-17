"""
Passo DSSP unificado: coleta todos os dados de hélices em uma única passagem.

Substitui as etapas 4, 5 e 6 separadas por um único loop DSSP por estrutura,
acumulando simultaneamente todos os dados necessários para todos os gráficos.
"""
from __future__ import annotations

import math
import pickle
import warnings
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

warnings.filterwarnings("ignore", category=UserWarning, message="parse error")

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from .config import (
    CODON_DEGENERACY,
    EISENBERG_SCALE,
    HELIX_PROPENSITY,
    HELIX_TYPES,
    HYDROPHOBIC_AA,
    ONE_TO_THREE,
    PROCESS_WORKERS,
    PROTEOME_FREQ,
    STANDARD_AMINO_ACIDS,
    THREE_TO_ONE,
    glob_pdb,
)

# Importacao no nivel de modulo para evitar re-importar a cada chamada
try:
    from Bio.Data.IUPACData import protein_letters_3to1 as _BIOPYTHON_3TO1
except ImportError:
    _BIOPYTHON_3TO1 = None

_HELIX_CODES = {"H", "G", "I"}


def _res3_to_1(res3: str) -> str | None:
    """Converte codigo de 3 letras para 1 letra."""
    if _BIOPYTHON_3TO1 is not None:
        result = _BIOPYTHON_3TO1.get(res3.upper())
        if result:
            return result
    return THREE_TO_ONE.get(res3.upper())


# ── Processamento por estrutura ───────────────────────────────────────────────

def _process_single_pdb_unified(pdb_path: Path) -> dict | None:
    """Roda DSSP uma vez e coleta TODOS os dados necessarios para todos os graficos.

    Returns:
        Dict com todos os dados acumulados, ou None em caso de falha.
    """
    from Bio.PDB import PDBParser, DSSP

    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_path.stem, str(pdb_path))
        model = structure[0]
        dssp = DSSP(model, str(pdb_path), dssp="mkdssp", file_type="PDB")
    except Exception:
        return None

    dssp_records = list(dssp)
    total_residues = len(dssp_records)
    if total_residues == 0:
        return None

    # Contadores para propensidade (3-letras, compativel com STANDARD_AMINO_ACIDS)
    aa_total_3: Counter = Counter()
    aa_in_H_3: Counter = Counter()

    # Contadores para tipos (1-letra)
    aa_in_H_1: Counter = Counter()
    aa_in_G_1: Counter = Counter()
    aa_in_I_1: Counter = Counter()

    # Posicoes N/mid/C dentro de helices H
    nterm_1: Counter = Counter()
    mid_1: Counter = Counter()
    cterm_1: Counter = Counter()

    # Contadores de tipos e comprimentos
    helix_type_counts: Counter = Counter()
    helix_lengths: list[int] = []

    # Dados evolutivos
    helix_sequences: list[list[str]] = []
    ncap_position: dict[int, Counter] = {0: Counter(), 1: Counter(), 2: Counter()}
    ccap_position: dict[int, Counter] = {-1: Counter(), -2: Counter(), -3: Counter()}
    per_helix_data: list[dict] = []
    transitions: Counter = Counter()
    helix_lengths_by_type: dict[str, list] = {"H": [], "G": [], "I": []}

    # Heptad sobre TODOS os residuos (indice global)
    heptad_local: dict[int, Counter] = {i: Counter() for i in range(7)}

    # Contagem de residuos helicais para helix_content
    helix_residue_count = 0

    # --- Primeiro loop: acumula aa_total_3 e heptad para todos os residuos ---
    for global_pos, record in enumerate(dssp_records):
        ss = record[2]
        aa1 = record[1]  # DSSP retorna codigo de 1 letra diretamente
        res3 = ONE_TO_THREE.get(aa1)

        if res3 and res3 in STANDARD_AMINO_ACIDS:
            aa_total_3[res3] += 1

        if aa1 and res3:
            heptad_local[global_pos % 7][aa1] += 1

        if ss in _HELIX_CODES:
            helix_residue_count += 1

    helix_content = helix_residue_count / total_residues if total_residues > 0 else 0.0

    # --- Segundo loop: detecta segmentos helicais e acumula por segmento ---
    i = 0
    n = len(dssp_records)
    prev_type: str | None = None

    while i < n:
        ss = dssp_records[i][2]
        if ss not in _HELIX_CODES:
            prev_type = None
            i += 1
            continue

        helix_type = ss
        segment_1: list[str] = []
        segment_3: list[str] = []

        # Coleta segmento continuo do mesmo tipo
        while i < n and dssp_records[i][2] == helix_type:
            aa1 = dssp_records[i][1]  # DSSP retorna codigo de 1 letra diretamente
            res3 = ONE_TO_THREE.get(aa1)
            if aa1 and res3:
                segment_1.append(aa1)
            if res3 and res3 in STANDARD_AMINO_ACIDS:
                segment_3.append(res3)
            i += 1

        seg_len = len(segment_1)

        # Transicoes entre tipos de helice
        if prev_type is not None and prev_type != helix_type:
            transitions[(prev_type, helix_type)] += 1
        prev_type = helix_type

        # Comprimentos por tipo
        if helix_type in helix_lengths_by_type:
            helix_lengths_by_type[helix_type].append(seg_len)
        helix_lengths.append(seg_len)
        helix_type_counts[helix_type] += 1

        # Acumula por tipo (1-letra)
        if helix_type == "H":
            for aa1 in segment_1:
                aa_in_H_1[aa1] += 1
            for aa3 in segment_3:
                aa_in_H_3[aa3] += 1
        elif helix_type == "G":
            for aa1 in segment_1:
                aa_in_G_1[aa1] += 1
        elif helix_type == "I":
            for aa1 in segment_1:
                aa_in_I_1[aa1] += 1

        # Per helix data (todos os tipos)
        if seg_len > 0:
            per_helix_data.append({"length": seg_len, "aa_counts": Counter(segment_1)})

        # Analise especifica de helices H
        if helix_type == "H":
            # Posicoes N/mid/C-terminal (tercos)
            third = max(1, seg_len // 3)
            for j, aa1 in enumerate(segment_1):
                if j < third:
                    nterm_1[aa1] += 1
                elif j >= seg_len - third:
                    cterm_1[aa1] += 1
                else:
                    mid_1[aa1] += 1

            # N-cap e C-cap (apenas H >= 4)
            if seg_len >= 4:
                helix_sequences.append(segment_1)
                for pos in range(min(3, seg_len)):
                    ncap_position[pos][segment_1[pos]] += 1
                for offset, pos_key in enumerate([-1, -2, -3]):
                    idx = seg_len - 1 - offset
                    if idx >= 0:
                        ccap_position[pos_key][segment_1[idx]] += 1

    return {
        # Propensidade (equivalente etapa 4)
        "aa_total_3": aa_total_3,
        "aa_in_H_3": aa_in_H_3,
        # Posicoes
        "nterm_1": nterm_1,
        "mid_1": mid_1,
        "cterm_1": cterm_1,
        # Tipos (equivalente etapa 5)
        "aa_in_H_1": aa_in_H_1,
        "aa_in_G_1": aa_in_G_1,
        "aa_in_I_1": aa_in_I_1,
        "helix_type_counts": helix_type_counts,
        "helix_lengths": helix_lengths,
        # Dados evolutivos (equivalente etapa 6)
        "helix_content": helix_content,
        "heptad_local": heptad_local,
        "helix_sequences": helix_sequences,
        "ncap_position": ncap_position,
        "ccap_position": ccap_position,
        "per_helix_data": per_helix_data,
        "transitions": transitions,
        "helix_lengths_by_type": helix_lengths_by_type,
    }


# ── Agregacao paralela ────────────────────────────────────────────────────────

def collect_all(clean_dir: Path, workers: int) -> dict:
    """Processa todas as estruturas em paralelo e agrega os resultados.

    Args:
        clean_dir: Diretorio com PDBs limpos.
        workers: Numero de threads paralelas.

    Returns:
        Dict agregado com todos os dados coletados.
    """
    pdb_files = sorted(glob_pdb(clean_dir))
    if not pdb_files:
        raise FileNotFoundError(f"Nenhum .pdb encontrado em {clean_dir}")

    print(f"  Processando {len(pdb_files):,} estruturas (workers: {workers})...")

    # Acumuladores agregados
    agg: dict = {
        "aa_total_3": Counter(),
        "aa_in_H_3": Counter(),
        "nterm_1": Counter(),
        "mid_1": Counter(),
        "cterm_1": Counter(),
        "aa_in_H_1": Counter(),
        "aa_in_G_1": Counter(),
        "aa_in_I_1": Counter(),
        "helix_type_counts": Counter(),
        "helix_lengths": [],
        "helix_content_per_structure": [],
        "heptad_aa_distribution": {i: Counter() for i in range(7)},
        "helix_sequences": [],
        "ncap_position": {0: Counter(), 1: Counter(), 2: Counter()},
        "ccap_position": {-1: Counter(), -2: Counter(), -3: Counter()},
        "per_helix_data": [],
        "transitions": Counter(),
        "helix_lengths_by_type": {"H": [], "G": [], "I": []},
        "ncap_residues": Counter(),
        "ccap_residues": Counter(),
        "structures_processed": 0,
        "structures_failed": 0,
    }

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(_process_single_pdb_unified, f): f for f in pdb_files}
        with tqdm(total=len(pdb_files), desc="DSSP unificado", unit="pdb") as pbar:
            for future in as_completed(futures):
                result = future.result()
                if result is None:
                    agg["structures_failed"] += 1
                    pbar.update(1)
                    continue

                agg["structures_processed"] += 1
                agg["aa_total_3"].update(result["aa_total_3"])
                agg["aa_in_H_3"].update(result["aa_in_H_3"])
                agg["nterm_1"].update(result["nterm_1"])
                agg["mid_1"].update(result["mid_1"])
                agg["cterm_1"].update(result["cterm_1"])
                agg["aa_in_H_1"].update(result["aa_in_H_1"])
                agg["aa_in_G_1"].update(result["aa_in_G_1"])
                agg["aa_in_I_1"].update(result["aa_in_I_1"])
                agg["helix_type_counts"].update(result["helix_type_counts"])
                agg["helix_lengths"].extend(result["helix_lengths"])
                agg["helix_content_per_structure"].append(result["helix_content"])

                for pos, cnt in result["heptad_local"].items():
                    agg["heptad_aa_distribution"][pos].update(cnt)

                agg["helix_sequences"].extend(result["helix_sequences"])

                for pos, cnt in result["ncap_position"].items():
                    agg["ncap_position"][pos].update(cnt)
                    agg["ncap_residues"].update(cnt)
                for pos, cnt in result["ccap_position"].items():
                    agg["ccap_position"][pos].update(cnt)
                    agg["ccap_residues"].update(cnt)

                agg["per_helix_data"].extend(result["per_helix_data"])
                agg["transitions"].update(result["transitions"])

                for ht, lens in result["helix_lengths_by_type"].items():
                    agg["helix_lengths_by_type"][ht].extend(lens)

                pbar.update(1)

    print(
        f"  Concluido: {agg['structures_processed']:,} sucessos, "
        f"{agg['structures_failed']:,} falhas"
    )
    return agg


# ── DataFrames derivados ──────────────────────────────────────────────────────

def _compute_propensity_df(agg: dict) -> pd.DataFrame:
    """Calcula propensao de helice por aminoacido.

    Columns: AA, Freq_Total (%), Freq_Helix (%), Propensity_Observed, Propensity_Theoretical
    """
    total_all = sum(agg["aa_total_3"].values()) or 1
    total_helix = sum(agg["aa_in_H_3"].values()) or 1

    rows = []
    for aa3 in sorted(STANDARD_AMINO_ACIDS):
        aa1 = THREE_TO_ONE.get(aa3)
        if aa1 is None:
            continue
        freq_total = agg["aa_total_3"].get(aa3, 0) / total_all * 100
        freq_helix = agg["aa_in_H_3"].get(aa3, 0) / total_helix * 100
        prop_obs = freq_helix / freq_total if freq_total > 0 else 0.0
        rows.append({
            "AA": aa1,
            "Freq_Total": freq_total,
            "Freq_Helix": freq_helix,
            "Propensity_Observed": prop_obs,
            "Propensity_Theoretical": HELIX_PROPENSITY.get(aa3, 1.0),
        })

    df = pd.DataFrame(rows)
    return df.sort_values("Propensity_Observed", ascending=False).reset_index(drop=True)


def _compute_positions_df(agg: dict) -> pd.DataFrame:
    """Frequencia de aminoacidos por posicao na helice (N-terminal/Middle/C-terminal).

    Columns: AA, Position, Frequency (%)
    """
    position_map = {
        "N-terminal": agg["nterm_1"],
        "Middle": agg["mid_1"],
        "C-terminal": agg["cterm_1"],
    }
    rows = []
    for pos_name, counter in position_map.items():
        total = sum(counter.values()) or 1
        for aa3 in sorted(STANDARD_AMINO_ACIDS):
            aa1 = THREE_TO_ONE.get(aa3)
            if aa1 is None:
                continue
            count = counter.get(aa1, 0)
            rows.append({
                "AA": aa1,
                "Position": pos_name,
                "Frequency": count / total * 100,
            })
    return pd.DataFrame(rows)


def _compute_type_dfs(agg: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Composicao de aminoacidos por tipo de helice (H/G/I).

    Returns:
        (composition_df, top_residues_df, stat_df)
    """
    type_counters = {
        "H": agg["aa_in_H_1"],
        "G": agg["aa_in_G_1"],
        "I": agg["aa_in_I_1"],
    }
    helix_names = {
        "H": "Alpha helix",
        "G": "3-10 helix",
        "I": "Pi helix",
    }

    comp_rows = []
    for h_type, counter in type_counters.items():
        total = sum(counter.values()) or 1
        h_name = helix_names[h_type]
        for aa3 in sorted(STANDARD_AMINO_ACIDS):
            aa1 = THREE_TO_ONE.get(aa3)
            if aa1 is None:
                continue
            count = counter.get(aa1, 0)
            comp_rows.append({
                "AA": aa1,
                "Helix_Type": h_type,
                "Helix_Name": h_name,
                "Count": count,
                "Frequency": count / total * 100,
            })

    composition_df = pd.DataFrame(comp_rows)

    def _add_rank(df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        df["Rank"] = df.groupby("Helix_Type")["Frequency"].rank(
            ascending=False, method="first"
        ).astype(int)
        return df

    composition_df = _add_rank(composition_df)

    # Top 10 por tipo
    top_residues_df = (
        composition_df[composition_df["Rank"] <= 10]
        .sort_values(["Helix_Type", "Rank"])
        .reset_index(drop=True)
    )

    # Comparacao estatistica H vs G
    types_present = composition_df["Helix_Type"].unique()
    stat_rows = []
    if "H" in types_present and "G" in types_present:
        for aa3 in sorted(STANDARD_AMINO_ACIDS):
            aa1 = THREE_TO_ONE.get(aa3)
            if aa1 is None:
                continue
            h_data = composition_df[
                (composition_df["AA"] == aa1) & (composition_df["Helix_Type"] == "H")
            ]["Frequency"].values
            g_data = composition_df[
                (composition_df["AA"] == aa1) & (composition_df["Helix_Type"] == "G")
            ]["Frequency"].values
            if h_data.size == 0 or g_data.size == 0:
                continue
            diff = float(h_data[0] - g_data[0])
            stat_rows.append({
                "AA": aa1,
                "Alpha_Freq": float(h_data[0]),
                "3-10_Freq": float(g_data[0]),
                "Difference": diff,
            })

    stat_df = pd.DataFrame(stat_rows)
    if not stat_df.empty:
        stat_df = stat_df.sort_values("Difference", key=abs, ascending=False).reset_index(drop=True)

    return composition_df, top_residues_df, stat_df


def _compute_heptad_df(agg: dict) -> pd.DataFrame:
    """Entropia de Shannon e AA mais frequente por posicao no heptad.

    Columns: Position (a-g), Entropy, Top_AA, Top_Freq
    """
    position_labels = list("abcdefg")
    rows = []
    for pos_idx in range(7):
        counter = agg["heptad_aa_distribution"].get(pos_idx, Counter())
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
        rows.append({
            "Position": position_labels[pos_idx],
            "Entropy": entropy,
            "Top_AA": top_aa,
            "Top_Freq": top_count / total * 100,
        })
    return pd.DataFrame(rows)


def _compute_hydrophobic_moments(helix_sequences: list) -> list[float]:
    """Calcula momento hidrofobico de Eisenberg para cada helice."""
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


def _compute_length_vs_composition(per_helix_data: list) -> pd.DataFrame:
    """Composicao por grupo fisico-quimico vs comprimento de helice."""
    groups_1 = {
        "Hidrofobico": {"A", "V", "I", "L", "M", "F", "W", "P"},
        "Polar":       {"S", "T", "C", "Y", "N", "Q"},
        "Carregado+":  {"K", "R", "H"},
        "Carregado-":  {"D", "E"},
        "Especial":    {"G", "P"},
    }

    bins: dict[str, list] = {
        "Curta (4-9)":   [],
        "Media (10-19)": [],
        "Longa (>=20)":  [],
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
            bin_key = "Media (10-19)"
        else:
            bin_key = "Longa (>=20)"
        bins[bin_key].append((counts, total))

    rows = []
    for bin_name, items in bins.items():
        if not items:
            for group in groups_1:
                rows.append({"Group": group, "Bin": bin_name, "Frequency": 0.0})
            continue
        for group_name, aa_set in groups_1.items():
            freqs = [
                sum(counts.get(aa, 0) for aa in aa_set) / total * 100
                for counts, total in items
            ]
            rows.append({
                "Group": group_name,
                "Bin": bin_name,
                "Frequency": float(np.mean(freqs)),
            })
    return pd.DataFrame(rows)


def _compute_codon_df(propensity_df: pd.DataFrame) -> pd.DataFrame:
    """Degenerescencia de codons vs propensao para helice."""
    rows = []
    for aa3, deg in sorted(CODON_DEGENERACY.items()):
        aa1 = THREE_TO_ONE.get(aa3)
        if aa1 is None:
            continue
        row_mask = propensity_df["AA"] == aa1
        if row_mask.any():
            prop_obs = float(propensity_df.loc[row_mask, "Propensity_Observed"].iloc[0])
            prop_theo = float(propensity_df.loc[row_mask, "Propensity_Theoretical"].iloc[0])
        else:
            prop_obs = 0.0
            prop_theo = HELIX_PROPENSITY.get(aa3, 1.0)
        rows.append({
            "AA": aa1,
            "Codon_Degeneracy": deg,
            "Propensity_Theoretical": prop_theo,
            "Propensity_Observed": prop_obs,
        })
    return pd.DataFrame(rows)


def _compute_proteome_df(propensity_df: pd.DataFrame) -> pd.DataFrame:
    """Enriquecimento de AAs em relacao ao proteoma humano."""
    rows = []
    for aa3, freq_prot in sorted(PROTEOME_FREQ.items()):
        aa1 = THREE_TO_ONE.get(aa3)
        if aa1 is None:
            continue
        row_mask = propensity_df["AA"] == aa1
        if row_mask.any():
            freq_obs = float(propensity_df.loc[row_mask, "Freq_Helix"].iloc[0])
        else:
            freq_obs = 0.0
        enrichment = freq_obs / freq_prot if freq_prot > 0 else 0.0
        rows.append({
            "AA": aa1,
            "Freq_Proteome": freq_prot,
            "Freq_Observed": freq_obs,
            "Enrichment": enrichment,
        })
    return pd.DataFrame(rows)


# ── Persistencia ──────────────────────────────────────────────────────────────

def save_all_results(agg: dict, analysis_path: Path) -> dict:
    """Calcula todos os DataFrames derivados e salva em disco.

    Args:
        agg: Dicionario agregado retornado por collect_all().
        analysis_path: Diretorio de saida.

    Returns:
        Dict com todos os DataFrames calculados.
    """
    analysis_path.mkdir(parents=True, exist_ok=True)

    print("\n  Calculando DataFrames...")

    propensity_df = _compute_propensity_df(agg)
    positions_df = _compute_positions_df(agg)
    composition_df, top_residues_df, stat_df = _compute_type_dfs(agg)
    heptad_df = _compute_heptad_df(agg)
    moments = _compute_hydrophobic_moments(agg["helix_sequences"])
    length_comp_df = _compute_length_vs_composition(agg["per_helix_data"])
    codon_df = _compute_codon_df(propensity_df)
    proteome_df = _compute_proteome_df(propensity_df)

    print("  Salvando CSVs...")

    def _save(df_or_list, filename: str, col: str | None = None) -> None:
        path = analysis_path / filename
        if isinstance(df_or_list, pd.DataFrame):
            df_or_list.to_csv(path, index=False)
        else:
            pd.DataFrame({col: df_or_list}).to_csv(path, index=False)
        print(f"    -> {filename}")

    _save(propensity_df, "helix_propensities.csv")
    _save(positions_df, "helix_positions.csv")
    _save(agg["helix_lengths"], "helix_lengths.csv", "Length")
    _save(composition_df, "helix_type_composition.csv")
    _save(top_residues_df, "helix_type_top_residues.csv")
    if not stat_df.empty:
        _save(stat_df, "helix_type_statistical_comparison.csv")
    _save(agg["helix_content_per_structure"], "helix_content_per_structure.csv", "Helix_Content")

    # helix_lengths_by_type em formato longo
    lbt_rows = []
    for ht, lens in agg["helix_lengths_by_type"].items():
        for l in lens:
            lbt_rows.append({"Helix_Type": ht, "Length": l})
    lbt_df = pd.DataFrame(lbt_rows)
    lbt_df.to_csv(analysis_path / "helix_lengths_by_type.csv", index=False)
    print("    -> helix_lengths_by_type.csv")

    _save(heptad_df, "shannon_entropy_heptad.csv")
    _save(moments, "hydrophobic_moments.csv", "Moment")
    _save(length_comp_df, "helix_length_vs_composition.csv")
    _save(codon_df, "codon_degeneracy.csv")
    _save(proteome_df, "proteome_comparison.csv")

    # evo_data.pkl (compativel com plots.ipynb)
    evo_data = {
        "helix_sequences": agg["helix_sequences"],
        "ncap_position": agg["ncap_position"],
        "ccap_position": agg["ccap_position"],
        "helix_content_per_structure": agg["helix_content_per_structure"],
        "per_helix_data": agg["per_helix_data"],
        "heptad_aa_distribution": agg["heptad_aa_distribution"],
        "transitions": agg["transitions"],
        "helix_lengths_by_type": agg["helix_lengths_by_type"],
        "ncap_residues": agg["ncap_residues"],
        "ccap_residues": agg["ccap_residues"],
    }
    pkl_path = analysis_path / "evo_data.pkl"
    with open(pkl_path, "wb") as f:
        pickle.dump(evo_data, f)
    print("    -> evo_data.pkl")

    # Relatorio texto
    _write_unified_report(agg, propensity_df, analysis_path)

    return {
        "propensity_df": propensity_df,
        "positions_df": positions_df,
        "composition_df": composition_df,
        "top_residues_df": top_residues_df,
        "stat_df": stat_df,
        "heptad_df": heptad_df,
        "evo_data": evo_data,
    }


def _write_unified_report(agg: dict, propensity_df: pd.DataFrame, analysis_path: Path) -> None:
    """Escreve relatorio de texto com estatisticas gerais."""
    report_path = analysis_path / "unified_report.txt"
    helix_counts = agg["helix_type_counts"]
    total_helices = sum(helix_counts.values())
    mean_content = (
        float(np.mean(agg["helix_content_per_structure"]))
        if agg["helix_content_per_structure"]
        else 0.0
    )

    top5 = propensity_df.head(5)

    with report_path.open("w") as f:
        f.write("=" * 60 + "\n")
        f.write("RELATORIO DSSP UNIFICADO\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Gerado em: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"Estruturas processadas: {agg['structures_processed']:,}\n")
        f.write(f"Estruturas com falha:   {agg['structures_failed']:,}\n")
        f.write(f"Total de residuos:      {sum(agg['aa_total_3'].values()):,}\n\n")
        f.write("### HELICES ###\n\n")
        f.write(f"  Total de helices:     {total_helices:,}\n")
        f.write(f"  Alpha (H):            {helix_counts.get('H', 0):,}\n")
        f.write(f"  3-10 (G):             {helix_counts.get('G', 0):,}\n")
        f.write(f"  Pi (I):               {helix_counts.get('I', 0):,}\n")
        f.write(f"  Conteudo medio (%):   {mean_content * 100:.1f}%\n\n")
        f.write("### TOP 5 PROPENSAO (Observada) ###\n\n")
        for _, row in top5.iterrows():
            f.write(f"  {row['AA']}: {row['Propensity_Observed']:.3f}\n")
        f.write("\n" + "=" * 60 + "\n")

    print("    -> unified_report.txt")


# ── Ponto de entrada ──────────────────────────────────────────────────────────

def run_unified_dssp(clean_dir: Path, analysis_path: Path, workers: int) -> dict:
    """Ponto de entrada principal do passo DSSP unificado.

    Args:
        clean_dir: Diretorio com PDBs limpos.
        analysis_path: Diretorio de saida para CSVs/pkl.
        workers: Numero de threads paralelas.

    Returns:
        Dict com todos os DataFrames calculados (vazio se ja existiam).
    """
    report_file = analysis_path / "unified_report.txt"
    if report_file.exists():
        # Deteccao de existencia tratada em main.py; aqui apenas retorna vazio
        return {}

    agg = collect_all(clean_dir, workers)
    results = save_all_results(agg, analysis_path)
    return results
