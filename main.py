#!/usr/bin/env python3
"""
Pipeline Evo-Hex - ponto de entrada principal.

Etapas
------
1. Download da lista de dominios CATH e estruturas PDB (filtro S40 por padrao)
2. Limpeza (cadeia A, sem heteroatomos)
3. Analise de frequencia de aminoacidos
4. Passo DSSP unificado (coleta dados para etapas 4, 5 e 6 de uma vez)
5. Tipos de helices (H/G/I) - dados ja coletados na etapa 4
6. Dados evolutivos - dados ja coletados na etapa 4

Uso
---
    python main.py                      # todas as etapas, modo interativo
    python main.py --explore            # pula download/limpeza, vai direto a analise
    python main.py --steps 3 4 5 6     # mesmo que --explore (manual)
    python main.py --steps 1 2         # apenas download e limpeza
    python main.py --force-download    # baixa mesmo se dados existirem
    python main.py --skip-download     # pula etapas 1-2 sem perguntar
    python main.py --no-interactive    # sem menus, executa tudo automaticamente
    python main.py --plots g3 g4       # pre-seleciona grupos no menu de graficos
    python main.py --list-plots        # lista todos os graficos disponiveis
    python main.py --full-dataset      # baixa dataset completo (~50k) sem filtro S40
"""

import argparse
import logging
import sys
import traceback
from pathlib import Path

from cath_analysis.config import (
    ANALYSIS_PATH,
    BASE_PATH,
    LOGS_PATH,
    STRUCTURES_CLEAN_PATH,
    STRUCTURES_PATH,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


# ── Etapas do pipeline ────────────────────────────────────────────────────────

def _show_existing_data_info(info: dict) -> None:
    """Imprime resumo dos dados ja existentes no disco."""
    print("\n  +- Dados existentes ---------------------------------------------------+")
    if info["has_cath_list"]:
        print(f"  |  Indice CATH:       atualizado em {info['cath_list_date']}")
    else:
        print("  |  Indice CATH:       nao encontrado")
    if info["raw_count"]:
        print(f"  |  Estruturas raw:    {info['raw_count']:,} PDBs  (ultimo: {info['raw_date']})")
    else:
        print("  |  Estruturas raw:    nenhuma")
    if info["clean_count"]:
        print(f"  |  Estruturas clean:  {info['clean_count']:,} PDBs  (ultimo: {info['clean_date']})")
    else:
        print("  |  Estruturas clean:  nenhuma")
    print("  +----------------------------------------------------------------------+")


def _ask_redownload() -> bool:
    """Pergunta ao usuario se deseja baixar novamente. Retorna True se sim."""
    while True:
        resp = input("\n  Baixar novamente? [s/N] ").strip().lower()
        if resp in ("s", "sim", "y", "yes"):
            return True
        if resp in ("", "n", "nao", "no"):
            return False
        print("  Responda 's' para sim ou 'n' para nao.")


def _ask_rerun(label: str) -> bool:
    """Pergunta ao usuario se deseja re-executar uma etapa. Retorna True se sim."""
    while True:
        resp = input(f"\n  Re-executar {label}? [s/N] ").strip().lower()
        if resp in ("s", "sim", "y", "yes"):
            return True
        if resp in ("", "n", "nao", "no"):
            return False
        print("  Responda 's' para sim ou 'n' para nao.")


def step1_download(force: bool = False, interactive: bool = True, full_dataset: bool = False) -> None:
    """Etapa 1: Download de estruturas PDB mainly-alpha.

    Por padrao aplica filtro S40 (~2-5k estruturas).
    Use full_dataset=True para baixar o conjunto completo (~50k).
    """
    from cath_analysis.downloader import (
        check_existing_data,
        download_cath_domain_list,
        download_s40_list,
        download_structures_parallel,
        parse_cath_domains,
        print_download_report,
        save_download_results,
        setup_directories,
    )

    print("\n" + "=" * 60)
    print("ETAPA 1 - DOWNLOAD DE ESTRUTURAS")
    if not full_dataset:
        print("  [Modo: S40 nao-redundante]")
    else:
        print("  [Modo: dataset completo]")
    print("=" * 60)

    info = check_existing_data()
    has_data = info["raw_count"] > 0

    if has_data and not force:
        _show_existing_data_info(info)
        if interactive:
            if not _ask_redownload():
                print("  Download ignorado. Usando dados existentes.")
                return
        else:
            print("\n  [INFO] Dados existentes encontrados. Use --force-download para re-baixar.")
            return

    print()
    base_path, structures_path, logs_path = setup_directories()
    cath_file = download_cath_domain_list(base_path)

    s40_codes = None if full_dataset else download_s40_list(base_path)
    pdb_codes = parse_cath_domains(cath_file, s40_codes=s40_codes)

    results = download_structures_parallel(pdb_codes, structures_path)
    save_download_results(base_path, pdb_codes, results)
    print_download_report(results, structures_path)


def step2_clean(force: bool = False, interactive: bool = True) -> None:
    """Etapa 2: Limpeza das estruturas PDB."""
    from cath_analysis.downloader import check_existing_data
    from cath_analysis.cleaner import (
        clean_all_structures,
        print_cleaning_summary,
        save_cleaning_results,
    )

    print("\n" + "=" * 60)
    print("ETAPA 2 - LIMPEZA DE ESTRUTURAS")
    print("=" * 60)

    info = check_existing_data()
    has_clean = info["clean_count"] > 0

    if has_clean and not force:
        _show_existing_data_info(info)
        if interactive:
            print(f"\n  Ja existem {info['clean_count']:,} estruturas limpas.")
            if not _ask_redownload():
                print("  Limpeza ignorada. Usando estruturas existentes.")
                return
        else:
            print(f"\n  [INFO] {info['clean_count']:,} estruturas limpas ja existem. Use --force-download para reprocessar.")
            return

    print()
    STRUCTURES_CLEAN_PATH.mkdir(parents=True, exist_ok=True)
    clean_results = clean_all_structures(STRUCTURES_PATH, STRUCTURES_CLEAN_PATH)
    if clean_results is not None:
        save_cleaning_results(LOGS_PATH, clean_results)
        print_cleaning_summary(clean_results, STRUCTURES_PATH, STRUCTURES_CLEAN_PATH)


def step3_frequency(cache: dict, interactive: bool = True) -> None:
    """Etapa 3: Analise de frequencia de aminoacidos."""
    from collections import Counter
    from cath_analysis.frequency_analysis import (
        analyze_amino_acid_frequency,
        derive_frequencies_from_unified,
        print_frequency_summary,
        save_frequency_results,
    )

    print("\n" + "=" * 60)
    print("ETAPA 3 - FREQUENCIA DE AMINOACIDOS")
    print("=" * 60 + "\n")

    freq_csv = ANALYSIS_PATH / "amino_acid_frequencies_per_structure.csv"
    if freq_csv.exists():
        mtime = freq_csv.stat().st_mtime
        from datetime import datetime
        date_str = datetime.fromtimestamp(mtime).strftime("%Y-%m-%d %H:%M")
        print(f"  Resultados existentes encontrados (gerados em {date_str}).")
        if interactive:
            if not _ask_rerun("analise de frequencia"):
                print("  Carregando resultados anteriores...")
                global_txt = ANALYSIS_PATH / "amino_acid_frequencies_global.txt"
                global_counter: Counter = Counter()
                if global_txt.exists():
                    for line in global_txt.read_text().splitlines():
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            try:
                                global_counter[parts[0]] = int(parts[1])
                            except ValueError:
                                pass
                cache["global_counter"] = global_counter
                cache["per_structure"] = {}
                return
        else:
            print("  [INFO] Use --force-download para re-executar.")
            cache["global_counter"] = Counter()
            cache["per_structure"] = {}
            return

    # Se o passo DSSP unificado ja foi executado, deriva frequencias instantaneamente
    unified_csv = ANALYSIS_PATH / "helix_propensities.csv"
    if unified_csv.exists():
        print("  Derivando frequencias do passo DSSP unificado...")
        global_counter, per_structure = derive_frequencies_from_unified(ANALYSIS_PATH)
        ANALYSIS_PATH.mkdir(parents=True, exist_ok=True)
        save_frequency_results(ANALYSIS_PATH, global_counter, per_structure)
        print_frequency_summary(global_counter, per_structure)
        cache["global_counter"] = global_counter
        cache["per_structure"] = per_structure
        return

    # Fallback: parseia todos os PDBs
    global_counter, per_structure = analyze_amino_acid_frequency(STRUCTURES_CLEAN_PATH)
    ANALYSIS_PATH.mkdir(parents=True, exist_ok=True)
    save_frequency_results(ANALYSIS_PATH, global_counter, per_structure)
    print_frequency_summary(global_counter, per_structure)

    cache["global_counter"] = global_counter
    cache["per_structure"] = per_structure


def step4_helix_analysis(cache: dict, interactive: bool = True) -> None:
    """Etapa 4: Passo DSSP unificado (coleta dados para etapas 4, 5 e 6 de uma vez)."""
    from cath_analysis.unified_dssp import run_unified_dssp
    from cath_analysis.config import PROCESS_WORKERS

    print("\n" + "=" * 60)
    print("ETAPA 4 - ANALISE DSSP UNIFICADA (etapas 4+5+6 em uma passagem)")
    print("=" * 60 + "\n")

    report_file = ANALYSIS_PATH / "unified_report.txt"
    if report_file.exists():
        from datetime import datetime
        date_str = datetime.fromtimestamp(report_file.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
        print(f"  Resultados existentes encontrados (gerados em {date_str}).")
        if interactive:
            if not _ask_rerun("analise DSSP unificada"):
                print("  Carregando resultados anteriores dos CSVs...")
                _load_unified_results(cache)
                return
        else:
            print("  [INFO] Use --force-download para re-executar.")
            _load_unified_results(cache)
            return

    results = run_unified_dssp(STRUCTURES_CLEAN_PATH, ANALYSIS_PATH, PROCESS_WORKERS)
    cache.update(results)


def _load_unified_results(cache: dict) -> None:
    """Carrega resultados do passo unificado dos arquivos salvos."""
    import pandas as pd
    import pickle

    def _try_csv(name: str):
        p = ANALYSIS_PATH / name
        return pd.read_csv(p) if p.exists() else None

    cache["propensity_df"]   = _try_csv("helix_propensities.csv")
    cache["positions_df"]    = _try_csv("helix_positions.csv")
    cache["composition_df"]  = _try_csv("helix_type_composition.csv")
    cache["top_residues_df"] = _try_csv("helix_type_top_residues.csv")
    cache["stat_df"]         = _try_csv("helix_type_statistical_comparison.csv")

    pkl = ANALYSIS_PATH / "evo_data.pkl"
    if pkl.exists():
        with open(pkl, "rb") as f:
            cache["evo_data"] = pickle.load(f)


def step5_helix_types(cache: dict, interactive: bool = True) -> None:
    """Etapa 5: No passo unificado, dados ja coletados na etapa 4."""
    print("\n" + "=" * 60)
    print("ETAPA 5 - TIPOS DE HELICES (dados do passo unificado)")
    print("=" * 60 + "\n")

    if cache.get("composition_df") is not None:
        print("  Dados de tipos de helices ja disponiveis (coletados na etapa 4).")
        return

    # Se etapa 4 foi pulada, tenta carregar dos CSVs
    _load_unified_results(cache)
    if cache.get("composition_df") is not None:
        print("  Dados carregados dos CSVs existentes.")
    else:
        print("  [AVISO] Dados nao encontrados. Execute a etapa 4 primeiro.")


def step6_evolutionary(cache: dict) -> None:
    """Etapa 6: No passo unificado, dados ja coletados na etapa 4."""
    print("\n" + "=" * 60)
    print("ETAPA 6 - ANALISE EVOLUTIVA (dados do passo unificado)")
    print("=" * 60 + "\n")

    if cache.get("evo_data") is not None:
        evo = cache["evo_data"]
        print(f"  Estruturas processadas: {len(evo.get('helix_content_per_structure', [])):,}")
        print(f"  Helices alpha (>=4 res): {len(evo.get('helix_sequences', [])):,}")
        print(f"  Pares de transicao:     {sum(evo.get('transitions', {}).values()):,}")
        return

    _load_unified_results(cache)
    if cache.get("evo_data") is not None:
        print("  Dados evolutivos carregados do pickle.")
    else:
        print("  [AVISO] evo_data.pkl nao encontrado. Execute a etapa 4 primeiro.")


# ── Execucao de graficos ──────────────────────────────────────────────────────

def run_plots(selected_ids: list, cache: dict, analysis_path: Path) -> None:
    """Despacha cada plot_id para a funcao de plotagem correspondente.

    Cada helper le de (1) cache flat ou (2) CSV em analysis_path.

    Args:
        selected_ids: Lista de plot_ids (strings) a plotar.
        cache: Dicionario com dados acumulados pelas etapas.
        analysis_path: Diretorio de saida para os graficos.
    """
    import pandas as pd
    from cath_analysis import plotting as P
    from cath_analysis import evolutionary_analysis as E

    analysis_path.mkdir(parents=True, exist_ok=True)

    evo_data: dict = cache.get("evo_data", {})
    global_counter = cache.get("global_counter")
    per_structure: dict = cache.get("per_structure") or {}

    def _get_propensity_df():
        if cache.get("propensity_df") is not None:
            return cache["propensity_df"]
        csv_path = analysis_path / "helix_propensities.csv"
        return pd.read_csv(csv_path) if csv_path.exists() else None

    def _get_positions_df():
        if cache.get("positions_df") is not None:
            return cache["positions_df"]
        csv_path = analysis_path / "helix_positions.csv"
        return pd.read_csv(csv_path) if csv_path.exists() else None

    def _get_helix_lengths():
        if cache.get("helix_lengths"):
            return cache["helix_lengths"]
        csv_path = analysis_path / "helix_lengths.csv"
        if csv_path.exists():
            df = pd.read_csv(csv_path)
            return df["Length"].tolist() if "Length" in df.columns else []
        return None

    def _get_heptad():
        if evo_data:
            from cath_analysis.config import HYDROPHOBIC_AA
            heptad_dist = evo_data.get("heptad_aa_distribution", {})
            heptad = {}
            for pos in range(7):
                cnts = heptad_dist.get(pos, {})
                total = sum(cnts.values())
                hydro = sum(v for aa, v in cnts.items() if aa in HYDROPHOBIC_AA)
                heptad[pos] = {"hydrophobic": hydro, "total": total}
            return heptad
        return None

    def _get_composition_df():
        if cache.get("composition_df") is not None:
            return cache["composition_df"]
        csv_path = analysis_path / "helix_type_composition.csv"
        return pd.read_csv(csv_path) if csv_path.exists() else None

    def _get_top_residues_df():
        if cache.get("top_residues_df") is not None:
            return cache["top_residues_df"]
        csv_path = analysis_path / "helix_type_top_residues.csv"
        return pd.read_csv(csv_path) if csv_path.exists() else None

    def _get_stat_df():
        if cache.get("stat_df") is not None:
            return cache["stat_df"]
        csv_path = analysis_path / "helix_type_statistical_comparison.csv"
        return pd.read_csv(csv_path) if csv_path.exists() else None

    def _get_helix_counts():
        comp_df = _get_composition_df()
        if comp_df is not None and not comp_df.empty:
            return comp_df.groupby("Helix_Type")["Count"].sum().to_dict()
        return None

    def _get_lengths_by_type():
        if evo_data:
            lbt = evo_data.get("helix_lengths_by_type")
            if lbt:
                return lbt
        csv_path = analysis_path / "helix_lengths_by_type.csv"
        if csv_path.exists():
            df = pd.read_csv(csv_path)
            result = {}
            for ht_type, grp in df.groupby("Helix_Type"):
                result[ht_type] = grp["Length"].tolist()
            return result
        return None

    for pid in selected_ids:
        print(f"\n-> Plotando: {pid}")
        try:
            # ── Grupo 1: Analise basica ───────────────────────────────────────
            if pid == "helix_propensities":
                df = _get_propensity_df()
                if df is not None:
                    P.plot_helix_propensities(df, analysis_path)
                else:
                    print("  [AVISO] Dados de propensao nao disponiveis. Execute a etapa 4.")

            elif pid == "helix_positions":
                df = _get_positions_df()
                if df is not None:
                    P.plot_helix_positions(df, analysis_path)
                else:
                    print("  [AVISO] Dados de posicoes nao disponiveis. Execute a etapa 4.")

            elif pid == "helix_lengths":
                lengths = _get_helix_lengths()
                if lengths:
                    P.plot_helix_length_distribution(lengths, analysis_path)
                else:
                    print("  [AVISO] Dados de comprimentos nao disponiveis. Execute a etapa 4.")

            elif pid == "heptad_pattern":
                heptad = _get_heptad()
                if heptad:
                    P.plot_hydrophobic_heptad(heptad, analysis_path)
                else:
                    print("  [AVISO] Dados de heptad nao disponiveis. Execute a etapa 4.")

            # ── Grupo 2: Tipos de helices ─────────────────────────────────────
            elif pid == "helix_type_distribution":
                counts = _get_helix_counts()
                if counts:
                    P.plot_helix_type_distribution(counts, analysis_path)
                else:
                    print("  [AVISO] Dados de tipos nao disponiveis. Execute a etapa 4.")

            elif pid == "helix_type_composition":
                df = _get_composition_df()
                if df is not None:
                    P.plot_helix_type_composition(df, analysis_path)
                else:
                    print("  [AVISO] Dados de composicao por tipo nao disponiveis. Execute a etapa 4.")

            elif pid == "helix_type_top_amino_acids":
                df = _get_top_residues_df()
                if df is not None:
                    P.plot_top_amino_acids_by_type(df, analysis_path)
                else:
                    print("  [AVISO] Dados de comparacao nao disponiveis. Execute a etapa 4.")

            elif pid == "helix_length_by_type":
                lengths_by_type = _get_lengths_by_type()
                if lengths_by_type:
                    P.plot_helix_length_comparison(lengths_by_type, analysis_path)
                else:
                    print("  [AVISO] Dados de comprimento por tipo nao disponiveis. Execute a etapa 4.")

            elif pid == "helix_type_statistical_comparison":
                df = _get_stat_df()
                if df is not None:
                    P.plot_statistical_differences(df, analysis_path)
                else:
                    print("  [AVISO] Dados estatisticos nao disponiveis. Execute a etapa 4.")

            # ── Grupo 3: Restricoes estruturais ──────────────────────────────
            elif pid == "helical_wheel_average":
                seqs = evo_data.get("helix_sequences") if evo_data else None
                if seqs:
                    P.plot_helical_wheel_average(seqs, analysis_path)
                else:
                    print("  [AVISO] Dados evolutivos nao disponiveis. Execute a etapa 4.")

            elif pid == "hydrophobic_moment_distribution":
                seqs = evo_data.get("helix_sequences") if evo_data else None
                if seqs:
                    moments = E.compute_hydrophobic_moments(seqs)
                    if moments:
                        P.plot_hydrophobic_moment_distribution(moments, analysis_path)
                    else:
                        print("  [AVISO] Nenhum momento calculado.")
                else:
                    print("  [AVISO] Dados evolutivos nao disponiveis. Execute a etapa 4.")

            elif pid == "ncap_ccap_preferences":
                ncap_pos = evo_data.get("ncap_position") if evo_data else None
                ccap_pos = evo_data.get("ccap_position") if evo_data else None
                if ncap_pos and ccap_pos:
                    P.plot_ncap_ccap(ncap_pos, ccap_pos, analysis_path)
                else:
                    print("  [AVISO] Dados de N-cap/C-cap nao disponiveis. Execute a etapa 4.")

            # ── Grupo 4: Composicao e variabilidade ───────────────────────────
            elif pid == "pca_aa_composition":
                ps = per_structure
                if not ps and evo_data:
                    ps = evo_data.get("per_structure") or {}
                if not ps and evo_data:
                    from collections import Counter as _Counter
                    per_helix = evo_data.get("per_helix_data", [])
                    ps = {str(i): dict(item["aa_counts"]) for i, item in enumerate(per_helix)}
                if ps:
                    P.plot_pca_composition(ps, analysis_path)
                else:
                    print("  [AVISO] Dados por estrutura nao disponiveis. Execute a etapa 3.")

            elif pid == "helix_content_distribution":
                contents = evo_data.get("helix_content_per_structure") if evo_data else None
                if contents:
                    P.plot_helix_content_distribution(contents, analysis_path)
                else:
                    print("  [AVISO] Dados de conteudo helicoidal nao disponiveis. Execute a etapa 4.")

            elif pid == "aa_cooccurrence":
                seqs = evo_data.get("helix_sequences") if evo_data else None
                if seqs:
                    cooc_df = E.compute_aa_cooccurrence(seqs)
                    P.plot_aa_cooccurrence(cooc_df, analysis_path)
                else:
                    print("  [AVISO] Dados evolutivos nao disponiveis. Execute a etapa 4.")

            elif pid == "helix_length_vs_composition":
                per_helix = evo_data.get("per_helix_data") if evo_data else None
                if per_helix:
                    length_df = E.compute_helix_length_composition(per_helix)
                    P.plot_helix_length_vs_composition(length_df, analysis_path)
                else:
                    csv_path = analysis_path / "helix_length_vs_composition.csv"
                    if csv_path.exists():
                        length_df = pd.read_csv(csv_path)
                        P.plot_helix_length_vs_composition(length_df, analysis_path)
                    else:
                        print("  [AVISO] Dados por helice nao disponiveis. Execute a etapa 4.")

            # ── Grupo 5: Transicoes e historia ───────────────────────────────
            elif pid == "helix_transition_matrix":
                trans = evo_data.get("transitions") if evo_data else None
                if trans:
                    P.plot_helix_transition_matrix(trans, analysis_path)
                else:
                    print("  [AVISO] Dados de transicao nao disponiveis. Execute a etapa 4.")

            elif pid == "g_ratio_by_length":
                lbt = _get_lengths_by_type()
                if lbt:
                    P.plot_g_ratio_by_length(lbt, analysis_path)
                else:
                    print("  [AVISO] Dados de comprimento por tipo nao disponiveis. Execute a etapa 4.")

            elif pid == "shannon_entropy_heptad":
                heptad_dist = evo_data.get("heptad_aa_distribution") if evo_data else None
                if heptad_dist:
                    entropy_df = E.compute_shannon_entropy_heptad(heptad_dist)
                    P.plot_shannon_entropy_heptad(entropy_df, analysis_path)
                else:
                    csv_path = analysis_path / "shannon_entropy_heptad.csv"
                    if csv_path.exists():
                        entropy_df = pd.read_csv(csv_path)
                        P.plot_shannon_entropy_heptad(entropy_df, analysis_path)
                    else:
                        print("  [AVISO] Dados de heptad nao disponiveis. Execute a etapa 4.")

            # ── Grupo 6: Vies do codigo genetico ─────────────────────────────
            elif pid == "codon_degeneracy_vs_propensity":
                csv_path = analysis_path / "codon_degeneracy.csv"
                if csv_path.exists():
                    codon_df = pd.read_csv(csv_path)
                    P.plot_codon_degeneracy_vs_propensity(codon_df, analysis_path)
                else:
                    print("  [AVISO] codon_degeneracy.csv nao encontrado. Execute a etapa 4.")

            elif pid == "proteome_comparison":
                csv_path = analysis_path / "proteome_comparison.csv"
                if csv_path.exists():
                    prot_df = pd.read_csv(csv_path)
                    P.plot_proteome_comparison(prot_df, analysis_path)
                else:
                    print("  [AVISO] proteome_comparison.csv nao encontrado. Execute a etapa 4.")

            else:
                print(f"  [AVISO] plot_id desconhecido: {pid}")

        except Exception as e:
            print(f"  [ERRO] {e}")
            logger.debug("Traceback completo:", exc_info=True)


# ── Mapa de etapas ────────────────────────────────────────────────────────────

STEPS_WITH_CACHE = {
    3: step3_frequency,
    4: step4_helix_analysis,
    5: step5_helix_types,
    6: step6_evolutionary,
}

ALL_STEPS = {1: step1_download, 2: step2_clean, **STEPS_WITH_CACHE}

_EXPLORE_STEPS = [3, 4, 5, 6]  # pula download e limpeza


def _resolve_plot_preselection(tokens: list[str]) -> set[str]:
    """Converte tokens --plots (ex: 'g3', 'helical_wheel_average') em set de plot_ids."""
    from cath_analysis.menu import PLOT_CATALOG
    groups: dict[str, list[str]] = {}
    for plot_id, info in PLOT_CATALOG.items():
        groups.setdefault(info["grupo"], []).append(plot_id)
    group_list = list(groups.keys())

    selected: set[str] = set()
    for token in tokens:
        t = token.lower()
        if t.startswith("g") and t[1:].isdigit():
            g_idx = int(t[1:]) - 1
            if 0 <= g_idx < len(group_list):
                selected.update(groups[group_list[g_idx]])
            else:
                print(f"  [AVISO] Grupo invalido: {token} (use g1-g{len(group_list)})")
        elif token in PLOT_CATALOG:
            selected.add(token)
        else:
            print(f"  [AVISO] Plot desconhecido: '{token}' (use --list-plots para ver opcoes)")
    return selected


def _list_plots() -> None:
    """Imprime o catalogo completo de graficos disponiveis."""
    from cath_analysis.menu import PLOT_CATALOG
    current_group = None
    print()
    for i, (plot_id, info) in enumerate(PLOT_CATALOG.items(), start=1):
        if info["grupo"] != current_group:
            current_group = info["grupo"]
            print(f"\n  {current_group}")
            print("  " + "-" * 56)
        tags = []
        if info["requer_dssp"]:
            tags.append("DSSP")
        if info["requer_evo"]:
            tags.append("EVO")
        tag_str = f"[{'+'.join(tags)}]" if tags else ""
        print(f"  {i:>2}. {tag_str:<10} {info['nome']}")
        print(f"       ID: {plot_id}")
    print()


# ── Ponto de entrada ──────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Pipeline Evo-Hex",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--steps",
        nargs="+",
        type=int,
        choices=list(ALL_STEPS),
        default=None,
        metavar="N",
        help="Etapas a executar (1-6). Padrao: todas.",
    )
    parser.add_argument(
        "--explore",
        action="store_true",
        default=False,
        help="Atalho para --steps 3 4 5 6: pula download e limpeza, vai direto a analise.",
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        default=False,
        help="Re-baixa e reprocessa mesmo que os dados ja existam (sem perguntar).",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        default=False,
        help="Pula etapas 1 e 2 sem perguntar, mesmo que os dados existam.",
    )
    parser.add_argument(
        "--no-interactive",
        action="store_true",
        default=False,
        help="Desativa todos os menus interativos (download e graficos).",
    )
    parser.add_argument(
        "--plots",
        nargs="+",
        default=None,
        metavar="GRUPO_OU_ID",
        help=(
            "Pre-seleciona graficos no menu. Use 'g1'..'g6' para grupos "
            "ou IDs individuais. Ex: --plots g3 g4 helical_wheel_average"
        ),
    )
    parser.add_argument(
        "--list-plots",
        action="store_true",
        default=False,
        help="Lista todos os graficos disponiveis e sai.",
    )
    parser.add_argument(
        "--full-dataset",
        action="store_true",
        default=False,
        help=(
            "Baixa o dataset completo (~50k estruturas) sem filtro S40. "
            "Por padrao o pipeline usa apenas o subconjunto nao-redundante S40."
        ),
    )
    args = parser.parse_args()

    if args.list_plots:
        _list_plots()
        return

    interactive = not args.no_interactive

    # Resolve quais etapas executar
    if args.steps:
        steps_to_run = sorted(args.steps)
    elif args.explore or args.skip_download:
        steps_to_run = _EXPLORE_STEPS
    else:
        steps_to_run = sorted(ALL_STEPS.keys())

    print("=" * 60)
    print("PIPELINE EVO-HEX")
    print(f"  Base:      {BASE_PATH}")
    print(f"  Etapas:    {steps_to_run}")
    print(f"  Interativo: {'Sim' if interactive else 'Nao'}")
    if args.force_download:
        print("  Modo:      forca re-download")
    elif args.skip_download or args.explore:
        print("  Modo:      exploracao (sem download)")
    if args.full_dataset:
        print("  Dataset:   completo (~50k estruturas)")
    else:
        print("  Dataset:   S40 nao-redundante")
    print("=" * 60)

    cache: dict = {}

    for step_num in steps_to_run:
        try:
            if step_num == 1:
                step1_download(
                    force=args.force_download,
                    interactive=interactive,
                    full_dataset=args.full_dataset,
                )
            elif step_num == 2:
                step2_clean(force=args.force_download, interactive=interactive)
            elif step_num in (3, 4, 5):
                STEPS_WITH_CACHE[step_num](cache, interactive=interactive)
            elif step_num == 6:
                STEPS_WITH_CACHE[step_num](cache)
            else:
                STEPS_WITH_CACHE[step_num](cache)
        except Exception:  # noqa: BLE001
            logger.error("Etapa %d falhou:\n%s", step_num, traceback.format_exc())
            sys.exit(1)

    print("\n" + "=" * 60)
    print("PIPELINE CONCLUIDO")
    print("=" * 60)

    # ── Menu interativo de graficos ───────────────────────────────────────────
    has_data = any(
        k in cache
        for k in ("global_counter", "evo_data", "propensity_df", "composition_df")
    )

    if not has_data:
        print("\n  Nenhum dado de analise em cache. Execute etapas 3-6 para plotar.")
    elif not interactive:
        from cath_analysis.menu import PLOT_CATALOG
        all_ids = list(PLOT_CATALOG.keys())
        print(f"\nModo nao interativo: plotando {len(all_ids)} graficos...")
        run_plots(all_ids, cache, ANALYSIS_PATH)
        print(f"Graficos salvos em: {ANALYSIS_PATH}")
    else:
        from cath_analysis.menu import select_plots

        # Pre-selecao via --plots
        preselected: set[str] = set()
        if args.plots:
            preselected = _resolve_plot_preselection(args.plots)
            if preselected:
                print(f"\n  {len(preselected)} grafico(s) pre-selecionados via --plots.")

        print("\nAbrindo seletor de graficos...")
        selected_ids = select_plots(preselected=preselected)

        if selected_ids:
            print(f"\nPlotando {len(selected_ids)} grafico(s)...")
            run_plots(selected_ids, cache, ANALYSIS_PATH)
            print(f"\nGraficos salvos em: {ANALYSIS_PATH}")
        else:
            print("\nNenhum grafico selecionado.")

    print("\n" + "=" * 60)
    print("ENCERRADO")
    print("=" * 60)


if __name__ == "__main__":
    main()
