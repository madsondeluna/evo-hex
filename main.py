#!/usr/bin/env python3
"""
Pipeline Evo-Hex – ponto de entrada principal.

Etapas
------
1. Download da lista de domínios CATH e estruturas PDB
2. Limpeza (cadeia A, sem heteroátomos)
3. Análise de frequência de aminoácidos
4. Análise avançada de hélices com DSSP (propensão, posições, comprimentos)
5. Classificação e comparação de tipos de hélices (H/G/I)
6. Coleta de dados evolutivos (hélice wheel, momentos hidrofóbicos, etc.)

Uso
---
    python main.py                      # todas as etapas, modo interativo
    python main.py --explore            # pula download/limpeza, vai direto à análise
    python main.py --steps 3 4 5 6     # mesmo que --explore (manual)
    python main.py --steps 1 2         # apenas download e limpeza
    python main.py --force-download    # baixa mesmo se dados existirem
    python main.py --skip-download     # pula etapas 1-2 sem perguntar
    python main.py --no-interactive    # sem menus, executa tudo automaticamente
    python main.py --plots g3 g4       # pré-seleciona grupos no menu de gráficos
    python main.py --list-plots        # lista todos os gráficos disponíveis
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
    """Imprime resumo dos dados já existentes no disco."""
    print("\n  ┌─ Dados existentes ───────────────────────────────────┐")
    if info["has_cath_list"]:
        print(f"  │  Índice CATH:       atualizado em {info['cath_list_date']}")
    else:
        print("  │  Índice CATH:       não encontrado")
    if info["raw_count"]:
        print(f"  │  Estruturas raw:    {info['raw_count']:,} PDBs  (último: {info['raw_date']})")
    else:
        print("  │  Estruturas raw:    nenhuma")
    if info["clean_count"]:
        print(f"  │  Estruturas clean:  {info['clean_count']:,} PDBs  (último: {info['clean_date']})")
    else:
        print("  │  Estruturas clean:  nenhuma")
    print("  └──────────────────────────────────────────────────────┘")


def _ask_redownload() -> bool:
    """Pergunta ao usuário se deseja baixar novamente. Retorna True se sim."""
    while True:
        resp = input("\n  Baixar novamente? [s/N] ").strip().lower()
        if resp in ("s", "sim", "y", "yes"):
            return True
        if resp in ("", "n", "não", "nao", "no"):
            return False
        print("  Responda 's' para sim ou 'n' para não.")


def step1_download(force: bool = False, interactive: bool = True) -> None:
    """Etapa 1: Download de estruturas PDB mainly-alpha.

    Args:
        force: Se True, baixa sem perguntar mesmo que os dados existam.
        interactive: Se True, pergunta ao usuário quando dados já existem.
    """
    from cath_analysis.downloader import (
        check_existing_data,
        download_cath_domain_list,
        download_structures_parallel,
        parse_cath_domains,
        print_download_report,
        save_download_results,
        setup_directories,
    )

    print("\n" + "=" * 60)
    print("ETAPA 1 – DOWNLOAD DE ESTRUTURAS")
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
    pdb_codes = parse_cath_domains(cath_file)
    results = download_structures_parallel(pdb_codes, structures_path)
    save_download_results(base_path, pdb_codes, results)
    print_download_report(results, structures_path)


def step2_clean(force: bool = False, interactive: bool = True) -> None:
    """Etapa 2: Limpeza das estruturas PDB.

    Args:
        force: Se True, reprocessa mesmo que já existam arquivos limpos.
        interactive: Se True, pergunta ao usuário quando dados já existem.
    """
    from cath_analysis.downloader import check_existing_data
    from cath_analysis.cleaner import (
        clean_all_structures,
        print_cleaning_summary,
        save_cleaning_results,
    )

    print("\n" + "=" * 60)
    print("ETAPA 2 – LIMPEZA DE ESTRUTURAS")
    print("=" * 60)

    info = check_existing_data()
    has_clean = info["clean_count"] > 0

    if has_clean and not force:
        _show_existing_data_info(info)
        if interactive:
            print(f"\n  Já existem {info['clean_count']:,} estruturas limpas.")
            if not _ask_redownload():
                print("  Limpeza ignorada. Usando estruturas existentes.")
                return
        else:
            print(f"\n  [INFO] {info['clean_count']:,} estruturas limpas já existem. Use --force-download para reprocessar.")
            return

    print()
    STRUCTURES_CLEAN_PATH.mkdir(parents=True, exist_ok=True)
    clean_results = clean_all_structures(STRUCTURES_PATH, STRUCTURES_CLEAN_PATH)
    if clean_results is not None:
        save_cleaning_results(LOGS_PATH, clean_results)
        print_cleaning_summary(clean_results, STRUCTURES_PATH, STRUCTURES_CLEAN_PATH)


def step3_frequency(cache: dict) -> None:
    """Etapa 3: Análise de frequência de aminoácidos."""
    from cath_analysis.frequency_analysis import (
        analyze_amino_acid_frequency,
        print_frequency_summary,
        save_frequency_results,
    )

    print("\n" + "=" * 60)
    print("ETAPA 3 – FREQUÊNCIA DE AMINOÁCIDOS")
    print("=" * 60 + "\n")

    global_counter, per_structure = analyze_amino_acid_frequency(STRUCTURES_CLEAN_PATH)
    ANALYSIS_PATH.mkdir(parents=True, exist_ok=True)
    save_frequency_results(ANALYSIS_PATH, global_counter, per_structure)
    print_frequency_summary(global_counter, per_structure)

    cache["global_counter"] = global_counter
    cache["per_structure"] = per_structure


def step4_helix_analysis(cache: dict) -> None:
    """Etapa 4: Análise avançada de hélices com DSSP."""
    from cath_analysis.helix_analysis import run_helix_analysis

    print("\n" + "=" * 60)
    print("ETAPA 4 – ANÁLISE AVANÇADA DE HÉLICES")
    print("=" * 60 + "\n")

    dssp_results = run_helix_analysis()
    cache["dssp_results"] = dssp_results


def step5_helix_types(cache: dict) -> None:
    """Etapa 5: Classificação de tipos de hélices (H/G/I)."""
    from cath_analysis.helix_types import run_helix_type_analysis

    print("\n" + "=" * 60)
    print("ETAPA 5 – TIPOS DE HÉLICES")
    print("=" * 60 + "\n")

    helix_type_results = run_helix_type_analysis()
    cache["helix_type_results"] = helix_type_results


def step6_evolutionary(cache: dict) -> None:
    """Etapa 6: Coleta de dados evolutivos em um único passo DSSP."""
    from cath_analysis.evolutionary_analysis import collect_evolutionary_data

    print("\n" + "=" * 60)
    print("ETAPA 6 – ANÁLISE EVOLUTIVA DE HÉLICES")
    print("=" * 60 + "\n")

    evo_data = collect_evolutionary_data(STRUCTURES_CLEAN_PATH)
    cache["evo_data"] = evo_data

    n_helices = len(evo_data.get("helix_sequences", []))
    n_structs = len(evo_data.get("helix_content_per_structure", []))
    print(f"  Estruturas processadas: {n_structs:,}")
    print(f"  Hélices alpha (≥4 res): {n_helices:,}")
    print(f"  Pares de transição:     {sum(evo_data['transitions'].values()):,}")


# ── Execução de gráficos ──────────────────────────────────────────────────────

def run_plots(selected_ids: list, cache: dict, analysis_path: Path) -> None:
    """Despacha cada plot_id para a função de plotagem correspondente.

    Args:
        selected_ids: Lista de plot_ids (strings) a plotar.
        cache: Dicionário com dados acumulados pelas etapas.
        analysis_path: Diretório de saída para os gráficos.
    """
    from cath_analysis import plotting as P
    from cath_analysis import evolutionary_analysis as E

    analysis_path.mkdir(parents=True, exist_ok=True)

    evo_data: dict = cache.get("evo_data", {})
    dssp_results: dict = cache.get("dssp_results") or {}
    helix_type_results: dict = cache.get("helix_type_results") or {}
    global_counter = cache.get("global_counter")
    per_structure: dict = cache.get("per_structure") or {}

    for pid in selected_ids:
        print(f"\n-> Plotando: {pid}")
        try:
            # ── Grupo 1: Análise básica ──────────────────────────────────────
            if pid == "helix_propensities":
                df = dssp_results.get("propensity_df") if dssp_results else None
                if df is not None:
                    P.plot_helix_propensities(df, analysis_path)
                else:
                    print("  [AVISO] Dados de propensão não disponíveis. Execute a etapa 4.")

            elif pid == "helix_positions":
                df = dssp_results.get("positions_df") if dssp_results else None
                if df is not None:
                    P.plot_helix_positions(df, analysis_path)
                else:
                    print("  [AVISO] Dados de posições não disponíveis. Execute a etapa 4.")

            elif pid == "helix_lengths":
                lengths = dssp_results.get("helix_lengths") if dssp_results else None
                if lengths:
                    P.plot_helix_length_distribution(lengths, analysis_path)
                else:
                    print("  [AVISO] Dados de comprimentos não disponíveis. Execute a etapa 4.")

            elif pid == "heptad_pattern":
                heptad = dssp_results.get("heptad_data") if dssp_results else None
                if heptad is None and evo_data:
                    # Converte heptad_aa_distribution para o formato esperado
                    from cath_analysis.config import HYDROPHOBIC_AA
                    heptad_dist = evo_data.get("heptad_aa_distribution", {})
                    heptad = {}
                    for pos in range(7):
                        cnts = heptad_dist.get(pos, {})
                        total = sum(cnts.values())
                        hydro = sum(v for aa, v in cnts.items() if aa in HYDROPHOBIC_AA)
                        heptad[pos] = {"hydrophobic": hydro, "total": total}
                if heptad:
                    P.plot_hydrophobic_heptad(heptad, analysis_path)
                else:
                    print("  [AVISO] Dados de heptad não disponíveis. Execute a etapa 4 ou 6.")

            # ── Grupo 2: Tipos de hélices ────────────────────────────────────
            elif pid == "helix_type_distribution":
                counts = helix_type_results.get("helix_counts") if helix_type_results else None
                if counts:
                    P.plot_helix_type_distribution(counts, analysis_path)
                else:
                    print("  [AVISO] Dados de tipos não disponíveis. Execute a etapa 5.")

            elif pid == "helix_type_composition":
                df = helix_type_results.get("composition_df") if helix_type_results else None
                if df is not None:
                    P.plot_helix_type_composition(df, analysis_path)
                else:
                    print("  [AVISO] Dados de composição por tipo não disponíveis. Execute a etapa 5.")

            elif pid == "helix_type_top_amino_acids":
                df = helix_type_results.get("comparison_df") if helix_type_results else None
                if df is not None:
                    P.plot_top_amino_acids_by_type(df, analysis_path)
                else:
                    print("  [AVISO] Dados de comparação não disponíveis. Execute a etapa 5.")

            elif pid == "helix_length_by_type":
                lengths_by_type = helix_type_results.get("helix_lengths") if helix_type_results else None
                if lengths_by_type is None and evo_data:
                    lengths_by_type = evo_data.get("helix_lengths_by_type")
                if lengths_by_type:
                    P.plot_helix_length_comparison(lengths_by_type, analysis_path)
                else:
                    print("  [AVISO] Dados de comprimento por tipo não disponíveis. Execute a etapa 5 ou 6.")

            elif pid == "helix_type_statistical_comparison":
                df = helix_type_results.get("stat_df") if helix_type_results else None
                if df is not None:
                    P.plot_statistical_differences(df, analysis_path)
                else:
                    print("  [AVISO] Dados estatísticos não disponíveis. Execute a etapa 5.")

            # ── Grupo 3: Restrições estruturais ─────────────────────────────
            elif pid == "helical_wheel_average":
                seqs = evo_data.get("helix_sequences") if evo_data else None
                if seqs:
                    P.plot_helical_wheel_average(seqs, analysis_path)
                else:
                    print("  [AVISO] Dados evolutivos não disponíveis. Execute a etapa 6.")

            elif pid == "hydrophobic_moment_distribution":
                seqs = evo_data.get("helix_sequences") if evo_data else None
                if seqs:
                    moments = E.compute_hydrophobic_moments(seqs)
                    if moments:
                        P.plot_hydrophobic_moment_distribution(moments, analysis_path)
                    else:
                        print("  [AVISO] Nenhum momento calculado.")
                else:
                    print("  [AVISO] Dados evolutivos não disponíveis. Execute a etapa 6.")

            elif pid == "ncap_ccap_preferences":
                ncap_pos = evo_data.get("ncap_position") if evo_data else None
                ccap_pos = evo_data.get("ccap_position") if evo_data else None
                if ncap_pos and ccap_pos:
                    P.plot_ncap_ccap(ncap_pos, ccap_pos, analysis_path)
                else:
                    print("  [AVISO] Dados de N-cap/C-cap não disponíveis. Execute a etapa 6.")

            # ── Grupo 4: Composição e variabilidade ──────────────────────────
            elif pid == "pca_aa_composition":
                ps = per_structure or (evo_data.get("per_structure") if evo_data else None)
                if ps:
                    P.plot_pca_composition(ps, analysis_path)
                else:
                    print("  [AVISO] Dados por estrutura não disponíveis. Execute a etapa 3.")

            elif pid == "helix_content_distribution":
                contents = evo_data.get("helix_content_per_structure") if evo_data else None
                if contents:
                    P.plot_helix_content_distribution(contents, analysis_path)
                else:
                    print("  [AVISO] Dados de conteúdo helicoidal não disponíveis. Execute a etapa 6.")

            elif pid == "aa_cooccurrence":
                seqs = evo_data.get("helix_sequences") if evo_data else None
                if seqs:
                    import pandas as pd
                    cooc_df = E.compute_aa_cooccurrence(seqs)
                    P.plot_aa_cooccurrence(cooc_df, analysis_path)
                else:
                    print("  [AVISO] Dados evolutivos não disponíveis. Execute a etapa 6.")

            elif pid == "helix_length_vs_composition":
                per_helix = evo_data.get("per_helix_data") if evo_data else None
                if per_helix:
                    length_df = E.compute_helix_length_composition(per_helix)
                    P.plot_helix_length_vs_composition(length_df, analysis_path)
                else:
                    print("  [AVISO] Dados por hélice não disponíveis. Execute a etapa 6.")

            # ── Grupo 5: Transições e história ──────────────────────────────
            elif pid == "helix_transition_matrix":
                trans = evo_data.get("transitions") if evo_data else None
                if trans is None:
                    trans = helix_type_results.get("transitions") if helix_type_results else None
                if trans:
                    P.plot_helix_transition_matrix(trans, analysis_path)
                else:
                    print("  [AVISO] Dados de transição não disponíveis. Execute a etapa 5 ou 6.")

            elif pid == "g_ratio_by_length":
                lbt = evo_data.get("helix_lengths_by_type") if evo_data else None
                if lbt is None:
                    lbt = helix_type_results.get("helix_lengths") if helix_type_results else None
                if lbt:
                    P.plot_g_ratio_by_length(lbt, analysis_path)
                else:
                    print("  [AVISO] Dados de comprimento por tipo não disponíveis. Execute a etapa 5 ou 6.")

            elif pid == "shannon_entropy_heptad":
                heptad_dist = evo_data.get("heptad_aa_distribution") if evo_data else None
                if heptad_dist:
                    entropy_df = E.compute_shannon_entropy_heptad(heptad_dist)
                    P.plot_shannon_entropy_heptad(entropy_df, analysis_path)
                else:
                    print("  [AVISO] Dados de heptad não disponíveis. Execute a etapa 6.")

            # ── Grupo 6: Viés do código genético ────────────────────────────
            elif pid == "codon_degeneracy_vs_propensity":
                if global_counter and evo_data:
                    from collections import Counter
                    # Reconstrói helix_residues a partir das sequências evolutivas
                    seqs = evo_data.get("helix_sequences", [])
                    # Mapeia 1-letra → 3-letras
                    _one_to_three = {
                        "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
                        "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
                        "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
                        "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
                    }
                    helix_residues: Counter = Counter()
                    for seq in seqs:
                        for aa1 in seq:
                            aa3 = _one_to_three.get(aa1)
                            if aa3:
                                helix_residues[aa3] += 1
                    codon_df = E.compute_codon_degeneracy_vs_propensity(helix_residues, global_counter)
                    P.plot_codon_degeneracy_vs_propensity(codon_df, analysis_path)
                else:
                    print("  [AVISO] Dados não disponíveis. Execute as etapas 3 e 6.")

            elif pid == "proteome_comparison":
                gc = global_counter
                if gc:
                    prot_df = E.compute_proteome_comparison(gc)
                    P.plot_proteome_comparison(prot_df, analysis_path)
                else:
                    print("  [AVISO] Dados globais não disponíveis. Execute a etapa 3.")

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
                print(f"  [AVISO] Grupo inválido: {token} (use g1–g{len(group_list)})")
        elif token in PLOT_CATALOG:
            selected.add(token)
        else:
            print(f"  [AVISO] Plot desconhecido: '{token}' (use --list-plots para ver opções)")
    return selected


def _list_plots() -> None:
    """Imprime o catálogo completo de gráficos disponíveis."""
    from cath_analysis.menu import PLOT_CATALOG
    current_group = None
    print()
    for i, (plot_id, info) in enumerate(PLOT_CATALOG.items(), start=1):
        if info["grupo"] != current_group:
            current_group = info["grupo"]
            print(f"\n  {current_group}")
            print("  " + "─" * 56)
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
        help="Etapas a executar (1-6). Padrão: todas.",
    )
    parser.add_argument(
        "--explore",
        action="store_true",
        default=False,
        help="Atalho para --steps 3 4 5 6: pula download e limpeza, vai direto à análise.",
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        default=False,
        help="Re-baixa e reprocessa mesmo que os dados já existam (sem perguntar).",
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
        help="Desativa todos os menus interativos (download e gráficos).",
    )
    parser.add_argument(
        "--plots",
        nargs="+",
        default=None,
        metavar="GRUPO_OU_ID",
        help=(
            "Pré-seleciona gráficos no menu. Use 'g1'..'g6' para grupos "
            "ou IDs individuais. Ex: --plots g3 g4 helical_wheel_average"
        ),
    )
    parser.add_argument(
        "--list-plots",
        action="store_true",
        default=False,
        help="Lista todos os gráficos disponíveis e sai.",
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
    print(f"  Interativo: {'Sim' if interactive else 'Não'}")
    if args.force_download:
        print("  Modo:      força re-download")
    elif args.skip_download or args.explore:
        print("  Modo:      exploração (sem download)")
    print("=" * 60)

    cache: dict = {}

    for step_num in steps_to_run:
        try:
            if step_num == 1:
                step1_download(force=args.force_download, interactive=interactive)
            elif step_num == 2:
                step2_clean(force=args.force_download, interactive=interactive)
            else:
                STEPS_WITH_CACHE[step_num](cache)
        except Exception:  # noqa: BLE001
            logger.error("Etapa %d falhou:\n%s", step_num, traceback.format_exc())
            sys.exit(1)

    print("\n" + "=" * 60)
    print("PIPELINE CONCLUÍDO")
    print("=" * 60)

    # ── Menu interativo de gráficos ──────────────────────────────────────────
    has_data = any(k in cache for k in ("global_counter", "dssp_results",
                                         "helix_type_results", "evo_data"))

    if not has_data:
        print("\n  Nenhum dado de análise em cache. Execute etapas 3-6 para plotar.")
    elif not interactive:
        from cath_analysis.menu import PLOT_CATALOG
        all_ids = list(PLOT_CATALOG.keys())
        print(f"\nModo não interativo: plotando {len(all_ids)} gráficos...")
        run_plots(all_ids, cache, ANALYSIS_PATH)
        print(f"Gráficos salvos em: {ANALYSIS_PATH}")
    else:
        from cath_analysis.menu import select_plots

        # Pré-seleção via --plots
        preselected: set[str] = set()
        if args.plots:
            preselected = _resolve_plot_preselection(args.plots)
            if preselected:
                print(f"\n  {len(preselected)} gráfico(s) pré-selecionados via --plots.")

        print("\nAbrindo seletor de gráficos...")
        selected_ids = select_plots(preselected=preselected)

        if selected_ids:
            print(f"\nPlotando {len(selected_ids)} gráfico(s)...")
            run_plots(selected_ids, cache, ANALYSIS_PATH)
            print(f"\nGráficos salvos em: {ANALYSIS_PATH}")
        else:
            print("\nNenhum gráfico selecionado.")

    print("\n" + "=" * 60)
    print("ENCERRADO")
    print("=" * 60)


if __name__ == "__main__":
    main()
