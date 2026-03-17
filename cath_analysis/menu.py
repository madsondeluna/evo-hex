"""
Menu interativo para seleção de gráficos a plotar.
"""

from pathlib import Path
from typing import Callable

# Definição completa de todos os gráficos disponíveis
# Estrutura: id -> {nome, grupo, descricao, requer_dssp, requer_evo}
PLOT_CATALOG = {
    # ── Grupo 1: Análise básica de hélices ───────────────────────────────────
    "helix_propensities": {
        "nome": "Propensão para hélices (observada vs Chou-Fasman)",
        "grupo": "1. Análise Básica de Hélices",
        "descricao": "Bar chart duplo + scatter: propensão de cada AA vs escala teórica",
        "requer_dssp": True,
        "requer_evo": False,
    },
    "helix_positions": {
        "nome": "Preferência por posição na hélice",
        "grupo": "1. Análise Básica de Hélices",
        "descricao": "Heatmap: frequência de AAs no N-terminal, meio e C-terminal",
        "requer_dssp": True,
        "requer_evo": False,
    },
    "helix_lengths": {
        "nome": "Distribuição de comprimentos de hélices",
        "grupo": "1. Análise Básica de Hélices",
        "descricao": "Histograma + boxplot dos comprimentos (em resíduos)",
        "requer_dssp": True,
        "requer_evo": False,
    },
    "heptad_pattern": {
        "nome": "Padrão heptad de hidrofobicidade",
        "grupo": "1. Análise Básica de Hélices",
        "descricao": "Barras: frequência hidrofóbica em cada posição a-g do heptad",
        "requer_dssp": False,
        "requer_evo": False,
    },
    # ── Grupo 2: Tipos de hélices ─────────────────────────────────────────────
    "helix_type_distribution": {
        "nome": "Distribuição de tipos (H/G/I)",
        "grupo": "2. Tipos de Hélices",
        "descricao": "Bar chart + pie: proporção de Alpha, 3-10 e Pi helices",
        "requer_dssp": True,
        "requer_evo": False,
    },
    "helix_type_composition": {
        "nome": "Composição por tipo de hélice (heatmap)",
        "grupo": "2. Tipos de Hélices",
        "descricao": "Heatmap 20×3: frequência de cada AA em H, G e I",
        "requer_dssp": True,
        "requer_evo": False,
    },
    "helix_type_top_amino_acids": {
        "nome": "Top 10 AAs por tipo de hélice",
        "grupo": "2. Tipos de Hélices",
        "descricao": "Barras horizontais: AAs mais frequentes em cada tipo",
        "requer_dssp": True,
        "requer_evo": False,
    },
    "helix_length_by_type": {
        "nome": "Comprimento por tipo de hélice",
        "grupo": "2. Tipos de Hélices",
        "descricao": "Boxplot + histograma comparando comprimentos H/G/I",
        "requer_dssp": True,
        "requer_evo": False,
    },
    "helix_type_statistical_comparison": {
        "nome": "Comparação estatística Alpha vs 3-10",
        "grupo": "2. Tipos de Hélices",
        "descricao": "Diferença de frequência e scatter direto entre H e G",
        "requer_dssp": True,
        "requer_evo": False,
    },
    # ── Grupo 3: Restrições estruturais ──────────────────────────────────────
    "helical_wheel_average": {
        "nome": "Helical wheel composto",
        "grupo": "3. Restrições Estruturais",
        "descricao": "Projeção polar: frequência hidrofóbica por ângulo na hélice (100°/resíduo)",
        "requer_dssp": True,
        "requer_evo": True,
    },
    "hydrophobic_moment_distribution": {
        "nome": "Momento hidrofóbico de Eisenberg",
        "grupo": "3. Restrições Estruturais",
        "descricao": "Distribuição do μH por hélice + classificação por anfipatiticidade",
        "requer_dssp": True,
        "requer_evo": True,
    },
    "ncap_ccap_preferences": {
        "nome": "N-cap e C-cap (3 posições cada)",
        "grupo": "3. Restrições Estruturais",
        "descricao": "Heatmap das preferências nas 3 primeiras/últimas posições da hélice",
        "requer_dssp": True,
        "requer_evo": True,
    },
    # ── Grupo 4: Composição e variabilidade ───────────────────────────────────
    "pca_aa_composition": {
        "nome": "PCA da composição por estrutura",
        "grupo": "4. Composição e Variabilidade",
        "descricao": "Scatter PC1×PC2: agrupa estruturas por assinatura de AA",
        "requer_dssp": False,
        "requer_evo": False,
    },
    "helix_content_distribution": {
        "nome": "Distribuição de conteúdo helicoidal",
        "grupo": "4. Composição e Variabilidade",
        "descricao": "Histograma + violin: % de resíduos helicais por proteína",
        "requer_dssp": True,
        "requer_evo": True,
    },
    "aa_cooccurrence": {
        "nome": "Co-ocorrência de AAs em hélices",
        "grupo": "4. Composição e Variabilidade",
        "descricao": "Heatmap 20×20: pares de AAs que co-ocorrem na mesma hélice",
        "requer_dssp": True,
        "requer_evo": True,
    },
    "helix_length_vs_composition": {
        "nome": "Comprimento vs composição",
        "grupo": "4. Composição e Variabilidade",
        "descricao": "Barras agrupadas: frequência por grupo (hidrofóbico/polar/etc) × comprimento",
        "requer_dssp": True,
        "requer_evo": True,
    },
    # ── Grupo 5: Transições e história ────────────────────────────────────────
    "helix_transition_matrix": {
        "nome": "Matriz de transição H→G→I",
        "grupo": "5. Transições e História Evolutiva",
        "descricao": "Heatmap normalizado: probabilidade de transição entre tipos de hélice",
        "requer_dssp": True,
        "requer_evo": False,
    },
    "g_ratio_by_length": {
        "nome": "Razão 3-10/(H+G+I) por comprimento",
        "grupo": "5. Transições e História Evolutiva",
        "descricao": "Linha + barras: hélices curtas são predominantemente 3-10?",
        "requer_dssp": True,
        "requer_evo": False,
    },
    "shannon_entropy_heptad": {
        "nome": "Entropia de Shannon por posição do heptad",
        "grupo": "5. Transições e História Evolutiva",
        "descricao": "Barras coloridas + scatter: conservação por posição a-g",
        "requer_dssp": False,
        "requer_evo": True,
    },
    # ── Grupo 6: Viés do código genético ──────────────────────────────────────
    "codon_degeneracy_vs_propensity": {
        "nome": "Degenerescência de códons vs propensão",
        "grupo": "6. Viés do Código Genético",
        "descricao": "Scatter: nº de códons vs propensão teórica e observada para hélice",
        "requer_dssp": True,
        "requer_evo": False,
    },
    "proteome_comparison": {
        "nome": "Enriquecimento vs proteoma humano",
        "grupo": "6. Viés do Código Genético",
        "descricao": "Barras de enriquecimento + scatter obs vs proteoma de referência",
        "requer_dssp": False,
        "requer_evo": False,
    },
}


def _clear() -> None:
    import os
    os.system('clear' if os.name == 'posix' else 'cls')


def _print_header() -> None:
    print("\n" + "=" * 65)
    print("  CATH ANALYSIS – SELETOR DE GRAFICOS EVOLUTIVOS")
    print("=" * 65)


def _print_catalog(selected: set) -> None:
    """Imprime o catálogo agrupado com status de seleção."""
    current_group = None
    items = list(PLOT_CATALOG.items())

    for i, (plot_id, info) in enumerate(items, start=1):
        if info["grupo"] != current_group:
            current_group = info["grupo"]
            print(f"\n  {current_group}")
            print("  " + "-" * 55)
        mark = "x" if plot_id in selected else " "
        dssp_tag = "[DSSP]" if info["requer_dssp"] else "      "
        evo_tag  = "[EVO]" if info["requer_evo"]  else "     "
        print(f"  [{mark}] {i:>2}. {dssp_tag}{evo_tag} {info['nome']}")


def _print_legend() -> None:
    print("\n  Legenda:")
    print("  [DSSP] = requer DSSP instalado (mkdssp)")
    print("  [EVO]  = requer etapa de coleta evolutiva (mais lenta)")


def _print_commands() -> None:
    print("\n" + "-" * 65)
    print("  Comandos:")
    print("  <numero>      Selecionar/deselecionar grafico")
    print("  g<numero>     Selecionar/deselecionar grupo inteiro (ex: g3)")
    print("  all           Selecionar todos")
    print("  none          Desmarcar todos")
    print("  info <num>    Ver descricao detalhada")
    print("  run           Gerar graficos selecionados")
    print("  q             Sair sem plotar")
    print("-" * 65)


def select_plots(
    preselect_all: bool = False,
    preselected: set | None = None,
) -> list:
    """Exibe o menu interativo e retorna lista de plot_ids selecionados.

    Args:
        preselect_all: Se True, inicia com todos selecionados.
        preselected: Conjunto de plot_ids a pré-selecionar (ex: via --plots).

    Returns:
        Lista de plot_ids escolhidos pelo usuário.
    """
    try:
        initial: set = set()
        if preselect_all:
            initial = set(PLOT_CATALOG.keys())
        elif preselected:
            initial = preselected & set(PLOT_CATALOG.keys())
        return _run_menu(initial)
    except KeyboardInterrupt:
        print("\n\n  Interrompido pelo usuário.")
        return []


def _run_menu(initial_selected: set | None = None) -> list:
    items = list(PLOT_CATALOG.items())
    selected: set = initial_selected.copy() if initial_selected else set()

    # Agrupa por grupo para comando 'g'
    groups: dict = {}
    for plot_id, info in PLOT_CATALOG.items():
        groups.setdefault(info["grupo"], []).append(plot_id)
    group_list = list(groups.keys())

    while True:
        _print_header()
        _print_catalog(selected)
        _print_legend()
        _print_commands()
        print(f"\n  Selecionados: {len(selected)}/{len(PLOT_CATALOG)}")
        print()

        try:
            cmd = input("  > ").strip().lower()
        except KeyboardInterrupt:
            return []

        if cmd == "q":
            return []
        elif cmd == "run":
            if not selected:
                print("\n  Nenhum grafico selecionado. Selecione ao menos um.")
                input("  [Enter para continuar]")
            else:
                return [pid for pid, _ in items if pid in selected]
        elif cmd == "all":
            selected = set(PLOT_CATALOG.keys())
        elif cmd == "none":
            selected = set()
        elif cmd.startswith("info "):
            try:
                n = int(cmd.split()[1]) - 1
                if 0 <= n < len(items):
                    pid, info = items[n]
                    print(f"\n  {'-'*55}")
                    print(f"  {info['nome']}")
                    print(f"  Grupo: {info['grupo']}")
                    print(f"  Descricao: {info['descricao']}")
                    print(f"  Requer DSSP: {'Sim' if info['requer_dssp'] else 'Nao'}")
                    print(f"  Requer coleta evolutiva: {'Sim' if info['requer_evo'] else 'Nao'}")
                    print(f"  {'-'*55}")
                    input("  [Enter para continuar]")
            except (ValueError, IndexError):
                pass
        elif cmd.startswith("g") and cmd[1:].isdigit():
            g_idx = int(cmd[1:]) - 1
            if 0 <= g_idx < len(group_list):
                group_ids = groups[group_list[g_idx]]
                if all(pid in selected for pid in group_ids):
                    selected -= set(group_ids)
                else:
                    selected |= set(group_ids)
        elif cmd.isdigit():
            n = int(cmd) - 1
            if 0 <= n < len(items):
                pid = items[n][0]
                if pid in selected:
                    selected.discard(pid)
                else:
                    selected.add(pid)
        # Suporte a múltiplos números: "1 3 5"
        elif all(part.isdigit() for part in cmd.split()):
            for part in cmd.split():
                n = int(part) - 1
                if 0 <= n < len(items):
                    pid = items[n][0]
                    if pid in selected:
                        selected.discard(pid)
                    else:
                        selected.add(pid)
