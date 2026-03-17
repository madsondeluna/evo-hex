"""
Funções de visualização para o pipeline CATH.

Todas as funções salvam o plot em disco (dpi=300) e fecham a figura,
evitando acúmulo de memória em processamentos longos.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

# Paleta padrão para tipos de hélices H/G/I
_HELIX_COLORS = ["steelblue", "coral", "lightgreen"]

sns.set_style("whitegrid")
plt.rcParams["figure.figsize"] = (12, 6)


def _save_and_close(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Gráfico salvo em: {path}")


# ── Análise avançada de hélices (helix_analysis.py) ───────────────────────────

def plot_helix_propensities(df: pd.DataFrame, analysis_path: Path) -> None:
    """Bar chart de propensão observada vs Chou-Fasman + scatter freq total×hélice.

    Args:
        df: Retorno de calculate_helix_propensities.
        analysis_path: Diretório de saída.
    """
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    x = np.arange(len(df))
    width = 0.35

    axes[0].bar(x - width / 2, df["Propensity_Observed"], width,
                label="Propensão Observada", color="steelblue", edgecolor="black")
    axes[0].bar(x + width / 2, df["Propensity_Theoretical"], width,
                label="Propensão Teórica (Chou-Fasman)", color="coral", edgecolor="black")
    axes[0].axhline(y=1.0, color="red", linestyle="--", alpha=0.5, label="Propensão neutra")
    axes[0].set_xlabel("Aminoácido")
    axes[0].set_ylabel("Propensão para Hélice")
    axes[0].set_title("Propensão de Aminoácidos para Hélices Alpha")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(df["AA"], rotation=0)
    axes[0].legend()
    axes[0].grid(axis="y", alpha=0.3)

    axes[1].scatter(df["Freq_Total"], df["Freq_Helix"], s=100, alpha=0.6, color="steelblue")
    for _, row in df.iterrows():
        axes[1].annotate(row["AA"], (row["Freq_Total"], row["Freq_Helix"]),
                         fontsize=9, ha="center")

    max_val = max(df["Freq_Total"].max(), df["Freq_Helix"].max())
    axes[1].plot([0, max_val], [0, max_val], "r--", alpha=0.5, label="Frequência igual")
    axes[1].set_xlabel("Frequência Total (%)")
    axes[1].set_ylabel("Frequência em Hélices (%)")
    axes[1].set_title("Frequência de Aminoácidos: Total vs Hélices")
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    fig.tight_layout()
    _save_and_close(fig, analysis_path / "helix_propensities.png")


def plot_helix_positions(df: pd.DataFrame, analysis_path: Path) -> None:
    """Heatmap de preferência de aminoácidos por posição na hélice.

    Args:
        df: Retorno de analyze_helix_positions.
        analysis_path: Diretório de saída.
    """
    pivot = df.pivot(index="AA", columns="Position", values="Frequency")
    pivot = pivot[["N-terminal", "Middle", "C-terminal"]]

    fig, ax = plt.subplots(figsize=(8, 10))
    sns.heatmap(pivot, annot=True, fmt=".1f", cmap="YlOrRd",
                cbar_kws={"label": "Frequência (%)"}, ax=ax)
    ax.set_title("Preferência de Aminoácidos por Posição na Hélice")
    ax.set_ylabel("Aminoácido")
    ax.set_xlabel("Posição")
    fig.tight_layout()
    _save_and_close(fig, analysis_path / "helix_positions.png")


def plot_helix_length_distribution(helix_lengths: list[int], analysis_path: Path) -> None:
    """Histograma e boxplot da distribuição de comprimentos de hélices.

    Args:
        helix_lengths: Lista de comprimentos em resíduos.
        analysis_path: Diretório de saída.
    """
    arr = np.array(helix_lengths)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].hist(arr, bins=50, color="steelblue", edgecolor="black", alpha=0.7)
    axes[0].set_xlabel("Comprimento (resíduos)")
    axes[0].set_ylabel("Frequência")
    axes[0].set_title("Distribuição de Comprimentos de Hélices")
    axes[0].grid(axis="y", alpha=0.3)

    axes[1].boxplot(arr, vert=True)
    axes[1].set_ylabel("Comprimento (resíduos)")
    axes[1].set_title("Box Plot – Comprimento de Hélices")
    axes[1].grid(axis="y", alpha=0.3)

    fig.tight_layout()
    _save_and_close(fig, analysis_path / "helix_lengths.png")

    print(f"\nEstatísticas de comprimento de hélices:")
    print(f"  Média:   {arr.mean():.1f} resíduos")
    print(f"  Mediana: {np.median(arr):.1f} resíduos")
    print(f"  DP:      {arr.std():.1f}")
    print(f"  Mín/Máx: {arr.min()} / {arr.max()}")


def plot_hydrophobic_heptad(heptad_data: dict, analysis_path: Path) -> None:
    """Frequência de resíduos hidrofóbicos por posição no heptad.

    Args:
        heptad_data: Retorno de analyze_hydrophobic_patterns.
        analysis_path: Diretório de saída.
    """
    positions = list("abcdefg")
    hydrophobic_freq = [
        heptad_data[i]["hydrophobic"] / heptad_data[i]["total"] * 100
        if heptad_data[i]["total"] > 0 else 0.0
        for i in range(7)
    ]
    mean_freq = np.mean(hydrophobic_freq)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(positions, hydrophobic_freq, color="steelblue", edgecolor="black")
    ax.axhline(y=mean_freq, color="red", linestyle="--", alpha=0.5,
               label=f"Média: {mean_freq:.1f}%")
    ax.set_xlabel("Posição no Heptad Repeat")
    ax.set_ylabel("Frequência de Aminoácidos Hidrofóbicos (%)")
    ax.set_title("Padrão Heptad: Distribuição de Resíduos Hidrofóbicos")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    _save_and_close(fig, analysis_path / "heptad_pattern.png")


# ── Tipos de hélices (helix_types.py) ─────────────────────────────────────────

def plot_helix_type_distribution(helix_counts: dict, analysis_path: Path) -> None:
    """Bar chart + pie chart da distribuição de tipos de hélices.

    Args:
        helix_counts: Dict {tipo: n_hélices}.
        analysis_path: Diretório de saída.
    """
    from .config import HELIX_TYPES  # noqa: PLC0415

    ordered = [(t, helix_counts[t]) for t in ("H", "G", "I") if t in helix_counts]
    names = [HELIX_TYPES[t] for t, _ in ordered]
    counts = [n for _, n in ordered]
    total = sum(counts)
    percentages = [c / total * 100 for c in counts]
    colors = _HELIX_COLORS[: len(ordered)]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    axes[0].bar(names, counts, color=colors, edgecolor="black")
    axes[0].set_ylabel("Número de Hélices")
    axes[0].set_title("Distribuição de Tipos de Hélices")
    axes[0].grid(axis="y", alpha=0.3)
    for i, (name, count, pct) in enumerate(zip(names, counts, percentages)):
        axes[0].text(i, count, f"{count:,}\n({pct:.1f}%)",
                     ha="center", va="bottom", fontweight="bold")

    axes[1].pie(counts, labels=names, autopct="%1.1f%%", colors=colors,
                startangle=90, wedgeprops={"edgecolor": "black"})
    axes[1].set_title("Proporção de Tipos de Hélices")

    fig.tight_layout()
    _save_and_close(fig, analysis_path / "helix_type_distribution.png")


def plot_helix_type_composition(composition_df: pd.DataFrame, analysis_path: Path) -> None:
    """Heatmap comparativo de composição por tipo de hélice.

    Args:
        composition_df: Retorno de analyze_helix_type_composition.
        analysis_path: Diretório de saída.
    """
    pivot = composition_df.pivot(index="AA", columns="Helix_Name", values="Frequency")
    if "Alpha helix" in pivot.columns:
        pivot = pivot.sort_values("Alpha helix", ascending=False)

    fig, ax = plt.subplots(figsize=(10, 12))
    sns.heatmap(pivot, annot=True, fmt=".1f", cmap="YlOrRd",
                cbar_kws={"label": "Frequência (%)"}, linewidths=0.5, ax=ax)
    ax.set_title("Composição de Aminoácidos por Tipo de Hélice", fontsize=14, fontweight="bold")
    ax.set_ylabel("Aminoácido", fontsize=12)
    ax.set_xlabel("Tipo de Hélice", fontsize=12)
    fig.tight_layout()
    _save_and_close(fig, analysis_path / "helix_type_composition_heatmap.png")


def plot_top_amino_acids_by_type(comparison_df: pd.DataFrame, analysis_path: Path) -> None:
    """Barras horizontais dos top-10 aminoácidos por tipo de hélice.

    Args:
        comparison_df: Retorno de compare_helix_types.
        analysis_path: Diretório de saída.
    """
    helix_names = comparison_df["Helix_Name"].unique()
    n_types = len(helix_names)
    colors = _HELIX_COLORS

    fig, axes = plt.subplots(1, n_types, figsize=(6 * n_types, 6))
    if n_types == 1:
        axes = [axes]

    for idx, h_name in enumerate(helix_names):
        data = comparison_df[comparison_df["Helix_Name"] == h_name].nsmallest(10, "Rank")
        axes[idx].barh(data["AA"], data["Frequency"],
                       color=colors[idx % len(colors)], edgecolor="black")
        axes[idx].set_xlabel("Frequência (%)")
        axes[idx].set_title(f"{h_name}\nTop 10 Aminoácidos", fontweight="bold")
        axes[idx].invert_yaxis()
        axes[idx].grid(axis="x", alpha=0.3)
        for aa, freq in zip(data["AA"], data["Frequency"]):
            axes[idx].text(freq, aa, f" {freq:.1f}%", va="center")

    fig.tight_layout()
    _save_and_close(fig, analysis_path / "helix_type_top_amino_acids.png")


def plot_helix_length_comparison(helix_lengths: dict, analysis_path: Path) -> None:
    """Boxplot e histograma comparativo de comprimentos por tipo.

    Args:
        helix_lengths: Dict {tipo: [comprimentos]}.
        analysis_path: Diretório de saída.
    """
    from .config import HELIX_TYPES  # noqa: PLC0415

    types_to_plot = [t for t in ("H", "G", "I") if t in helix_lengths]
    data_to_plot = [helix_lengths[t] for t in types_to_plot]
    labels = [HELIX_TYPES[t] for t in types_to_plot]
    colors = _HELIX_COLORS[: len(types_to_plot)]

    fig, axes = plt.subplots(2, 1, figsize=(12, 10))

    bp = axes[0].boxplot(data_to_plot, labels=labels, patch_artist=True)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
    axes[0].set_ylabel("Comprimento (resíduos)")
    axes[0].set_title("Distribuição de Comprimentos por Tipo de Hélice")
    axes[0].grid(axis="y", alpha=0.3)

    for h_type, color in zip(types_to_plot, colors):
        axes[1].hist(helix_lengths[h_type], bins=30, alpha=0.5,
                     label=HELIX_TYPES[h_type], color=color, edgecolor="black")
    axes[1].set_xlabel("Comprimento (resíduos)")
    axes[1].set_ylabel("Frequência")
    axes[1].set_title("Histograma de Comprimentos")
    axes[1].legend()
    axes[1].grid(axis="y", alpha=0.3)

    fig.tight_layout()
    _save_and_close(fig, analysis_path / "helix_length_by_type.png")


def plot_statistical_differences(stat_df: pd.DataFrame, analysis_path: Path) -> None:
    """Diferenças de frequência entre Alpha e 3-10 helix.

    Args:
        stat_df: Retorno de statistical_comparison.
        analysis_path: Diretório de saída.
    """
    if stat_df.empty:
        return

    top_diff = stat_df.head(10)
    colors = ["steelblue" if x > 0 else "coral" for x in top_diff["Difference"]]

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    axes[0].barh(top_diff["AA"], top_diff["Difference"], color=colors, edgecolor="black")
    axes[0].axvline(x=0, color="black", linewidth=0.8)
    axes[0].set_xlabel("Diferença de Frequência (%)")
    axes[0].set_title("Diferenças Alpha helix vs 3-10 helix", fontweight="bold")
    axes[0].invert_yaxis()
    axes[0].grid(axis="x", alpha=0.3)
    legend_elements = [
        Patch(facecolor="steelblue", label="Enriquecido em Alpha"),
        Patch(facecolor="coral", label="Enriquecido em 3-10"),
    ]
    axes[0].legend(handles=legend_elements)

    axes[1].scatter(top_diff["Alpha_Freq"], top_diff["3-10_Freq"],
                    s=200, alpha=0.6, c=colors, edgecolor="black", linewidth=2)
    for _, row in top_diff.iterrows():
        axes[1].annotate(row["AA"], (row["Alpha_Freq"], row["3-10_Freq"]),
                         fontsize=10, ha="center", fontweight="bold")

    max_val = max(top_diff["Alpha_Freq"].max(), top_diff["3-10_Freq"].max())
    axes[1].plot([0, max_val], [0, max_val], "k--", alpha=0.3, label="Frequência igual")
    axes[1].set_xlabel("Frequência em Alpha helix (%)")
    axes[1].set_ylabel("Frequência em 3-10 helix (%)")
    axes[1].set_title("Comparação Direta: Alpha vs 3-10", fontweight="bold")
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    fig.tight_layout()
    _save_and_close(fig, analysis_path / "helix_type_statistical_comparison.png")


# ── Gráficos evolutivos (evolutionary_analysis.py) ─────────────────────────

def plot_helical_wheel_average(helix_sequences: list, analysis_path: Path) -> None:
    """
    Helical wheel composite: frequência de resíduos hidrofóbicos por posição angular.
    Usa projeção polar com 36 bins de 10° cada. Alpha helix: 100°/resíduo.
    """
    from .config import HYDROPHOBIC_AA

    # Para cada resíduo em cada hélice, calcula ângulo = (i * 100) % 360
    angle_hydrophobic = np.zeros(36)
    angle_total = np.zeros(36)

    for seq in helix_sequences:
        for i, aa in enumerate(seq):
            angle_deg = (i * 100) % 360
            bin_idx = int(angle_deg / 10)
            angle_total[bin_idx] += 1
            if aa in HYDROPHOBIC_AA:
                angle_hydrophobic[bin_idx] += 1

    with np.errstate(divide='ignore', invalid='ignore'):
        freq = np.where(angle_total > 0, angle_hydrophobic / angle_total, 0)

    angles = np.linspace(0, 2 * np.pi, 36, endpoint=False)

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': 'polar'})
    bars = ax.bar(angles, freq, width=2 * np.pi / 36, bottom=0.0,
                  color=plt.cm.RdYlBu_r(freq), edgecolor='black', linewidth=0.5)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(['25%', '50%', '75%', '100%'], fontsize=8)
    ax.set_title('Helical Wheel: Frequência de Resíduos\nHidrofóbicos por Posição Angular',
                 fontsize=12, fontweight='bold', pad=20)

    sm = plt.cm.ScalarMappable(cmap='RdYlBu_r', norm=plt.Normalize(0, 1))
    sm.set_array([])
    plt.colorbar(sm, ax=ax, pad=0.1, shrink=0.7, label='Frequência Hidrofóbica')

    fig.tight_layout()
    _save_and_close(fig, analysis_path / 'helical_wheel_average.png')


def plot_hydrophobic_moment_distribution(moments: list, analysis_path: Path) -> None:
    """Distribuição do momento hidrofóbico de Eisenberg por hélice."""
    arr = np.array(moments)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].hist(arr, bins=50, color='steelblue', edgecolor='black', alpha=0.8)
    axes[0].axvline(arr.mean(), color='red', linestyle='--', label=f'Média: {arr.mean():.3f}')
    axes[0].axvline(np.median(arr), color='orange', linestyle='--', label=f'Mediana: {np.median(arr):.3f}')
    axes[0].set_xlabel('Momento Hidrofóbico (μH)')
    axes[0].set_ylabel('Frequência')
    axes[0].set_title('Distribuição do Momento Hidrofóbico\n(Escala de Eisenberg)')
    axes[0].legend()
    axes[0].grid(axis='y', alpha=0.3)

    # Categorias: baixo < 0.2, médio 0.2-0.4, alto > 0.4
    low = np.sum(arr < 0.2)
    med = np.sum((arr >= 0.2) & (arr < 0.4))
    high = np.sum(arr >= 0.4)
    total = len(arr)

    axes[1].bar(['Baixo\n(<0.2)', 'Médio\n(0.2-0.4)', 'Alto\n(≥0.4)'],
                [low/total*100, med/total*100, high/total*100],
                color=['lightblue', 'steelblue', 'darkblue'], edgecolor='black')
    axes[1].set_ylabel('% de Hélices')
    axes[1].set_title('Classificação por Anfipatiticidade')
    axes[1].grid(axis='y', alpha=0.3)
    for i, v in enumerate([low/total*100, med/total*100, high/total*100]):
        axes[1].text(i, v + 0.5, f'{v:.1f}%', ha='center', fontweight='bold')

    fig.suptitle(f'Momento Hidrofóbico – {total:,} hélices alpha', fontsize=13, fontweight='bold')
    fig.tight_layout()
    _save_and_close(fig, analysis_path / 'hydrophobic_moment_distribution.png')


def plot_ncap_ccap(ncap_pos: dict, ccap_pos: dict, analysis_path: Path) -> None:
    """Heatmap de preferência de AAs nas 3 posições do N-cap e C-cap."""
    from .config import STANDARD_AMINO_ACIDS

    all_aa = sorted(STANDARD_AMINO_ACIDS)

    # Normaliza cada posição
    def _to_freq(pos_dict, positions):
        rows = {}
        for p in positions:
            cnts = pos_dict.get(p, {})
            total = sum(cnts.values()) or 1
            rows[str(p)] = {aa: cnts.get(aa, 0) / total * 100 for aa in all_aa}
        return pd.DataFrame(rows, index=all_aa).T

    ncap_df = _to_freq(ncap_pos, [0, 1, 2])
    ccap_df = _to_freq(ccap_pos, [-3, -2, -1])
    ncap_df.index = ['N1', 'N2', 'N3']
    ccap_df.index = ['C3', 'C2', 'C1']

    fig, axes = plt.subplots(1, 2, figsize=(18, 5))

    sns.heatmap(ncap_df, annot=True, fmt='.1f', cmap='Blues',
                cbar_kws={'label': 'Frequência (%)'}, linewidths=0.3, ax=axes[0])
    axes[0].set_title('N-cap: Preferência de Aminoácidos\n(posições N1, N2, N3)', fontsize=12, fontweight='bold')
    axes[0].set_xlabel('Aminoácido')
    axes[0].set_ylabel('Posição')

    sns.heatmap(ccap_df, annot=True, fmt='.1f', cmap='Oranges',
                cbar_kws={'label': 'Frequência (%)'}, linewidths=0.3, ax=axes[1])
    axes[1].set_title('C-cap: Preferência de Aminoácidos\n(posições C3, C2, C1)', fontsize=12, fontweight='bold')
    axes[1].set_xlabel('Aminoácido')
    axes[1].set_ylabel('Posição')

    fig.tight_layout()
    _save_and_close(fig, analysis_path / 'ncap_ccap_preferences.png')


def plot_pca_composition(per_structure: dict, analysis_path: Path) -> None:
    """PCA da composição de AAs por estrutura. Colorido por % do AA dominante."""
    try:
        from sklearn.preprocessing import StandardScaler
        from sklearn.decomposition import PCA
    except ImportError:
        print("  [AVISO] sklearn não disponível. Instale scikit-learn para este gráfico.")
        return

    from .config import STANDARD_AMINO_ACIDS

    all_aa = sorted(STANDARD_AMINO_ACIDS)
    codes = sorted(per_structure.keys())

    if not codes:
        print("  [AVISO] Nenhuma estrutura disponível para PCA.")
        return

    # Matriz N_structures × 20
    matrix = []
    for code in codes:
        counts = per_structure[code]
        total = sum(counts.values()) or 1
        matrix.append([counts.get(aa, 0) / total for aa in all_aa])

    X = np.array(matrix)
    X_scaled = StandardScaler().fit_transform(X)

    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(X_scaled)

    # Cor = AA dominante
    dominant_idx = X.argmax(axis=1)

    fig, ax = plt.subplots(figsize=(10, 8))
    scatter = ax.scatter(coords[:, 0], coords[:, 1], c=dominant_idx,
                         cmap='tab20', alpha=0.6, s=20, edgecolors='none')
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)')
    ax.set_title(f'PCA – Composição de Aminoácidos por Estrutura\n({len(codes):,} estruturas)')
    ax.grid(alpha=0.3)

    # Legenda compacta: top 10 AAs dominantes
    unique, cnts_dom = np.unique(dominant_idx, return_counts=True)
    top_idx = unique[np.argsort(-cnts_dom)[:10]]
    handles = [plt.Line2D([0], [0], marker='o', color='w',
                           markerfacecolor=plt.cm.tab20(i / 20), markersize=8,
                           label=all_aa[i]) for i in top_idx]
    ax.legend(handles=handles, title='AA dominante', loc='best', fontsize=8)

    fig.tight_layout()
    _save_and_close(fig, analysis_path / 'pca_aa_composition.png')


def plot_helix_content_distribution(helix_contents: list, analysis_path: Path) -> None:
    """Distribuição do % de conteúdo helicoidal por estrutura."""
    arr = np.array(helix_contents) * 100

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].hist(arr, bins=40, color='steelblue', edgecolor='black', alpha=0.8)
    axes[0].axvline(arr.mean(), color='red', linestyle='--', label=f'Média: {arr.mean():.1f}%')
    axes[0].axvline(np.median(arr), color='orange', linestyle='--', label=f'Mediana: {np.median(arr):.1f}%')
    axes[0].set_xlabel('Conteúdo Helicoidal (%)')
    axes[0].set_ylabel('Número de Estruturas')
    axes[0].set_title('Distribuição de Conteúdo Helicoidal\npor Estrutura')
    axes[0].legend()
    axes[0].grid(axis='y', alpha=0.3)

    # Violin plot
    axes[1].violinplot(arr, showmedians=True, showextrema=True)
    axes[1].set_xticks([1])
    axes[1].set_xticklabels(['Estruturas'])
    axes[1].set_ylabel('Conteúdo Helicoidal (%)')
    axes[1].set_title('Violin Plot – Conteúdo Helicoidal')
    axes[1].grid(axis='y', alpha=0.3)

    # Anota percentis
    for p, label in [(25, 'Q1'), (50, 'Med'), (75, 'Q3')]:
        val = np.percentile(arr, p)
        axes[1].axhline(val, color='gray', linestyle=':', alpha=0.7)
        axes[1].text(1.3, val, f'{label}: {val:.1f}%', va='center', fontsize=9)

    fig.suptitle(f'Conteúdo Helicoidal – {len(arr):,} estruturas', fontsize=13, fontweight='bold')
    fig.tight_layout()
    _save_and_close(fig, analysis_path / 'helix_content_distribution.png')


def plot_aa_cooccurrence(cooc_df: pd.DataFrame, analysis_path: Path) -> None:
    """Heatmap de co-ocorrência de AAs dentro das hélices."""
    fig, ax = plt.subplots(figsize=(12, 10))

    mask = np.eye(len(cooc_df), dtype=bool)  # mascara diagonal
    sns.heatmap(cooc_df * 100, annot=False, cmap='YlOrRd', mask=mask,
                cbar_kws={'label': 'Frequência de Co-ocorrência (×100)'}, ax=ax,
                linewidths=0.2)
    ax.set_title('Co-ocorrência de Aminoácidos em Hélices Alpha\n'
                 '(frequência relativa de pares na mesma hélice)',
                 fontsize=12, fontweight='bold')
    ax.set_xlabel('Aminoácido')
    ax.set_ylabel('Aminoácido')
    fig.tight_layout()
    _save_and_close(fig, analysis_path / 'aa_cooccurrence.png')


def plot_helix_length_vs_composition(length_comp_df: pd.DataFrame, analysis_path: Path) -> None:
    """Composição por propriedade físico-química em função do comprimento da hélice."""
    bins = ['Curta (4-9)', 'Média (10-19)', 'Longa (≥20)']
    groups = length_comp_df['Group'].unique()
    colors_map = {
        'Hidrofóbico': 'steelblue', 'Polar': 'lightgreen',
        'Carregado+': 'coral', 'Carregado-': 'gold', 'Especial': 'violet'
    }

    x = np.arange(len(bins))
    width = 0.15

    fig, ax = plt.subplots(figsize=(12, 6))

    for i, group in enumerate(groups):
        data = length_comp_df[length_comp_df['Group'] == group]
        freqs = [data[data['Bin'] == b]['Frequency'].values[0]
                 if len(data[data['Bin'] == b]) > 0 else 0
                 for b in bins]
        offset = (i - len(groups)/2 + 0.5) * width
        ax.bar(x + offset, freqs, width, label=group,
               color=colors_map.get(group, 'gray'), edgecolor='black', alpha=0.85)

    ax.set_xlabel('Comprimento da Hélice')
    ax.set_ylabel('Frequência (%)')
    ax.set_title('Composição de Aminoácidos por Comprimento de Hélice')
    ax.set_xticks(x)
    ax.set_xticklabels(bins)
    ax.legend(title='Grupo')
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    _save_and_close(fig, analysis_path / 'helix_length_vs_composition.png')


def plot_helix_transition_matrix(transitions, analysis_path: Path) -> None:
    """Heatmap da matriz de transição entre tipos de hélice (H↔G↔I)."""
    from .config import HELIX_TYPES

    all_types = ['H', 'G', 'I']

    # Monta matriz 3×3
    matrix = pd.DataFrame(0.0, index=all_types, columns=all_types)
    for (src, dst), count in transitions.items():
        if src in all_types and dst in all_types:
            matrix.loc[src, dst] = count

    # Normaliza por linha
    row_sums = matrix.sum(axis=1)
    matrix_norm = matrix.div(row_sums.replace(0, 1), axis=0) * 100
    matrix_norm.index = [HELIX_TYPES[t] for t in all_types]
    matrix_norm.columns = [HELIX_TYPES[t] for t in all_types]

    fig, ax = plt.subplots(figsize=(7, 6))
    sns.heatmap(matrix_norm, annot=True, fmt='.1f', cmap='Blues',
                cbar_kws={'label': '% de Transições'}, linewidths=1,
                annot_kws={'size': 13, 'weight': 'bold'}, ax=ax)
    ax.set_title('Matriz de Transição entre Tipos de Hélice\n'
                 '(% de transições por tipo de origem)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Para')
    ax.set_ylabel('De')
    fig.tight_layout()
    _save_and_close(fig, analysis_path / 'helix_transition_matrix.png')


def plot_g_ratio_by_length(helix_lengths_by_type: dict, analysis_path: Path) -> None:
    """Razão 3-10/(Alpha+3-10+Pi) em função do comprimento da hélice."""
    h_len = helix_lengths_by_type.get('H', [])
    g_len = helix_lengths_by_type.get('G', [])
    i_len = helix_lengths_by_type.get('I', [])

    all_lengths = list(h_len) + list(g_len) + list(i_len)
    if not all_lengths:
        return

    max_len = min(max(all_lengths), 50)
    bins = list(range(3, max_len + 2))

    h_counts = np.histogram(h_len, bins=bins)[0] if h_len else np.zeros(len(bins)-1)
    g_counts = np.histogram(g_len, bins=bins)[0] if g_len else np.zeros(len(bins)-1)
    i_counts = np.histogram(i_len, bins=bins)[0] if i_len else np.zeros(len(bins)-1)
    total = h_counts + g_counts + i_counts

    with np.errstate(divide='ignore', invalid='ignore'):
        g_ratio = np.where(total > 0, g_counts / total * 100, np.nan)

    bin_centers = [(bins[i] + bins[i+1]) / 2 for i in range(len(bins)-1)]

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(bin_centers, g_ratio, 'o-', color='coral', linewidth=2, markersize=4)
    ax.fill_between(bin_centers, 0, g_ratio, alpha=0.2, color='coral')
    ax.set_xlabel('Comprimento da Hélice (resíduos)')
    ax.set_ylabel('% de Hélices 3-10 G/(H+G+I)')
    ax.set_title('Razão de Hélices 3-10 por Comprimento\n'
                 '(hélices curtas tendem a ser 3-10)')
    ax.grid(alpha=0.3)
    ax.set_xlim(bins[0], bins[-1])

    # Barra de contagem abaixo
    ax2 = ax.twinx()
    ax2.bar(bin_centers, total, width=0.8, color='lightgray', alpha=0.4, label='Total de hélices')
    ax2.set_ylabel('Total de Hélices', color='gray')
    ax2.tick_params(axis='y', labelcolor='gray')

    fig.tight_layout()
    _save_and_close(fig, analysis_path / 'g_ratio_by_length.png')


def plot_shannon_entropy_heptad(entropy_df: pd.DataFrame, analysis_path: Path) -> None:
    """Entropia de Shannon por posição do heptad repeat."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    colors = plt.cm.RdYlGn_r(entropy_df['Entropy'] / entropy_df['Entropy'].max())

    bars = axes[0].bar(entropy_df['Position'], entropy_df['Entropy'],
                       color=colors, edgecolor='black')
    axes[0].set_xlabel('Posição no Heptad (a-g)')
    axes[0].set_ylabel('Entropia de Shannon (bits)')
    axes[0].set_title('Conservação por Posição no Heptad\n(baixa entropia = maior conservação)')
    axes[0].grid(axis='y', alpha=0.3)

    # Anota AA mais frequente
    for bar, (_, row) in zip(bars, entropy_df.iterrows()):
        axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                     f"{row['Top_AA']}\n{row['Top_Freq']:.0f}%",
                     ha='center', va='bottom', fontsize=9, fontweight='bold')

    # Scatter: posição × entropia com tamanho = top_freq
    axes[1].scatter(entropy_df['Position'], entropy_df['Entropy'],
                    s=entropy_df['Top_Freq'] * 5, c=entropy_df['Entropy'],
                    cmap='RdYlGn_r', edgecolors='black', linewidth=1.5)
    for _, row in entropy_df.iterrows():
        axes[1].annotate(row['Top_AA'], (row['Position'], row['Entropy']),
                        fontsize=11, ha='center', va='bottom', fontweight='bold')
    axes[1].set_xlabel('Posição no Heptad (a-g)')
    axes[1].set_ylabel('Entropia de Shannon (bits)')
    axes[1].set_title('AA mais frequente × Entropia\n(tamanho do ponto = frequência do AA)')
    axes[1].grid(alpha=0.3)

    fig.tight_layout()
    _save_and_close(fig, analysis_path / 'shannon_entropy_heptad.png')


def plot_codon_degeneracy_vs_propensity(codon_df: pd.DataFrame, analysis_path: Path) -> None:
    """Scatter: degenerescência de códons × propensão para hélice."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    colors = plt.cm.Set2(np.linspace(0, 1, len(codon_df['Codon_Degeneracy'].unique())))
    color_map = {d: colors[i] for i, d in enumerate(sorted(codon_df['Codon_Degeneracy'].unique()))}

    for _, row in codon_df.iterrows():
        c = color_map[row['Codon_Degeneracy']]
        axes[0].scatter(row['Codon_Degeneracy'], row['Propensity_Theoretical'],
                       color=c, s=100, edgecolors='black', zorder=5)
        axes[0].annotate(row['AA'], (row['Codon_Degeneracy'], row['Propensity_Theoretical']),
                        fontsize=8, ha='left', va='bottom')

    axes[0].axhline(1.0, color='red', linestyle='--', alpha=0.5, label='Propensão neutra')
    axes[0].set_xlabel('Degenerescência de Códons (nº de códons)')
    axes[0].set_ylabel('Propensão Teórica (Chou-Fasman)')
    axes[0].set_title('Degenerescência de Códons\nvs Propensão para Hélice (Teórica)')
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    for _, row in codon_df.iterrows():
        c = color_map[row['Codon_Degeneracy']]
        axes[1].scatter(row['Codon_Degeneracy'], row['Propensity_Observed'],
                       color=c, s=100, edgecolors='black', zorder=5)
        axes[1].annotate(row['AA'], (row['Codon_Degeneracy'], row['Propensity_Observed']),
                        fontsize=8, ha='left', va='bottom')

    axes[1].axhline(1.0, color='red', linestyle='--', alpha=0.5, label='Propensão neutra')
    axes[1].set_xlabel('Degenerescência de Códons')
    axes[1].set_ylabel('Propensão Observada')
    axes[1].set_title('Degenerescência de Códons\nvs Propensão para Hélice (Observada)')
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    fig.suptitle('Viés do Código Genético na Seleção de Aminoácidos em Hélices',
                 fontsize=13, fontweight='bold')
    fig.tight_layout()
    _save_and_close(fig, analysis_path / 'codon_degeneracy_vs_propensity.png')


def plot_proteome_comparison(proteome_df: pd.DataFrame, analysis_path: Path) -> None:
    """Enriquecimento de AAs em hélices vs proteoma humano de referência."""
    df = proteome_df.sort_values('Enrichment', ascending=False)

    colors = ['steelblue' if x >= 1 else 'coral' for x in df['Enrichment']]

    fig, axes = plt.subplots(2, 1, figsize=(14, 10))

    # Barras de enriquecimento
    axes[0].bar(df['AA'], df['Enrichment'], color=colors, edgecolor='black')
    axes[0].axhline(1.0, color='black', linestyle='--', linewidth=1.5, label='Sem enriquecimento')
    axes[0].set_xlabel('Aminoácido')
    axes[0].set_ylabel('Enriquecimento (Obs/Proteoma)')
    axes[0].set_title('Enriquecimento de Aminoácidos em Hélices\nrelativo ao Proteoma Humano de Referência')
    axes[0].legend()
    axes[0].grid(axis='y', alpha=0.3)

    legend_elements = [Patch(facecolor='steelblue', label='Enriquecido em hélices'),
                       Patch(facecolor='coral', label='Depleto em hélices')]
    axes[0].legend(handles=legend_elements)

    # Scatter: freq proteoma vs freq observada
    axes[1].scatter(df['Freq_Proteome'], df['Freq_Observed'], s=100,
                   c=colors, edgecolors='black', zorder=5)
    for _, row in df.iterrows():
        axes[1].annotate(row['AA'], (row['Freq_Proteome'], row['Freq_Observed']),
                        fontsize=9, ha='left', va='bottom')

    max_val = max(df['Freq_Proteome'].max(), df['Freq_Observed'].max()) * 1.05
    axes[1].plot([0, max_val], [0, max_val], 'k--', alpha=0.3, label='Freq. igual')
    axes[1].set_xlabel('Frequência no Proteoma (%)')
    axes[1].set_ylabel('Frequência Observada em Hélices (%)')
    axes[1].set_title('Frequência Observada vs Proteoma de Referência')
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    fig.tight_layout()
    _save_and_close(fig, analysis_path / 'proteome_comparison.png')
