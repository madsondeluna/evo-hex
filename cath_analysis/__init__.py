"""
cath_analysis – Pipeline de análise de estruturas Evo-Hex.

Módulos
-------
config             : Constantes e caminhos (edite BASE_PATH aqui)
downloader         : Download do índice CATH e estruturas PDB
cleaner            : Limpeza de PDBs (cadeia A, apenas AA padrão)
frequency_analysis : Frequência global de aminoácidos
unified_dssp       : Passo DSSP unificado (etapas 4+5+6 em uma passagem)
evolutionary_analysis : Funções compute_* para gráficos evolutivos
plotting           : Visualizações (matplotlib/seaborn)
"""
