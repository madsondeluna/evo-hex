"""
cath_analysis – Pipeline de análise de estruturas CATH Mainly-Alpha.

Módulos
-------
config          : Constantes e caminhos (edite BASE_PATH aqui)
downloader      : Download do índice CATH e estruturas PDB
cleaner         : Limpeza de PDBs (cadeia A, apenas AA padrão)
frequency_analysis : Frequência global de aminoácidos
helix_analysis  : Análise avançada de hélices com DSSP
helix_types     : Classificação por tipo de hélice (H/G/I)
plotting        : Visualizações (matplotlib/seaborn)
"""
