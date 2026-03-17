# Evo-Hex — Pipeline de Análise Estrutural e Evolutiva

Pipeline modular em Python para download, limpeza e análise de estruturas proteicas
da classe **Mainly Alpha** do banco de dados CATH. Gera 21 visualizações cobrindo
desde composição de aminoácidos em hélices até análises evolutivas com base em
restrições estruturais e viés do código genético.

---

## Índice

1. [O que é o CATH e por que Mainly-Alpha?](#1-o-que-é-o-cath-e-por-que-mainly-alpha)
2. [Instalação](#2-instalação)
3. [Configuração](#3-configuração)
4. [Como executar](#4-como-executar)
5. [Etapas do pipeline](#5-etapas-do-pipeline)
6. [Catálogo de análises e gráficos](#6-catálogo-de-análises-e-gráficos)
7. [Estrutura de arquivos gerados](#7-estrutura-de-arquivos-gerados)
8. [Estrutura do código](#8-estrutura-do-código)
9. [Dependências externas](#9-dependências-externas)
10. [Perguntas frequentes](#10-perguntas-frequentes)

---

## 1. O que é o CATH e por que Mainly-Alpha?

**CATH** (Class, Architecture, Topology, Homology) é um banco de dados hierárquico
de domínios proteicos com estrutura 3D determinada experimentalmente. Cada domínio é
classificado pela composição de sua estrutura secundária.

A classe **1 — Mainly Alpha** agrupa domínios cuja estrutura é dominada por hélices
alpha (≥60% de resíduos em conformação helicoidal). Essa classe é biologicamente
relevante por incluir:

- Proteínas transportadoras de oxigênio (globinas)
- Receptores acoplados à proteína G (GPCRs)
- Canais iônicos e proteínas de membrana
- Fatores de transcrição com motivos helix-turn-helix
- Proteínas de coiled-coil envolvidas em sinalização

Estudar essa classe em larga escala permite responder perguntas sobre **restrições
evolutivas** na composição de hélices: quais aminoácidos são selecionados pela
estrutura? Qual o papel do código genético nessa seleção? Como diferentes tipos de
hélice evoluem?

---

## 2. Instalação

### Pré-requisitos

- Python 3.10 ou superior
- `mkdssp` (obrigatório para análises de estrutura secundária)

### Instalar dependências Python

```bash
pip install -r requirements.txt
```

Para o gráfico de PCA também é necessário:

```bash
pip install scikit-learn
```

### Instalar DSSP

O DSSP é um programa externo que determina a estrutura secundária a partir de
coordenadas atômicas. É necessário para as análises das etapas 4, 5 e 6.

**macOS:**
```bash
brew install brewsci/bio/dssp
```

**Ubuntu/Debian:**
```bash
sudo apt install dssp
```

**Manual:** https://swift.cmbi.umcn.nl/gv/dssp/

Verifique a instalação:
```bash
which mkdssp   # deve retornar o caminho do executável
```

---

## 3. Configuração

Toda a configuração está centralizada em `cath_analysis/config.py`. O único
parâmetro que você **precisa editar** é o diretório base:

```python
# cath_analysis/config.py
BASE_PATH = Path("/seu/caminho/para/cath")
```

Os demais caminhos são derivados automaticamente:

| Constante | Caminho padrão | Conteúdo |
|---|---|---|
| `STRUCTURES_PATH` | `BASE_PATH/structures/` | PDBs brutos baixados |
| `STRUCTURES_CLEAN_PATH` | `BASE_PATH/structures_clean/` | PDBs limpos (cadeia A) |
| `LOGS_PATH` | `BASE_PATH/logs/` | Relatórios e logs de erros |
| `ANALYSIS_PATH` | `BASE_PATH/analysis/` | Gráficos e CSVs gerados |

---

## 4. Como executar

### Fluxo básico

```bash
# Executa todas as 6 etapas com menus interativos
python main.py
```

### Modo exploração (já tem a base baixada)

Se você já baixou e limpou as estruturas, pule diretamente para a análise:

```bash
python main.py --explore
```

O pipeline detecta automaticamente dados existentes. Ao rodar com as etapas 1 ou 2,
você verá:

```
  ┌─ Dados existentes ───────────────────────────────────┐
  │  Índice CATH:       atualizado em 2025-11-08 13:30
  │  Estruturas raw:    8.421 PDBs  (último: 2025-11-08 14:12)
  │  Estruturas clean:  7.890 PDBs  (último: 2025-11-08 15:44)
  └──────────────────────────────────────────────────────┘

  Baixar novamente? [s/N]
```

Responda `N` (ou Enter) para pular e usar os dados existentes.

### Referência completa de flags

| Comando | Descrição |
|---|---|
| `python main.py` | Pipeline completo, modo interativo |
| `python main.py --explore` | Pula etapas 1-2, vai direto à análise |
| `python main.py --steps 3 4 5` | Executa apenas as etapas especificadas |
| `python main.py --force-download` | Re-baixa e relimpa sem perguntar |
| `python main.py --skip-download` | Pula etapas 1-2 sem perguntar |
| `python main.py --no-interactive` | Sem menus — executa tudo automaticamente |
| `python main.py --plots g3 g4` | Pré-seleciona grupos de gráficos no menu |
| `python main.py --plots helical_wheel_average pca_aa_composition` | Pré-seleciona gráficos individuais |
| `python main.py --list-plots` | Lista todos os 21 gráficos disponíveis e sai |

### Exemplos práticos

```bash
# Primeiro uso: baixa tudo e analisa
python main.py

# Sessão de re-análise (dados já existem)
python main.py --explore

# Gerar só gráficos evolutivos e de código genético
python main.py --explore --plots g3 g5 g6

# Pipeline completo em servidor (sem interação)
python main.py --no-interactive --force-download

# Só baixar e limpar
python main.py --steps 1 2

# Ver o que pode plotar
python main.py --list-plots
```

### Menu interativo de gráficos

Após as etapas de análise, o menu aparece automaticamente:

```
═══════════════════════════════════════════════════════════════════
  CATH ANALYSIS – SELETOR DE GRÁFICOS EVOLUTIVOS
═══════════════════════════════════════════════════════════════════

  1. Análise Básica de Hélices
  ───────────────────────────────────────────────────────
  [ ]  1. [DSSP]       Propensão para hélices (observada vs Chou-Fasman)
  [ ]  2. [DSSP]       Preferência por posição na hélice
  ...

  Comandos:
  <número>      Selecionar/deselecionar gráfico
  g<número>     Selecionar/deselecionar grupo inteiro (ex: g3)
  all           Selecionar todos
  none          Desmarcar todos
  info <num>    Ver descrição detalhada
  run           Gerar gráficos selecionados
  q             Sair sem plotar
```

---

## 5. Etapas do pipeline

### Etapa 1 — Download de estruturas

Baixa o índice de domínios do CATH e faz download paralelo das estruturas PDB
da classe Mainly-Alpha (classe 1) diretamente do RCSB.

- Fonte do índice: `http://download.cathdb.info/...`
- Fonte dos PDBs: `https://files.rcsb.org/download/{}.pdb`
- Download paralelo com até 20 workers simultâneos
- Retoma automaticamente (pula arquivos já baixados)
- Gera: `logs/download_report.txt`, `mainly_alpha_pdb_codes.txt`

**Saída esperada:** milhares de arquivos `.pdb` em `structures/`

### Etapa 2 — Limpeza de estruturas

Processa cada PDB para remover ruído e padronizar o formato:

- Mantém apenas o **modelo 0** (primeiro modelo em NMR)
- Mantém apenas a **cadeia A**
- Remove **heteroátomos** (HETATM): água, íons, ligantes, cofatores
- Remove resíduos não-padrão (mantém apenas os 20 aminoácidos canônicos)
- Estruturas sem cadeia A ou com 0 resíduos após limpeza são descartadas
- Processamento paralelo com `ThreadPoolExecutor`

**Saída esperada:** arquivos `.pdb` limpos em `structures_clean/`

### Etapa 3 — Frequência de aminoácidos

Conta a ocorrência de cada aminoácido em todas as estruturas limpas.

- Contagem global e por estrutura
- Gera: `amino_acid_frequencies_global.txt`, `amino_acid_frequencies_per_structure.csv`

### Etapa 4 — Análise avançada de hélices (DSSP)

Usa o DSSP para determinar a estrutura secundária de cada resíduo e analisa
a composição das hélices em profundidade:

- Propensão de cada AA para hélices (observada vs. Chou-Fasman)
- Preferência por posição: N-terminal, meio, C-terminal
- Distribuição de comprimentos de hélices
- Padrão heptad de hidrofobicidade

### Etapa 5 — Tipos de hélices (H/G/I)

Classifica cada hélice pelo tipo DSSP e compara composição:

| Tipo DSSP | Nome | Ligação de H | Resíduos/volta |
|---|---|---|---|
| H | Alpha helix | i → i+4 | 3,6 |
| G | 3-10 helix | i → i+3 | 3,0 |
| I | Pi helix | i → i+5 | 4,4 |

### Etapa 6 — Coleta de dados evolutivos

Passo DSSP único e abrangente que coleta dados para todas as análises evolutivas:

- Sequências individuais de cada hélice (para moment hidrofóbico e helical wheel)
- Preferências de N-cap e C-cap (primeiras/últimas 3 posições)
- Distribuição completa de AAs por posição no heptad
- Pares de transição entre tipos de hélice
- Fração helicoidal por estrutura

---

## 6. Catálogo de análises e gráficos

### Grupo 1 — Análise Básica de Hélices

#### 1.1 Propensão para hélices
**Arquivo:** `helix_propensities.png`

Dois subplots:
- **Barras duplas:** compara a propensão observada (calculada dos dados) com a
  propensão teórica da escala de Chou-Fasman para cada aminoácido.
- **Scatter:** frequência global do AA vs. frequência dentro de hélices, com linha
  diagonal indicando ausência de enriquecimento.

**Valor:** Revela quais AAs estão sobre- ou sub-representados em hélices.
ALA, GLU, LEU com propensão >1.2 são os grandes formadores; GLY e PRO (<0.6) são
disruptores. Desvios entre observado e Chou-Fasman indicam pressões específicas
do conjunto CATH.

---

#### 1.2 Preferência por posição na hélice
**Arquivo:** `helix_positions.png`

Heatmap 20×3 com a frequência de cada aminoácido no N-terminal, meio e C-terminal
da hélice.

**Valor:** Certas posições são estruturalmente mais restritivas. AAs com capacidade
de fazer N-cap (ASN, ASP) aparecem enriquecidos no N-terminal. GLY é tolerado no
C-terminal mas raramente no meio de hélices longas. Essa análise é base para
identificar restrições evolutivas posicionais.

---

#### 1.3 Distribuição de comprimentos
**Arquivo:** `helix_lengths.png`

Histograma + boxplot dos comprimentos de hélices em resíduos.

**Valor:** Hélices alpha estáveis tendem a ter 10–20 resíduos. A distribuição
de comprimentos reflete a pressão seletiva para manter a hélice coesa: hélices
muito curtas (4–6 res.) têm menor estabilidade térmica.

---

#### 1.4 Padrão heptad
**Arquivo:** `heptad_pattern.png`

Frequência de resíduos hidrofóbicos em cada posição do heptad repeat (a–g).

**Valor:** Em coiled-coils, as posições *a* e *d* do heptad são hidrofóbicas para
formar a interface entre as cadeias. Esse gráfico revela se o conjunto CATH
Mainly-Alpha exibe essa periodicidade, sugerindo presença de coiled-coils ou
de pressão seletiva para periodicidade hidrofóbica.

---

### Grupo 2 — Tipos de Hélices

#### 2.1 Distribuição de tipos H/G/I
**Arquivo:** `helix_type_distribution.png`

Bar chart e pie chart com a contagem e proporção de hélices Alpha (H), 3-10 (G)
e Pi (I).

**Valor:** Em proteínas globulares, >90% das hélices são do tipo H. 3-10 helices
aparecem tipicamente nos terminais de hélices maiores ou em loops curtos. Pi helices
são raras (<1%) e geralmente associadas a regiões funcionais. A proporção informa
o grau de diversidade estrutural do conjunto.

---

#### 2.2 Composição por tipo (heatmap)
**Arquivo:** `helix_type_composition_heatmap.png`

Heatmap 20×3 comparando a frequência de cada AA em H, G e I helices
simultaneamente.

**Valor:** Permite identificar AAs diferencialmente enriquecidos por tipo. PRO,
por exemplo, tende a ocorrer em 3-10 helices mais do que em alpha helices. Essas
diferenças refletem restrições geométricas distintas de cada tipo.

---

#### 2.3 Top-10 AAs por tipo
**Arquivo:** `helix_type_top_amino_acids.png`

Barras horizontais com os 10 AAs mais frequentes em cada tipo de hélice.

---

#### 2.4 Comprimento por tipo
**Arquivo:** `helix_length_by_type.png`

Boxplot + histograma sobrepostos comparando a distribuição de comprimentos entre H,
G e I.

**Valor:** 3-10 helices são estruturalmente limitadas a 3–6 resíduos pela geometria
da ligação de hidrogênio i→i+3. Alpha helices são estáveis de 5 a >30 resíduos.
O gráfico quantifica essa distribuição no conjunto real.

---

#### 2.5 Comparação estatística Alpha vs. 3-10
**Arquivo:** `helix_type_statistical_comparison.png`

Barras de diferença de frequência e scatter direto Alpha vs. 3-10 para cada AA.

**Valor:** Destaca os AAs com maior diferenciação entre os dois tipos. AAs
enriquecidos em 3-10 (como PRO, GLY) são candidatos a estudos de transição
evolutiva entre tipos de hélice.

---

### Grupo 3 — Restrições Estruturais

#### 3.1 Helical wheel composto
**Arquivo:** `helical_wheel_average.png`

Projeção polar: frequência de resíduos hidrofóbicos por posição angular na hélice.
Cada resíduo ocupa 100° em relação ao anterior, completando ~3,6 resíduos por volta.
A escala de cores vai de azul (polar) a vermelho (hidrofóbico).

**Valor:** Revela a **anfipatiticidade** do conjunto. Em hélices anfipáticas
funcionais (e.g., hélices de membrana, hélices de ligação a lipídeos), um setor
angular claro com alta frequência hidrofóbica deve emergir. A ausência desse padrão
indica que o conjunto é dominado por hélices de interior proteico com distribuição
uniforme. É uma das análises mais diretas de restrição estrutural evolutiva.

---

#### 3.2 Momento hidrofóbico de Eisenberg
**Arquivo:** `hydrophobic_moment_distribution.png`

Distribuição do momento hidrofóbico (μH) por hélice, calculado pela escala de
Eisenberg. Hélices são classificadas em anfipáticas baixas (<0,2), médias (0,2–0,4)
e altas (≥0,4).

**Valor:** O momento hidrofóbico é um preditor de função. Hélices com μH alto são
associadas a: inserção em membranas, interação com lipídeos, atividade antimicrobiana.
A distribuição do μH no conjunto CATH informa a proporção de hélices funcionalmente
anfipáticas vs. estruturalmente internas.

---

#### 3.3 N-cap e C-cap (3 posições)
**Arquivo:** `ncap_ccap_preferences.png`

Heatmaps das frequências de AAs nas 3 primeiras posições da hélice (N1, N2, N3)
e nas 3 últimas (C3, C2, C1).

**Valor:** As posições de *cap* têm preferências fortíssimas sob seleção purificadora:
- **N-cap:** ASN e ASP são enriquecidos porque fazem ligações de hidrogênio com
  os NH da hélice não satisfeitos pelas posições iniciais.
- **C-cap:** GLY é favorecido pela flexibilidade que permite a curvatura da cadeia
  ao sair da hélice.

Desvios dessas preferências podem indicar hélices com funções especiais ou pressões
evolutivas distintas (e.g., hélices de proteínas termófilas).

---

### Grupo 4 — Composição e Variabilidade

#### 4.1 PCA da composição por estrutura
**Arquivo:** `pca_aa_composition.png`

Cada proteína é representada como um vetor de 20 frequências (uma por AA), reduzido
para 2 dimensões por PCA. Pontos são coloridos pelo AA dominante.

**Valor:** Agrupa proteínas por "assinatura" de composição. Estruturas funcionalmente
relacionadas devem agrupar. Outliers revelam proteínas com composição incomum.
PC1 tipicamente captura a variação hidrofóbico/hidrofílico; PC2 captura outros
gradientes físico-químicos. Requer `scikit-learn`.

---

#### 4.2 Distribuição de conteúdo helicoidal
**Arquivo:** `helix_content_distribution.png`

Histograma + violin plot da fração de resíduos em hélice por proteína.

**Valor:** Embora todas as proteínas sejam da classe Mainly-Alpha, a distribuição
de conteúdo helicoidal varia significativamente. Proteínas com >80% são casos
extremos (e.g., hélices de membrana). A distribuição real quantifica a heterogeneidade
da classe e pode ser comparada com outras classes CATH como referência evolutiva.

---

#### 4.3 Co-ocorrência de AAs em hélices
**Arquivo:** `aa_cooccurrence.png`

Heatmap 20×20 onde cada célula representa a frequência relativa de pares de AAs
que aparecem juntos na mesma hélice.

**Valor:** Pares com alta co-ocorrência podem refletir:
- **Interações físicas:** GLU–LYS e ASP–ARG formam salt bridges (i→i+4) que
  estabilizam hélices.
- **Perfis evolutivos compartilhados:** AAs com propensão similar tendem a ser
  trocados por mutação e co-ocorrem mais.
- **Viés composicional:** LEU e ALA co-ocorrem simplesmente por serem os mais
  frequentes.

A análise deve ser interpretada com baseline de frequência esperada (produto das
frequências marginais).

---

#### 4.4 Comprimento vs. composição
**Arquivo:** `helix_length_vs_composition.png`

Barras agrupadas: frequência de cada grupo físico-químico (hidrofóbico, polar,
carregado+, carregado–, especial) em hélices curtas (4–9 res.), médias (10–19)
e longas (≥20).

**Valor:** Hélices longas tendem a ser mais ricas em AAs de alta propensão (ALA,
LEU, GLU). Hélices curtas toleram mais GLY e PRO. Essa relação comprimento-composição
é uma evidência direta de seleção estrutural: AAs disruptores são tolerados apenas
em contextos curtos onde a hélice não precisa ser mantida por muitos resíduos.

---

### Grupo 5 — Transições e História Evolutiva

#### 5.1 Matriz de transição H→G→I
**Arquivo:** `helix_transition_matrix.png`

Heatmap normalizado com a probabilidade de cada tipo de hélice ser seguido por
outro tipo dentro da mesma proteína.

**Valor:** Revela se existe uma hierarquia de transições. A hipótese evolutiva é
que 3-10 helices são estados ancestrais ou transicionais de alpha helices — proteínas
que "experimentam" hélices mais curtas antes de estabilizar em alpha. Uma alta taxa
de H→G seguida de G→H indicaria coexistência de ambos os tipos no mesmo contexto
estrutural.

---

#### 5.2 Razão 3-10/(H+G+I) por comprimento
**Arquivo:** `g_ratio_by_length.png`

Linha com a fração de hélices 3-10 em cada comprimento. Barra cinza de fundo mostra
o total de hélices por comprimento.

**Valor:** Testa diretamente a hipótese de que hélices curtas são predominantemente
3-10. Se a curva cai monotonicamente com o comprimento, confirma a restrição
geométrica: hélices longas simplesmente não conseguem manter a geometria 3-10.
Picos em comprimentos específicos podem indicar contextos estruturais especiais.

---

#### 5.3 Entropia de Shannon por posição do heptad
**Arquivo:** `shannon_entropy_heptad.png`

Entropia de Shannon calculada a partir da distribuição de AAs em cada posição (a–g)
do heptad. Baixa entropia = maior conservação = maior pressão seletiva.

**Valor:** Em coiled-coils clássicos, posições *a* e *d* têm baixa entropia
(dominadas por LEU, ILE, VAL) enquanto posições *b*, *c*, *e*, *f*, *g* têm alta
entropia (qualquer AA). A assinatura heptad de entropia é um indicador da proporção
de coiled-coils no conjunto e da magnitude da seleção nessas posições.

---

### Grupo 6 — Viés do Código Genético

#### 6.1 Degenerescência de códons vs. propensão
**Arquivo:** `codon_degeneracy_vs_propensity.png`

Dois scatters lado a lado:
- Número de códons do AA (degenerescência) vs. propensão teórica (Chou-Fasman)
- Número de códons vs. propensão observada no conjunto CATH

**Valor:** Permite separar dois mecanismos distintos de composição:
- **Seleção estrutural:** AAs com alta propensão são favorecidos pela geometria da hélice.
- **Deriva neutra / viés mutacional:** AAs com mais códons (LEU=6, SER=6) ocorrem
  mais por pura probabilidade mutacional.

Se a correlação propensão×códons é fraca, a composição das hélices é dominada por
seleção. Se é forte, o viés do código genético tem papel relevante — uma questão
aberta em evolução molecular.

---

#### 6.2 Enriquecimento vs. proteoma humano
**Arquivo:** `proteome_comparison.png`

Dois subplots:
- **Barras:** razão Observado/Proteoma (>1 = enriquecido em hélices, <1 = depleto)
- **Scatter:** frequência no proteoma humano de referência vs. frequência observada

Frequências de referência: proteoma humano médio (UniProt/Swiss-Prot).

**Valor:** Normaliza a composição das hélices pelo background proteômico. Um AA
pode ser frequente nas hélices simplesmente porque é frequente em proteínas em geral
(caso do LEU). O enriquecimento real — controlando pelo background — revela a
seleção verdadeira. ALA com enriquecimento >1.3 e PRO com <0.5 são esperados;
desvios desses valores indicam particularidades da classe Mainly-Alpha.

---

## 7. Estrutura de arquivos gerados

```
BASE_PATH/
├── cath-domain-list.txt              # Índice CATH baixado
├── mainly_alpha_pdb_codes.txt        # Códigos PDB da classe 1
│
├── structures/                       # Etapa 1: PDBs brutos
│   ├── 1abc.pdb
│   └── ...
│
├── structures_clean/                 # Etapa 2: PDBs limpos
│   ├── 1abc.pdb                      # Apenas cadeia A, sem heteroátomos
│   └── ...
│
├── logs/
│   ├── download_report.txt
│   ├── failed_downloads.txt
│   ├── cleaning_report.txt
│   ├── no_chain_a.txt
│   └── cleaning_errors.txt
│
└── analysis/
    ├── amino_acid_frequencies_global.txt
    ├── amino_acid_frequencies_per_structure.csv
    │
    ├── helix_analysis_comprehensive_report.txt
    ├── helix_propensities.csv
    ├── helix_positions.csv
    ├── aa_properties_distribution.csv
    │
    ├── helix_types_comprehensive_report.txt
    ├── helix_type_composition.csv
    ├── helix_type_top_residues.csv
    ├── helix_type_statistical_comparison.csv
    │
    └── *.png                         # 21 gráficos gerados
```

---

## 8. Estrutura do código

```
cath/
├── main.py                           # Ponto de entrada (CLI)
├── requirements.txt
│
└── cath_analysis/
    ├── __init__.py
    ├── config.py                     # Constantes e caminhos
    ├── downloader.py                 # Etapa 1: download CATH + PDB
    ├── cleaner.py                    # Etapa 2: limpeza PDB
    ├── frequency_analysis.py         # Etapa 3: frequência de AAs
    ├── helix_analysis.py             # Etapa 4: propensão, posições, heptad
    ├── helix_types.py                # Etapa 5: classificação H/G/I
    ├── evolutionary_analysis.py      # Etapa 6: dados evolutivos (DSSP único)
    ├── plotting.py                   # Todas as 21 funções de visualização
    └── menu.py                       # Menu interativo de seleção de gráficos
```

### Dependências entre módulos

```
config.py ←── todos os módulos
downloader.py ──→ Etapa 1
cleaner.py ──────→ Etapa 2
frequency_analysis.py ──→ Etapa 3 → PCA (plotting)
helix_analysis.py ──────→ Etapa 4 → plotting
helix_types.py ─────────→ Etapa 5 → plotting
evolutionary_analysis.py → Etapa 6 → plotting
menu.py ──────────────────────────────→ main.py
```

---

## 9. Dependências externas

| Pacote | Versão mínima | Uso |
|---|---|---|
| `biopython` | 1.81 | Parser PDB, DSSP, seleção de cadeias |
| `numpy` | 1.24 | Cálculos numéricos |
| `pandas` | 2.0 | DataFrames de resultados |
| `matplotlib` | 3.7 | Visualizações |
| `seaborn` | 0.12 | Heatmaps e estilos |
| `requests` | 2.28 | Download de PDBs |
| `tqdm` | 4.65 | Barras de progresso |
| `scipy` | 1.11 | (reservado para extensões) |
| `scikit-learn` | — | Opcional: PCA (etapa 4.1) |
| `mkdssp` | — | Externo: estrutura secundária |

---

## 10. Perguntas frequentes

**Q: O download trava em alguns PDBs. É normal?**

Sim. Alguns códigos PDB foram obsoletados ou renomeados no RCSB. O pipeline marca
esses casos como `error` e continua. Os códigos que falharam ficam em
`logs/failed_downloads.txt`.

---

**Q: DSSP falha em várias estruturas. Por quê?**

DSSP requer que os átomos backbone (N, Cα, C, O) estejam presentes e com distâncias
razoáveis. Estruturas com resolução muito baixa, modelos incompletos ou erros de
modelagem podem falhar. O pipeline registra e pula essas estruturas.

---

**Q: Posso usar estruturas de outra fonte (não CATH)?**

Sim. Basta colocar os PDBs limpos (cadeia A, sem heteroátomos) em
`structures_clean/` e executar a partir da etapa 3:

```bash
python main.py --steps 3 4 5 6
```

---

**Q: O PCA não está disponível.**

Instale o `scikit-learn`:
```bash
pip install scikit-learn
```

---

**Q: Como citar o CATH?**

> Sillitoe I. et al. (2021). CATH: expanding the horizons of structure-based
> functional annotations for genome sequences. *Nucleic Acids Research*, 49(D1),
> D377–D385. https://doi.org/10.1093/nar/gkaa1079
