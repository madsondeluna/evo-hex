"""
Configurações centralizadas do pipeline CATH.

Altere BASE_PATH para o diretório onde os dados serão armazenados.
"""
import multiprocessing
from pathlib import Path

# ── Caminhos ──────────────────────────────────────────────────────────────────
BASE_PATH = Path("/Volumes/promethion/cath")

STRUCTURES_PATH = BASE_PATH / "structures"
STRUCTURES_CLEAN_PATH = BASE_PATH / "structures_clean"
LOGS_PATH = BASE_PATH / "logs"
ANALYSIS_PATH = BASE_PATH / "analysis"


def glob_pdb(directory: Path) -> list[Path]:
    """Retorna PDBs válidos de um diretório, excluindo resource forks do macOS (._*)."""
    return [p for p in directory.glob("*.pdb") if not p.name.startswith("._")]


# ── URLs ──────────────────────────────────────────────────────────────────────
CATH_DOMAIN_LIST_URL = (
    "http://download.cathdb.info/cath/releases/latest-release"
    "/cath-classification-data/cath-domain-list.txt"
)
PDB_DOWNLOAD_URL = "https://files.rcsb.org/download/{}.pdb"
CATH_S40_URL = (
    "http://download.cathdb.info/cath/releases/latest-release"
    "/non-redundant-data-sets/cath-dataset-nonredundant-S40.list"
)

# ── Parâmetros de download ────────────────────────────────────────────────────
MAINLY_ALPHA_CLASS = "1"   # Class 1 no CATH = Mainly Alpha
CHUNK_SIZE = 1024 * 1024   # 1 MB para streaming
REQUEST_TIMEOUT = 30        # segundos

# ── Paralelismo ───────────────────────────────────────────────────────────────
# Downloads: mais workers (I/O bound)
DOWNLOAD_WORKERS = min(20, multiprocessing.cpu_count() * 2)
# Processamento local: deixa 1 CPU livre para o SO
PROCESS_WORKERS = max(1, multiprocessing.cpu_count() - 1)

# ── Aminoácidos padrão ────────────────────────────────────────────────────────
STANDARD_AMINO_ACIDS = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
})

# ── Propriedades físico-químicas ──────────────────────────────────────────────
AA_PROPERTIES: dict[str, list[str]] = {
    "Hidrofóbico":         ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO"],
    "Polar":               ["SER", "THR", "CYS", "TYR", "ASN", "GLN"],
    "Carregado positivo":  ["LYS", "ARG", "HIS"],
    "Carregado negativo":  ["ASP", "GLU"],
    "Especial":            ["GLY", "PRO"],
}

# Subconjunto estritamente hidrofóbico (usado na análise de heptad)
HYDROPHOBIC_AA = frozenset({"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP"})

# ── Mapeamentos de aminoácidos ─────────────────────────────────────────────────
THREE_TO_ONE: dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

ONE_TO_THREE: dict[str, str] = {v: k for k, v in THREE_TO_ONE.items()}

# ── Escala de Eisenberg ────────────────────────────────────────────────────────
EISENBERG_SCALE: dict[str, float] = {
    "ALA":  0.62, "ARG": -2.53, "ASN": -0.78, "ASP": -0.90,
    "CYS":  0.29, "GLN": -0.85, "GLU": -0.74, "GLY":  0.48,
    "HIS": -0.40, "ILE":  1.38, "LEU":  1.06, "LYS": -1.50,
    "MET":  0.64, "PHE":  1.19, "PRO":  0.12, "SER": -0.18,
    "THR": -0.05, "TRP":  0.81, "TYR":  0.26, "VAL":  1.08,
}

# ── Degenerescência de códons ─────────────────────────────────────────────────
CODON_DEGENERACY: dict[str, int] = {
    "ALA": 4, "ARG": 6, "ASN": 2, "ASP": 2, "CYS": 2,
    "GLN": 2, "GLU": 2, "GLY": 4, "HIS": 2, "ILE": 3,
    "LEU": 6, "LYS": 2, "MET": 1, "PHE": 2, "PRO": 4,
    "SER": 6, "THR": 4, "TRP": 1, "TYR": 2, "VAL": 4,
}

# ── Frequência no proteoma humano (%) ─────────────────────────────────────────
PROTEOME_FREQ: dict[str, float] = {
    "ALA": 6.97, "ARG": 5.53, "ASN": 4.06, "ASP": 5.25, "CYS": 2.27,
    "GLN": 3.93, "GLU": 6.75, "GLY": 6.87, "HIS": 2.29, "ILE": 5.49,
    "LEU": 9.68, "LYS": 5.19, "MET": 2.32, "PHE": 3.87, "PRO": 5.02,
    "SER": 7.14, "THR": 5.57, "TRP": 1.33, "TYR": 3.21, "VAL": 6.47,
}

# ── Escala de Chou-Fasman ─────────────────────────────────────────────────────
HELIX_PROPENSITY: dict[str, float] = {
    "ALA": 1.42, "GLU": 1.51, "LEU": 1.21, "MET": 1.45,
    "GLN": 1.11, "LYS": 1.16, "ARG": 0.98, "HIS": 1.00,
    "VAL": 1.06, "ILE": 1.08, "TYR": 0.69, "PHE": 1.13,
    "TRP": 1.08, "THR": 0.83, "SER": 0.77, "CYS": 0.70,
    "ASP": 1.01, "ASN": 0.67, "GLY": 0.57, "PRO": 0.57,
}

# ── Tipos de hélices (notação DSSP) ───────────────────────────────────────────
HELIX_TYPES: dict[str, str] = {
    "H": "Alpha helix",
    "G": "3-10 helix",
    "I": "Pi helix",
}

HELIX_CHARACTERISTICS: dict[str, dict] = {
    "H": {
        "name": "Alpha helix",
        "residues_per_turn": 3.6,
        "rise_per_residue": 1.5,
        "hydrogen_bond": "i to i+4",
        "typical_length": "10-20 residues",
    },
    "G": {
        "name": "3-10 helix",
        "residues_per_turn": 3.0,
        "rise_per_residue": 2.0,
        "hydrogen_bond": "i to i+3",
        "typical_length": "3-5 residues",
    },
    "I": {
        "name": "Pi helix",
        "residues_per_turn": 4.4,
        "rise_per_residue": 1.15,
        "hydrogen_bond": "i to i+5",
        "typical_length": "7-10 residues",
    },
}
