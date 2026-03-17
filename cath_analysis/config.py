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
STRUCTURES_CLEAN_S40_PATH = BASE_PATH / "structures_clean_s40"
LOGS_PATH = BASE_PATH / "logs"
ANALYSIS_PATH = BASE_PATH / "analysis"


def glob_pdb(directory: Path) -> list[Path]:
    """Retorna PDBs válidos de um diretório, excluindo resource forks do macOS (._*)."""
    return [p for p in directory.glob("*.pdb") if not p.name.startswith("._") and p.exists()]


# ── URLs ──────────────────────────────────────────────────────────────────────
CATH_DOMAIN_LIST_URL = (
    "http://download.cathdb.info/cath/releases/latest-release"
    "/cath-classification-data/cath-domain-list.txt"
)
PDB_DOWNLOAD_URL = "https://files.rcsb.org/download/{}.pdb"

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
