"""Filtragem S40: cria subconjunto nao-redundante com symlinks."""

import urllib.request
from pathlib import Path

CATH_S40_URL = (
    "http://download.cathdb.info/cath/releases/latest-release"
    "/non-redundant-data-sets/cath-dataset-nonredundant-S40.list"
)


def download_s40_list(base_path: Path) -> list[str]:
    """Baixa a lista S40 e retorna domain IDs."""
    s40_file = base_path / "cath-s40-domains.txt"
    print(f"  Baixando lista S40 de {CATH_S40_URL} ...")
    urllib.request.urlretrieve(CATH_S40_URL, s40_file)
    domain_ids = []
    for line in s40_file.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            domain_ids.append(line.split()[0].lower())
    print(f"  Lista S40: {len(domain_ids):,} dominios")
    return domain_ids


def create_s40_directory(clean_dir: Path, s40_dir: Path, domain_ids: list[str]) -> dict:
    """Cria s40_dir com symlinks para os PDBs representantes S40."""
    s40_dir.mkdir(parents=True, exist_ok=True)

    s40_set = set(domain_ids)
    found = 0
    not_found = 0
    skipped_mac = 0

    for pdb_file in clean_dir.glob("*.pdb"):
        if pdb_file.name.startswith("._"):
            skipped_mac += 1
            continue
        domain_id = pdb_file.stem.lower()
        if domain_id in s40_set:
            target = s40_dir / pdb_file.name
            if not target.exists():
                target.symlink_to(pdb_file.resolve())
            found += 1
        else:
            not_found += 1

    return {
        "total_s40": len(domain_ids),
        "found": found,
        "not_found_in_s40": not_found,
        "skipped_mac": skipped_mac,
    }


def run_s40_filter(base_path: Path, clean_dir: Path) -> Path:
    """Entry point: baixa S40, cria diretorio filtrado. Retorna caminho do s40_dir."""
    s40_dir = base_path / "structures_clean_s40"

    # Se ja existe com arquivos, nao refaz
    existing = sum(1 for _ in s40_dir.glob("*.pdb")) if s40_dir.exists() else 0
    if existing > 0:
        print(f"  Diretorio S40 ja existe: {existing:,} estruturas em {s40_dir}")
        return s40_dir

    domain_ids = download_s40_list(base_path)
    stats = create_s40_directory(clean_dir, s40_dir, domain_ids)

    print(f"\n  Filtro S40 concluido:")
    print(f"    Dominios na lista S40:   {stats['total_s40']:,}")
    print(f"    Encontrados no clean/:   {stats['found']:,}")
    print(f"    Nao encontrados (S40 somente): {stats['not_found_in_s40']:,}")
    print(f"    Symlinks criados em:     {s40_dir}")

    return s40_dir
