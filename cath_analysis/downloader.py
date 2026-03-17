"""
Download de estruturas PDB mainly-alpha do banco CATH.

Fluxo principal:
    1. setup_directories()          – garante hierarquia de pastas
    2. download_cath_domain_list()  – baixa o índice CATH
    3. parse_cath_domains()         – extrai códigos PDB da classe 1
    4. download_structures_parallel() – baixa PDBs em paralelo
    5. save_download_results()      – persiste listas e relatório
"""

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests
from tqdm.auto import tqdm

from .config import (
    BASE_PATH,
    CATH_DOMAIN_LIST_URL,
    CHUNK_SIZE,
    DOWNLOAD_WORKERS,
    LOGS_PATH,
    MAINLY_ALPHA_CLASS,
    PDB_DOWNLOAD_URL,
    REQUEST_TIMEOUT,
    STRUCTURES_PATH,
    glob_pdb,
)

logger = logging.getLogger(__name__)


def check_existing_data() -> dict:
    """Inspeciona o que já existe no BASE_PATH e retorna um resumo.

    Returns:
        Dict com chaves:
          - has_cath_list   : bool
          - cath_list_date  : str (data de modificação formatada) ou ""
          - raw_count       : int (nº de PDBs em structures/)
          - raw_date        : str (data do PDB mais recente) ou ""
          - clean_count     : int (nº de PDBs em structures_clean/)
          - clean_date      : str (data do PDB mais recente) ou ""
    """
    from datetime import datetime

    def _fmt(path: Path) -> str:
        return datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y-%m-%d %H:%M")

    cath_file = BASE_PATH / "cath-domain-list.txt"
    raw_dir = STRUCTURES_PATH
    clean_dir = BASE_PATH / "structures_clean"

    # Usa mtime do diretório para evitar varrer todos os arquivos (lento em drives externos)
    raw_count = sum(1 for _ in raw_dir.glob("*.pdb") if not _.name.startswith("._")) if raw_dir.exists() else 0
    clean_count = sum(1 for _ in clean_dir.glob("*.pdb") if not _.name.startswith("._")) if clean_dir.exists() else 0
    raw_date = _fmt(raw_dir) if raw_dir.exists() and raw_count else ""
    clean_date = _fmt(clean_dir) if clean_dir.exists() and clean_count else ""

    return {
        "has_cath_list": cath_file.exists(),
        "cath_list_date": _fmt(cath_file) if cath_file.exists() else "",
        "raw_count": raw_count,
        "raw_date": raw_date,
        "clean_count": clean_count,
        "clean_date": clean_date,
    }


def setup_directories() -> tuple[Path, Path, Path]:
    """Cria e valida a hierarquia de diretórios necessária.

    Returns:
        Tupla (base_path, structures_path, logs_path).

    Raises:
        RuntimeError: Se algum diretório não puder ser criado.
    """
    for directory in (BASE_PATH, STRUCTURES_PATH, LOGS_PATH):
        directory.mkdir(parents=True, exist_ok=True)

    missing = [d for d in (BASE_PATH, STRUCTURES_PATH, LOGS_PATH) if not d.exists()]
    if missing:
        raise RuntimeError(f"Falha ao criar diretórios: {missing}")

    logger.info("Diretórios verificados: %s, %s, %s", BASE_PATH, STRUCTURES_PATH, LOGS_PATH)
    return BASE_PATH, STRUCTURES_PATH, LOGS_PATH


def download_cath_domain_list(base_path: Path) -> Path:
    """Baixa o arquivo de índice de domínios CATH (skip se já existir).

    Args:
        base_path: Diretório base onde o arquivo será salvo.

    Returns:
        Caminho para o arquivo de lista CATH.
    """
    cath_file = base_path / "cath-domain-list.txt"

    if cath_file.exists():
        logger.info("Arquivo CATH já existe: %s", cath_file)
        return cath_file

    logger.info("Baixando lista de domínios CATH...")
    response = requests.get(CATH_DOMAIN_LIST_URL, timeout=60)
    response.raise_for_status()

    cath_file.write_text(response.text, encoding="utf-8")
    logger.info("Arquivo CATH salvo: %s", cath_file)
    return cath_file


def parse_cath_domains(cath_file: Path) -> set[str]:
    """Extrai códigos PDB de 4 letras para estruturas da classe mainly-alpha.

    Args:
        cath_file: Caminho para cath-domain-list.txt.

    Returns:
        Conjunto de códigos PDB em minúsculas.
    """
    pdb_codes: set[str] = set()

    with cath_file.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            domain_id, cath_class = parts[0], parts[1]
            if cath_class == MAINLY_ALPHA_CLASS:
                pdb_codes.add(domain_id[:4].lower())

    logger.info("Total de estruturas mainly-alpha: %d", len(pdb_codes))
    return pdb_codes


def _download_single_pdb(pdb_code: str, output_dir: Path) -> dict:
    """Baixa uma única estrutura PDB.

    Args:
        pdb_code: Código PDB de 4 letras.
        output_dir: Diretório de destino.

    Returns:
        Dict com chaves 'pdb_code', 'status' e 'message'.
    """
    output_file = output_dir / f"{pdb_code}.pdb"

    if output_file.exists():
        return {"pdb_code": pdb_code, "status": "exists", "message": "Arquivo já existe"}

    try:
        url = PDB_DOWNLOAD_URL.format(pdb_code.upper())
        response = requests.get(url, timeout=REQUEST_TIMEOUT, stream=True)
        response.raise_for_status()

        with output_file.open("wb") as f:
            for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
                if chunk:
                    f.write(chunk)

        return {"pdb_code": pdb_code, "status": "success", "message": "Download concluído"}

    except requests.exceptions.HTTPError as exc:
        return {
            "pdb_code": pdb_code,
            "status": "error",
            "message": f"HTTP {exc.response.status_code}",
        }
    except Exception as exc:  # noqa: BLE001
        return {"pdb_code": pdb_code, "status": "error", "message": str(exc)}


def download_structures_parallel(
    pdb_codes: set[str],
    output_dir: Path,
    max_workers: int = DOWNLOAD_WORKERS,
) -> dict:
    """Baixa estruturas PDB em paralelo com barra de progresso.

    Args:
        pdb_codes: Conjunto de códigos PDB.
        output_dir: Diretório de destino.
        max_workers: Número de threads paralelas.

    Returns:
        Dict com contadores 'total', 'success', 'exists', 'error' e 'failed_codes'.
    """
    pdb_list = sorted(pdb_codes)
    results: dict = {
        "total": len(pdb_list),
        "success": 0,
        "exists": 0,
        "error": 0,
        "failed_codes": [],
    }

    print(f"\nIniciando download de {results['total']} estruturas...")
    print(f"Workers paralelos: {max_workers}")
    print(f"Salvando em: {output_dir}\n")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(_download_single_pdb, code, output_dir): code
            for code in pdb_list
        }
        with tqdm(total=results["total"], desc="Downloading PDB structures") as pbar:
            for future in as_completed(futures):
                result = future.result()
                status = result["status"]
                if status == "success":
                    results["success"] += 1
                elif status == "exists":
                    results["exists"] += 1
                else:
                    results["error"] += 1
                    results["failed_codes"].append(result["pdb_code"])

                pbar.update(1)
                pbar.set_postfix(
                    Success=results["success"],
                    Exists=results["exists"],
                    Error=results["error"],
                )

    return results


def save_download_results(base_path: Path, pdb_codes: set[str], results: dict) -> None:
    """Persiste listas de códigos e relatório de download.

    Args:
        base_path: Diretório base do projeto.
        pdb_codes: Conjunto completo de códigos mainly-alpha.
        results: Dict retornado por download_structures_parallel.
    """
    logs_path = base_path / "logs"
    logs_path.mkdir(parents=True, exist_ok=True)

    # Lista completa
    codes_file = base_path / "mainly_alpha_pdb_codes.txt"
    codes_file.write_text("\n".join(sorted(pdb_codes)) + "\n")
    print(f"Lista de códigos PDB salva em: {codes_file}")

    # Códigos que falharam
    if results["failed_codes"]:
        failed_file = logs_path / "failed_downloads.txt"
        failed_file.write_text("\n".join(sorted(results["failed_codes"])) + "\n")
        print(f"Códigos com falha salvos em: {failed_file}")

    # Relatório
    total = results["total"]
    success_rate = (results["success"] + results["exists"]) / total * 100 if total else 0.0
    report = (
        f"{'='*60}\n"
        f"RELATÓRIO DE DOWNLOAD - EVO-HEX\n"
        f"{'='*60}\n\n"
        f"Total de estruturas:      {total}\n"
        f"Downloads bem-sucedidos:  {results['success']}\n"
        f"Já existentes:            {results['exists']}\n"
        f"Falhas:                   {results['error']}\n"
        f"Taxa de sucesso:          {success_rate:.1f}%\n\n"
        f"{'='*60}\n"
    )
    report_file = logs_path / "download_report.txt"
    report_file.write_text(report)
    print(f"Relatório salvo em: {report_file}")


def print_download_report(results: dict, structures_path: Path) -> None:
    """Exibe relatório final de download no console.

    Args:
        results: Dict retornado por download_structures_parallel.
        structures_path: Caminho onde as estruturas foram salvas.
    """
    total = results["total"]
    success_rate = (results["success"] + results["exists"]) / total * 100 if total else 0.0

    print(f"\n{'='*60}")
    print("RELATÓRIO FINAL DE DOWNLOAD")
    print(f"{'='*60}")
    print(f"Total de estruturas:           {total}")
    print(f"Downloads bem-sucedidos:       {results['success']}")
    print(f"Já existentes (não baixados):  {results['exists']}")
    print(f"Falhas:                        {results['error']}")
    print(f"Taxa de sucesso:               {success_rate:.1f}%")
    print(f"\nEstruturas salvas em: {structures_path}")
    print(f"{'='*60}")
