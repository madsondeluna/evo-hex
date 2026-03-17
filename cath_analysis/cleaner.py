"""
Limpeza de estruturas PDB.

Remove heteroátomos, água, íons e mantém apenas a cadeia A do modelo 0.
Processamento paralelo via ThreadPoolExecutor (tarefa I/O-bound).
"""

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

from Bio import PDB
from Bio.PDB import PDBIO, Select
from tqdm.auto import tqdm

from .config import PROCESS_WORKERS, STANDARD_AMINO_ACIDS

logger = logging.getLogger(__name__)


class ChainAOnlySelect(Select):
    """Seletor BioPython: mantém apenas cadeia A, modelo 0 e resíduos padrão."""

    def accept_model(self, model) -> bool:
        return model.id == 0

    def accept_chain(self, chain) -> bool:
        return chain.id == "A"

    def accept_residue(self, residue) -> bool:
        # hetfield != ' ' → HETATM (água, íon, ligante)
        if residue.id[0] != " ":
            return False
        return residue.resname.strip() in STANDARD_AMINO_ACIDS

    def accept_atom(self, atom) -> bool:  # noqa: ARG002
        return True


def clean_single_structure(input_file: Path, output_file: Path) -> dict:
    """Limpa uma estrutura PDB individual.

    Args:
        input_file: Arquivo PDB bruto.
        output_file: Destino do arquivo limpo.

    Returns:
        Dict com 'pdb_code', 'status' e 'residues'.
        Status possíveis: 'exists', 'no_chain_a', 'empty', 'success', 'error'.
    """
    pdb_code = input_file.stem

    if output_file.exists():
        return {"pdb_code": pdb_code, "status": "exists", "residues": 0}

    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_code, str(input_file))

        has_chain_a = any("A" in {c.id for c in m} for m in structure)
        if not has_chain_a:
            return {"pdb_code": pdb_code, "status": "no_chain_a", "residues": 0}

        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_file), ChainAOnlySelect())

        # Conta resíduos padrão após limpeza
        clean_structure = parser.get_structure(pdb_code, str(output_file))
        residue_count = sum(
            1
            for m in clean_structure
            for c in m
            for r in c
            if r.id[0] == " "
        )

        if residue_count == 0:
            output_file.unlink(missing_ok=True)
            return {"pdb_code": pdb_code, "status": "empty", "residues": 0}

        return {"pdb_code": pdb_code, "status": "success", "residues": residue_count}

    except Exception as exc:  # noqa: BLE001
        logger.debug("Erro ao limpar %s: %s", pdb_code, exc)
        return {"pdb_code": pdb_code, "status": "error", "residues": 0, "message": str(exc)}


def clean_all_structures(
    raw_dir: Path,
    clean_dir: Path,
    max_workers: int = PROCESS_WORKERS,
) -> Optional[dict]:
    """Processa todas as estruturas PDB em paralelo.

    Args:
        raw_dir: Diretório com PDBs brutos.
        clean_dir: Diretório de saída.
        max_workers: Threads para processamento paralelo.

    Returns:
        Dict de resultados, ou None se não houver arquivos.
    """
    pdb_files = list(raw_dir.glob("*.pdb"))
    if not pdb_files:
        logger.error("Nenhum arquivo .pdb encontrado em %s", raw_dir)
        return None

    total = len(pdb_files)
    existing_clean = len(list(clean_dir.glob("*.pdb")))

    print(f"\n{'='*60}")
    print("LIMPEZA DE ESTRUTURAS PDB")
    print(f"{'='*60}")
    print(f"Total de estruturas brutas: {total:,}")
    print(f"Já processadas:             {existing_clean:,}")
    print(f"Restantes:                  {total - existing_clean:,}")
    print(f"Workers:                    {max_workers}")

    results: dict = {
        "total": total,
        "success": 0,
        "exists": 0,
        "no_chain_a": 0,
        "empty": 0,
        "error": 0,
        "total_residues": 0,
        "no_chain_a_codes": [],
        "error_codes": [],
    }

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(clean_single_structure, f, clean_dir / f.name): f
            for f in pdb_files
        }
        with tqdm(total=total, desc="Cleaning structures") as pbar:
            for future in as_completed(futures):
                result = future.result()
                status = result["status"]

                if status == "success":
                    results["success"] += 1
                    results["total_residues"] += result["residues"]
                elif status == "exists":
                    results["exists"] += 1
                elif status == "no_chain_a":
                    results["no_chain_a"] += 1
                    results["no_chain_a_codes"].append(result["pdb_code"])
                elif status == "empty":
                    results["empty"] += 1
                elif status == "error":
                    results["error"] += 1
                    results["error_codes"].append(result["pdb_code"])

                pbar.update(1)
                pbar.set_postfix(
                    OK=results["success"],
                    Exist=results["exists"],
                    NoChainA=results["no_chain_a"],
                    Error=results["error"],
                )

    return results


def save_cleaning_results(logs_path: Path, results: dict) -> None:
    """Persiste logs e relatório da etapa de limpeza.

    Args:
        logs_path: Diretório de logs.
        results: Dict retornado por clean_all_structures.
    """
    logs_path.mkdir(parents=True, exist_ok=True)

    if results["no_chain_a_codes"]:
        no_chain_file = logs_path / "no_chain_a.txt"
        no_chain_file.write_text("\n".join(sorted(results["no_chain_a_codes"])) + "\n")
        print(f"Estruturas sem cadeia A: {no_chain_file}")

    if results["error_codes"]:
        error_file = logs_path / "cleaning_errors.txt"
        error_file.write_text("\n".join(sorted(results["error_codes"])) + "\n")
        print(f"Erros de processamento: {error_file}")

    success = results["success"]
    avg = results["total_residues"] / success if success > 0 else 0.0
    report_lines = [
        f"{'='*60}",
        "RELATÓRIO DE LIMPEZA - EVO-HEX",
        f"{'='*60}\n",
        f"Total de estruturas:      {results['total']:,}",
        f"Limpas com sucesso:       {results['success']:,}",
        f"Já existiam:              {results['exists']:,}",
        f"Sem cadeia A:             {results['no_chain_a']:,}",
        f"Vazias após limpeza:      {results['empty']:,}",
        f"Erros:                    {results['error']:,}",
        f"Total de resíduos:        {results['total_residues']:,}",
        f"Média resíduos/estrutura: {avg:.1f}",
        f"\n{'='*60}",
    ]
    report_file = logs_path / "cleaning_report.txt"
    report_file.write_text("\n".join(report_lines))
    print(f"Relatório de limpeza: {report_file}")


def print_cleaning_summary(results: dict, raw_dir: Path, clean_dir: Path) -> None:
    """Exibe resumo da limpeza no console.

    Args:
        results: Dict retornado por clean_all_structures.
        raw_dir: Diretório de PDBs brutos.
        clean_dir: Diretório de PDBs limpos.
    """
    total_ok = results["success"] + results["exists"]
    print(f"\n{'='*60}")
    print("RESUMO DA LIMPEZA")
    print(f"{'='*60}")
    print(f"\n  Raw (originais): {raw_dir}")
    print(f"  Clean (limpas):  {clean_dir}")
    print(f"\n  Total processadas:     {results['total']:,}")
    print(f"  Sucesso:               {results['success']:,}")
    print(f"  Já existiam:           {results['exists']:,}")
    print(f"  Sem cadeia A:          {results['no_chain_a']:,}")
    print(f"  Vazias:                {results['empty']:,}")
    print(f"  Erros:                 {results['error']:,}")
    print(f"\n  Estruturas disponíveis: {total_ok:,}")
    if results["total_residues"] > 0:
        print(f"  Total de resíduos:      {results['total_residues']:,}")
    print(f"{'='*60}\n")
