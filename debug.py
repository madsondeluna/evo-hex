#!/usr/bin/env python3
"""
Sessao de debug do pipeline Evo-Hex.

Amostra N PDBs do structures_clean/, roda DSSP unificado e salva tudo em
debug/analysis/ para validacao rapida de graficos sem re-processar a base inteira.

Uso
---
    python debug.py              # 100 PDBs (padrao)
    python debug.py --n 50       # 50 PDBs
    python debug.py --n 200      # 200 PDBs
    python debug.py --clean      # apaga debug/ e recria
"""

import argparse
import logging
import random
import shutil
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT))

from cath_analysis.config import STRUCTURES_CLEAN_PATH, PROCESS_WORKERS, glob_pdb

DEBUG_ROOT        = ROOT / "debug"
DEBUG_CLEAN_PATH  = DEBUG_ROOT / "structures_clean"
DEBUG_ANALYSIS    = DEBUG_ROOT / "analysis"


def sample_pdbs(n: int, seed: int = 42) -> list[Path]:
    all_pdbs = glob_pdb(STRUCTURES_CLEAN_PATH)
    if not all_pdbs:
        raise FileNotFoundError(
            f"Nenhum PDB encontrado em {STRUCTURES_CLEAN_PATH}.\n"
            "Execute as etapas 1 e 2 do pipeline principal primeiro."
        )
    if n >= len(all_pdbs):
        logger.warning(
            "N=%d >= total disponivel (%d). Usando todos.", n, len(all_pdbs)
        )
        return all_pdbs
    random.seed(seed)
    return random.sample(all_pdbs, n)


def setup_debug_clean(pdbs: list[Path], force: bool = False) -> None:
    if DEBUG_CLEAN_PATH.exists() and not force:
        existing = len(glob_pdb(DEBUG_CLEAN_PATH))
        print(f"\n  debug/structures_clean/ ja existe com {existing} PDBs.")
        resp = input("  Recriar? [s/N] ").strip().lower()
        if resp not in ("s", "sim", "y", "yes"):
            print("  Usando PDBs existentes.")
            return

    if DEBUG_CLEAN_PATH.exists():
        shutil.rmtree(DEBUG_CLEAN_PATH)
    DEBUG_CLEAN_PATH.mkdir(parents=True)

    print(f"\n  Copiando {len(pdbs)} PDBs para debug/structures_clean/...")
    for pdb in pdbs:
        shutil.copy2(pdb, DEBUG_CLEAN_PATH / pdb.name)
    print(f"  Pronto: {len(pdbs)} PDBs copiados.")


def run_debug_dssp() -> dict:
    from cath_analysis.unified_dssp import collect_all, save_all_results

    report = DEBUG_ANALYSIS / "unified_report.txt"
    if report.exists():
        print(f"\n  Resultados existentes em debug/analysis/.")
        resp = input("  Re-executar DSSP? [s/N] ").strip().lower()
        if resp not in ("s", "sim", "y", "yes"):
            print("  Carregando resultados anteriores...")
            import pickle, pandas as pd
            cache = {}
            for name in (
                "helix_propensities.csv",
                "helix_positions.csv",
                "helix_type_composition.csv",
                "helix_type_top_residues.csv",
                "helix_type_statistical_comparison.csv",
            ):
                p = DEBUG_ANALYSIS / name
                if p.exists():
                    cache[name.replace(".csv", "")] = pd.read_csv(p)
            pkl = DEBUG_ANALYSIS / "evo_data.pkl"
            if pkl.exists():
                with open(pkl, "rb") as f:
                    cache["evo_data"] = pickle.load(f)
            return cache

        shutil.rmtree(DEBUG_ANALYSIS, ignore_errors=True)

    DEBUG_ANALYSIS.mkdir(parents=True, exist_ok=True)

    print(f"\n  Rodando DSSP em {len(glob_pdb(DEBUG_CLEAN_PATH))} PDBs...")
    agg = collect_all(DEBUG_CLEAN_PATH, workers=min(4, PROCESS_WORKERS))
    results = save_all_results(agg, DEBUG_ANALYSIS)
    print(f"\n  Resultados salvos em: {DEBUG_ANALYSIS}")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Debug Evo-Hex com subset de PDBs")
    parser.add_argument("--n",     type=int, default=100, help="Numero de PDBs (padrao: 100)")
    parser.add_argument("--seed",  type=int, default=42,  help="Semente aleatoria (padrao: 42)")
    parser.add_argument("--clean", action="store_true",   help="Apaga debug/ e recria do zero")
    args = parser.parse_args()

    if args.clean and DEBUG_ROOT.exists():
        shutil.rmtree(DEBUG_ROOT)
        print(f"  debug/ removido.")

    print("=" * 60)
    print("DEBUG EVO-HEX")
    print(f"  N PDBs:   {args.n}")
    print(f"  Seed:     {args.seed}")
    print(f"  Saida:    {DEBUG_ROOT}")
    print("=" * 60)

    pdbs = sample_pdbs(args.n, seed=args.seed)
    setup_debug_clean(pdbs, force=args.clean)
    run_debug_dssp()

    print("\n" + "=" * 60)
    print("DEBUG CONCLUIDO")
    print(f"  Abra o plots.ipynb e mude ANALYSIS para:")
    print(f"  ANALYSIS = ROOT / 'debug' / 'analysis'")
    print("=" * 60)


if __name__ == "__main__":
    main()
