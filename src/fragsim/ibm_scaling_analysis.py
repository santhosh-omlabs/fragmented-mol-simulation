"""Step 7c: right-electron-count fraction and matched-size error against two-qubit gate count.

Run:  python -m fragsim.ibm_scaling_analysis   (writes data/ibm_scaling_analysis.json)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from pyscf import fci, scf

from fragsim.embedding import build_atom_basis, density_in_basis
from fragsim.matched import _error, noisy_order
from fragsim.molecule import build_mol
from fragsim.sampling import has_correct_electron_count
from fragsim.solver import build_cluster_problem
from fragsim.sqd_circuit import canonicalize

DATA = Path(__file__).resolve().parents[2] / "data"
KS = [4, 6, 8, 12, 16]


def main() -> None:
    mol = build_mol()
    mf = scf.RHF(mol).run()
    basis = build_atom_basis(mol, mf)
    raw = build_cluster_problem(mol, mf, basis, density_in_basis(basis, mf), [0], n_virtual=3)
    problem, _ = canonicalize(raw)
    e_fci, _ = fci.direct_spin0.FCI(mol).kernel(problem.h1, problem.eri, 7, problem.nelec)
    run = json.loads((DATA / "ibm_scaling_run.json").read_text(encoding="utf-8"))
    sets = dict(zip(run["labels"], zip(run["cz"], [run["counts"][k] for k in run["labels"]])))
    sets["all-to-all 2 layers (first job)"] = (
        1107,
        json.loads((DATA / "ibm_counts_n_fragment.json").read_text(encoding="utf-8")),
    )
    rows = {}
    for name, (cz, counts) in sets.items():
        total = sum(counts.values())
        right = (
            sum(
                c
                for b, c in counts.items()
                if has_correct_electron_count(b, problem.n_orbitals, problem.nelec)
            )
            / total
        )
        order = noisy_order(problem, counts)
        rows[name] = {
            "cz": cz,
            "right_count_fraction": right,
            "matched": {str(k): _error(problem, order[:k], e_fci) for k in KS},
        }
        m = rows[name]["matched"]
        print(
            f"{name:32s} CZ {cz:5d}  right count {right:.3f}  "
            + "  ".join(f"k={k}: {m[str(k)]:5.1f}" for k in KS)
        )
    hh = [v for k, v in rows.items() if k.startswith("heavy-hex")]
    summary = {
        "heavy_hex_right_count_mean": float(np.mean([v["right_count_fraction"] for v in hh])),
        "heavy_hex_right_count_std": float(np.std([v["right_count_fraction"] for v in hh], ddof=1)),
        **{f"heavy_hex_k{k}_mean": float(np.mean([v["matched"][str(k)] for v in hh])) for k in KS},
        **{
            f"heavy_hex_k{k}_std": float(np.std([v["matched"][str(k)] for v in hh], ddof=1))
            for k in KS
        },
    }
    print(summary)
    (DATA / "ibm_scaling_analysis.json").write_text(
        json.dumps({"rows": rows, "heavy_hex_summary": summary}, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
