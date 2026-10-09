"""Step 6c: analyse the IBM hardware counts against the noisy simulation, a random control and classical orders.

Run:  python -m fragsim.ibm_analysis   (writes data/ibm_analysis.json)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from pyscf import fci, scf

from fragsim.embedding import build_atom_basis, density_in_basis
from fragsim.matched import HARTREE_TO_KCAL, KS, NOISY_COUNTS, _error, noisy_order
from fragsim.molecule import build_mol
from fragsim.sampling import has_correct_electron_count
from fragsim.solver import build_cluster_problem
from fragsim.sqd import run_sqd
from fragsim.sqd_circuit import canonicalize

DATA = Path(__file__).resolve().parents[2] / "data"


def _right_fraction(counts: dict[str, int], problem) -> float:
    total = sum(counts.values())
    good = sum(
        c
        for b, c in counts.items()
        if has_correct_electron_count(b, problem.n_orbitals, problem.nelec)
    )
    return good / total


def main() -> None:
    mol = build_mol()
    mf = scf.RHF(mol).run()
    basis = build_atom_basis(mol, mf)
    raw = build_cluster_problem(mol, mf, basis, density_in_basis(basis, mf), [0], n_virtual=3)
    problem, _ = canonicalize(raw)
    e_fci, _ = fci.direct_spin0.FCI(mol).kernel(problem.h1, problem.eri, 7, problem.nelec)
    hardware = json.loads((DATA / "ibm_counts_n_fragment.json").read_text(encoding="utf-8"))
    simulated = json.loads(NOISY_COUNTS.read_text(encoding="utf-8"))
    rng = np.random.default_rng(1)
    uniform = {}
    for _ in range(4000):
        b = "".join(map(str, rng.integers(0, 2, 2 * problem.n_orbitals)))
        uniform[b] = uniform.get(b, 0) + 1
    out = {
        "what": "IBM hardware counts vs noisy simulation vs random bits (14-qubit N fragment)",
        "sources": {},
    }
    for name, counts in (
        ("ibm_fez", hardware),
        ("noisy_simulation_0.5pct", simulated),
        ("random_bits", uniform),
    ):
        try:
            best, _ = run_sqd(problem, counts, iterations=5, samples_per_batch=200, num_batches=3)
        except ValueError as exc:  # recovery on pure random bits can hand back malformed batches
            print(f"{name}: SQD loop failed ({exc})")
            best = None
        try:
            order = noisy_order(problem, counts)
        except ValueError:  # same failure mode as above
            order = []
        out["sources"][name] = {
            "shots": sum(counts.values()),
            "fraction_right_electron_count": _right_fraction(counts, problem),
            "sqd_error_kcal_mol": (best.energy - e_fci) * HARTREE_TO_KCAL if best else None,
            "sqd_strings": best.n_strings if best else None,
            "matched": {str(k): _error(problem, order[:k], e_fci) for k in KS if k <= len(order)},
        }
        s = out["sources"][name]
        print(
            f"{name:26s} right count {s['fraction_right_electron_count']:.3f}  SQD {s['sqd_error_kcal_mol']} kcal/mol "
            f"({s['sqd_strings']} strings)  matched k=8 {s['matched'].get('8')} k=12 {s['matched'].get('12')}"
        )
    (DATA / "ibm_analysis.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
