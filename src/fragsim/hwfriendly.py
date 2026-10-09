"""Step 7a: does a heavy-hex-friendly LUCJ circuit keep the physics and cut the hardware cost?

For the 14-qubit N fragment compare the all-to-all circuit with one restricted to heavy_hex_pairs:
  ideal quality   matched-size SQD error of the ideal (noiseless) circuit's string ranking
  hardware cost   two-qubit gates after routing to a real IBM backend (read-only, nothing submitted)

Run:  python -m fragsim.hwfriendly [--backend ibm_fez]   (writes data/hwfriendly_study.json)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import ffsim
import numpy as np
from pyscf import fci, scf
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_ibm_runtime import QiskitRuntimeService

from fragsim.embedding import build_atom_basis, density_in_basis
from fragsim.matched import _error, circuit_order
from fragsim.molecule import build_mol
from fragsim.solver import build_cluster_problem
from fragsim.sqd_circuit import (
    canonicalize,
    final_state,
    heavy_hex_pairs,
    lucj_circuit,
    lucj_operator,
)

DATA = Path(__file__).resolve().parents[2] / "data"
KS = [4, 6, 8, 12, 16, 20]


def route(circuit, backend, seed: int):
    """ISA circuit using ffsim's pre-init passes (they merge and decompose orbital rotations well)."""
    manager = generate_preset_pass_manager(
        optimization_level=3, backend=backend, seed_transpiler=seed
    )
    manager.pre_init = ffsim.qiskit.PRE_INIT
    return manager.run(circuit)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--backend", default="ibm_fez")
    backend_name = parser.parse_args().backend
    backend = QiskitRuntimeService().backend(backend_name)

    mol = build_mol()
    mf = scf.RHF(mol).run()
    basis = build_atom_basis(mol, mf)
    raw = build_cluster_problem(mol, mf, basis, density_in_basis(basis, mf), [0], n_virtual=3)
    problem, _ = canonicalize(raw)
    e_fci, _ = fci.direct_spin0.FCI(mol).kernel(problem.h1, problem.eri, 7, problem.nelec)
    out = {
        "what": "all-to-all vs heavy-hex-restricted LUCJ, 14-qubit N fragment",
        "backend": backend_name,
        "variants": {},
    }
    for name, pairs in (
        ("all-to-all", None),
        ("heavy-hex pairs", heavy_hex_pairs(problem.n_orbitals)),
    ):
        for reps in (1, 2):
            op = lucj_operator(problem, reps, interaction_pairs=pairs)
            order = circuit_order(problem, final_state(problem, op))
            quality = {str(k): _error(problem, order[:k], e_fci) for k in KS}
            circuit = lucj_circuit(problem, op)
            routed = [route(circuit, backend, s) for s in range(5)]
            cz = [int(r.count_ops().get("cz", 0)) for r in routed]
            depth = [r.depth() for r in routed]
            key = f"{name}, {reps} layer(s)"
            out["variants"][key] = {
                "ideal_matched_error_kcal_mol": quality,
                "cz_per_seed": cz,
                "best_cz": min(cz),
                "depth_at_best": depth[int(np.argmin(cz))],
            }
            print(
                f"{key:28s} CZ best {min(cz):5d} (seeds {cz})  ideal k=8 {quality['8']:.2f} k=12 {quality['12']:.2f} k=16 {quality['16']:.2f}",
                flush=True,
            )
    (DATA / "hwfriendly_study.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
