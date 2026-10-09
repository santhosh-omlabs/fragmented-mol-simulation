"""Step 6a: price the 14-qubit N-fragment job on a real IBM backend WITHOUT submitting anything.

Transpiles to the device's connectivity and native gates, then reports routed two-qubit gate count,
scheduled circuit duration, expected QPU seconds for a given shot count, and the clean-run fraction
implied by the device's own reported two-qubit errors. Only reads backend properties (no job, no QPU time).

Run:  python -m fragsim.ibm_dry_run [--backend ibm_fez] [--shots 4000]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from pyscf import scf
from qiskit import transpile
from qiskit_ibm_runtime import QiskitRuntimeService

from fragsim.embedding import build_atom_basis, density_in_basis
from fragsim.molecule import build_mol
from fragsim.solver import build_cluster_problem
from fragsim.sqd_circuit import canonicalize, lucj_circuit, lucj_operator

RECORD = Path(__file__).resolve().parents[2] / "data" / "ibm_dry_run.json"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--backend", default="ibm_fez")
    parser.add_argument("--shots", type=int, default=4000)
    args = parser.parse_args()

    mol = build_mol()
    mf = scf.RHF(mol).run()
    basis = build_atom_basis(mol, mf)
    raw = build_cluster_problem(mol, mf, basis, density_in_basis(basis, mf), [0], n_virtual=3)
    problem, _ = canonicalize(raw)
    circuit = lucj_circuit(problem, lucj_operator(problem, 2))

    backend = QiskitRuntimeService().backend(args.backend)
    target = backend.target
    two_q = next(n for n in target.operation_names if n in ("cz", "ecr", "cx"))
    errors = [p.error for p in target[two_q].values() if p is not None and p.error is not None]
    median_error = float(np.median(errors))

    results = {}
    for seed in range(5):  # routing is stochastic: report the spread
        isa = transpile(
            circuit,
            backend=backend,
            optimization_level=3,
            seed_transpiler=seed,
            scheduling_method="alap",
        )
        n2 = int(isa.count_ops().get(two_q, 0))
        results[seed] = (n2, isa.depth(), float(isa.duration * target.dt))
    best = min(results, key=lambda s: results[s][0])
    n2, _, seconds_per_shot_circuit = results[best]
    rep_delay = getattr(backend, "default_rep_delay", None) or 250e-6
    per_shot = seconds_per_shot_circuit + rep_delay
    qpu_seconds = args.shots * per_shot
    clean = float(np.prod([1 - median_error] * n2))
    print(f"backend {args.backend} ({backend.num_qubits} qubits), native 2-qubit gate: {two_q}")
    print(f"median {two_q} error reported by device: {median_error:.4f}")
    for s, (a, d, t) in results.items():
        print(f"  routing seed {s}: {a} {two_q}, depth {d}, circuit {t * 1e6:.0f} us")
    print(
        f"best: {n2} {two_q}; per shot {per_shot * 1e6:.0f} us (rep delay {rep_delay * 1e6:.0f} us)"
    )
    print(f"{args.shots} shots -> about {qpu_seconds:.1f} s of quantum execution time")
    print(
        f"clean-run fraction from device errors: {clean:.2e}  (expect {clean * args.shots:.0f} clean shots)"
    )

    RECORD.write_text(
        json.dumps(
            {
                "what": "Dry run: 14-qubit N-fragment LUCJ circuit routed for a real IBM backend; nothing submitted",
                "backend": args.backend,
                "native_two_qubit_gate": two_q,
                "median_two_qubit_error": median_error,
                "routing_seeds": {
                    str(s): {"two_qubit_gates": a, "depth": d, "circuit_seconds": t}
                    for s, (a, d, t) in results.items()
                },
                "best_two_qubit_gates": n2,
                "shots": args.shots,
                "estimated_quantum_seconds": qpu_seconds,
                "clean_fraction_estimate": clean,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"wrote {RECORD}")


if __name__ == "__main__":
    main()
