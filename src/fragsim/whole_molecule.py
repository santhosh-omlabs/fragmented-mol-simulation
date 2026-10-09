"""Step 8: the whole molecule (28 qubits, no fragmentation) with the same SQD pipeline.

Closes the three-way comparison from the README: exact FCI, SQD on the whole molecule, SQD on fragments.

  python -m fragsim.whole_molecule dry      ideal circuit quality + routed gate count (no job)
  python -m fragsim.whole_molecule submit   ONE IBM job, writes data/whole_molecule_job_id.txt
  python -m fragsim.whole_molecule analyse  [--job-id id]  hardware counts -> data/whole_molecule.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from pyscf import fci, scf
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2

from fragsim.embedding import build_atom_basis, density_in_basis
from fragsim.hwfriendly import route
from fragsim.matched import (
    HARTREE_TO_KCAL,
    _error,
    ccsd_order,
    circuit_order,
    excitation_order,
    noisy_order,
)
from fragsim.molecule import build_mol
from fragsim.sampling import has_correct_electron_count
from fragsim.solver import build_cluster_problem
from fragsim.sqd import run_sqd
from fragsim.sqd_circuit import (
    canonicalize,
    final_state,
    heavy_hex_pairs,
    lucj_circuit,
    lucj_operator,
)

DATA = Path(__file__).resolve().parents[2] / "data"
RECORD = DATA / "whole_molecule.json"
COUNTS = DATA / "whole_molecule_counts.json"
KS = [8, 16, 32, 64, 128, 256]
SHOTS = 10_000
SEED = 1


def _setup():
    mol = build_mol()
    mf = scf.RHF(mol).run()
    basis = build_atom_basis(mol, mf)
    raw = build_cluster_problem(mol, mf, basis, density_in_basis(basis, mf), [0, 1, 2, 3])
    problem, eps = canonicalize(raw)
    e_fci, _ = fci.direct_spin0.FCI(mol).kernel(
        problem.h1, problem.eri, problem.n_orbitals, problem.nelec
    )
    op = lucj_operator(problem, 1, interaction_pairs=heavy_hex_pairs(problem.n_orbitals))
    return problem, eps, e_fci, op


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["dry", "submit", "analyse"])
    parser.add_argument("--backend", default="ibm_fez")
    parser.add_argument("--job-id", default=None)
    args = parser.parse_args()
    problem, eps, e_fci, op = _setup()
    print(
        f"{problem.n_qubits} qubits, {problem.n_orbitals} orbitals, nelec {problem.nelec}",
        flush=True,
    )

    if args.mode == "dry":
        state = final_state(problem, op)
        orders = {
            "ideal heavy-hex circuit": circuit_order(problem, state),
            "excitation order": excitation_order(problem, eps),
            "CCSD-ranked": ccsd_order(problem, eps),
        }
        ideal = {n: {str(k): _error(problem, o[:k], e_fci) for k in KS} for n, o in orders.items()}
        for n, row in ideal.items():
            print(f"{n:26s}", "  ".join(f"k={k}: {v:6.2f}" for k, v in row.items()), flush=True)
        backend = QiskitRuntimeService().backend(args.backend)
        cz = [
            int(route(lucj_circuit(problem, op), backend, s).count_ops().get("cz", 0))
            for s in range(5)
        ]
        print("routed CZ per seed:", cz)
        (DATA / "whole_molecule_dry.json").write_text(
            json.dumps(
                {"ideal_matched_error_kcal_mol": ideal, "routed_cz": cz, "k_values": KS}, indent=2
            )
            + "\n",
            encoding="utf-8",
        )
    elif args.mode == "submit":
        service = QiskitRuntimeService()
        backend = service.backend(args.backend)
        isa = route(lucj_circuit(problem, op), backend, SEED)
        job = SamplerV2(mode=backend).run([isa], shots=SHOTS)
        (DATA / "whole_molecule_job_id.txt").write_text(job.job_id() + "\n", encoding="utf-8")
        print(
            f"submitted {job.job_id()}: {isa.count_ops().get('cz', 0)} CZ, {SHOTS} shots",
            flush=True,
        )
        COUNTS.write_text(
            json.dumps(job.result()[0].data.meas.get_counts()) + "\n", encoding="utf-8"
        )
        print("usage:", job.usage())
    else:
        service = QiskitRuntimeService()
        if args.job_id:
            COUNTS.write_text(
                json.dumps(service.job(args.job_id).result()[0].data.meas.get_counts()) + "\n",
                encoding="utf-8",
            )
        counts = json.loads(COUNTS.read_text(encoding="utf-8"))
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
        matched = {str(k): _error(problem, order[:k], e_fci) for k in KS if k <= len(order)}
        best, history = run_sqd(problem, counts, iterations=5, samples_per_batch=300, num_batches=3)
        out = {
            "what": "whole molecule, 28 qubits, heavy-hex LUCJ 1 layer, IBM hardware raw",
            "shots": total,
            "right_count_fraction": right,
            "distinct_strings_after_recovery": len(order),
            "matched_error_kcal_mol": matched,
            "sqd_loop": {
                "energy_error_kcal_mol": (best.energy - e_fci) * HARTREE_TO_KCAL,
                "strings": best.n_strings,
                "determinants": best.dimension,
                "history_kcal_mol": [(h.energy - e_fci) * HARTREE_TO_KCAL for h in history],
            },
        }
        # controls at the SAME subspace size the hardware loop reached
        k = best.n_strings
        rng = np.random.default_rng(1)
        random_counts: dict[str, int] = {}
        for _ in range(total):
            b = "".join(map(str, rng.integers(0, 2, 2 * problem.n_orbitals)))
            random_counts[b] = random_counts.get(b, 0) + 1
        try:
            rnd, _ = run_sqd(
                problem, random_counts, iterations=5, samples_per_batch=300, num_batches=3
            )
            random_error = (rnd.energy - e_fci) * HARTREE_TO_KCAL
        except ValueError:
            random_error = None
        out["controls_at_k"] = {
            "k": k,
            "excitation_order": _error(problem, excitation_order(problem, eps)[:k], e_fci),
            "ccsd_ranked": _error(problem, ccsd_order(problem, eps)[:k], e_fci),
            "random_bits_sqd_loop": random_error,
            "random_bits_right_count_fraction": sum(
                c
                for b, c in random_counts.items()
                if has_correct_electron_count(b, problem.n_orbitals, problem.nelec)
            )
            / total,
        }
        RECORD.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
