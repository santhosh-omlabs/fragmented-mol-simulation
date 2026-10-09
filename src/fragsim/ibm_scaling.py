"""Step 7b: ONE IBM job, four circuits: noise against two-qubit gate count.

  heavy-hex pairs, 1 layer, three routing seeds (error bars)    ~160 CZ
  all-to-all, 1 layer                                           ~485 CZ
(the earlier all-to-all 2-layer job, ~1,107 CZ, is the third point)

Run:  python -m fragsim.ibm_scaling [--backend ibm_fez] [--shots 4000]
      python -m fragsim.ibm_scaling --job-id <id>
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from pyscf import scf
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2

from fragsim.embedding import build_atom_basis, density_in_basis
from fragsim.hwfriendly import route
from fragsim.molecule import build_mol
from fragsim.solver import build_cluster_problem
from fragsim.sqd_circuit import canonicalize, heavy_hex_pairs, lucj_circuit, lucj_operator

DATA = Path(__file__).resolve().parents[2] / "data"
RECORD = DATA / "ibm_scaling_run.json"
LABELS = ["heavy-hex seed 1", "heavy-hex seed 4", "heavy-hex seed 0", "all-to-all 1 layer"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--backend", default="ibm_fez")
    parser.add_argument("--shots", type=int, default=4000)
    parser.add_argument("--job-id", default=None)
    args = parser.parse_args()
    service = QiskitRuntimeService()
    cz = None
    if args.job_id:
        job = service.job(args.job_id)
    else:
        mol = build_mol()
        mf = scf.RHF(mol).run()
        basis = build_atom_basis(mol, mf)
        raw = build_cluster_problem(mol, mf, basis, density_in_basis(basis, mf), [0], n_virtual=3)
        problem, _ = canonicalize(raw)
        backend = service.backend(args.backend)
        hh = lucj_circuit(problem, lucj_operator(problem, 1, interaction_pairs=heavy_hex_pairs(7)))
        a2a = lucj_circuit(problem, lucj_operator(problem, 1))
        isa = [
            route(hh, backend, 1),
            route(hh, backend, 4),
            route(hh, backend, 0),
            route(a2a, backend, 3),
        ]
        cz = [int(c.count_ops().get("cz", 0)) for c in isa]
        job = SamplerV2(mode=backend).run(isa, shots=args.shots)
        (DATA / "ibm_scaling_job_id.txt").write_text(job.job_id() + "\n", encoding="utf-8")
        print(f"submitted {job.job_id()} on {args.backend}: CZ per circuit {cz}", flush=True)
    result = job.result()
    counts = {LABELS[i]: result[i].data.meas.get_counts() for i in range(len(LABELS))}
    try:
        usage = job.usage()
    except Exception as exc:  # noqa: BLE001 - usage can lag behind completion
        usage = f"unavailable: {exc}"
    RECORD.write_text(
        json.dumps(
            {
                "what": "noise vs two-qubit gate count, 14-qubit N fragment, raw, no mitigation",
                "job_id": job.job_id(),
                "backend": job.backend().name,
                "shots_each": args.shots,
                "labels": LABELS,
                "cz": cz,
                "quantum_seconds_billed": usage,
                "counts": counts,
            },
            indent=1,
        )
        + "\n",
        encoding="utf-8",
    )
    print("usage:", usage)


if __name__ == "__main__":
    main()
