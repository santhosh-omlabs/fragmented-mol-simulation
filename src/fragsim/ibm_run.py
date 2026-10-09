"""Step 6b: run the 14-qubit N-fragment circuit once on IBM hardware (raw, no error mitigation).

Submits ONE Sampler job, writes the job id to disk before waiting, then saves the raw counts.
Only run this after reading data/ibm_dry_run.json.

Run:  python -m fragsim.ibm_run [--backend ibm_fez] [--shots 4000] [--seed 1]
      python -m fragsim.ibm_run --job-id <id>      (fetch a job already submitted)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from pyscf import scf
from qiskit import transpile
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2

from fragsim.embedding import build_atom_basis, density_in_basis
from fragsim.molecule import build_mol
from fragsim.solver import build_cluster_problem
from fragsim.sqd_circuit import canonicalize, lucj_circuit, lucj_operator

DATA = Path(__file__).resolve().parents[2] / "data"
RECORD = DATA / "ibm_run.json"
COUNTS = DATA / "ibm_counts_n_fragment.json"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--backend", default="ibm_fez")
    parser.add_argument("--shots", type=int, default=4000)
    parser.add_argument(
        "--seed", type=int, default=1, help="routing seed (1 was best in the dry run)"
    )
    parser.add_argument("--job-id", default=None)
    args = parser.parse_args()
    service = QiskitRuntimeService()

    if args.job_id:
        job = service.job(args.job_id)
        n2 = None
    else:
        mol = build_mol()
        mf = scf.RHF(mol).run()
        basis = build_atom_basis(mol, mf)
        raw = build_cluster_problem(mol, mf, basis, density_in_basis(basis, mf), [0], n_virtual=3)
        problem, _ = canonicalize(raw)
        circuit = lucj_circuit(problem, lucj_operator(problem, 2))
        backend = service.backend(args.backend)
        isa = transpile(circuit, backend=backend, optimization_level=3, seed_transpiler=args.seed)
        n2 = int(isa.count_ops().get("cz", 0))
        sampler = SamplerV2(mode=backend)  # defaults: no dynamical decoupling, no twirling
        job = sampler.run([isa], shots=args.shots)
        DATA.mkdir(parents=True, exist_ok=True)
        (DATA / "ibm_job_id.txt").write_text(job.job_id() + "\n", encoding="utf-8")
        print(
            f"submitted job {job.job_id()} on {args.backend}: {n2} CZ, {args.shots} shots",
            flush=True,
        )

    result = job.result()
    counts = result[0].data.meas.get_counts()
    COUNTS.write_text(json.dumps(counts, indent=1) + "\n", encoding="utf-8")
    usage = None
    try:
        usage = job.usage()
    except Exception as exc:  # noqa: BLE001 - usage can lag behind completion
        print("usage not available yet:", exc)
    record = {
        "what": "14-qubit N-fragment LUCJ circuit on IBM hardware, raw, no error mitigation",
        "job_id": job.job_id(),
        "backend": job.backend().name,
        "shots": sum(counts.values()),
        "two_qubit_gates_cz": n2,
        "quantum_seconds_billed": usage,
        "counts_file": COUNTS.name,
    }
    RECORD.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
