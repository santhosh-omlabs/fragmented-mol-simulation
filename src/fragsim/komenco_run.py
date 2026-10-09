"""Step 5: run the 14-qubit N-fragment circuit through the Komenco gateway (pipeline demonstration).

One POST sends the LUCJ circuit; the gateway returns bitstring probabilities. We then

  1. compare them with our own exact simulation of the same circuit (total variation distance),
  2. draw shots from them and run the normal SQD loop,
  3. rank strings by gateway probability and compare with the local ideal ranking.

Komenco is a third-party platform provided by Automatski (https://automatski.com/platform.html).
This is a plumbing test of "circuit out, samples in, SQD on top". The gateway is treated as a
classical emulator, so the result is not evidence about quantum hardware.

Run:  python -m fragsim.komenco_run --dry-run   (build and size the request, no network)
      python -m fragsim.komenco_run             (one POST; writes data/komenco_run.json)
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from pyscf import fci, scf
from pyscf.fci import cistring

from fragsim.backends.komenco import KomencoClient, circuit_request
from fragsim.embedding import build_atom_basis, density_in_basis
from fragsim.matched import HARTREE_TO_KCAL, _error, circuit_order
from fragsim.molecule import build_mol
from fragsim.sampling import has_correct_electron_count
from fragsim.solver import build_cluster_problem
from fragsim.sqd import run_sqd
from fragsim.sqd_circuit import canonicalize, final_state, lucj_circuit, lucj_operator

DATA = Path(__file__).resolve().parents[2] / "data"
RECORD = DATA / "komenco_run.json"
SHOTS = 10_000
KS = [4, 6, 8, 10, 12, 16, 20, 35]


def _local_probabilities(problem, state: np.ndarray) -> dict[str, float]:
    """Exact bitstring probabilities of the ideal circuit, Qiskit order [beta][alpha]."""
    norb = problem.n_orbitals
    strings = cistring.make_strings(range(norb), problem.nelec[0])
    probs = np.abs(state) ** 2
    out = {}
    for row, col in zip(*np.nonzero(probs > 1e-15)):
        key = format(int(strings[col]), f"0{norb}b") + format(int(strings[row]), f"0{norb}b")
        out[key] = float(probs[row, col])
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true", help="build the request, make no call")
    dry = parser.parse_args().dry_run

    mol = build_mol()
    mf = scf.RHF(mol).run()
    basis = build_atom_basis(mol, mf)
    raw = build_cluster_problem(mol, mf, basis, density_in_basis(basis, mf), [0], n_virtual=3)
    problem, _ = canonicalize(raw)
    op = lucj_operator(problem, 2)
    circuit = lucj_circuit(problem, op)
    body = circuit_request(circuit, shots=SHOTS, top_k=1 << circuit.num_qubits, api_key="open")
    n_cx = sum(1 for g in body["operations"] if g["gate"] == "cx")
    size_kb = len(json.dumps(body)) / 1024
    print(
        f"request: {body['num_qubits']} qubits, {len(body['operations'])} gates ({n_cx} CX), "
        f"{size_kb:.0f} kB"
    )
    if dry:
        return

    client = KomencoClient(shots=SHOTS)
    started = time.perf_counter()
    remote = client.run_probabilities(circuit)
    seconds = time.perf_counter() - started
    print(f"gateway answered in {seconds:.1f} s with {len(remote)} bitstrings")

    state = final_state(problem, op)
    local = _local_probabilities(problem, state)
    keys = set(local) | set(remote)
    tv = 0.5 * sum(abs(local.get(k, 0.0) - remote.get(k, 0.0)) for k in keys)
    print(f"total variation distance, gateway vs local exact: {tv:.2e}")

    counts = client.sample_counts(circuit, SHOTS, seed=0)
    best, _ = run_sqd(problem, counts, iterations=5, samples_per_batch=200, num_batches=3)
    e_fci, _ = fci.direct_spin0.FCI(mol).kernel(problem.h1, problem.eri, 7, problem.nelec)
    sqd_error = (best.energy - e_fci) * HARTREE_TO_KCAL
    print(f"SQD on gateway samples: {best.n_strings} strings, error {sqd_error:.3f} kcal/mol")

    # rank strings by gateway probability (both spin halves), compare with local ideal ranking
    norb = problem.n_orbitals
    mass: dict[int, float] = {}
    for bits, p in remote.items():
        if not has_correct_electron_count(bits, norb, problem.nelec):
            continue  # the gateway returns all 2^14 bitstrings; wrong-count ones are ~zero
        for half in (bits[:norb], bits[norb:]):
            mass[int(half, 2)] = mass.get(int(half, 2), 0.0) + p
    remote_order = sorted(mass, key=lambda s: -mass[s])
    local_order = circuit_order(problem, state)
    rows = []
    for k in KS:
        rows.append(
            {
                "k": k,
                "gateway_ranked": _error(problem, remote_order[:k], e_fci),
                "local_ideal_ranked": _error(problem, local_order[:k], e_fci),
            }
        )
        print(rows[-1])

    record = {
        "what": "14-qubit N-fragment LUCJ circuit through the Komenco gateway (pipeline demonstration)",
        "gateway": f"http://{client.host}:{client.port}/api/komenco (vendor trial, treated as a classical emulator)",
        "requests": client.n_requests,
        "qubits": body["num_qubits"],
        "gates": len(body["operations"]),
        "cx_gates": n_cx,
        "seconds": seconds,
        "bitstrings_returned": len(remote),
        "total_variation_vs_local_exact": tv,
        "shots_drawn": SHOTS,
        "sqd_strings": best.n_strings,
        "sqd_error_kcal_mol": sqd_error,
        "rows": rows,
    }
    RECORD.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {RECORD}")


if __name__ == "__main__":
    main()
