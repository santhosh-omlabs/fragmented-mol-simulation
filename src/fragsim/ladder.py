"""Fragment ladder: exact (FCI) cluster solves for growing fragments, truncated virtuals.

Run:  python -m fragsim.ladder          (about 4 minutes, writes data/ladder_exact_fragments.json)
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from pyscf import scf

from fragsim.embedding import build_atom_basis, density_in_basis
from fragsim.molecule import BASIS, build_mol, run_reference
from fragsim.solver import frozen_core_energy, solve_cluster

HARTREE_TO_KCAL = 627.5095
RECORD = Path(__file__).resolve().parents[2] / "data" / "ladder_exact_fragments.json"

# atom 0 = N, atoms 1-3 = H. Every rung is a partition of the four atoms.
RUNGS: dict[str, list[list[int]]] = {
    "1: N | H1 | H2 | H3": [[0], [1], [2], [3]],
    "2: N+H1 | H2 | H3": [[0, 1], [2], [3]],
    "3: N+H1+H2 | H3": [[0, 1, 2], [3]],
}
VIRTUALS_KEPT: list[int | None] = [0, 1, 2, 3, 4, 5, None]  # None = keep all empty orbitals


def run_ladder() -> dict:
    """Total energy and per-fragment qubit counts for every rung and truncation level."""
    mol = build_mol()
    mf = scf.RHF(mol).run()
    basis = build_atom_basis(mol, mf)
    d = density_in_basis(basis, mf)
    e_core = frozen_core_energy(mol, mf)
    ref = run_reference(mol)
    rows = []
    for rung, fragments in RUNGS.items():
        for k in VIRTUALS_KEPT:
            t0 = time.perf_counter()
            parts = [solve_cluster(mol, mf, basis, d, atoms, n_virtual=k) for atoms in fragments]
            energy = e_core + sum(p.e_frag for p in parts)
            rows.append(
                {
                    "rung": rung,
                    "virtuals_kept": k,
                    "energy_hartree": energy,
                    "error_kcal_mol": (energy - ref.e_fci) * HARTREE_TO_KCAL,
                    "qubits_per_fragment": [p.n_qubits for p in parts],
                    "seconds": round(time.perf_counter() - t0, 1),
                }
            )
            print(
                f"{rung:22s} k={k!s:5s} error {rows[-1]['error_kcal_mol']:8.2f} kcal/mol",
                flush=True,
            )
    return {
        "system": "NH3 pyramidal, N-H 1.012 A, H-N-H 106.7 deg",
        "basis": BASIS,
        "frozen_core": "N 1s",
        "e_hf": ref.e_hf,
        "e_fci": ref.e_fci,
        "solver": "FCI on fragment + bath; MP2 natural virtuals; projected-amplitude energy",
        "rows": rows,
    }


def main() -> None:
    record = run_ladder()
    RECORD.parent.mkdir(parents=True, exist_ok=True)
    RECORD.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {RECORD}")


if __name__ == "__main__":
    main()
