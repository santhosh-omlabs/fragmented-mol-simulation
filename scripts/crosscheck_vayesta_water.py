"""Cross-check fragsim's embedding against Vayesta (EWF, FCI solver, DMET bath) on water clusters.

Both sides run all-electron in STO-3G with exact FCI on every cluster and no truncation. STO-3G is used
because all-electron 6-31G dimer clusters (17 orbitals, 14 electrons, 3.8e8 determinants) do not fit in
7 GB; the check is of the algorithm on molecules without symmetry, not of 6-31G numbers.
The frozen two-oxygen core of the 6-31G runs is covered by the identities in tests/test_water_embedding.py.

Needs a Vayesta source checkout:  export VAYESTA_PATH=<folder>   (see scripts/crosscheck_vayesta.py).

Run:  python scripts/crosscheck_vayesta_water.py   (a few minutes, writes data/crosscheck_vayesta_water.json)
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

from pyscf import scf

from fragsim.embedding import build_atom_basis, density_in_basis, make_cluster
from fragsim.molecule import build_mol
from fragsim.solver import frozen_core_energy, solve_cluster
from fragsim.water import dimer_atoms, flatten, ring_atoms, water_fragments

HARTREE_TO_KCAL = 627.5095
BASIS = "sto-3g"
RECORD = Path(__file__).resolve().parents[1] / "data" / "crosscheck_vayesta_water.json"
SYSTEMS = {"dimer": dimer_atoms(), "trimer ring": ring_atoms(3)}


def fragsim_energy(mol, mf, partition) -> tuple[float, list[int]]:
    """Our one-shot embedding, all-electron (n_frozen=0). Also returns the qubits per cluster."""
    basis = build_atom_basis(mol, mf, n_frozen=0)
    d = density_in_basis(basis, mf, n_frozen=0)
    parts = [solve_cluster(mol, mf, basis, d, atoms, n_frozen=0) for atoms in partition]
    qubits = [2 * make_cluster(basis, d, atoms).n_orbitals for atoms in partition]
    return frozen_core_energy(mol, mf, n_frozen=0) + sum(p.e_frag for p in parts), qubits


def vayesta_energy(mf, partition) -> float:
    import vayesta.ewf

    emb = vayesta.ewf.EWF(mf, solver="FCI", bath_options={"bathtype": "dmet"})
    with emb.iaopao_fragmentation() as frag:
        for atoms in partition:
            frag.add_atomic_fragment(atoms)
    emb.kernel()
    return float(emb.e_tot)


def main() -> None:
    path = os.environ.get("VAYESTA_PATH")
    if not path:
        sys.exit("Set VAYESTA_PATH to a Vayesta source checkout.")
    sys.path.insert(0, path)
    rows = []
    print(
        f"{'system':14s} {'fragsim':>14s} {'Vayesta':>14s} {'diff (kcal/mol)':>16s}  qubits per cluster"
    )
    for name, waters in SYSTEMS.items():
        mol = build_mol(flatten(waters), BASIS)
        mf = scf.RHF(mol).run()
        partition = water_fragments(len(waters))
        t0 = time.perf_counter()
        ours, qubits = fragsim_energy(mol, mf, partition)
        theirs = vayesta_energy(mf, partition)
        rows.append(
            {
                "system": name,
                "partition": partition,
                "e_hf": float(mf.e_tot),
                "fragsim_hartree": ours,
                "vayesta_hartree": theirs,
                "difference_kcal_mol": (ours - theirs) * HARTREE_TO_KCAL,
                "qubits_per_cluster": qubits,
                "seconds": round(time.perf_counter() - t0, 1),
            }
        )
        print(
            f"{name:14s} {ours:14.8f} {theirs:14.8f} {rows[-1]['difference_kcal_mol']:16.3f}  {qubits}",
            flush=True,
        )
    RECORD.write_text(
        json.dumps(
            {
                "what": "fragsim vs Vayesta EWF (FCI, DMET bath), all-electron STO-3G, one fragment per water",
                "vayesta_commit": "7f1639d",
                "rows": rows,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"wrote {RECORD}")


if __name__ == "__main__":
    main()
