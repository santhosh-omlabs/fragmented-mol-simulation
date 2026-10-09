"""Cross-check fragsim's embedding against Vayesta (EWF, FCI solver, DMET bath).

Both sides run all-electron (nothing frozen), because Vayesta's "frozen" option means something
else. Same molecule, same basis, exact FCI on every cluster, no truncation of empty orbitals.

Vayesta is not on PyPI and needs a BLAS dev library to build, so it is used from a source
checkout:  git clone https://github.com/BoothGroup/Vayesta  and  export VAYESTA_PATH=<that folder>.

Run:  python scripts/crosscheck_vayesta.py
"""

from __future__ import annotations

import os
import sys

from pyscf import scf

from fragsim.embedding import build_atom_basis, density_in_basis
from fragsim.molecule import build_mol
from fragsim.solver import frozen_core_energy, solve_cluster

HARTREE_TO_KCAL = 627.5095
PARTITIONS = {
    "N | H1 | H2 | H3": [[0], [1], [2], [3]],
    "N+H1 | H2 | H3": [[0, 1], [2], [3]],
    "N+H1+H2 | H3": [[0, 1, 2], [3]],
}


def fragsim_energy(mol, mf, partition) -> float:
    """Our one-shot embedding, all-electron (n_frozen=0)."""
    basis = build_atom_basis(mol, mf, n_frozen=0)
    d = density_in_basis(basis, mf, n_frozen=0)
    parts = [solve_cluster(mol, mf, basis, d, atoms, n_frozen=0) for atoms in partition]
    return frozen_core_energy(mol, mf, n_frozen=0) + sum(p.e_frag for p in parts)


def vayesta_energy(mf, partition) -> float:
    """Vayesta EWF with atomic IAO+PAO fragments, DMET bath, FCI solver."""
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
        sys.exit("Set VAYESTA_PATH to a Vayesta source checkout (see the module docstring).")
    sys.path.insert(0, path)
    mol = build_mol()
    mf = scf.RHF(mol).run()
    print(f"E_HF = {mf.e_tot:.8f}")
    print(f"{'partition':20s} {'fragsim':>14s} {'Vayesta':>14s} {'diff (kcal/mol)':>16s}")
    for name, partition in PARTITIONS.items():
        ours = fragsim_energy(mol, mf, partition)
        theirs = vayesta_energy(mf, partition)
        print(
            f"{name:20s} {ours:14.8f} {theirs:14.8f} {(ours - theirs) * HARTREE_TO_KCAL:16.3f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
