"""Lesson 7.4 side check: the UNTRUNCATED cut error of the water dimer in STO-3G, against exact FCI.

STO-3G is small enough that the whole dimer (12 active orbitals, frozen O 1s) and every fragment cluster can be solved
exactly, so the cut error is measured with no truncation. It says what a one-shot DMET-style cut costs when there are no
empty functions on the neighbouring water. It does not transfer to 6-31G.

Run:  python scripts/water_sto3g_cut.py   (about 1 minute, writes data/water_sto3g_cut.json)
"""

from __future__ import annotations

import json
from pathlib import Path

from pyscf import mcscf, scf
from pyscf.fci import direct_spin0

from fragsim.embedding import build_atom_basis, density_in_basis, make_cluster
from fragsim.molecule import build_mol
from fragsim.solver import frozen_core_energy, solve_cluster
from fragsim.water import HARTREE_TO_KCAL, dimer_atoms, flatten, monomer_atoms, water_fragments

RECORD = Path(__file__).resolve().parents[1] / "data" / "water_sto3g_cut.json"


def fci_energy(atoms, n_core: int):
    mol = build_mol(atoms, "sto-3g")
    mf = scf.RHF(mol).run()
    cas = mcscf.CASCI(mf, mol.nao - n_core, mol.nelectron - 2 * n_core)
    cas.fcisolver = direct_spin0.FCI(mol)
    return float(cas.kernel()[0]), mol, mf


def main() -> None:
    e_dimer, mol, mf = fci_energy(flatten(dimer_atoms()), 2)
    e_water, _, _ = fci_energy(monomer_atoms(), 1)
    basis = build_atom_basis(mol, mf, n_frozen=2)
    d = density_in_basis(basis, mf, n_frozen=2)
    qubits = [2 * make_cluster(basis, d, a).n_orbitals for a in water_fragments(2)]
    parts = [solve_cluster(mol, mf, basis, d, a, n_frozen=2) for a in water_fragments(2)]
    e_frag = frozen_core_energy(mol, mf, n_frozen=2) + sum(p.e_frag for p in parts)
    out = {
        "what": "water dimer, STO-3G, frozen O 1s: exact FCI against two-water fragments with exact clusters, no truncation",
        "qubits_per_cluster": qubits,
        "e_int_exact_kcal_mol": (e_dimer - 2 * e_water) * HARTREE_TO_KCAL,
        "e_int_fragments_kcal_mol": (e_frag - 2 * e_water) * HARTREE_TO_KCAL,
        "cut_error_kcal_mol": (e_frag - e_dimer) * HARTREE_TO_KCAL,
    }
    RECORD.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
