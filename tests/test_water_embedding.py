"""Lesson 7.2: the embedding pipeline on water (two equivalent oxygen cores, no C3v symmetry)."""

import numpy as np
import pytest
from pyscf import scf

from fragsim.embedding import build_atom_basis, density_in_basis
from fragsim.molecule import BASIS, build_mol, run_reference
from fragsim.solver import frozen_core_energy, solve_cluster
from fragsim.water import dimer_atoms, flatten, monomer_atoms, water_fragments


@pytest.fixture(scope="module")
def dimer():
    mol = build_mol(flatten(dimer_atoms()), BASIS)
    mf = scf.RHF(mol).run()
    basis = build_atom_basis(mol, mf, n_frozen=2)
    return mol, mf, basis, density_in_basis(basis, mf, n_frozen=2)


def test_water_fragments_group_atoms_by_molecule():
    assert water_fragments(1) == [[0, 1, 2]]
    assert water_fragments(3) == [[0, 1, 2], [3, 4, 5], [6, 7, 8]]


def test_dimer_basis_is_orthonormal_and_gives_each_water_twelve_functions(dimer):
    mol, mf, basis, _ = dimer
    c, s = basis.coeff, mf.get_ovlp()
    assert c.shape[1] == mol.nao - 2 == 24
    assert np.allclose(c.T @ s @ c, np.eye(24), atol=1e-8)
    for water in water_fragments(2):
        assert int(np.isin(basis.atom, water).sum()) == 12


def test_both_oxygen_core_functions_are_removed_not_one_twice(dimer):
    _, _, basis, _ = dimer
    # each oxygen keeps 8 of its 9 functions (the 1s core is gone); a double removal on one atom would give 7 and 9
    for oxygen in (0, 3):
        assert int((basis.atom == oxygen).sum()) == 8


def test_one_water_as_a_single_fragment_reproduces_fci():
    mol = build_mol(monomer_atoms(), BASIS)
    mf = scf.RHF(mol).run()
    basis = build_atom_basis(mol, mf)
    d = density_in_basis(basis, mf)
    part = solve_cluster(mol, mf, basis, d, [0, 1, 2])
    total = frozen_core_energy(mol, mf) + part.e_frag
    # agreement to 1e-7 hartree (6e-5 kcal/mol): the two FCI paths stop at slightly different convergence thresholds
    assert total == pytest.approx(run_reference(mol).e_fci, abs=1e-7)


def test_dimer_with_no_empty_orbitals_reproduces_hartree_fock(dimer):
    mol, mf, basis, d = dimer
    e_core = frozen_core_energy(mol, mf, n_frozen=2)
    parts = [
        solve_cluster(mol, mf, basis, d, atoms, n_virtual=0, n_frozen=2)
        for atoms in water_fragments(2)
    ]
    assert e_core + sum(p.e_frag for p in parts) == pytest.approx(mf.e_tot, abs=1e-8)


def test_dimer_cluster_with_few_empty_orbitals_is_a_proper_correlated_energy(dimer):
    mol, mf, basis, d = dimer
    e_core = frozen_core_energy(mol, mf, n_frozen=2)
    parts = [
        solve_cluster(mol, mf, basis, d, atoms, n_virtual=2, n_frozen=2)
        for atoms in water_fragments(2)
    ]
    energy = e_core + sum(p.e_frag for p in parts)
    assert energy < mf.e_tot - 0.01  # recovers correlation energy
    assert energy > mf.e_tot - 0.60  # and nothing absurd (dimer correlation is about 0.5 hartree)
    assert all(p.n_qubits <= 24 for p in parts)


def test_atom_basis_does_not_depend_on_how_the_two_core_orbitals_are_mixed(dimer):
    """A rotation inside the frozen-core pair is an equally valid SCF output; the basis must not change."""
    import copy

    mol, mf, basis, _ = dimer
    t = np.radians(45.0)
    rotated = copy.copy(mf)
    rotated.mo_coeff = mf.mo_coeff.copy()
    rotated.mo_coeff[:, :2] = mf.mo_coeff[:, :2] @ np.array(
        [[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]]
    )
    other = build_atom_basis(mol, rotated, n_frozen=2)
    s = mf.get_ovlp()
    assert [int((other.atom == a).sum()) for a in range(6)] == [
        int((basis.atom == a).sum()) for a in range(6)
    ]
    p_a = basis.coeff @ basis.coeff.T @ s
    p_b = other.coeff @ other.coeff.T @ s
    assert np.allclose(p_a, p_b, atol=1e-7)  # same 24-dimensional active space
