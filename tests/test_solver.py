"""Cluster solver: exactness limits and the headline untruncated fragment error."""

import pytest
from pyscf import scf

from fragsim.embedding import build_atom_basis, density_in_basis
from fragsim.molecule import build_mol, run_reference
from fragsim.solver import frozen_core_energy, solve_cluster

HARTREE_TO_KCAL = 627.5095


@pytest.fixture(scope="module")
def system():
    mol = build_mol()
    mf = scf.RHF(mol).run()
    basis = build_atom_basis(mol, mf)
    return mol, mf, basis, density_in_basis(basis, mf), frozen_core_energy(mol, mf)


@pytest.fixture(scope="module")
def e_fci():
    return run_reference().e_fci


def test_whole_molecule_cluster_equals_fci(system, e_fci):
    mol, mf, basis, d, e_core = system
    result = solve_cluster(mol, mf, basis, d, [0, 1, 2, 3])
    assert e_core + result.e_frag == pytest.approx(e_fci, abs=1e-7)
    assert result.n_qubits == 28


@pytest.mark.slow
def test_fragment_shares_add_to_fci_when_cluster_is_the_whole_space(system, e_fci):
    mol, mf, basis, d, e_core = system
    parts = [solve_cluster(mol, mf, basis, d, f, entire_space=True) for f in ([0, 1], [2, 3])]
    assert e_core + sum(p.e_frag for p in parts) == pytest.approx(e_fci, abs=1e-7)


def test_no_empty_orbitals_gives_hartree_fock(system):
    mol, mf, basis, d, e_core = system
    parts = [solve_cluster(mol, mf, basis, d, [a], n_virtual=0) for a in range(4)]
    assert e_core + sum(p.e_frag for p in parts) == pytest.approx(mf.e_tot, abs=1e-7)


def test_rung1_untruncated_fragmentation_error(system, e_fci):
    mol, mf, basis, d, e_core = system
    parts = [solve_cluster(mol, mf, basis, d, [a]) for a in range(4)]
    assert [p.n_qubits for p in parts] == [22, 6, 6, 6]
    error = (e_core + sum(p.e_frag for p in parts) - e_fci) * HARTREE_TO_KCAL
    assert error == pytest.approx(17.27, abs=0.05)


def test_n_fragment_with_three_empty_orbitals_is_a_14_qubit_problem(system):
    """The Hamiltonian handed to SQD in Lesson 4: 7 orbitals, 4 up + 4 down electrons."""
    import math

    import numpy as np
    from pyscf import fci

    from fragsim.solver import build_cluster_problem, fragment_energy

    mol, mf, basis, d, _ = system
    problem = build_cluster_problem(mol, mf, basis, d, [0], n_virtual=3)
    assert (problem.n_orbitals, problem.n_qubits, problem.nelec) == (7, 14, (4, 4))
    assert math.comb(7, 4) ** 2 == 1225
    assert np.allclose(problem.h1, problem.h1.T)
    _, ci = fci.direct_spin0.FCI(mol).kernel(problem.h1, problem.eri, 7, problem.nelec)
    assert np.shape(ci) == (35, 35)
    assert ci[0, 0] ** 2 == pytest.approx(0.971, abs=0.001)  # Hartree-Fock determinant weight
    assert fragment_energy(problem, ci) == pytest.approx(-14.41254, abs=1e-4)
