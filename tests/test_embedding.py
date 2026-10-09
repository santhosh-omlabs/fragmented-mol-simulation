"""Step 2 of embedding: atom basis, fragment/environment split, bath and cluster sizes."""

import numpy as np
import pytest
from pyscf import scf

from fragsim.embedding import build_atom_basis, density_in_basis, make_cluster
from fragsim.molecule import build_mol


@pytest.fixture(scope="module")
def system():
    mol = build_mol()
    mf = scf.RHF(mol).run()
    basis = build_atom_basis(mol, mf)
    return mol, mf, basis, density_in_basis(basis, mf)


def test_basis_is_complete_orthonormal_and_atom_tagged(system):
    mol, mf, basis, _ = system
    s = mf.get_ovlp()
    assert np.allclose(basis.coeff.T @ s @ basis.coeff, np.eye(mol.nao - 1), atol=1e-8)
    assert basis.n_iao == 7  # 8 IAOs minus the one that is the frozen N 1s core
    # N owns 8 functions (4 IAO + 4 virtual), each H owns 2 (1 IAO + 1 virtual)
    assert np.bincount(basis.atom).tolist() == [8, 2, 2, 2]
    assert np.abs(mf.mo_coeff[:, :1].T @ s @ basis.coeff).max() < 1e-8  # orthogonal to the core


def test_density_is_a_projector_with_eight_electrons(system):
    d = system[3]
    assert np.allclose(d @ d, d, atol=1e-8)
    assert 2 * np.trace(d) == pytest.approx(8.0)


@pytest.mark.parametrize("atoms", [[0], [1], [0, 1]])
def test_bath_strength_is_x_times_one_minus_x(system, atoms):
    """sigma^2 = x (1 - x), x = eigenvalue of the fragment block of D. So sigma <= 0.5."""
    _, _, basis, d = system
    cluster = make_cluster(basis, d, atoms)
    frag = cluster.fragment
    x = np.linalg.eigvalsh(d[np.ix_(frag, frag)])
    expected = np.sort(x * (1 - x))[::-1][: len(cluster.singular_values)]
    assert np.allclose(np.sort(cluster.singular_values**2)[::-1], expected, atol=1e-8)


@pytest.mark.parametrize(
    "atoms, n_frag, n_bath, n_orb, n_elec",
    [
        ([0], 8, 3, 11, 8),  # N fragment: lone pair + 3 bonds seen through the bath (22 qubits)
        ([1], 2, 1, 3, 2),  # one H: its bond partner on N is the single bath orbital
        ([0, 1], 10, 2, 12, 8),
        ([0, 1, 2], 12, 1, 13, 8),
    ],
)
def test_cluster_sizes(system, atoms, n_frag, n_bath, n_orb, n_elec):
    cluster = make_cluster(system[2], system[3], atoms)
    assert len(cluster.fragment) == n_frag
    assert cluster.bath.shape[1] == n_bath
    assert cluster.n_orbitals == n_orb
    assert cluster.n_electrons == n_elec
    assert np.all(cluster.singular_values <= 0.5 + 1e-9)  # one-spin bound
