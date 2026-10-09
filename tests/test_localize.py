"""Step 1 of embedding: localized occupied orbitals are core, lone pair and three N-H bonds."""

import numpy as np
import pytest
from pyscf import scf

from fragsim.localize import atom_weights, localized_occupied
from fragsim.molecule import build_mol


@pytest.fixture(scope="module")
def setup():
    mol = build_mol()
    mf = scf.RHF(mol).run()
    return mol, mf, localized_occupied(mf)


def test_rotation_keeps_the_same_occupied_space(setup):
    _, mf, loc = setup
    occ = mf.mo_coeff[:, mf.mo_occ > 0]
    assert np.allclose(occ @ occ.T, loc @ loc.T, atol=1e-10)
    assert np.allclose(loc.T @ mf.get_ovlp() @ loc, np.eye(5), atol=1e-8)


def test_orbitals_are_core_lone_pair_and_three_bonds(setup):
    mol, mf, loc = setup
    w = atom_weights(mol, mf, loc)
    on_n_only = np.sum(w[:, 0] > 0.99)
    bonds = [row for row in w if row[0] < 0.99]
    assert on_n_only == 2  # N 1s core + lone pair
    assert len(bonds) == 3
    for row in bonds:
        assert row[0] == pytest.approx(0.622, abs=0.01)  # N share of an N-H bond
        assert np.sum(row[1:] > 0.3) == 1  # exactly one H carries the rest
