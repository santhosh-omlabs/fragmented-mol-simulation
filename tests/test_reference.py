"""Headline facts of the classical baseline. These fail if the README numbers drift."""

import math

import numpy as np
import pytest

from fragsim.molecule import n_determinants, nh3_atoms, run_reference

HARTREE_TO_KCAL = 627.5095


def test_geometry_is_c3v_with_requested_bond_and_angle():
    atoms = nh3_atoms(1.012, 106.7)
    n = np.array(atoms[0][1])
    h = [np.array(a[1]) - n for a in atoms[1:]]
    for v in h:
        assert np.linalg.norm(v) == pytest.approx(1.012)
    for i, j in ((0, 1), (0, 2), (1, 2)):
        cos_a = h[i] @ h[j] / (np.linalg.norm(h[i]) * np.linalg.norm(h[j]))
        assert math.degrees(math.acos(cos_a)) == pytest.approx(106.7, abs=1e-6)


def test_active_space_size():
    assert n_determinants(14, 8) == 1_002_001
    assert 2 * 14 == 28


@pytest.fixture(scope="module")
def ref():
    return run_reference()


def test_whole_molecule_is_28_qubits(ref):
    assert (ref.n_orbitals_total, ref.n_active_orbitals, ref.n_active_electrons) == (15, 14, 8)
    assert ref.n_qubits == 28


def test_energies_6_31g(ref):
    assert ref.e_hf == pytest.approx(-56.161063, abs=1e-5)
    assert ref.e_fci == pytest.approx(-56.291514, abs=1e-5)
    assert ref.correlation_energy * HARTREE_TO_KCAL == pytest.approx(-81.86, abs=0.05)
