"""Step 4.2: LUCJ circuit for the 14-qubit N fragment (3 empty orbitals kept)."""

import numpy as np
import pytest
from pyscf import fci, scf
from qiskit import transpile

from fragsim.embedding import build_atom_basis, density_in_basis
from fragsim.molecule import build_mol
from fragsim.solver import build_cluster_problem, fragment_energy
from fragsim.sqd_circuit import (
    canonicalize,
    final_state,
    lucj_circuit,
    lucj_operator,
    state_energy,
)

HARTREE_TO_KCAL = 627.5095


@pytest.fixture(scope="module")
def fragment():
    mol = build_mol()
    mf = scf.RHF(mol).run()
    basis = build_atom_basis(mol, mf)
    raw = build_cluster_problem(mol, mf, basis, density_in_basis(basis, mf), [0], n_virtual=3)
    problem, eps = canonicalize(raw)
    e_fci, ci = fci.direct_spin0.FCI(mol).kernel(problem.h1, problem.eri, 7, problem.nelec)
    return raw, problem, eps, e_fci, ci, mol


def test_canonical_orbitals_do_not_change_the_physics(fragment):
    raw, problem, _, e_fci, ci, mol = fragment
    e_raw, ci_raw = fci.direct_spin0.FCI(mol).kernel(raw.h1, raw.eri, 7, raw.nelec)
    assert e_raw == pytest.approx(e_fci, abs=1e-9)
    assert fragment_energy(raw, ci_raw) == pytest.approx(fragment_energy(problem, ci), abs=1e-5)


def test_hartree_fock_determinant_sits_43_kcal_above_exact(fragment):
    import ffsim

    _, problem, _, e_fci, _, _ = fragment
    hf = ffsim.hartree_fock_state(7, problem.nelec).reshape(35, 35)
    assert (state_energy(problem, hf) - e_fci) * HARTREE_TO_KCAL == pytest.approx(43.28, abs=0.05)


@pytest.mark.parametrize("n_reps, expected_kcal", [(1, 13.57), (2, 10.67), (3, 8.45)])
def test_lucj_state_energy_improves_with_layers(fragment, n_reps, expected_kcal):
    _, problem, _, e_fci, ci, _ = fragment
    state = final_state(problem, lucj_operator(problem, n_reps))
    assert np.linalg.norm(state) == pytest.approx(1.0)
    assert abs(np.vdot(ci, state)) ** 2 > 0.989
    energy_error = (state_energy(problem, state) - e_fci) * HARTREE_TO_KCAL
    assert energy_error == pytest.approx(expected_kcal, abs=0.1)


def test_circuit_cost_grows_with_layers(fragment):
    _, problem, *_ = fragment
    cx = []
    for n_reps in (1, 2):
        circuit = lucj_circuit(problem, lucj_operator(problem, n_reps))
        assert circuit.num_qubits == 14
        compiled = transpile(circuit, basis_gates=["rz", "sx", "x", "cx"], seed_transpiler=7)
        cx.append(compiled.count_ops()["cx"])
    assert cx[0] > 300 and cx[1] > cx[0]  # hundreds of two-qubit gates for a 14-qubit fragment
