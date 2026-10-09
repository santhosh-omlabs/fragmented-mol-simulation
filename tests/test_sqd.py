"""Steps 4.4-4.5: subspace diagonalization, configuration recovery, and the random control."""

import itertools

import numpy as np
import pytest
from pyscf import fci, scf
from qiskit_addon_sqd.counts import counts_to_arrays

from fragsim.embedding import build_atom_basis, density_in_basis
from fragsim.molecule import build_mol
from fragsim.sampling import sample_ideal, sample_noisy
from fragsim.solver import build_cluster_problem, fragment_energy
from fragsim.sqd import run_sqd, solve_subspace, uniform_random_counts
from fragsim.sqd_circuit import canonicalize, final_state, lucj_circuit, lucj_operator

HARTREE_TO_KCAL = 627.5095


@pytest.fixture(scope="module")
def fragment():
    mol = build_mol()
    mf = scf.RHF(mol).run()
    basis = build_atom_basis(mol, mf)
    raw = build_cluster_problem(mol, mf, basis, density_in_basis(basis, mf), [0], n_virtual=3)
    problem, _ = canonicalize(raw)
    e_fci, ci = fci.direct_spin0.FCI(mol).kernel(problem.h1, problem.eri, 7, problem.nelec)
    return problem, e_fci, ci


def _error_kcal(result, e_fci):
    return (result.energy - e_fci) * HARTREE_TO_KCAL


def test_subspace_with_every_string_reproduces_fci_and_the_fragment_share(fragment):
    problem, e_fci, ci = fragment
    patterns = [
        "".join("1" if i in c else "0" for i in range(6, -1, -1))
        for c in itertools.combinations(range(7), 4)
    ]
    matrix = np.array([[ch == "1" for ch in down + up] for down in patterns for up in patterns])
    result = solve_subspace(problem, matrix)
    assert result.dimension == 1225
    assert result.energy == pytest.approx(e_fci, abs=1e-9)
    assert fragment_energy(problem, result.ci) == pytest.approx(
        fragment_energy(problem, ci), abs=1e-5
    )


def test_ideal_circuit_samples_give_a_few_kcal_error(fragment):
    problem, e_fci, _ = fragment
    state = final_state(problem, lucj_operator(problem, 2))
    matrix, _ = counts_to_arrays(sample_ideal(state, problem, 1000, seed=11))
    result = solve_subspace(problem, matrix)
    assert result.dimension < 100  # of 1225
    assert _error_kcal(result, e_fci) == pytest.approx(4.58, abs=0.05)


def test_random_bitstrings_make_this_fragment_trivial(fragment):
    """100 random right-count bitstrings already contain all 35 up strings: exact answer."""
    problem, e_fci, _ = fragment
    matrix, _ = counts_to_arrays(uniform_random_counts(problem, 100, seed=3))
    result = solve_subspace(problem, matrix)
    assert result.dimension == 1225
    assert abs(_error_kcal(result, e_fci)) < 1e-6


@pytest.mark.slow
def test_recovery_loop_runs_on_noisy_samples(fragment):
    problem, e_fci, _ = fragment
    circuit = lucj_circuit(problem, lucj_operator(problem, 2))
    counts = sample_noisy(circuit, shots=200, p_cx=0.005, p_1q=0.0005, seed=5)
    best, history = run_sqd(problem, counts, iterations=3, samples_per_batch=100, num_batches=2)
    assert [h.index for h in history] == [0, 1, 2]
    assert history[0].n_repaired > 0.5 * len(counts)  # most noisy bitstrings have a wrong count
    assert _error_kcal(best, e_fci) >= -1e-6  # a subspace can never beat the exact energy
    assert _error_kcal(best, e_fci) < 10.0
