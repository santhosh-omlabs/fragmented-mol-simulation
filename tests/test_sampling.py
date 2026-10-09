"""Step 4.3: bitstring convention and sampling statistics for the 14-qubit N fragment."""

import ffsim
import pytest
from pyscf import fci, scf

from fragsim.embedding import build_atom_basis, density_in_basis
from fragsim.molecule import build_mol
from fragsim.sampling import (
    determinant_index,
    has_correct_electron_count,
    sample_ideal,
    sample_noisy,
    sample_statistics,
    split_spins,
)
from fragsim.solver import build_cluster_problem
from fragsim.sqd_circuit import canonicalize, final_state, lucj_circuit, lucj_operator

HF_BITSTRING = "0001111" * 2  # [beta][alpha], orbitals 0-3 occupied, orbital 0 = rightmost


@pytest.fixture(scope="module")
def fragment():
    mol = build_mol()
    mf = scf.RHF(mol).run()
    basis = build_atom_basis(mol, mf)
    raw = build_cluster_problem(mol, mf, basis, density_in_basis(basis, mf), [0], n_virtual=3)
    problem, _ = canonicalize(raw)
    _, ci = fci.direct_spin0.FCI(mol).kernel(problem.h1, problem.eri, 7, problem.nelec)
    return problem, ci


def test_bitstring_convention():
    assert split_spins("0000001" + "0001111", 7) == ("0001111", "0000001")
    assert has_correct_electron_count(HF_BITSTRING, 7, (4, 4))
    assert not has_correct_electron_count("0001111" + "0000111", 7, (4, 4))


def test_hartree_fock_bitstring_is_the_first_determinant(fragment):
    problem, _ = fragment
    state = ffsim.hartree_fock_state(7, problem.nelec).reshape(35, 35)
    counts = sample_ideal(state, problem, shots=5, seed=1)
    assert counts == {HF_BITSTRING: 5}
    assert determinant_index(HF_BITSTRING, 7, problem.nelec) == (0, 0)


def test_ideal_sampling_stays_in_the_right_sector_and_covers_the_ground_state(fragment):
    problem, ci = fragment
    state = final_state(problem, lucj_operator(problem, 2))
    stats = sample_statistics(sample_ideal(state, problem, 1000, seed=11), problem, ci)
    assert stats["fraction_correct_electron_count"] == 1.0
    assert 40 <= stats["subspace_size"] <= 80  # of 1225 determinants
    assert stats["exact_weight_in_subspace"] > 0.99


@pytest.mark.slow
def test_noise_pushes_samples_out_of_the_right_electron_sector(fragment):
    problem, ci = fragment
    circuit = lucj_circuit(problem, lucj_operator(problem, 2))
    counts = sample_noisy(circuit, shots=100, p_cx=0.005, p_1q=0.0005, seed=5)
    stats = sample_statistics(counts, problem, ci)
    assert stats["fraction_correct_electron_count"] < 0.5
    assert counts.get(HF_BITSTRING, 0) / 100 < 0.3  # ideal circuit puts ~97% here
