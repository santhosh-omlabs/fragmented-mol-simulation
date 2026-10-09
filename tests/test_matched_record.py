"""The committed matched-size study (data/matched_subspace_study.json) backs the README table."""

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from fragsim.matched import all_strings, excitation_order

RECORD = Path(__file__).resolve().parents[1] / "data" / "matched_subspace_study.json"


@pytest.fixture(scope="module")
def by_k():
    rows = json.loads(RECORD.read_text(encoding="utf-8"))["rows"]
    return {r["k"]: r for r in rows}


def test_every_method_is_exact_when_all_35_strings_are_used(by_k):
    row = by_k[35]
    for key in ("oracle", "circuit_infinite_shots", "ccsd_ranked_classical", "excitation_ordered"):
        assert abs(row[key]) < 1e-6
    assert abs(row["random_median"]) < 1e-6


def test_random_strings_are_far_worse_than_any_informed_choice(by_k):
    assert (
        by_k[8]["random_median"] > 500
    )  # without Hartree-Fock the energy is thousands of kcal/mol off
    assert by_k[12]["random_with_hf_median"] > 25
    assert by_k[12]["random_with_hf_median"] > 5 * by_k[12]["circuit_infinite_shots"]


def test_classical_ccsd_ranking_matches_the_best_possible_choice(by_k):
    for k in (1, 2, 3, 4, 5, 6, 8, 12, 16):
        assert by_k[k]["ccsd_ranked_classical"] == pytest.approx(by_k[k]["oracle"], abs=0.01)


def test_ideal_circuit_beats_simple_excitation_ordering_at_small_k(by_k):
    for k in (4, 5, 6, 8, 10, 12, 16):
        assert by_k[k]["circuit_infinite_shots"] < by_k[k]["excitation_ordered"]


def test_noise_removes_the_circuit_advantage(by_k):
    for k in (4, 5, 6, 8, 10, 12, 16):
        assert by_k[k]["noisy_circuit_recovered"] > by_k[k]["circuit_infinite_shots"]
    # with noise the circuit is no better than the classical excitation ordering at k = 8..12
    for k in (8, 10, 12):
        gap = by_k[k]["noisy_circuit_recovered"] - by_k[k]["excitation_ordered"]
        assert abs(gap) < 1.5


def test_finite_shots_cannot_reach_large_subspaces(by_k):
    assert 300 < by_k[8]["circuit_shots_to_reach_k_median"] < 600
    assert by_k[12]["circuit_shots_to_reach_k_median"] > 5000
    assert by_k[16]["circuit_finite_seeds_reaching_k"] == 0  # not reached in 20,000 shots


def test_excitation_order_puts_hartree_fock_first_then_single_excitations():
    problem = SimpleNamespace(n_orbitals=7, n_occ=4)
    order = excitation_order(problem, np.arange(7, dtype=float))
    assert len(order) == len(all_strings(problem)) == 35
    assert order[0] == 0b0001111  # Hartree-Fock pattern
    hf = 0b0001111
    assert all((s & ~hf).bit_count() == 1 for s in order[1:13])  # 4 x 3 single excitations
    assert all((s & ~hf).bit_count() >= 2 for s in order[13:])
