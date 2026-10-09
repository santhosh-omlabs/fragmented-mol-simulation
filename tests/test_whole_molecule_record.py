"""Committed whole-molecule IBM run (data/whole_molecule*.json)."""

import json
from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parents[1] / "data"


@pytest.fixture(scope="module")
def rec():
    return json.loads((DATA / "whole_molecule.json").read_text(encoding="utf-8"))


def test_hardware_samples_carry_signal_over_random_bits_in_electron_count(rec):
    assert (
        rec["right_count_fraction"] > 5 * rec["controls_at_k"]["random_bits_right_count_fraction"]
    )


def test_hardware_sqd_is_barely_better_than_random_bits_through_the_same_loop(rec):
    hw = rec["sqd_loop"]["energy_error_kcal_mol"]
    rnd = rec["controls_at_k"]["random_bits_sqd_loop"]
    assert hw < rnd and rnd - hw < 2.0


def test_free_classical_ordering_beats_the_hardware_at_equal_subspace_size(rec):
    c = rec["controls_at_k"]
    assert c["excitation_order"] < rec["sqd_loop"]["energy_error_kcal_mol"] / 2
    assert c["ccsd_ranked"] <= c["excitation_order"] + 1e-6


def test_ideal_whole_molecule_circuit_is_worse_than_classical_orderings_at_large_k():
    path = DATA / "whole_molecule_dry.json"
    dry = json.loads(path.read_text(encoding="utf-8"))["ideal_matched_error_kcal_mol"]
    for k in ("64", "128", "256"):
        assert (
            dry["ideal heavy-hex circuit"][k] > dry["excitation order"][k] > dry["CCSD-ranked"][k]
        )
