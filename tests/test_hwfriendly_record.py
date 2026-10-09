"""Committed heavy-hex study and IBM gate-count scaling (data/hwfriendly_study.json, ibm_scaling_analysis.json)."""

import json
from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parents[1] / "data"


@pytest.fixture(scope="module")
def study():
    return json.loads((DATA / "hwfriendly_study.json").read_text(encoding="utf-8"))["variants"]


@pytest.fixture(scope="module")
def scaling():
    return json.loads((DATA / "ibm_scaling_analysis.json").read_text(encoding="utf-8"))


def test_heavy_hex_pairs_cut_routed_gates_several_fold(study):
    assert study["heavy-hex pairs, 1 layer(s)"]["best_cz"] < 200
    assert (
        study["all-to-all, 2 layer(s)"]["best_cz"]
        > 4 * study["heavy-hex pairs, 1 layer(s)"]["best_cz"]
    )


def test_heavy_hex_keeps_most_of_the_ideal_quality_at_large_k(study):
    hh = study["heavy-hex pairs, 1 layer(s)"]["ideal_matched_error_kcal_mol"]
    a2a = study["all-to-all, 1 layer(s)"]["ideal_matched_error_kcal_mol"]
    assert hh["12"] == pytest.approx(a2a["12"], abs=0.1)
    assert hh["8"] > a2a["8"]  # restriction costs something at small k


def test_hardware_sample_quality_falls_monotonically_with_gate_count(scaling):
    rows = scaling["rows"]
    heavy = scaling["heavy_hex_summary"]["heavy_hex_right_count_mean"]
    mid = rows["all-to-all 1 layer"]["right_count_fraction"]
    big = rows["all-to-all 2 layers (first job)"]["right_count_fraction"]
    assert heavy > mid > big
    assert heavy > 0.3 and big < 0.1


def test_better_samples_did_not_give_better_subspaces_than_the_classical_ordering(scaling):
    s = scaling["heavy_hex_summary"]
    excitation = {8: 10.5, 12: 8.4}  # data/matched_subspace_study.json
    assert s["heavy_hex_k8_mean"] > excitation[8] - 1.0
    assert s["heavy_hex_k12_mean"] > excitation[12] - 1.0
