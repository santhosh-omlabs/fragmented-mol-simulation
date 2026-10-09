"""The committed Komenco run (data/komenco_run.json) backs the README statement."""

import json
from pathlib import Path

import pytest

RECORD = Path(__file__).resolve().parents[1] / "data" / "komenco_run.json"


@pytest.fixture(scope="module")
def record():
    return json.loads(RECORD.read_text(encoding="utf-8"))


def test_circuit_fits_under_the_open_gateway_cap(record):
    assert record["qubits"] == 14
    assert record["qubits"] <= 15


def test_gateway_distribution_equals_our_exact_simulation(record):
    # exact probabilities, no shot noise: the gateway behaves as a classical emulator
    assert record["total_variation_vs_local_exact"] < 1e-9


def test_gateway_ranking_reproduces_local_ideal_ranking(record):
    for row in record["rows"]:
        assert row["gateway_ranked"] == pytest.approx(row["local_ideal_ranked"], abs=1e-6)


def test_sqd_on_gateway_samples_matches_finite_shot_expectation(record):
    # 10,000 shots reach about 11 of the 35 strings; error is close to the k = 12 finite-shot value
    assert 8 <= record["sqd_strings"] <= 14
    assert 3.0 < record["sqd_error_kcal_mol"] < 6.0
