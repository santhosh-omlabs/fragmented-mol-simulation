"""The committed IBM run (data/ibm_*.json) backs the README statement. No network."""

import json
from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parents[1] / "data"


@pytest.fixture(scope="module")
def run():
    return json.loads((DATA / "ibm_run.json").read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def src():
    return json.loads((DATA / "ibm_analysis.json").read_text(encoding="utf-8"))["sources"]


def test_run_record_is_a_single_4000_shot_job(run):
    assert (
        run["shots"] == 4000 and run["backend"] == "ibm_fez" and run["quantum_seconds_billed"] <= 10
    )


def test_hardware_is_close_to_random_bits_in_electron_count(src):
    # random 14-bit strings have (35/128)^2 = 7.5% chance of the right count in both halves
    assert src["random_bits"]["fraction_right_electron_count"] == pytest.approx(0.075, abs=0.01)
    assert src["ibm_fez"]["fraction_right_electron_count"] < 0.10
    assert (
        src["ibm_fez"]["fraction_right_electron_count"]
        < src["noisy_simulation_0.5pct"]["fraction_right_electron_count"]
    )


def test_hardware_keeps_some_signal_at_small_subspaces_and_loses_it_at_large(src):
    hw, rnd = src["ibm_fez"]["matched"], src["random_bits"]["matched"]
    assert hw["4"] < rnd["4"] - 10 and hw["8"] < rnd["8"] - 8
    assert abs(hw["12"] - rnd["12"]) < 1.5  # indistinguishable from random bits plus recovery
