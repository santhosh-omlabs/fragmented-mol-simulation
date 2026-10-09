"""The committed scheduling study (data/schedule_study.json) backs the README table."""

import json
from pathlib import Path

import pytest

RECORD = Path(__file__).resolve().parents[1] / "data" / "schedule_study.json"


@pytest.fixture(scope="module")
def rows():
    return {r["config"]: r for r in json.loads(RECORD.read_text(encoding="utf-8"))["rows"]}


def _shots(row):
    return max(j["shots_needed"] for j in row["jobs"])


def test_fragmenting_cuts_required_shots_by_orders_of_magnitude(rows):
    assert _shots(rows["whole molecule"]) > 1e6
    assert _shots(rows["rung 1, 3 empty kept"]) < 5e3
    assert _shots(rows["whole molecule"]) > 1000 * _shots(rows["rung 1, 3 empty kept"])


def test_whole_molecule_does_not_fit_a_ten_minute_budget_on_one_qpu(rows):
    assert rows["whole molecule"]["total_qpu_seconds"] > 600
    assert rows["rung 1, 3 empty kept"]["total_qpu_seconds"] < 60


def test_fragment_parallel_speedup_is_capped_by_the_largest_fragment(rows):
    row = rows["rung 1, all empty kept"]
    biggest = max(j["seconds"] for j in row["jobs"])
    assert row["wall_seconds"]["4"] == pytest.approx(biggest)
    assert row["wall_seconds"]["1"] / row["wall_seconds"]["4"] < 1.2


def test_shot_parallel_speedup_is_below_linear(rows):
    sp = rows["whole molecule"]["shot_parallel_wall_seconds"]
    assert 3.5 < sp["1"] / sp["4"] < 4.0


def test_bigger_fragments_cost_more_but_err_less(rows):
    order = ["rung 1, 3 empty kept", "rung 1, all empty kept", "rung 3, all empty kept"]
    errors = [rows[k]["error_kcal_mol"] for k in order]
    costs = [rows[k]["total_qpu_seconds"] for k in order]
    assert errors == sorted(errors, reverse=True) and costs == sorted(costs)
