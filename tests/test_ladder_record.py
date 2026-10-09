"""The committed ladder record (data/ladder_exact_fragments.json) backs the README table."""

import json
from pathlib import Path

import pytest

RECORD = Path(__file__).resolve().parents[1] / "data" / "ladder_exact_fragments.json"
HARTREE_TO_KCAL = 627.5095


@pytest.fixture(scope="module")
def rows():
    return json.loads(RECORD.read_text(encoding="utf-8"))["rows"]


def _rung(rows, prefix):
    return [r for r in rows if r["rung"].startswith(prefix)]


def test_no_virtuals_equals_hartree_fock(rows):
    record = json.loads(RECORD.read_text(encoding="utf-8"))
    expected = (record["e_hf"] - record["e_fci"]) * HARTREE_TO_KCAL
    for r in (r for r in rows if r["virtuals_kept"] == 0):
        assert r["error_kcal_mol"] == pytest.approx(expected, abs=1e-3)


@pytest.mark.parametrize("prefix", ["1:", "2:", "3:"])
def test_more_virtuals_never_increase_the_error(rows, prefix):
    errors = [r["error_kcal_mol"] for r in _rung(rows, prefix)]  # k = 0..5, then all
    assert errors == sorted(errors, reverse=True)
    assert all(e > 0 for e in errors)  # fragments recover less correlation than exact


def test_untruncated_fragmentation_error_by_rung(rows):
    full = {r["rung"][:2]: r for r in rows if r["virtuals_kept"] is None}
    assert full["1:"]["error_kcal_mol"] == pytest.approx(17.27, abs=0.05)
    assert full["2:"]["error_kcal_mol"] == pytest.approx(11.69, abs=0.05)
    assert full["3:"]["error_kcal_mol"] == pytest.approx(5.94, abs=0.05)
    assert max(full["1:"]["qubits_per_fragment"]) == 22


def test_largest_fragment_under_komenco_cap_is_far_from_chemical_accuracy(rows):
    fits = [r for r in _rung(rows, "1:") if max(r["qubits_per_fragment"]) <= 15]
    best = min(r["error_kcal_mol"] for r in fits)
    assert best == pytest.approx(38.82, abs=0.05)  # 14 qubits, 3 empty orbitals kept
