"""The committed Vayesta cross-check (data/crosscheck_vayesta.json) stays within 0.5 kcal/mol."""

import json
from pathlib import Path

RECORD = Path(__file__).resolve().parents[1] / "data" / "crosscheck_vayesta.json"
HARTREE_TO_KCAL = 627.5095


def test_fragsim_matches_vayesta_within_half_a_kcal():
    rows = json.loads(RECORD.read_text(encoding="utf-8"))["rows"]
    assert len(rows) == 3
    for r in rows:
        diff = (r["fragsim_hartree"] - r["vayesta_hartree"]) * HARTREE_TO_KCAL
        assert abs(diff) < 0.5, r["partition"]
        assert diff > 0  # fragsim recovers slightly less correlation in every case


def test_rows_are_ordered_by_fragment_size():
    rows = json.loads(RECORD.read_text(encoding="utf-8"))["rows"]
    ours = [r["fragsim_hartree"] for r in rows]
    assert ours == sorted(ours, reverse=True)  # bigger fragments: lower energy
