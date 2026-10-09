"""The committed water cross-check against Vayesta (data/crosscheck_vayesta_water.json)."""

import json
from pathlib import Path

RECORD = Path(__file__).resolve().parents[1] / "data" / "crosscheck_vayesta_water.json"


def _rows():
    return {r["system"]: r for r in json.loads(RECORD.read_text(encoding="utf-8"))["rows"]}


def test_fragsim_matches_vayesta_to_a_tenth_of_a_kcal_on_water():
    for name, row in _rows().items():
        assert abs(row["difference_kcal_mol"]) < 0.1, name


def test_fragsim_recovers_slightly_less_correlation_as_for_ammonia():
    for name, row in _rows().items():
        assert row["difference_kcal_mol"] > 0, name


def test_one_fragment_per_water():
    rows = _rows()
    assert rows["dimer"]["partition"] == [[0, 1, 2], [3, 4, 5]]
    assert len(rows["trimer ring"]["partition"]) == 3
    assert all(q <= 26 for r in rows.values() for q in r["qubits_per_cluster"])
