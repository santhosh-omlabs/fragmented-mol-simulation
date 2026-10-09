"""Committed water-dimer cut-error curve (data/water_dimer_cut.json, data/water_sto3g_cut.json)."""

import json
from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parents[1] / "data"


def _rows():
    return json.loads((DATA / "water_dimer_cut.json").read_text(encoding="utf-8"))["rows"]


def test_curve_covers_k_zero_to_nine_with_the_expected_cluster_sizes():
    rows = [r for r in _rows() if r["k"] <= 9]
    assert [r["k"] for r in rows] == list(range(10))
    for r in rows:
        assert (
            r["qubits_per_fragment"] == [12 + 2 * r["k"]] * 2
        )  # 6 occupied cluster orbitals + k empty


def test_primary_error_falls_steadily_as_more_empty_orbitals_are_kept():
    errors = [r["primary_error_kcal_mol"] for r in sorted(_rows(), key=lambda r: r["k"])]
    assert errors == sorted(errors, reverse=True)
    assert (
        errors[0] > 100
    )  # k = 0 is Hartree-Fock for the dimer against exact monomers: not a cut error


def test_primary_and_total_errors_differ_by_the_monomer_method_offset():
    # the offset is 2 x (E_CCSD(T) - E_FCI) of one water, 0.66 kcal/mol, the same at every k
    for r in _rows():
        assert r["primary_error_kcal_mol"] - r[
            "total_energy_error_vs_ccsd_t_kcal_mol"
        ] == pytest.approx(0.658, abs=0.01)


def test_tier_one_is_not_met_at_any_size_up_to_30_qubits():
    for r in _rows():
        if r["k"] <= 9:
            assert r["primary_error_kcal_mol"] > 1.0, r["k"]
    assert (
        min(r["primary_error_kcal_mol"] for r in _rows() if r["k"] <= 9) < 6.0
    )  # but it is getting close to the size of the interaction


def test_the_cut_itself_is_cheap_in_a_minimal_basis():
    rec = json.loads((DATA / "water_sto3g_cut.json").read_text(encoding="utf-8"))
    assert abs(rec["cut_error_kcal_mol"]) < 0.2
    assert rec["e_int_exact_kcal_mol"] == pytest.approx(-4.99, abs=0.01)
    assert rec["qubits_per_cluster"] == [20, 20]


def _error_by_largest_qubits(name):
    rows = json.loads((DATA / name).read_text(encoding="utf-8"))["rows"]
    return {max(r["qubits_per_fragment"]): r["primary_error_kcal_mol"] for r in rows}


def test_moving_the_cut_onto_a_covalent_bond_is_worse_at_equal_largest_cluster_size():
    whole = _error_by_largest_qubits("water_dimer_cut.json")
    shifted = _error_by_largest_qubits("water_dimer_cut_unequal.json")
    common = sorted(set(whole) & set(shifted))
    assert len(common) >= 6  # 20 to 28 qubits at least
    for qubits in common:
        if qubits >= 22:  # below that both are dominated by missing correlation
            assert shifted[qubits] > whole[qubits], qubits


def test_unequal_split_cluster_sizes():
    rows = json.loads((DATA / "water_dimer_cut_unequal.json").read_text(encoding="utf-8"))
    assert rows["partition"] == [[3, 4, 5, 1], [0, 2]]
    for r in rows["rows"]:
        assert r["qubits_per_fragment"] == [14 + 2 * r["k"], 10 + 2 * r["k"]]
