"""Water cluster geometries and the committed CCSD(T) reference (Lesson 7.1)."""

import json
import math
from pathlib import Path

import numpy as np
import pytest

from fragsim.molecule import n_determinants
from fragsim.water import (
    ANGLE_HOH,
    R_OH,
    dimer_atoms,
    flatten,
    min_intermolecular_distance,
    ring_atoms,
    sizes,
)

RECORD = Path(__file__).resolve().parents[1] / "data" / "water_reference.json"
GEOMETRIES = {"dimer": dimer_atoms(), "trimer ring": ring_atoms(3), "tetramer ring": ring_atoms(4)}


@pytest.mark.parametrize("name", list(GEOMETRIES))
def test_every_water_keeps_the_monomer_geometry(name):
    for water in GEOMETRIES[name]:
        o, h1, h2 = (np.array(c) for _, c in water)
        assert np.linalg.norm(h1 - o) == pytest.approx(R_OH, abs=1e-9)
        assert np.linalg.norm(h2 - o) == pytest.approx(R_OH, abs=1e-9)
        cos = np.dot(h1 - o, h2 - o) / R_OH**2
        assert math.acos(cos) == pytest.approx(ANGLE_HOH, abs=1e-9)


@pytest.mark.parametrize("name", list(GEOMETRIES))
def test_no_atoms_from_different_waters_overlap(name):
    # a hydrogen bond puts H and O 1.9-2.0 A apart; anything under 1.5 A would be a clash
    assert min_intermolecular_distance(GEOMETRIES[name]) > 1.5


def test_sizes_match_the_design_note():
    one, two = sizes(1), sizes(2)
    assert (one["active_orbitals"], one["active_electrons"], one["qubits"]) == (12, 8, 24)
    assert (two["active_orbitals"], two["active_electrons"], two["qubits"]) == (24, 16, 48)
    assert 10 ** one["log10_determinants"] == pytest.approx(n_determinants(12, 8), rel=1e-9)
    assert n_determinants(12, 8) == 245025


def test_committed_geometry_is_what_the_code_builds():
    rec = json.loads(RECORD.read_text(encoding="utf-8"))["systems"]
    for name, waters in GEOMETRIES.items():
        built = [(s, list(c)) for s, c in flatten(waters)]
        stored = [(s, list(c)) for s, c in rec[name]["atoms"]]
        assert [s for s, _ in built] == [s for s, _ in stored]
        assert np.allclose([c for _, c in built], [c for _, c in stored], atol=1e-9)


def test_dimer_interaction_energy_reference():
    e = json.loads(RECORD.read_text(encoding="utf-8"))["systems"]["dimer"]["interaction_kcal_mol"]
    assert e["ccsd_t"]["raw"] == pytest.approx(-6.59, abs=0.02)
    assert e["ccsd_t"]["counterpoise"] == pytest.approx(-4.62, abs=0.02)
    assert e["ccsd_t"]["bsse"] < -1.5  # superposition error is ~30% of the binding in 6-31G


def test_mean_field_carries_most_of_every_interaction():
    systems = json.loads(RECORD.read_text(encoding="utf-8"))["systems"]
    for name, system in systems.items():
        e = system["interaction_kcal_mol"]
        assert e["ccsd_t"]["raw"] < e["hf"]["raw"] < 0, name  # bound, correlation adds binding
        assert e["hf"]["raw"] / e["ccsd_t"]["raw"] > 0.80, name  # observed 83% to 93%


def test_rings_are_bound_and_binding_grows_with_size():
    systems = json.loads(RECORD.read_text(encoding="utf-8"))["systems"]
    raw = {n: s["interaction_kcal_mol"]["ccsd_t"]["raw"] for n, s in systems.items()}
    assert raw["tetramer ring"] < raw["trimer ring"] < raw["dimer"] < 0
    for s in systems.values():
        assert s["interaction_kcal_mol"]["ccsd_t"]["counterpoise"] < 0


def test_ccsd_t_is_accurate_for_the_interaction_energy_where_fci_is_possible():
    rec = json.loads(RECORD.read_text(encoding="utf-8"))
    sto = rec["sto3g_reference_check"]
    assert abs(sto["ccsd_t_minus_fci_kcal_mol"]) < 0.05
    # monomer FCI sits below CCSD(T) in 6-31G by a third of a kcal/mol: a constant that must cancel, not be subtracted
    assert -0.45 < rec["monomer_fci_check"]["fci_minus_ccsd_t_kcal_mol"] < -0.2
