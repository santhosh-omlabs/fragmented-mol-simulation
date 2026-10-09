"""Lesson 7.1: water clusters, model geometries and the CCSD(T) reference.

Geometries are MODEL structures, not optimised: every water has the experimental monomer geometry
(O-H 0.9572 A, H-O-H 104.52 deg) and O...O = 2.91 A along each hydrogen bond. Only fragment-versus-
reference differences at the same geometry matter in Lesson 7.

  dimer   linear O-H...O, acceptor bisector tilted 55 degrees from the O...O axis
  rings   n waters on a regular polygon (n = 3, 4); each donates one H straight to the next oxygen,
          its free H points up or down alternately (an odd ring has one repeat; it is a model)

Interaction energies (kcal/mol, negative = bound), all in the 6-31G basis with frozen oxygen 1s cores:

  raw  E_int = E(cluster) - sum E(water in its own basis)
  cp   same, with every water computed in the whole cluster's basis (counterpoise, for context only)

Run:  python -m fragsim.water fit   (about 3 minutes; prints RING_PARAMS to paste above)
      python -m fragsim.water       (about 2 minutes, writes data/water_reference.json)
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
from pyscf import cc, gto, scf

from fragsim.molecule import BASIS, build_mol, run_reference

HARTREE_TO_KCAL = 627.5095
R_OH = 0.9572
ANGLE_HOH = math.radians(104.52)
R_OO = 2.91
DATA = Path(__file__).resolve().parents[2] / "data"
RECORD = DATA / "water_reference.json"

Atom = tuple[str, tuple[float, float, float]]
Monomer = list[Atom]


def _water(o: np.ndarray, h_bonded: np.ndarray, h_free: np.ndarray) -> Monomer:
    return [("O", tuple(o)), ("H", tuple(o + R_OH * h_bonded)), ("H", tuple(o + R_OH * h_free))]


def monomer_atoms() -> Monomer:
    """One water in its own frame."""
    z = np.array([0.0, 0.0, 1.0])
    free = np.array([math.sin(ANGLE_HOH), 0.0, math.cos(ANGLE_HOH)])
    return _water(np.zeros(3), z, free)


def dimer_atoms(r_oo: float = R_OO) -> list[Monomer]:
    """Donor along +z from the origin; acceptor oxygen at (0, 0, r_oo)."""
    z, x = np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0])
    donor = _water(np.zeros(3), z, np.array([math.sin(ANGLE_HOH), 0.0, math.cos(ANGLE_HOH)]))
    o2 = r_oo * z
    bisector = np.array([0.0, math.sin(math.radians(55)), math.cos(math.radians(55))])
    half = ANGLE_HOH / 2
    acceptor = [
        ("O", tuple(o2)),
        ("H", tuple(o2 + R_OH * (math.cos(half) * bisector + math.sin(half) * x))),
        ("H", tuple(o2 + R_OH * (math.cos(half) * bisector - math.sin(half) * x))),
    ]
    return [donor, acceptor]


# Ring parameters fitted by `python -m fragsim.water fit` (rigid monomers, Hartree-Fock/6-31G energy minimum):
# O...O distance (angstrom), in-plane bend of each donated hydrogen away from the O...O line (degrees),
# and tilt of the free hydrogen from straight up/down toward the ring's outside (degrees).
RING_PARAMS: dict[int, tuple[float, float, float]] = {
    3: (2.931, -21.499, 0.021),
    4: (2.873, -10.118, 0.006),
}


def ring_atoms(n: int, params: tuple[float, float, float] | None = None) -> list[Monomer]:
    """n waters on a regular polygon, each donating one hydrogen to the next oxygen.

    params = (r_oo, bend_deg, tilt_deg); defaults to RING_PARAMS[n]. The free hydrogens alternate up/down.
    """
    r_oo, bend, tilt = params if params is not None else RING_PARAMS[n]
    radius = r_oo / (2.0 * math.sin(math.pi / n))
    oxygens = [
        np.array(
            [radius * math.cos(2 * math.pi * i / n), radius * math.sin(2 * math.pi * i / n), 0.0]
        )
        for i in range(n)
    ]
    b, tl = math.radians(bend), math.radians(tilt)
    waters = []
    for i, o in enumerate(oxygens):
        toward_next = oxygens[(i + 1) % n] - o
        toward_next /= np.linalg.norm(toward_next)
        # rotate about z by `bend` (positive = toward the ring's inside)
        bonded = np.array(
            [
                toward_next[0] * math.cos(b) - toward_next[1] * math.sin(b),
                toward_next[0] * math.sin(b) + toward_next[1] * math.cos(b),
                0.0,
            ]
        )
        outward = np.array([o[0], o[1], 0.0]) / np.linalg.norm(o)
        outward -= np.dot(outward, bonded) * bonded  # in-plane, perpendicular to the bonded O-H
        outward /= np.linalg.norm(outward)
        up = np.array([0.0, 0.0, 1.0 if i % 2 == 0 else -1.0])
        free = math.cos(ANGLE_HOH) * bonded + math.sin(ANGLE_HOH) * (
            math.cos(tl) * up + math.sin(tl) * outward
        )
        waters.append(_water(o, bonded, free))
    return waters


def fit_ring(n: int) -> dict:
    """Minimise the Hartree-Fock energy of an n-ring over (r_oo, bend, tilt) with rigid monomers."""
    from scipy.optimize import minimize

    def energy(x: np.ndarray) -> float:
        if x[0] < 2.4:
            return 1e3
        mol = gto.M(atom=flatten(ring_atoms(n, tuple(x))), basis=BASIS, unit="Angstrom", verbose=0)
        return float(scf.RHF(mol).run().e_tot)

    start = np.array([2.85, 10.0, 0.0])
    result = minimize(
        energy, start, method="Nelder-Mead", options={"xatol": 0.02, "fatol": 2e-5, "maxfev": 120}
    )
    return {
        "n": n,
        "params": [float(v) for v in result.x],
        "e_hf_start": energy(start),
        "e_hf_fit": float(result.fun),
        "evaluations": int(result.nfev),
    }


def water_fragments(n_waters: int) -> list[list[int]]:
    """Atom indices of each water in `flatten(...)` order: one fragment per molecule."""
    return [[3 * i, 3 * i + 1, 3 * i + 2] for i in range(n_waters)]


def flatten(waters: list[Monomer]) -> list[Atom]:
    return [atom for water in waters for atom in water]


def min_intermolecular_distance(waters: list[Monomer]) -> float:
    """Closest pair of atoms on different waters (angstrom)."""
    best = math.inf
    for a in range(len(waters)):
        for b in range(a + 1, len(waters)):
            for _, p in waters[a]:
                for _, q in waters[b]:
                    best = min(best, float(np.linalg.norm(np.array(p) - np.array(q))))
    return best


def sizes(n_waters: int) -> dict:
    """Active-space size of n waters (frozen O 1s)."""
    orbitals, electrons = 12 * n_waters, 8 * n_waters
    return {
        "n_waters": n_waters,
        "active_orbitals": orbitals,
        "active_electrons": electrons,
        "qubits": 2 * orbitals,
        "log10_determinants": 2 * math.log10(math.comb(orbitals, electrons // 2)),
    }


def _energies(atoms: list[Atom], n_core: int, ghost: set[int] = frozenset()) -> dict:
    """RHF, frozen-core CCSD and CCSD(T) total energies (hartree). Ghost atoms carry basis only."""
    spec = [(f"X-{s}" if i in ghost else s, c) for i, (s, c) in enumerate(atoms)]
    mol = gto.M(atom=spec, basis=BASIS, unit="Angstrom", verbose=0)
    t0 = time.perf_counter()
    mf = scf.RHF(mol).run()
    solver = cc.CCSD(mf, frozen=list(range(n_core))).run()
    if not (mf.converged and solver.converged):
        raise RuntimeError("SCF or CCSD did not converge.")
    e_t = solver.ccsd_t()
    return {
        "hf": float(mf.e_tot),
        "ccsd": float(solver.e_tot),
        "ccsd_t": float(solver.e_tot + e_t),
        "seconds": round(time.perf_counter() - t0, 1),
    }


def reference(waters: list[Monomer]) -> dict:
    """Cluster and monomer energies, raw and counterpoise interaction energies at three levels."""
    atoms = flatten(waters)
    n = len(waters)
    cluster = _energies(atoms, n_core=n)
    own = _energies(waters[0], n_core=1)  # all waters share one internal geometry
    cp = []
    for i in range(n):
        ghost = {j for j in range(len(atoms)) if j // 3 != i}
        cp.append(_energies(atoms, n_core=1, ghost=ghost))
    out = {
        "energies_hartree": {"cluster": cluster, "monomer_own_basis": own},
        "interaction_kcal_mol": {},
    }
    for level in ("hf", "ccsd", "ccsd_t"):
        raw = (cluster[level] - n * own[level]) * HARTREE_TO_KCAL
        counterpoise = (cluster[level] - sum(c[level] for c in cp)) * HARTREE_TO_KCAL
        out["interaction_kcal_mol"][level] = {
            "raw": raw,
            "counterpoise": counterpoise,
            "bsse": raw - counterpoise,
            "raw_per_hydrogen_bond": raw / n,
        }
    return out


def monomer_fci_check() -> dict:
    """One water: exact FCI against CCSD(T) in the same space (decides which monomer reference to use)."""
    ref = run_reference(build_mol(monomer_atoms()), n_frozen=1)
    mono = _energies(monomer_atoms(), n_core=1)
    return {
        "e_fci": ref.e_fci,
        "e_ccsd_t": mono["ccsd_t"],
        "fci_minus_ccsd_t_kcal_mol": (ref.e_fci - mono["ccsd_t"]) * HARTREE_TO_KCAL,
        "fci_seconds": round(ref.fci_seconds, 1),
        "qubits": ref.n_qubits,
    }


def sto3g_reference_check() -> dict:
    """How far CCSD(T) is from exact in the dimer interaction energy, in a basis small enough for FCI.

    STO-3G dimer: 12 active orbitals, 6 up electrons (853,776 determinants). A bound on CCSD(T)'s own
    error in E_int for this geometry; it does not transfer exactly to 6-31G.
    """
    from pyscf import mcscf
    from pyscf.fci import direct_spin0

    def energies(waters, n_core):
        mol = gto.M(atom=flatten(waters), basis="sto-3g", unit="Angstrom", verbose=0)
        mf = scf.RHF(mol).run()
        solver = cc.CCSD(mf, frozen=list(range(n_core))).run()
        e_t = solver.e_tot + solver.ccsd_t()
        cas = mcscf.CASCI(mf, mol.nao - n_core, mol.nelectron - 2 * n_core)
        cas.fcisolver = direct_spin0.FCI(mol)
        return float(e_t), float(cas.kernel()[0])

    dimer = dimer_atoms()
    cc_d, fci_d = energies(dimer, 2)
    cc_m, fci_m = energies([dimer[0]], 1)
    e_cc = (cc_d - 2 * cc_m) * HARTREE_TO_KCAL
    e_fci = (fci_d - 2 * fci_m) * HARTREE_TO_KCAL
    return {
        "basis": "sto-3g",
        "e_int_ccsd_t_kcal_mol": e_cc,
        "e_int_fci_kcal_mol": e_fci,
        "ccsd_t_minus_fci_kcal_mol": e_cc - e_fci,
        "monomer_fci_minus_ccsd_t_kcal_mol": (fci_m - cc_m) * HARTREE_TO_KCAL,
    }


def fit_main() -> None:
    for n in (3, 4):
        r = fit_ring(n)
        print(
            f"ring {n}: params (r_oo, bend, tilt) = {tuple(round(v, 3) for v in r['params'])}, "
            f"HF {r['e_hf_start']:.6f} -> {r['e_hf_fit']:.6f} Ha ({r['evaluations']} evaluations)",
            flush=True,
        )


def main() -> None:
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "fit":
        return fit_main()
    DATA.mkdir(parents=True, exist_ok=True)
    record: dict = {
        "what": "Water clusters, 6-31G, frozen O 1s, model geometries; CCSD(T) reference (Lesson 7.1)",
        "geometry": {"r_oh": R_OH, "angle_hoh_deg": math.degrees(ANGLE_HOH), "r_oo": R_OO},
        "monomer_fci_check": monomer_fci_check(),
        "sto3g_reference_check": sto3g_reference_check(),
        "systems": {},
    }
    m = record["monomer_fci_check"]
    print(
        f"monomer: FCI - CCSD(T) = {m['fci_minus_ccsd_t_kcal_mol']:+.3f} kcal/mol ({m['qubits']} qubits)",
        flush=True,
    )
    s = record["sto3g_reference_check"]
    print(
        f"STO-3G dimer: E_int CCSD(T) {s['e_int_ccsd_t_kcal_mol']:.3f}, FCI {s['e_int_fci_kcal_mol']:.3f}, "
        f"difference {s['ccsd_t_minus_fci_kcal_mol']:+.3f} kcal/mol (monomer FCI-CCSD(T) {s['monomer_fci_minus_ccsd_t_kcal_mol']:+.3f})",
        flush=True,
    )
    for name, waters in (
        ("dimer", dimer_atoms()),
        ("trimer ring", ring_atoms(3)),
        ("tetramer ring", ring_atoms(4)),
    ):
        result = reference(waters)
        result["sizes"] = sizes(len(waters))
        result["min_intermolecular_distance"] = min_intermolecular_distance(waters)
        result["atoms"] = flatten(waters)
        record["systems"][name] = result
        e = result["interaction_kcal_mol"]
        print(
            f"{name:14s} E_int raw: HF {e['hf']['raw']:7.2f}  CCSD {e['ccsd']['raw']:7.2f}  CCSD(T) {e['ccsd_t']['raw']:7.2f}"
            f"  | cp CCSD(T) {e['ccsd_t']['counterpoise']:7.2f}  | per H-bond {e['ccsd_t']['raw_per_hydrogen_bond']:6.2f}",
            flush=True,
        )
    RECORD.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {RECORD}")


if __name__ == "__main__":
    main()
