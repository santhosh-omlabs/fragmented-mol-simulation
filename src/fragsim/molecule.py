"""Ammonia geometry and the classical references (RHF, frozen-core FCI)."""

from __future__ import annotations

import math
import time
from dataclasses import dataclass

from pyscf import fci, gto, mcscf, scf

BASIS = "6-31g"


def nh3_atoms(
    r_nh: float = 1.012, angle_hnh_deg: float = 106.7
) -> list[tuple[str, tuple[float, float, float]]]:
    """Pyramidal NH3 (C3v): N at the origin, three H at distance r_nh (angstrom).

    For H-N-H angle a and polar angle b of each N-H bond from the C3 axis:
    cos(a) = cos(b)^2 - 0.5 sin(b)^2, so cos(b)^2 = (cos(a) + 0.5) / 1.5.
    """
    cos_b = math.sqrt((math.cos(math.radians(angle_hnh_deg)) + 0.5) / 1.5)
    sin_b = math.sqrt(1.0 - cos_b**2)
    atoms: list[tuple[str, tuple[float, float, float]]] = [("N", (0.0, 0.0, 0.0))]
    for k in range(3):
        phi = 2.0 * math.pi * k / 3.0
        atoms.append(
            ("H", (r_nh * sin_b * math.cos(phi), r_nh * sin_b * math.sin(phi), -r_nh * cos_b))
        )
    return atoms


def build_mol(atoms=None, basis: str = BASIS) -> gto.Mole:
    """PySCF molecule, singlet, neutral."""
    return gto.M(
        atom=atoms or nh3_atoms(), basis=basis, spin=0, charge=0, unit="Angstrom", verbose=0
    )


@dataclass(frozen=True)
class Reference:
    """Classical numbers every quantum result is judged against (hartree)."""

    n_orbitals_total: int
    n_active_orbitals: int
    n_active_electrons: int
    e_nuc: float
    e_hf: float
    e_fci: float
    fci_seconds: float

    @property
    def n_qubits(self) -> int:
        return 2 * self.n_active_orbitals

    @property
    def correlation_energy(self) -> float:
        return self.e_fci - self.e_hf


def run_reference(mol: gto.Mole | None = None, n_frozen: int = 1) -> Reference:
    """RHF, then FCI in the frozen-core space (N 1s frozen). The exact answer in this basis."""
    mol = mol or build_mol()
    mf = scf.RHF(mol).run()
    if not mf.converged:
        raise RuntimeError("RHF did not converge.")
    n_act = mol.nao - n_frozen
    n_elec = mol.nelectron - 2 * n_frozen
    t0 = time.perf_counter()
    cas = mcscf.CASCI(mf, n_act, n_elec)
    cas.fcisolver = fci.direct_spin0.FCI(mol)
    e_fci = cas.kernel()[0]
    return Reference(
        n_orbitals_total=mol.nao,
        n_active_orbitals=n_act,
        n_active_electrons=n_elec,
        e_nuc=float(mol.energy_nuc()),
        e_hf=float(mf.e_tot),
        e_fci=float(e_fci),
        fci_seconds=time.perf_counter() - t0,
    )


def n_determinants(n_orbitals: int, n_electrons: int) -> int:
    """Size of the FCI space for a singlet: (n_orbitals choose n_alpha)^2."""
    return math.comb(n_orbitals, n_electrons // 2) ** 2


__all__ = ["BASIS", "Reference", "build_mol", "n_determinants", "nh3_atoms", "run_reference"]
