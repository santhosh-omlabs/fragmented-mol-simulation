"""Step 2 of embedding: atom-tagged orthonormal basis, fragment/environment split, bath.

Notation (all in an orthonormal basis B whose columns are tagged to atoms):
  D   one-spin Hartree-Fock density matrix, eigenvalues 1 (occupied) and 0 (empty)
  F   indices of basis functions on the fragment's atoms
  E   all other indices (the environment)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from pyscf import gto, lo, scf


@dataclass(frozen=True)
class AtomBasis:
    """Orthonormal basis of the full AO space with an atom label per column."""

    coeff: np.ndarray  # AO x n, columns are orthonormal in the AO overlap metric
    atom: np.ndarray  # length n, atom index owning each column
    n_iao: int  # first n_iao columns are IAOs, the rest are virtual (PAO-like)


def build_atom_basis(mol: gto.Mole, mf: scf.hf.RHF, n_frozen: int = 1) -> AtomBasis:
    """Orthonormal, atom-tagged basis of the ACTIVE space (frozen core removed).

    IAOs carry the occupied space; the IAO that is the frozen core is dropped.
    The remaining empty functions are atom-local virtuals. Columns: n_ao - n_frozen.
    """
    s = mf.get_ovlp()
    n_occ = int(np.sum(mf.mo_occ > 0))
    core = mf.mo_coeff[:, :n_frozen]
    iao = lo.vec_lowdin(lo.iao.iao(mol, mf.mo_coeff[:, :n_occ]), s)
    pmol = lo.iao.reference_mol(mol)
    iao_atom = np.array([int(lbl[0]) for lbl in pmol.ao_labels(fmt=None)])

    # Drop the IAOs that carry the frozen core. With several equivalent cores (two oxygens) the
    # core MOs are mixed combinations (1s_A +- 1s_B), so match the core SPACE, not orbital by orbital:
    # take the n_frozen distinct IAOs with the largest weight inside it.
    weight_in_core = np.sum((core.T @ s @ iao) ** 2, axis=0)
    drop = {int(i) for i in np.argsort(weight_in_core)[::-1][:n_frozen]}
    keep = [i for i in range(iao.shape[1]) if i not in drop]
    core_atoms = [int(iao_atom[i]) for i in sorted(drop)]
    iao = iao[:, keep] - core @ (core.T @ s @ iao[:, keep])
    iao = lo.vec_lowdin(iao, s)
    iao_atom = iao_atom[keep]

    # atom-local virtuals: AOs of one atom with IAO and core directions projected out
    proj = np.eye(mol.nao) - iao @ iao.T @ s - core @ core.T @ s
    sliced = mol.aoslice_by_atom()
    vir_cols, vir_atom = [], []
    for atom in range(mol.natm):
        a0, a1 = sliced[atom, 2], sliced[atom, 3]
        n_keep = (a1 - a0) - int(np.sum(iao_atom == atom)) - core_atoms.count(atom)
        u, _, _ = np.linalg.svd(s @ proj[:, a0:a1], full_matrices=False)
        vir_cols.append(np.linalg.solve(s, u[:, :n_keep]))
        vir_atom += [atom] * n_keep
    vir = lo.vec_lowdin(np.hstack(vir_cols), s)  # make atom blocks mutually orthogonal
    coeff = np.hstack([iao, vir])
    return AtomBasis(coeff, np.concatenate([iao_atom, vir_atom]), iao.shape[1])


def density_in_basis(basis: AtomBasis, mf: scf.hf.RHF, n_frozen: int = 1) -> np.ndarray:
    """One-spin density of the active occupied RHF orbitals, written in the atom basis."""
    s = mf.get_ovlp()
    occ_in_b = basis.coeff.T @ s @ mf.mo_coeff[:, n_frozen : int(np.sum(mf.mo_occ > 0))]
    return occ_in_b @ occ_in_b.T


@dataclass(frozen=True)
class Cluster:
    """A fragment plus its bath, expressed in the atom basis."""

    fragment: np.ndarray  # basis indices (columns of the atom basis)
    bath: np.ndarray  # (n_basis x n_bath) bath orbitals, in the atom basis
    singular_values: np.ndarray  # strength of fragment-environment entanglement
    n_electrons: int  # electrons living in fragment + bath

    @property
    def n_orbitals(self) -> int:
        return len(self.fragment) + self.bath.shape[1]


def make_cluster(basis: AtomBasis, d: np.ndarray, atoms: list[int], tol: float = 1e-6) -> Cluster:
    """Bath from the SVD of the fragment-environment block of D."""
    n = d.shape[0]
    frag = np.where(np.isin(basis.atom, atoms))[0]
    env = np.setdiff1d(np.arange(n), frag)
    # D_EF = U S V^T : columns of U (with S > tol) are the environment orbitals
    # that the fragment can see. Everything else in E is untouched core or empty.
    u, sv, _ = np.linalg.svd(d[np.ix_(env, frag)], full_matrices=False)
    keep = sv > tol
    bath = np.zeros((n, int(keep.sum())))
    bath[env, :] = u[:, keep]
    cluster_basis = np.zeros((n, len(frag) + bath.shape[1]))
    cluster_basis[frag, np.arange(len(frag))] = 1.0
    cluster_basis[:, len(frag) :] = bath
    n_elec = 2 * float(np.trace(cluster_basis.T @ d @ cluster_basis))
    return Cluster(frag, bath, sv, round(n_elec))
