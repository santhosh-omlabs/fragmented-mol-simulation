"""Step 1 of embedding: rotate occupied orbitals so each sits on an atom or bond."""

from __future__ import annotations

import numpy as np
from pyscf import gto, lo, scf


def localized_occupied(mf: scf.hf.RHF) -> np.ndarray:
    """Intrinsic bond orbitals from all occupied RHF orbitals (AO x n_occ).

    A unitary rotation of the occupied space: same HF energy, same density.
    """
    occ = mf.mo_coeff[:, mf.mo_occ > 0]
    return lo.ibo.ibo(mf.mol, occ, locmethod="IBO")


def atom_weights(mol: gto.Mole, mf: scf.hf.RHF, orbitals: np.ndarray) -> np.ndarray:
    """Fraction (rows sum to 1) of each orbital that lives on each atom.

    Measured in the intrinsic atomic orbital (IAO) basis, which is a minimal
    atom-labelled basis built from the occupied space.
    """
    occ = mf.mo_coeff[:, mf.mo_occ > 0]
    iao = lo.iao.iao(mol, occ)
    iao = lo.vec_lowdin(iao, mf.get_ovlp())
    # orbital coefficients in the IAO basis: C_iao = IAO^T S C
    c_iao = iao.T @ mf.get_ovlp() @ orbitals
    # which atom owns each IAO: the minimal-basis atom labels
    pmol = lo.iao.reference_mol(mol)
    owner = np.array([int(lbl[0]) for lbl in pmol.ao_labels(fmt=None)])
    weights = np.zeros((orbitals.shape[1], mol.natm))
    for atom in range(mol.natm):
        weights[:, atom] = (c_iao[owner == atom, :] ** 2).sum(axis=0)
    return weights / weights.sum(axis=1, keepdims=True)
