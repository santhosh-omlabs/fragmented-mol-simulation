"""Step 3 of embedding: solve a cluster exactly and read off the fragment's energy share.

Energy bookkeeping (frozen core c, environment core e; p, q, r, s cluster orbitals;
i, j occupied and a, b empty cluster orbitals):

  E_total = E_nuc + E_c + sum over fragments of ( E_HF,frag + E_corr,frag )

  E_HF,frag   = sum_{i in F, k} [ 2 (h + v_c + 1/2 v_e)_ki + sum_j ( 2 (ki|jj) - (kj|ji) ) ]
                the Hartree-Fock (single determinant) energy, weight on the first occupied index
  E_corr,frag = sum_{i in F, j, a, b} c_ij^ab ( 2 (ia|jb) - (ib|ja) )

c_ij^ab are the double-excitation amplitudes of the cluster FCI wavefunction relative to its
Hartree-Fock determinant. The fragment weight touches only the OCCUPIED index i, so throwing
away empty orbitals cannot unbalance the energy split between fragments.
h = bare one-electron operator, v_c / v_e = mean field of the frozen core / environment core.
The 1/2 on v_e stops two fragments from both counting the environment's interaction.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from pyscf import ao2mo, fci, gto, scf
from pyscf.ci import cisd

from fragsim.embedding import AtomBasis, Cluster, make_cluster


@dataclass(frozen=True)
class ClusterResult:
    """Outcome of one cluster solve."""

    atoms: tuple[int, ...]
    n_orbitals: int  # orbitals the solver actually used
    n_electrons: int
    n_virtual_full: int  # empty orbitals available before truncation
    e_frag: float  # this fragment's share, hartree (without E_nuc and E_c)

    @property
    def n_qubits(self) -> int:
        return 2 * self.n_orbitals


def frozen_core_energy(mol: gto.Mole, mf: scf.hf.RHF, n_frozen: int = 1) -> float:
    """E_nuc + energy of the frozen core electrons alone (closed shell)."""
    core = mf.mo_coeff[:, :n_frozen]
    d_core = 2.0 * core @ core.T
    return float(
        mol.energy_nuc()
        + np.einsum("ij,ji->", d_core, mf.get_hcore())
        + 0.5 * np.einsum("ij,ji->", d_core, mf.get_veff(mol, d_core))
    )


def _generic_dipole(mol: gto.Mole) -> np.ndarray:
    """Dipole operator along a fixed, irrational-looking direction; splits symmetry degeneracies."""
    direction = np.array([1.0, 0.6180339887, 0.4142135623])
    return np.einsum("k,kpq->pq", direction, mol.intor("int1e_r", comp=3))


def _mp2_virtual_natural_orbitals(mol, mf, c_occ, c_vir):
    """Cluster MP2 natural virtual orbitals: columns of c_vir rotated, strongest first."""
    fock = mf.get_fock()
    e_o, u_o = np.linalg.eigh(c_occ.T @ fock @ c_occ)
    e_v, u_v = np.linalg.eigh(c_vir.T @ fock @ c_vir)
    co, cv = c_occ @ u_o, c_vir @ u_v
    no, nv = co.shape[1], cv.shape[1]
    ovov = ao2mo.general(mol, (co, cv, co, cv), compact=False).reshape(no, nv, no, nv)
    denom = (
        e_o[:, None, None, None]
        + e_o[None, None, :, None]
        - e_v[None, :, None, None]
        - e_v[None, None, None, :]
    )
    t = ovov.transpose(0, 2, 1, 3) / denom.transpose(0, 2, 1, 3)  # t[i, j, a, b]
    dm_vv = 2.0 * np.einsum("ijac,ijbc->ab", t, t) - np.einsum("ijac,jibc->ab", t, t)
    occupation, vecs = np.linalg.eigh(dm_vv)
    order = np.argsort(occupation)[::-1]
    return cv @ vecs[:, order], occupation[order]


@dataclass(frozen=True)
class ClusterProblem:
    """One fragment's embedded Hamiltonian: what any solver (FCI, SQD, ...) needs."""

    atoms: tuple[int, ...]
    h1: np.ndarray  # one-electron integrals incl. frozen-core and environment mean field
    eri: np.ndarray  # two-electron integrals, chemists' (pq|rs), 4-index
    n_orbitals: int
    n_occ: int  # doubly occupied orbitals in the cluster Hartree-Fock determinant
    n_virtual_full: int  # empty orbitals available before truncation
    h_energy: np.ndarray  # h1 with half the environment field, for the energy share
    q: np.ndarray  # weight of the fragment functions inside each cluster orbital
    mo_coeff: np.ndarray  # AO x n_orbitals coefficients of the cluster orbitals
    gauge: np.ndarray  # dipole along a fixed generic direction, in cluster orbitals (tie-breaker)

    @property
    def nelec(self) -> tuple[int, int]:
        return (self.n_occ, self.n_occ)

    @property
    def n_qubits(self) -> int:
        return 2 * self.n_orbitals


def build_cluster_problem(
    mol: gto.Mole,
    mf: scf.hf.RHF,
    basis: AtomBasis,
    d: np.ndarray,
    atoms: list[int],
    n_virtual: int | None = None,
    n_frozen: int = 1,
    entire_space: bool = False,
) -> ClusterProblem:
    """Fragment + bath, keeping the `n_virtual` strongest MP2 natural virtuals.

    `n_virtual=None` keeps every empty cluster orbital (no truncation).
    `entire_space=True` is a diagnostic: the cluster is every orbital, so the fragment
    shares of any partition must add up to the exact FCI energy.
    """
    cluster: Cluster = make_cluster(basis, d, atoms)
    n = d.shape[0]
    if entire_space:
        z = np.eye(n)
    else:
        z = np.zeros((n, cluster.n_orbitals))
        z[cluster.fragment, np.arange(len(cluster.fragment))] = 1.0
        z[:, len(cluster.fragment) :] = cluster.bath

    # cluster density splits into exactly occupied (1) and exactly empty (0) orbitals
    occupation, vecs = np.linalg.eigh(z.T @ d @ z)
    is_occ = occupation > 0.5
    z_occ, z_vir = z @ vecs[:, is_occ], z @ vecs[:, ~is_occ]
    c = basis.coeff
    c_occ, c_vir = c @ z_occ, c @ z_vir
    n_virtual_full = c_vir.shape[1]

    if n_virtual is not None and n_virtual < n_virtual_full:
        c_vir, _ = _mp2_virtual_natural_orbitals(mol, mf, c_occ, c_vir)
        c_vir = c_vir[:, :n_virtual]
    c_cl = np.hstack([c_occ, c_vir])
    m = c_cl.shape[1]

    # mean fields: frozen core and environment core (occupied, outside the cluster)
    core = mf.mo_coeff[:, :n_frozen]
    d_env = d - z @ (z.T @ d @ z) @ z.T
    d_env_ao = c @ d_env @ c.T
    v_core = mf.get_veff(mol, 2.0 * (core @ core.T))
    v_env = mf.get_veff(mol, 2.0 * d_env_ao)
    hcore = mf.get_hcore()

    # weight of the fragment functions inside the cluster orbitals
    z_cl = c.T @ mf.get_ovlp() @ c_cl
    p_frag = np.zeros((n, n))
    p_frag[cluster.fragment, cluster.fragment] = 1.0
    return ClusterProblem(
        atoms=tuple(atoms),
        h1=c_cl.T @ (hcore + v_core + v_env) @ c_cl,
        eri=ao2mo.restore(1, ao2mo.kernel(mol, c_cl), m),
        n_orbitals=m,
        n_occ=c_occ.shape[1],
        n_virtual_full=n_virtual_full,
        h_energy=c_cl.T @ (hcore + v_core + 0.5 * v_env) @ c_cl,
        q=z_cl.T @ p_frag @ z_cl,
        mo_coeff=c_cl,
        gauge=c_cl.T @ _generic_dipole(mol) @ c_cl,
    )


def fragment_energy(problem: ClusterProblem, ci: np.ndarray) -> float:
    """Fragment's energy share from a (possibly approximate) cluster wavefunction `ci`.

    `ci` is a full-CI vector of shape (C(m, n_occ), C(m, n_occ)) in the cluster orbitals,
    occupied orbitals first. A solver that only fills part of the space (SQD) passes zeros elsewhere.
    """
    m, n_occ, eri, q = problem.n_orbitals, problem.n_occ, problem.eri, problem.q

    # Hartree-Fock part of the share, written with occupied indices only:
    # E_HF = sum_i 2 h_ii + sum_ij ( 2 (ii|jj) - (ij|ji) ), fragment weight on the first i.
    q_oo = q[:n_occ, :n_occ]
    h_oo = problem.h_energy[:n_occ, :n_occ]
    eri_o = eri[:n_occ, :n_occ, :n_occ, :n_occ]
    e_hf = (
        2.0 * np.einsum("ik,ki->", q_oo, h_oo)
        + 2.0 * np.einsum("ik,kijj->", q_oo, eri_o)
        - np.einsum("ik,kjji->", q_oo, eri_o)
    )

    # correlation part: doubles amplitudes of the wavefunction, projected on occupied i
    e_corr = 0.0
    if m > n_occ:
        civec = cisd.from_fcivec(ci, m, problem.nelec)
        c0, _, c2 = cisd.cisdvec_to_amplitudes(civec, m, n_occ)
        c2 = c2 / c0
        ovov = eri[:n_occ, n_occ:, :n_occ, n_occ:]
        pair = 2.0 * ovov.transpose(0, 2, 1, 3) - ovov.transpose(0, 2, 3, 1)  # [i, j, a, b]
        e_corr = np.einsum("ik,kjab,ijab->", q_oo, c2, pair)
    return float(e_hf + e_corr)


def solve_cluster(
    mol: gto.Mole,
    mf: scf.hf.RHF,
    basis: AtomBasis,
    d: np.ndarray,
    atoms: list[int],
    n_virtual: int | None = None,
    n_frozen: int = 1,
    entire_space: bool = False,
) -> ClusterResult:
    """Exact (FCI) solve of one fragment cluster; see `build_cluster_problem` for the options."""
    problem = build_cluster_problem(mol, mf, basis, d, atoms, n_virtual, n_frozen, entire_space)
    _, ci = fci.direct_spin0.FCI(mol).kernel(
        problem.h1, problem.eri, problem.n_orbitals, problem.nelec
    )
    return ClusterResult(
        tuple(atoms),
        problem.n_orbitals,
        2 * problem.n_occ,
        problem.n_virtual_full,
        fragment_energy(problem, ci),
    )
