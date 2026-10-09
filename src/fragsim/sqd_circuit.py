"""Step 4.2: the LUCJ circuit that proposes which determinants matter.

Pipeline for one fragment:
  cluster problem -> canonical orbitals -> CCSD amplitudes (t1, t2) -> LUCJ operator -> circuit

LUCJ = local unitary cluster Jastrow. The circuit starts from the Hartree-Fock bitstring and
applies a few layers of orbital rotations and number-number phase gates. Its parameters come
from the CCSD doubles amplitudes, so there is no variational training loop.
"""

from __future__ import annotations

from dataclasses import replace

import ffsim
import numpy as np
from pyscf import ao2mo, cc, gto, scf
from qiskit import QuantumCircuit, QuantumRegister

from fragsim.solver import ClusterProblem

GAUGE_STRENGTH = 1e-6  # hartree per a.u. of dipole; only breaks exact degeneracies


def canonicalize(problem: ClusterProblem) -> tuple[ClusterProblem, np.ndarray]:
    """Rotate occupied and empty cluster orbitals so the Fock matrix is (almost) diagonal.

    The rotation mixes occupied with occupied and empty with empty only, so the Hartree-Fock
    state and `fragment_energy` are unchanged. Returns the rotated problem and orbital energies.

    The gauge is fixed so that results repeat from run to run: a tiny dipole term splits the
    exactly degenerate orbitals of a symmetric molecule, and each orbital's sign is chosen so its
    largest atomic-orbital coefficient is positive.
    """
    n_occ, eri = problem.n_occ, problem.eri
    j = np.einsum("pqii->pq", eri[:, :, :n_occ, :n_occ])
    k = np.einsum("piiq->pq", eri[:, :n_occ, :n_occ, :])
    fock = problem.h1 + 2.0 * j - k
    # exactly zero in theory; ~1e-7 in practice because bath orbitals are cut at singular value 1e-6
    if np.abs(fock[:n_occ, n_occ:]).max() > 1e-5:
        raise RuntimeError("Cluster Fock matrix couples occupied and empty orbitals.")
    u = np.zeros_like(fock)
    blocks = ((slice(0, n_occ), slice(0, n_occ)), (slice(n_occ, None), slice(n_occ, None)))
    for block, _ in blocks:
        tilted = fock[block, block] + GAUGE_STRENGTH * problem.gauge[block, block]
        _, vecs = np.linalg.eigh(tilted)
        ao = problem.mo_coeff[:, block] @ vecs
        vecs = vecs * np.sign(ao[np.argmax(np.abs(ao), axis=0), np.arange(ao.shape[1])])
        u[block, block] = vecs
    rotated = replace(
        problem,
        h1=u.T @ problem.h1 @ u,
        eri=np.einsum("pqrs,pa,qb,rc,sd->abcd", eri, u, u, u, u, optimize=True),
        h_energy=u.T @ problem.h_energy @ u,
        q=u.T @ problem.q @ u,
        mo_coeff=problem.mo_coeff @ u,
        gauge=u.T @ problem.gauge @ u,
    )
    return rotated, np.diag(u.T @ fock @ u).copy()


def ccsd_amplitudes(problem: ClusterProblem, eps: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """CCSD (t1, t2) of a canonical cluster problem, using its own integrals."""
    m, n_occ = problem.n_orbitals, problem.n_occ
    mol = gto.M(verbose=0)
    mol.nelectron = 2 * n_occ
    mf = scf.RHF(mol)
    mf.get_hcore = lambda *_: problem.h1
    mf.get_ovlp = lambda *_: np.eye(m)
    mf._eri = ao2mo.restore(8, problem.eri, m)
    mf.mo_coeff = np.eye(m)
    mf.mo_occ = np.array([2.0] * n_occ + [0.0] * (m - n_occ))
    mf.mo_energy = eps
    solver = cc.CCSD(mf)
    solver.kernel()
    if not solver.converged:
        raise RuntimeError("Cluster CCSD did not converge.")
    return solver.t1, solver.t2


def heavy_hex_pairs(n_orbitals: int) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """Interaction pairs that fit a heavy-hex chip: same-spin neighbours, and opposite-spin links every 4th orbital."""
    same_spin = [(p, p + 1) for p in range(n_orbitals - 1)]
    opposite_spin = [(p, p) for p in range(0, n_orbitals, 4)]
    return same_spin, opposite_spin


def lucj_operator(
    problem: ClusterProblem, n_reps: int, use_t1: bool = True, interaction_pairs=None
) -> ffsim.UCJOpSpinBalanced:
    """LUCJ operator from the cluster's CCSD amplitudes (problem must be canonical).

    `interaction_pairs=None` allows every orbital pair (all-to-all). `heavy_hex_pairs(n)` restricts
    the Jastrow part to a sparse set that maps onto real hardware with few swaps.
    """
    _, eps = canonicalize(problem)
    t1, t2 = ccsd_amplitudes(problem, eps)
    return ffsim.UCJOpSpinBalanced.from_t_amplitudes(
        t2, t1=t1 if use_t1 else None, n_reps=n_reps, interaction_pairs=interaction_pairs
    )


def lucj_circuit(problem: ClusterProblem, op: ffsim.UCJOpSpinBalanced) -> QuantumCircuit:
    """Hartree-Fock bitstring, then the LUCJ layers, then measure every qubit."""
    qubits = QuantumRegister(2 * problem.n_orbitals)
    circuit = QuantumCircuit(qubits)
    circuit.append(ffsim.qiskit.PrepareHartreeFockJW(problem.n_orbitals, problem.nelec), qubits)
    circuit.append(ffsim.qiskit.UCJOpSpinBalancedJW(op), qubits)
    circuit.measure_all()
    return circuit


def final_state(problem: ClusterProblem, op: ffsim.UCJOpSpinBalanced) -> np.ndarray:
    """Exact state the circuit prepares, as a (strings x strings) array in pyscf's layout."""
    m = problem.n_orbitals
    vec = ffsim.apply_unitary(ffsim.hartree_fock_state(m, problem.nelec), op, m, problem.nelec)
    side = round(np.sqrt(vec.size))
    return vec.reshape(side, side)


def state_energy(problem: ClusterProblem, state: np.ndarray) -> float:
    """<state|H|state> for a normalized, possibly complex state in pyscf's layout (no constants)."""
    from pyscf.fci import direct_spin1

    m, nelec = problem.n_orbitals, problem.nelec
    h2e = direct_spin1.absorb_h1e(problem.h1, problem.eri, m, nelec, 0.5)
    h_state = direct_spin1.contract_2e(h2e, np.ascontiguousarray(state.real), m, nelec)
    if np.iscomplexobj(state):
        h_state = h_state + 1j * direct_spin1.contract_2e(
            h2e, np.ascontiguousarray(state.imag), m, nelec
        )
    return float(np.real(np.vdot(state, h_state)))
