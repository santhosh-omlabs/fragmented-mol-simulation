"""Steps 4.4-4.5: diagonalize inside the sampled subspace, and repair noisy samples.

SQD in one paragraph: take the bitstrings a quantum computer returned, keep the distinct
up-spin and down-spin electron patterns they contain, and write the Hamiltonian only on the
determinants that pair those patterns. Diagonalize that small matrix classically. The lowest
eigenvalue approximates the ground-state energy. Samples with the wrong electron count are
repaired by flipping bits toward the orbital occupancies of the previous solution.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from pyscf.fci import cistring
from qiskit_addon_sqd.configuration_recovery import recover_configurations
from qiskit_addon_sqd.counts import counts_to_arrays
from qiskit_addon_sqd.fermion import solve_fermion
from qiskit_addon_sqd.subsampling import postselect_by_hamming_right_and_left, subsample

from fragsim.solver import ClusterProblem


@dataclass(frozen=True)
class SubspaceResult:
    """Ground-state estimate from one subspace diagonalization."""

    energy: float  # cluster energy, hartree (same convention as FCI on problem.h1 / eri)
    ci: np.ndarray  # full-CI-shaped amplitudes, zero outside the subspace
    occupancies: tuple[np.ndarray, np.ndarray]  # mean occupation of each up / down orbital
    n_strings: int  # distinct electron patterns per spin
    dimension: int  # determinants in the subspace


def solve_subspace(problem: ClusterProblem, bitstring_matrix: np.ndarray) -> SubspaceResult:
    """Lowest eigenvalue of the cluster Hamiltonian restricted to the sampled subspace."""
    solved = solve_fermion(bitstring_matrix, problem.h1, problem.eri, open_shell=False)
    return _finish(problem, solved)


def solve_strings(problem: ClusterProblem, strings: list[int]) -> SubspaceResult:
    """Same, for an explicit set of electron patterns (integers, bit p = orbital p) for both spins."""
    chosen = sorted({int(s) for s in strings})
    solved = solve_fermion((chosen, chosen), problem.h1, problem.eri, open_shell=True)
    return _finish(problem, solved)


def _finish(problem: ClusterProblem, solved) -> SubspaceResult:
    """Embed a selected-CI result into a full-CI-shaped array so `fragment_energy` can use it."""
    energy, state, occupancies, _ = solved
    norb, n_occ = problem.n_orbitals, problem.n_occ
    side = len(cistring.make_strings(range(norb), n_occ))
    ci = np.zeros((side, side))
    rows = [int(cistring.str2addr(norb, n_occ, int(s))) for s in state.ci_strs_a]
    cols = [int(cistring.str2addr(norb, n_occ, int(s))) for s in state.ci_strs_b]
    ci[np.ix_(rows, cols)] = state.amplitudes
    return SubspaceResult(
        float(energy),
        ci,
        occupancies,
        len(state.ci_strs_a),
        len(state.ci_strs_a) * len(state.ci_strs_b),
    )


@dataclass(frozen=True)
class Iteration:
    """One round of the recovery loop."""

    index: int
    energy: float
    dimension: int
    n_unique_bitstrings: int
    n_repaired: int  # distinct input bitstrings that had a wrong electron count


def run_sqd(
    problem: ClusterProblem,
    counts: dict[str, int],
    iterations: int = 5,
    samples_per_batch: int = 200,
    num_batches: int = 3,
    seed: int = 0,
) -> tuple[SubspaceResult, list[Iteration]]:
    """Iterative SQD with configuration recovery. Returns the best result and the history.

    Round 0 keeps only samples with the right electron count. Each later round repairs the
    wrong-count samples using the previous round's occupancies, then solves `num_batches`
    random batches of `samples_per_batch` bitstrings and keeps the lowest energy.
    """
    norb, (na, nb) = problem.n_orbitals, problem.nelec
    matrix, probs = counts_to_arrays(counts)
    n_wrong = int(np.sum(~_right_count(matrix, norb, na, nb)))
    occupancies = None
    best: SubspaceResult | None = None
    history: list[Iteration] = []
    for index in range(iterations):
        if occupancies is None:
            mat, pr = postselect_by_hamming_right_and_left(
                matrix, probs, hamming_right=na, hamming_left=nb
            )
        else:
            mat, pr = recover_configurations(
                matrix, probs, occupancies, na, nb, rand_seed=seed + index
            )
        batches = subsample(mat, pr, samples_per_batch, num_batches, rand_seed=seed + index)
        results = [solve_subspace(problem, batch) for batch in batches]
        round_best = min(results, key=lambda r: r.energy)
        occupancies = (
            np.mean([r.occupancies[0] for r in results], axis=0),
            np.mean([r.occupancies[1] for r in results], axis=0),
        )
        if best is None or round_best.energy < best.energy:
            best = round_best
        history.append(Iteration(index, round_best.energy, round_best.dimension, len(mat), n_wrong))
    assert best is not None
    return best, history


def _right_count(matrix: np.ndarray, norb: int, na: int, nb: int) -> np.ndarray:
    """Boolean per row: right half has `na` ones (up), left half has `nb` ones (down)."""
    return (matrix[:, norb:].sum(axis=1) == na) & (matrix[:, :norb].sum(axis=1) == nb)


def uniform_random_counts(problem: ClusterProblem, shots: int, seed: int) -> dict[str, int]:
    """Control: bitstrings with the right electron count chosen uniformly at random."""
    rng = np.random.default_rng(seed)
    norb, (na, nb) = problem.n_orbitals, problem.nelec

    def pattern(n_elec: int) -> str:
        bits = np.zeros(norb, dtype=int)
        bits[rng.choice(norb, size=n_elec, replace=False)] = 1
        return "".join(map(str, bits))

    counts: dict[str, int] = {}
    for _ in range(shots):
        bitstring = pattern(nb) + pattern(na)
        counts[bitstring] = counts.get(bitstring, 0) + 1
    return counts
