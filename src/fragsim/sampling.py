"""Step 4.3: sample bitstrings from the LUCJ circuit and measure how useful they are.

Bitstring convention (Qiskit): 2*norb characters, rightmost = qubit 0.
Qubits 0..norb-1 hold the up (alpha) orbitals, norb..2*norb-1 the down (beta) orbitals, so
    bitstring = [beta bits][alpha bits]     (left half beta, right half alpha)
and orbital p of a spin is the character p places from the right end of its half.
"""

from __future__ import annotations

import numpy as np
from pyscf.fci import cistring
from qiskit import QuantumCircuit, transpile

from fragsim.solver import ClusterProblem

BASIS_GATES = ["rz", "sx", "x", "cx"]


def split_spins(bitstring: str, norb: int) -> tuple[str, str]:
    """(alpha, beta) halves of a Qiskit bitstring."""
    return bitstring[norb:], bitstring[:norb]


def has_correct_electron_count(bitstring: str, norb: int, nelec: tuple[int, int]) -> bool:
    """True when the up half has nelec[0] ones and the down half has nelec[1] ones."""
    alpha, beta = split_spins(bitstring, norb)
    return alpha.count("1") == nelec[0] and beta.count("1") == nelec[1]


def sample_ideal(
    state: np.ndarray, problem: ClusterProblem, shots: int, seed: int
) -> dict[str, int]:
    """Draw `shots` bitstrings from an exact state (ideal, noiseless circuit); reproducible by seed."""
    norb = problem.n_orbitals
    if problem.nelec[0] != problem.nelec[1]:
        raise ValueError("sample_ideal assumes a closed-shell cluster.")
    strings = cistring.make_strings(range(norb), problem.nelec[0])
    probs = np.abs(state).ravel() ** 2
    draws = np.random.default_rng(seed).multinomial(shots, probs / probs.sum())
    counts: dict[str, int] = {}
    for flat in np.nonzero(draws)[0]:
        row, col = divmod(int(flat), state.shape[1])
        alpha, beta = format(int(strings[row]), f"0{norb}b"), format(int(strings[col]), f"0{norb}b")
        counts[beta + alpha] = int(draws[flat])
    return counts


def sample_noisy(
    circuit: QuantumCircuit, shots: int, p_cx: float, p_1q: float, seed: int
) -> dict[str, int]:
    """Sample the compiled circuit with depolarizing noise on every gate (illustrative model).

    Statevector trajectories: each shot is an independent noisy run.
    """
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, depolarizing_error

    noise = NoiseModel(basis_gates=BASIS_GATES)
    noise.add_all_qubit_quantum_error(depolarizing_error(p_cx, 2), ["cx"])
    noise.add_all_qubit_quantum_error(depolarizing_error(p_1q, 1), ["sx", "x"])
    backend = AerSimulator(method="statevector", noise_model=noise, seed_simulator=seed)
    compiled = transpile(circuit, backend, basis_gates=BASIS_GATES, seed_transpiler=7)
    return dict(backend.run(compiled, shots=shots).result().get_counts())


def determinant_index(bitstring: str, norb: int, nelec: tuple[int, int]) -> tuple[int, int]:
    """Row (up string) and column (down string) of this bitstring in a pyscf FCI array."""
    alpha, beta = split_spins(bitstring, norb)
    return (
        int(cistring.str2addr(norb, nelec[0], int(alpha, 2))),
        int(cistring.str2addr(norb, nelec[1], int(beta, 2))),
    )


def sample_statistics(
    counts: dict[str, int], problem: ClusterProblem, exact_ci: np.ndarray
) -> dict[str, float]:
    """How good a sample is, measured against the exact cluster ground state `exact_ci`.

    Keys: shots, distinct bitstrings, fraction of shots with the right electron count,
    distinct up and down strings among right-count samples, subspace size (up x down),
    and the exact ground-state probability that lies inside that subspace.
    """
    norb, nelec = problem.n_orbitals, problem.nelec
    shots = sum(counts.values())
    good = {b: c for b, c in counts.items() if has_correct_electron_count(b, norb, nelec)}
    alphas = {split_spins(b, norb)[0] for b in good}
    betas = {split_spins(b, norb)[1] for b in good}
    rows = sorted({determinant_index(b, norb, nelec)[0] for b in good})
    cols = sorted({determinant_index(b, norb, nelec)[1] for b in good})
    inside = float(np.sum(np.abs(exact_ci[np.ix_(rows, cols)]) ** 2)) if rows and cols else 0.0
    return {
        "shots": shots,
        "distinct_bitstrings": len(counts),
        "fraction_correct_electron_count": sum(good.values()) / shots,
        "distinct_up_strings": len(alphas),
        "distinct_down_strings": len(betas),
        "subspace_size": len(alphas) * len(betas),
        "exact_weight_in_subspace": inside,
    }
