"""Step 4.5: fair comparison at matched subspace size.

Every method below picks k electron patterns (strings) per spin; the SQD subspace is all pairs
of them, k*k determinants. Same k, same Hamiltonian, so the only difference is WHICH k strings.

  oracle        the k strings with the most weight in the exact ground state (best possible)
  circuit       ideal LUCJ circuit: k most probable strings (infinite shots) and
                the first k strings seen while sampling (finite shots)
  noisy circuit noisy samples, repaired by configuration recovery: k strings with most mass
  excitation    classical ordering: fewest electrons moved out of the Hartree-Fock pattern first,
                then lowest sum of orbital energies
  ccsd ranked   classical, no circuit: strings ranked by weight in the CCSD wavefunction
                (the same classical calculation that sets the circuit's angles)
  random        k strings chosen uniformly (median over trials); also with Hartree-Fock forced in

Run:  python -m fragsim.matched   (a few minutes, writes data/matched_subspace_study.json)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from pyscf import fci, scf
from pyscf.ci import cisd
from pyscf.fci import cistring
from qiskit_addon_sqd.configuration_recovery import recover_configurations
from qiskit_addon_sqd.counts import counts_to_arrays
from qiskit_addon_sqd.subsampling import postselect_by_hamming_right_and_left

from fragsim.embedding import build_atom_basis, density_in_basis
from fragsim.molecule import build_mol
from fragsim.sampling import sample_noisy
from fragsim.solver import ClusterProblem, build_cluster_problem
from fragsim.sqd import solve_strings, solve_subspace
from fragsim.sqd_circuit import (
    canonicalize,
    ccsd_amplitudes,
    final_state,
    lucj_circuit,
    lucj_operator,
)

HARTREE_TO_KCAL = 627.5095
DATA = Path(__file__).resolve().parents[2] / "data"
RECORD = DATA / "matched_subspace_study.json"
NOISY_COUNTS = DATA / "noisy_counts_n_fragment_2layers.json"
KS = [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 35]
NOISE = {"p_cx": 0.005, "p_1q": 0.0005, "shots": 1000, "seed": 5}


def all_strings(problem: ClusterProblem) -> np.ndarray:
    """Every electron pattern for one spin, as integers (bit p set = orbital p occupied)."""
    return cistring.make_strings(range(problem.n_orbitals), problem.n_occ)


def _error(problem: ClusterProblem, strings, e_fci: float) -> float:
    return (solve_strings(problem, list(strings)).energy - e_fci) * HARTREE_TO_KCAL


def oracle_order(problem: ClusterProblem, ci: np.ndarray) -> list[int]:
    """Strings ranked by their weight in the exact ground state."""
    strings = all_strings(problem)
    weight = np.sum(ci**2, axis=1) + np.sum(ci**2, axis=0)
    return [int(strings[i]) for i in np.argsort(-weight, kind="stable")]


def circuit_order(problem: ClusterProblem, state: np.ndarray) -> list[int]:
    """Strings ranked by probability under the ideal circuit (the infinite-shot limit)."""
    strings = all_strings(problem)
    mass = np.sum(np.abs(state) ** 2, axis=1) + np.sum(np.abs(state) ** 2, axis=0)
    return [int(strings[i]) for i in np.argsort(-mass, kind="stable")]


def circuit_first_seen(
    problem: ClusterProblem, state: np.ndarray, seed: int, shots: int = 20000
) -> tuple[list[int], list[int]]:
    """Strings in order of first appearance while sampling, with the shot index of each."""
    strings = all_strings(problem)
    probs = np.abs(state).ravel() ** 2
    draws = np.random.default_rng(seed).choice(probs.size, size=shots, p=probs / probs.sum())
    seen: dict[int, int] = {}
    for shot, flat in enumerate(draws, start=1):
        row, col = divmod(int(flat), state.shape[1])
        for index in (row, col):
            seen.setdefault(int(strings[index]), shot)
    ordered = sorted(seen.items(), key=lambda kv: kv[1])
    return [s for s, _ in ordered], [n for _, n in ordered]


def excitation_order(problem: ClusterProblem, eps: np.ndarray) -> list[int]:
    """Classical baseline: fewest excitations out of the Hartree-Fock pattern, then lowest energy."""
    hf = (1 << problem.n_occ) - 1
    strings = [int(s) for s in all_strings(problem)]
    level = {s: (s & ~hf).bit_count() for s in strings}
    energy = {s: sum(eps[p] for p in range(problem.n_orbitals) if s >> p & 1) for s in strings}
    return sorted(strings, key=lambda s: (level[s], energy[s]))


def ccsd_order(problem: ClusterProblem, eps: np.ndarray) -> list[int]:
    """Classical, no circuit: strings ranked by their weight in the CCSD wavefunction."""
    t1, t2 = ccsd_amplitudes(problem, eps)
    vec = cisd.to_fcivec(cisd.amplitudes_to_cisdvec(1.0, t1, t2), problem.n_orbitals, problem.nelec)
    return circuit_order(problem, np.asarray(vec).reshape(len(all_strings(problem)), -1))


def noisy_order(problem: ClusterProblem, counts: dict[str, int], seed: int = 0) -> list[int]:
    """Strings ranked by mass after one round of configuration recovery on noisy samples."""
    norb, (na, nb) = problem.n_orbitals, problem.nelec
    matrix, probs = counts_to_arrays(counts)
    kept, _ = postselect_by_hamming_right_and_left(matrix, probs, hamming_right=na, hamming_left=nb)
    occupancies = solve_subspace(problem, kept).occupancies
    repaired, weights = recover_configurations(matrix, probs, occupancies, na, nb, rand_seed=seed)
    weight_of: dict[int, float] = {}
    for bits, w in zip(repaired, weights):
        for half in (bits[:norb], bits[norb:]):
            value = int("".join("1" if b else "0" for b in half), 2)
            weight_of[value] = weight_of.get(value, 0.0) + float(w)
    return sorted(weight_of, key=lambda s: -weight_of[s])


def _noisy_counts(problem: ClusterProblem) -> dict[str, int]:
    if NOISY_COUNTS.exists():
        return json.loads(NOISY_COUNTS.read_text(encoding="utf-8"))
    circuit = lucj_circuit(problem, lucj_operator(problem, 2))
    counts = sample_noisy(circuit, NOISE["shots"], NOISE["p_cx"], NOISE["p_1q"], seed=NOISE["seed"])
    NOISY_COUNTS.parent.mkdir(parents=True, exist_ok=True)
    NOISY_COUNTS.write_text(json.dumps(counts, indent=1) + "\n", encoding="utf-8")
    return counts


def _quartiles(values: list[float]) -> list[float]:
    return [float(np.quantile(values, 0.25)), float(np.quantile(values, 0.75))]


def run_study(n_random: int = 100, n_circuit_seeds: int = 50) -> dict:
    """Energy error (kcal/mol) against exact cluster FCI for every method and every k."""
    mol = build_mol()
    mf = scf.RHF(mol).run()
    basis = build_atom_basis(mol, mf)
    raw = build_cluster_problem(mol, mf, basis, density_in_basis(basis, mf), [0], n_virtual=3)
    problem, eps = canonicalize(raw)
    e_fci, ci = fci.direct_spin0.FCI(mol).kernel(problem.h1, problem.eri, 7, problem.nelec)
    state = final_state(problem, lucj_operator(problem, 2))
    hf = (1 << problem.n_occ) - 1
    pool = [int(s) for s in all_strings(problem)]
    orders = {
        "oracle": oracle_order(problem, ci),
        "circuit_infinite_shots": circuit_order(problem, state),
        "excitation_ordered": excitation_order(problem, eps),
        "ccsd_ranked_classical": ccsd_order(problem, eps),
        "noisy_circuit_recovered": noisy_order(problem, _noisy_counts(problem)),
    }
    sequences = [circuit_first_seen(problem, state, seed) for seed in range(n_circuit_seeds)]
    rng = np.random.default_rng(0)
    rows = []
    for k in KS:
        row: dict = {"k": k, "determinants": k * k}
        for name, order in orders.items():
            row[name] = _error(problem, order[:k], e_fci)
        reached = [(seq, seen) for seq, seen in sequences if len(seq) >= k]
        finite = [_error(problem, seq[:k], e_fci) for seq, _ in reached]
        # only report a finite-shot value when most seeds actually produced k distinct strings
        enough = len(reached) >= len(sequences) // 2
        row["circuit_finite_shots_median"] = float(np.median(finite)) if enough else None
        row["circuit_finite_shots_q25_q75"] = _quartiles(finite) if enough else None
        row["circuit_finite_seeds_reaching_k"] = len(reached)
        shots = [seen[k - 1] for _, seen in reached]
        row["circuit_shots_to_reach_k_median"] = float(np.median(shots)) if enough else None
        plain, with_hf = [], []
        for _ in range(n_random):
            plain.append(_error(problem, rng.choice(pool, size=k, replace=False), e_fci))
            rest = [s for s in pool if s != hf]
            pick = [hf, *rng.choice(rest, size=k - 1, replace=False)]
            with_hf.append(_error(problem, pick, e_fci))
        row["random_median"] = float(np.median(plain))
        row["random_q25_q75"] = _quartiles(plain)
        row["random_with_hf_median"] = float(np.median(with_hf))
        row["random_with_hf_q25_q75"] = _quartiles(with_hf)
        rows.append(row)
        print(
            f"k={k:2d} oracle {row['oracle']:7.2f} | circuit(inf) {row['circuit_infinite_shots']:7.2f} "
            f"circuit(finite) {row['circuit_finite_shots_median'] or float('nan'):7.2f} "
            f"| ccsd-ranked {row['ccsd_ranked_classical']:7.2f} "
            f"| noisy {row['noisy_circuit_recovered']:7.2f} | excitation {row['excitation_ordered']:7.2f} "
            f"| random {row['random_median']:7.2f} random+HF {row['random_with_hf_median']:7.2f}",
            flush=True,
        )
    return {
        "what": "SQD energy error (kcal/mol) vs exact cluster FCI at matched number of strings per spin",
        "fragment": "N fragment of rung 1, 3 empty orbitals kept: 7 orbitals, 14 qubits, 1,225 determinants",
        "circuit": "LUCJ from CCSD amplitudes, 2 layers",
        "noise_model": "depolarizing, 0.5% per CNOT, 0.05% per one-qubit gate, 1000 shots (illustrative)",
        "n_random_trials": n_random,
        "n_circuit_seeds": n_circuit_seeds,
        "rows": rows,
    }


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    RECORD.write_text(json.dumps(run_study(), indent=2) + "\n", encoding="utf-8")
    print(f"wrote {RECORD}")


if __name__ == "__main__":
    main()
