"""Lesson 6: what fragmenting buys on hardware, and how to schedule the pieces on several QPUs.

Two separate effects, kept separate on purpose:

  1. Shots needed.  A circuit with G two-qubit gates, each failing with probability e, runs with no
     error at all with probability p = (1-e)^G (validated against a noisy simulation in Lesson 4).
     To see N_CLEAN error-free runs you need about N_CLEAN / p shots. Fewer gates -> exponentially
     fewer shots. This is the real, physical benefit of fragmenting.
  2. Wall-clock time.  Given the jobs (one per fragment), a scheduler places them on M QPUs.
     Finish time (makespan) is set by the busiest QPU. One big fragment plus three tiny ones
     cannot be balanced, so fragment-parallel speed-up is capped by the largest fragment.

Every hardware number below (e, rep delay, per-job overhead, gate layer time) is an ASSUMPTION
for illustration. Replace them with your device's numbers. Nothing here was measured on a QPU.

Run:  python -m fragsim.schedule   (about 2 minutes, writes data/schedule_study.json)
"""

from __future__ import annotations

import heapq
import json
import math
from dataclasses import dataclass
from pathlib import Path

from pyscf import scf
from qiskit import transpile

from fragsim.embedding import build_atom_basis, density_in_basis
from fragsim.ladder import RUNGS
from fragsim.molecule import build_mol
from fragsim.sampling import BASIS_GATES
from fragsim.solver import build_cluster_problem
from fragsim.sqd_circuit import canonicalize, lucj_circuit, lucj_operator

DATA = Path(__file__).resolve().parents[2] / "data"
RECORD = DATA / "schedule_study.json"
N_CLEAN = 100  # error-free runs wanted per fragment


@dataclass(frozen=True)
class CostModel:
    """Illustrative hardware assumptions (all adjustable)."""

    error_per_cx: float = 0.005  # same value as the Lesson 4 noise model
    rep_delay_s: float = 250e-6  # wait between shots (typical superconducting default)
    layer_s: float = 2e-7  # time per circuit layer
    overhead_s: float = 3.0  # per-job fixed cost: load, compile, classical hand-off


@dataclass(frozen=True)
class Job:
    name: str
    qubits: int
    cx: int  # two-qubit gates, all-to-all connectivity (no routing)
    depth: int


def clean_fraction(job: Job, model: CostModel) -> float:
    """Probability that one shot suffers no two-qubit gate error."""
    return (1.0 - model.error_per_cx) ** job.cx


def shots_needed(job: Job, model: CostModel, n_clean: int = N_CLEAN) -> int:
    return math.ceil(n_clean / clean_fraction(job, model))


def job_seconds(job: Job, shots: int, model: CostModel) -> float:
    return model.overhead_s + shots * (model.rep_delay_s + job.depth * model.layer_s)


def makespan(durations: list[float], n_qpus: int) -> float:
    """Longest-processing-time-first list scheduling: finish time of the busiest QPU."""
    loads = [0.0] * n_qpus
    heapq.heapify(loads)
    for d in sorted(durations, reverse=True):
        heapq.heappush(loads, heapq.heappop(loads) + d)
    return max(loads)


def shot_parallel_seconds(job: Job, shots: int, n_qpus: int, model: CostModel) -> float:
    """One circuit, its shots split evenly over QPUs; each QPU pays the per-job overhead."""
    return job_seconds(job, math.ceil(shots / n_qpus), model)


def _jobs_for(mol, mf, basis, d, fragments, n_virtual) -> list[Job]:
    jobs = []
    for atoms in fragments:
        problem, _ = canonicalize(
            build_cluster_problem(mol, mf, basis, d, atoms, n_virtual=n_virtual)
        )
        circuit = lucj_circuit(problem, lucj_operator(problem, 2))
        ops = transpile(circuit, basis_gates=BASIS_GATES, optimization_level=1)
        jobs.append(
            Job(f"atoms{atoms}", problem.n_qubits, int(ops.count_ops().get("cx", 0)), ops.depth())
        )
    return jobs


DEFAULT_MODEL = CostModel()


def run_study(model: CostModel = DEFAULT_MODEL) -> dict:
    mol = build_mol()
    mf = scf.RHF(mol).run()
    basis = build_atom_basis(mol, mf)
    d = density_in_basis(basis, mf)
    ladder = {
        (r["rung"], r["virtuals_kept"]): r["error_kcal_mol"]
        for r in json.loads((DATA / "ladder_exact_fragments.json").read_text(encoding="utf-8"))[
            "rows"
        ]
    }
    configs = [
        ("whole molecule", [[0, 1, 2, 3]], None, None),
        ("rung 1, 3 empty kept", RUNGS["1: N | H1 | H2 | H3"], 3, ("1: N | H1 | H2 | H3", 3)),
        ("rung 3, 3 empty kept", RUNGS["3: N+H1+H2 | H3"], 3, ("3: N+H1+H2 | H3", 3)),
        (
            "rung 1, all empty kept",
            RUNGS["1: N | H1 | H2 | H3"],
            None,
            ("1: N | H1 | H2 | H3", None),
        ),
        ("rung 3, all empty kept", RUNGS["3: N+H1+H2 | H3"], None, ("3: N+H1+H2 | H3", None)),
    ]
    rows = []
    for label, fragments, n_virtual, key in configs:
        jobs = _jobs_for(mol, mf, basis, d, fragments, n_virtual)
        shots = [shots_needed(j, model) for j in jobs]
        secs = [job_seconds(j, s, model) for j, s in zip(jobs, shots)]
        row = {
            "config": label,
            "error_kcal_mol": ladder[key] if key else 0.0,
            "jobs": [
                {
                    "qubits": j.qubits,
                    "cx": j.cx,
                    "depth": j.depth,
                    "clean_fraction": clean_fraction(j, model),
                    "shots_needed": s,
                    "seconds": t,
                }
                for j, s, t in zip(jobs, shots, secs)
            ],
            "total_qpu_seconds": sum(secs),
            "wall_seconds": {str(m): makespan(secs, m) for m in (1, 2, 4)},
        }
        if len(jobs) == 1:  # a single circuit can only be sped up by splitting its shots
            row["shot_parallel_wall_seconds"] = {
                str(m): shot_parallel_seconds(jobs[0], shots[0], m, model) for m in (1, 2, 4)
            }
        rows.append(row)
        top = max(jobs, key=lambda j: j.cx)
        print(
            f"{label:24s} largest {top.qubits:2d}q {top.cx:5d} CX | shots for {N_CLEAN} clean runs: "
            f"{max(shots):.2e} | QPU time {row['total_qpu_seconds']:.3g} s",
            flush=True,
        )
    return {
        "what": "Shots needed and QPU schedule for NH3 fragment sets (illustrative cost model, no hardware)",
        "cost_model": model.__dict__,
        "n_clean_target": N_CLEAN,
        "rows": rows,
    }


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    RECORD.write_text(json.dumps(run_study(), indent=2) + "\n", encoding="utf-8")
    print(f"wrote {RECORD}")


if __name__ == "__main__":
    main()
