# Roadmap

| Lesson | Topic | Output | Status |
|---|---|---|---|
| 0 | Why fragment: cost scaling | docs/INTUITION.md | written |
| 1 | Exact baseline: HF and FCI for NH3 in 6-31G | src/fragsim/molecule.py, tests | done |
| 2 | Cut the molecule: localized orbitals, fragments, bath | localize.py, embedding.py | done |
| 3 | Fragmentation error with exact fragment solvers, virtual truncation | solver.py, ladder.py, data/ladder_exact_fragments.json | done |
| 4 | SQD: LUCJ circuit, sampling, configuration recovery, matched-size test | sqd_circuit.py, sampling.py, sqd.py, matched.py, komenco_run.py, data/*.json | done (classical simulation, Komenco run) |
| 5 | IBM Runtime: budget estimate, dry run, jobs on 14- and 28-qubit circuits | ibm_dry_run.py, ibm_run.py, ibm_scaling.py, whole_molecule.py, hwfriendly.py | done (3 jobs, 14 of 600 quantum-seconds) |
| 6 | Shots and scheduling: shot-parallel vs fragment-parallel | schedule.py, data/schedule_study.json | done (cost model, no hardware) |
| 7 | Fragmenting a water cluster: interaction energy with weak cuts, CCSD(T) reference | docs/LESSON7_DESIGN.md | closed, partial and negative result: design written; 7.1 reference data done (`python -m fragsim.water`); 7.2 generalised pipeline done; 7.3 Vayesta cross-check done (STO-3G); 7.4 cut-error curve through k = 9: tier T1 not met (4.8 kcal/mol at 30 qubits; untruncated point not obtained); unequal split tested and worse; 7.5 onward not run (hard stop). Closed as a partial and negative result |

IBM allowance: 10 minutes per month. Target 4 minutes, hard cap 10. Estimate with a dry run before every submission,
and ask before each real submission.
