# Fragment-based sample-based quantum diagonalization on ammonia (6-31G)

**Santhosh Reddy** · tutorial · MIT · v0.2.0 (2026)

**Status:** complete as a study (Lessons 0-7). Not peer-reviewed. Classical simulation throughout, plus one gateway run and three IBM hardware jobs (14 quantum-seconds of a 600-second monthly allowance). Several results are negative and are reported as such.

This repository teaches fragment-based quantum chemistry by doing it, end to end, on small molecules. The question: when a molecule is cut into fragments so that each piece needs fewer qubits and a shallower circuit, how much accuracy does the cut cost, how much hardware error does it avoid, and does a quantum circuit help choose which electron arrangements matter?

On ammonia (NH3, 6-31G, 28 qubits whole) we compare exact full configuration interaction (FCI), sample-based quantum diagonalization (SQD) on the whole molecule, and embedded fragments, using simulation, a classical gateway emulator and IBM hardware, and we test the circuit against free classical baselines at matched subspace size. A final lesson asks the question fragment methods were designed for: the interaction energy of a hydrogen-bonded water dimer.

The method follows two papers (see [docs/LITERATURE.md](docs/LITERATURE.md)). This is a tutorial, not a recreation: it matches neither paper's instance and does not claim to reproduce their numbers.

**Keywords:** ammonia, water dimer, 6-31G, fragmentation, embedding, SQD, FCI, Qiskit, IBM Quantum Runtime.

## Headline findings

All errors are against the exact (or, for water, CCSD(T)) energy in the same basis, in kcal/mol. One kcal/mol is chemical accuracy.

| # | Finding | Evidence |
|---|---|---|
| 1 | Cutting NH3 saves few qubits: each 2 extra qubits per fragment buys about 6 kcal/mol (38.8 at 14 qubits, 17.3 at 22, 11.7 at 24, 5.9 at 26; Hartree-Fock alone is 81.9). | [Results](docs/RESULTS.md), `data/ladder_exact_fragments.json` |
| 2 | Our embedding agrees with an independent code (Vayesta) to 0.17-0.34 kcal/mol on NH3 and 0.01-0.06 on water clusters. | `data/crosscheck_vayesta*.json` |
| 3 | At equal subspace size, an ideal LUCJ circuit beats random and simple choices but not a free classical ranking built from the same CCSD calculation (k = 8 strings: random + Hartree-Fock 36.0, simple ordering 10.5, CCSD-ranked 5.7, ideal circuit 4.8). | `data/matched_subspace_study.json` |
| 4 | An earlier comparison made SQD look exact; that was an artifact (100 random samples cover all 35 strings). Matched subspace size is the fair test. | [Results](docs/RESULTS.md) |
| 5 | On IBM hardware the share of shots with a valid electron count rises from 6.7% to 37% as the routed two-qubit gate count falls from 1,107 to 161, as the (1 - e)^G rule predicts. | `data/ibm_scaling_analysis.json` |
| 6 | Better samples did not give better subspaces: on the whole 28-qubit molecule, hardware SQD gives 3.4 kcal/mol, random bits through the same loop 4.5, a plain classical ordering 0.8. | `data/whole_molecule.json` |
| 7 | What fragmenting does buy is exponentially fewer shots (cost model: about 7.5 million for the whole molecule against about 3,100 for the 14-qubit fragment); a second QPU helps little when one fragment dominates (1.09x on 4 QPUs against about 4x for splitting shots). | `data/schedule_study.json` |
| 8 | Water dimer, one fragment per water, exact solvers: the interaction-energy error falls to 4.8 kcal/mol at 30 qubits, against a 1 kcal/mol target and a 6.6 kcal/mol interaction. The cut itself is cheap in a minimal basis (0.13). The pre-registered target was not met, so no quantum step was run on water. | `data/water_dimer_cut.json` |

Not claimed anywhere: quantum advantage, a match to either paper's numbers, or any result beyond the sizes and one geometry per system studied here.

## Where to start

| You are | Read |
|---|---|
| New to the idea | [docs/INTUITION.md](docs/INTUITION.md), then the findings above |
| Wanting the lessons, not the tables | [docs/TAKEAWAYS.md](docs/TAKEAWAYS.md) |
| Wanting every number and its source | [docs/RESULTS.md](docs/RESULTS.md) |
| Wanting the short version | [docs/SHORT.md](docs/SHORT.md) |
| Wanting the hardware details (job ids, circuits, billed time) | [docs/IBM_HARDWARE_JOBS.md](docs/IBM_HARDWARE_JOBS.md) |
| Wanting to reproduce it | [GUIDE.md](GUIDE.md) (install, verify, regenerate each result) |
| Checking how a result was decided | [docs/RESEARCH_NOTE.md](docs/RESEARCH_NOTE.md) (lab book: what ran, what failed) and [docs/LESSON7_DESIGN.md](docs/LESSON7_DESIGN.md) (pass/fail criteria written before the run) |

## Method

```text
geometry -> RHF (6-31G) -> freeze core 1s -> orbitals tagged to atoms (IAO + atom-local virtuals)
  for each fragment:   bath = SVD of the fragment-environment density block
                       cluster = fragment + bath (+ MP2 natural virtual orbitals, truncated to k)
                       cluster Hamiltonian with frozen-core and environment mean fields
  solver per cluster:  exact FCI                                   (classical reference)
                       SQD: LUCJ circuit from CCSD amplitudes -> sample bitstrings (simulator, Komenco
                            gateway, or IBM hardware) -> configuration recovery -> diagonalize in the
                            subspace spanned by the sampled electron patterns
  energy:              sum over fragments of (Hartree-Fock share + correlation share), occupied indices only
  compare:             exact FCI (NH3) or frozen-core CCSD(T) (water) in the same basis
```

Say out loud: the quantum computer only proposes which electron arrangements matter; a classical computer then diagonalizes the Hamiltonian inside that set. The fragment
energies are not variational, so errors of either sign are possible in principle (all observed errors are positive).

## Reproducibility

Tested with Python 3.12 on Ubuntu under WSL (`pyproject.toml` allows 3.10+, untested); Linux, macOS, or Windows with WSL (PySCF has no Windows wheels). Full instructions, expected output for every command and known failure modes are in [GUIDE.md](GUIDE.md).

```text
python3 -m venv ~/venvs/fragsim
~/venvs/fragsim/bin/pip install -e ".[dev,quantum,ibm]"
~/venvs/fragsim/bin/python -m pytest -q          # 109 tests, about 2 minutes
```

The tests assert the headline numbers from the committed records in `data/`; regenerating a record takes between 1 and 45 minutes per command (listed in the guide).
No IBM job and no gateway call runs unless you run the command that submits it. IBM credentials are never read from the repository.

## Limitations

- **Small molecules, one geometry each.** NH3 is easy classically (FCI of the whole molecule takes about 12 s), so nothing here says anything about molecules too large for classical methods. Water geometries are model structures (monomers at experimental geometry, O...O 2.91 angstrom; rings fitted at Hartree-Fock level), not optimized structures.
- **Small basis.** 6-31G gives absolute energies that differ from experiment; basis-set superposition error is about 30% of the water dimer's raw interaction energy. Interaction energies here are raw, in the same basis, to be compared with the reference only.
- **Embedding is one-shot** (bath from the Hartree-Fock density, no self-consistency) with a hand-written cluster solver. It agrees with Vayesta where both can run (NH3 all-electron; water in STO-3G); the frozen-core path and the virtual truncation in 6-31G have no independent check beyond exact identities in the tests.
- **The SQD circuit's angles come from classical CCSD**, so it cannot hold information CCSD lacks; the finding that a CCSD-ranked classical choice matches it follows from that. Circuits not derived from CCSD were not tried.
- **Hardware evidence is thin:** three jobs on one device (`ibm_fez`) on one day, 3 routing layouts, no repeats of the whole-molecule job, no error mitigation. The Komenco gateway returns exact probabilities and is treated as a classical emulator.
- **The shot and scheduling numbers are a cost model** with assumed hardware parameters (0.5% error per two-qubit gate, 250 microseconds between shots, 3 s per job), not measurements. "Error-free run" is a pessimistic proxy because SQD tolerates some noisy shots.
- **Water was not carried through.** The untruncated 32-qubit point of the equal-waters scheme and a 30-qubit point of the unequal split could not be run on a 7 GB machine (a memory kill and a silent exit); the verdict rests on the curve through 30 qubits. The hypothesis for why 6-31G fails where STO-3G does not is untested.
- Third-party code and services are listed in the next section; none of them is part of this work.

## Third-party code and services

| Item | Provider | How it is used here |
|---|---|---|
| Komenco gateway and Python client | Automatski, https://automatski.com/platform.html (a provider separate from this repository and its author) | The client in `src/fragsim/backends/komenco.py` is adapted from the Python client supplied for the platform. One shared trial endpoint (15-qubit cap, open key, plain HTTP) was called twice, to run one circuit. Treated as a classical emulator because it returns exact probabilities. Check the provider's terms before redistributing the adapted client. |
| IBM Quantum hardware (`ibm_fez`) and Qiskit Runtime | IBM | Three jobs, 14 quantum-seconds; details in [docs/IBM_HARDWARE_JOBS.md](docs/IBM_HARDWARE_JOBS.md). |
| Vayesta | Booth group, https://github.com/BoothGroup/Vayesta | Independent cross-check, used from a source checkout (commit 7f1639d); not included. |
| PySCF, Qiskit, Qiskit Aer, ffsim, qiskit-addon-sqd | their authors | Integrals and reference energies; circuits and simulation; LUCJ circuits; SQD subspace solves. Versions in the lab book. |

## Layout

```text
README.md  GUIDE.md  LICENSE  CITATION.cff  pyproject.toml  verify.sh  verify.ps1
src/fragsim/
  molecule.py localize.py embedding.py solver.py     NH3 reference; orbitals, bath, cluster, energy share (Lessons 1-3)
  ladder.py                                          python -m fragsim.ladder          -> data/ladder_exact_fragments.json
  sqd_circuit.py sampling.py sqd.py                  LUCJ circuit, sampling, configuration recovery, subspace solve (Lesson 4)
  matched.py                                         python -m fragsim.matched         -> data/matched_subspace_study.json
  komenco_run.py backends/komenco.py                 gateway client and run             -> data/komenco_run.json
  schedule.py                                        shots and multi-QPU scheduling model (Lesson 6)
  hwfriendly.py ibm_dry_run.py ibm_run.py ibm_analysis.py ibm_scaling*.py whole_molecule.py   IBM steps (Lesson 5)
  water.py water_cut.py                              water clusters, CCSD(T) reference, cut-error curve (Lesson 7)
scripts/                                             Vayesta cross-checks, STO-3G water check
tests/                                               headline numbers and committed records
data/                                                small JSON records, each with the command that made it (data/README.md)
docs/                                                SHORT, INTUITION, RESULTS, TAKEAWAYS, IBM_HARDWARE_JOBS, ROADMAP, LITERATURE, RESEARCH_NOTE, LESSON7_DESIGN
```

## Citation

```text
@software{reddy2026fragsim,
  author  = {Santhosh Reddy},
  title   = {Fragment-based sample-based quantum diagonalization on ammonia (6-31G)},
  year    = {2026},
  version = {0.2.0},
  license = {MIT},
  note    = {Tutorial with simulation, a gateway run and three IBM jobs; includes negative results.}
}
```

See also [CITATION.cff](CITATION.cff). Cite the two papers in [docs/LITERATURE.md](docs/LITERATURE.md) for the methods.
