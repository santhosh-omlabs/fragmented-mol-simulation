# Research note (lab book)

In date order: what ran, what failed, what I decided and why. Written in the first person.

## 2026-10-09 - setup and baseline

Environment: Windows 11, WSL Ubuntu (Python 3.12.3), 12 cores, 7 GB RAM. PySCF 2.14.0, Qiskit 2.5.2,
qiskit-addon-sqd 0.14.0, ffsim 0.0.84, qiskit-ibm-runtime 0.50.0.

What failed: the first venv, created inside the repo on the mounted `D:` drive, was corrupted during pip install
(`No module named 'pip._vendor.chardet.langturkishmodel'`). I rebuilt it in `~/venvs/fragsim` on the WSL disk.

What ran: RHF and frozen-core FCI for NH3 / 6-31G at N-H = 1.012 A, H-N-H = 106.7 degrees.
E_HF = -56.161063, E_FCI = -56.291514 hartree, 1,002,001 determinants, FCI in about 12 s.

Komenco gateway (a third-party platform provided by Automatski, https://automatski.com/platform.html): the client I adapted from the one supplied for it returns exact probabilities with no shot noise, and the
vendor's own note calls the 15-qubit instance a shared trial simulator. An earlier zero-error match is what an ideal
simulator gives, and a physical QPU without error correction would show noise. This repo therefore treats Komenco as a
classical emulator. At the time of this entry no request to the gateway had been made yet and the client was tested only with a fake transport (superseded: see the Komenco entry dated later).

## 2026-10-09 - embedding and the fragment ladder (Lessons 2-3)

Built localization (IBO), an atom-tagged orthonormal basis with the N 1s core frozen out (14 functions), the bath from the SVD of the
fragment-environment block of the density matrix, and an exact cluster solver.

Two energy formulas failed before the third worked, and the failures are worth keeping:
1. Weighting the fragment on every orbital index of the cluster 1- and 2-body density matrices. With whole-molecule and full-space
   checks passing, but empty orbitals truncated, the total came out 300+ kcal/mol above Hartree-Fock (an impossible value).
   Cause: the fragment weight touched empty orbitals, so dropping some of them broke the cancellation between fragments.
2. Projected amplitudes for the correlation part but the old formula for the Hartree-Fock part: still wrong for the same reason.
3. Final: both parts written with the fragment weight on occupied indices only. Checks that pass: a whole-molecule cluster equals FCI
   (7e-9 hartree), a full-space cluster partition adds up to FCI (7e-9), zero kept empty orbitals gives exactly Hartree-Fock.

Result: untruncated errors of 17.3, 11.7 and 5.9 kcal/mol for three rungs of the ladder, at 22, 24 and 26 qubits for the largest
fragment (28 for the whole). Truncating to the 15-qubit Komenco cap leaves about 39 kcal/mol.

## 2026-10-09 - Vayesta cross-check

`pip install vayesta` finds nothing on PyPI; `pip install git+https://github.com/BoothGroup/Vayesta` fails at cmake and then at
`Could NOT find BLAS` (needs a BLAS development package, i.e. sudo apt). Running from a source checkout via `PYTHONPATH` works
(commit 7f1639d); only the optional `dyson` module is missing.

Vayesta's `frozen` options mean environment orbitals inside a cluster, not a frozen molecular core, so the comparison ran all-electron
on both sides (fragsim with `n_frozen=0`). Result: differences of 0.34, 0.28 and 0.17 kcal/mol for the three partitions, fragsim
always slightly higher. I did not track down the cause. Each Vayesta run of the larger partitions takes about 3 minutes (FCI on
12 to 13 orbitals, single thread).

## 2026-10-09 - SQD on the 14-qubit N fragment (Lesson 4, steps 4.1-4.4)

Problem: N fragment of rung 1, 3 empty orbitals kept: 7 orbitals, 4 up + 4 down electrons, 14 qubits, 1,225 determinants.

Reproducibility bug (found, fixed): the same code gave different circuits in different processes. NH3 is symmetric, so some cluster
orbitals are exactly degenerate; the diagonalizer returns an arbitrary rotation inside a degenerate pair, multithreaded floating-point noise
picks a different one each run, and the LUCJ circuit (which keeps only the 2 leading terms of a factorization of the CCSD amplitudes)
depends on that gauge. The integrals differed by 1.5 between runs and the circuit probabilities by 1e-3. Fix: `canonicalize` adds a 1e-6
dipole term along a fixed generic direction to split the degeneracies and fixes every orbital's sign from its largest AO coefficient;
now the integrals repeat to 5e-10 and the state to 3e-11. An earlier `ffsim.sample_state_vector` call was also replaced by a seeded multinomial
sampler. Numbers quoted before this fix (the first Step 4.3 table) came from one non-reproducible draw and are superseded.

Ideal LUCJ circuit, 2 layers (96.7% of the probability sits on the Hartree-Fock bitstring), SQD error against exact cluster FCI,
20 seeds, mean +- std: 30 shots 34.8 +- 7.6 kcal/mol; 100 shots 21.0 +- 11.2; 300 shots 8.9 +- 3.3; 1000 shots 5.2 +- 1.0;
3000 shots 4.3 +- 0.1; 10000 shots 4.2 +- 0.05. The error plateaus near 4.2 kcal/mol: the circuit rarely produces some determinants that carry energy.

Noisy sampling (illustrative depolarizing noise, 0.5% per CNOT, 0.05% per one-qubit gate, 1,000 shots, Aer statevector trajectories):
fraction of shots with the right electron count 33% (1 layer), 17% (2 layers); the simple estimate (1-e)^G predicted 17% and 4.5% clean runs
and the observed Hartree-Fock-bitstring share relative to ideal was 20% and 5%.

Control that changes the interpretation: 100 uniformly random bitstrings with the right electron count already contain all 35 possible
up-spin patterns, so the subspace is all 1,225 determinants and SQD returns the exact energy (error < 1e-6 kcal/mol). With noisy circuit samples
plus configuration recovery the subspace also grows to 900-1,225 determinants and the error is 0.0-0.2 kcal/mol. Both are exact for the same
trivial reason. At this size SQD does not test whether the quantum circuit is useful; a fair test needs a matched subspace size (next step).

## 2026-10-09 - matched-size comparison (step 4.5)

Question: at the same number k of electron patterns per spin, does the LUCJ circuit pick better patterns than other methods?
Answer for the 14-qubit N fragment, ideal (noise-free) circuit unless stated:
- Much better than random patterns (with Hartree-Fock forced in, 33 kcal/mol at k = 12 against 4.2) and than a simple
  excitation ordering at k = 4 to 16 (for example 4.8 against 10.5 at k = 8).
- Not better than a classical ranking by the CCSD wavefunction, which needs no circuit: CCSD-ranked equals the best-possible
  ordering (ranking by exact ground-state weight) at k = 1-6, 8, 12, 16. The circuit wins narrowly only at k = 8 and 10 (4.8 and 4.4 against 5.7 and 4.5);
  the "oracle" ordering is not a strict lower bound on energy, which is why the circuit can beat it at k = 8.
- With 0.5% two-qubit noise plus configuration recovery the circuit is no better than the simple excitation ordering at k = 8-12 (gap under 1.5 kcal/mol).
- Finite shots: 50 seeds, 20,000 shots each. Median shots to see k distinct strings: 56 (k=4), 166 (k=6), 439 (k=8), 1,373 (k=10),
  10,148 (k=12, only 26 of 50 seeds got there). No seed reached k = 16. The remaining energy sits in strings with probability around 1e-4 or less.

Mistakes avoided on the way: the first version of the study compared the circuit only with uniformly random strings, which looks impressive and
means nothing. The CCSD-ranked baseline was added because the circuit's angles come from CCSD, so any gain over CCSD must come from the circuit itself.
A first run also reported finite-shot values at k > 13 although the circuit never produced that many distinct strings; those are now reported as not reached.

## 2026-10-09 - Komenco run of the 14-qubit N fragment (Lesson 4, step 4.6)

Pipeline demonstration, not evidence. One POST of the LUCJ circuit (14 qubits, 1,428 gates of which 420 CX, 95 kB) to the open
gateway; answer in about 2 s. Two POSTs were made in total because my first ranking code crashed after the answer arrived
(the gateway returns all 2^14 bitstrings, including wrong-electron-count ones with probability near zero; fixed by ignoring them).
Record: `data/komenco_run.json`.

- Gateway probabilities equal our own exact simulation of the same circuit: total variation distance 3e-14. The bitstring convention
  (qubit 0 rightmost) and gate translation are therefore right, and the gateway behaves as an exact classical emulator, as assumed.
- Ranking strings by gateway probability gives the same energies as the local ideal ranking at every k (identical to 1e-6 kcal/mol).
- 10,000 shots drawn from the gateway probabilities, run through the SQD loop, reach 11 distinct strings and 4.2 kcal/mol error,
  the same as the finite-shot simulation at k = 12 in step 4.5.
- Nothing here says anything about quantum hardware: no noise, no queue, no device connectivity.

## 2026-10-09 - Shots and scheduling model (Lesson 6)

Cost model, not a measurement. Gate counts come from real LUCJ circuits transpiled to a CNOT basis with all-to-all connectivity;
hardware routing would raise them (the Komenco request counts 420 CX for the 14-qubit N fragment because its gate set includes
controlled-phase gates; the CNOT-basis count is 684). Hardware numbers (e = 0.5% per CNOT, 250 us rep delay, 3 s per job) are assumptions.
Record: `data/schedule_study.json`.

- Shots for 100 error-free runs: 7.5 million for the 28-qubit whole molecule, 3,084 for the 14-qubit N fragment (about 2,400x fewer).
- Fragment-parallel scheduling is capped by the largest fragment: N / H / H / H with all empty orbitals kept finishes in 103 s on 4 QPUs
  against 112 s on one. Splitting one circuit's shots scales almost linearly (4.0x less the fixed cost), which is why the multi-QPU paper pairs both.
- Limitation: "error-free run" is a pessimistic proxy because configuration recovery uses some noisy shots. Lesson 4 validated the proxy at
  the 0.5% level only against a simulated depolarizing model.

## 2026-10-09 - IBM hardware run of the 14-qubit N fragment (step 4.7)

Dry run first (`data/ibm_dry_run.json`): routed to `ibm_fez`, 1,107 CZ (best of 5 routing seeds; 1,107-1,148), median CZ error 0.28%,
estimate 1.2 s of execution for 4,000 shots, expected clean fraction 4.5% from the device's own errors (optimistic: ignores readout and crosstalk).
Real run: one Sampler job, 4,000 shots, no mitigation, billed 3 quantum-seconds of the 10-minute monthly allowance. Job db4bg7klf4us73c2h6ig.
IBM's job metadata (UTC): created 09:47:10, running from 09:48:17, finished 10:01:23, i.e. about 1 minute queued and 13 minutes between the running and finished timestamps, of which 3 seconds were billed. Details of all three jobs: `docs/IBM_HARDWARE_JOBS.md`.

- Right-electron-count fraction 6.7% against 7.5% for random bits and 19.9% for the 0.5%-per-CNOT simulation. The simulation under-estimated the noise:
  it used the all-to-all gate count (684) rather than the routed one (1,107), and ignored readout error.
- Matched-size errors: k=4 24.4, k=8 10.6, k=12 5.8 kcal/mol, against 41.8, 23.3, 5.8 for random bits with recovery. Signal at small k, none at k>=10.
- Bug fixed on the way: my first random control used 28-bit strings for a 14-qubit problem and gave a meaningless 0.2% right-count rate.
  The SQD loop also crashes on fully random bits (configuration recovery returns malformed batches); handled by recording "failed" for that control only.
- Not done: no repeat, no second device, no error mitigation, no comparison against a hardware-efficient (heavy-hex) LUCJ layout.

## 2026-10-09 - Heavy-hex circuit, gate-count scaling, whole molecule (steps 4.8-4.9)

Jobs: db4bs3qmb58s738904i0 (four circuits x 4,000 shots, 6 quantum-seconds), db4c58svf2bc73cude20 (whole molecule, 10,000 shots, 5 quantum-seconds).
Total spent across all three jobs: 3 + 6 + 5 = 14 of 600 quantum-seconds.

- Heavy-hex pairs plus ffsim's pre-init passes cut routed CZ from 1,107 to 161 for the 14-qubit fragment; ideal quality barely changes beyond k = 8.
- Right-count fraction: 36.8% (161 CZ), 10.9% (485), 6.7% (1,107). The monotone rise is the one clear hardware gain.
- It did not translate into better subspaces: hardware heavy-hex errors at k = 8, 12, 16 (10.1, 8.7, 6.6) are about the simple classical ordering (10.5, 8.4, 3.6).
  Configuration recovery pulls samples toward Hartree-Fock-like strings, so a classical low-excitation list does about as well.
- Whole molecule, 28 qubits: hardware SQD 3.4 kcal/mol at 314 strings, random bits through the same loop 4.5, classical ordering 0.8. The ideal circuit is itself
  worse than the classical orderings at large k, so noise is not the only problem: with one layer built from CCSD amplitudes, the circuit has no information a CCSD ranking lacks.
- Mistakes on the way: a detached launch died with its shell and submitted nothing (verified: no job id was written) and was rerun properly;
  the random control first had the wrong bit length (28 for a 14-qubit problem); two of my shell commands failed to parse and were redone from script files.
- Not done: error mitigation, a second device, repeats of the whole-molecule job, a 22-qubit fragment, circuits that do not come from CCSD.

## 2026-10-09 - Lesson 7.1: water cluster references

Design note: `docs/LESSON7_DESIGN.md`. Code: `src/fragsim/water.py`. Record: `data/water_reference.json`.

- Dimer interaction energy (CCSD(T), 6-31G, raw) -6.59 kcal/mol; counterpoise -4.62; Hartree-Fock gives 93% of the raw value.
- First ring geometry was unbound (trimer HF +0.35); refitted (rigid monomers, HF energy minimum) to trimer -18.82 and tetramer -32.64 kcal/mol (CCSD(T), raw).
- Monomer FCI minus CCSD(T) = -0.329 kcal/mol triggered a pre-written rule that would have biased the test; amended before any fragment result (design note, section 11).
  STO-3G check: CCSD(T) vs FCI differ by 0.008 kcal/mol in the dimer interaction energy.
- Not done: any fragment calculation; the code risks listed in the design note (core-orbital localization, NH3-specific partitions) are still open for 7.2.

## 2026-10-09 - Lesson 7.2: pipeline generalised to water

- `build_atom_basis` matches the frozen-core space instead of orbital by orbital (needed in principle for two equivalent oxygen cores). I had predicted this was a bug;
  a direct test (core pair rotated by 0, 20, 45 degrees) showed the old rule still worked on this dimer, so it is recorded as hardening, not a fix.
- New identities verified: single water = FCI (1.6e-8 hartree), dimer with zero empty orbitals = Hartree-Fock (1e-8), core-pair rotation leaves the basis unchanged.
- One test tolerance was loosened from 1e-8 to 1e-7 hartree because the two FCI code paths stop at different convergence thresholds; the observed gap was 1.6e-8.
- All 80 NH3 tests unchanged. Total now 99 tests.
- Not done: independent cross-check of the dimer (7.3), any cut-error number (7.4).

## 2026-10-09 - Lesson 7.3: Vayesta cross-check on water

- Probing cluster sizes first showed that all-electron 6-31G dimer clusters are 34 qubits (3.8e8 determinants, 3 GB per vector) and frozen-core ones 32 qubits:
  the NH3-style exact cross-check is impossible here. Done in STO-3G instead (22 to 26 qubit clusters): dimer 0.056 and trimer ring 0.012 kcal/mol from Vayesta.
- This also corrects the design note: a water fragment with its bath is already 32 qubits untruncated, not 24. Tier T1 will be measured on the truncation curve.
- The STO-3G gaps are smaller than the NH3 ones (0.17-0.34); consistent with, but not proof of, the earlier guess that the NH3 gap comes from atom-local virtual functions.
- Not done: any 6-31G cross-check, the cut-error curve (7.4).

## 2026-10-09 - Lesson 7.4: water dimer cut-error curve (interim; superseded by the closing entry below)

- Primary E_int error vs CCSD(T): 171 (k=0), 25 (k=4), 15 (k=6), 6.3 (k=8), 4.8 kcal/mol (k=9, 30-qubit fragments). Tier T1 (<= 1 kcal/mol) not met; hard stop in force.
  The untruncated k=10 point (32 qubits) is running; the stated expectation is 3 to 4 kcal/mol.
- The secondary (matched-k monomer) reading was added before the run to remove the correlation mismatch; it did not behave (non-monotone, 12 kcal/mol at k=6) because the
  two sides have different numbers of occupied orbitals at equal k. Recorded as an approach that did not work.
- STO-3G control: untruncated cut error 0.13 kcal/mol, so the energy partition and bath are sound; the 6-31G error is a feature of the scheme in a basis with empty functions (hypothesis untested).
- Timing: k=8 179 s, k=9 649 s for two clusters; k=10 expected about 45 minutes.
- Not done: k=10, any SQD/hardware for water (blocked by the stop), the options listed in section 15.

## 2026-10-09 - Lesson 7 closed: option 1 (larger fragments) and what could not be run

- Overlapping fragments collapse to the whole dimer for two molecules, so option 1 was run as an unequal split (water + donated H against the rest). Predicted worse; it was:
  16.5 against 9.95 kcal/mol at 26 qubits, 12.9 against 6.3 at 28.
- Two jobs did not complete. The equal-waters k = 10 run (32 qubits) was stopped by the system for low memory; the unequal split's k = 8 exited with no output and no row (probably out of
  memory). Neither was retried. Consequence: tier T1 was never measured at its pre-registered untruncated point; the verdict rests on the curve through 30 qubits.
- Left running a second large job while the first was alive on a 7 GB machine; the memory kill followed. Lesson: run one large FCI at a time here.
- Hard stop respected: no SQD, no hardware, no quantum seconds used in Lesson 7.
- Not done: k = 10, unequal split at 30 qubits, trimer and tetramer cut errors, a correlated bath, FMO-style pair corrections, a direct test of the empty-orbital hypothesis.
