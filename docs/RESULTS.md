# Results

Every number in this study, with the command or data file that produced it. The short version is in [../README.md](../README.md); what to take away
from it is in [TAKEAWAYS.md](TAKEAWAYS.md). Energies in hartree unless a table says kcal/mol (chemical accuracy is 1 kcal/mol).
All results are from classical simulation unless a section says IBM hardware. Records are the JSON files in [../data/](../data/).

## Lesson 1: the classical baseline

NH3 pyramidal, N-H = 1.012 angstrom, H-N-H = 106.7 degrees, 6-31G, N 1s frozen. Energies in hartree.

| Quantity | Value | Meaning |
|---|---:|---|
| Spatial orbitals (total / active) | 15 / 14 | basis functions; the N 1s core is frozen |
| Active electrons | 8 | |
| Qubits for the whole molecule | 28 | 2 per spatial orbital (spin up and down) |
| FCI determinants | 1,002,001 | ways to place 4 up and 4 down electrons in 14 orbitals |
| E_HF | -56.161063 | mean-field energy |
| E_FCI | -56.291514 | exact in this basis and frozen-core space (classical, about 12 s) |
| E_FCI - E_HF | -0.130451 (-81.9 kcal/mol) | correlation energy: what a quantum or fragment method must recover |

`python -m pytest -q` asserts these values.

### Fragmentation error with exact fragment solvers (Lesson 3, classical only)

Each fragment plus its bath is solved by exact FCI (no quantum step yet). "Empty orbitals kept" is how many MP2 natural virtual
orbitals each fragment keeps; error is the total energy minus E_FCI, in kcal/mol (chemical accuracy is 1). The largest
fragment's qubit count is 2 x (cluster orbitals). Record: `data/ladder_exact_fragments.json`.

| Fragments | Empty orbitals kept | Largest fragment (qubits) | Error vs E_FCI (kcal/mol) |
|---|---|---:|---:|
| none (Hartree-Fock) | 0 | - | 81.9 |
| N / H / H / H | 3 | 14 | 38.8 |
| N / H / H / H | 5 | 18 | 21.3 |
| N / H / H / H | all | 22 | 17.3 |
| N+H1 / H / H | all | 24 | 11.7 |
| N+H1+H2 / H | all | 26 | 5.9 |
| whole molecule | all | 28 | 0 (by definition) |

Regenerate: `python -m fragsim.ladder` (about 4 minutes). With every fragment at or under Komenco's 15-qubit cap, the best total error
in the table is 38.8 kcal/mol, about half the Hartree-Fock error. Errors are positive: fragments recover less correlation than exact.

### Independent cross-check against Vayesta

The same three partitions, all-electron (Vayesta cannot freeze the N 1s core the way this repo does), exact FCI on every
cluster, DMET bath, no truncation. Energies in hartree. Record: `data/crosscheck_vayesta.json`, Vayesta commit 7f1639d.

| Fragments | fragsim | Vayesta (EWF) | difference (kcal/mol) |
|---|---:|---:|---:|
| N / H / H / H | -56.26404822 | -56.26459310 | 0.34 |
| N+H1 / H / H | -56.27326027 | -56.27371184 | 0.28 |
| N+H1+H2 / H | -56.28276578 | -56.28303551 | 0.17 |

The two agree to about 0.3 kcal/mol, which is 2% of the fragmentation error itself. fragsim recovers slightly less correlation in
every case. The cause was not investigated; a likely candidate is that the two codes build the empty fragment functions (PAOs)
slightly differently. The frozen-core path and the truncation are not covered by this check; they are covered by the
whole-molecule and zero-virtual tests.

### SQD on the 14-qubit N fragment: does the circuit beat classical choices? (Lesson 4, classical simulation only)

The N fragment of rung 1 with 3 empty orbitals kept is a 14-qubit problem with 1,225 determinants (7 orbitals, 4 up and 4 down
electrons). SQD builds its subspace from electron patterns per spin ("strings"; there are 35 possible); with k strings the subspace has
k x k determinants. The table fixes k, so every method gets the same subspace size and differs only in WHICH strings it picks.
Entries are the SQD energy error against exact FCI of this cluster, in kcal/mol. Record: `data/matched_subspace_study.json`.

| k (determinants) | random + HF | excitation order | CCSD-ranked (classical) | ideal circuit | circuit, finite shots (shots to reach k) | noisy circuit + recovery |
|---|---:|---:|---:|---:|---:|---:|
| 4 (16) | 42.1 | 38.2 | **12.3** | **12.3** | 24.3 (56) | 24.5 |
| 6 (36) | 38.4 | 26.4 | **7.5** | 9.5 | 9.5 (166) | 11.4 |
| 8 (64) | 36.0 | 10.5 | 5.7 | **4.8** | 6.2 (439) | 10.0 |
| 12 (144) | 32.8 | 8.4 | **2.9** | 4.2 | 4.2 (10,148; 26 of 50 seeds) | 8.6 |
| 16 (256) | 27.7 | 3.6 | **0.7** | 2.7 | not reached in 20,000 shots | 6.8 |
| 20 (400) | 19.6 | 0.7 | **0.3** | 0.7 | not reached | 0.7 |
| 35 (1,225) | 0 | 0 | 0 | 0 | 0 | 0 |

Random patterns without the Hartree-Fock pattern are off by thousands of kcal/mol and are not shown. "CCSD-ranked" needs no circuit:
it ranks strings by their weight in the same classical coupled-cluster wavefunction that sets the circuit's angles. The noisy column uses
an illustrative depolarizing model (0.5% per two-qubit gate, 1,000 shots), not a real device.

What this shows, for this fragment: the ideal circuit is much better than random choice and than a simple excitation ordering at k = 4 to 16,
but it is not better than the classical CCSD-ranked choice, which matches the best possible selection at most k. Noise removes the
circuit's edge over the simple classical ordering (k = 8 to 12). The ideal circuit also needs many shots to produce rare strings
(about 440 for k = 8, about 10,000 for k = 12) and cannot reach k = 16 in 20,000 shots. An earlier, simpler comparison against uniformly
random strings (without the matched-size control) made SQD look exact; that was only because 100 random samples cover all 35 strings.

### Komenco gateway run (pipeline demonstration)

The same 14-qubit circuit (1,428 gates, 420 of them CX) was sent to the Komenco gateway, a third-party platform provided by Automatski (https://automatski.com/platform.html). Its probabilities equal our own exact
simulation to 3e-14 (total variation distance), so the gateway acts as an exact classical emulator and the client, gate translation and
bit ordering are verified. 10,000 shots drawn from them and run through SQD reach 11 of the 35 strings and 4.2 kcal/mol error, the same as the
finite-shot simulation above. This demonstrates the workflow, not hardware behaviour. Record: `data/komenco_run.json`.

### Why fragment on hardware, and how much does a second QPU help? (Lesson 6, cost model only)

A circuit with G two-qubit gates, each failing with probability e = 0.5%, runs error-free with probability (1-e)^G. Seeing 100 error-free runs
therefore takes about 100 / (1-e)^G shots. Gate counts are measured on our real LUCJ circuits (CNOT basis, all-to-all connectivity, 2 layers,
no routing, which would add more). Timing assumes 250 us between shots and 3 s fixed cost per job. These are assumptions, not measurements.
Record: `data/schedule_study.json`.

| Fragment set | Largest job (qubits / CX) | Shots for 100 clean runs (largest job) | Total QPU time | Wall time, 4 QPUs | Energy error (kcal/mol) |
|---|---|---:|---:|---:|---:|
| whole molecule | 28 / 2,240 | 7.5 million | 52 min | 13 min (shots split 4 ways) | 0 by definition (not an SQD result) |
| N / H / H / H, 3 empty kept | 14 / 684 | 3,084 | 13 s | 4 s | 38.8 |
| N / H / H / H, all empty kept | 22 / 1,568 | 259,066 | 112 s | 103 s | 17.3 |
| N+H1+H2 / H, all empty kept | 26 / 2,000 | 2.3 million | 15 min | 15 min | 5.9 |

What it shows: cutting the gate count cuts the shots needed exponentially, which is the real benefit of fragmenting (the whole molecule would not fit a
10-minute budget on one QPU in this model). A second QPU helps little when one fragment dominates: with N / H / H / H the finish time is set by the
N job, so 4 QPUs give only 1.09x (all empty kept) while splitting one circuit's shots over 4 QPUs gives 4.0x minus the per-job overhead. The three H fragments are
identical by C3v symmetry, so one run could stand in for all three. The "clean run" criterion is pessimistic: SQD tolerates some noisy shots through
configuration recovery (Lesson 4), so real shot needs are lower, but by an amount not measured here.

### IBM hardware run: one job, 14-qubit N fragment

One Sampler job on `ibm_fez`, 4,000 shots, raw (no error mitigation, no dynamical decoupling, no twirling), billed 3 quantum-seconds. The circuit is the
same LUCJ circuit as above, routed to the chip: 1,107 CZ gates (684 CNOT with all-to-all connectivity). Job `db4bg7klf4us73c2h6ig`. Records:
`data/ibm_run.json`, `data/ibm_counts_n_fragment.json`, `data/ibm_analysis.json`; dry-run estimate in `data/ibm_dry_run.json`.

| Source of 4,000 bitstrings | Shots with the right electron count | Matched-size error at k = 4 | k = 8 | k = 12 |
|---|---:|---:|---:|---:|
| ideal circuit (exact) | 100% | 12.3 | 4.8 | 4.2 |
| noisy simulation (0.5% per CNOT, 1,000 shots) | 19.9% | 24.4 | 10.0 | 8.6 |
| **IBM `ibm_fez`** | **6.7%** | **24.4** | **10.6** | **5.8** |
| uniformly random bits (control) | 7.4% | 41.8 | 23.3 | 5.8 |

(errors in kcal/mol against exact FCI of the cluster; the last three rows use configuration recovery on the samples.)

What it shows, and does not show:

- The hardware output is mostly noise. Only 6.7% of shots have the right electron count, about what random bits give (7.5%), and well below the
  noisy simulation (19.9%). My 0.5% per-CNOT simulation was too optimistic for this routed circuit: it had 684 gates, the chip ran 1,107.
- Some signal survives. At k = 4 to 8 the hardware-derived strings are clearly better than random bits (24.4 against 41.8 at k = 4; 10.6 against 23.3 at k = 8)
  and match the noisy simulation. At k >= 10 the hardware is indistinguishable from random bits plus configuration recovery (5.8 against 5.8 at k = 12).
- The SQD energy from the full loop is exact (error about 1e-9) for hardware, simulation and random bits alike, because each gets all 35 strings.
  That number says nothing about the hardware; it is the same trap as in Lesson 4. The matched-size columns are the fair comparison.
- A single job on one device on one day. No error bars, no repeat. It says nothing about quantum advantage.

### A circuit built for the chip: fewer gates, cleaner samples (second IBM job)

Restricting the Jastrow part of the circuit to heavy-hex-friendly pairs (`heavy_hex_pairs`) and using one layer cuts the routed gate count from 1,107 to 161 CZ
(`data/hwfriendly_study.json`). The ideal circuit loses little: same error at k = 12 (4.2), 6.2 instead of 4.8 at k = 8. One job on `ibm_fez` with four circuits,
4,000 shots each, billed 6 quantum-seconds (`data/ibm_scaling_run.json`, `data/ibm_scaling_analysis.json`):

| Circuit | CZ gates | Shots with the right electron count | Error at k = 8 / 12 / 16 (kcal/mol) |
|---|---:|---:|---|
| heavy-hex, 1 layer (3 layouts) | 161-167 | 36.8% +- 6.1 (std of 3) | 10.1 / 8.7 / 6.6 (std 0.3 / 0.7 / 2.2) |
| all-to-all, 1 layer | 485 | 10.9% | 9.8 / 7.9 / 4.1 |
| all-to-all, 2 layers (first job) | 1,107 | 6.7% | 10.6 / 5.8 / 2.4 |
| simple classical ordering (no circuit) | - | - | 10.5 / 8.4 / 3.6 |

Sample quality rises steeply as the gate count falls (6.7% to 37%), as the (1-e)^G rule predicts. But better samples did NOT give better subspaces:
with configuration recovery, the strings chosen from hardware samples are no better than the simple classical ordering at any k tested.
Only 3 layouts, one device, one day.

### Closing the three-way comparison: the whole 28-qubit molecule on IBM (third job)

The whole molecule (28 qubits, no fragmentation) with the same pipeline: heavy-hex LUCJ, 1 layer, 668 CZ, 10,000 shots, 5 quantum-seconds
(`data/whole_molecule.json`, job `db4c58svf2bc73cude20`). Right electron count in 4.0% of shots (uniform random bits: 0.4%). SQD with configuration
recovery on the hardware samples ends at 314 strings per spin, 98,596 of 1,002,001 determinants. Errors against exact FCI (-56.291514), in kcal/mol:

| Method (whole molecule, 28 qubits, 314 strings per spin) | Error |
|---|---:|
| SQD on IBM hardware samples | 3.4 |
| the same SQD loop on uniformly random bits | 4.5 |
| simple classical ordering, no circuit | **0.8** |
| CCSD-ranked, no circuit | **0.8** |

The hardware is barely better than random bits and four times worse than a free classical ordering. Even the ideal noiseless circuit would be worse than
both classical orderings at k >= 64 (25.5 against 10.6 and 4.0 at k = 64). The whole molecule is also easy classically (FCI in about 12 s), so this is a pipeline
stamp, not a result about hardware usefulness.

### Lesson 7: a weakly bound system, the water dimer (negative result, no quantum time used)

NH3 was cut through covalent bonds. Water dimers are held by a hydrogen bond, which should be a cheaper cut. Design and pass/fail criteria were written before any fragment
calculation (`docs/LESSON7_DESIGN.md`). Setting: 6-31G, frozen O 1s, model geometry, one fragment per water, exact FCI on every cluster. Reference: frozen-core CCSD(T),
raw interaction energy -6.59 kcal/mol (counterpoise-corrected -4.62; Hartree-Fock gives 93% of the raw value). Errors are in the interaction energy, kcal/mol:

| Largest cluster (qubits) | 20 | 22 | 24 | 26 | 28 | 30 |
|---|---:|---:|---:|---:|---:|---:|
| two whole waters | 25.0 | 17.5 | 15.2 | 9.95 | 6.3 | **4.8** |
| water + donated H against the rest (covalent cut) | 70.4 | 43.3 | 28.0 | 16.5 | 12.9 | not obtained |

Records: `data/water_reference.json`, `data/water_dimer_cut.json`, `data/water_dimer_cut_unequal.json`, `data/water_sto3g_cut.json`, `data/crosscheck_vayesta_water.json`.

- The pre-registered target (1 kcal/mol with exact solvers) is not met at any computed size, so by the rule fixed beforehand no SQD and no hardware were run on water.
  The untruncated point (32 qubits) could not be measured: the job was stopped for low memory on a 7 GB machine.
- The cut itself is not the problem in a minimal basis: in STO-3G it costs 0.13 kcal/mol against an exact interaction of -4.99.
  A hypothesis, untested: in 6-31G the neighbouring water's empty orbitals carry correlation that a bath built from the occupied density cannot represent.
- Moving the cut onto a covalent O-H bond is worse at every matched size.
- Our embedding agrees with Vayesta to 0.06 kcal/mol (dimer) and 0.01 (trimer ring) in STO-3G, all-electron, exact clusters.
- This is not evidence against fragment methods for water clusters in general: FMO-type methods with pair corrections are designed for them. It shows that one-shot embedding with
  exact clusters of up to 30 qubits does not reach chemical accuracy here.

### The original question, answered for this molecule

| Route to E_FCI | Largest circuit | Error vs E_FCI (kcal/mol) |
|---|---|---:|
| fragments, exact solvers (N / H / H / H, 3 empty kept) | 14 qubits | 38.8 |
| fragments, exact solvers (N / H / H / H, all kept) | 22 qubits | 17.3 |
| fragments, exact solvers (N+H1+H2 / H, all kept) | 26 qubits | 5.9 |
| whole molecule, SQD on IBM | 28 qubits | 3.4 |
| whole molecule, classical ordering (no QPU) | 28 qubits | 0.8 |

Fragmenting NH3 saves little (each +2 qubits buys about 6 kcal/mol), the circuit adds nothing a classical ordering lacks, and today's noise removes most of the
rest. What fragmenting does buy is exponentially fewer shots (Lesson 6). These are the conclusions for this size; they say nothing about
molecules too large for classical methods. The fragment rows use exact solvers, not SQD, so they show the cost of the cut alone.
