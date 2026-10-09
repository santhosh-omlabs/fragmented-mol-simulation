# IBM hardware jobs

Every job submitted to IBM hardware in this study, and every read-only call made to the service. Three jobs, all on `ibm_fez`, all on 2026-10-09, **14 quantum-seconds billed out of a
600-second monthly allowance (2.3%)**. All numbers below come from the committed records in [../data/](../data/) and from IBM's job metadata (timestamps in UTC).
The jobs are described as hardware stamps: one device, one day, no repeats, no error mitigation. They are not a hardware result in the usual sense.

## Setup common to all jobs

| Item | Value |
|---|---|
| Backend | `ibm_fez`, 156 qubits, native two-qubit gate CZ (read from the backend target) |
| Median CZ error reported by the device on the day | 0.28% (`data/ibm_dry_run.json`) |
| Primitive | `SamplerV2` in job mode, default options: no dynamical decoupling, no gate twirling, no error mitigation |
| Software | Qiskit 2.5.2, qiskit-ibm-runtime 0.50.0, ffsim 0.0.84, qiskit-addon-sqd 0.14.0, PySCF 2.14.0 |
| Circuit family | LUCJ ansatz built from classical CCSD amplitudes of the cluster; Hartree-Fock state, then the LUCJ layers, then measure every qubit |
| Credentials | Saved once in the local Qiskit account store; never in this repository |

## Job 1 - the 14-qubit nitrogen fragment

| Item | Value |
|---|---|
| Job id | `db4bg7klf4us73c2h6ig` |
| Submitted by | `python -m fragsim.ibm_run` |
| Circuit | N fragment of NH3 (rung 1, 3 empty orbitals kept): 7 orbitals, 4 + 4 electrons, 14 qubits; LUCJ with 2 layers, all-to-all interactions |
| Compilation | `transpile`, optimization level 3, routing seed 1 (the best of five in the dry run: 1,139, 1,107, 1,148, 1,137, 1,114 CZ for seeds 0-4) |
| Routed size | 1,107 CZ (684 CNOT with all-to-all connectivity) |
| Shots | 4,000 |
| Timestamps (UTC) | created 09:47:10, running from 09:48:17, finished 10:01:23 |
| Billed | 3 quantum-seconds |
| Dry-run estimate | about 1.2 s of execution; billing is higher and I did not investigate why |
| Records | `data/ibm_run.json`, `data/ibm_counts_n_fragment.json`, `data/ibm_dry_run.json`, `data/ibm_analysis.json` |

Outcome. 3,527 distinct bitstrings in 4,000 shots. Shots with the right electron count in both spin halves: **6.7%**, against 7.5% expected for uniformly random bits
(`(35/128)^2`) and 19.9% in the 0.5%-per-CNOT noise simulation, which used 684 gates and no readout error, so it under-estimated the noise. Matched-size SQD error (kcal/mol,
strings chosen from the samples after configuration recovery): k = 4: 24.4, k = 8: 10.6, k = 12: 5.8, against 41.8, 23.3, 5.8 for random bits. Some signal at small k, none at k of 10 or more.

## Job 2 - gate count against sample quality

| Item | Value |
|---|---|
| Job id | `db4bs3qmb58s738904i0` |
| Submitted by | `python -m fragsim.ibm_scaling` |
| Circuits | four circuits in one job, 4,000 shots each, same 14-qubit fragment: heavy-hex-restricted LUCJ, 1 layer, three routing seeds (1, 4, 0); all-to-all LUCJ, 1 layer, routing seed 3 |
| Compilation | preset pass manager, optimization level 3, with ffsim's pre-initialization passes |
| Routed size (CZ) | 161, 161, 167 (heavy-hex); 485 (all-to-all) |
| Timestamps (UTC) | created 10:12:31, running from 10:17:17, finished 10:17:39 (about 4 minutes 46 seconds queued) |
| Billed | 6 quantum-seconds |
| Records | `data/ibm_scaling_run.json`, `data/ibm_scaling_analysis.json`, `data/hwfriendly_study.json` |

| Circuit | CZ | Shots with the right electron count | Error at k = 8 / 12 / 16 (kcal/mol) |
|---|---:|---:|---|
| heavy-hex, seed 1 | 161 | 41.7% | 10.3 / 9.2 / 7.9 |
| heavy-hex, seed 4 | 161 | 38.7% | 9.8 / 7.9 / 4.1 |
| heavy-hex, seed 0 | 167 | 30.0% | 10.3 / 8.9 / 7.9 |
| all-to-all, 1 layer | 485 | 10.9% | 9.8 / 7.9 / 4.1 |

Heavy-hex mean 36.8% with a standard deviation of 6.1 (three layouts). Together with job 1 (1,107 CZ, 6.7%) this gives the gate-count scaling quoted in the README.
Cleaner samples did not give better subspaces: the errors above are about those of a simple classical ordering (10.5 / 8.4 / 3.6).

## Job 3 - the whole 28-qubit molecule

| Item | Value |
|---|---|
| Job id | `db4c58svf2bc73cude20` |
| Submitted by | `python -m fragsim.whole_molecule submit` |
| Circuit | whole NH3: 14 orbitals, 4 + 4 electrons, 28 qubits; heavy-hex-restricted LUCJ, 1 layer |
| Compilation | preset pass manager, optimization level 3, ffsim pre-initialization passes, routing seed 1 |
| Routed size | 668 CZ. Seed 1 is fixed in the script; the dry run showed 596, 668, 608, 571, 648 CZ for seeds 0-4, so the job did not use the best seed (571). I did not revisit this. |
| Shots | 10,000 |
| Timestamps (UTC) | created 10:32:03, running from 10:32:04, finished 10:32:11 |
| Billed | 5 quantum-seconds |
| Records | `data/whole_molecule.json`, `data/whole_molecule_counts.json`, `data/whole_molecule_dry.json` |

Outcome. Right electron count in 4.0% of shots (uniform random bits: 0.38%). After configuration recovery 985 distinct up/down strings remained; the SQD loop converged to
**314 strings per spin (98,596 of 1,002,001 determinants) and an error of 3.43 kcal/mol** against exact FCI (rounds: 10.67, 4.45, 3.43, 3.97, 3.86).
Controls at the same 314 strings: the same loop on random bits 4.47 kcal/mol; a plain classical excitation ordering 0.78; a CCSD-ranked choice 0.78.
The ideal noiseless circuit would itself be worse than both classical orderings at k of 64 or more (25.5 against 10.6 and 4.0 at k = 64), so noise is not the only limit.

## Totals

| | Billed (s) | CZ gates (routed) | Shots |
|---|---:|---:|---:|
| Job 1 | 3 | 1,107 | 4,000 |
| Job 2 | 6 | 161 / 161 / 167 / 485 | 4 x 4,000 |
| Job 3 | 5 | 668 | 10,000 |
| **Total** | **14** | | 30,000 |

IBM also reported resource-unit charges of 0.0481, 0.0995 and 0.0736 for the three jobs; I did not look into how those relate to quantum-seconds.

## Calls that used no quantum time

- Listing operational backends (`ibm_fez`, `ibm_kingston`, `ibm_marrakesh`) and their pending-job counts, to choose a device.
- Reading backend properties to transpile for the device: `python -m fragsim.ibm_dry_run`, `python -m fragsim.hwfriendly`, `python -m fragsim.whole_molecule dry`.
- Reading job metadata afterwards (the timestamps above).

## One launch that submitted nothing

The first attempt to start job 3 used a detached shell that ended with its parent; no job id was written, and I confirmed that no job existed before running it again. Job 1's command was
started in the foreground, exceeded my local 10-minute timeout, and kept running; it is one job, not two.

## What these jobs do not show

No repeats of any job; one device on one day; no error mitigation; circuits not chosen to be best for the chip beyond the heavy-hex restriction; for job 3 not the best of five routing seeds.
Nothing here is evidence of quantum advantage, and the hardware SQD result on the whole molecule is worse than a free classical ordering.
To run any of it yourself, see [../GUIDE.md](../GUIDE.md) section 6. A rerun will not reproduce these numbers exactly.
