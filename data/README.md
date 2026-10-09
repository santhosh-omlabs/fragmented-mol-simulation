# data/

Small JSON records, each the output of one command from [../GUIDE.md](../GUIDE.md). The README and [../docs/RESULTS.md](../docs/RESULTS.md) quote these files; the tests assert their headline values.
Records from IBM hardware depend on the device, its calibration and the day, so rerunning will not reproduce them exactly. Large regenerable files are gitignored.

| File | Produced by | What it holds |
|---|---|---|
| `ladder_exact_fragments.json` | `python -m fragsim.ladder` | NH3 fragment ladder: error against exact FCI for every partition and truncation |
| `crosscheck_vayesta.json` | `scripts/crosscheck_vayesta.py` | fragsim against Vayesta, NH3, all-electron |
| `matched_subspace_study.json` | `python -m fragsim.matched` | SQD error at matched subspace size for six ways of choosing electron patterns |
| `noisy_counts_n_fragment_2layers.json` | `python -m fragsim.matched` (cached input) | 1,000 simulated noisy shots of the 14-qubit fragment circuit |
| `komenco_run.json` | `python -m fragsim.komenco_run` | one run of the 14-qubit circuit on the Komenco gateway (Automatski, https://automatski.com/platform.html) |
| `schedule_study.json` | `python -m fragsim.schedule` | shots and QPU time per fragment set (cost model with assumed hardware numbers) |
| `ibm_dry_run.json` | `python -m fragsim.ibm_dry_run` | routed gate count and time estimate for the 14-qubit fragment (no job) |
| `ibm_run.json`, `ibm_counts_n_fragment.json`, `ibm_job_id.txt` | `python -m fragsim.ibm_run` | IBM job 1: raw counts and billed time |
| `ibm_analysis.json` | `python -m fragsim.ibm_analysis` | job 1 against the noisy simulation and random bits |
| `hwfriendly_study.json` | `python -m fragsim.hwfriendly` | all-to-all against heavy-hex-restricted circuits: routed gates and ideal quality |
| `ibm_scaling_run.json`, `ibm_scaling_job_id.txt` | `python -m fragsim.ibm_scaling` | IBM job 2: four circuits, raw counts |
| `ibm_scaling_analysis.json` | `python -m fragsim.ibm_scaling_analysis` | valid-shot fraction and matched-size errors against gate count |
| `whole_molecule_dry.json` | `python -m fragsim.whole_molecule dry` | ideal 28-qubit circuit quality and routed gates (no job) |
| `whole_molecule_counts.json`, `whole_molecule_job_id.txt`, `whole_molecule.json` | `python -m fragsim.whole_molecule submit`, then `analyse` | IBM job 3: raw counts, SQD result and controls |
| `water_reference.json` | `python -m fragsim.water` | water dimer, trimer ring, tetramer ring: CCSD(T) interaction energies, geometries, monomer FCI check, STO-3G check |
| `crosscheck_vayesta_water.json` | `scripts/crosscheck_vayesta_water.py` | fragsim against Vayesta on water, STO-3G |
| `water_sto3g_cut.json` | `scripts/water_sto3g_cut.py` | untruncated cut error of the water dimer in STO-3G against exact FCI |
| `water_dimer_cut.json` | `python -m fragsim.water_cut --kmax=9` | interaction-energy error against clusters of 12 to 30 qubits (two whole waters) |
| `water_dimer_cut_unequal.json` | `python -m fragsim.water_cut --partition=hbond-shift --kmax=7` | the same for the unequal split (water plus donated hydrogen) |

IBM job ids are included as identifiers for the hardware runs; every job is described in [../docs/IBM_HARDWARE_JOBS.md](../docs/IBM_HARDWARE_JOBS.md). No credentials are stored here or anywhere in the repository.
