# Guide - reproducing the study

Every command below was run on the author's machine (Windows 11, WSL Ubuntu, Python 3.12, 12 cores, 7 GB RAM). The text after "Expected" is what that run printed.
Commands that submit work to an outside service are marked **EXTERNAL**; nothing in this guide runs one unless you type that command.

## 0. Conventions

- Commands assume a Linux or macOS shell, or the Ubuntu shell inside WSL on Windows (PySCF has no Windows wheels). Run them from the repository root.
- `PY` stands for `~/venvs/fragsim/bin/python`. On Windows PowerShell, prefix any command with `wsl --cd $PWD --`, for example `wsl --cd $PWD -- ~/venvs/fragsim/bin/python -m pytest -q`.
- Keep the virtual environment in your Linux home folder, not on a mounted Windows drive: one on `/mnt/d` was corrupted during install on the author's machine.
- Flags use the `--name=value` form.
- Timings are for the author's machine. "Records" are the JSON files in `data/` that the README and [docs/RESULTS.md](docs/RESULTS.md) quote.

## 1. Install (about 5 minutes)

```text
python3 -m venv ~/venvs/fragsim
~/venvs/fragsim/bin/pip install -e ".[dev,quantum,ibm]"
```

Expected: ends without error. To keep pip caches off a small system drive: `PIP_CACHE_DIR=/mnt/d/pip-cache TMPDIR=/mnt/d/tmp` (adjust the paths).

## 2. Verify (about 2 minutes)

```text
PY -m pytest -q
```

Expected: `109 passed` (106 with `-m "not slow"`). The tests recompute RHF and exact FCI for NH3, check the orbital, bath, cluster and energy-share steps against exact identities
(a whole-molecule cluster reproduces FCI, zero empty orbitals reproduces Hartree-Fock), check the Komenco client against a fake transport, and assert the headline numbers in
the committed records. They need no network and no IBM account.

## 3. Regenerate the classical results (no external service)

| Result | Command | Time | Expected | Record |
|---|---|---|---|---|
| NH3 fragment ladder (Lesson 3) | `PY -m fragsim.ladder` | 4 min | 21 lines such as `1: N \| H1 \| H2 \| H3    k=None  error 17.27 kcal/mol` | `data/ladder_exact_fragments.json` |
| Matched-size SQD study (Lesson 4) | `PY -W ignore -m fragsim.matched` | 5 min | 13 lines `k= 1 ...` to `k=35 ...` | `data/matched_subspace_study.json` |
| Shots and scheduling model (Lesson 6) | `PY -W ignore -m fragsim.schedule` | 2 min | 5 lines, first `whole molecule  largest 28q 2240 CX \| shots for 100 clean runs: 7.52e+06` | `data/schedule_study.json` |
| Water references (Lesson 7.1) | `PY -W ignore -m fragsim.water` | 3 min | dimer -6.59, trimer ring -18.82, tetramer ring -32.64 kcal/mol (CCSD(T), raw) | `data/water_reference.json` |
| Water STO-3G cut check | `PY -W ignore scripts/water_sto3g_cut.py` | 1 min | `cut_error_kcal_mol` near 0.13 | `data/water_sto3g_cut.json` |
| Water dimer cut curve (Lesson 7.4) | `PY -W ignore -m fragsim.water_cut --kmax=9` | 15 min | error falls from about 171 to 4.8 kcal/mol | `data/water_dimer_cut.json` |
| Unequal split (Lesson 7) | `PY -W ignore -m fragsim.water_cut --partition=hbond-shift --kmax=7` | 6 min | 12.9 kcal/mol at k = 7 | `data/water_dimer_cut_unequal.json` |

Notes. `matched` reuses `data/noisy_counts_n_fragment_2layers.json` (1,000 simulated noisy shots); delete it to resample (about 1 minute more). The scheduling hardware numbers are assumptions in
`CostModel`; edit them to model your own device. `python -m fragsim.water fit` refits the ring geometries (3 minutes) and prints numbers to paste into `RING_PARAMS`.

**Memory warning.** The water cut curve at k = 10 (32-qubit clusters, about 45 minutes) and the unequal split at k = 8 were killed or exited silently on the author's 7 GB machine, and
are not in the records. Do not expect them to fit in 7 GB; run one large FCI at a time.

## 4. Cross-checks against Vayesta (optional, no external service beyond a git clone)

Vayesta is not on PyPI, and installing it from GitHub needs a BLAS development library. It runs from a source checkout:

```text
git clone --depth 1 https://github.com/BoothGroup/Vayesta ~/src/Vayesta
VAYESTA_PATH=~/src/Vayesta PY scripts/crosscheck_vayesta.py            # NH3, about 10 minutes
VAYESTA_PATH=~/src/Vayesta PY -W ignore scripts/crosscheck_vayesta_water.py   # water STO-3G, about 8 minutes
```

Expected: NH3 differences 0.17 to 0.34 kcal/mol (the clone should be at or near commit 7f1639d); water `dimer ... 0.056` and `trimer ring ... 0.012`. Records: `data/crosscheck_vayesta.json`,
`data/crosscheck_vayesta_water.json`. The line `No module named 'dyson'` is harmless. Vayesta prints a long log. All-electron 6-31G water is not attempted: its clusters do not fit in 7 GB.

## 5. EXTERNAL: the Komenco gateway (optional)

Komenco is a platform provided by Automatski (https://automatski.com/platform.html), a provider separate from this repository. This sends one 14-qubit circuit (about 95 kB) over plain HTTP to a shared trial server (15-qubit cap, open key). Do not send anything you would not post publicly.
The endpoint and key can be overridden with `KOMENCO_HOST`, `KOMENCO_PORT`, `KOMENCO_API_KEY`; see `.env.example`.

```text
PY -W ignore -m fragsim.komenco_run --dry-run     # builds and sizes the request; no network
PY -W ignore -m fragsim.komenco_run               # one POST
```

Expected: `request: 14 qubits, 1428 gates (420 CX)`; the gateway answers in a couple of seconds; `total variation distance ... ~1e-14`;
`SQD on gateway samples: 11 strings, error 4.207 kcal/mol`; then `wrote .../data/komenco_run.json`. The gateway returns exact probabilities, so it is treated as a classical emulator.

## 6. EXTERNAL: IBM hardware (optional, uses your own runtime allowance)

You need an IBM Quantum account. Save credentials once, outside the repository, in your own terminal (never paste a token into a file in this repo):

```text
PY -c "from qiskit_ibm_runtime import QiskitRuntimeService as S; S.save_account(channel='ibm_quantum_platform', token='YOUR_TOKEN', instance='YOUR_INSTANCE', set_as_default=True)"
```

Commands that only read backend properties (no job, no quantum time): `fragsim.ibm_dry_run`, `fragsim.hwfriendly`, `fragsim.whole_molecule dry`. They print routed gate counts and an estimate.
Commands that submit exactly one job each:

| Command | What it submits | Billed in the author's run |
|---|---|---|
| `PY -W ignore -m fragsim.ibm_run` | 14-qubit fragment, 4,000 shots (run `fragsim.ibm_dry_run` first to see the estimate) | 3 s |
| `PY -W ignore -m fragsim.ibm_scaling` | four circuits, 4,000 shots each | 6 s |
| `PY -W ignore -m fragsim.whole_molecule submit` | 28-qubit whole molecule, 10,000 shots | 5 s |

Analysis commands (no job): `fragsim.ibm_analysis`, `fragsim.ibm_scaling_analysis`, `fragsim.whole_molecule analyse`. Each submit command writes the job id to `data/` before waiting;
if interrupted, fetch results with `--job-id=<id>` (`ibm_run`, `ibm_scaling`) or `analyse --job-id=<id>` (`whole_molecule`). Queue time is not billed. Hardware results differ run to run
(noise, calibration, layout), so a rerun will not reproduce the committed numbers exactly; the tests check the committed records. The scripts target the backend `ibm_fez`
(`--backend=<name>` to change it); it may be retired or renamed by the time you read this.

## 7. Known problems

- `Failed to start the systemd user session` printed by WSL is harmless.
- If WSL stops starting after a memory kill ("distribution failed to start"), run `wsl --shutdown` from Windows and retry; the virtual environment survives.
- Do not commit `.env` files or tokens. `.env.example` lists only optional Komenco overrides.
