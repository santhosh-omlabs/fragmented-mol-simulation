# Literature

## The two papers that motivated the study

| Source | Taken from it | Not taken |
|---|---|---|
| Herrmann et al., "Demonstration of Parallel Multi-QPU Execution for Fragment-Based Quantum Chemistry Using On-Premises Hardware", arXiv:2610.07702 | Idea: fragment molecular orbital (FMO) fragments; shot-parallel vs fragment-parallel execution; efficiency as a metric | Their Quoll hardware, He clusters, QWFS solver, their numbers |
| Das et al., "Quantum Computations on Fusion Blanket Molten Salts", arXiv:2606.30402 | Idea: embedded-wavefunction (EWF) fragments; SQD per fragment; fragmentation error vs solver error | Their FLiBe data, IBM devices, error values |

Both summaries come from the arXiv abstract pages; the papers' full texts were not used to set any parameter in this repository. Their reported errors (fragment solves within
0.7 kcal/mol of FCI; fragmentation error about 12 to 110 kcal/mol) are quoted from the abstracts as motivation, not as results of this work.

## Method sources

| Source | Used for |
|---|---|
| Robledo-Moreno et al., "Chemistry beyond the scale of exact diagonalization on a quantum-centric supercomputer", Science Advances 11, eadu9991 (2025), arXiv:2405.05068 | Sample-based quantum diagonalization (SQD) and configuration recovery |
| Motta, Sung, Whaley, Head-Gordon, Shee, "Bridging physical intuition and hardware efficiency for correlated electronic states: the local unitary cluster Jastrow ansatz for electronic structure", Chemical Science 14, 11213 (2023), doi:10.1039/d3sc02516k | The LUCJ circuit, including its hardware-friendly restricted form |
| Knizia and Chan, "Density matrix embedding: a simple alternative to dynamical mean-field theory", Physical Review Letters 109, 186404 (2012) | The idea of a bath from the fragment-environment block of the density matrix |
| Knizia, "Intrinsic atomic orbitals: an unbiased bridge between quantum theory and chemical concepts", Journal of Chemical Theory and Computation 9, 4834 (2013) | Intrinsic atomic and bond orbitals (localization and the atom-tagged basis) |

The last two references are cited from memory of well-known papers and were not re-checked against the journals; verify page numbers before formal citation. The first two were checked against their arXiv and journal records.

## Software and services

| Source | Used for |
|---|---|
| Vayesta (Booth group), https://github.com/BoothGroup/Vayesta, commit 7f1639d | Independent cross-check of the embedding energies (`scripts/crosscheck_vayesta.py`, `scripts/crosscheck_vayesta_water.py`). Not a dependency of the package. |
| PySCF | Integrals, Hartree-Fock, FCI, CCSD(T), intrinsic atomic and bond orbitals |
| Qiskit, Qiskit Aer, ffsim, qiskit-addon-sqd | Circuits and simulation; LUCJ circuits; subspace diagonalization and configuration recovery |
| IBM Quantum, Qiskit Runtime, backend `ibm_fez` | Three hardware jobs, see [IBM_HARDWARE_JOBS.md](IBM_HARDWARE_JOBS.md) |
| Komenco, provided by Automatski, https://automatski.com/platform.html | A trial gateway that returns exact bitstring probabilities; used as a classical emulator. The client in `src/fragsim/backends/komenco.py` is adapted from the Python client supplied for it. A third-party provider, separate from this repository. |
