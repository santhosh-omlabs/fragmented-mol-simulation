# Takeaways - what this study taught about fragmentation

**What this is:** the lessons, in plain language, with the evidence for each. For the numbers behind them read [RESULTS.md](RESULTS.md); for the claim and its limits, [../README.md](../README.md).

**What this is not:** a claim about quantum advantage or about molecules too large for classical methods. Everything here was learned on small molecules we could check exactly.

## 1. The bottom line

Fragmenting a molecule cuts the quantum cost (gates, shots, qubits per piece) and pays for it in accuracy. On ammonia and on a water dimer the accuracy paid was larger than
the quantum cost saved was worth, and a quantum circuit built from classical CCSD amplitudes never beat the free classical choice it was built from. That is not a failure of the
idea; it is what a careful test on molecules this easy should find. The value of the study is the method for testing the idea honestly, and the specific traps that
look like successes.

## 2. What we learned about fragmentation itself

| Question | Answer from this study | Evidence |
|---|---|---|
| Does cutting save qubits? | Yes, but each extra 2 qubits per fragment bought only about 6 kcal/mol on NH3 (38.8 at 14 qubits down to 5.9 at 26; the whole molecule is 28). | ladder, [RESULTS](RESULTS.md) (Lesson 3 section) |
| Is the cut error real or an artifact of the solver? | Real. Exact solvers on every fragment still leave 5.9 to 17.3 kcal/mol on NH3 before any truncation. Cutting through covalent bonds is expensive. | ladder |
| Is a weak cut cheap? | Cheaper per cut (about 4 against 17 kcal/mol), but the thing you want to measure, an interaction energy, is itself small (6.6 kcal/mol), so the cut error can be as large as the signal. In a minimal basis the cut cost only 0.13. | Lesson 7 |
| Does moving the cut to a better place help? | Not by putting it on a covalent bond: that was worse at every matched size. | Lesson 7, unequal split |
| What does fragmenting reliably buy? | Exponentially fewer shots. Two-qubit gate count G sets the chance of an error-free shot, (1 - e)^G, so cutting G by 3x cuts shots by orders of magnitude (cost model: 7.5 million against about 3,100). | schedule, Lesson 6 |
| Does a second QPU speed things up? | Only if the work divides evenly. With one big fragment and three small ones, four QPUs gave 1.09x; splitting one circuit's shots gave close to 4x. | schedule |
| What do fragments cost that is easy to miss? | The bath carries electrons, so a cluster is bigger than its fragment: a water (24 qubits alone) plus its bath is 32 qubits untruncated. Clusters, not fragments, set the qubit count. | Lesson 7 |

## 3. What we learned about the quantum step

- SQD asks the quantum computer which rows of a huge matrix matter, then solves the small matrix classically. The energy always comes from the classical solve.
- With the same number of electron patterns, the ideal circuit beat random choices (36.0 against 4.8 kcal/mol at 8 patterns) and a simple excitation ordering (10.5), but
  not a ranking by the CCSD wavefunction (5.7), which needs no quantum computer. This follows from where the circuit's angles come from: they are computed from CCSD.
- A circuit can only add value if it carries information the classical baseline lacks. Ours could not. Circuits not derived from CCSD were not tried.
- On the whole 28-qubit molecule, hardware SQD gave 3.4 kcal/mol, random bits through the same loop 4.5, and a plain classical ordering 0.8.

## 4. What we learned from the hardware

- A real chip returns mostly noise at these gate counts: 6.7% of shots had a valid electron count at 1,107 routed gates, about what random bits give (7.5%).
- Circuit design matters more than anything else we varied. Restricting the circuit to hardware-friendly interactions cut the routed gates to 161 and raised valid shots to 37%.
- Cleaner samples did not translate into better answers here, because the classical post-processing (configuration recovery) already pulls samples toward low-excitation patterns, which a simple classical list also produces.
- All of this is three jobs, one device, one day, 14 quantum-seconds. It is a stamp, not a hardware result.

## 5. Ten working rules, each learned from something that went wrong or nearly did

1. **Compare at equal budget.** We first fed 100 random bitstrings to SQD and got the exact answer. The problem had only 35 electron patterns, so 100 samples cover all of them. Fix the subspace size and compare only which patterns were chosen.
2. **Run the strongest free classical baseline before the quantum method.** The CCSD-ranked choice was added late; it changed the conclusion.
3. **Control with random.** A hardware number means nothing until you know what random bits give. For a 14-qubit problem that is 7.5% valid shots, and hardware at 1,107 gates was 6.7%.
4. **A number that is the same for every method is a warning.** The full SQD loop gave an error of 1e-9 for hardware, simulation and random bits alike, because all of them covered all 35 patterns.
5. **Define the target and the pass/fail rule before you run.** For water we fixed the interaction-energy definition, the reference and the thresholds in advance. One of my own pre-written rules turned out to be wrong; it was caught before any result existed and the amendment is recorded. A rule written afterwards could not have been trusted.
6. **Obey the stop rule.** The water target was missed, so no quantum step was run on water. Running it anyway would have produced numbers that look like results.
7. **Test the code with identities, not just with outputs.** Whole molecule as one fragment must equal FCI; zero empty orbitals must equal Hartree-Fock; two independent codes must agree. These caught real energy-partition bugs (errors of 300+ kcal/mol) early.
8. **Reproducibility is a bug class.** The same code gave different circuits in different runs because of a free rotation among exactly degenerate orbitals and a multithreaded sampler. It silently invalidated an early table until fixed.
9. **Cost models are hypotheses.** Label every assumed hardware number, and say which direction the proxy errs ("error-free run" is pessimistic).
10. **Write down what you did not do.** Two water jobs ran out of memory and are reported as not obtained, with the effect on the verdict stated.

## 6. Where the real value is

- A working, tested, end-to-end fragment-and-SQD pipeline with independent checks, that you understand line by line.
- The vocabulary and instincts of the field: strings against determinants, bath against fragment, cut error against solver error, valid-shot fraction, matched subspace size, routed against logical gates.
- A set of honest numbers on where the idea stands for small molecules, including the parts that do not flatter it.
- A demonstrated habit that matters more than any single result: predict, test, control, and report what failed.

## 7. What it would take to turn this into a real research result

1. A molecule where the classical baselines actually fail (strong correlation, or size beyond coupled cluster), because only there can a circuit add information.
2. Circuits that do not inherit their angles from CCSD (variational or learned ansatz), compared at equal subspace size against the CCSD-ranked baseline.
3. Fragment schemes that handle weak interactions properly: a bath that includes correlated information, or FMO-style pair corrections, tested on water clusters with more memory than 7 GB.
4. Hardware with error mitigation, repeats for error bars, more than one device, and hardware-compatible compilation from the start.
5. Reference energies by a method that scales (DMRG or selected CI), so the comparison does not stop at the size where exact diagonalisation does.

## 8. Cheat sheet

```text
determinants          C(orbitals, up electrons)^2        NH3: C(14,4)^2 = 1,002,001
qubits                2 x spatial orbitals               14 orbitals -> 28 qubits
SQD subspace          k up strings x k down strings      k strings -> k^2 determinants
clean-shot chance     (1 - e)^G                          e = gate error, G = two-qubit gates
shots for N clean     N / (1 - e)^G
chemical accuracy     1 kcal/mol                          1 hartree = 627.5 kcal/mol
interaction energy    E(cluster) - sum E(monomers)        small difference of large numbers
cut error             total error with exact solvers, no truncation
```
