# Intuition - Fragment-based SQD on ammonia

**What this is:** pedagogy. Analogies and mental models for Lessons 0-7, in reading order.

**What this is not:** the claim. For numbers and limits, read [../README.md](../README.md).

## How to read the symbols

| Say it | Written | Means |
|---|---|---|
| n choose k | C(n, k) | number of ways to pick k items out of n |
| E-FCI | E_FCI | exact energy in this basis, found by classical brute force |
| E-HF | E_HF | energy if each electron only feels the average of the others |
| dets | determinants | one specific arrangement of electrons over orbitals |

## 1. What the problem is

An energy of a molecule is the lowest eigenvalue of a huge matrix (the Hamiltonian). Each row and column of that matrix is one
arrangement of the electrons over the orbitals (a determinant). "Solve the molecule" means "find the lowest eigenvalue of that matrix".

Analogy: a seating plan. N guests (electrons) pick seats (orbitals), and every complete seating is one row of the matrix.
Where the analogy breaks: electrons are identical and cannot be told apart, so the matrix counts seatings, not named guests.

## 2. Why the matrix explodes

For NH3 in 6-31G we keep 14 orbitals and 8 electrons: 4 spin-up and 4 spin-down. Up and down electrons are placed independently.

```text
number of dets = C(orbitals, up electrons) ^ 2
NH3, 14 orbitals, 4 up:   C(14, 4)^2 = 1001^2 = 1,002,001
```

Say out loud: pick 4 of 14 seats for the up electrons, pick 4 of 14 for the down electrons, and multiply.
Our laptop handles one million rows (12 seconds). The growth is what hurts:

```text
orbitals   electrons   dets
   7          4              441        (C(7,2)^2)
  14          8        1,002,001        (C(14,4)^2)
  28         16   ~ 9.7e12              (C(28,8)^2)
```

Doubling the orbitals (and electrons) multiplies the matrix by roughly ten million at this size. This is the whole reason people look at quantum
computers and at fragmenting.

## 3. Qubits and circuit depth

One spatial orbital needs 2 qubits (spin up, spin down). So 14 orbitals need 28 qubits, and 7 orbitals need 14 qubits.
Fewer qubits also means fewer gates, and gates are where hardware errors come from:

```text
chance a circuit runs without any error  ~  (1 - e) ^ G
e = error rate of one gate, G = number of gates
```

Say out loud: each gate is a small chance to fail, and the chances multiply. Cutting G by a factor of about four makes a clean run
far more likely. (Illustration only: for e = 0.5% and G = 600, (1 - e)^G is about 5%; for G = 150 it is about 47%.)

## 4. What fragmenting does

Cut the molecule into pieces, solve each piece with its own small matrix, and add the answers up.

Analogy: planning a wedding seating chart table by table instead of for the whole hall at once. Each table is cheap to plan.
Where the analogy breaks: electrons in different fragments still feel each other. Two guests at different tables may be feuding.
Ignoring that costs accuracy. That cost is called the fragmentation error.

To soften it, each fragment is solved together with a small "bath": a few extra orbitals that stand in for the rest of the
molecule (this is the idea behind embedding methods such as the embedded wavefunction method).

## 5. The trade, in one line

```text
total error  =  error of solving each fragment (noise, sampling)  +  error of cutting the molecule
```

Fragmenting trades the first term for the second. It pays only if the second term stays small. In the fusion-salt paper
(arXiv:2606.30402, from its abstract) the cutting term was the larger one, which is why this repository measures it directly.

## 6. Where the quantum computer fits (SQD)

Sample-based quantum diagonalization (SQD) does not read the answer off the qubits.

```text
1. run a circuit on the QPU, measure many times -> a list of bitstrings (electron arrangements)
2. keep the arrangements that appear (and repair obviously wrong ones)
3. build the small matrix using only those arrangements   (classical)
4. find its lowest eigenvalue                              (classical)
```

Say out loud: the quantum computer suggests which rows of the huge matrix matter. A classical computer then solves the small matrix.
For NH3 the whole matrix has about a million rows, so a classical computer can do it directly. SQD here is practice for
molecules where that is no longer true.

## 7. Lesson 1 - give every orbital an address

Out of the box, the orbitals a computer finds are spread over the whole molecule, like a rumour that everyone has partly heard.
To cut a molecule into fragments we first need orbitals that belong to a place: "this one lives on nitrogen", "this one sits on hydrogen 2".
The recipe is called localization (here: intrinsic bond orbitals), and the result is a set of orbitals each tagged with an atom.

```text
before:   orbital = a little bit everywhere        (cannot be assigned to a fragment)
after:    orbital 3 -> N     orbital 5 -> H1     ...   (each has an address)
```

Analogy: sorting a pile of letters by street. Where it breaks: some orbitals (the N-H bonds) honestly live between two atoms, so the
sorting is a convention, not a fact of nature. Different conventions give slightly different fragment energies.

Two housekeeping facts. The nitrogen 1s orbital is a tight inner shell that barely takes part in bonding, so we freeze it (that is why 15 orbitals
become 14). And NH3 has threefold symmetry, so some orbitals have exactly the same energy and any rotation among them is equally valid. That arbitrary
choice once changed our circuits from run to run; the code now fixes it (see the lab book entry on reproducibility).

## 8. Lesson 2 - a fragment is not alone: the bath

Take the nitrogen fragment. It is entangled with the three hydrogens. If we solved it by itself we would pretend nitrogen has no neighbours.
The fix, from a family of methods called embedding (we follow the DMET-style bath), is to add a few extra orbitals that stand in for the neighbours.

```text
fragment orbitals  --- coupled to ---  environment orbitals
                       (one matrix, D_EF, records how strongly)
SVD of that matrix -> the few combinations of environment orbitals that
                      actually touch the fragment = the BATH
cluster = fragment + bath          <- this is what we solve
```

Say out loud: out of all the environment, keep only the part the fragment can see. Everything else is replaced by a smooth average field.
A number worth knowing: for a closed-shell molecule each bath orbital's coupling is at most 0.5, and our bath keeps those above 1e-6.

Analogy: asking one table at the wedding to plan around "the people who will actually walk past us", not all 200 guests.
Where it breaks: the bath is built once from a mean-field picture and never updated. A self-consistent method would iterate; ours does not.

## 9. Lesson 3 - adding the pieces up without double counting

Each cluster (fragment + bath) is solved exactly, but the bath orbitals belong to the neighbours too. If we added all cluster energies we would count the
bath twice. The rule that worked: each fragment only claims the energy tied to the occupied orbitals that live on its own atoms.

```text
E_total = E_nuclear + E_core + sum over fragments ( E_HF share + E_correlation share )
```

We trust this because of three tests: if the fragment is the whole molecule we get exact FCI (to 1e-9 hartree); if all empty orbitals are kept and the fragments
partition the molecule, the pieces sum to FCI; and with no empty orbitals we get plain Hartree-Fock. An independent code (Vayesta) agrees to about 0.3 kcal/mol.

The result is a ladder, from small circuits and big errors to big circuits and small errors:

```text
fragments                  largest circuit    error vs exact (kcal/mol)
N | H | H | H, 3 empty        14 qubits              38.8
N | H | H | H, all empty      22 qubits              17.3
N+H1 | H | H                  24 qubits              11.7
N+H1+H2 | H                    26 qubits               5.9
whole molecule                 28 qubits               0
```

One kcal/mol is "chemical accuracy". Each rung buys roughly 6 kcal/mol for 2 qubits. For a molecule as small as NH3 that is a poor deal, and
the lab says so. It is a teaching exercise: the same ladder is what you would climb on a molecule too large for the exact solver.

## 10. Lesson 4 - what SQD actually does, step by step

A "string" is one spin's electron pattern, for example "orbitals 0, 1, 2, 3 are occupied" out of 7 orbitals. There are C(7,4) = 35 of them for our
14-qubit nitrogen cluster. A determinant is a pair (up string, down string), so the full matrix has 35 x 35 = 1,225 rows.

```text
SQD in four moves
 1. circuit    prepare a state close to the ground state (here: LUCJ, its angles come from a cheap classical CCSD calculation)
 2. measure    each shot gives a bitstring = one up string + one down string
 3. collect    keep the distinct strings you saw: say k up strings and k down strings
 4. solve      keep ALL k x k pairs, write the Hamiltonian only on those rows, find the lowest eigenvalue
```

Say out loud: the quantum computer is asked "which rows matter?", never "what is the energy?". The energy always comes from the classical solve in move 4.
Two consequences you can check by hand. First, k strings give k squared determinants, so the matrix grows fast with k. Second, the answer can only
get better as strings are added, and with all 35 strings it is exactly FCI, for any source of strings at all.

Noisy hardware returns bitstrings with the wrong number of electrons. "Configuration recovery" flips bits toward the orbital occupancies of the previous
solution until the count is right, then repeats the solve a few rounds.

## 11. Lesson 4.5 - the fair test (and a trap we fell into)

Our first comparison fed 100 random bitstrings to SQD and got the exact answer. That looked wonderful and meant nothing: with only 35 strings in
total, 100 random draws cover all of them. The lesson generalises.

```text
A comparison is only fair if every method gets the SAME subspace size k.
Then the only question is: which k strings did you pick?
```

At k = 8 strings (64 determinants), errors in kcal/mol on the nitrogen cluster:

```text
random strings + Hartree-Fock    36.0
simple classical ordering        10.5     (fewest excitations from Hartree-Fock first)
ideal circuit                     4.8
CCSD-ranked (classical)           5.7     (rank strings by weight in the classical CCSD state)
```

The circuit beats naive choices, but a free classical ranking built from the same CCSD calculation does about as well, and at most k it equals the best possible
choice. That is not surprising once you know the circuit's angles come from CCSD: it cannot contain information CCSD lacks. Finite shots hurt more: about 440 shots to see
8 strings, about 10,000 for 12, and 16 strings were never reached in 20,000 shots, because the energy that is left sits in strings with probability near 1e-4.

## 12. Komenco versus real hardware

The Komenco gateway (a third-party platform provided by Automatski, https://automatski.com/platform.html) takes a circuit and returns exact bitstring probabilities, with no shot noise and no hardware errors. We treat it as a classical emulator and
verified that: its answers match our own exact simulation of the same 14-qubit circuit to 3e-14. It is useful as a plumbing test (circuit out, samples in, SQD on top)
and as a check that bit ordering and gate translation are right. It says nothing about what a real chip would do.

## 13. Lesson 5 - real hardware: what noise does

Run the same circuit on IBM's `ibm_fez` and look at one simple health number: the fraction of shots with the right electron count.
Pure random bits would pass this check 7.5% of the time for our 14-qubit problem (35 over 128, squared).

```text
two-qubit gates (CZ, after routing)    shots with the right electron count
        161   (heavy-hex circuit)               37 %
        485                                       11 %
      1,107                                        7 %   <- random bits give 7.5 %
```

Three ideas sit behind that table:

1. Routing. A chip connects each qubit to only a few neighbours, so a circuit that wants distant qubits to talk needs extra swap gates. Our circuit grew from 684
   all-to-all gates to 1,107 on the chip.
2. Circuit design. Restricting which orbitals are allowed to interact (heavy-hex-friendly pairs) cut the routed count to 161 while the ideal quality changed little beyond 8 strings.
3. The (1 - e)^G rule from section 3 predicts the direction exactly: fewer gates, more clean shots.

The surprise is on the other side of the table. Clean shots rose from 7% to 37%, but the subspaces chosen from those shots were no better than the simple classical ordering.
On the whole 28-qubit molecule, at the same 314 strings per spin: hardware 3.4 kcal/mol, random bits through the same loop 4.5, a plain classical ordering 0.8.
Hardware samples carried real signal (4.0% valid against 0.4% for random bits), but not enough to matter after the classical ordering is taken into account.

## 14. Lesson 6 - what fragmenting really buys, and what a second QPU does

Fragmenting shortens circuits, and shortening a circuit helps exponentially. To see 100 error-free runs you need about 100 / (1 - e)^G shots.
In our cost model (e = 0.5%, an assumption, not a measurement):

```text
whole molecule      2,240 gates     about 7.5 million shots    52 minutes of QPU time
N fragment (14q)      684 gates     about 3,100 shots          13 seconds
```

That is the real, physical benefit, and it matches what the multi-QPU paper is after. Now suppose you have several QPUs. There are two ways to use them:

```text
fragment-parallel:  one fragment per QPU       finish time = the slowest fragment
shot-parallel:      split one circuit's shots  finish time ~ shots / number of QPUs, plus a fixed cost per job
```

Analogy: a kitchen with several cooks. Four dishes of equal size finish in a quarter of the time. One huge casserole and three salads do not: everyone waits for the casserole.
For N / H / H / H the nitrogen fragment is the casserole, so four QPUs give only 1.09 times the speed, while splitting that one circuit's shots gives close to 4 times.
Bonus: the three hydrogen fragments are identical by symmetry, so you only need to run one of them.

## 15. The honest verdict for ammonia

- The cut costs real accuracy: about 6 kcal/mol for each 2 qubits you give back.
- A circuit built from CCSD amplitudes does not know more than CCSD. A free classical ranking matches it, and the best classical choice beat today's hardware at equal subspace size.
- Today's noise matters: 7% of shots valid at 1,107 gates, 37% at 161.
- Fragmenting still pays where it should: far fewer shots, shorter circuits.
- NH3 is easy classically (FCI in about 12 seconds). None of this says anything about molecules too big for classical methods. It says how to run the experiment honestly.

## 16. Lesson 7 - a weak cut is not automatically a cheap one

Ammonia's pieces are glued by covalent bonds, the strongest kind. Water molecules in a cluster are held by hydrogen bonds, roughly ten times weaker. That is the case
fragment methods were invented for, so Lesson 7 asks: cut a water dimer into its two waters, and how well do the pieces reproduce the energy of the pair sticking together?

The quantity of interest is the interaction energy:

```text
E_int = E(dimer) - E(water A) - E(water B)         about -6.6 kcal/mol in our basis (negative = they attract)
```

Say out loud: how much lower the pair sits than the two molecules apart. Two facts make this a hard target.

1. It is a small difference of large numbers. Each water carries about 130 kcal/mol of correlation energy, so a 5% error on one water is as big as the whole interaction.
   Errors must cancel between the dimer and the monomers, and that only happens if both are treated by the same procedure.
2. In a small basis, each molecule "borrows" the other's basis functions to improve itself (basis-set superposition error). In our basis that is about 30% of the raw binding.
   So we fix beforehand which version we compare (raw, same basis) and only compare against a reference computed the same way.

Analogy: weighing a captain by weighing the ship with and without them aboard. If your scale is off by a tenth of a percent, the captain disappears into the error. Where it breaks:
errors on the two weighings are partly correlated, which is exactly the cancellation we hope for and have to check.

What we found, reading the curve (error in the interaction energy against the largest cluster, exact solvers):

```text
20 qubits   25.0 kcal/mol          28 qubits   6.3
24 qubits   15.2                   30 qubits   4.8     <- target was 1.0
```

Why it falls slowly: a cluster is the fragment plus its bath, and the bath carries electrons, so one water plus its bath is already 32 qubits. We cannot keep every empty orbital, so
we keep the most important ones (natural virtual orbitals) and the error shrinks as we keep more.

Two controls tell us where the error does and does not come from:

- In a minimal basis (STO-3G), where we can solve everything exactly with no truncation, the cut error is 0.13 kcal/mol. So the machinery is sound.
- Moving the cut onto a covalent O-H bond (a bigger fragment, a stronger cut) was worse at every size. A bigger fragment is not automatically better if you cut something stronger.

Our guess for why 6-31G is harder, which we did not test: each water has empty orbitals that the other water's electrons correlate with, and a bath built from the occupied orbitals
alone cannot represent that. A minimal basis has no such empty orbitals.

The lesson about doing science: we wrote the pass/fail rule before running (1 kcal/mol with exact solvers). It failed, so we stopped and did not run the quantum step on water,
because that would only have produced numbers that look like results. A negative result, reported with what was not run and why, is a result.

## 17. Check yourself

1. Why does 14 orbitals mean 28 qubits?
2. Why is the number of arrangements C(14, 4) squared, not C(14, 8)?
3. In the trade-off line of section 5, which term do hardware errors belong to?
4. What does the bath contain, and why is it not just "all the other atoms"?
5. A subspace of 20 strings per spin has how many determinants? What is the subspace when you use all 35 strings?
6. Why did 100 random bitstrings recover the exact energy of a 35-string problem, and why is that not impressive?
7. Why can the CCSD-ranked classical ordering never be much worse than the circuit built from CCSD angles?
8. With 161 gates at 0.5% error each, about what fraction of shots is error-free? (Check with section 3.)
9. Why can four QPUs fail to speed up the N / H / H / H run, and what could you change?
10. Why is the water interaction energy harder to get right than the total energy of one water?
11. The STO-3G cut error was 0.13 kcal/mol but the 6-31G error was 4.8. What does that tell you, and what does it not tell you?

Answers.
1. Two qubits per spatial orbital (spin up and spin down).
2. Up and down electrons are placed independently: 4 up among 14 orbitals and 4 down among 14, so 1001 x 1001. Choosing 8 of 14 would not respect spin.
3. The first term, the error of solving each fragment.
4. The few combinations of environment orbitals that couple to the fragment (from an SVD of the fragment-environment block of the density matrix). The rest of the environment is replaced by an average field.
5. 400 determinants; 1,225, which is exactly FCI for that cluster.
6. Once every one of the 35 strings is present, the solve is exact whatever produced the strings. Matched subspace size is what makes a comparison meaningful.
7. The ranking comes from the same CCSD state that sets the circuit's angles, so the circuit has no extra information; it can only blur it with finite shots and noise.
8. (0.995)^161 is about 0.45, so roughly 45% in the cost model. On the chip, with the device's own errors, we measured 37% valid electron counts.
9. The nitrogen job dominates the finish time and cannot be split. Split its shots across the QPUs (shot-parallel) or give it a smaller bath.
10. It is a small difference of large numbers (about -6.6 against about 130 kcal/mol of correlation per water), so any error that does not cancel between the dimer and the monomers is as large as the effect.
11. It shows the energy partition and the bath are sound, because the cut itself cost little when everything could be solved exactly. It does not show why 6-31G is harder; our explanation (empty orbitals on the neighbour carry correlation the bath misses) was not tested.
