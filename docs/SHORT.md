# In short

What if you could shrink a quantum chemistry problem just by cutting the molecule into pieces?

That was our question. Smaller pieces mean fewer qubits, shorter circuits, fewer errors. We took ammonia, split it into fragments, and built the whole pipeline ourselves: orbitals, embedding, exact solvers, a quantum sampling circuit, and a method called SQD. We ran it in simulation, on a classical emulator, and on a real IBM quantum computer, using fourteen seconds of quantum time.

Then reality pushed back. Our first result looked perfect, until we noticed that a hundred random guesses did equally well. The quantum circuit beat random choices but lost to a free classical shortcut. On the real chip, most samples came back as noise. Our water experiment missed its target, so we stopped on purpose.

But look what we found. Cutting circuits shrinks the shots needed from millions to thousands, in our model. Redesigning the circuit for the chip lifted valid results from seven percent to thirty-seven. Two independent codes agreed to a fraction of a kilocalorie.

Most of all, we learned how to test these ideas honestly: fair comparisons, strong baselines, random controls. That is exactly the skill real quantum chemistry needs next.
