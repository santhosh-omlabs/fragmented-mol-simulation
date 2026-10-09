"""Lesson 7.4: how much does cutting a water dimer into two waters cost, as a function of cluster size?

For k = 0, 1, 2, ... empty orbitals kept per cluster, with exact FCI on every cluster:

  E_frag(dimer, k)         sum of the two water fragments' energy shares (+ frozen-core energy)
  primary   E_int = E_frag(dimer, k) - 2 E_FCI(water)                 pre-registered, tiers T1 and T2
  secondary E_int = E_frag(dimer, k) - 2 E_cluster(water alone, k)    exploratory, same procedure on both sides
  error     E_int - E_int(CCSD(T), raw) = -6.59 kcal/mol from data/water_reference.json

Run:  python -m fragsim.water_cut [--kmax 9] [--partition waters|hbond-shift]   (minutes; writes data/water_dimer_cut.json after every k)
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from pyscf import scf

from fragsim.embedding import build_atom_basis, density_in_basis
from fragsim.molecule import BASIS, build_mol, run_reference
from fragsim.solver import frozen_core_energy, solve_cluster
from fragsim.water import HARTREE_TO_KCAL, dimer_atoms, flatten, monomer_atoms, water_fragments

DATA = Path(__file__).resolve().parents[2] / "data"
RECORD = DATA / "water_dimer_cut.json"

# atoms of flatten(dimer_atoms()): donor water A = 0 (O), 1 (H donated to B's oxygen), 2 (free H); acceptor water B = 3, 4, 5
PARTITIONS: dict[str, tuple[list[list[int]], Path]] = {
    "waters": (water_fragments(2), RECORD),
    "hbond-shift": ([[3, 4, 5, 1], [0, 2]], DATA / "water_dimer_cut_unequal.json"),
}


def monomer_cluster_energy(k: int | None) -> tuple[float, int]:
    """One water as its own (whole-molecule) cluster with k empty orbitals kept. Returns (hartree, qubits)."""
    mol = build_mol(monomer_atoms(), BASIS)
    mf = scf.RHF(mol).run()
    basis = build_atom_basis(mol, mf)
    d = density_in_basis(basis, mf)
    part = solve_cluster(mol, mf, basis, d, [0, 1, 2], n_virtual=k)
    return frozen_core_energy(mol, mf) + part.e_frag, part.n_qubits


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kmax", type=int, default=9)
    parser.add_argument("--kmin", type=int, default=0)
    parser.add_argument("--partition", choices=list(PARTITIONS), default="waters")
    args = parser.parse_args()
    partition, record_path = PARTITIONS[args.partition]

    ref = json.loads((DATA / "water_reference.json").read_text(encoding="utf-8"))["systems"][
        "dimer"
    ]
    e_int_ref = ref["interaction_kcal_mol"]["ccsd_t"]["raw"]
    e_dimer_ref = ref["energies_hartree"]["cluster"]["ccsd_t"]
    e_fci_water = run_reference(build_mol(monomer_atoms(), BASIS), n_frozen=1).e_fci

    mol = build_mol(flatten(dimer_atoms()), BASIS)
    mf = scf.RHF(mol).run()
    basis = build_atom_basis(mol, mf, n_frozen=2)
    d = density_in_basis(basis, mf, n_frozen=2)
    e_core = frozen_core_energy(mol, mf, n_frozen=2)

    rows = (
        json.loads(record_path.read_text(encoding="utf-8"))["rows"]
        if record_path.exists() and args.kmin > 0
        else []
    )
    for k in range(args.kmin, args.kmax + 1):
        t0 = time.perf_counter()
        parts = [
            solve_cluster(mol, mf, basis, d, atoms, n_virtual=k, n_frozen=2) for atoms in partition
        ]
        e_dimer = e_core + sum(p.e_frag for p in parts)
        e_mono_k, mono_qubits = (
            monomer_cluster_energy(k) if (k <= 8 and args.partition == "waters") else (None, None)
        )
        row = {
            "k": k,
            "qubits_per_fragment": [p.n_qubits for p in parts],
            "e_dimer_fragments_hartree": e_dimer,
            "total_energy_error_vs_ccsd_t_kcal_mol": (e_dimer - e_dimer_ref) * HARTREE_TO_KCAL,
            "primary_e_int_kcal_mol": (e_dimer - 2 * e_fci_water) * HARTREE_TO_KCAL,
            "monomer_cluster_qubits": mono_qubits,
            "secondary_e_int_kcal_mol": None
            if e_mono_k is None
            else (e_dimer - 2 * e_mono_k) * HARTREE_TO_KCAL,
            "seconds": round(time.perf_counter() - t0, 1),
        }
        row["primary_error_kcal_mol"] = row["primary_e_int_kcal_mol"] - e_int_ref
        row["secondary_error_kcal_mol"] = (
            None
            if row["secondary_e_int_kcal_mol"] is None
            else row["secondary_e_int_kcal_mol"] - e_int_ref
        )
        rows.append(row)
        sec = row["secondary_error_kcal_mol"]
        print(
            f"k={k:2d} qubits {row['qubits_per_fragment']}  total err {row['total_energy_error_vs_ccsd_t_kcal_mol']:8.2f}  "
            f"E_int primary err {row['primary_error_kcal_mol']:8.2f}  secondary err "
            f"{'n/a' if sec is None else format(sec, '8.2f')}  ({row['seconds']} s)",
            flush=True,
        )
        record_path.write_text(
            json.dumps(
                {
                    "what": "water dimer, exact FCI per cluster, interaction-energy error vs CCSD(T) raw",
                    "partition_name": args.partition,
                    "partition": partition,
                    "e_int_reference_kcal_mol": e_int_ref,
                    "e_fci_water_hartree": e_fci_water,
                    "rows": rows,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    print(f"wrote {record_path}")


if __name__ == "__main__":
    main()
