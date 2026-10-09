"""Komenco gateway client: circuits in, bitstring probabilities out.

Third-party code. Komenco is a platform provided by Automatski (https://automatski.com/platform.html),
a provider separate from this repository and its author. This module is adapted from the Python
client (`komenco.py`) supplied for that platform. Check the provider's terms before redistributing it.

Sends a measured circuit to http://<host>:<port>/api/komenco. `KomencoEstimator`
turns the returned probabilities into <psi|P|psi>, one POST per group of
qubit-wise commuting Paulis. `KomencoClient.sample_counts` draws bitstring
samples from the returned distribution, which is the input SQD needs.

The wire format is JSON with num_qubits, operations, measurements, shots,
topK, and apiKey. The gateway is plain HTTP on a shared trial instance: do not
send anything you would not post publicly. Override the target with the
KOMENCO_HOST, KOMENCO_PORT and KOMENCO_API_KEY environment variables.

The trial instance is capped at 15 qubits. A 16-qubit circuit raises here,
before any HTTP call. The vendor describes the instance as a trial simulator and
it returns exact probabilities with no shot noise, so this repo treats it as a
classical emulator.
"""

from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.request
from collections.abc import Callable, Iterable
from typing import Any

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.primitives import BaseEstimatorV2
from qiskit.primitives.containers import DataBin, EstimatorPubLike, PrimitiveResult, PubResult
from qiskit.primitives.containers.estimator_pub import EstimatorPub
from qiskit.primitives.primitive_job import PrimitiveJob
from qiskit.quantum_info import SparsePauliOp

DEFAULT_HOST = os.environ.get("KOMENCO_HOST", "168.220.234.162")
DEFAULT_PORT = int(os.environ.get("KOMENCO_PORT", "80"))
DEFAULT_API_KEY = os.environ.get("KOMENCO_API_KEY", "open")  # vendor's public trial key
MAX_QUBITS = 15

# Gates the Komenco samples accept. Local transpile only; not a remote call.
_BASIS = [
    "ccx",
    "ccz",
    "cp",
    "crz",
    "cs",
    "csdg",
    "cswap",
    "cu",
    "cx",
    "cy",
    "cz",
    "h",
    "id",
    "measure",
    "p",
    "rx",
    "ry",
    "rz",
    "s",
    "sdg",
    "swap",
    "sx",
    "sxdg",
    "t",
    "tdg",
    "u",
    "x",
    "y",
    "z",
]
_GATE_ALIAS = {
    "cnot": "cx",
    "cu1": "cp",
    "phase": "p",
    "cphase": "cp",
    "toffoli": "ccx",
    "fredkin": "cswap",
}

Transport = Callable[[dict[str, Any]], dict[str, Any]]


def fits_open_gateway(n_qubits: int) -> bool:
    """True when the open key's 15-qubit cap allows this register."""
    return int(n_qubits) <= MAX_QUBITS


def count_measurement_groups(operator: SparsePauliOp) -> int:
    """How many gateway circuits one energy evaluation needs."""
    return len(_measurement_groups(operator))


class KomencoClient:
    """POST one circuit. `transport` replaces HTTP in tests."""

    def __init__(
        self,
        host: str = DEFAULT_HOST,
        port: int = DEFAULT_PORT,
        api_key: str = DEFAULT_API_KEY,
        *,
        shots: int = 4096,
        timeout_s: float = 300.0,
        transport: Transport | None = None,
    ) -> None:
        self.host = host
        self.port = int(port)
        self.api_key = api_key
        self.shots = int(shots)
        self.timeout_s = float(timeout_s)
        self.transport = transport
        self.n_requests = 0

    def run_probabilities(
        self, circuit: QuantumCircuit, *, top_k: int | None = None
    ) -> dict[str, float]:
        """Return bitstring -> probability. Qubit 0 is the rightmost character."""
        n = int(circuit.num_qubits)
        if n > MAX_QUBITS:
            raise RuntimeError(
                f"Komenco open gateway allows {MAX_QUBITS} qubits; this circuit has {n}."
            )
        if top_k is None:
            top_k = 1 << n
        body = circuit_request(
            circuit,
            shots=self.shots,
            top_k=int(top_k),
            api_key=self.api_key,
        )
        payload = self._post(body)
        if payload.get("error"):
            raise RuntimeError(str(payload["error"]))
        raw = payload.get("measurements")
        if not isinstance(raw, dict) or not raw:
            raise RuntimeError("Komenco response has no measurements.")
        out: dict[str, float] = {}
        for key, value in raw.items():
            bits = str(key).replace(" ", "")
            bits = bits.removeprefix("0b")
            bits = bits.zfill(n)
            if len(bits) != n or any(ch not in "01" for ch in bits):
                raise RuntimeError(f"Komenco bitstring is not {n} bits: {key!r}")
            prob = float(value)
            if prob == 0.0:
                continue
            out[bits] = out.get(bits, 0.0) + prob
        if not out:
            raise RuntimeError("Komenco measurements are all zero.")
        return out

    def sample_counts(
        self, circuit: QuantumCircuit, shots: int, *, seed: int | None = None
    ) -> dict[str, int]:
        """Draw `shots` bitstrings from the gateway's probabilities (one POST)."""
        probs = self.run_probabilities(circuit)
        keys = list(probs)
        weights = np.array([probs[k] for k in keys], dtype=float)
        weights /= weights.sum()
        rng = np.random.default_rng(seed)
        draws = rng.multinomial(int(shots), weights)
        return {k: int(c) for k, c in zip(keys, draws) if c}

    def _post(self, body: dict[str, Any]) -> dict[str, Any]:
        self.n_requests += 1
        if self.transport is not None:
            return self.transport(body)
        url = f"http://{self.host}:{self.port}/api/komenco"
        data = json.dumps(body).encode("utf-8")
        request = urllib.request.Request(
            url,
            data=data,
            headers={
                "Content-Type": "application/json; charset=utf-8",
                "Accept": "application/json",
            },
            method="POST",
        )
        try:
            with urllib.request.urlopen(request, timeout=self.timeout_s) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            detail = exc.read().decode("utf-8", errors="replace")
            raise RuntimeError(f"Komenco HTTP {exc.code}: {detail}") from exc
        except urllib.error.URLError as exc:
            raise RuntimeError(f"Komenco gateway unreachable: {exc.reason}") from exc


def circuit_request(
    circuit: QuantumCircuit,
    *,
    shots: int,
    top_k: int,
    api_key: str,
) -> dict[str, Any]:
    """JSON body for /api/komenco. Local transpile only."""
    if circuit.num_parameters:
        raise ValueError("Bind circuit parameters before sending to Komenco.")
    measured = transpile(circuit, basis_gates=_BASIS, optimization_level=1)
    operations: list[dict[str, Any]] = []
    measurements: list[int] = []
    for inst in measured.data:
        op = inst.operation
        name = _GATE_ALIAS.get(op.name.lower(), op.name.lower())
        qubits = [measured.find_bit(q).index for q in inst.qubits]
        if name == "barrier":
            continue
        if name == "measure":
            measurements.extend(qubits)
            continue
        if name in ("reset", "delay", "initialize", "store"):
            raise RuntimeError(f"Komenco client refuses gate {name!r}.")
        operations.append(
            {
                "gate": name,
                "params": [float(p) for p in op.params],
                "qubits": qubits,
            }
        )
    if not measurements:
        raise RuntimeError("Komenco circuit has no measurements.")
    return {
        "num_qubits": int(measured.num_qubits),
        "operations": operations,
        "measurements": measurements,
        "shots": int(shots),
        "topK": int(top_k),
        "apiKey": api_key,
    }


def pauli_expectation(probabilities: dict[str, float], label: str) -> float:
    """<P> from Z-basis probabilities after P has been rotated onto Z.

    `label` uses Qiskit order: the rightmost character is qubit 0, matching
    the rightmost bit of each bitstring.
    """
    n = len(label)
    acc = 0.0
    for bits, prob in probabilities.items():
        if len(bits) != n:
            raise ValueError(f"Bitstring length {len(bits)} != Pauli length {n}.")
        sign = 1
        for q, char in enumerate(reversed(label)):
            if char == "I":
                continue
            if char not in "XYZ":
                raise ValueError(f"Bad Pauli character {char!r} in {label!r}.")
            if bits[-(q + 1)] == "1":
                sign = -sign
        acc += float(prob) * sign
    return acc


def measurement_circuit(state: QuantumCircuit, axes: list[str]) -> QuantumCircuit:
    """Append the basis change that turns each axis into a Z measurement."""
    qc = state.copy()
    qc.remove_final_measurements(inplace=True)
    for q, axis in enumerate(axes):
        if axis == "I" or axis == "Z":
            continue
        if axis == "X":
            qc.h(q)
        elif axis == "Y":
            qc.sdg(q)
            qc.h(q)
        else:
            raise ValueError(f"Unknown Pauli axis {axis!r}.")
    qc.measure_all()
    return qc


def expectation_value(
    client: KomencoClient,
    state: QuantumCircuit,
    operator: SparsePauliOp,
    *,
    progress: bool = False,
) -> float:
    """<state|operator|state> using one POST per qubit-wise commuting group."""
    if int(state.num_qubits) != int(operator.num_qubits):
        raise ValueError("Circuit and operator qubit counts differ.")
    groups = list(_measurement_groups(operator))
    total = 0.0
    for index, (group, axes) in enumerate(groups, start=1):
        if all(axis == "I" for axis in axes):
            total += float(np.real(np.sum(group.coeffs)))
            continue
        measured = measurement_circuit(state, axes)
        if progress:
            print(
                f"     Komenco circuit {index}/{len(groups)}  gates={measured.size()}",
                flush=True,
            )
        t0 = time.perf_counter()
        probs = client.run_probabilities(measured)
        if progress:
            print(
                f"     returned in {time.perf_counter() - t0:.1f}s",
                flush=True,
            )
        for pauli, coeff in zip(group.paulis, group.coeffs):
            term = SparsePauliOp([pauli], [coeff]).simplify()
            label = term.paulis[0].to_label()
            weight = complex(term.coeffs[0])
            total += float(np.real(weight * pauli_expectation(probs, label)))
    return total


class KomencoEstimator(BaseEstimatorV2):
    """Estimator V2. VQE calls run([(circuit, observable, parameters)])."""

    def __init__(self, client: KomencoClient, *, default_precision: float = 0.0) -> None:
        self.client = client
        self._default_precision = float(default_precision)

    @property
    def default_precision(self) -> float:
        return self._default_precision

    def run(
        self,
        pubs: Iterable[EstimatorPubLike],
        *,
        precision: float | None = None,
    ) -> PrimitiveJob[PrimitiveResult[PubResult]]:
        if precision is None:
            precision = self._default_precision
        coerced = [EstimatorPub.coerce(pub, precision) for pub in pubs]
        job = PrimitiveJob(self._run, coerced)
        job._submit()
        return job

    def _run(self, pubs: list[EstimatorPub]) -> PrimitiveResult[PubResult]:
        return PrimitiveResult([self._run_pub(pub) for pub in pubs], metadata={"version": 2})

    def _run_pub(self, pub: EstimatorPub) -> PubResult:
        bound = pub.parameter_values.bind_all(pub.circuit)
        circuits, observables = np.broadcast_arrays(bound, pub.observables)
        evs = np.zeros_like(circuits, dtype=np.float64)
        stds = np.zeros_like(circuits, dtype=np.float64)
        for index in np.ndindex(*circuits.shape):
            evs[index] = expectation_value(
                self.client,
                circuits[index],
                _as_sparse_pauli(observables[index]),
            )
        data = DataBin(evs=evs, stds=stds, shape=evs.shape)
        return PubResult(
            data,
            metadata={
                "target_precision": pub.precision,
                "backend": "komenco",
                "requests": self.client.n_requests,
            },
        )


class KomencoSession:
    """Client plus estimator, with the metadata a run record stores."""

    def __init__(self, client: KomencoClient, estimator: KomencoEstimator) -> None:
        self.client = client
        self.estimator = estimator
        self.sampler = None
        self.device = "komenco"
        self.use_gpu = False
        self.backend_info: dict[str, Any] = {
            "device": "komenco",
            "host": client.host,
            "port": client.port,
            "api_key": client.api_key,
            "shots": client.shots,
            "max_qubits": MAX_QUBITS,
            "primitive": "KomencoEstimator (grouped Pauli probabilities)",
        }
        self.profile: dict[str, float] = {}

    def report(self) -> None:
        print("----- runtime split -----")
        print(f"     device           {self.device}  (use_gpu={self.use_gpu})")
        print(f"     gateway requests {self.client.n_requests}")
        print("     note             probabilities from /api/komenco, not Aer")
        print("-----")


def make_komenco_session(
    *,
    host: str = DEFAULT_HOST,
    port: int = DEFAULT_PORT,
    api_key: str = DEFAULT_API_KEY,
    shots: int = 4096,
    transport: Transport | None = None,
) -> KomencoSession:
    client = KomencoClient(
        host,
        port,
        api_key,
        shots=shots,
        transport=transport,
    )
    print(
        f"[ok] backend      Komenco  http://{host}:{port}/api/komenco  "
        f"apiKey={api_key!r}  shots={shots}"
    )
    return KomencoSession(client, KomencoEstimator(client))


def _as_sparse_pauli(observable: Any) -> SparsePauliOp:
    if isinstance(observable, SparsePauliOp):
        return observable
    try:
        paulis, coeffs = zip(*observable.items())
    except Exception as exc:
        raise TypeError(f"Komenco estimator needs a Pauli sum, got {type(observable)!r}.") from exc
    return SparsePauliOp(list(paulis), list(coeffs))


def _measurement_groups(operator: SparsePauliOp) -> list[tuple[SparsePauliOp, list[str]]]:
    """Qubit-wise commuting groups, each with one measurement axis per qubit."""
    simplified = SparsePauliOp(operator.paulis, operator.coeffs).simplify()
    if len(simplified) == 0:
        return []
    groups = simplified.group_commuting(qubit_wise=True)
    out: list[tuple[SparsePauliOp, list[str]]] = []
    for group in groups:
        group = group.simplify()
        axes = ["I"] * int(group.num_qubits)
        for pauli in group.paulis:
            for q, char in enumerate(reversed(pauli.to_label())):
                if char == "I":
                    continue
                if axes[q] not in ("I", char):
                    raise RuntimeError(f"Qubit {q} has both {axes[q]} and {char} in one group.")
                axes[q] = char
        out.append((group, axes))
    return out
