"""Komenco client against a fake transport: no network."""

import pytest
from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp

from fragsim.backends.komenco import (
    MAX_QUBITS,
    KomencoClient,
    KomencoEstimator,
    circuit_request,
    expectation_value,
    fits_open_gateway,
)


def _bell() -> QuantumCircuit:
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cx(0, 1)
    qc.measure_all()
    return qc


def _fake(body):
    assert body["num_qubits"] == 2 and body["measurements"] == [0, 1]
    return {"measurements": {"00": 0.5, "11": 0.5}}


def test_cap_matches_vendor_trial_limit():
    assert fits_open_gateway(12) and fits_open_gateway(15) and not fits_open_gateway(16)


def test_sixteen_qubits_rejected_before_any_call():
    client = KomencoClient(transport=lambda b: pytest.fail("must not post"))
    qc = QuantumCircuit(MAX_QUBITS + 1)
    qc.measure_all()
    with pytest.raises(RuntimeError, match="15 qubits"):
        client.run_probabilities(qc)


def test_request_wire_format():
    body = circuit_request(_bell(), shots=10, top_k=4, api_key="open")
    assert {"num_qubits", "operations", "measurements", "shots", "topK", "apiKey"} <= body.keys()


def test_bell_expectations():
    client = KomencoClient(transport=_fake)
    state = QuantumCircuit(2)
    state.h(0)
    state.cx(0, 1)
    zz = expectation_value(client, state, SparsePauliOp("ZZ"))
    assert zz == pytest.approx(1.0)


def test_estimator_runs_through_qiskit_primitive_api():
    client = KomencoClient(transport=lambda b: {"measurements": {"00": 0.5, "11": 0.5}})
    state = QuantumCircuit(2)
    state.h(0)
    state.cx(0, 1)
    result = KomencoEstimator(client).run([(state, SparsePauliOp("ZZ"))]).result()
    assert float(result[0].data.evs) == pytest.approx(1.0)


def test_sample_counts_sum_to_shots_and_are_seeded():
    client = KomencoClient(transport=_fake)
    a = client.sample_counts(_bell(), 1000, seed=7)
    b = client.sample_counts(_bell(), 1000, seed=7)
    assert sum(a.values()) == 1000 and a == b and set(a) <= {"00", "11"}
