"""Scheduling arithmetic on toy numbers (no chemistry, fast)."""

import math

import pytest

from fragsim.schedule import (
    CostModel,
    Job,
    clean_fraction,
    job_seconds,
    makespan,
    shot_parallel_seconds,
    shots_needed,
)

MODEL = CostModel(error_per_cx=0.01, rep_delay_s=1e-3, layer_s=0.0, overhead_s=2.0)


def test_clean_fraction_and_shots_follow_the_exponential_rule():
    job = Job("x", 4, 100, 50)
    assert clean_fraction(job, MODEL) == pytest.approx(0.99**100)
    assert shots_needed(job, MODEL, n_clean=100) == math.ceil(100 / 0.99**100)


def test_halving_the_gates_squares_the_clean_fraction():
    big, small = Job("b", 8, 200, 0), Job("s", 4, 100, 0)
    assert clean_fraction(big, MODEL) == pytest.approx(clean_fraction(small, MODEL) ** 2)


def test_job_seconds_is_overhead_plus_shots_times_shot_time():
    assert job_seconds(Job("x", 2, 1, 0), 1000, MODEL) == pytest.approx(2.0 + 1000 * 1e-3)


def test_makespan_is_capped_by_the_largest_job():
    assert makespan([10.0, 1.0, 1.0, 1.0], 4) == 10.0  # one big fragment cannot be split
    assert makespan([10.0, 1.0, 1.0, 1.0], 1) == 13.0
    assert makespan([4.0, 4.0, 4.0, 4.0], 2) == 8.0  # balanced jobs scale perfectly


def test_shot_parallel_pays_overhead_on_every_qpu():
    job = Job("x", 2, 1, 0)
    one = shot_parallel_seconds(job, 4000, 1, MODEL)
    four = shot_parallel_seconds(job, 4000, 4, MODEL)
    assert one == pytest.approx(2.0 + 4.0)
    assert four == pytest.approx(2.0 + 1.0)
    assert one / four < 4  # overhead keeps speed-up below the number of QPUs
