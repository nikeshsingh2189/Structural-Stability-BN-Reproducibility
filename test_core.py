"""Sanity tests: the fixed-graph characterization, the sparse algorithm and the returned deletion."""
import itertools, random, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from bn_radius import (d_separated, deletion_radius, deletion_radius_bruteforce, radius_moral_graph)


def random_dag(rng, n, p):
    return {v: [u for u in range(v) if rng.random() < p] for v in range(n)}


def test_proposition1_example():
    # X -> M <- Z -> Y  (labels 0..3 = X, Z, M, Y); C = {M, Z}
    X, Z, M, Y = 0, 1, 2, 3
    par = {X: [], Z: [], M: [X, Z], Y: [Z]}
    assert d_separated(par, X, Y, {M, Z})
    assert not d_separated(par, X, Y, {M})          # R1 = {Z} breaks
    assert d_separated(par, X, Y, set())            # R2 = {M, Z} restores separation
    rho, R = deletion_radius(par, X, Y, {M, Z})
    assert rho == 1 and R == frozenset({Z})


def test_random_unit_and_weighted():
    rng = random.Random(1)
    checked = 0
    for _ in range(1500):
        n = rng.randint(5, 10)
        par = random_dag(rng, n, rng.choice([0.2, 0.35, 0.5]))
        X, Y = rng.sample(range(n), 2)
        rest = [v for v in range(n) if v not in (X, Y)]
        for _ in range(3):
            C = rng.sample(rest, rng.randint(1, min(6, len(rest))))
            if not d_separated(par, X, Y, C):
                continue
            checked += 1
            w = {c: rng.choice([0, 0, 1, 2, 5]) for c in C}
            for weights in (None, w):
                brute, _ = deletion_radius_bruteforce(par, X, Y, C, weights)
                assert radius_moral_graph(par, X, Y, C, weights) == brute          # Theorem 1 / Cor. 1
                rho, R = deletion_radius(par, X, Y, C, weights)                    # Algorithm 1
                assert rho == brute
                if rho != float("inf"):
                    assert not d_separated(par, X, Y, set(C) - set(R))            # R* really breaks
                    cost = len(R) if weights is None else sum(weights[c] for c in R)
                    assert cost == brute
    assert checked > 500


def test_invalid_input_rejected():
    par = {0: [], 1: [0]}
    try:
        deletion_radius(par, 0, 1, set())
    except ValueError:
        return
    raise AssertionError("expected ValueError for a d-connected input")
