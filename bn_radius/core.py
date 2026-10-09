"""Core routines: d-separation, deletion radius (Theorem 1 / Corollary 1), Algorithm 1.

A DAG is a dict  parents[v] -> list of parent vertices  (vertices are hashable, usually int).
"""
from __future__ import annotations
import heapq
import itertools
from collections import deque
from typing import Dict, Hashable, Iterable, List, Optional, Sequence, Set, Tuple

Vertex = Hashable
Parents = Dict[Vertex, List[Vertex]]
INF = float("inf")


def ancestors(par: Parents, S: Iterable[Vertex]) -> Set[Vertex]:
    """Ancestral closure An_G(S) (includes S)."""
    seen = set(S)
    stack = list(seen)
    while stack:
        v = stack.pop()
        for p in par[v]:
            if p not in seen:
                seen.add(p)
                stack.append(p)
    return seen


def moral_graph(par: Parents, A: Set[Vertex]) -> Dict[Vertex, Set[Vertex]]:
    """moral(G[A]): skeleton edges plus edges between co-parents, undirected."""
    adj: Dict[Vertex, Set[Vertex]] = {v: set() for v in A}
    for v in A:
        ps = [p for p in par[v] if p in A]
        for p in ps:
            adj[v].add(p)
            adj[p].add(v)
        for a, b in itertools.combinations(ps, 2):
            adj[a].add(b)
            adj[b].add(a)
    return adj


def d_separated(par: Parents, X: Vertex, Y: Vertex, Z: Iterable[Vertex]) -> bool:
    """True iff X and Y are d-separated given Z (ancestral-moral criterion)."""
    Z = set(Z)
    adj = moral_graph(par, ancestors(par, {X, Y} | Z))
    seen = {X}
    queue = deque([X])
    while queue:
        v = queue.popleft()
        if v == Y:
            return False
        for w in adj[v]:
            if w not in seen and w not in Z:
                seen.add(w)
                queue.append(w)
    return True


# ---------------------------------------------------------------- reference (exhaustive) radius
def deletion_radius_bruteforce(par: Parents, X: Vertex, Y: Vertex, C: Sequence[Vertex],
                               weights: Optional[Dict[Vertex, float]] = None,
                               protected: Iterable[Vertex] = ()):
    """Definition 2 / 3 by enumeration of all deletions R subset of C (minus `protected`).

    Returns (radius, [all minimum-cost breaking deletions]).  radius = inf if none exists.
    """
    C = list(C)
    prot = set(protected)
    free = [c for c in C if c not in prot]
    w = (lambda c: 1.0) if weights is None else (lambda c: weights[c])
    best, best_sets = INF, []
    for r in range(len(free) + 1):
        for R in itertools.combinations(free, r):
            cost = sum(w(c) for c in R)
            if cost > best:
                continue
            if not d_separated(par, X, Y, set(C) - set(R)):
                if cost < best:
                    best, best_sets = cost, [frozenset(R)]
                elif cost == best:
                    best_sets.append(frozenset(R))
    return best, best_sets


# ---------------------------------------------------------------- fixed-graph shortest path
def radius_moral_graph(par: Parents, X: Vertex, Y: Vertex, C: Iterable[Vertex],
                       weights: Optional[Dict[Vertex, float]] = None,
                       protected: Iterable[Vertex] = ()) -> float:
    """Right-hand side of Theorem 1 / Corollary 1 on the explicit moral graph H_C.

    `protected` items get infinite deletion cost (used for the protection experiments)."""
    C = set(C)
    prot = set(protected)
    adj = moral_graph(par, ancestors(par, {X, Y} | C))

    def cost(v):
        if v not in C:
            return 0.0
        if v in prot:
            return INF
        return 1.0 if weights is None else weights[v]

    dist = {X: cost(X)}
    heap = [(dist[X], 0, X)]
    tick = 1
    while heap:
        d, _, v = heapq.heappop(heap)
        if d > dist.get(v, INF):
            continue
        for w in adj[v]:
            nd = d + cost(w)
            if nd < dist.get(w, INF):
                dist[w] = nd
                heapq.heappush(heap, (nd, tick, w))
                tick += 1
    return dist.get(Y, INF)


# ---------------------------------------------------------------- sparse auxiliary graph J_C
def sparse_auxiliary_graph(par: Parents, A: Set[Vertex]):
    """J_C: original vertices + skeleton edges; for each child v with >= 2 parents in A an
    auxiliary vertex ('aux', v) joined to every parent of v (no explicit moral cliques)."""
    adj: Dict[Vertex, Set[Vertex]] = {v: set() for v in A}
    for v in A:
        ps = [p for p in par[v] if p in A]
        for p in ps:
            adj[v].add(p)
            adj[p].add(v)
        if len(ps) >= 2:
            m = ("aux", v)
            adj[m] = set(ps)
            for p in ps:
                adj[p].add(m)
    return adj


def deletion_radius(par: Parents, X: Vertex, Y: Vertex, C: Iterable[Vertex],
                    weights: Optional[Dict[Vertex, float]] = None,
                    check_valid: bool = True) -> Tuple[float, Optional[frozenset]]:
    """Algorithm 1.  Returns (radius, R*) where R* is a minimum breaking deletion, or (inf, None).

    Shortest path on J_C with lexicographic label (conditioned-vertex cost, number of original
    vertices): among minimum-cost paths this picks one with the fewest original vertices,
    which is the tie-break required by Lemma 3.
    """
    C = set(C)
    if check_valid and not d_separated(par, X, Y, C):
        raise ValueError("input is not a valid d-separation X _||_ Y | C")
    A = ancestors(par, {X, Y} | C)
    J = sparse_auxiliary_graph(par, A)

    def step(v):  # cost of entering v
        if isinstance(v, tuple) and len(v) == 2 and v[0] == "aux":
            return (0.0, 0)
        a = (1.0 if weights is None else weights[v]) if v in C else 0.0
        return (a, 1)

    start = (step(X)[0], step(X)[1])
    label = {X: start}
    pred = {X: None}
    heap = [(start, 0, X)]
    tick = 1
    while heap:
        lab, _, v = heapq.heappop(heap)
        if lab > label.get(v, (INF, INF)):
            continue
        if v == Y:
            break
        for w in J[v]:
            s = step(w)
            nl = (lab[0] + s[0], lab[1] + s[1])
            if nl < label.get(w, (INF, INF)):
                label[w] = nl
                pred[w] = v
                heapq.heappush(heap, (nl, tick, w))
                tick += 1
    if Y not in label:
        return INF, None
    path, v = [], Y
    while v is not None:
        path.append(v)
        v = pred[v]
    R = frozenset(v for v in path if not (isinstance(v, tuple)) and v in C)
    return label[Y][0], R
