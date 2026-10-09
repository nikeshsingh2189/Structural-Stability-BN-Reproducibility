import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from bn_radius.pipeline_fit import PAPER_EDGES, common_relation_partition, enumerate_relations


def test_identical_graphs_give_no_partial_patterns():
    rows, s = common_relation_partition([PAPER_EDGES] * 3)
    assert (s["common"], s["all"], s["some"], s["none"]) == (556, 190, 0, 366)
    assert s["rates"] == [100 * 190 / 556] * 3


def test_partition_is_consistent_with_single_graph_counts():
    sparse = {("SO", "UR"), ("SO", "SPO"), ("SPO", "L"), ("L", "W"), ("UR", "UE"), ("UE", "US"), ("W", "SIPO"), ("L", "MD")}
    dense = PAPER_EDGES | {("SO", "W"), ("UR", "MD")}
    rows, s = common_relation_partition([sparse, PAPER_EDGES, dense])
    assert s["all"] + s["some"] + s["none"] == s["common"] == len(rows)
    # the per-graph break counts over the common set must equal the column sums of the pattern table
    for j in range(3):
        assert round(s["rates"][j] * s["common"] / 100) == sum(r[3][j] for r in rows)
    # cross-check one graph against a direct enumeration
    rel = enumerate_relations(PAPER_EDGES)
    common = {(r[0], r[1], frozenset(r[2])) for r in rows}
    assert all(k in rel and rel[k][1] is not None for k in common)
