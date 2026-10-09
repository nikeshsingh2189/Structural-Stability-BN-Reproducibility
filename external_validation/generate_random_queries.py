#!/usr/bin/env python3
"""Regenerate the 1,000 random valid relations used in the reported experiment."""
from __future__ import annotations
import math, random
from pathlib import Path
import networkx as nx
import pandas as pd
from bn_radius.core import d_separated, deletion_radius

ROOT=Path(__file__).resolve().parent
SEED=2026092201
N=1000
ed=pd.read_csv(ROOT/'data'/'smartds_sfo_p3r_edge_pairs.csv')
G=nx.DiGraph(); G.add_edges_from(ed[['source','target']].itertuples(index=False,name=None))
par={v:list(G.predecessors(v)) for v in sorted(G.nodes())}; nodes=sorted(G.nodes())
rng=random.Random(SEED); rows=[]; tries=0
while len(rows)<N:
    tries+=1
    x,y=rng.sample(nodes,2); k=rng.randint(1,6)
    C=tuple(sorted(rng.sample([v for v in nodes if v not in (x,y)],k)))
    if not d_separated(par,x,y,C): continue
    radius,_=deletion_radius(par,x,y,C)
    if math.isinf(radius): continue
    rows.append((x,y,';'.join(map(str,C)),len(C)))
out=pd.DataFrame(rows,columns=['X','Y','C','C_size'])
out.to_csv(ROOT/'data'/'random_query_set_1000_regenerated.csv',index=False)
ref=pd.read_csv(ROOT/'data'/'random_query_set_1000.csv')
print(f'candidate draws: {tries}')
print(f'exact query-set match: {out.equals(ref)}')
if not out.equals(ref): raise SystemExit(1)
