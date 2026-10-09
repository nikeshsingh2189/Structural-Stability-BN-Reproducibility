#!/usr/bin/env python3
"""Independent deterministic verification of the reported CSV outputs.

This script ignores stored agreement flags and recomputes every radius/cost from
`bn_radius.core`. It exits nonzero if any stored result is inconsistent.
"""
from __future__ import annotations
import math
from pathlib import Path
import networkx as nx
import pandas as pd
from bn_radius.core import d_separated, deletion_radius, deletion_radius_bruteforce, radius_moral_graph

ROOT=Path(__file__).resolve().parent

def Cparse(s): return tuple(int(x) for x in str(s).split(';')) if str(s).strip() else tuple()
def Wparse(s):
    d={}
    for z in str(s).split(';'):
        k,v=z.split(':'); d[int(k)]=float(v)
    return d

def eq(a,b,tol=1e-6):
    return (math.isinf(float(a)) and math.isinf(float(b))) or abs(float(a)-float(b))<=tol

ed=pd.read_csv(ROOT/'data'/'smartds_sfo_p3r_edge_pairs.csv')
G=nx.DiGraph(); G.add_edges_from(ed[['source','target']].itertuples(index=False,name=None))
assert nx.is_directed_acyclic_graph(G)
par={v:list(G.predecessors(v)) for v in sorted(G.nodes())}

fail=[]

r=pd.read_csv(ROOT/'results'/'random_valid_relation_validation_1000.csv')
for i,row in r.iterrows():
    C=Cparse(row.C); x=int(row.X); y=int(row.Y)
    a,_=deletion_radius(par,x,y,C); m=radius_moral_graph(par,x,y,C); b,_=deletion_radius_bruteforce(par,x,y,C)
    if not d_separated(par,x,y,C) or not (eq(a,row.radius) and eq(m,row.radius_moral) and eq(b,row.radius_bruteforce) and eq(a,m) and eq(a,b)):
        fail.append(('random',i))

s=pd.read_csv(ROOT/'results'/'stratified_radius_validation_1000.csv')
for i,row in s.iterrows():
    C=Cparse(row.C); x=int(row.X); y=int(row.Y)
    a,_=deletion_radius(par,x,y,C); m=radius_moral_graph(par,x,y,C); b,_=deletion_radius_bruteforce(par,x,y,C)
    if not d_separated(par,x,y,C) or not (eq(a,row.radius_algorithm) and eq(m,row.radius_explicit_moral) and eq(b,row.radius_bruteforce) and eq(a,row.target_radius) and eq(a,m) and eq(a,b)):
        fail.append(('stratified',i))

w=pd.read_csv(ROOT/'results'/'weighted_validation_500.csv')
for i,row in w.iterrows():
    C=Cparse(row.C); weights=Wparse(row.weights); x=int(row.X); y=int(row.Y)
    a,_=deletion_radius(par,x,y,C,weights=weights); m=radius_moral_graph(par,x,y,C,weights=weights); b,_=deletion_radius_bruteforce(par,x,y,C,weights=weights)
    if not d_separated(par,x,y,C) or not (eq(a,row.cost_algorithm) and eq(m,row.cost_moral) and eq(b,row.cost_bruteforce) and eq(a,m) and eq(a,b)):
        fail.append(('weighted',i))

if fail:
    print('FAILED:', fail[:20]); raise SystemExit(1)
print('Verified 1000 random + 1000 stratified + 500 weighted cases: 2500/2500 exact agreement.')
