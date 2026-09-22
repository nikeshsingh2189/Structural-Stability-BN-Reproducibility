"""The five ALARM monitoring relations used in the paper (Supplement C, Table S: exact inputs)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from bn_radius import load_network

CASES = [
    ("INTUBATION", "HYPOVOLEMIA", ["ARTCO2", "BP", "CVP", "HRBP", "PRESS", "SAO2"]),
    ("MINVOLSET", "SAO2", ["CO", "ERRCAUTER", "VENTMACH", "VENTTUBE"]),
    ("CO", "PRESS", ["CATECHOL", "CVP", "HR"]),
    ("HRSAT", "VENTMACH", ["ERRCAUTER", "FIO2", "INTUBATION", "LVFAILURE", "SHUNT", "VENTALV", "VENTLUNG", "VENTTUBE"]),
    ("HYPOVOLEMIA", "ERRLOWOUTPUT", ["BP", "MINVOL"]),
]


def load_cases():
    names, idx, par = load_network("alarm")
    rel = [(idx[x], idx[y], [idx[c] for c in C]) for x, y, C in CASES]
    return names, idx, par, rel
