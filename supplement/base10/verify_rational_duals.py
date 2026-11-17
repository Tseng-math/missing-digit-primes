"""Exact verification of rational dual and Farkas certificates.

The stored nonnegative coefficients are multiplied into the constraint
matrix, and every coefficient identity and strict sign is checked exactly.
The fixed mathematical systems are reconstructed by verify_type_ii.py.
"""
from fractions import Fraction as F

def verify_one(item):
    system=item['system']
    rows=system['constraints']
    D=len(system['names'])
    coefficients=[F(0) for _ in range(D)]
    bound=F(0)
    entries=item['dual'] if item['feasible'] else item['farkas']
    for entry in entries:
        y=F(entry['coefficient'])
        assert y>=0
        row=rows[entry['row']]
        assert entry['label']==row['label']
        coefficients=[x+y*F(a) for x,a in zip(coefficients,row['coefficients'])]
        bound+=y*F(row['rhs'])
    if item['feasible']:
        obj=system['objective']
        assert coefficients==[F(a) for a in obj['coefficients']]
        bound+=F(obj['constant'])
        assert bound==F(item['exact_bound'])
        assert bound<0
    else:
        assert coefficients==[F(0)]*D
        assert bound==F(item['exact_contradiction'])
        assert bound<0
    return bound
