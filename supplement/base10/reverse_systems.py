"""Exact rational systems for the reverse-carry Type II estimates.

The low-frequency mass is bounded by the minimum of its second moment
and the affine Hellinger estimate. The depth parameter retains the
contribution from slopes congruent to one modulo a large base power.
"""
from fractions import Fraction as F
import forward_systems as h

CLASSICAL_P=F(235,154)
CLASSICAL_G=F(59,433)
h.NAMES += ['ell0','deep']
h.D=len(h.NAMES)

def build(args,branch=None,deep_branch=0,generic=True):
    common=dict(nu_min=F(args.nu_min),nu_max=F(args.nu_max),
                generic_margin=F(args.generic_margin) if generic else None,
                arc_margin=F(args.arc_margin),width=F(args.width))
    rows,labels,m=h.setup(**common,branch=branch)
    raw_ell=h.NAMES.index('ell');ell0i=h.NAMES.index('ell0')
    for row,label in zip(rows,labels):
        if label in ['low layer nonempty','low second moment','low p moment','low fourth moment']:
            row[ell0i]+=row[raw_ell]
            row[raw_ell]=F(0)
    def le(a,b,label):
        rows.append(h.sub(a,b));labels.append(label)
    ell0=h.ex(ell0=1);ell=h.ex(ell=1);deep=h.ex(deep=1)
    x=m['x'];s=m['s'];blo=m['blo'];theta=m['theta'];nu=m['nu']
    le(ell,ell0,'good overlap bounded by low L2')
    le(ell,h.sub(h.scale(blo,h.EPS),h.scale(h.sub(blo,deep),F(args.delta))),
       'good overlap affine Hellinger')
    # Retain the classical fractional moment together with the 3/2 moment.
    le(ell0,h.add(h.scale(blo,CLASSICAL_G),h.scale(h.ex(a0=1),CLASSICAL_P-2)),
       'low classical p moment')
    le(s,h.add(h.scale(theta,CLASSICAL_G),h.scale(x,CLASSICAL_P-1)),
       'high classical p moment')
    excess=h.add(ell0,s,h.scale(x,-1))
    if deep_branch==0:
        le(excess,h.ex(),'deep zero branch')
        desired=h.ex(F(args.deep_margin))
    else:
        le(h.ex(),excess,'deep positive branch')
        desired=h.add(excess,h.ex(F(args.deep_margin)))
    le(deep,desired,'chosen deep exponent upper')
    le(desired,deep,'chosen deep exponent lower')
    m.update({'ell0':ell0,'deep':deep})
    return rows,labels,m
