"""Exact rational systems for the forward-carry Type II estimates.

Every exponent is measured relative to q. The high-frequency scale is
Q=q**theta with theta=10*nu/3. The functions below reconstruct each
constraint and objective in the paper from fixed rational parameters.
"""
from fractions import Fraction as F

P=F(235,154)
G=F(59,433)
GAMMA=F(27,77)
EPS=F(23,500)
G4=F(3,5000)
KAPPA=F(23,80)
NAMES=['nu','u','a0','ell','s','j','a','b','c','d','e','f','s1','s2','total']
D=len(NAMES)

def ex(const=0,**kw):
    a=[F(0) for _ in range(D+1)]
    a[-1]=F(const)
    for k,v in kw.items():a[NAMES.index(k)]=F(v)
    return a

def add(*args):
    return [sum((v[i] for v in args),F(0)) for i in range(D+1)]

def scale(v,t):
    return [F(t)*a for a in v]

def sub(a,b):return add(a,scale(b,-1))


def setup(nu_min=F(257,1000), nu_max=F(29,100), generic_margin=F(1,1000000),
          arc_margin=F(1,1000000), width=F(1,1000000), branch=None):
    rows=[];labels=[]
    def le(a,b,label):
        rows.append(sub(a,b));labels.append(label)
    nu=ex(nu=1);theta=scale(nu,F(10,3));blo=sub(ex(1),theta)
    u=ex(u=1);a0=ex(a0=1);x=sub(u,a0);ell=ex(ell=1);s=ex(s=1)
    le(ex(nu_min),nu,'nu lower')
    le(nu,ex(nu_max),'nu upper')
    le(add(nu,ex(-F(7,100))),u,'retained full amplitude')
    le(u,ex(KAPPA),'full amplitude upper')
    le(ex(),a0,'low amplitude nonnegative')
    le(a0,u,'high cap exponent nonnegative')
    # Every nonempty amplitude layer contains at least one atom.
    le(scale(a0,-2),ell,'low layer nonempty')
    le(scale(x,-1),s,'high layer nonempty')
    le(ell,scale(blo,EPS),'low second moment')
    le(ell,add(scale(blo,G),scale(a0,P-2)),'low p moment')
    le(ell,add(scale(blo,G4),scale(a0,2)),'low fourth moment')
    le(s,scale(theta,GAMMA),'high L1')
    le(s,add(scale(theta,G),scale(x,P-1)),'high p')
    le(s,add(scale(theta,EPS),x),'high second moment')
    le(s,add(scale(theta,G4),scale(x,3)),'high fourth moment')
    le(s,sub(theta,x),'high count cap')
    if generic_margin is not None:
        le(ex(-generic_margin),add(ell,s,scale(x,-1)),
           'not already controlled by nonzero generic input')
    meta={'theta':theta,'blo':blo,'x':x,'ell':ell,'s':s,'nu':nu}
    if branch is None:
        return rows,labels,meta

    j=ex(j=1);lam=add(scale(nu,2),scale(j,-1))
    aa=ex(a=1);bb=ex(b=1);cc=ex(c=1);dd=ex(d=1);ee=ex(e=1);ff=ex(f=1)
    rho=add(aa,bb,cc,dd,ee)
    E=add(lam,scale(rho,-1),ex(width))
    z1=add(lam,bb,scale(cc,-1),ex(width))
    z2=add(lam,aa,cc,scale(bb,-1),scale(ee,-1),ex(width))
    zarc=add(lam,aa,bb,scale(cc,-1),ex(width))
    for name in ['j','a','b','c','d','e','f']:
        le(ex(),ex(**{name:1}),name+' nonnegative')
    le(j,scale(nu,2),'count trivial maximum')
    le(rho,lam,'common denominator height')
    le(add(bb,ee),cc,'gcd orientation')
    le(ff,cc,'smooth part of g2')
    le(ex(-arc_margin),add(ell,scale(zarc,G),scale(x,-(2-P))),
       'not already source-arc deleted')
    inner_main=add(scale(add(lam,scale(add(aa,bb,cc),-1),ex(width)),GAMMA),
                   scale(add(aa,scale(bb,2)),F(1,21)))
    inner_err=add(scale(add(lam,ex(width)),F(5,6)),scale(aa,F(1,6)),
                  scale(bb,F(7,6)),scale(cc,-F(5,6)),scale(add(dd,ee),F(2,3)),
                  scale(theta,-F(10,21)))
    outer_main=add(scale(add(lam,scale(add(aa,bb,cc,ee),-1),ff,ex(width)),GAMMA),
                   scale(add(aa,cc,scale(ff,-1)),F(2,21)))
    outer_err=add(scale(add(lam,ex(width)),F(5,6)),scale(aa,F(7,6)),
                  scale(bb,-F(5,6)),scale(cc,F(7,6)),scale(dd,F(2,3)),
                  scale(ee,-F(5,6)),scale(ff,-F(1,2)),scale(theta,-F(10,21)))
    row1_base=add(dd,ee,aa,scale(bb,2),E,scale(x,-(1-P/2)))
    row1=[sub(row1_base,scale(add(aa,bb),1-G)),
          sub(row1_base,scale(sub(theta,E),(1-G)/2))]
    row2_base=add(dd,scale(aa,2),scale(cc,2),scale(ff,-1),E,scale(x,-(1-P/2)))
    row2=[sub(row2_base,scale(add(aa,cc,scale(ff,-1)),1-G)),
          sub(row2_base,scale(sub(theta,E),(1-G)/2))]
    choices=[[inner_main,inner_err],[scale(z2,GAMMA),sub(z2,scale(theta,1-GAMMA))],
             [outer_main,outer_err],row1,row2]
    for k,(pair,bit) in enumerate(zip(choices,branch)):
        le(pair[1-bit],pair[bit],f'branch {k} dominance')
    s1=ex(s1=1);s2=ex(s2=1);total=ex(total=1)
    le(s1,scale(z1,GAMMA),'inner basic hybrid')
    le(s1,sub(z1,x),'inner cap')
    le(s1,add(scale(z1,G),scale(x,P-1)),'inner cap p hybrid')
    le(s1,choices[0][branch[0]],'inner coprime hybrid')
    le(s1,choices[3][branch[3]],'inner adapted row')
    le(s2,choices[1][branch[1]],'outer basic hybrid')
    le(s2,choices[2][branch[2]],'outer coprime hybrid')
    le(s2,choices[4][branch[4]],'outer adapted row')
    le(s2,sub(z2,x),'outer cap')
    overlap=ex() if branch[1]==0 else sub(z2,theta)
    le(s2,add(s,overlap),'all-outer-representations mass')
    outer_p=scale(z2,G) if branch[1]==0 else sub(z2,scale(theta,1-G))
    le(s2,add(scale(x,P-1),outer_p),'outer cap p hybrid')
    le(total,add(s1,s2),'rank2 product')
    le(total,scale(s,2),'full source mass squared')
    meta.update({'lam':lam,'E':E,'z1':z1,'z2':z2,'zarc':zarc,
                 'objective':add(ell,total,scale(lam,-1))})
    return rows,labels,meta


def encoded_system(rows,labels,obj=None):
    result={'names':NAMES,'constraints':[
        {'label':label,'coefficients':[str(v) for v in row[:-1]],'rhs':str(-row[-1])}
        for row,label in zip(rows,labels)]}
    if obj is not None:
        result['objective']={'coefficients':[str(v) for v in obj[:-1]],
                             'constant':str(obj[-1])}
    return result
