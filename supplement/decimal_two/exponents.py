"""Rational endpoint exponents for all 45 omission pairs."""

from fractions import Fraction as Q
from itertools import combinations, product

from mask_mellin import check_cells
from affine import (
    add, const, delta, intersection, kap2_range, kap5_range, log_bounds, scale, sub, value,
)

a,f,a2,a5,f2,f5 = map(Q,('.428','.258','.028','.2241','.019','.137'))
m_min,m_max = Q('.4195'),Q('.4801')
h_low,h_high,sep,taper = map(Q,('.0918','.09693','.00003','.306'))
delta_lower = Q('.096910')
coarse,support,mask_floor,density_loss = map(Q,('.249','.0477','.117','.0474'))


def check_high():
    m, s, d = ((Q(0), Q(1), Q(0), Q(0)),
               (Q(0), Q(0), Q(1), Q(0)),
               (Q(0), Q(0), Q(0), Q(1)))
    n = sub(const(1),m)
    low = add(m,const(h_low+sep))
    high = add(m,const(h_high+sep))
    v = add(scale(Q(1,2),high),const('.00001'))
    w = sub(sub(low,v),const('.00015'))
    end = add(v,w)
    base_a = Q(3,4)*(1-delta)-a
    for is_high, (kr,kp,ap,fp) in product(
        (False,True),
        ((kap2_range,kap5_range,a2,f2),(kap5_range,kap2_range,a5,f5)),
    ):
        h,lam = (h_high,high) if is_high else (h_low,low)
        xi_lo = add(scale(kr[0],s),scale(kp[0],d))
        xi_hi = add(scale(kr[1],s),scale(kp[1],d))
        activity = add(m,const('.020'))
        mismatch = sub(s,d)
        base_twist = sub(s,m)
        principal = (sub(m,const(m_min)), sub(const(m_max),m),
                     sub(s,const('.29')), sub(lam,s), d, mismatch)
        narrow = (sub(const('.05'),mismatch),) if kr==kap2_range else ()
        wide = (sub(mismatch,const('.05')),)
        cases = [
            ('inactive ordinary', (*narrow,sub(activity,xi_lo)), False, True, base_twist, None,Q('-.00349')),
            ('active affine', (*narrow,sub(xi_hi,activity)), True, True, base_twist, None,Q('-.00068')),
        ]
        if kr==kap2_range:
            cases.extend((
                ('wide imbalance', wide, False,False,sub(base_twist,const(kp[0]*Q('.024'))),None,Q('-.00345')),
                ('wide short twist', wide,False,False,base_twist,Q('.075'),Q('-.00104')),
                ('wide averaged tail', wide,False,False,base_twist,None,Q('.0134')),
            ))
        for label,extra,active,use_filter,twist,cap,wanted in cases:
            floors=(*principal,*extra)
            inc=sub(xi_hi,n)
            planes=list(floors)+[inc,twist]
            if cap is not None:
                planes.append(sub(twist,const(cap)))
            if active:
                length = add(const(h),sub(xi_lo,m))
                planes.append(sub(scale(Q(1,2),d),scale(2,length)))
            if use_filter and is_high:
                planes += [sub(d,v),sub(d,end),sub(s,v),sub(s,end)]
                balanced=(const(0),sub(d,v),w)
                pure=(const(0),sub(s,v),sub(s,d),w,sub(end,d))
                planes += [
                    add(scale(a-f,b),scale((ap-fp)*kr[0],p),scale(-Q('.4')*taper,w))
                    for b,p in product(balanced,pure)
                ]
            exponent_max=None
            for equations in combinations(set(planes),3):
                point=intersection(equations)
                if point is None or any(value(x,point)<0 for x in floors):
                    continue
                mm,ss,dd=point
                tw=max(0,value(twist,point))
                if cap is not None:
                    tw=min(tw,cap)
                exponent=(delta*(1-ss)-base_a*dd-(1-2*ap)*kr[0]*(ss-dd)
                          +tw+max(0,value(inc,point))/2)
                if active:
                    assert h+value(xi_hi,point)-mm-dd/2<Q('-.010')
                    exponent -= (delta_lower*max(0,value(inc,point))+
                                 delta_lower/2*min(dd/2,2*value(length,point)))
                if is_high and use_filter:
                    vv,ww=value(v,point),value(w,point)
                    balance=max(0,min(vv+ww,dd)-vv)
                    pure=max(0,min(vv+ww,ss)-max(vv,dd))
                    exponent -= 2*max(0,(a-f)*balance+(ap-fp)*kr[0]*pure-Q('.4')*taper*ww)
                exponent_max=exponent if exponent_max is None else max(exponent_max,exponent)
            assert exponent_max is not None and exponent_max < wanted, (label,is_high,kr,exponent_max)
            print(label, 'high cut' if is_high else 'low cut', kr, float(exponent_max))


def check_sides():
    lo2, hi2 = log_bounds(2)
    lo54, hi54 = log_bounds(5, 4)
    low_kappa2 = lo2/(3*hi2+hi54)
    high_kappa2 = hi2/(3*lo2+lo54)
    assert kap2_range[0] < low_kappa2 < high_kappa2 < kap2_range[1]
    assert delta_lower < 1-3*high_kappa2 < 1-3*low_kappa2 < delta
    b = support/Q('.806')
    assert density_loss+sep+Q('.00015')<support
    assert b<Q('.060')
    assert max(mask_floor+b,support/coarse)<Q('.192')
    assert mask_floor+4*b<Q('.354')<m_min-Q('.0001')
    assert mask_floor*2*kap2_range[0]>Q('.070')
    assert a*b<Q('.026')
    assert b-m_min+Q('.0001')<Q('-.360')
    assert Q(3,4)*(1-delta)-a>coarse
    assert (1-2*a2)*kap2_range[0]>coarse
    assert (1-2*a5)*kap5_range[0]>coarse
    assert density_loss-coarse*Q('.234')+Q('.097')*sep<0
    assert density_loss-coarse*Q('.29')-Q('.403')*(m_min-Q('.0001')-Q('.29'))+Q('.151')*Q('.29')<0
    assert Q('.097')+Q('.655')*Q('.29')-m_min+Q('.0001')<0
    assert density_loss-coarse*Q('.192')<0
    assert delta<h_high
    assert Q('.6')*a2+Q('.4')*Q('.004')<f2
    assert Q('.6')*a5+Q('.4')*Q('.0053')<f5
    average_gain=((2*a2-1)*kap2_range[1]+kap5_range[0]-delta)*Q('.05')
    assert average_gain>Q('.0158')
    assert Q('.075')-kap5_range[1]*Q('.05')>Q('.0158')
    assert Q('.024')<Q('.05')
    assert (1-delta-a)/2>Q('.2375')
    for kap, other, ap in ((kap2_range,kap5_range,a2),(kap5_range,kap2_range,a5)):
        assert (1-delta-ap)*kap[0]/2>Q('.201')*(kap[1]+other[1]/2)
    assert (1-delta-a)/2>Q('.201')
    assert Q('.201')-Q('.799')*kap2_range[0]<0
    for m in (m_min,m_max):
        low,high=m+h_low+sep,m+h_high+sep
        v=high/2+Q('.00001')
        w=low-v-Q('.00015')
        assert Q('.258')<v and Q('.252')<w
        assert delta*(1-low)<density_loss
        assert (m-Q('.0001')-kap5_range[1]*high)/kap2_range[1]>Q('.193')
        assert (m-Q('.0001')-kap2_range[1]*high)/kap5_range[1]>Q('.37')
        assert max(kap2_range[1]*(high-v),Q('.192')-kap5_range[0]*v)<Q('.087')
        assert Q('.097')*Q('.192')+f*Q('.087')-Q('.6')*taper*w<Q('-.005')
        assert delta+(f-delta-Q('.6')*taper)*w-h_low<Q('-.0005')
        assert Q('.2375')*(1-high)>Q('.100')
        assert (1-a)*(1-high-Q('.230'))>Q('.100')
        for xi in (2*mask_floor,low,high):
            assert density_loss+a*Q('.192')-(1-a)*xi+max(0,xi-m)+max(0,xi-(1-m))<0
    gam,z,y=Q('.192'),Q('.4836'),Q('.234')
    D=Q('2.4')*y
    assert -y+gam*(D-Q('.201'))+Q('.230')*D<Q('-.009')
    assert -y+gam*(D-Q('.201'))+z*(D-Q('.2375'))<Q('-.008')
    assert -y+(Q('.0485')+D/2)*(gam+z)<Q('-.011')
    assert Q('1.2')*(gam+z)<1
    assert Q('1.5')*(gam+z)/Q('1.25')**2<1
    assert -Q('.5')+(Q('.5')+Q('.0485'))*(gam+z)<-Q('.10')
    assert (-Q('.5')+Q('.5')*(Q('.750')+Q('.097')*Q('.230'))
            + Q('.0485')*(Q('.192')+b+Q('.48867'))+Q('.194')*b)<-Q('.06')


def verify():
    check_sides()
    check_high()
    vertices=check_cells(
        a=a,budget=support,min_short=m_min,order_slope=coarse,
        upper_decay=Q('.2375'),min_orbit=mask_floor,linear_saving=Q('.201'),
        cut_at_eta_one=h_low+sep,low_break=Q('.180'),high_break=Q('.233'),
        height_ceiling=Q('.750'),
    )
    assert vertices==110
    print('All pairs: Mellin certificate passed at',vertices,'vertices.')
