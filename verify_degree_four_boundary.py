"""Exact finite-dimensional boundary-moment certificates for C4.

Standard library only; all intervals are outward-rounded
dyadic enclosures. Analytic complement bound: S_tail <= D/3 at cutoff128.
No floating eigensolver or mixed-spectrum approximation enters acceptance.
"""
from fractions import Fraction as Q
from math import isqrt
import json

SCALE = 2**48

def need(test, message):
    if not test:
        raise ArithmeticError(message)

def ceildiv(a,b):
    return -((-a)//b)

class I:
    __slots__=('lo','hi')
    def __init__(self,lo,hi):
        need(lo<=hi,'interval order')
        self.lo,self.hi=lo,hi
    def __add__(self,other):
        if not isinstance(other,I): other=rat(Q(other))
        return I(self.lo+other.lo,self.hi+other.hi)
    __radd__=__add__
    def __neg__(self): return I(-self.hi,-self.lo)
    def __sub__(self,other): return self+-other
    def __mul__(self,other):
        if not isinstance(other,I): other=rat(Q(other))
        vals=(self.lo*other.lo,self.lo*other.hi,self.hi*other.lo,self.hi*other.hi)
        return I(min(vals)//SCALE,ceildiv(max(vals),SCALE))
    __rmul__=__mul__
    def __truediv__(self,other):
        if not isinstance(other,I): other=rat(Q(other))
        need(other.lo>0,'positive denominator')
        vals=[(x*SCALE,y) for x in (self.lo,self.hi) for y in (other.lo,other.hi)]
        return I(min(x//y for x,y in vals),max(ceildiv(x,y) for x,y in vals))

def rat(q):
    q=Q(q)
    return I((q.numerator*SCALE)//q.denominator,
             ceildiv(q.numerator*SCALE,q.denominator))

ZERO=rat(0)

def atan_sum(q,n):
    return sum((Q((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(n)),Q(0))

PI_LO=16*atan_sum(5,12)-4*atan_sum(239,5)
PI_HI=16*atan_sum(5,13)-4*atan_sum(239,4)
need(PI_LO>Q(103993,33102) and PI_HI<Q(22,7),'Machin enclosure')
TAU=I((SCALE*(1/PI_HI).numerator)//(1/PI_HI).denominator,
      ceildiv(SCALE*(1/PI_LO).numerator,(1/PI_LO).denominator))

def J(a,b):
    """Integral over [-1,-1/2] union [1/2,1] of cosine products.

    The cosine coordinate is t=(x+1)/2. Opposite parity vanishes.
    Same parity gives a rational plus a rational times 1/pi.
    """
    if (a+b)%2: return ZERO
    const=Q(0); invpi=Q(0)
    for k in (a-b,a+b):
        if k==0: const+=Q(1,2)
        else:
            need(k%2==0,'even frequency difference')
            sine=(0,1,0,-1)[(k//2)%4]
            invpi+=Q(2*sine,k)
    return rat(const)+rat(invpi)*TAU

def boundary_pair(a,b,m,n,jcache):
    return (1+(-1)**(a+m))*jcache[b][n]+(1+(-1)**(b+n))*jcache[a][m]

# Endpoint P, lowest included shell m²+n². All relaxed modes between that
# shell and the endpoint are included, including entire degenerate shells.
# The special symmetric row keeps phi02+phi20 but discards its antisymmetric
# partner by one extra L² orthogonality condition.
ROWS=[
    ('ground',Q(2),0),('p4',Q(11,2),1),
    ('p20_symmetric',Q(41,2),4),('p25',Q(55,2),5),
    ('p32',Q(65,2),8),('p41',Q(85,2),10),
    ('p68',Q(70),17),('p73',Q(75),18),('p80',Q(85),20),
    ('p100',Q(102),25),('p104',Q(211,2),26),
    ('p160',Q(323,2),40),('p164',Q(337,2),40),
    ('p169',Q(343,2),41),('p185',Q(371,2),45),
    ('p212',Q(214),53),('p272',Q(551,2),68),
    ('p292',Q(294),73),('p296',Q(595,2),74),
    ('p340',Q(683,2),85),('p425',Q(428),106),
    ('p436',Q(438),109),('p468',Q(469),117),
    ('p601',Q(603),149),('p833',Q(834),208),
    ('p1108',Q(2221,2),277),('p1361',Q(1365),340),
    ('p1700',Q(3413,2),425),
]

def ldl_positive(matrix):
    """Interval symmetric elimination; every exact pivot is positive."""
    n=len(matrix); pivots=[]
    a=[[x for x in row] for row in matrix]
    for k in range(n):
        pivot=a[k][k]
        need(pivot.lo>0,'positive LDL pivot '+str(k))
        pivots.append(Q(pivot.lo,SCALE))
        for i in range(k+1,n):
            for j in range(i,n):
                val=a[i][j]-(a[i][k]*a[j][k])/pivot
                a[i][j]=a[j][i]=val
    return min(pivots)

def verify_row(name,end_p,start,cutoff=128):
    e=end_p/16
    low=[]
    for m in range(isqrt((4*e).numerator//(4*e).denominator)+1):
        for n in range(isqrt((4*e).numerator//(4*e).denominator)+1):
            need(Q(m*m+n*n,4)!=e,'endpoint avoids every NN pole')
            if Q(m*m+n*n,4)<e: low.append((m,n))
    block=[[(m,n,Q(1))] for m,n in low if m*m+n*n>=start]
    if name=='p20_symmetric':
        need(block==[[(0,2,Q(1))],[(1,2,Q(1))],[(2,0,Q(1))],[(2,1,Q(1))]],
             'symmetric row untransformed block')
        block=[[(0,2,Q(1)),(2,0,Q(1))],[(1,2,Q(1))],[(2,1,Q(1))]]
    k=len(block); n_low=len(low)-k
    max_index=max(max(a,b) for column in block for a,b,_ in column)
    jcache=[[J(a,b) for b in range(cutoff+1)] for a in range(max_index+1)]
    def row_value(j,m,n):
        return sum((rat(c)*boundary_pair(a,b,m,n,jcache) for a,b,c in block[j]),ZERO)
    d=[[sum((rat(c)*row_value(i,a,b) for a,b,c in column),ZERO)
        for column in block] for i in range(k)]
    w=[]
    for col in block:
        q0=Q(col[0][0]**2+col[0][1]**2,4)
        need(all(Q(a*a+b*b,4)==q0 for a,b,_ in col),'one eigenvalue per column')
        norm=sum((c*c*(2 if a==0 else 1)*(2 if b==0 else 1) for a,b,c in col),Q(0))
        w.append(rat((e-q0)*norm))
    gram=[[ZERO for _ in range(k)] for _ in range(k)]
    for m in range(cutoff+1):
        for n in range(cutoff+1):
            q=Q(m*m+n*n,4)
            if q<=e: continue
            denominator=(q-e)*(2 if m==0 else 1)*(2 if n==0 else 1)
            bv=[row_value(j,m,n) for j in range(k)]
            nonzero=[j for j,v in enumerate(bv) if v.lo or v.hi]
            den=rat(denominator)
            for i in nonzero:
                for j in nonzero:
                    if j>=i:
                        gram[i][j]=gram[i][j]+bv[i]*bv[j]/den
    for i in range(k):
        for j in range(i): gram[i][j]=gram[j][i]
    cert=[[sum((d[i][t]*d[j][t]/w[t] for t in range(k)),ZERO)
           -gram[i][j]-d[i][j]/3 for j in range(k)] for i in range(k)]
    pivot=ldl_positive(cert)
    need(pivot>Q(3,4),'common three-quarter lower pivot margin')
    data=dict(name=name,end_p=str(end_p),upper_count=n_low,
              block_size=k,smallest_pivot_lower=str(pivot))
    print(json.dumps(data,sort_keys=True),flush=True)
    return end_p,n_low

def verify_tail_constant():
    b=128; emax=Q(3413,32)
    beta=(b+1)**2-4*emax
    need(beta>127**2,'square root lower bound for boundary tail')
    eps=4*(4/(Q(b)*(1-4*emax/(b+1)**2))+2/beta+2*Q(22,7)/127)
    need(eps<Q(1,3),'uniform one-third boundary Gram tail')
    print(json.dumps(dict(section='tail',epsilon_upper=str(eps),passed=True)),flush=True)

def main():
    verify_tail_constant()
    caps=[verify_row(*row) for row in ROWS]
    verify_coverage(caps)
    print(json.dumps(dict(section='matrix_certificates',passed=True,
                          caps=[(str(p),n) for p,n in caps])),flush=True)

def verify_coverage(caps):
    """All actual AND reset events to1707, conditional on scalar odd33 proof.

    The odd33 bound is independently proved by verify_degree_four_odd33.py;
    it is added here only after the exact sector count interface is checked.
    The all-energy coarse range P>=1707 comes from the separate tail proof.
    """
    # At P=163/2 the p80 matrix removes the two even-even modes (2,4),(4,2),
    # and the odd33 certificate removes the odd-odd mode(3,3).
    special=Q(163,2)
    all_low=[(m,n) for m in range(6) for n in range(6)
             if 4*(m*m+n*n)<=special]
    need(len(all_low)==22,'all NN modes at the special cap')
    odd=[(m,n) for m,n in all_low if m%2==n%2==1]
    need(odd==[(1,1),(1,3),(3,1),(3,3)],'odd-sector inclusive cap')
    even=[(m,n) for m,n in all_low if m%2==n%2==0]
    need((2,4) in even and (4,2) in even,'two distinct even-sector modes')
    caps=caps+[(special,19)]
    stop=1707
    ch=[0]*(stop+1); rh=[0]*(stop+1)
    for m in range(isqrt(stop)//2+1):
        for n in range(isqrt(stop)//2+1):
            if 4*(m*m+n*n)<=stop: ch[4*(m*m+n*n)]+=1
    for m in range(isqrt(stop)+1):
        for n in range(1,isqrt(stop)//4+1):
            if m*m+16*n*n<=stop: rh[m*m+16*n*n]+=1
    events={Q(p) for p in range(stop+1) if ch[p] or rh[p]}
    events.update(p for p,n in caps); events.add(Q(stop))
    hc=rc=0; cumulative=[]
    for p in range(stop+1):
        hc+=ch[p];rc+=rh[p];cumulative.append((hc,rc))
    minimum=None
    for p in sorted(events):
        h,r=cumulative[p.numerator//p.denominator]
        # Exclude each cap at its own endpoint. The resulting upper jump
        # is checked even though the proved chamber bound is inclusive.
        c=min([h]+[n for end,n in caps if p<end])
        margin=3*103993*p-16*33102*(c+2*r)
        need(margin>=0,'degree-four corrected inclusive event '+str(p))
        if p>0 and (minimum is None or margin<minimum[0]):minimum=(margin,p,c,r)
    need(len(caps)==29 and len(events)==300, 'complete cap/reset event totals')
    need(minimum==(Q(291),Q(73),19,12), 'minimum corrected event margin')
    print(json.dumps(dict(section='finite_boundary_coverage',passed=True,
        events=len(events),reset_caps=len(caps),minimum_positive_energy_margin=str(minimum[0]),
        minimizing_p=str(minimum[1]),chamber_bound=minimum[2],corridor_count=minimum[3],
        depends_on='separate exact odd33 and tail certificates')),flush=True)

if __name__=='__main__': main()
