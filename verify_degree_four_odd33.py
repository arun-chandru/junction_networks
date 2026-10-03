"""Exact full-complement certificate: lambda_4(odd,odd)>163*pi^2/32.

Standard library only. Run with python -B -O. All checks remain active.
No floating arithmetic or approximate mixed eigenvalue enters this script.
Analytic interpretation and complete tail derivation are in the manuscript.
"""
from fractions import Fraction as F
from math import isqrt
import json


def require(test, message):
    if not test:
        raise ArithmeticError(message)


def floor(x):
    return x.numerator // x.denominator


def ceil(x):
    return -((-x.numerator) // x.denominator)


def atan_bounds(q, n):
    a = sum((F((-1)**j, (2*j+1)*q**(2*j+1)) for j in range(n)), F(0))
    b = a + F((-1)**n, (2*n+1)*q**(2*n+1))
    return min(a,b), max(a,b)


def sin_even(k):
    require(k % 2 == 0, 'only integer sine values occur')
    return {0:0, 2:1, 4:0, 6:-1}[k % 8]


def j_coeff(a, n):
    """Exact J_an = c + s/pi for odd nonnegative indices a,n."""
    require(a % 2 == n % 2 == 1, 'odd tangential modes')
    constant = F(0)
    coefficient = F(0)
    for k in (a-n, a+n):
        if k == 0:
            constant += F(1,2)
        else:
            coefficient += F(2*sin_even(k), k)
    return constant, coefficient


def main():
    lo5, hi5 = atan_bounds(5, 12)
    lo239, hi239 = atan_bounds(239, 4)
    pi_lo = 16*lo5-4*hi239
    pi_hi = 16*hi5-4*lo239
    require(F(103993,33102)<pi_lo<pi_hi<F(22,7), 'Machin pi enclosure')
    reciprocal_lo, reciprocal_hi = 1/pi_hi, 1/pi_lo
    scale = 10**14

    def enclosure(c, s):
        endpoints = (c+s*reciprocal_lo, c+s*reciprocal_hi)
        return floor(scale*min(endpoints)), ceil(scale*max(endpoints))

    coeff = {1:F(-35,100), 3:F(1), 5:F(-53,100)}
    cutoff = 256
    e = F(163,32)
    delta = e-F(9,2)
    require(delta == F(19,32), 'target spectral displacement')
    fn = {}
    for n in range(1, cutoff+1, 2):
        c = s = F(0)
        for j,a in coeff.items():
            cj,sj = j_coeff(j,n)
            c += a*cj
            s += a*sj
        fn[n] = enclosure(c,s)

    # D_f = integral over all four Dirichlet face portions of |f|^2.
    dc = ds = F(0)
    for i,ai in coeff.items():
        for j,aj in coeff.items():
            cj,sj = j_coeff(i,j)
            dc += 4*ai*aj*cj
            ds += 4*ai*aj*sj
    d_lo,d_hi = enclosure(dc,ds)
    require(d_lo>0, 'positive boundary norm')
    beta = (cutoff+1)**2-4*e
    sqrt_beta_lower = isqrt(floor(beta))
    require(sqrt_beta_lower>0 and sqrt_beta_lower**2<=beta, 'tail square root')
    eps = 4*(F(4,cutoff)/(1-4*e/(cutoff+1)**2)
             +2/beta+2*F(22,7)/sqrt_beta_lower)
    tail_upper = ceil(eps*d_hi)  # Units 1/scale.

    finite_upper = 0
    terms = 0
    omitted_low = []
    for m in range(1,cutoff+1,2):
        for n in range(1,cutoff+1,2):
            den = 8*(m*m+n*n)-163
            if den <= 0:
                omitted_low.append((m,n))
                continue
            lo = fn[m][0]+fn[n][0]
            hi = fn[m][1]+fn[n][1]
            square_upper = max(lo*lo,hi*hi)
            # L(psi_mn)=2(F_m+F_n); q_mn-e=den/32.
            finite_upper += (128*square_upper+den*scale-1)//(den*scale)
            terms += 1
    require(omitted_low == [(1,1),(1,3),(3,1),(3,3)], 'complete low odd block')
    require(terms == 16380, 'complete finite complement size')

    # L(phi_33)=4F_3, positive for this witness.
    b_lo = 4*fn[3][0]
    b_hi = 4*fn[3][1]
    require(b_lo>0, 'positive low-mode boundary moment')
    required_lower = F(b_lo*b_lo,scale*scale)/delta
    gram_upper = F(finite_upper+tail_upper,scale)
    margin = required_lower-gram_upper
    require(margin>F(1,50), 'strict scalar Gram certificate with uniform margin')
    print(json.dumps(dict(passed=True, cutoff=cutoff, scale=scale,
        coefficients={str(k):str(v) for k,v in coeff.items()},
        target_q=str(e), target_p='163/2', odd_sector_index=4,
        finite_high_modes=terms, boundary_norm_upper=str(F(d_hi,scale)),
        low_moment_interval=[str(F(b_lo,scale)),str(F(b_hi,scale))],
        finite_gram_upper=str(F(finite_upper,scale)),
        tail_epsilon_upper=str(eps), tail_gram_upper=str(F(tail_upper,scale)),
        gram_upper=str(gram_upper), lower_certificate_margin=str(margin),
        certified_margin='>1/50'),sort_keys=True,indent=2))


if __name__=='__main__':
    main()
