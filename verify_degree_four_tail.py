#!/usr/bin/env python3
"""Exact degree-four separated-count tail and complete finite defect catalogue.

Run with Python 3.8+: python -B -O verify_degree_four_tail.py
Standard library only; no floating-point arithmetic and no removable asserts.
The analytic centered-block estimate is proved in the manuscript.
The output concerns the fully relaxed NN chamber, not the actual C4 spectrum.
"""
from array import array
from fractions import Fraction as Q
from math import isqrt
import json


def require(ok, detail):
    if not ok:
        raise ArithmeticError(detail)


def emit(section, **data):
    print(json.dumps(dict(section=section, **data), sort_keys=True), flush=True)


def atan_bounds(q, terms):
    value = sum((Q((-1)**j, (2*j+1)*q**(2*j+1)) for j in range(terms)), Q(0))
    other = value + Q((-1)**terms, (2*terms+1)*q**(2*terms+1))
    return min(value,other), max(value,other)


def integral_lower(a, size, maximum, grid=32768, scale=10**12):
    lip = size*(size-1)//2
    previous = 0
    total = 0
    for j in range(grid*maximum):
        t = j+1
        following = isqrt(scale*scale*t*t//(a*a*grid*grid+t*t))
        m = min(sum(((s-k)*j) % grid for k in range(size)) for s in range(size))
        total += max(0,m-lip)*max(0,following-previous-1)
        previous = following
    return Q(total,size*grid*scale), total


def verify_tail():
    h,htotal = integral_lower(1,13,3)
    r,rtotal = integral_lower(4,7,8)
    require(h > Q(353,1000), "chamber lower integral")
    require(r > Q(279,1000), "corridor lower integral")
    require(2*Q(353,1000)+2*Q(279,1000)-1 == Q(33,125), "linear gain")
    chamber_error = Q(13**2*3,8)+Q(13,2)
    corridor_error = Q(7**2*8,8)+Q(7,2)
    constant = 1+chamber_error+2*corridor_error
    require(constant == Q(1407,8), "centered-block constant")
    require(Q(2,27) < Q(9,121) and 2 < Q(100,49), "radical bounds")
    require(Q(3,11)*(Q(10,7)+8) == Q(18,7), "square-root coefficient")
    require(1024 == 32**2, "endpoint square root")
    endpoint = -Q(33,125)*1024+Q(18,7)*32+constant
    require(endpoint < 0, "degree-four tail endpoint")
    emit("tail", passed=True, chamber_integral_lower=str(h), corridor_integral_lower=str(r),
         chamber_integer_sum=htotal, corridor_integer_sum=rtotal, grid=32768, root_scale=10**12,
         chamber_block_error=str(chamber_error), corridor_block_error=str(corridor_error),
         linear_gain="33/125", total_constant=str(constant), radius_cutoff=1024,
         endpoint_upper=str(endpoint))


def histogram(limit, first_x, first_y, step_x, step_y):
    result = array("I",[0])*(limit+1)
    require(result.itemsize >= 4, "unsigned histogram capacity")
    root = isqrt(limit)
    require((root+1)**2 < 2**32, "total representation bound")
    for x in range(first_x,root+1,step_x):
        xx=x*x
        for y in range(first_y,isqrt(limit-xx)+1,step_y):
            result[xx+y*y] += 1
    return result


def verify_events():
    lo5,hi5=atan_bounds(5,12)
    lo239,hi239=atan_bounds(239,4)
    pi_lo=16*lo5-4*hi239
    pi_hi=16*hi5-4*lo239
    require(pi_lo > Q(103993,33102), "finite-scan pi lower bound")
    limit=16*1024**2
    chamber=histogram(limit,0,0,2,2)
    corridor=histogram(limit,0,4,1,4)
    ch=rc=event_count=0
    defects=[]
    expected=[0,4,20,25,32,41,68,73,80,100,104,160,164,169,185,
              212,272,292,296,340,425,436,468,596,601,833,1108,1360,1361,1700]
    minimum_nondefect=None
    for p in range(limit+1):
        ch+=chamber[p]
        rc+=corridor[p]
        if not (chamber[p] or corridor[p]):
            continue
        event_count+=1
        total=ch+2*rc
        # Compact budget = 3*pi*P/16; integer scan uses strict lower pi.
        margin=3*103993*p-16*33102*total
        if margin < 0:
            # Check these are failures for exact pi as well. No spectral
            # assertion about the physical mixed chamber is being made.
            require(Q(total) > Q(3*p,16)*pi_hi, "genuine separated-count defect")
            needed=(16*33102*total-3*103993*p+16*33102-1)//(16*33102)
            require(Q(total-needed+1) > Q(3*p,16)*pi_hi,
                    "exact minimal integer removal at the event")
            catchup=Q(16*33102*total,3*103993)
            require(catchup < 1707, "all separated defects paid by P=1707")
            defects.append(dict(p=p,chamber_count=ch,corridor_count=rc,
                                separated_count=total,required_integer_removal=needed,
                                rational_margin=margin,
                                budget_catchup_upper=str(catchup)))
        elif minimum_nondefect is None or margin < minimum_nondefect[0]:
            minimum_nondefect=(margin,p,total)
    require([item['p'] for item in defects] == expected, "complete relaxed defect catalogue")
    require(event_count == 1657038, "complete relaxed event total")
    require((ch, rc) == (3296231, 3292602), "finite-cutoff multiplicities")
    require(minimum_nondefect == (7323, 17, 10), "minimum nondeficient event")
    require(corridor[limit] > 0, "finite-tail join is included")
    emit("finite_events", passed=True, scope="relaxed NN chamber plus twice corridor",
         actual_union_events=event_count, cutoff_p=limit, defect_count=len(defects),
         proved_all_energy_from_p=1707,
         cutoff_chamber_count=ch, cutoff_corridor_count=rc,
         minimum_nondefect=minimum_nondefect, defects=defects)


if __name__=='__main__':
    verify_tail()
    verify_events()
    emit("all_checks",passed=True)
