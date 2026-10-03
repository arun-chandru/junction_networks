#!/usr/bin/env python3
"""Exact scalar constants for the nonuniform corridor extension.

Standard library only; run python -B -O verify_nonuniform.py.
The proof's continuum inequalities are derived in the manuscript;
this script verifies their rational endpoints and radical comparisons.
"""
from fractions import Fraction as Q
import json


def require(ok,message):
    if not ok:
        raise ArithmeticError(message)


def atan_sum(q,n):
    return sum((Q((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(n)),Q(0))


def main():
    pi_lower=16*atan_sum(5,12)-4*atan_sum(239,5)
    pi_upper=16*atan_sum(5,13)-4*atan_sum(239,4)
    require(pi_lower>Q(25,8) and pi_upper<Q(22,7),'Machin pi enclosure')
    require(Q(2,27)<Q(1,4),'c < 1/2')
    require(2<Q(9,4),'sqrt(2) < 3/2')
    require(5<Q(81,16),'sqrt(5) < 9/4')

    one=Q(25,32)-Q(8,25)
    two=Q(125,64)-Q(32,25)
    weighted_two=Q(75,32)-Q(36,25)
    require(one==Q(369,800) and one>Q(1,3),'one-row strip gap')
    require(two==Q(1077,1600) and two>Q(2,3),'two-row strip gap')
    require(weighted_two==Q(723,800),'weighted two-row gap')
    require(Q(2,27)*Q(8,3)==(Q(8,3)/6)**2,'large-r strip endpoint')
    require(Q(16,27)<2,'c sqrt(2T) <= T/2 for T >= 2')

    first_margin=9*one-4
    midpoint_gap=Q(225,128)-Q(9,8)
    second_margin=9*midpoint_gap-5
    third_margin=9*weighted_two-8
    require(first_margin==Q(121,800) and first_margin>0,'first active interval')
    require(midpoint_gap==Q(81,128),'one-row value at T=3/2')
    require(second_margin==Q(89,128) and second_margin>0,'second active interval')
    require(third_margin==Q(107,800) and third_margin>0,'third active interval')
    require(Q(3)+2*Q(3,2)/3==4,'d1 first-band radical comparison')
    require(Q(15,4)+Q(3,2)/6==4,'d4 first-band radical comparison')
    require(Q(13,2)+2*Q(3,2)/9<7,'two-row maximum')
    require(9*Q(1,3)-1==2,'active high-energy comparison')

    require(Q(9)*Q(1,2)>=4,'width-one-half length-nine example area')
    require((1-Q(1,2))/2==Q(1,4),'example port containment')
    require(Q(17,64)*Q(22,7)==Q(187,224)<1,'near-four increment obstruction')

    print(json.dumps(dict(passed=True,strip_aspect_threshold=3,
                         active_length_threshold=9,minimum_corridor_area=4,
                         first_band_margin=str(first_margin),
                         second_band_margin=str(second_margin),
                         third_band_margin=str(third_margin),
                         one_row_lower=str(one),two_row_lower=str(two),
                         weighted_two_row_lower=str(weighted_two)),sort_keys=True))


if __name__ == '__main__':
    main()

