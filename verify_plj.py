#!/usr/bin/env python3
"""Exact arithmetic certificates for Pólya inequalities for junction networks.

Run: python -B -O verify_plj.py

Python 3.8 or later; standard library only; no input files or network access.
The script verifies all finite counts, fractional-block integral lower bounds,
tail constants, and rational boundary-lift comparisons used in the paper. It
does not approximate mixed-boundary eigenvalues: those estimates follow from
the analytic form, trace, and min--max arguments in the manuscript.

Every acceptance check raises an exception on failure, including under -O.
Histogram entries use unsigned storage; all count totals and every inequality
are evaluated as unbounded Python integers or exact fractions. No floating
point arithmetic is used. The largest calculation uses about 128 MiB for its
two histograms, plus the interpreter. Each type is released before the next.
"""

from array import array
from fractions import Fraction as Q
from math import isqrt
import json


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def report(section, **data):
    print(json.dumps(dict(section=section, passed=True, **data), sort_keys=True), flush=True)


def atan_bounds(reciprocal, terms):
    """Consecutive partial sums enclose atan(1/reciprocal), strictly."""
    partial = sum((Q((-1) ** j, (2 * j + 1) * reciprocal ** (2 * j + 1))
                   for j in range(terms)), Q(0))
    adjacent = partial + Q((-1) ** terms,
                           (2 * terms + 1) * reciprocal ** (2 * terms + 1))
    return min(partial, adjacent), max(partial, adjacent)


def verify_pi_and_boundary_constants():
    low5, high5 = atan_bounds(5, 12)
    low239, high239 = atan_bounds(239, 4)
    # Machin: pi = 16 atan(1/5) - 4 atan(1/239).
    lower = 16 * low5 - 4 * high239
    upper = 16 * high5 - 4 * low239
    require(lower > Q(103993, 33102), "strong rational pi lower bound")
    require(Q(103993, 33102) > Q(333, 106) > Q(25, 8),
            "nested rational pi lower bounds")
    require(upper < Q(22, 7) and upper * upper < 10,
            "pi upper bounds")

    # Compact ground lift: 2 pi^2/(pi^2+2) > pi^2/6.
    require(upper * upper + 2 < 12, "compact ground lift")
    # Compact opposite-port lift: 2/(1+pi sqrt(13)) > 1/32.
    require(13 < 4 ** 2 and 1 + Q(22, 7) * 4 < 64,
            "opposite-port lift")

    # Entire complementary-spectrum Gram bound for the degree-three cluster.
    require(Q(2, 3) < Q(25, 36), "sqrt(2/3) < 5/6")
    require(Q(2, 3) + Q(22, 7) * Q(5, 6) == Q(23, 7),
            "high transverse-mode sum bound")
    require(Q(16, 5) * (Q(10, 6) - 1) == Q(32, 15),
            "zeroth transverse-mode sum bound")
    gram_a = Q(64, 625) ** 2 * Q(32, 15) + Q(64, 625) * Q(23, 7) * Q(3, 20)
    require(gram_a < Q(73, 1000), "complete Gram complement bound")
    require(Q(1, 12) + Q(73, 1000) == Q(469, 3000), "first Gram diagonal")
    require(2 * (Q(73, 1000) + Q(3, 8) * Q(64, 625)) == Q(557, 2500),
            "second Gram diagonal")
    normalized_trace = 4 * Q(469, 3000) + Q(557, 2500)
    require(normalized_trace == Q(6361, 7500), "normalized Gram trace")
    matrix_margin = Q(17, 20) - normalized_trace
    require(matrix_margin == Q(7, 3750) and matrix_margin > 0, "Gram trace margin")
    # PSD bounds its largest eigenvalue by its trace, hence G <= (17/20)D^2.
    # ||u||^2 <= (17/20 + 1/pi^2) Q(h); pi^2 < 10 then suffices.
    require(Q(20, 17 * 10 + 20) == Q(2, 19), "uniform cluster increment")
    rational_lift = Q(25, 16) + Q(2, 19)
    require(rational_lift == Q(507, 304), "strengthened cluster lift")
    require(rational_lift - Q(5, 3) == Q(1, 912), "stronger intermediate lift")
    require(Q(5, 3) - Q(53, 32) == Q(1, 96), "uniform lift endpoint gap")
    require(Q(25, 16) + Q(3, 32) == Q(53, 32), "cluster lift endpoint")
    budget_margin = Q(53, 32) * Q(333, 106) - Q(26, 5)
    require(budget_margin > 0, "cluster endpoint area comparison")

    # All possible indices below the endpoint have m,n <= 2. The larger
    # auxiliary lists also check the next levels and the cluster multiplicity.
    even = sorted(16 * m * m + (2 * n + 1) ** 2
                  for m in range(3) for n in range(4))
    odd = sorted(4 * (2 * m + 1) ** 2 + (2 * n + 1) ** 2
                 for m in range(3) for n in range(4))
    require(even[:6] == [1, 9, 17, 25, 25, 41], "even-sector ordering")
    require(odd[:3] == [5, 13, 29], "odd-sector ordering")
    require([v for v in even if 2 * v <= 53] == [1, 9, 17, 25, 25],
            "even-sector inclusive threshold")
    require([v for v in odd if 2 * v <= 53] == [5, 13],
            "odd-sector inclusive threshold")

    # Offset-chain trace masses: g1 > 16/75 and g2 > 4/3.
    require(3 < Q(7, 4) ** 2, "sqrt(3) < 7/4")
    require(Q(4, 3) - 2 * Q(7, 4) / Q(25, 8) == Q(16, 75),
            "first chain trace mass")
    require(63 * Q(16, 75) > Q(22, 7) * 4, "first chain spectral lift")
    require(10 < 4 ** 2 and 19 * Q(4, 3) > 2 * Q(22, 7) * 4,
            "second chain spectral lift")
    report("pi_and_boundary_constants", pi_lower="103993/33102",
           pi_upper="22/7", pi_squared_upper="10",
           gram_majorant=str(gram_a), gram_majorant_margin=str(Q(73, 1000) - gram_a),
           normalized_gram_trace=str(normalized_trace), matrix_margin=str(matrix_margin),
           common_gram_coefficient="17/20", strengthened_cluster_lift=str(rational_lift),
           strengthened_endpoint_gap=str(rational_lift - Q(53, 32)),
           cluster_endpoint_budget_margin=str(budget_margin),
           even_sector_first_six=even[:6], odd_sector_first_three=odd[:3])


def block_numerator(j, denominator, size):
    return min(sum(((ell - k) * j) % denominator for k in range(size))
               for ell in range(size))


def integral_lower(a, size, slope_limit, grid=8192, root_scale=10 ** 9):
    """Integer lower sum for (1/size) integral M_size(g_a(x)) dx.

    b_j = floor(root_scale * j/sqrt(a^2*grid^2+j^2)) is obtained
    exactly by isqrt of the integer quotient. Width and height are
    clipped below at zero before multiplication; no rounding is hidden.
    """
    lipschitz = size * (size - 1) // 2
    total = 0
    previous = 0
    for j in range(grid * slope_limit):
        following_j = j + 1
        following = isqrt(root_scale ** 2 * following_j ** 2 //
                          (a * a * grid * grid + following_j ** 2))
        height = max(0, block_numerator(j, grid, size) - lipschitz)
        width = max(0, following - previous - 1)
        total += height * width
        previous = following
    return Q(total, size * grid * root_scale)


def verify_compact_tail():
    chamber = integral_lower(1, 6, 2)
    corridor = integral_lower(4, 4, 4)
    require(chamber == Q(860001557153, 3276800000000), "chamber lower-sum value")
    require(corridor == Q(138748241337, 819200000000), "corridor lower-sum value")
    require(chamber > Q(13, 50), "chamber integral threshold")
    require(corridor > Q(169, 1000), "corridor integral threshold")
    require(Q(2, 27) < Q(9, 121), "quarter-circle remainder constant")
    require(2 < Q(100, 49), "sqrt(2) bound")
    require(Q(3, 11) * (Q(10, 7) + 6) == Q(156, 77), "tail root coefficient")
    require(Q(831 * 3, 2000) - Q(38, 25) == -Q(547, 2000),
            "degree-three linear coefficient")
    require(Q(6**2 * 2, 8) + Q(6, 2) == 12, "centered chamber block error")
    require(Q(4**2 * 4, 8) + Q(4, 2) == 10, "centered corridor block error")
    require(13 + 5 * 3 == 28, "degree-three tail constant")
    require(800 < 29 ** 2, "sqrt(800) bound")
    endpoint = -Q(547, 2000) * 800 + Q(156, 77) * 29 + 28
    require(endpoint == -Q(50838, 385) and endpoint < 0, "compact tail endpoint")
    # Dividing the displayed analytic majorant by r leaves a negative
    # constant plus positive multiples of r^(-1/2) and r^(-1).
    report("compact_tail", chamber_integral_lower=str(chamber),
           corridor_integral_lower=str(corridor), radius_cutoff=800,
           endpoint_upper=str(endpoint))


def verify_degree_four_block_tail():
    chamber = integral_lower(1, 13, 3, grid=32768, root_scale=10**12)
    corridor = integral_lower(4, 7, 8, grid=32768, root_scale=10**12)
    require(chamber == Q(37632015538794073, 106496000000000000),
            "degree-four chamber lower-sum value")
    require(corridor == Q(32022499598134509, 114688000000000000),
            "degree-four corridor lower-sum value")
    require(chamber > Q(353, 1000) and corridor > Q(279, 1000),
            "degree-four block integral thresholds")
    chamber_error = Q(13**2 * 3, 8) + Q(13, 2)
    corridor_error = Q(7**2 * 8, 8) + Q(7, 2)
    require(chamber_error == Q(559, 8) and corridor_error == Q(105, 2),
            "degree-four centered-block errors")
    require(2*Q(353, 1000) + 2*Q(279, 1000) - 1 == Q(33, 125),
            "degree-four linear gain")
    require(Q(3, 11) * (Q(10, 7) + 8) == Q(18, 7),
            "degree-four root coefficient")
    constant = 1 + chamber_error + 2*corridor_error
    require(constant == Q(1407, 8), "degree-four total constant")
    require(1024 == 32**2, "degree-four endpoint root")
    endpoint = -Q(33, 125)*1024 + Q(18, 7)*32 + constant
    require(endpoint == -Q(85227, 7000) and endpoint < 0,
            "degree-four tail endpoint")
    report("degree_four_block_tail", chamber_integral_lower=str(chamber),
           corridor_integral_lower=str(corridor), radius_cutoff=1024,
           endpoint_upper=str(endpoint))


COMPACT_TYPES = (
    # name, first x, first y, degree, event count, min margin, argmin, raw failures
    ("d0", 2, 2, 0, 1028343, 604624, 8, []),
    ("d1", 2, 1, 1, 1390540, 500631, 5, []),
    ("d2_opposite", 2, 0, 2, 1028657, 440576, 80, [20]),
    ("d2_adjacent", 1, 1, 2, 1402109, 222632, 26, [2]),
    ("d3", 1, 0, 3, 1391119, 198901, 89, [1, 25]),
)


def histogram(limit, first_x, first_y, step_x=2, step_y=2):
    """All events x^2+y^2 <= limit, retaining ordered-pair multiplicity."""
    hist = array("I", [0]) * (limit + 1)
    require(hist.itemsize >= 4, "histogram storage must have at least 32 bits")
    # At the stated cutoff even the containing rectangle has < 2^32
    # pairs, so no entry can overflow. Totals below use Python integers.
    root = isqrt(limit)
    rectangular_bound = ((root - first_x) // step_x + 1) * ((root - first_y) // step_y + 1)
    require(0 <= rectangular_bound < 2 ** 32, "histogram entry capacity")
    for x in range(first_x, root + 1, step_x):
        xx = x * x
        for y in range(first_y, isqrt(limit - xx) + 1, step_y):
            hist[xx + y * y] += 1
    return hist


def compact_correction(name, p, chamber_count):
    correction = 0
    if name in ("d2_adjacent", "d3") and 3 * p < 8 and chamber_count > 0:
        correction += 1
    if name == "d2_opposite" and 20 <= p and 2 * p < 41:
        correction += 1
    if name == "d3" and 25 <= p and 2 * p < 53:
        correction += 2
    return correction


def verify_compact_events():
    radius = 800
    limit = 16 * radius * radius
    passage_hist = histogram(limit, 0, 4, 1, 4)
    switch_spec = {
        "d0": [], "d1": [],
        "d2_opposite": [(Q(41, 2), Q(53314), 16)],
        "d2_adjacent": [(Q(8, 3), Q(149984, 3), 2)],
        "d3": [(Q(8, 3), Q(981928, 3), 2), (Q(53, 2), Q(17281, 2), 26)],
    }
    for name, sx, sy, degree, expected_events, expected_minimum, expected_p, expected_raw in COMPACT_TYPES:
        chamber_hist = histogram(limit, sx, sy)
        chamber_count = passage_count = event_count = 0
        minimum = None
        raw_failures = []
        # These three snapshots suffice for every noninteger switch.
        snapshots = {}
        for p in range(limit + 1):
            ch, rh = chamber_hist[p], passage_hist[p]
            chamber_count += ch
            passage_count += rh
            if p in (2, 20, 26):
                snapshots[p] = (chamber_count, passage_count)
            if not (ch or rh):
                continue
            event_count += 1
            count2 = 2 * chamber_count + degree * passage_count
            rhs = 103993 * (2 + degree) * p
            if rhs < 16 * 33102 * count2:
                raw_failures.append(p)
            correction = compact_correction(name, p, chamber_count)
            require(chamber_count >= correction, name + ": nonnegative chamber count")
            margin = rhs - 16 * 33102 * (count2 - 2 * correction)
            require(margin >= 0, name + ": actual event at P=" + str(p))
            if minimum is None or margin < minimum[0]:
                minimum = (margin, p)
        require(event_count == expected_events, name + ": event count")
        require(minimum == (expected_minimum, expected_p), name + ": minimum margin")
        require(raw_failures == expected_raw, name + ": complete coarse failure list")
        switches = []
        for p, expected_margin, expected_count2 in switch_spec[name]:
            cc, rc = snapshots[p.numerator // p.denominator]
            correction = compact_correction(name, p, cc)
            require(correction == 0, name + ": correction reset at " + str(p))
            count2 = 2 * (cc - correction) + degree * rc
            margin = 103993 * (2 + degree) * p - 16 * 33102 * count2
            require(margin >= 0 and margin == expected_margin and count2 == expected_count2,
                    name + ": artificial event at P=" + str(p))
            switches.append(dict(p=str(p), twice_upper_count=count2, margin=str(margin)))
        # The common finite/tail endpoint is an actual corridor event (a=0,b=800).
        require(passage_hist[limit] > 0, "compact finite/tail endpoint included")
        report("compact_events", type=name, cutoff_p=limit,
               actual_events=event_count, minimum_integer_margin=minimum[0],
               minimizing_p=minimum[1], coarse_failing_p=raw_failures,
               switches=switches)
        del chamber_hist


CHAIN_WINDOWS = {
    2: [(125, 351), (530, 658), (748, 808), (852, 882), (905, 921),
        (935, 944), (953, 958), (964, 968), (972, 974), (978, 979)],
    4: [(63, 184), (299, 400), (491, 566), (631, 683), (729, 764),
        (796, 820), (843, 860), (876, 888), (900, 909), (918, 925)],
}


def verify_chain_tail():
    total_width = 0
    window_count = 0
    for a, windows in CHAIN_WINDOWS.items():
        previous = 0
        for j, (p, q) in enumerate(windows):
            require(previous < p < q < 1000, "chain window order")
            low, high = 4 * j + 1, 4 * j + 3
            require((16 * a * a + low * low) * p * p >= low * low * 10 ** 6,
                    "chain window lower derivative")
            require((16 * a * a + high * high) * q * q <= high * high * 10 ** 6,
                    "chain window upper derivative")
            total_width += q - p
            window_count += 1
            previous = q
    require(window_count == 20 and total_width == 934, "chain window totals")
    require(Q(8, 3) < Q(25, 9), "chain root remainder coefficient")
    require(17 ** 2 <= 300 and 5 <= Q(300, 60), "chain endpoint bounds")
    margin = Q(934, 8000) - Q(1, 60) - Q(5, 51)
    require(margin > 0, "chain tail coefficient margin")
    # Independent exact endpoint/monotonicity check of the same tail:
    slope = Q(total_width, 8000)
    affine = 300 * slope - 5
    require(affine > 0 and affine ** 2 > Q(8 * 300, 3), "chain tail at radius 300")
    require(2 * slope * affine > Q(8, 3), "chain squared-tail monotonicity")
    report("chain_tail", derivative_windows=window_count,
           total_window_width=str(Q(total_width, 1000)), radius_cutoff=300,
           coefficient_margin=str(margin), endpoint_affine_lower=str(affine))


def verify_chain_events():
    radius = 300
    limit = 32 * radius * radius
    hist = array("I", [0]) * (limit + 1)
    require(hist.itemsize >= 4, "chain histogram storage")
    require(600 * 301 + 300 * 1201 < 2 ** 32, "chain histogram entry capacity")
    # Chamber spectrum X=8m^2+32n^2, m>=1,n>=0.
    for m in range(1, isqrt(limit // 8) + 1):
        mm = 8 * m * m
        for n in range(isqrt((limit - mm) // 32) + 1):
            hist[mm + 32 * n * n] += 1
    # Corridor spectrum X=32m^2+2n^2, m>=1,n>=0.
    for m in range(1, isqrt(limit // 32) + 1):
        mm = 32 * m * m
        for n in range(isqrt((limit - mm) // 2) + 1):
            hist[mm + 2 * n * n] += 1
    count = actual_events = total_events = 0
    minimum = None
    raw_failures = []
    switches = []
    for x in range(limit + 1):
        multiplicity = hist[x]
        count += multiplicity
        if multiplicity:
            actual_events += 1
        if not multiplicity and x not in (41, 168):
            continue
        total_events += 1
        correction = int(40 <= x < 41 or 160 <= x < 168)
        require(count >= correction, "chain nonnegative upper count")
        corrected = count - correction
        margin = 999 * x - 6784 * corrected
        require(margin >= 0, "chain event at X=" + str(x))
        if minimum is None or margin < minimum[0]:
            minimum = (margin, x, corrected)
        if 999 * x < 6784 * count:
            raw_failures.append(x)
        if x in (41, 168):
            require(correction == 0, "chain reset is inclusive")
            switches.append(dict(x=x, upper_count=corrected, margin=margin))
    require(actual_events == 135424, "chain actual event count")
    require(total_events == 135426, "chain total event count")
    require(count == 423767, "chain cutoff multiplicity")
    require(raw_failures == [40, 160], "chain coarse failure list")
    require(minimum == (255, 41, 6), "chain minimum margin")
    require(hist[limit] > 0, "chain finite/tail endpoint included")
    report("chain_events", cutoff_x=limit, actual_jump_levels=actual_events,
           total_events_including_switches=total_events, coarse_count_at_cutoff=count,
           coarse_failing_x=raw_failures, minimum_integer_margin=minimum[0],
           minimizing_x=minimum[1], minimizing_upper_count=minimum[2], switches=switches)


def main():
    import verify_degree_four_odd33
    import verify_degree_four_boundary
    import verify_degree_four_tail
    import verify_nonuniform

    verify_pi_and_boundary_constants()
    verify_compact_tail()
    verify_degree_four_block_tail()
    verify_compact_events()
    verify_degree_four_odd33.main()
    verify_degree_four_boundary.main()
    verify_degree_four_tail.verify_events()
    verify_nonuniform.main()
    verify_chain_tail()
    verify_chain_events()
    report("all_certificates")


if __name__ == "__main__":
    main()
