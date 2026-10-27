from math import sqrt
vname = "validity"
# The voter's proof that its commitment holds a valid vote (validity.c),
# for LaZer's lin-codegen.sage. Over R_p, p = MODP = 2^39 - 19, degree 1024:
#
#   c1 = sum_j B1_j r_j
#   c2 = sum_j b2_j r_j + m
#   0  = m - u
#
# with witness w = (r_0, ..., r_3, m, u): r short, m binary and ||u|| <= 1, so
# that m is X^j or 0, a vote for one candidate or none.
deg   = 1024
mod   = 549755813869
dim   = (3, 6)
wpart = [[0, 1, 2, 3], [4], [5]]
# r is uniform ternary: E ||r||^2 = 4 * 1024 * 2/3, and 1.2 times that norm
# is about 40 standard deviations above it.
wl2   = [1.2*sqrt(4*1024*2/3), 0, 1]
wbin  = [0, 1, 0]
wlinf = 1
# LaZer's default target, a root Hermite factor of 1.0044, is about 102 bits
# core-SVP (lattice-estimator); 1.0036 gives 128 or more.
DELTA_TARGET = 1.0036
