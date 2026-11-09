from math import sqrt
vname = "validity"
# The voter's proof that its commitment holds a valid vote (validity.c), for
# LaZer's lin-codegen.sage, stated over LaZer's own ring R_64 = Z_p[Y]/(Y^64 + 1),
# p = MODP = 2^39 - 19. R_p of degree 1024 is R_64[X]/(X^16 - Y), an element a
# being sum_k X^k a_k(Y); a vote for candidate j is X^(16 j) = Y^j, in piece 0
# alone. With witness w = (the 16 pieces of each of r_0, ..., r_3, m, m', u):
#
#   c1_k = (sum_j B1_j r_j)_k                      for k < 16
#   c2_k = (sum_j b2_j r_j)_k + [k = 0] m          for k < 16
#   0    = m - u
#   0    = Y^(64 - NCAND) m - m'
#
# with r short, m and m' binary and ||u|| <= 1, so that m is Y^j or 0, a vote
# for one candidate or none. The last row bounds j: R_64 is negacyclic, so
# Y^(64 - NCAND) Y^j is Y^(64 - NCAND + j), binary, when j < NCAND, and
# -Y^(j - NCAND), not binary, when j >= NCAND. NCAND is a coefficient of A, set
# at run time; the parameters here do not depend on it.
deg   = 64
mod   = 549755813869
dim   = (34, 67)
wpart = [list(range(64)), [64, 65], [66]]
# r is uniform ternary: E ||r||^2 = 4 * 1024 * 2/3, and 1.2 times that norm
# is about 40 standard deviations above it.
wl2   = [1.2*sqrt(4*1024*2/3), 0, 1]
wbin  = [0, 1, 0]
# The witness is ternary, so 1 would do; 2 makes the proof modulus 2^64 + 885
# rather than 2^63 + 1701, on which LaZer's prover does not terminate (q above
# 2^63 in one 64-bit limb). A looser bound only widens the modulus.
wlinf = 2
# LaZer's default target, a root Hermite factor of 1.0044, is about 102 bits
# core-SVP (lattice-estimator); 1.0036 gives 128 or more.
DELTA_TARGET = 1.0036
