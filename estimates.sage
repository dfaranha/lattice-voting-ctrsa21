# Security estimates for every branch, as the README's table quotes them.
#
# Run with the lattice-estimator (https://github.com/malb/lattice-estimator)
# importable, for instance from a checkout next to this file:
#
#   git clone https://github.com/malb/lattice-estimator
#   sage estimates.sage
#
# or with LATTICE_ESTIMATOR pointing at a checkout elsewhere. Figures are
# core-SVP bits (the ADPS16 model); the README's were last produced with
# lattice-estimator 53da598.
#
# The three problems, per branch:
#
# - Binding: MSIS over R_p with HEIGHT rows, WIDTH columns and bound
#   16 * sigma * sqrt(nu * N), the extractable opening of the paper's section
#   3.1, with sigma the width of the masks the extractor sees: SIGMA_C, or on
#   lnp SIGMA_B, since one rejection test there spans the whole batch. The
#   last column is the same MSIS at 2 ||delta||_1 * 4 sqrt(N) * sigma, with
#   ||delta||_1 <= 2 nu, the bound Section 5 of SOUNDNESS.md on the patched
#   branches uses; it is larger by a factor sqrt(nu), and which of the two the
#   extraction needs has not been settled.
# - Hiding: MLWE over R_p of rank WIDTH - HEIGHT - 1 with HEIGHT + 1 samples,
#   secret and error uniform in {-1, 0, 1}, as commit_sample_short draws them.
# - Encryption: MLWE over R_q of rank DIM with DIM samples, secret and error
#   centred binomial with parameter 1, as encrypt_sample_short draws them.

import os
import sys
from math import log2, sqrt

sys.path.insert(0, os.environ.get("LATTICE_ESTIMATOR", "lattice-estimator"))
from estimator import LWE, SIS, ND, RC

TERNARY = ND.Uniform(-1, 1)
BINOMIAL = ND.CenteredBinomial(1)


def lwe(n, q, m, Xs, Xe):
    params = LWE.Parameters(n=n, q=q, Xs=Xs, Xe=Xe, m=m)
    costs = LWE.estimate.rough(params, quiet=True)
    return float(log2(min(c["rop"] for c in costs.values())))


def sis(N, k, n, q, beta):
    params = SIS.Parameters(n=n * N, q=q, length_bound=beta, m=k * N, norm=2)
    return float(log2(SIS.lattice(params, red_cost_model=RC.ADPS16)["rop"]))


# DEGREE, MODP, WIDTH, HEIGHT, Q, DIM, NONZERO and the mask width the binding
# bound follows, as each branch's param.h, commit.h and encrypt.h set them.
BRANCHES = {
    "master":          (1024, 3906450253, 3, 1, 72057594037928893, 2, 36, 54000),
    "fix-pkc":         (1024, 3906450253, 3, 1, 72057594037928893, 2, 36, 54000),
    "balanced-params": (1024, 137438953501, 4, 1, 1986862199687640941, 3, 20,
                        27501),
    "lnp":             (2048, 1099511627917, 3, 1, 18446744073709551557, 3, 18,
                        295465),
}

print(f"{'branch':16s} {'binding':>8s} {'hiding':>8s} {'encryption':>10s} "
      f"{'level':>7s} {'binding at 2|delta|_1 4 sqrt(N) sigma':>38s}", flush=True)
for name, (N, p, k, n, q, dim, nu, sigma) in BRANCHES.items():
    binding = sis(N, k, n, p, 16 * sigma * sqrt(nu * N))
    hiding = lwe((k - n - 1) * N, p, (n + 1) * N, TERNARY, TERNARY)
    encryption = lwe(dim * N, q, dim * N, BINOMIAL, BINOMIAL)
    rigorous = sis(N, k, n, p, 2 * (2 * nu) * 4 * sqrt(N) * sigma)
    print(f"{name:16s} {binding:8.1f} {hiding:8.1f} {encryption:10.1f} "
          f"{min(binding, hiding, encryption):7.1f} {rigorous:38.1f}", flush=True)
