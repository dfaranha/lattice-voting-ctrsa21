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
# - Binding: MSIS over R_p with HEIGHT rows and WIDTH columns, at the norm
#   of the vector the extractor builds from two relaxed openings,
#   delta' zbar - delta zbar'. Each challenge difference has
#   ||delta||_1 <= 2 nu, each of the WIDTH elements of zbar has l2-norm at
#   most 4 sigma sqrt(N), and multiplying by delta stretches a vector by at
#   most ||delta||_1, so the bound is 2 * 2 nu * 4 sigma sqrt(WIDTH * N)
#   = 16 nu sigma sqrt(WIDTH * N). sigma is the width of the masks the
#   extractor sees: SIGMA_C, or on lnp SIGMA_B, since one rejection test there
#   spans the whole batch. The paper's 16 sigma sqrt(nu N) is smaller by
#   sqrt(nu WIDTH): it takes ||delta||_2 for ||delta||_1, which does not bound
#   a product, and the norm of one element for that of the vector. It is
#   printed for reference and does not enter the level.
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
    "balanced-params": (1024, 549755813869, 4, 1, 1152921504606846869, 3, 20,
                        27501),
    "lnp":             (2048, 1099511627581, 3, 1, 18446744073709551557, 3, 18,
                        295465),
}

print(f"{'branch':16s} {'binding':>8s} {'hiding':>8s} {'encryption':>10s} "
      f"{'level':>7s} {'binding at the paper bound':>27s}", flush=True)
for name, (N, p, k, n, q, dim, nu, sigma) in BRANCHES.items():
    binding = sis(N, k, n, p, 16 * nu * sigma * sqrt(k * N))
    hiding = lwe((k - n - 1) * N, p, (n + 1) * N, TERNARY, TERNARY)
    encryption = lwe(dim * N, q, dim * N, BINOMIAL, BINOMIAL)
    paper = sis(N, k, n, p, 16 * sigma * sqrt(nu * N))
    print(f"{name:16s} {binding:8.1f} {hiding:8.1f} {encryption:10.1f} "
          f"{min(binding, hiding, encryption):7.1f} {paper:27.1f}", flush=True)
