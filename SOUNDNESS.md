# Soundness of the relaxed proof of `sigma_i in D`

## Status of this document

This is a **draft argument, not a reviewed result.** It was written to accompany
the implementation in `shuffle.c` and records why that implementation is
believed to be sound. It restates and modifies Lemma 5 of Bootle, Lyubashevsky
and Merino-Gallardo, *Efficient Verifiable Mixnets from Lattices, Revisited*
(ePrint 2025/658), and it has not been checked by anyone else. Sections 9 and 10
list what is asserted here without proof and what a full write-up would still
have to establish.

On the `lnp` branch the relaxed proof has been replaced by `is_bin`, which
puts `sigma_i` in a binary `D` as Protocol 1 does; see `LNP.md`. The figures
below are for the parameters this argument was written at, except where
marked as the `lnp` branch's.

## 1. The gap this closes

Lemma 5 of that paper proves that

```
prod_{i=1..N} (a_i + g(i) * X1 - X2) = prod_{i=1..N} (b_i + sigma_i * X1 - X2)     (*)
```

over `R[X1, X2]` implies `(a_1, ..., a_N) ~P (b_1, ..., b_N)`, **provided** each
`sigma_i` lies in a set `D` whose pairwise differences are invertible. Protocol 1
of that paper enforces `sigma_i in D` with a sub-proof of `is_bin(sigma_i)`
delegated to a general-purpose proof system.

Without that sub-proof the shuffle proof is broken, and concretely so: a prover
that swaps the first CRT component of two messages and applies the *same* swap
to the two index encodings obtains `sigma_i` satisfying (\*) in each CRT
component separately, and the proof goes through. This is implemented as the
test `shuffle proof rejects the CRT-mixing attack on sigma` in `shuffle.c`, and
it succeeded against the version of the code that lacked the countermeasure
below.

The observation exploited here is that `is_bin` is stronger than Lemma 5 needs.
The lemma only needs `sigma_i - g(j)` never to be a nonzero zero-divisor, and
Lemma 1 supplies that for any element of small norm. So it is enough to prove
that `sigma_i` is **short**, which is what Fiat-Shamir with aborts proves
natively. The price is that the resulting guarantee is *relaxed*, which is what
Section 7 has to accommodate.

## 2. Notation and parameters

| symbol | meaning | value |
| --- | --- | --- |
| `q` | modulus (`MODP`) | 3906450253, `= 5 mod 8` |
| `n` | ring degree (`DEGREE`) | 1024 |
| `R` | `Z_q[x]/(x^n + 1)` | |
| `k` | number of CRT factors | 2 |
| `p_1, p_2` | irreducible factors of `x^n + 1` | degree `n/2 = 512` |
| `sigma_C` | Gaussian width for commitment randomness (`SIGMA_C`) | 54000 |
| `sigma_S` | Gaussian width for the committed `sigma_i` (`SIGMA_S`) | 256 |
| `C_lin` | challenge set of the linear proof | `NONZERO = 36` coefficients of +-1, the rest zero |
| `C` | challenge set for `alpha, beta, gamma, lambda, tau, rho` | `{p in R : deg p < n/2}` |

By Lemma 1 with `k = 2`, `R = R_1 x R_2` with `R_l = Z_q[x]/(p_l)` a field, and

> **(Inv)** any `y` with `0 < ||y||_2 < q^(1/2) ~ 62501.6` is invertible in `R`.

Two norm facts used repeatedly:

* For `d in C_lin`: coefficients in `{-1,0,1}`, exactly 36 nonzero, so
  `||d||_2 = 6` and `||d||_1 = 36`.
* For a difference `delta = d - d'` of two challenges: coefficients in
  `{-2,...,2}`, so `||delta||_oo <= 2`, `||delta||_1 <= 72` and
  `||delta||_2 <= sqrt(||delta||_oo ||delta||_1) <= 12`. By (Inv), **every
  nonzero challenge difference is invertible**, since `12 << 62501`.

(These figures are for the challenge space of the paper. The code first drew
`d` as a difference `c_0 - c_1` of two binary challenges, for which they were
`sqrt(72)`, 72, 24 and 144. At the current `NONZERO = 18` they are
`sqrt(18) ~ 4.24`, 18, `2 sqrt(18) ~ 8.49` and 36.)

For `c, c'` distinct in `C`, `c - c'` is nonzero of degree `< n/2 = deg p_l`,
hence nonzero modulo each `p_l`, hence invertible. This is why the challenge
space was narrowed to degree `< n/2`.

## 3. Commitments and relaxed openings

The commitment is BDLOP-style: `Com(m; r) = (B1 r, <b2, r> + m)` with `r` short.

**Definition (relaxed opening).** A *relaxed opening* of `c = (c1, c2)` is a
triple `(mbar, zbar, delta)` where `delta` is a nonzero difference of two
challenges in `C_lin`, `||zbar||` is bounded as in Section 4, and

```
delta * c1 = B1 zbar ,      delta * c2 = <b2, zbar> + mbar .
```

Since `delta` is invertible, define the **effective message** of the opening as
`m* := delta^(-1) mbar in R`.

## 4. What two accepting transcripts yield

The linear proof of `shuffle.c` has first message `(t, t_p, t', v_sigma, u)`,
challenge `d`, and response `(z, z_p, z', z_sigma)`. The verifier checks

```
(1)  B1 z    = t     + d * x.c1
(2)  B1 z_p  = t_p   + d * p.c1
(3)  B1 z'   = t'    + d * _x.c1
(4)  <b2, z_p> + z_sigma = v_sigma + d * p.c2
(5)  alpha <b2,z> + gamma <b2,z_p> - <b2,z'>
         = u + d * (alpha * x.c2 + gamma * p.c2 - _x.c2 + beta)
```

together with `||z||, ||z_p||, ||z'|| <= 2 sqrt(n) sigma_C` and
`||z_sigma|| <= 2 sqrt(n) sigma_S = 16384`.

Two accepting transcripts sharing the first message with challenges `d != d'`
give, writing `delta = d - d'` and `zbar = z - z'` etc.,

```
(1')  B1 zbar   = delta * x.c1
(2')  B1 zbar_p = delta * p.c1
(3')  B1 zbar'  = delta * _x.c1
(4')  <b2, zbar_p> + zbar_sigma = delta * p.c2
(5')  alpha <b2,zbar> + gamma <b2,zbar_p> - <b2,zbar'>
          = delta * (alpha * x.c2 + gamma * p.c2 - _x.c2 + beta)
```

with `||zbar[i]|| <= 4 sqrt(n) sigma_C` for each of the `WIDTH` elements of
`zbar`, so `||zbar|| <= 4 sqrt(WIDTH n) sigma_C`, and `||zbar_sigma|| <= 4 sqrt(n) sigma_S
= 32768`.

So `(mbar_x, zbar, delta)` with `mbar_x := delta * x.c2 - <b2, zbar>` is a
relaxed opening of `x`, and similarly for `p` and `_x`. Line (4') says something
stronger about `p`, and that is the whole point:

> **(Short)** The relaxed opening of `p` extracted this way has
> `mbar_sigma = zbar_sigma`, hence `||delta * sigma*|| = ||zbar_sigma||
> <= 4 sqrt(n) sigma_S = 32768`, where `sigma* = delta^(-1) zbar_sigma` is the
> effective message of `p`.

That is the relaxed form of `sigma_i in D`: not `sigma*` itself is short, but
`delta * sigma*` is, for an invertible `delta` the extractor holds.

## 5. Effective messages are well defined

The extractor visits many leaves and obtains relaxed openings with *different*
`delta`. For the argument below to speak of "the" message of a commitment, the
effective message must not depend on which opening was used.

**Claim.** Under MSIS, if `(mbar, zbar, delta)` and `(mbar', zbar', delta')` are
relaxed openings of the same `c`, then `delta^(-1) mbar = delta'^(-1) mbar'`.

*Proof.* From `delta c1 = B1 zbar` and `delta' c1 = B1 zbar'`, multiply the
first by `delta'`, the second by `delta`, and subtract:
`B1 (delta' zbar - delta zbar') = 0`. Multiplying by a ring element `a`
stretches an `l2` norm by at most `||a||_1`, so the vector
`delta' zbar - delta zbar'` has norm at most
`2 * 2 NONZERO * 4 sqrt(WIDTH n) sigma_C = 16 NONZERO sigma_C sqrt(WIDTH n)`,
and under MSIS at that bound it is zero, i.e. `delta' zbar = delta zbar'`.
(The paper's `16 sigma_C sqrt(nu n)` is smaller by `sqrt(NONZERO WIDTH)` and
does not bound this vector: it charges a challenge its `l2` norm, which does
not bound the stretch of a product -- `(1 + x)` times the all-ones element
already exceeds it -- and charges `zbar` the norm of one element rather than
of all `WIDTH`.) Then

```
delta' delta c2 = delta' <b2,zbar> + delta' mbar = delta <b2,zbar'> + delta mbar' ,
```

and cancelling the equal inner-product terms gives `delta' mbar = delta mbar'`.
Both `delta, delta'` are invertible, so `delta^(-1) mbar = delta'^(-1) mbar'`. QED

Write `m*` for this common value. Note the MSIS bound is driven by the
`sigma_C`-sized responses; the `sigma_S`-sized response is much smaller, so
**the countermeasure does not weaken the binding assumption.**

## 6. The linear relation holds exactly for effective messages

Substituting `delta * x.c2 = <b2,zbar> + mbar_x` and its analogues into (5'),
the inner-product terms cancel and

```
alpha * mbar_x + gamma * mbar_sigma - mbar__x + delta * beta = 0 .
```

Dividing by the invertible `delta`,

```
alpha * m*_x + gamma * sigma* - m*__x + beta = 0 .
```

This is worth stating explicitly: **the relaxation does not degrade the linear
relation.** Every term is divided by the same `delta`, so the relation the
extractor obtains over effective messages is exact, not approximate. The
determinant argument of the paper's Theorem 5, and the Schwartz-Zippel steps
that follow it, therefore apply verbatim to effective messages and yield the
product identity (\*) over `R[X1, X2, X3]`, with `a_i, b_i` built from effective
messages.

## 7. Relaxed Lemma 5

**Lemma 5' (relaxed).** Let `R = R_1 x ... x R_t` be a product of integral
domains, `N in N`, and `g : [N] -> R` injective. Let
`a_1..a_N, b_1..b_N, sigma_1..sigma_N in R` satisfy (\*) over `R[X1, X2]`.
Suppose there is an invertible `delta in R` such that for all `i, j in [N]`,

```
delta * (sigma_i - g(j))  is either zero or invertible.                      (H)
```

Then `(a_1,...,a_N) ~P (b_1,...,b_N)`, and `sigma_i = g(pi(i))` for the
permutation `pi` with `b_i = a_{pi(i)}`.

*Proof.* Fix `j in [N]`. Reading (\*) as polynomials in `X2`, the element
`a_j + g(j) X1 in R[X1]` is a root of the left-hand side, hence of the right,
so `prod_i ((sigma_i - g(j)) X1 + (b_i - a_j)) = 0`. Fix `l in [t]`. Since
`R_l[X1]` is an integral domain, some factor vanishes over `R_l`: there is
`i_{j,l}` with

```
sigma_{i_{j,l}} = g(j)  and  b_{i_{j,l}} = a_j     over R_l .
```

Write `i = i_{j,l}` and suppose `sigma_i - g(j) != 0` over `R`. Since `delta` is
invertible, `delta (sigma_i - g(j)) != 0`, so by (H) it is invertible, hence
nonzero in *every* component, contradicting `sigma_i = g(j)` over `R_l`.
Therefore `sigma_i = g(j)` over `R`.

The rest is the original proof: this holds for every `l`, and `g` is injective,
so `i_{j,1} = ... = i_{j,t} =: i_j`; injectivity makes `j -> i_j` a permutation;
and `b_{i_j} = a_j` over every component, hence over `R`. QED

The only change from Lemma 5 is hypothesis (H). The original assumes
`sigma_i in D` with `D` having invertible pairwise differences, which is (H)
with `delta = 1`. Allowing a general invertible `delta` is what accommodates a
relaxed proof.

**Discharging (H).** By (Short), `||delta * sigma*_i|| <= 32768`, and
`g(j) = x^j` is a monomial, so `||delta * g(j)||_2 = ||delta||_2 <= 12`
(multiplication by a monomial is a signed rotation and preserves the l2-norm
exactly). Hence

```
|| delta * (sigma*_i - g(j)) ||_2  <=  32768 + 12  =  32780  <  62501.6 ,
```

so by (Inv) it is zero or invertible, which is (H).

## 8. Parameter budget

```
ceiling            q^(1/2)                          =  62501.6
extracted norm     4 sqrt(n) sigma_S + ||delta||_2  =  32780
margin                                                 1.91x
```

`sigma_S` is constrained from both sides:

| bound | value | source |
| --- | --- | --- |
| `sigma_S <= (q^(1/2) - 12) / (4 sqrt n)` | 488 | (H) via (Inv) |
| `sigma_S >= 21.5 * ||d * sigma_i||` | 129 | rejection sampling at `M = 1.75` |

The feasible window is roughly `[129, 488]`, a factor of 3.8 wide (it was
`[183, 488]` while `d` was a difference of two binary challenges); the
implementation uses 256. **This window is the fragile part of the construction.**
It is narrow because `k = 2`, so (Inv) only reaches `q^(1/2)`, which is small
against `n = 1024`-dimensional Gaussian responses.

*On the `lnp` branch (2026-10-09).* The figures above are for the parameters
this argument was written at. With `q = 2^40 - 195`, `n = 2048` and
`NONZERO = 18`, `||d * sigma_i|| = sqrt(18) ~ 4.24` and `||delta||_2 <= 8.49`,
so the window is

| bound | value | source |
| --- | --- | --- |
| `sigma_S <= (q^(1/2) - 8.49) / (4 sqrt n)` | 5792 | (H) via (Inv) |
| `sigma_S >= 21.485 * ||d * sigma_i||` | 91.2 | rejection sampling at `M = 1.921` |

and `SIGMA_S` was narrowed from 256 to 92, the lower edge rounded up. The
extracted norm is then `4 sqrt(n) sigma_S + ||delta||_2 = 16663` against a
ceiling of 1048576. The constant is 21.485 rather than 21.5: it is the least
`alpha` with `exp(12/alpha + 1/(2 alpha^2)) <= 1.75`, the code's `M` for
statistical distance `2^-100`; the code now uses `M = 1.921`, which keeps
`2^-128` at the same `alpha` (`exp(14/alpha + 1/(2 alpha^2))`, LNP22,
Lemma 2.14). No proof on that branch
samples at `SIGMA_S` any more, since `is_bin` replaced the relaxed proof, so
the window is recorded for completeness rather than load-bearing there.

This is also why `g` maps to monomials rather than to binary representations of
indices, as in the paper. For a binary `g`, `||d * sigma_i||_2 <= ||d||_2 *
||sigma_i||_1` grows with the bit-length of the index: at `N = 25` the window
still exists but is about a factor of 1.3 wide, and it closes entirely for
larger `N`. With monomials, `||d * sigma_i||_2 = ||d||_2` independently of `N`,
and `g` stays injective up to `N = n`.

## 9. Completeness and zero-knowledge

*Completeness.* An honest `z_sigma` has `||z_sigma|| ~ sigma_S sqrt(n) = 8192`
against a bound of 16384, so the added check fails with negligible probability.
The added rejection sampling runs at `alpha = sigma_S / ||d sigma_i|| ~ 30`,
giving `M ~ 1.49` for `2^-100`, within the `M = 1.75` used by the code. At
the `lnp` branch's `SIGMA_S = 92` and `NONZERO = 18` it is `alpha ~ 21.7` and
`M ~ 1.741`, just inside 1.75 by construction; that branch now uses
`M = 1.921`, which keeps `2^-128` at any `alpha >= 21.485` (`M ~ 1.908` at
21.7).

*Zero-knowledge.* The proof carries one extra masked value. The simulator picks
`z_p` and `z_sigma` from their distributions and sets

```
v_sigma = <b2, z_p> + z_sigma - d * p.c2 ,
```

which is the same move as for the existing masked values. **This has been
sketched, not written out**, and the SHVZK proof of the paper's Theorem 4 would
need to be extended accordingly.

## 10. What is not established here

1. **Not reviewed.** Nothing here has been checked by a second person.
2. **Theorem 5 is not reproven.** Section 6 asserts that the tree-of-transcripts
   bookkeeping of the paper's Theorem 5 carries over to effective messages. The
   argument for why it should is given, but the knowledge-error accounting has
   not been redone.
3. **MSIS at the bound Section 5 needs is weaker than the paper's.**
   Section 5 appeals to MSIS at `16 NONZERO sigma sqrt(WIDTH n)`. On this
   branch the masks the extractor sees are at `SIGMA_B = 295465`, and at
   `NONZERO = 18`, `WIDTH = 3`, `DEGREE = 2048` that bound is `2^32.64`,
   which the lattice-estimator (53da598, core-SVP) puts at **256.1 bits**
   (278.0 MATZOV). The paper's `16 sigma sqrt(nu N)` would give 327.0, but
   Section 5 explains why it does not bound the vector extraction builds;
   earlier versions of this note left the choice open, and quoted 273.6 at a
   bound that had the `l1` norm but not `sqrt(WIDTH)`. Binding is far above
   the hiding that sets this branch's level either way.
4. **Zero-knowledge sketched only**, as noted in Section 9.
5. **The batching aggregations.** The implementation proves `is_bin` for
   every message at once, aggregating the constant-coefficient and range
   relations by random scalars and the quadratic relations by random ring
   elements `rho_l`. The scalar aggregation contributes `p^-LNP_LAMBDA`. The
   `rho` aggregation contributes at most `1/|C|` per hash query:

   *Claim.* Let `F_1, ..., F_MSGS` in `R_p` be fixed before the `rho_l` are
   drawn, and not all zero. Drawing each `rho_l` independently and uniformly
   from `C`, `Pr[sum_l rho_l F_l = 0] <= 1/|C|`.

   *Proof.* Take `j` with `F_j != 0` and fix every `rho_l` with `l != j`; the
   sum vanishes only if `rho_j F_j = a` for one fixed `a`. Since
   `p = 5 mod 8`, `x^N + 1` splits into two irreducible factors and `R_p` is
   a product of two fields, and `F_j` is nonzero in at least one of them, say
   modulo `f`. There `rho_j F_j = a` has exactly one solution `rho_j mod f`.
   For distinct `rho, rho'` in `C`, `rho - rho'` has `||.||_oo <= 2`, below
   the `sqrt(p/2) = 741455` of Lemma 1, so it is invertible, hence nonzero
   modulo `f`: reducing `C` modulo `f` is injective, and at most one element
   of `C` satisfies the equation. QED

   This applies because every message's quadratic relation enters the
   aggregated check as `rho_l` times a ring element `F_l` that is zero exactly
   when the relation holds -- the garbage terms carry `rho_l` and
   `sigma(rho_l)` too, but only in the lower-order coefficients -- and the
   `F_l` are fixed by the commitments `rho_hash` absorbs before it derives the
   `rho_l`. With the `rho_l` derived by hashing, a prover who makes `Q` hash
   queries finds a set that cancels a false relation with probability at most
   `Q/|C|`; at `NONZERO = 18`, `DEGREE = 2048`, `|C| = 2^18 C(2048, 18)`,
   about `2^163.4`. What this does not do is replace the analysis of the
   challenge `d` that follows, which treats the aggregated relation as one
   statement; it bounds only the step from many relations to one. Section 5b
   of LNP-PARAMS.md has the construction.
6. **The window in Section 8 has no slack to spare.** Any change to `n`, `q`,
   `NONZERO`, `M` or the encoding `g` invalidates the budget and must be
   rechecked.

The implementation is tested against both attacks it is meant to stop, and the
verifier was instrumented to confirm each is rejected by the intended check
rather than incidentally; but tests cannot establish any of the six points
above.
