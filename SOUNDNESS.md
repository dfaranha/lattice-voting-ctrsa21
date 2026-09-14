# Soundness of the relaxed proof of `sigma_i in D`

## Status of this document

This is a **draft argument, not a reviewed result.** It was written to accompany
the implementation in `shuffle.c` and records why that implementation is
believed to be sound. It restates and modifies Lemma 5 of Bootle, Lyubashevsky
and Merino-Gallardo, *Efficient Verifiable Mixnets from Lattices, Revisited*
(ePrint 2025/658), and it has not been checked by anyone else. Sections 9 and 10
list what is asserted here without proof and what a full write-up would still
have to establish.

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
| `C_lin` | challenge set of the linear proof | `c_0 - c_1`, each `c_j` with `NONZERO = 36` ones |
| `C` | challenge set for `alpha, beta, gamma, lambda, tau, rho` | `{p in R : deg p < n/2}` |

By Lemma 1 with `k = 2`, `R = R_1 x R_2` with `R_l = Z_q[x]/(p_l)` a field, and

> **(Inv)** any `y` with `0 < ||y||_2 < q^(1/2) ~ 62501.6` is invertible in `R`.

Two norm facts used repeatedly:

* For `d in C_lin`: coefficients in `{-1,0,1}`, at most 72 nonzero, so
  `||d||_2 <= sqrt(72) ~ 8.49` and `||d||_1 <= 72`.
* For a difference `delta = d - d'` of two challenges: coefficients in
  `{-2,...,2}`, at most 144 nonzero, so `||delta||_2 <= 24` and
  `||delta||_1 <= 144`. By (Inv), **every nonzero challenge difference is
  invertible**, since `24 << 62501`.

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

with `||zbar|| <= 4 sqrt(n) sigma_C` and `||zbar_sigma|| <= 4 sqrt(n) sigma_S
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
`B1 (delta' zbar - delta zbar') = 0`. The vector `delta' zbar - delta zbar'` has
norm at most `2 ||delta||_1 * 4 sqrt(n) sigma_C`, so under MSIS at that bound it
is zero, i.e. `delta' zbar = delta zbar'`. Then

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
`g(j) = x^j` is a monomial, so `||delta * g(j)||_2 = ||delta||_2 <= 24`
(multiplication by a monomial is a signed rotation and preserves the l2-norm
exactly). Hence

```
|| delta * (sigma*_i - g(j)) ||_2  <=  32768 + 24  =  32792  <  62501.6 ,
```

so by (Inv) it is zero or invertible, which is (H).

## 8. Parameter budget

```
ceiling            q^(1/2)                          =  62501.6
extracted norm     4 sqrt(n) sigma_S + ||delta||_2  =  32792
margin                                                 1.91x
```

`sigma_S` is constrained from both sides:

| bound | value | source |
| --- | --- | --- |
| `sigma_S <= (q^(1/2) - 24) / (4 sqrt n)` | 488 | (H) via (Inv) |
| `sigma_S >= 21.5 * ||d * sigma_i||` | 183 | rejection sampling at `M = 1.75` |

The feasible window is roughly `[183, 488]`, a factor of 2.7 wide; the
implementation uses 256. **This window is the fragile part of the construction.**
It is narrow because `k = 2`, so (Inv) only reaches `q^(1/2)`, which is small
against `n = 1024`-dimensional Gaussian responses.

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
giving `M ~ 1.49`, within the `M = 1.75` used by the code.

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
3. **MSIS parameters not re-derived.** Section 5 appeals to MSIS at a bound
   driven by `||delta||_1 * 4 sqrt(n) sigma_C`. This is the same regime as the
   original proof, and the new response is smaller, but the concrete hardness
   has not been recomputed.
4. **Zero-knowledge sketched only**, as noted in Section 9.
5. **The window in Section 8 has no slack to spare.** Any change to `n`, `q`,
   `NONZERO`, `M` or the encoding `g` invalidates the budget and must be
   rechecked.

The implementation is tested against both attacks it is meant to stop, and the
verifier was instrumented to confirm each is rejected by the intended check
rather than incidentally; but tests cannot establish any of the five points
above.
