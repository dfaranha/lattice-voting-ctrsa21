# Parameter selection

This note records how the parameters in `param.h` were chosen, so that they can
be re-derived rather than copied. All figures are classical core-SVP
(`0.292 * beta`, the ADPS16 model) unless stated otherwise; the MATZOV
gate-count figures are given alongside because the two models disagree by
roughly 25 to 29 bits here and the published claim depends on which one is
used.

Estimates come from the reference
[lattice-estimator](https://github.com/malb/lattice-estimator), which was
checked first against the published figures for Kyber (118 / 181 / 254) and for
Dilithium's MSIS (123 / 186 / 265) and reproduces both.

## What the published parameters actually give

The paper states "at least 100 bits of security". Resolving that:

| Problem | Protects | core-SVP | MATZOV |
| --- | --- | --- | --- |
| MLWE, commitment hiding | ballot privacy, ZK of the shuffle | **72** | 100 |
| MSIS, commitment binding | soundness of the shuffle | 117 | 144 |
| MLWE, encryption | ballot secrecy | 87 | 116 |

So the claim holds under MATZOV, where the weakest link is 100.5, and does not
hold under core-SVP, which is the model Kyber and Dilithium quote. The binding
constraint is hiding, not soundness.

Two things make hiding weak. The commitment randomness is ternary
(`BETA = 1`), and the MLWE rank is `WIDTH - HEIGHT - 1`, which at `WIDTH = 3`
is 1. A rank-1 instance of dimension 1024 modulo a 32-bit prime has a noise
rate of about `2^-32`; Kyber-512 has half the dimension but a noise rate of
`2^-11.4`, and is far harder. The commitment key is also fixed and reused
across every ballot and every election, so the one-time cost of reducing that
lattice is amortised over all of them.

## The constraints

Writing `N = DEGREE`, `k = WIDTH`, `n = HEIGHT`, `nu = NONZERO`:

1. **Hiding** is MLWE of rank `k - n - 1` over `R_p`, with `n + 1` samples and
   ternary secret and error.
2. **Binding** is MSIS over `R_p` with `n` rows, `k` columns and bound
   `16 * SIGMA_C * sqrt(nu * N)` (paper, section 3.1). Note this is the bound
   the reduction gives, which is about 9 times larger than the `2 * sqrt(k) *
   2 * sqrt(N) * SIGMA_C` that the verifier's own norm check suggests; using
   the latter overstates binding security by some 71 bits (239.7 against
   168.2 at the current `p`; it was 191.0 against 132.9 at the previous one).
   Section 5 of SOUNDNESS.md appeals to MSIS at the larger
   `2 * ||delta||_1 * 4 * sqrt(N) * SIGMA_C`, with `||delta||_1 <= 2 * nu`,
   which is `sqrt(nu)` times the paper's bound: `2^28.07` against `2^25.91`.
   Which of the two the extraction needs has not been settled, so `p` is
   sized to clear 128 bits at both.
3. **Invertibility (Lemma 1)**: the extracted `||(d - d') * sigma_i||` must
   stay below `sqrt(p / 2)`, which floors `p` at `32 * N * SIGMA_S^2`
   `= 2^28.20`. The current `p` clears it by a factor of 446 (it was 8.3 at
   `2^31.25`): `sqrt(p / 2) = 262144`, against an extracted norm of at most
   12416.
4. **Decryption correctness**: `q > 2p(2 * DIM * N^2 * BETA^2 + N + 1)`
   (paper, Table 1). The factor is 12584962, about `2^23.6`, so `q` scales
   linearly with `p`: at the current `p` the bound is `2^60.585`, and `q` is
   `2^60.785`, 1.15 times it.
5. **Quasi-unique responses**: `q > 24 * SIGMA_E^2 = 2^36`. Never binding.
6. **Challenge space**: Theorem 1 has a `4 * tau * t` term, so `|C|` should
   exceed about `2^135`. The sampler in `shuffle.c` draws challenges with `nu`
   coefficients of `+-1`, so `|C| = 2^nu * C(N, nu)`, which is `2^158.7` at
   `nu = 20`.
7. Both `p` and `q` must be prime and congruent to 5 mod 8, so that
   `x^N + 1` splits into exactly `NCRT = 2` factors.

The interesting interaction is that constraints 1 and 2 pull `p` in opposite
directions, and 3 and 4 pin it from below and tie `q` to it. Raising `p`
strengthens binding but weakens hiding, and through constraint 4 it raises `q`
and so weakens encryption as well; it also widens the uniform coefficients of
the proof. Hiding is no longer what holds `p` down: at rank 2 it has surplus
(151.9 bits at the current `p`). What sets `p` is binding at Section 5's
conservative bound, which `p = 2^37` lifts to 134.3, and what caps it is
encryption, which `q` pulls down to 136.7. At the paper's bound binding has
surplus (168.2) and the level is set by encryption; at the conservative bound
it is set by binding. Raising `p` further would buy binding with encryption.

## What changed

| | published | here |
| --- | --- | --- |
| `WIDTH` | 3 | 4 |
| `NONZERO` | 36 | 20 |
| `MODP` | 3906450253 (2^31.86) | 137438953501 (2^37.00) |
| `SIGMA_C` | 54000 | 27501 |
| `SIGMA_S` | 256 | 97 |
| `DIM` | 2 | 3 (now derived as `WIDTH - HEIGHT`) |
| `VECTOR` | 3 | 4 (now derived as `WIDTH`) |
| `Q` | 2^56 | 1986862199687640941 (2^60.79) |

`MODP` and `Q` were 2553802997 (2^31.25) and 36028797018964429 (2^55) until
they were raised for binding; see "Raising p" below.

`WIDTH = 4` is what takes the MLWE rank from 1 to 2. There is nothing in
between: rank 1 gives 72 bits and rank 2 gave 186 at the 2^31.25 modulus
(151.9 at the current 2^37), so hiding ends up over-provisioned, and the only way to trade that surplus back would be a
smaller degree with a higher rank, which cuts against raising `DEGREE` later.

`NONZERO = 20` pays for most of `WIDTH = 4`. The MSIS bound grows as
`nu^1.5` (once through `SIGMA_C`, which is linear in `nu`, and again through
the explicit `sqrt(nu)`), so cutting `nu` from 36 improves binding and shrinks
the responses at the same time. At 36 the challenge space was `2^257.0`, far more
than soundness needs.

Both Gaussian widths are the smallest that rejection sampling at the code's
`M = 1.75` allows. Lyubashevsky's bound keeps the statistical distance at
`2^-100` when `M >= exp(12 / alpha + 1 / (2 alpha^2))` for `alpha = sigma /
T`, where `T` bounds the norm of the shift being masked; at `M = 1.75` that
is `alpha >= 21.485`. The challenge `d` has `nu` coefficients of `+-1`, so:

- `SIGMA_C = ceil(21.485 * nu * sqrt(k * N)) = 27501`, since `||d * r|| <=
  nu * BETA * sqrt(k * N) = 1280` for ternary `r`. This replaces the 11 of
  Table 1, which is not enough for `M = 1.75`, and the 27.06 the code used to
  carry over by rescaling the paper's 54000, which left 21 % of slack.
- `SIGMA_S = ceil(21.485 * sqrt(nu)) = 97`, since the `sigma_i` are monomials
  and `||d * sigma_i|| = ||d|| = sqrt(nu)`. The Lemma 1 floor on `p` is then
  `32 * N * SIGMA_S^2 = 2^28.20`, well below `p = 2^37`.

Narrowing the widths does not change the prover's running time. A rejection
test whose ratio stays below `M` accepts with probability `1 / M` whatever
`alpha` is, and the linear proof of the shuffle runs four of them (the three
`SIGMA_C` masks and the `SIGMA_S` one) and repeats until all four accept, so
it expects `M^4 = 9.38` repetitions either way. Instrumenting `lin_prover`
over 150 linear proofs gave 10.05 repetitions and a per-test acceptance rate
of 0.576 at these widths, against 9.09 and 0.563 at the previous 34641 and
190, both consistent with `1 / M = 0.571`.

`SIGMA_E` is left at 54000. Table 1 would allow 28160 at these dimensions, but
the decryption bound in `vericrypt.c` has not been re-derived for a smaller
value, and 54000 is the conservative direction.

## Result

| Problem | published | here | |
| --- | --- | --- | --- |
| | core-SVP | core-SVP | MATZOV |
| MLWE, commitment hiding | 72 | **151.9** | 177.9 |
| MSIS, binding at the paper's `16 * SIGMA_C * sqrt(nu * N)` | 117 | **168.2** | 193.1 |
| MSIS, binding at Section 5's `2 * ||delta||_1 * 4 * sqrt(N) * SIGMA_C` | 89.9 | **134.3** | 160.7 |
| MLWE, encryption | 87 | **136.7** | 164.6 |

The level is 136.7 core-SVP at the paper's bound, set by encryption, and 134.3
at the conservative one, set by binding. At the previous `p = 2^31.25` the
same columns read 186.1 (213.9) for hiding, 132.9 (159.2) and 105.4 for
binding, and 155.4 (183.6) for encryption, so the level was 132.9 at the
paper's bound and 105.4 at the conservative one. Binding at the paper's bound
was 128.2 (154.7) at the earlier `SIGMA_C = 34641`: the narrower width shrinks
that bound from `2^26.24` to `2^25.91`.

Cost: the shuffle proof is 50.0 KB per message at `MSGS = 25` (1280192 bytes),
measured by the serialisation test, at 16.26 bits a Gaussian coefficient
against 18.38 for a flat encoding at the bound and 16.17 of entropy. It was
46.3 KB at the 2^31.25 modulus: the Gaussian coefficients do not depend on
`p`, and only the uniform ones widen, from 32 to 38 bits. Before that it was
46.4 KB (1188970 bytes, at 16.25 bits a Gaussian coefficient) when it still
carried the unused last `s_i`, against 47.1 KB at the previous widths; a
ciphertext grows by one ring element over the published parameters, and its
coefficients, which live modulo `q`, widen with `q` from 56 to 61 bits.
Prover time is unchanged: the linear proof expects `M^4 = 9.38` repetitions
whatever the widths are.

## Raising p (2026-10-10)

At `p = 2553802997` (2^31.25), binding was 132.9 bits at the paper's bound
`16 * SIGMA_C * sqrt(nu * N) = 2^25.91`, but SOUNDNESS.md, Section 5, proves
effective messages well defined only under MSIS at
`2 * ||delta||_1 * 4 * sqrt(N) * SIGMA_C = 2^28.07`, larger by `sqrt(nu)`,
where it was 105.4. So this branch reached 128 bits only if the paper's bound
is the one extraction needs, which has not been settled.

The levers scanned, by core-SVP binding at the conservative bound:

| change | binding (conservative) | cost |
| --- | --- | --- |
| `p = 2^33` | 114.2 | |
| `p = 2^35` | 124.1 | |
| `p = 2^35`, `nu = 17` | 130.2 | `|C| = 2^138.5`, thin against the `2^135` of constraint 6 |
| `p = 2^37` | **134.3** | uniform coefficients 32 -> 38 bits, `q` to `2^60.79` |
| `HEIGHT = 2`, `WIDTH = 5` | 155.3 | about +14 KB a message, and the code assumes `HEIGHT = 1` |

Raising `p` is the cheapest: `SIGMA_C` and `SIGMA_S` do not depend on it, so
the masked openings, which are most of the proof, do not grow. The new `p` is
137438953501, prime and 5 mod 8, with square roots of -1 `P0 = 11860615386`
and `P1 = 125578338115`. `q` follows it through constraint 4:
1986862199687640941, about `2^60.785`, against the `2^60.585` required, also 5
mod 8, with `Q0 = 403703418285289875` and `Q1 = 1583158781402351066`, and far
above `24 * SIGMA_E^2 = 2^36.03`. The challenge space is unchanged at
`2^158.7`. Coefficients now hash in 5 bytes, the fewest that hold `p - 1`.

Side effects: hiding falls from 186.1 to 151.9 and encryption from 155.4 to
136.7, both still above binding at the conservative bound; the Lemma 1 margin
of constraint 3 grows from 8.3 to 446, and the `SIGMA_S` window of
SOUNDNESS.md, Section 8, widens from a factor of 4.1 to 30.1; the proof grows
from 46.3 to 50.0 KB a message.

## Caveats

- The decryption bound in constraint 4 uses the worst-case ring product bound
  `||ab||_inf <= N ||a||_inf ||b||_inf` applied twice, hence the `N^2`. A
  statistical bound would be much smaller and would let `q` drop further, which
  would shrink ciphertexts and raise encryption security at the same time. That
  needs a failure-probability analysis, which has not been done here.
- These are estimates against the best currently known attacks, in a model that
  counts one SVP call. They are not proofs.
