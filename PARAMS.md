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
   `16 * nu * SIGMA_C * sqrt(k * N)`, `2^29.07`. That is the norm of the vector
   `delta' zbar - delta zbar'` that Section 5 of SOUNDNESS.md builds from two
   relaxed openings: each of the `k` elements of `zbar` has norm at most
   `4 * sqrt(N) * SIGMA_C`, so the vector has at most `sqrt(k)` times that,
   and multiplying by a challenge difference stretches it by at most
   `||delta||_1 <= 2 * nu`. The paper's `16 * SIGMA_C * sqrt(nu * N)`,
   `2^25.91`, is smaller by `sqrt(nu * k)`: it charges the challenge
   `||delta||_2`, which does not bound a product (`(1 + x)` times the
   all-ones vector already exceeds it), and an element where the vector is
   meant. Earlier versions of this note quoted binding at the paper's bound and
   at `2 * ||delta||_1 * 4 * sqrt(N) * SIGMA_C`, which has the first
   correction but not the second, and left open which one extraction needs.
   It needs this one. The paper's figure is still given, for comparison.
3. **Invertibility (Lemma 1)**: the extracted `||(d - d') * sigma_i||` must
   stay below `sqrt(p / 2)`, which floors `p` at `32 * N * SIGMA_S^2`
   `= 2^28.20`. The current `p` clears it by a factor of 1783 (446 at `2^37`,
   8.3 at `2^31.25`): `sqrt(p / 2) = 524288`, against an extracted norm of at
   most 12416.
4. **Decryption correctness**: `q >= 8 * nu * p * (2 * DIM * N + 1.5)`, which
   is `2^58.91` at the current `p`; `q` is `2^60`, twice it. This is what the
   code's decryption checks, bounded in the worst case: it computes
   `c_bar (w - <v, s_1>)` for `c_bar = c - c'`, `||c_bar||_1 <= 2 nu`, and
   accepts below `q / 4`, and for an honest ciphertext that is
   `c_bar (p x + mu)` with `||x||_inf <= 2 * DIM * N + 1`, since every secret
   and every piece of randomness is CB(1). The paper's Table 1 asks instead for
   `q > 2p(2 * DIM * N^2 * BETA^2 + N + 1)`, `2^62.585` here, which `q` was
   sized against until "Lowering q" below: its `N^2` is loose, a product of
   two ternary elements having coefficients at most `N`, and it leaves out
   `c_bar`. Either way `q` scales linearly with `p`.
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
(142.2 bits at the current `p`). What sets `p` is binding, which `p = 2^39`
lifts to 131.1, and what caps it is encryption, which `q` pulls down to 139.0,
and then hiding. Raising `p` further would buy binding with both.

## What changed

| | published | here |
| --- | --- | --- |
| `WIDTH` | 3 | 4 |
| `NONZERO` | 36 | 20 |
| `MODP` | 3906450253 (2^31.86) | 549755813869 (2^39) |
| `SIGMA_C` | 54000 | 27501 |
| `SIGMA_S` | 256 | 97 |
| `DIM` | 2 | 3 (now derived as `WIDTH - HEIGHT`) |
| `VECTOR` | 3 | 4 (now derived as `WIDTH`) |
| `Q` | 2^56 | 1152921504606846869 (2^60) |

`MODP` and `Q` were 2553802997 (2^31.25) and 36028797018964429 (2^55) until
they were raised for binding; see "Raising p" below. They were then
137438953349 (2^37) and 1986862199687640941 (2^60.79), then 2^37 and
288230376151711717 (2^58) once `q` was sized by what decryption needs ("Lowering
q"), and are now 2^39 and 2^60, for binding at the bound that extraction
actually needs ("Raising p again").

`WIDTH = 4` is what takes the MLWE rank from 1 to 2. There is nothing in
between: rank 1 gives 72 bits and rank 2 gave 186 at the 2^31.25 modulus
(142.2 at the current 2^39), so hiding ends up over-provisioned, and the only way to trade that surplus back would be a
smaller degree with a higher rank, which cuts against raising `DEGREE` later.

`NONZERO = 20` pays for most of `WIDTH = 4`. The MSIS bound grows as
`nu^2` (once through `SIGMA_C`, which is linear in `nu`, and again through
`||delta||_1`), so cutting `nu` from 36 improves binding and shrinks
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
  `32 * N * SIGMA_S^2 = 2^28.20`, well below `p = 2^39`.

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
| MLWE, commitment hiding | 72 | **142.2** | 168.8 |
| MSIS, binding at `16 * nu * SIGMA_C * sqrt(k * N)` | 83.2 | **131.1** | 157.7 |
| MSIS, binding at the paper's `16 * SIGMA_C * sqrt(nu * N)`, for comparison | 117 | 181.0 | 205.4 |
| MLWE, encryption | 87 | **139.0** | 166.7 |

The level is 131.1 core-SVP, set by binding. At `p = 2^37` and `q = 2^58` the
columns read 151.9 (177.9) for hiding, 121.8 for binding, 168.2 (193.1) at the
paper's bound, and 145.4 (172.5) for encryption, so the level was 121.8; at
`q = 2^60.79` encryption was 136.7 (164.6). At `p = 2^31.25` they read 186.1
(213.9) for hiding, 132.9 (159.2) at the paper's bound, and 155.4 (183.6) for
encryption. Binding at the paper's bound was 128.2 (154.7) at the earlier
`SIGMA_C = 34641`: the narrower width shrinks that bound from `2^26.24` to
`2^25.91`.

Cost: the shuffle proof is 50.6 KB per message at `MSGS = 25` (1295902 bytes),
measured by the serialisation test, at 16.25 bits a Gaussian coefficient
against 18.38 for a flat encoding at the bound and 16.17 of entropy. It was
49.4 KB at `2^37` and 46.3 KB at the 2^31.25 modulus: the Gaussian
coefficients do not depend on `p`, and only the uniform ones widen, from 32 to
37 to 39 bits. Before that it was
46.4 KB (1188970 bytes, at 16.25 bits a Gaussian coefficient) when it still
carried the unused last `s_i`, against 47.1 KB at the previous widths; a
ciphertext grows by one ring element over the published parameters, and its
coefficients, which live modulo `q`, widen with `q` from 56 to 60 bits.
Prover time is unchanged: the linear proof expects `M^4 = 9.38` repetitions
whatever the widths are.

## Raising p (2026-10-10)

*This section and the next are history. Their "conservative bound" is
`2 * ||delta||_1 * 4 * sqrt(N) * SIGMA_C`, which leaves out the `sqrt(k)` of
constraint 2; "Raising p again" below redoes the choice at the correct bound.*

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
| `p = 2^37` | **134.3** | uniform coefficients 32 -> 37 bits, `q` to `2^60.79` |
| `HEIGHT = 2`, `WIDTH = 5` | 155.3 | about +14 KB a message, and the code assumes `HEIGHT = 1` |

Raising `p` is the cheapest: `SIGMA_C` and `SIGMA_S` do not depend on it, so
the masked openings, which are most of the proof, do not grow. The new `p` is
137438953349, `2^37 - 123`, the largest prime below `2^37` that is 5 mod 8, so
that `p - 1` fits in 37 bits; the first choice, `2^37 + 29`, took 38 and cost
0.6 KB a message for no measurable security. Its square roots of -1 are
`P0 = 6781345549` and `P1 = 130657607800`. `q` follows it through constraint 4:
1986862199687640941, about `2^60.785`, against the `2^60.585` required, also 5
mod 8, with `Q0 = 403703418285289875` and `Q1 = 1583158781402351066`, and far
above `24 * SIGMA_E^2 = 2^36.03`. The challenge space is unchanged at
`2^158.7`. Coefficients now hash in 5 bytes, the fewest that hold `p - 1`.

Side effects: hiding falls from 186.1 to 151.9 and encryption from 155.4 to
136.7, both still above binding at the conservative bound; the Lemma 1 margin
of constraint 3 grows from 8.3 to 446, and the `SIGMA_S` window of
SOUNDNESS.md, Section 8, widens from a factor of 4.1 to 30.1; the proof grows
from 46.3 to 49.4 KB a message.

## Lowering q (2026-10-10)

`q` had been sized against the paper's decryption condition, whose `N^2`
asked for `2^60.585`. Bounding what the code's decryption actually checks, as
constraint 4 now does, needs `2^56.91` in the worst case, with no failure
probability to account for. The new `q` is 288230376151711717, `2^58 - 27`,
the largest prime below `2^58` that is 5 mod 8, with `Q0 = 142220625783425069`
and `Q1 = 146009750368286648`; `x^512 +- Q0` are irreducible and multiply to
`x^1024 + 1`. It clears the bound by a factor of 2.1, and its coefficients take
58 bits. Going to `2^57 - 195` would add 3.1 bits of encryption for a margin of
6 per cent; one bit of margin seemed the better trade.

Encryption rises from 136.7 to **145.4** core-SVP (164.6 to 172.5 MATZOV), so
the level at the paper's binding bound rises from 136.7 to 145.4, still set by
encryption, with hiding at 151.9 above it. At the conservative bound it stays
134.3, set by binding. A ciphertext coefficient shrinks from 61 to 58 bits. The
shuffle proof is unaffected, at 49.4 KB a message, since it lives modulo `p`.

What gets worse is a malicious voter's lot. A ballot whose encryption
randomness is uniform in `[-R, R]` stops decrypting at a smaller `R`, since
decryption's ceiling `q / 4` fell, while the verifier's tail bound, which
decides where proofs stop verifying, does not depend on `q`; the window of
ballots that verify but do not decrypt therefore widens downwards. `master`'s
README has the measurement.

## Raising p again (2026-10-10)

Settling which binding bound extraction needs (constraint 2) put binding at
121.8 bits at `p = 2^37`: Section 5's vector has `k` elements, and the earlier
"conservative" bound had counted one. The levers, at the correct bound:

| change | binding | hiding | encryption | cost |
| --- | --- | --- | --- | --- |
| `p = 2^37`, `q = 2^58` | 121.8 | 151.9 | 145.4 | |
| `p = 2^38`, `q = 2^59` | 126.4 | 147.2 | 142.2 | |
| `p = 2^39`, `q = 2^60` | **131.1** | 142.2 | 139.0 | uniform coefficients 37 -> 39 bits |
| `p = 2^40`, `q = 2^61` | 136.1 | 137.6 | 135.9 | |
| `p = 2^37`, challenges with operator norm at most 12.47 | 130.2 | 151.9 | 145.4 | a challenge sampler that rejects by a floating-point FFT |

The last row bounds the stretch of a challenge by its largest singular value
rather than its `l1` norm, as LNP22 does, keeping the 90 per cent of `C_lin`
whose singular values are at most 12.47 (`|C_lin|` loses 0.15 bits); it was
passed over for a parameter change that needs no new code.

The new `p` is 549755813869, `2^39 - 19`, the largest prime below `2^39` that
is 5 mod 8, with `P0 = 248826153816` and `P1 = 300929660053`. `q` follows it
through constraint 4, which asks for `2^58.91`: 1152921504606846869,
`2^60 - 107`, the largest prime below `2^60` that is 5 mod 8, a factor of 2.1
above it, with `Q0 = 383824447308448458` and `Q1 = 769097057298398411`. For
both, `x^512 +- P0` are irreducible and multiply to `x^1024 + 1`.
Coefficients still hash in 5 bytes. The proof grows from 49.4 to 50.6 KB a
message. The Lemma 1 margin of constraint 3 grows from 446 to 1783, and the
`SIGMA_S` window of SOUNDNESS.md, Section 8, from a factor of 30.1 to 60.3.

## Caveats

- The decryption bound in constraint 4 is worst case. A statistical bound
  would be smaller still and would let `q` drop further, at the price of a
  decryption-failure probability that would have to be analysed and accepted.
- These are estimates against the best currently known attacks, in a model that
  counts one SVP call. They are not proofs.
