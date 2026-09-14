# Parameter selection

This note records how the parameters in `param.h` were chosen, so that they can
be re-derived rather than copied. All figures are classical core-SVP
(`0.292 * beta`, the ADPS16 model) unless stated otherwise; the MATZOV
gate-count figures are given alongside because the two models disagree by
roughly 25 bits here and the published claim depends on which one is used.

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
| MLWE, encryption | ballot secrecy | 89 | 117 |

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
   the reduction gives, which is about 14 times larger than the `2 * sqrt(k) *
   2 * sqrt(N) * SIGMA_C` that the verifier's own norm check suggests; using
   the latter overstates binding security by some 60 bits.
3. **Invertibility (Lemma 1)**: the extracted `||(d - d') * sigma_i||` must
   stay below `sqrt(p / 2)`, which floors `p` at `32 * N * SIGMA_S^2`.
4. **Decryption correctness**: `q > 2p(2 * DIM * N^2 * BETA^2 + N + 1)`
   (paper, Table 1). The factor is `2^22`, so `q` scales linearly with `p`.
5. **Quasi-unique responses**: `q > 24 * SIGMA_E^2 = 2^36`. Never binding.
6. **Challenge space**: Theorem 1 has a `4 * tau * t` term, so `|C|` should
   exceed about `2^135`. The sampler in `shuffle.c` draws binary challenges, so
   `|C| = C(N, nu)`, which is `2^139` at `nu = 20`.
7. Both `p` and `q` must be prime and congruent to 5 mod 8, so that
   `x^N + 1` splits into exactly `NCRT = 2` factors.

The interesting interaction is that constraints 1 and 2 pull `p` in opposite
directions, and 3 and 4 pin it from below and tie `q` to it. Lowering `p`
improves hiding and shrinks the proof, but weakens binding; since binding is
the tight constraint after `WIDTH` grows, `p` stays roughly where it was.

## What changed

| | published | here |
| --- | --- | --- |
| `WIDTH` | 3 | 4 |
| `NONZERO` | 36 | 20 |
| `MODP` | 3906450253 (2^31.86) | 2553802997 (2^31.25) |
| `SIGMA_C` | 54000 | 34641 |
| `SIGMA_S` | 256 | 190 |
| `DIM` | 2 | 3 (now derived as `WIDTH - HEIGHT`) |
| `VECTOR` | 3 | 4 (now derived as `WIDTH`) |
| `Q` | 2^56 | 36028797018964429 (2^55) |

`WIDTH = 4` is what takes the MLWE rank from 1 to 2. There is nothing in
between: rank 1 gives 72 bits and rank 2 gives 186, so hiding ends up
over-provisioned, and the only way to trade that surplus back would be a
smaller degree with a higher rank, which cuts against raising `DEGREE` later.

`NONZERO = 20` pays for most of `WIDTH = 4`. The MSIS bound grows as
`nu^1.5` (once through `SIGMA_C`, which is linear in `nu`, and again through
the explicit `sqrt(nu)`), so cutting `nu` from 36 improves binding and shrinks
the responses at the same time. At 36 the challenge space was `2^222`, far more
than soundness needs.

`SIGMA_S = 190` rather than the size-optimal 76: the Lemma 1 floor at this `p`
is `2^30.14`, well below `p`, so the headroom is free, and the wider Gaussian
keeps the rejection-sampling cost at today's 3.07 expected repetitions instead
of 5.6.

`SIGMA_E` is left at 54000. Table 1 would allow 28160 at these dimensions, but
the decryption bound in `vericrypt.c` has not been re-derived for a smaller
value, and 54000 is the conservative direction.

## Result

| Problem | published | here | |
| --- | --- | --- | --- |
| | core-SVP | core-SVP | MATZOV |
| MLWE, commitment hiding | 72 | **186** | 210 |
| MSIS, commitment binding | 117 | **128** | 155 |
| MLWE, encryption | 89 | **157** | 184 |

Cost: the shuffle proof grows from 65.5 KB to 71.7 KB per message, about 9 %,
and a ciphertext grows by one ring element. Prover time is unchanged, since the
rejection-sampling rate is held at 3.07.

## Caveats

- The decryption bound in constraint 4 uses the worst-case ring product bound
  `||ab||_inf <= N ||a||_inf ||b||_inf` applied twice, hence the `N^2`. A
  statistical bound would be much smaller and would let `q` drop further, which
  would shrink ciphertexts and raise encryption security at the same time. That
  needs a failure-probability analysis, which has not been done here.
- These are estimates against the best currently known attacks, in a model that
  counts one SVP call. They are not proofs.
