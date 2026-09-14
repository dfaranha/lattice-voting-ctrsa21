# Stage 0 of the LNP plan: what was built and what it found

[ISBIN-PLAN.md](ISBIN-PLAN.md) section 6 says not to start at B1 with the
intent of finishing B6, but to run a spike first: B1, B2 against a toy
statement, and enough of B7 to see whether the parameters close. This is the
report of that spike, and of the constant-coefficient proof that it identified
as the one piece with real design content, which was then built as well.

`is_bin` now holds end to end for a committed witness, on the hypothesis that
the witness is short. Supplying that hypothesis is B5, which is not built;
section 4 says what else is missing.

## 1. What is implemented

**B1, automorphisms** (`lnp_auto`, `lnp_auto_crt`). `sigma_k : X -> X^k` in
coefficient representation and in CRT representation. The feasibility note
recorded that `sigma_{-1}` exchanges the two CRT factors rather than fixing
them; that is now a test rather than a claim, and the general rule turns out to
be that `sigma_k` fixes the components when `k = 1 mod 4` and exchanges them
when `k = 3 mod 4`. The CRT version is tested against the coefficient version
for `k = 1, 3, 5, 7`, and the coefficient version is tested against the
definition of `sigma_{-1}`, against the homomorphism property, against being an
involution, and against the identity that the rest depends on: the constant
coefficient of `sigma_{-1}(a) * b` is the inner product `<a, b>`.

**B2, quadratic proof with garbage terms** (`lnp_quad_prover`,
`lnp_quad_verifier`). Proves `m0 * m1 = m2` for committed `m0, m1, m2` over a
BDLOP commitment with four message slots. Tested for completeness, for
rejecting a product that is off by one, and for rejecting a tampered opening.

**The `is_bin` product relation** (`lnp_isbin_prover`, `lnp_isbin_verifier`).
Proves that a committed `f` equals `sigma_{-1}(s) * (s - ones)` for a committed
`s`, where `ones` is `1 + X + ... + X^(n-1)`. The constant coefficient of that
product is exactly the sum over `j` of `s_j (s_j - 1)`, which is the binary
constraint. This is the algebraic heart of `is_bin`, and it is where the
CRT-aware automorphism earns its place: `sigma` is applied to `u_0`, to the
challenge, and to `u_3`, and each of those exchanges the two CRT components.

The mixed challenge terms are what make this more than a relabelling of B2.
Because `sigma(u_0) = sigma(v_0) - sigma(c) sigma(s)`, the expansion carries
both `c` and `sigma(c)`, and the two garbage slots have to absorb them
separately: slot 2 commits the coefficient of `c`, and slot 3 commits *sigma
of* the coefficient of `sigma(c)`, so that applying `sigma` to its `u` recovers
a term multiplied by `sigma(c)` rather than by `c`.

**B3, the constant coefficient** (`lnp_ct_prover`, `lnp_ct_verifier`). Proves
that the constant coefficient of a committed value is zero, which is what turns
the product relation above into `is_bin`. For each of `LNP_LAMBDA` masks `g_i`,
sampled uniformly subject to `ct(g_i) = 0`, the prover publishes
`h_i = g_i + mu_i * f` for a scalar `mu_i`, and the verifier checks
`ct(h_i) = 0`. Tying `h_i` back to the committed values is linear, so it needs
no garbage terms: with `B_i = b2[SLOT_G + i] + mu_i * b2[SLOT_F]`, the value
`T_i - h_i` is a commitment to zero under `B_i`.

Two orderings carry the whole argument. The scalars are derived from the
commitment alone, so the masks are fixed before them; and `h` is absorbed into
the opening challenge, so it is fixed before that. Both are tested: one test
plays a prover that reads off the scalars it would face and then rewrites its
masks to cancel a non-zero constant coefficient, and fails, because committing
the rewritten masks changes the commitment and so changes the scalars.

**The two halves composed.** `is_bin` is the product proof and the
constant-coefficient proof over one commitment, the first filling the garbage
slots and the second reading slot 1. Tested to accept a binary witness and to
reject one carrying a single coefficient of 2 or of -1.

## 2. B7: the parameters close

The plan expected pressure towards a larger modulus. For the machinery built
here there is none. At `DEGREE = 1024` and `MODP` about `2^31.86`:

| constraint | value | comfortable? |
| --- | --- | --- |
| no wraparound in `sum_j s_j (s_j - 1) < p` | `\|\|s\|\|_inf <= 1953` | yes, against 1 for honest `s` |
| MLWE hiding, `LNP_RANK = 2` | 187 bits core-SVP | yes |
| MSIS binding, bound `16 sigma_C sqrt(nu N)` | 118 bits core-SVP | same as the base scheme |
| `\|\|c - c'\|\| < sqrt(p/2)` for Lemma 1 | 8.5 against 44195 | yes |
| soundness of the scalar aggregation, `p^-LNP_LAMBDA` | `2^-127.5` at `LNP_LAMBDA = 4` | just |

MSIS binding is flat in the commitment width, so widening the randomness to
make room for the message slots costs nothing there. Hiding is the constraint
that dictates the width: rank 1 is 73 bits and rank 2 is 187, which is why
`LNP_WIDTH` is `HEIGHT + SLOTS + 2` rather than `HEIGHT + SLOTS + 1`.

The aggregation is the one place where the modulus is close to binding. Each
mask contributes a factor `1/MODP`, so four give `2^-127.5` here, which clears
128 bits only just. On the smaller modulus proposed on the `balanced-params`
branch, `2^31.25`, four masks give `2^-125.0` and a fifth would be needed. A
mask costs one commitment slot and one published ring element, so this is
cheap to fix but is worth noting: it is a constraint the base scheme does not
have, and it couples `LNP_LAMBDA` to any future change in the modulus.

So the verdict the spike was meant to produce is: **the quadratic layer closes
at the existing modulus.** The modulus pressure the plan anticipated belongs to
B5, the approximate range proof, which is not implemented and not assessed
here.

## 3. The obstacle the spike found, and how it was removed

The spike stopped at `ct(f) = 0` for a committed `f`, and the reason is worth
keeping, because it dictates the shape of the solution.

The obvious construction does not work. Commit a mask `g` with `ct(g) = 0`,
send `h = g + f` in the clear, and have the verifier check `ct(h) = 0`. This is
unsound: nothing forces `ct(g) = 0`, so a prover with `ct(f) != 0` simply
commits `g` with `ct(g) = -ct(f)`. The quadratic relation `h = g + f` still
holds and the constant coefficient still vanishes.

The repair has to stop the prover from knowing what to cancel, which means the
mask must be committed before a challenge, and the challenge must act on the
constant coefficient linearly. That last requirement is what forces the
challenge to be a **scalar**: for a ring challenge `gamma`, `ct(gamma * f)` is
not `gamma * ct(f)`. It is a linear form in the coefficients of `gamma` applied
to `f`, so requiring it to vanish for a random ring `gamma` would force
`f = 0`, which is false for an honest binary `s`, whose `f` is non-zero in
every coefficient but the constant one. Multiplying by a ring challenge
destroys exactly the structure the statement is about.

With scalars, `ct(h_i) = ct(g_i) + mu_i * ct(f)` really is linear, a prover
facing `ct(f) != 0` needs `ct(g_i) = -mu_i * ct(f)` for a scalar it cannot
predict, and each mask independently catches it with probability `1 - 1/p`.
This is also why the construction wants to be used on many claims at once: the
`LNP_LAMBDA` masks are paid for once regardless of how many constant
coefficients are being proven zero, so proving `is_bin` for all `MSGS`
permutation elements together costs the same four masks as proving it for one.

## 4. B5, the approximate range proof, and the modulus it forced

`is_bin` needs the witness to be short, and B5 is what supplies that. Costing
it first, before writing it, is what set the modulus.

The wraparound condition allows `||s|| <= sqrt(p/2)`, against about 22.6 for an
honest binary witness. That looks like enormous slack and is not, once the
proof has to hide `s`. Taking the bound from an ABDLOP opening certifies about
`2^22` and would need `p > 2^53.5`. The JL projection does better, but its mask
is driven by `||R s||`, so the certified bound grows linearly in the number of
projected coordinates:

| PROJ | tau | certifies | needs |
| --- | --- | --- | --- |
| 64 | 6 | 69511 | `p > 2^33.2` |
| 256 | 6 | 278046 | `p > 2^37.2` |
| 256 | 9 | 417069 | `p > 2^38.4` |
| 256 | 12 | 556091 | `p > 2^39.2` |

`PROJ = 256` is what gives the projection lemma its `2^-128`. At the old
modulus of `2^31.86` none of these close: B5 would have certified 278046
against a relaxed norm proof already achieving 32768, so it would have been
about eight times weaker than the proof it replaces.

`TAU_PROJ = 9` is the chosen row. It costs about 3.8 rejection-sampling
repetitions and leaves the certified bound a factor 1.78 inside the ceiling.

### What the modulus change dragged with it

`p` at `2^40` is not a local change.

* **Hiding collapses at `WIDTH = 3`.** The commitment's hiding is MLWE of rank
  `k - n - 1`, and rank 1 at `2^40` is 49.6 bits. `WIDTH = 4` is therefore not
  optional here; it is what the modulus forces.
* **`q` follows `p` linearly** through `q > 2p(2 * DIM * N^2 + N + 1)`, and `q`
  has to stay inside 64 bits. That ceiling is what fixes `p`: at `2^41` even
  `DIM = 2` needs `2^64`. At `p = 2^40` the largest workable encryption rank is
  `DIM = 3`, with `q` saturated just under `2^64`.
* **`DIM = 3` needed the vericrypt indexing fix**, since three places walked
  CRT-indexed arrays with `DIM`. Those were harmless only while `DIM` and
  `NCRT` were both 2.

### Final parameters

| | was | now |
| --- | --- | --- |
| `MODP` | 3906450253 (`2^31.86`) | 1099511627917 (`2^40`) |
| `WIDTH` | 3 | 4 |
| `DIM` | 2 | 3 |
| `Q` | `2^56` | 18446744073709551557 (`2^64 - 59`) |
| `SIGMA_P` | -- | 3258 |

| problem | core-SVP | MATZOV |
| --- | --- | --- |
| MLWE commitment hiding | 137.6 | 164.7 |
| MSIS commitment binding | 162.1 | 187.3 |
| MLWE encryption | 129.1 | 157.3 |
| scalar aggregation | `2^-160` | |
| JL certified bound against the ceiling | 1.78x | |

Encryption at 129.1 is the thinnest, and is what saturating `q` under 64 bits
costs.

### Two latent bugs the modulus change exposed

Neither was reachable while `MODP` stayed under `2^32`, and both are now fixed.

* `commit_sample_short` held `MODP - 1 + d` in a `uint32_t`. At `2^40` that
  truncates, the randomness stops being ternary and becomes coefficients near
  141, and rejection sampling then never accepts: the prover hangs rather than
  failing. This is the kind of bug that only a parameter change finds.
* `commit_norm2_leq` bailed out when a coefficient exceeded the bound, but the
  bound is a *squared* norm, so a coefficient below it could still square past
  `2^64` and wrap. A malicious response with coefficients around `2^32` would
  have passed the norm check. It now compares by dividing instead.

## 5. What this changes about the decision

The plan's section 7 framed the choice as fidelity to Protocol 1 against
1500-2000 lines of new proof-system code. Three things have moved.

The parameter risk is gone: the modulus does not have to change, so the
cascade the plan worried about, into the CRT constants and the ciphertext
sizes, does not happen. The one new coupling is `LNP_LAMBDA` against the
modulus, and it is cheap.

The cost estimate came down. B1, B2 and B3 together are about 900 lines, well
under the plan's figure for the whole of Track B, and the piece that looked
hardest turned out to have a short answer once the reason it was hard was
clear.

B1, B2, B3 and B5 are built. What is left is B4, which the JL route does not
need, since the range proof bounds a BDLOP-committed witness directly, and B6,
the wiring into the shuffle. Until B6, nothing here replaces the relaxed norm
proof in `shuffle.c`; the two are independent.

The cost is no longer the 1500-2000 lines the plan estimated. It is the
modulus, and everything the modulus drags with it: `WIDTH`, `DIM`, `q`, the CRT
constants, and the two 32-bit assumptions that had been sitting unreachable in
the arithmetic. Track A would have paid the same price, since LaZer's range
proof obeys the same arithmetic.

Finally, none of this is a soundness proof. The extraction argument for the
product relation divides by `(c - c')(sigma(c) - sigma(c'))` and yields a
relaxed opening, and what that relaxation does to the binding bound has not
been worked out. The projection lemma is used with a rough constant of 2; the
exact one would move the 1.78x margin by a little. The tests establish
completeness and that specific cheating strategies fail. They do not establish
that no strategy succeeds.
