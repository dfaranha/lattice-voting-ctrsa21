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

## 4. B5, and why it does not pay at this modulus

`is_bin` as built holds for a committed witness **given that the witness is
short**. Supplying that hypothesis is B5, the approximate range proof. Working
out what B5 could actually certify here is what the rest of this section is,
and the answer is that it should not be built at the current modulus.

The argument needs `sum_j s_j (s_j - 1)` not to wrap, which in the 2-norm means
`||s|| <= sqrt(p/2) = 44195`, against about 22.6 for an honest binary witness.
That looks like enormous slack. It is not, once the proof has to hide `s`.

Two routes were costed.

**The opening bound on its own.** Put the witness in the Ajtai part of an
ABDLOP commitment, which is B4, and take the bound from the verifier's norm
check on the response. The mask has to hide `c * s`, so its width is
`tau * ||c||_1 * ||s||`, and extraction doubles the verifier's bound twice
over. That certifies about `2^22`, and needs `p > 2^53.5`.

**The JL projection, which is B5 proper.** Project onto `PROJ` coordinates with
a public sign matrix, publish the projection masked, and bound `||s||` from the
bound on the projection. The mask width is now driven by `||R s||`, which grows
as `sqrt(PROJ)`, and the published vector has `PROJ` entries, so the certified
bound grows linearly in `PROJ`:

| PROJ | tau | certifies `\|\|s\|\|` | needs |
| --- | --- | --- | --- |
| 64 | 6 | 69511 | `p > 2^33.2` |
| 128 | 6 | 139023 | `p > 2^35.2` |
| 256 | 6 | 278046 | `p > 2^37.2` |
| 256 | 12 | 556091 | `p > 2^39.2` |

`PROJ = 256` is what gives the projection lemma its `2^-128`, and `tau = 6`
already costs about seven rejection-sampling repetitions, so the honest target
is the third row: **`p` of about `2^37.2`**, against `2^31.86` now.

The comparison that settles it is not against the ceiling but against what the
branch already has. The relaxed norm proof in `shuffle.c` bounds the extracted
`||(d - d') sigma||` by `2 * 2 sqrt(n) SIGMA_S = 32768`, comfortably inside the
Lemma 1 ceiling of 44195. At this modulus B5 would certify 278046. **The range
proof would be about eight times weaker than the proof it is meant to
replace**, so building it here would not close the gap; it would widen it.

The factor of 2 used for the projection lemma is rough, and the exact constant
from LNP would move these numbers by a bit or so. It would not move them by the
three to four bits that separate `2^33.2` from `2^37.2`, nor by the six that
separate the certified bound from the ceiling.

So B5 is not blocked on code. It is blocked on a modulus of roughly `2^37` to
`2^39`. At `2^40` the other constraints stay comfortable: hiding at rank 2 is
140 bits, MSIS improves as `p` grows, the scalar aggregation improves to
`2^-160`, and the Lemma 1 floor on `SIGMA_S` relaxes. What it costs is the
cascade the original plan predicted, into the CRT constants and, through
`q > 2p(2 l N^2 + N + 1)`, into the ciphertexts.

B4 and B6 remain unbuilt as well, but neither is on the critical path until the
modulus question is settled: B4 only matters for the route that was costed at
`p > 2^53.5`, and B6 has nothing to wire in until B5 exists.

Finally, none of this is a soundness proof. The extraction argument for the
product relation divides by `(c - c')(sigma(c) - sigma(c'))` and yields a
relaxed opening, and what that relaxation does to the binding bound has not
been worked out. The tests establish completeness and that specific cheating
strategies fail. They do not establish that no strategy succeeds.

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

What is left is B4 and B5, and section 4 turns that from a coding question into
a parameter one. B5 needs `p` of about `2^37` to `2^39`; below that it
certifies a weaker bound than the relaxed proof it would replace. So the Track
A against Track B decision is now downstream of a modulus decision, and the
same modulus decision would be forced on Track A, since LaZer's range proof
obeys the same arithmetic.

That also means the honest comparison is no longer is_bin against nothing. It
is is_bin at a larger modulus, with the cascade that implies, against the
relaxed norm proof that already closes the gap at the modulus in use. The
first buys fidelity to Protocol 1 as published and removes the need to review
the restatement in SOUNDNESS.md. Whether that is worth a modulus change is a
judgement about the artifact, not about the cryptography.
