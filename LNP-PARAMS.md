# Stage 0 of the LNP plan: what was built and what it found

[ISBIN-PLAN.md](ISBIN-PLAN.md) section 6 says not to start at B1 with the
intent of finishing B6, but to run a spike first: B1, B2 against a toy
statement, and enough of B7 to see whether the parameters close. This is the
report of that spike. It is not a complete `is_bin` proof, and section 3 says
exactly what is missing.

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

## 2. B7: the parameters close

The plan expected pressure towards a larger modulus. For the machinery built
here there is none. At `DEGREE = 1024` and `MODP` about `2^31.86`:

| constraint | value | comfortable? |
| --- | --- | --- |
| no wraparound in `sum_j s_j (s_j - 1) < p` | `\|\|s\|\|_inf <= 1953` | yes, against 1 for honest `s` |
| MLWE hiding, `LNP_RANK = 2` | 187 bits core-SVP | yes |
| MSIS binding, bound `16 sigma_C sqrt(nu N)` | 118 bits core-SVP | same as the base scheme |
| `\|\|c - c'\|\| < sqrt(p/2)` for Lemma 1 | 8.5 against 44195 | yes |

MSIS binding is flat in the commitment width, so widening the randomness from
3 to 7 to make room for the message slots costs nothing there. Hiding is the
constraint that dictates the width: rank 1 is 73 bits and rank 2 is 187, which
is why `LNP_WIDTH` is `HEIGHT + SLOTS + 2` rather than `HEIGHT + SLOTS + 1`.

So the verdict the spike was meant to produce is: **the quadratic layer closes
at the existing modulus.** The modulus pressure the plan anticipated belongs to
B5, the approximate range proof, which is not implemented and not assessed
here.

## 3. What is missing, and the one real obstacle

`is_bin` needs two things: the product relation, which is implemented, and
`ct(f) = 0` for the committed `f`, which is not. Proving that the constant
coefficient of a *committed* polynomial is zero turned out to be the hard part,
and the spike pinned down why.

The obvious construction does not work. Commit a mask `g` with `ct(g) = 0`,
send `h = g + f` in the clear, and have the verifier check `ct(h) = 0`. This is
unsound: nothing forces `ct(g) = 0`, so a prover with `ct(f) != 0` simply
commits `g` with `ct(g) = -ct(f)`. The quadratic relation `h = g + f` still
holds and the constant coefficient still vanishes.

The repair is to make the mask commit before the challenge and to aggregate
with a challenge that acts on the constant coefficient linearly, so that the
prover cannot precompute the offset. That forces a **scalar** challenge, and
this is the obstacle: for a ring challenge `gamma`, `ct(gamma * f)` is not
`gamma * ct(f)`. It is a linear form in the coefficients of `gamma` applied to
`f`, so requiring it to vanish for a random ring `gamma` forces `f = 0`, which
is false for an honest binary `s`, whose `f` is non-zero in every coefficient
but the constant one. Multiplying by a ring challenge destroys exactly the
structure the statement is about.

So the missing piece is not more of the same algebra. It is a genuine
subroutine: either scalar-challenge aggregation of many constant-coefficient
claims at once, which is what LNP does and which only pays off when there are
many claims to batch, or a direct proof that a committed polynomial has zero
constant coefficient via the automorphism trace, which needs the whole
automorphism group rather than `sigma_{-1}` alone.

Also absent, as in the original plan: B4 (the Ajtai part of ABDLOP, so that the
committed `s` carries a norm bound at all), B5 (the approximate range proof
that would supply `\|\|s\|\|_inf <= 1953`), and B6 (wiring into the shuffle).
Without B4 and B5 the norm hypothesis is assumed, not proven, so nothing here
yet replaces the relaxed norm proof in `shuffle.c`.

## 4. What this changes about the decision

The plan's section 7 framed the choice as fidelity to Protocol 1 against
1500-2000 lines of new proof-system code. The spike moves two things.

It removes the parameter risk: the modulus does not have to change for the
quadratic layer, so the cascade the plan worried about, into the CRT constants
and the ciphertext sizes, does not happen at this stage.

It sharpens the cost. The automorphism layer and the quadratic layer came to
about 600 lines and behave. The constant-coefficient subroutine is the piece
with real design content, and it is the piece LaZer already has. If the
artifact is going to take Track A, this is the natural place to stop.
