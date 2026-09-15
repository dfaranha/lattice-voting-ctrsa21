# The LNP machinery in this repository

What `lnp.c` builds, why each piece is shaped the way it is, and how the
pieces combine into the `is_bin` argument that the shuffle needs. This is the
approach document. `ISBIN-PLAN.md` is the plan it was built from,
`LNP-PARAMS.md` records the parameters and the measurements, and
`SOUNDNESS.md` treats the relaxed Lemma 5 that sits above all of it.

The style is Lyubashevsky-Nguyen-Plancon: commit once to everything, prove
relations about the committed values with automorphisms and garbage terms,
and aggregate many claims into a few by random scalars.

## 1. The problem it solves

Protocol 1 of ePrint 2025/658 patches Neff's proof of shuffle over a ring that
is not a field. The patch is Lemma 5, and Lemma 5 needs each committed
permutation element `sigma_i` to lie in a set `D` on which the product
argument is sound.

`D` here is the set of **binary** ring elements: those whose coefficients are
all 0 or 1. A norm bound is not enough. `R_p` splits into two CRT components,
and an element can be short while its two components encode different
permutations; that is exactly the attack in Section 4.1 of the paper, and
`shuffle.c` has a test that mounts it. So membership in `D` has to be proven
coefficient-wise, not by a norm check, and `is_bin` is that proof.

## 2. The ring, and the one identity everything rests on

`R_p = Z_p[X]/(X^n + 1)` with `n = DEGREE = 1024` and `p = MODP` chosen so that
`p = 5 mod 8`. That congruence makes `X^n + 1` split into exactly two
irreducible factors, `X^(n/2) + P0` and `X^(n/2) - P0`. Everything in the code
carries elements in this two-component CRT representation, `pcrt_poly_t`,
because multiplication is componentwise there.

The automorphisms are `sigma_k : X -> X^k` for odd `k`. The one that matters is
`sigma_{-1}`, written in the code as `sigma_k` with `k = 2n - 1`. Two facts
about it drive the whole construction:

**It exchanges the CRT components.** `X^(n/2)` squares to `-1`, so raising it
to an odd power `k` gives `X^(n/2)` when `k = 1 mod 4` and `-X^(n/2)` when
`k = 3 mod 4`; in the second case one factor is sent onto the other. So
`lnp_auto_crt` is *not* componentwise, and `lnp_auto_swaps` is what decides.
Getting this wrong is silent, so `lnp_auto_crt` is tested against a
coefficient-representation reference, `lnp_auto`.

**It turns products into inner products.** For any `a`, `b`,

```
ct( sigma_{-1}(a) * b ) = <a, b>
```

where `ct` is the constant coefficient and `<.,.>` the inner product of
coefficient vectors. This is the hinge. It converts every statement about
coefficient vectors into a statement about one coefficient of one ring
element, which is a thing a commitment scheme can prove.

## 3. The commitment

A BDLOP commitment in Hermite normal form, generalised from `commit.h` to carry
several message slots under one randomness vector.

```
lnpkey_t:   B1[HEIGHT][LNP_WIDTH]      Ajtai part
            b2[SLOTS][LNP_WIDTH]       one row per message slot

lnpcom_t:   c1[i] = <B1[i], r>
            c2[i] = <b2[i], r> + m_i
```

`LNP_WIDTH = HEIGHT + SLOTS + LNP_RANK`. The Hermite form means `b2[i]` is the
unit vector `e_{HEIGHT+i}` plus a random tail on the last `LNP_RANK`
coordinates, so each slot needs its own randomness component on top of the
Ajtai part. Binding is MSIS on `[B1; b2]`; hiding is MLWE of rank `LNP_RANK`
in the trailing components, which is why `LNP_RANK` is 2: rank 1 does not
reach the target. `LNP-PARAMS.md` has the figures.

The slots:

| slot | holds |
| --- | --- |
| `SLOT_S` = 0 | the witness `s`, claimed binary |
| `SLOT_F` = 1 | `f`, the product whose constant coefficient is claimed zero |
| `SLOT_W` = 2 | the projection mask, packed into the first `PROJ` coefficients |

The garbage terms of the quadratic proof used to occupy two more slots here.
They are batched now, and live in a commitment of their own; section 6 has
the reason and the shape.

An opening is a vector `z = y + d r` for a challenge `d` and mask `y`, and the
verifier checks `<B1[i], z> = w_i + d c1[i]` against the Ajtai first message
`w_i`, plus a norm bound on `z`. Those two together are the binding.

## 4. The three claims

`s` is binary if and only if `s_j (s_j - 1) = 0` for every `j`. Over the
integers that is equivalent to `sum_j s_j (s_j - 1) = 0`, since every term is a
product of two integers differing by one and is therefore non-negative. The
sum is a constant coefficient by the identity of section 2:

```
f = sigma_{-1}(s) * (s - 1)        ct(f) = sum_j s_j (s_j - 1)
```

which is `lnp_isbin_product`.

That equivalence holds over the integers, and the coefficients live in `Z_p`,
so it needs `sum_j s_j (s_j - 1)` not to wrap around `p`. **This is why all
three claims below are needed and why none can be dropped.** The
constant-coefficient claim says the sum is zero mod `p`; only the shortness
claim makes that imply it is zero over the integers; and only the product
claim says the committed `f` is the thing whose constant coefficient is being
read. Each one is load-bearing for the next.

### 4.1 The product relation

That `slot SLOT_F` really holds `sigma_{-1}(s) * (s - 1)` for the `s` in
`SLOT_S`. This is quadratic in the committed values, which is what needs the
garbage terms.

Write `v_i = <b2[i], y>` for the mask's contribution to slot `i`, and let
`u_i = <b2[i], z> - d c2[i]`. Expanding an honest opening gives `u_i = v_i - d
m_i`: mask minus challenge times message. Substituting that into the quadratic
form produces terms of degree 0, 1 and 2 in the challenge `d`. The degree-2
terms cancel on their own because the relation holds. The degree-1 terms do
not, and that is what slots 2 and 3 are for: the prover puts

```
g1 = -sigma(v_0) * (s - 1)
g2 = sigma( v_1 - v_0 * sigma(s) )
```

chosen so that their contribution cancels the challenge-linear part exactly.
What survives is a challenge-free term

```
t = v_2 + sigma(v_3) + sigma(v_0) * v_0
```

against which the verifier checks the whole expansion recomputed from
`u_0 .. u_3`. Section 6 describes how `g1`, `g2` and `t` are shared across the
batch rather than paid for per message; the shape of the identity is the same
either way.

Note that the garbage terms depend on the mask, so they are part of the proof
rather than of the statement: they are computed when the mask is drawn, after
the rest of the commitment is fixed.

### 4.2 The constant coefficient is zero

That `ct(f) = 0`. The trick is aggregation by **scalars**.

For a scalar `nu` and a ring element `f`, `ct(nu * f) = nu * ct(f)`. This is
false for a general ring multiplier, which is the whole reason the construction
looks the way it does: `ct(gamma * f)` is not `gamma * ct(f)` for ring `gamma`,
so only scalar aggregation is available.

The prover samples `LNP_LAMBDA` masks `g_j` uniformly subject to `ct(g_j) = 0`,
commits them, and publishes

```
h_j = g_j + nu_j * f
```

The verifier checks `ct(h_j) = 0`. If `ct(f) = 0` this holds for any `nu_j`. If
`ct(f) != 0` the prover needs `ct(g_j) = -nu_j ct(f)`, and the `nu_j` are not
known when the `g_j` are committed, so each mask independently catches a
cheating prover with probability `1 - 1/p`. Four masks at this modulus give
about `2^-160`.

Tying `h_j` back to the commitment is linear, so it costs no garbage terms: with
`B_j` the appropriate combination of key rows, `T_j - h_j` is a commitment to
zero under `B_j`, and the verifier checks that in the same opening.

Two orderings carry this and are the fragile part. The scalars must be derived
*after* the masks are committed, and `h_j` must be fixed *before* the challenge.
Both are enforced by what goes into which hash, and both are tested.

### 4.3 The witness is short

An approximate range proof by Johnson-Lindenstrauss projection. `PROJ = 256`
rows `r_i` with entries in `{-1, 0, 1}` are derived from a seed. The prover
publishes

```
z_i = <r_i, s> + w_i
```

where `w` is a Gaussian mask packed into the first `PROJ` coefficients of
`SLOT_W`, and the verifier checks that the published `z` is short. The
projection lemma then says a short projection certifies a short `s`, up to a
constant.

The published `z` has to be tied back to the committed `s` and `w`, which is
again done by scalar aggregation, and this is where the pieces join. With
scalars `mu_ji`, define

```
P_j = sigma( sum_i mu_ji r_i )        multiplier for the witness slot
M_j = sigma( sum_i mu_ji X^i )        multiplier for the packed mask slot
Z_j = sum_i mu_ji z_i                 a scalar
```

Then by the identity of section 2,

```
ct( P_j * s ) = sum_i mu_ji <r_i, s>
ct( M_j * W ) = sum_i mu_ji w_i
```

so `ct( P_j s + M_j W - Z_j ) = sum_i mu_ji ( <r_i,s> + w_i - z_i )`, which is
zero exactly when the published projection is the honest one. Building `P_j`
as the automorphism of a *row combination* rather than one automorphism per row
is what keeps this to `LNP_LAMBDA` automorphisms instead of `PROJ` of them,
and works because `sigma` is linear.

The mask `w` cannot be reused: the prover rejection-samples it, and a rejection
means recommitting, because the projection is fixed by the commitment. That is
the `for (tries...)` loop around `lnp_bin_setup` in the caller. `TAU_PROJ`
trades the tightness of the certified bound against the number of repetitions;
9 gives about 3.8 repetitions.

## 5. Merging the three

All three claims speak about one commitment under one randomness, so running
them as separate protocols meant masking the same randomness three times and
sending three openings. They are merged instead:

- one mask `y`, one challenge `d`, one opening `z`;
- the product relation contributes its challenge-free term `t`;
- the constant-coefficient claim and the projection-consistency claim are
  aggregated into the *same* values `h_j`, because both are "this constant
  coefficient is zero" claims and scalars compose:

```
h_j = g_j + nu_j * f + P_j * s + M_j * W - Z_j
```

This merge was not only an optimisation. Run separately with the *same* masks
`g_j`, the difference `h_j(ct) - h_j(range)` cancels `g_j` and reveals a linear
function of the witness. Merging removed that.

The linear proof of the shuffle opens this same commitment under this same
randomness, so it shares the opening too: `lin_first` in `shuffle.c` drives
`lnp_bin_first` with its own mask, and `lin_verifier` hands its opening to
`lnp_bin_check`.

## 6. Batching across messages

The shuffle proves `is_bin` for all `MSGS` permutation elements at once, and
the aggregation does not care whether the claims being aggregated belong to
one message or many. Four things follow.

**One challenge.** `batch_hash` derives a single challenge over every message's
first messages. This is forced, not merely convenient: the verifier can only
form the challenge-weighted sum of per-message terms, while a shared aggregated
value carries the unweighted sum, and that value must be fixed before the
challenge.

**One rejection test.** One challenge means one rejection test, over the
concatenation of every response; testing each message separately would multiply
the abort probabilities. The masked term is therefore `sqrt(MSGS)` times longer
than for a single message, so the masks widen from `SIGMA_C` to
`SIGMA_B = SIGMA_C sqrt(MSGS)`, and the verifier's norm bounds widen with them.

**One set of garbage terms.** The quadratic proof's two garbage terms are also
paid for once. Weighting message `l` by a batching challenge `rho_l` and
summing, the per-message identities close as one:

```
sum_l rho_l [ sigma(u_l0)(u_l0 + d) + sigma(d) u_l1 ] + U_2 + sigma(U_3) = T
```

where `U_2`, `U_3` open a single garbage commitment holding

```
G1 = sum_l rho_l       g1_l
G2 = sum_l sigma(rho_l) g2_l
T  = sum_l rho_l sigma(v_l0) v_l0  +  V_2 + sigma(V_3)
```

The `sigma(rho_l)` on the second accumulator is the one detail that is easy to
get wrong: the verifier's identity applies `sigma` to that row, and
`sigma(sigma(rho_l) g2_l)` is `rho_l sigma(g2_l)`, which is what the expansion
needs. Weighting both by `rho_l` does not close.

`rho` has to be drawn after every message's commitment, since the garbage
terms are weighted by it, and before those terms are committed, since a prover
that knew it could adapt them. That is a separate hash from the one producing
the opening challenge, and it is why the garbage terms cannot share the mask
commitment: that one is fixed earlier still, before the projection scalars are
derived.

Batching them takes `SLOTS` from 5 to 3 and `LNP_WIDTH` from 8 to 6, and drops
the per-message `t`, at the cost of one commitment for the batch.

**One set of masks.** The `LNP_LAMBDA` masks `g_j` are paid for once for the
whole batch, not once per message. They live in their own commitment,
`lnpmaskcom_t`, under their own key and randomness of width `MASK_WIDTH`, and
the aggregated values become

```
h_j = g_j + sum_l [ nu_lj f_l + P_lj s_l + M_lj W_l - Z_lj ]
```

with each message contributing its share in `lnp_bin_setup` (for `h`) and
`lnp_bin_first` (for `v`), and `lnp_batch_check` settling the comparison once
after `lnp_bin_check` has accumulated every message's side into `acc`.

Two consequences are worth stating plainly.

*The mask commitment is a commitment, and has to be opened like one.* The batch
rows are `LNP_LAMBDA` equations in `MASK_WIDTH` unknowns, so without an Ajtai
opening equation and a norm bound on `z_mask` a prover can choose any `h` with
zero constant coefficient and solve for `z_mask`. `lnp_batch_check` enforces
both.

*`MSGS` became a security parameter.* Since `SIGMA_B` scales with
`sqrt(MSGS)`, so does the extractable opening, and MSIS binding falls with the
batch size: 162 bits core-SVP at one message, 129 at 25, 100 at 1000.

There is no way around this by blocking, because `MSGS` is the anonymity set
and not merely a batch size: the proof says the output list is a permutation
of the input list, so shuffling in blocks would prove only that each output
block permutes its own input block, and would reveal the partition. The way to
shuffle more messages at a given security level is a larger ring. At
`DEGREE = 2048` the MSIS lattice dimension doubles, which dominates the single
bit the bound gains, and `MSGS = 1000` sits at 243 bits instead of 100.
`shuffle.c` asserts both halves of this at compile time, and section 5c of
`LNP-PARAMS.md` has the tables.

## 7. The orderings

Most of the soundness is carried by what is fixed before what. In order:

1. The commitment `com` and the mask commitment `mcom` are formed.
2. `proj_seed` hashes both, and the projection rows and the scalars `mu`,
   `nu` are derived from that seed. So the masks are fixed before the scalars
   that will be applied to them.
3. The prover publishes the projection `z_p`, and `proj_public` derives
   `P_j`, `M_j`, `Z_j` from it.
4. `rho_hash` draws the batching challenge over every commitment. So every
   message's commitment is fixed before the weights applied to it.
5. The masks are drawn, the garbage terms are formed with those weights and
   committed, and the aggregated values `h_j`, the batched term `T` and the
   first messages are hashed by `batch_hash` into the challenge `d`. So the
   garbage is fixed after `rho` and before `d`.
6. Responses are computed, one rejection test is run over all of them, and the
   verifier checks norms, the Ajtai equations of all three commitments, the
   batched product identity, the constant coefficients and the batch rows.

Step 2 is tested by a test that plays a prover reading the scalars and then
adapting its masks to them. Step 5 is tested by breaking one message's product
in a batch of 25 while leaving its constant coefficient alone, which only the
batched quadratic relation can catch.

## 8. Where it lives

| function | role |
| --- | --- |
| `lnp_auto`, `lnp_auto_crt`, `lnp_auto_swaps` | the automorphisms, in both representations |
| `lnp_keyinit`, `lnp_keygen`, `lnp_commit` | the multi-slot commitment |
| `lnp_isbin_product`, `lnp_ones` | `f = sigma(s)(s-1)` and the all-ones element |
| `lnp_sample_proj_mask`, `proj_seed`, `proj_pass`, `proj_scalars`, `proj_public` | the projection and its public multipliers |
| `lnp_bin_setup` | projection, scalars, and this message's share of `h` |
| `lnp_bin_first` | garbage terms, the challenge-free term `t`, and this message's share of `v` |
| `lnp_bin_public` | what the verifier rebuilds rather than receives |
| `lnp_bin_check` | per-message checks, and accumulation into `acc` |
| `lnp_maskkey_*`, `lnp_mask_commit`, `lnp_batch_*` | the batch-wide masks and the final rows |
| `lnp_garbkey_*`, `lnp_garb_commit`, `lnp_garb_first` | the batch-wide garbage terms |
| `rho_hash` in `shuffle.c` | the batching challenge |

The split into `setup` / `first` / `check` exists so that the caller owns the
mask and the challenge: `shuffle.c` runs every message's setup, then every
message's first message, then derives one challenge, then collects every
response.

## 9. Parameters

| | | why |
| --- | --- | --- |
| `DEGREE` | 1024 | ring degree |
| `MODP` | `2^40 + 141` | set by the range proof's no-wraparound condition |
| `LNP_LAMBDA` | 4 | `p^-4` is about `2^-160` |
| `SLOTS` | 3 | witness, product, projection mask |
| `LNP_RANK` | 2 | rank 1 hiding is only 73 bits |
| `PROJ` | 256 | the projection lemma's requirement for `2^-128` |
| `TAU_PROJ` | 9 | about 3.8 repetitions, bound 1.78x inside the ceiling |
| `SIGMA_B` | `SIGMA_C sqrt(MSGS)` | one rejection test spans the batch |

The modulus is the parameter that drags everything: it forced `WIDTH` up, `DIM`
up, `q` against its 64-bit ceiling, and exposed two latent 32-bit assumptions
in arithmetic that had been correct for years. `LNP-PARAMS.md` section 4 has
that history.

## 10. What is not established

The tests establish completeness, that two specific attacks fail by the checks
meant to catch them, and that several specific forgeries against the
machinery are rejected. They do not establish that no attack succeeds. Specifically:

- The extraction argument for the product relation divides by
  `(c - c')(sigma(c) - sigma(c'))` and yields a *relaxed* opening. What that
  relaxation does to the binding bound has not been worked out, and the MSIS
  figures quoted here are for exact openings.
- The projection lemma is used with a rough constant of 2.
- Zero-knowledge is sketched, not proven. `SOUNDNESS.md` section 9.
- Nothing here has been reviewed by a second person.
