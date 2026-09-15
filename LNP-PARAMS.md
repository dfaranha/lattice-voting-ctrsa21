# Stage 0 of the LNP plan: what was built and what it found

[ISBIN-PLAN.md](ISBIN-PLAN.md) section 6 says not to start at B1 with the
intent of finishing B6, but to run a spike first: B1, B2 against a toy
statement, and enough of B7 to see whether the parameters close. This is the
report of that spike, and of the constant-coefficient proof that it identified
as the one piece with real design content, which was then built as well.

`is_bin` now holds end to end, and is wired into the shuffle, where it replaces
the relaxed norm proof. Sections 4 and 5 record what that cost.

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

**B3, the constant coefficient** (now `lnp_bin_setup`, `lnp_bin_check` and
`lnp_batch_check`, after the merges of section 5b). Proves
that the constant coefficient of a committed value is zero, which is what turns
the product relation above into `is_bin`. For each of `LNP_LAMBDA` masks `g_i`,
sampled uniformly subject to `ct(g_i) = 0`, the prover publishes
`h_i = g_i + mu_i * f` for a scalar `mu_i`, and the verifier checks
`ct(h_i) = 0`. Tying `h_i` back to the committed values is linear, so it needs
no garbage terms: with `B_i` the sum of the mask commitment's `b2[i]` and
`mu_i * b2[SLOT_F]`, the value `T_i - h_i` is a commitment to zero under
`B_i`.

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
Section 5b is that observation collected.

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

**The binding row no longer holds.** It was computed with the masks at
`SIGMA_C`, and since `d3491d0` they are at `SIGMA_B`. Section 5c has the
recomputation: 128.8, not 162.1, and encryption is no longer the constraint.
The other two rows are unaffected and stand.

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

## 5. B6: wiring into the shuffle

The permutation elements now live in the LNP commitment rather than in a
separate BDLOP one. That is the decision that shaped the rest: sharing a single
commitment between the linear proof and `is_bin` means the two speak about the
same object, so no linking proof between two commitments is needed. The linear
proof's `p` argument became an `lnpcom_t`, its randomness LNP_WIDTH wide, and
its message row `b2[SLOT_S]`.

The `SIGMA_S` machinery is gone. `lin_prover` and `lin_verifier` no longer
carry `ys`, `vs` or `sig`, and the check that tied a narrow `z_sigma` to the
committed message is deleted. Membership in the set `D` that Lemma 5 needs is
now established by `is_bin` against the same commitment, with `D` the binary
ring elements. The monomials `g(i) = X^i` are themselves binary, so the
encoding did not have to change.

Two orderings are load-bearing. The shuffle hash absorbs only the Ajtai part of
the commitment and slot `SLOT_S`, which leaves the garbage and mask slots free
to be rewritten; without that, the range proof's rejection sampling would
change `beta` on every retry. And within the prover the product proof runs
first, since it fills the garbage slots, then the range proof, then the
constant-coefficient proof, because the last two derive their challenges from
the finished commitment.

### It rejects for the right reason

The attack tests pass, but passing is not the same as passing for the right
reason, so the verifier was instrumented to report each check separately:

| | lin | isbin | ct | range |
| --- | --- | --- | --- | --- |
| honest witness | 1 | 1 | 1 | 1 |
| CRT-mixed sigma | **1** | **1** | **0** | **0** |

`lin = 1` on the attack confirms the old `SIGMA_S` guard really is gone rather
than still firing by accident. `isbin = 1` is correct: the product relation
holds, because the prover computed it honestly; what fails is that the constant
coefficient of that product is not zero, which is `ct`. The mixed element is
also far too large, which `range` catches independently.

### What it costs

Measured on the same machine, three full proofs at `MSGS = 25`:

| | per proof | per message |
| --- | --- | --- |
| before B6 | 2.37 s | 0.095 s |
| after B6 | 17.1 s | 0.68 s |

**7.2 times slower.** Three sub-proofs per message, and the range proof's
rejection sampling retries about 3.8 times, each retry recommitting and
re-running the product proof. Nothing here had been optimised at this point;
sections 5a and 5b supersede these figures, and the number below is 111 ms.

For scale, the paper reports 33 ms per vote.

## 5a. Measured against fix-pkc

Time is the invariant-TSC cycle count the benchmark harness reports, and the
harness is asked for the **mean over nine proofs** rather than the time of one.
That matters here. Since the transcript was batched the prover runs a single
rejection test over every message's response, so the cost of one proof is one
geometric draw; timing single proofs and taking a median across runs gave
per-round ratios anywhere from 0.82 to 3.35, and an earlier version of this
section reported 1.2x from exactly that mistake. `fix-pkc`, with twenty-five
independent rejection tests that average out, does not have the problem, which
is why the two cannot be measured the same careless way.

Runs are also bracketed: `fix-pkc` is measured immediately before and
immediately after each variant, and the mean of the two brackets is the
baseline, so drift in the machine's thermal state over a pair cancels. That
this is necessary is visible in the result: the three baselines come out at
2.17, 2.20 and 2.21 seconds and the three ratios at 1.25, 1.28 and 1.30,
where the same quantity measured without bracketing ranged from 1.22 to 1.46
depending on what else had been running.

The run length matters too, and shorter is better rather than worse. Asking
the harness for sixteen proofs instead of nine makes each run long enough to
heat the machine within itself, and since the TSC ticks at a fixed rate
whatever the core does, a throttled run genuinely costs more ticks. That
penalises the slower binary more, and inflated this ratio to 1.40.

Size is counted from the transmitted structures, since nothing here serialises
a proof: a uniform ring element costs `DEGREE` times `ceil(log2 p)` bits and a
Gaussian one `DEGREE` times `ceil(log2 12 sigma)`. Note that `MODP` is 141
above `2^40`, so a uniform element costs 41 bits per coefficient and not 40.

| | fix-pkc | lnp at `a63d9b3` | lnp at `d3491d0` | lnp now |
| --- | --- | --- | --- | --- |
| modulus, `WIDTH` | `2^31.86`, 3 | `2^40`, 4 | `2^40`, 4 | `2^40`, 4 |
| prover, per proof | 2.19 s | 14.99 s | 5.68 s | 2.77 s |
| prover, per message | 88 ms | 600 ms | 227 ms | 111 ms |
| proof, per message | 64.0 KB | 330.1 KB | 188.8 KB | 119.9 KB |
| proof, 25 messages | 1.56 MB | 8.06 MB | 4.61 MB | 2.93 MB |
| against fix-pkc | | 6.8x, 5.2x | 2.6x, 3.0x | **1.28x, 1.9x** |

`is_bin` as first wired in cost 6.8 times the prover and 5.2 times the
proof. Batching it, which is section 5b, brought that to **1.28 times the
prover and 1.9 times the proof**. Almost all of the gap was redundancy rather
than the argument: at `a63d9b3` the same commitment was opened four times,
once by the shuffle's linear proof and once by each of the three sub-proofs,
and every one of those openings carried its own mask, its own challenge and
its own rejection sampling loop.

The last column includes the mask-commitment opening added in `472edfc`, which
costs 1.04x on its own, within the noise of these measurements.

Where the 119.9 KB goes now, per message:

| | |
| --- | --- |
| sigma commitment, one Ajtai part and five message slots | 30.8 KB |
| the one masked opening, `WIDTH` twice and `LNP_WIDTH` once | 44.0 KB |
| first messages | 25.6 KB |
| product commitment and partial product | 15.4 KB |
| published projection | 0.5 KB |
| share of the batch-wide mask commitment, values and opening | 3.6 KB |

## 5b. Batching it down

`is_bin` as first wired in was three separate sub-proofs per message, each
with its own mask, its own challenge and its own opening of the same
commitment under the same randomness, and the shuffle's linear proof opened
that commitment a fourth time. Five changes removed the redundancy:

1. `bb3f90e` merged the constant-coefficient and range proofs. They were being
   handed the *same* masks, so the difference of their two aggregated values
   cancelled `g_j` and revealed a linear function of the witness. Merging
   them fixed that, and removed one of the two.
2. `358d91c` put the whole `is_bin` argument under a single opening.
3. `5ce854b` shared that opening with the shuffle's linear proof, which opens
   the same commitment under the same randomness.
4. `d3491d0` derived one challenge over every message, so one rejection test
   covers the whole batch.
5. This change gives the batch one set of constant-coefficient masks.

The last one is the only one that exploits the batch rather than the message.
The aggregation is over scalars: `h_j` is `g_j` plus a scalar combination of
the statements, and the scalars are drawn after `g_j` is committed. Nothing in
that requires the statements to belong to one message, so one set of
`LNP_LAMBDA` masks covers all 25. They moved out of the per-message commitment
into a commitment of their own, which takes `SLOTS` from 9 to 5, `LNP_WIDTH`
from 12 to 8, and the aggregated values `h` and `v` from one pair per message
to one pair for the batch. The mask commitment, the aggregated values and the
mask opening cost 3.6 KB per message once spread over 25, and remove 72.5, so
the proof falls by 68.9 KB per message.

The prover gains more than the size does, 2.6x against 1.6x, because the
inner products in `lnp_bin_first` and `lnp_bin_check` cost `SLOTS` times
`LNP_WIDTH` multiplications, which falls from 108 to 40.

### The mask commitment has to be opened, and that was missed

Moving the masks into a commitment of their own means that commitment has to
be opened like any other, and the first version of this change did not do it.
Neither half of the opening was checked: there was no Ajtai first message and
so no equation tying `z_mask` to `mcom`, and `z_mask` was never required to be
short.

That is not a small omission. The batch rows are `LNP_LAMBDA` equations in
`MASK_WIDTH` unknowns, and the key is in Hermite normal form, so `b2[j]` is
`e_{1+j}` plus a tail on the last two coordinates. A prover could therefore
set that tail to zero, choose *any* aggregated values `h` with zero constant
coefficient, and read off the `z_mask` that satisfies every row by a single
assignment per row. No lattice problem stands in the way. Since `ct(h_j) = 0`
is the entire content of the constant-coefficient argument, that argument
certified nothing.

It was checked by instrumenting the verifier to replace `h` with zero and
solve the rows for `z_mask`: the consistency test still passed. With the Ajtai
first message added to `lnpbatch_t` and to the challenge hash, and with the
norm bound applied to `z_mask`, the same forgery is rejected by both checks
independently. Two tests now cover it. The first perturbs coordinate 0 of the
opening, which the Ajtai row reaches and the message rows do not, so it breaks
only the commitment equation and isolates it. The second supplies a long
opening; that one is rejected by both checks, and isolating the norm bound the
way the first isolates the equation would mean exhibiting a long vector in the
kernel of the whole key, which is the MSIS problem the binding rests on.

The two attack tests did not catch this, because the CRT-mixing attacks are
rejected by the per-message half of `is_bin`, which was never affected. A test
suite that covers the attacks you thought of does not cover the ones you
introduced.

### The rejection test got noisier

One challenge over every message means one rejection test over the
concatenation of their responses, so the number of repetitions a proof needs
is a single geometric draw instead of 25 independent ones that average out.
The mean fell and the variance rose with it: timing one proof at a time gave
ratios against `fix-pkc` ranging from 0.82 to 3.35.

This is a property of the scheme and not only of the measurement, and it is
worth stating as such. A prover that is 1.3 times `fix-pkc` on average is not
1.3 times on every proof, and the tail is one-sided: rejection can only add
repetitions. `fix-pkc` has no equivalent, because averaging 25 independent
tests is what keeps its own spread to a few per cent. Section 5a explains what
measuring the mean rather than the median costs in method.

### What is left

The proof is now 14 uniform ring elements and 16 Gaussian ones per message.
The largest single item is the commitment `p_l`, at six uniform elements or
30.8 of the 119.9 KB. Two of its five slots hold the garbage terms of the quadratic
proof. Those are filled at first-message time and depend on that message's own
mask, so they look per message; but every message now answers one challenge,
which is the condition under which LNP22 batches garbage terms across
statements. Whether that applies here has not been checked, and it is the
obvious next thing to look at.

## 5c. What batching cost the parameters

The table in section 4 quotes MSIS binding at 162.1 bits core-SVP, computed
with the extractable opening `16 sigma_C sqrt(nu N)`. That was right when it
was written and is not right now. Since `d3491d0` one rejection test spans the
concatenation of every message's response, so every mask is drawn at
`SIGMA_B`, five times `SIGMA_C`, and the verifier's norm bounds moved with
them. The bound the extractor meets moved too, and had not been recomputed.

Recomputed with the `lattice-estimator`, taking the `SIGMA_C` row first as a
control against the number already recorded:

| MSIS binding | `beta` | core-SVP | MATZOV |
| --- | --- | --- | --- |
| as documented, at `SIGMA_C` | `2^27.31` | 162.1 | 187.3 |
| as built, at `SIGMA_B` | `2^29.63` | **128.8** | **155.5** |

The control reproduces 162.1 and 187.3 exactly, so the 33-bit drop is real and
not a change of method. MSIS is flat in the commitment width, so the same
figure covers all three commitments: the shuffle's at `WIDTH` 4, the
per-message one at `LNP_WIDTH` 8 and the mask one at `MASK_WIDTH` 7.

### MSGS became a security parameter

`SIGMA_B` is `SIGMA_C sqrt(MSGS)`, because the term one rejection test has to
mask is `sqrt(MSGS)` times longer than a single message's. So the batch size
now feeds straight into the binding bound, which it did not before:

| `MSGS` | `SIGMA_B` | `beta` | core-SVP | MATZOV |
| --- | --- | --- | --- | --- |
| 1 | 54000 | `2^27.31` | 162.1 | 187.3 |
| 25 | 270000 | `2^29.63` | 128.8 | 155.5 |
| 100 | 540000 | `2^30.63` | 117.1 | 144.3 |
| 1000 | 1707630 | `2^32.29` | 100.4 | 128.5 |
| 10000 | 5400000 | `2^33.95` | 86.4 | 115.1 |

This is the part worth carrying into any writeup. The speed in section 5a is
quoted at `MSGS = 25`, and it is not scale-free: an electorate of ten thousand
shuffled as one batch would sit at 86 bits, not 129.

An earlier version of this section said to shuffle a larger electorate in
blocks. **That was wrong.** `MSGS` is the anonymity set, not merely a batch
size: the proof says the output list is a permutation of the input list, so
blocking proves only that each output block permutes its own input block and
reveals the partition. It buys the security back by giving up the anonymity
the shuffle exists to provide.

The remedy is a larger ring. `SIGMA_C` scales as `sqrt(k N)` and the bound as
`16 sigma sqrt(nu N)`, so `beta` grows linearly in `N` while the MSIS lattice
dimension does too, and the dimension wins comfortably:

| `DEGREE` | `MSGS` 25 | `MSGS` 100 | `MSGS` 1000 |
| --- | --- | --- | --- |
| 1024 | 128.8 | 117.1 | 100.4 |
| 2048 | 303.1 | 278.6 | 242.9 |

The same near-doubling appears between Dilithium2 and Dilithium5, at 1024 and
2048 MSIS rows and 123 and 265 bits, which is one of the points this estimator
was validated against. So `DEGREE = 2048` covers any realistic electorate with
room to spare, at roughly twice the proof size and rather more than twice the
prover. Whether that trade is worth making depends on the target electorate,
and it is a parameter decision rather than a code one.

`shuffle.c` now asserts both halves of this at compile time: that `SIGMA_B` is
wide enough for `MSGS`, which is a completeness requirement and fails loudly,
and that `MSGS` does not exceed the size these numbers were derived at, which
is a soundness requirement and would otherwise fail silently.

### What did not move

Hiding and encryption depend on the modulus and the randomness distribution,
neither of which batching touched, so they are unchanged. Recomputed anyway
rather than assumed:

| | core-SVP | MATZOV |
| --- | --- | --- |
| MLWE commitment hiding, rank 2 | 137.6 | 164.7 |
| MLWE encryption, `DIM` 3 | 129.1 | 157.3 |

Those two are the figures already in section 4, and they stand. Recomputing
them here needed the attack set restricted to primal uSVP and dual, because
the exhaustive sweep does not terminate in reasonable time at rank 2 with
`q = 2^40`; that gives 140.2 and 130.5, which sit 2.6 and 1.4 bits above the
recorded numbers. The difference is the attacks the restricted set leaves out,
not a change in the parameters, so the lower recorded figures are the ones to
quote.

So the scheme as it stands is **128.8 bits core-SVP**, set by MSIS binding,
with encryption just behind it at 129.1. Before batching the same two were
162.1 and 129.1 and encryption was the constraint; MSIS binding has overtaken
it. The paper's "at least 100 bits" is still true at `MSGS = 25`, and it is
the right kind of claim to make precise: it is 129, it is set by the batch
size, and at `MSGS = 1000` it would be exactly 100.

## 6. What this changes about the decision

Track B is complete: B1, B2, B3, B5 and B6 are built and tested, and B4 is not
needed on this route because the projection bounds a BDLOP-committed witness
directly.

The plan sized Track B at 1500-2000 lines and framed the choice as fidelity to
Protocol 1 against that. Both halves of that framing moved. The code came to
roughly 1400 lines, and the piece that looked hardest, the constant
coefficient, had a short answer once the reason it was hard was clear. But the
real price was never the code. It was the modulus, and what the modulus drags
with it: `WIDTH`, `DIM`, `q` against its 64-bit ceiling, the CRT constants, and
two latent 32-bit assumptions in arithmetic that had been correct for years.
Track A would have paid the same price, since LaZer's range proof obeys the
same arithmetic.

The prover was the other half of that price, and it is no longer. `is_bin`
started at 6.8 times `fix-pkc` and batching brought it to 1.3 times, with the
proof 1.9 times larger; sections 5a and 5b have the numbers. What that cost
turned out to measure was redundancy in how the sub-proofs were wired, not the
argument itself.

So the decision is no longer Track A against Track B, and it is no longer a
question of a 7x prover either. It is whether fidelity to Protocol 1 as
published, and not having to review the restatement of Lemma 5 in
SOUNDNESS.md, is worth a modulus change. That is a judgement about what the
artifact is for.

Finally, none of this is a soundness proof. The extraction argument for the
product relation divides by `(c - c')(sigma(c) - sigma(c'))` and yields a
relaxed opening, and what that relaxation does to the binding bound has not
been worked out. The projection lemma is used with a rough constant of 2. The
tests establish completeness, and that two specific attacks fail, by the checks
that should catch them. They do not establish that no attack succeeds.
