# lattice-voting-ctrsa21

Code accompannying the paper "Lattice-Based Proof of Shuffle and Applications to Electronic Voting" by Diego F. Aranha, Carsten Baum, Kristian Gjøsteen,
Tjerand Silde, and Thor Tunge accepted at CT-RSA 2021. The code has been revised post-publication for quality and performance.

## WARNING

This is an academic proof of concept and has not received code review.
This implementation is NOT ready for any type of production use, in particular because of the attack below.

Bootle, Lyubashevsky and Merino-Gallardo, ["Efficient Verifiable Mixnets from
Lattices, Revisited"](https://eprint.iacr.org/2025/658), showed that our proof
is not sound, and mounted a working attack against this implementation.

**This branch does not contain the fix.** Work on the countermeasure of Lemma 5
of that paper lives on the `fix-pkc` branch, along with a draft of the
corresponding soundness argument.

**A ballot that verifies need not decrypt.** Verifiable encryption bounds a
voter's encryption randomness only through its response, `||z||_oo < 6 sigma_E`,
while decryption needs the resulting noise times `p` to stay below `q / 4`. A
malicious voter who draws that randomness from `[-R, R]` rather than ternary,
and skips the rejection sampling, produces a proof that verifies and a
ciphertext that decryption refuses. Measured, six trials per value:

| branch | decrypts | verifies, does not decrypt | does not verify |
| --- | --- | --- | --- |
| `master`, `fix-pkc` | `R <= 6000` | `7000 <= R <= 20000` | `R >= 25000` |
| `balanced-params` | `R <= 700` | `1000 <= R <= 25000` | `R >= 30000` |
| `lnp` | `R <= 4000` | `5000 <= R <= 25000` | `R >= 30000` |

Decryption never returned a wrong message: it refuses such a ciphertext. So
this costs neither privacy nor integrity, only the malicious voter's own ballot,
and the implementation assumes honest voters or accepts that risk. Closing it
would take `q` around `2^63` here and `2^72` on `lnp`, against `2^56` and
`2^64`, which would weaken encryption; the paper's constraints on `q` cover
honest ciphertexts only.

Honest ballots always decrypt, in the worst case and not merely with high
probability. Decryption multiplies by a difference of two challenges and
accepts below `q / 4`, and with every secret and every piece of randomness
CB(1) that needs `q >= 8 nu p (2 DIM N + 1.5)`: `2^52.0` on `master` and
`fix-pkc`, `2^58.9` on `balanced-params` and `2^60.8` on `lnp`, each at least
one bit below its `q`. The paper's condition, `q > 2p(2 DIM N^2 + N + 1)`, is
looser by its `N^2`, and `balanced-params` was sized against it until its `q`
was set by this bound instead, which put it four bits lower relative to `p`;
that is why its window starts lower than the others.

## Building

Dependencies are the GMP and FLINT libraries. FLINT 3 or later is recommended;
the sources use the context-based `fmpz_mod_poly` interface, so releases older
than FLINT 2.8 will not build.

For building the code, run `make` inside the source directory. This will build the binaries for `commit`, `vericrypt` and `shuffle` to test and benchmark different modules of the code.

Each binary runs its tests and then its benchmarks. Since the benchmarks
dominate the runtime by two orders of magnitude, either phase can be selected
on its own:

    make test     # every binary, tests only
    make bench    # every binary, benchmarks only
    ./shuffle test

The tests are also run under AddressSanitizer and UndefinedBehaviorSanitizer in
CI, which is what keeps the allocation behaviour honest:

    make CFLAGS="-O1 -g -march=native -pthread \
      -fsanitize=address,undefined -fno-sanitize-recover=all" test

## Third-party code

`vcl/` vendors Agner Fog's Vector Class Library, used by the constant-time
discrete Gaussian sampler in `gaussian_ct.cpp` (by Raymond K. Zhao).

## Branches, security levels and proof sizes

The repository carries several parameter sets and several versions of the
proof. This table is the reference for what each costs and what each is worth.

| branch | modulus | `DEGREE` | `WIDTH` | `NONZERO` | proof, per message | MSIS binding | MLWE hiding | MLWE encryption | security level |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `master` | `2^31.86` | 1024 | 3 | 36 | 25.3 KB | 83.2 | 72.4 | 87.0 | **72.4** |
| `fix-pkc` | `2^31.86` | 1024 | 3 | 36 | 41.2 KB | 83.2 | 72.4 | 87.0 | **72.4** |
| `balanced-params` | `2^39` | 1024 | 4 | 20 | 50.6 KB | 131.1 | 142.2 | 139.0 | **131.1** |
| `lnp` | `2^40` | 2048 | 3 | 18 | 132.8 KB | 256.1 | 137.6 | 316.2 | **137.6** |

Security is core-SVP bits, in the sense of Alkim-Ducas-Pöppelmann-Schwabe,
computed with the `lattice-estimator`. The security level of a branch is the
minimum of its three columns, since an attacker picks the weakest.
`estimates.sage` reproduces every column and says which problem each one is:
hiding with the commitment randomness uniform in `{-1, 0, 1}` and encryption
with its secrets centred binomial, as the code samples them, which puts the
hiding and encryption columns a few bits below figures that modelled them
otherwise.

**Binding is quoted at the bound extraction needs, not the paper's.** The
extractor turns two relaxed openings into `delta' zbar - delta zbar'`
(`SOUNDNESS.md`, Section 5, on the patched branches), and the column is MSIS
at its norm, `16 nu sigma sqrt(WIDTH N)`: a challenge difference stretches a
vector by up to its `l1` norm, `2 nu`, and `zbar` has `WIDTH` elements of norm
at most `4 sigma sqrt(N)` each. The paper's `16 sigma sqrt(nu N)` is smaller by
`sqrt(nu WIDTH)` because it charges the challenge its `l2` norm, which does
not bound a product, and charges `zbar` one element. At the paper's bound
binding would read 117.4 on `master` and `fix-pkc`, 181.0 on `balanced-params`
and 327.0 on `lnp`; `estimates.sage` prints it for comparison. Only
`balanced-params` is moved by the difference: at `p = 2^37` it had 121.8, and
`p` was raised to `2^39`, and `q` with it, to reach 131.1. The masked openings
do not depend on `p`, so only the uniform coefficients widen, 49.4 KB a
message becoming 50.6, while hiding and encryption, which a larger modulus
weakens, fall to 142.2 and 139.0. Raising `p` from `2^31.25` to `2^37` earlier
had been for the same reason, at a bound that missed the `sqrt(WIDTH)`.

Sizes are at `MSGS = 25` and are **measured, not estimated**. Every branch
carries a `serial.c` that packs the values the prover sends, at
`ceil(log2 p)` bits per uniform coefficient and, for a Gaussian one, against
the distribution rather than its bound, and its test suite round-trips an
honest proof through it and verifies the *decoded* values before reporting the
byte count. So each figure is the size of something that actually verifies
rather than a count of struct fields. What is counted is what the prover sends: the statement, which is the
input commitments and the shuffled list, is not included, and neither are the
first messages, which the verifier recovers rather than receiving.

Six things the table is worth reading for.

**The published parameters give about 72 bits as per latest estimates, not the 100 claimed in the
paper at the time.** `master` and `fix-pkc` are both limited by the hiding property
of the commitment, which is MLWE of rank `WIDTH - HEIGHT - 1`, and at
`WIDTH = 3` that rank is 1. Binding and encryption are comfortable; hiding is
not. `balanced-params` exists to fix this, and does, at `WIDTH = 4`.

**`master` is the unpatched protocol.** Its proofs are the
smallest in the table and establish the least: the attack of ePrint 2025/658
breaks them, and the test suites on the patched branches demonstrate it
succeeding. The 25.3 KB is a reference point, not a recommendation.

**The sizes are measured, and the masked openings are entropy-coded.** A
Gaussian coefficient used to be written at the width of its bound,
`ceil(log2 12 sigma)`; `serial.c` now zigzags it and Golomb-Rice codes it about
its own width instead, and each branch's own test reports what that comes to
and brackets the byte count by it -- below a flat encoding, above the entropy
`log2(sigma sqrt(2 pi e))` of what it codes.

| branch | achieved | flat | entropy | was | is |
| --- | --- | --- | --- | --- | --- |
| `master` | 17.86 | 20 | 17.77 | 27.0 KB | 25.3 KB |
| `fix-pkc` | 17.09 | 19.20 | 17.00 | 44.0 KB | 41.2 KB |
| `balanced-params` | 16.25 | 18.38 | 16.17 | 50.0 KB | 50.6 KB |
| `lnp` | 20.26 | 21.93 | 20.15 | 137.0 KB | 132.8 KB |

The first three columns are bits a Gaussian coefficient, averaged over the
widths a branch uses where it has more than one. Rice lands within about a
tenth of a bit of the entropy throughout, so what the change is worth per branch is
essentially the gap between the bound and the distribution, and how much of the
proof is Gaussian rather than uniform: close to 6 per cent on every branch but
`lnp`, where the projection and the batch masks make it 3. The "was" column
is at the widths each branch had when the coding landed; `balanced-params` and
`lnp` have since moved theirs, below, and both have moved `p`: `balanced-params`
to `2^39 - 19`, whose uniform coefficients take 39 bits rather than 37, and
`lnp` to `2^40 - 195`, whose take 40 rather than 41. The "is" column also no longer sends the
last `s_i`, which no linear proof reads: one ring element fewer, about 4 KB a
proof at `DEGREE = 1024` and 10 KB on `lnp`.

**`fix-pkc` is what the patch costs.** Lemma 5 needs the permutation element
committed and tied to the linear proof, which is one extra commitment and one
extra opening per message: 25.3 KB becomes 41.2 KB, for the same parameters and
therefore the same security.

**`lnp` is the same patch with the membership proof actually proven.** Where
`fix-pkc` establishes the set membership Lemma 5 needs by a norm bound, `lnp`
proves the committed element is binary, which is the statement the lemma really
requires. That needs the LNP machinery, a larger modulus, and a larger ring,
and it costs 132.8 KB.

**`DEGREE` is 2048 on `lnp` because the batch size became a security
parameter.** One rejection test spans every message's response, so the masks
widen as `sqrt(MSGS)` and the extractable opening widens with them. At
`DEGREE = 1024` that puts binding at 107 bits for 25 messages and 84 for a
thousand; at 2048 the same two are 256 and 208. Shuffling in blocks is not a
way around this, because `MSGS` is the anonymity set, not merely a batch size.

Every branch omits the first messages from the proof and carries the 32-byte
Fiat-Shamir digest instead, the verifier recovering each first message from the
equation that used to check it. That is worth 25 per cent on `lnp` and 30 to 33
on the others, and is already included above.

The CRT splitting of the commitment ring buys no speed. A `nocrt` branch did
the same arithmetic on whole polynomials, and once master stopped hashing each
CRT component as if it were a full polynomial, the two did the same work per
prover attempt, within 2 per cent of instructions retired. FLINT multiplies in
quasi-linear time, so two products of half the degree cost about what one of
the full degree does, and the conversions between the representations take the
rest. The branch is gone; its last state is the tag `archive/nocrt`.

**The mask widths are the minimum rejection sampling needs.** At the code's
`M = 1.75`, rejection sampling keeps its `2^-100` bound once `sigma` is 21.485
times the longest vector it masks, `nu * beta * sqrt(K * N)` for `K` ternary
ring elements. `balanced-params` had inherited a width scaled from the paper's
54000 and now sits at that minimum, which took 0.7 KB off a message and raised
binding from 128.2 to 132.9 bits at the modulus of the time and the paper's
bound. `lnp` went the other way: its batched test
masks `K = 285` ring elements, not the 75 its width was sized for, so `SIGMA_B`
had to grow from 165340 to 295465, at 2.2 KB a message and 25 bits of binding
at the paper's bound that its hiding-limited security level does not notice. `master` and `fix-pkc`
keep the paper's width, with about 26 per cent to spare.

**Timings.** Median of three interleaved rounds, in millions of cycles of the
time-stamp counter, on one core of an Intel Core i7-10510U, each branch's
benchmarks built with the Makefile's flags at `MSGS = 25`. The shuffle row is a
whole run, proving and verifying all 25 messages; the others are per call.

| | `master` | `fix-pkc` | `balanced-params` | `lnp` |
| --- | --- | --- | --- | --- |
| shuffle, prove and verify (25 messages) | 793 | 3044 | 4256 | 5647 |
| per message | 31.7 | 121.8 | 170.3 | 225.9 |
| verifiable encryption | 220 | 233 | 479 | 1376 |
| its verification | 46.1 | 44.8 | 94.4 | 290 |
| its decryption | 54.5 | 55.5 | 115.8 | 337 |
| commitment | 0.77 | 0.76 | 1.68 | 2.17 |
| opening | 2.08 | 2.01 | 3.03 | 5.79 |

The rounds agree within about 10 per cent except one round on `lnp` that ran
twice as slow; the median sets it aside. The shuffle carries the spread of
rejection sampling, since the prover's attempts are random: `lnp`'s one batched
test makes it the widest, as section 5b of its `LNP-PARAMS.md` explains. What
each branch pays for is in the sizes above: `fix-pkc` is the patched protocol
at the published parameters, and costs four times `master` to prove, mostly in
the linear proofs that now relate three commitments under four rejection
tests.

`LNP-PARAMS.md` on the `lnp` branch has the derivations and the caveats.
