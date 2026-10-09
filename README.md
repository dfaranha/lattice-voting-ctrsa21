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
| `master` | `2^31.86` | 1024 | 3 | 36 | 25.3 KB | 117.4 | 72.7 | 88.8 | **72.7** |
| `fix-pkc` | `2^31.86` | 1024 | 3 | 36 | 41.2 KB | 117.4 | 72.7 | 88.8 | **72.7** |
| `balanced-params` | `2^31.25` | 1024 | 4 | 20 | 46.3 KB | 132.9 | 191.6 | 159.7 | **132.9** |
| `lnp` | `2^40` | 2048 | 3 | 18 | 134.6 KB | 327.0 | 140.2 | 326.2 | **140.2** |

Security is core-SVP bits, in the sense of Alkim-Ducas-Pöppelmann-Schwabe,
computed with the `lattice-estimator`. The security level of a branch is the
minimum of its three columns, since an attacker picks the weakest.

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

**The published parameters give about 73 bits as per latest estimates, not the 100 claimed in the
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
| `balanced-params` | 16.25 | 18.38 | 16.17 | 50.0 KB | 46.3 KB |
| `lnp` | 20.26 | 21.93 | 20.15 | 137.0 KB | 134.6 KB |

The first three columns are bits a Gaussian coefficient, averaged over the
widths a branch uses where it has more than one. Rice lands within about a
tenth of a bit of the entropy throughout, so what the change is worth per branch is
essentially the gap between the bound and the distribution, and how much of the
proof is Gaussian rather than uniform: close to 6 per cent on every branch but
`lnp`, where the projection and the batch masks make it 3. The "was" column
is at the widths each branch had when the coding landed; `balanced-params` and
`lnp` have since moved theirs, below. The "is" column also no longer sends the
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
and it costs 134.6 KB.

**`DEGREE` is 2048 on `lnp` because the batch size became a security
parameter.** One rejection test spans every message's response, so the masks
widen as `sqrt(MSGS)` and the extractable opening widens with them. At
`DEGREE = 1024` that puts binding at 140 bits for 25 messages and 109 for a
thousand; at 2048 the same two are 327 and 261. Shuffling in blocks is not a
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
54000 and now sits at that minimum, which takes 0.7 KB off a message and raises
binding from 128.2 to 132.9 bits. `lnp` went the other way: its batched test
masks `K = 285` ring elements, not the 75 its width was sized for, so `SIGMA_B`
had to grow from 165340 to 295465, at 2.2 KB a message and 25 bits of binding
that its hiding-limited security level does not notice. `master` and `fix-pkc`
keep the paper's width, with about 26 per cent to spare.

`LNP-PARAMS.md` on the `lnp` branch has the derivations and the caveats.
