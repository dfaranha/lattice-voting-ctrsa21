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
| `master` | `2^31.86` | 1024 | 3 | 36 | 27.0 KB | 117.4 | 72.7 | 88.8 | **72.7** |
| `nocrt` | `2^31.86` | 1024 | 3 | 36 | 27.0 KB | 117.4 | 72.7 | 88.8 | **72.7** |
| `fix-pkc` | `2^31.86` | 1024 | 3 | 36 | 44.0 KB | 117.4 | 72.7 | 88.8 | **72.7** |
| `balanced-params` | `2^31.25` | 1024 | 4 | 20 | 50.0 KB | 128.2 | 191.6 | 159.7 | **128.2** |
| `lnp` | `2^40` | 2048 | 3 | 18 | 137.0 KB | 352.2 | 140.2 | 326.2 | **140.2** |

Security is core-SVP bits, in the sense of Alkim-Ducas-Pöppelmann-Schwabe,
computed with the `lattice-estimator`. The security level of a branch is the
minimum of its three columns, since an attacker picks the weakest. Sizes are at
`MSGS = 25` and count what the prover sends the verifier: the statement, which
is the input commitments and the shuffled list, is not included.

Five things the table is worth reading for.

**The published parameters give about 73 bits, not the 100 claimed in the
paper.** `master`, `nocrt` and `fix-pkc` are all limited by the hiding property
of the commitment, which is MLWE of rank `WIDTH - HEIGHT - 1`, and at
`WIDTH = 3` that rank is 1. Binding and encryption are comfortable; hiding is
not. `balanced-params` exists to fix this, and does, at `WIDTH = 4`.

**`master` and `nocrt` are the unpatched protocol.** Their proofs are the
smallest in the table and establish the least: the attack of ePrint 2025/658
breaks them, and the test suites on the patched branches demonstrate it
succeeding. The 27.0 KB is a reference point, not a recommendation.

**`fix-pkc` is what the patch costs.** Lemma 5 needs the permutation element
committed and tied to the linear proof, which is one extra commitment and one
extra opening per message: 27.0 KB becomes 44.0 KB, for the same parameters and
therefore the same security.

**`lnp` is the same patch with the membership proof actually proven.** Where
`fix-pkc` establishes the set membership Lemma 5 needs by a norm bound, `lnp`
proves the committed element is binary, which is the statement the lemma really
requires. That needs the LNP machinery, a larger modulus, and a larger ring,
and it costs 137.0 KB. It is the only branch whose size is *measured* rather
than modelled: `serial.c` packs a real proof and the verifier checks the
decoded values before the byte count is reported.

**`DEGREE` is 2048 on `lnp` because the batch size became a security
parameter.** One rejection test spans every message's response, so the masks
widen as `sqrt(MSGS)` and the extractable opening widens with them. At
`DEGREE = 1024` that puts binding at 129 bits for 25 messages and 100 for a
thousand; at 2048 the same two are 303 and 243. Shuffling in blocks is not a
way around this, because `MSGS` is the anonymity set, not merely a batch size.

Every branch omits the first messages from the proof and carries the 32-byte
Fiat-Shamir digest instead, the verifier recovering each first message from the
equation that used to check it. That is worth 25 to 31 per cent and is already
included above.

`LNP-PARAMS.md` on the `lnp` branch has the derivations, the measurements and
the caveats, including which of these numbers are modelled and which are
measured.
