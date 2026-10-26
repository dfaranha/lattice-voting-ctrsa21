# lattice-voting-ctrsa21

Code accompannying the paper "Lattice-Based Proof of Shuffle and Applications to Electronic Voting" by Diego F. Aranha, Carsten Baum, Kristian Gjøsteen,
Tjerand Silde, and Thor Tunge accepted at CT-RSA 2021.

## Soundness fix in the proof of shuffle

Bootle, Lyubashevsky and Merino-Gallardo, ["Efficient Verifiable Mixnets from
Lattices, Revisited"](https://eprint.iacr.org/2025/658), identified a soundness
issue in the proof of shuffle of the CT-RSA 2021 paper, and mounted a working
attack against this implementation.

The proof followed Neff's paradigm and checked the product identity

    \prod (a_i - X) = \prod (b_i - X) ,

which implies that the two lists are related by a permutation only when the
elements live in a field. The ring `R_p = Z_p[x]/(x^N + 1)` used here is never a
field: `x^N + 1` always factors, and with `p = 3906450253 = 5 (mod 8)` it splits
into the two factors that the CRT representation of `commit.c` works with. The
identity therefore only implies that the lists are permuted *within each CRT
component*, possibly by two different permutations. An adversary can take two
messages, exchange just their first CRT component, and obtain a list that is not
a permutation of the original one but still satisfies the identity.

`shuffle.c` now uses the two-variable product of Lemma 5 of that paper, which
ties every message to its index:

    \prod (a_i + g(i) * X1 - X2) = \prod (b_i + sigma_i * X1 - X2) ,

where `g` maps an index to a monomial `x^i` and `sigma_i = g(pi(i))` encodes the
secret permutation. The prover commits to the `sigma_i` *before* the challenges
`X1` and `X2` are drawn. Because the differences of distinct elements in the
image of `g` are short, they are invertible, and the per-component permutations
are forced to agree.

Concretely, this changed the following:

* The prover sends an extra list of commitments to the `sigma_i`, and the
  verifier draws two challenges (`tau` for `X1`, `rho` for `X2`) instead of one.
* The linear proof relates three commitments instead of two. The factor
  `b_i = _m_i + sigma_i * tau - rho` is no longer public: its public part is
  folded into the additive term, and the committed `sigma_i` enters with the
  public coefficient `tau`.
* Challenges are now sampled of degree below `N/2`, so that the difference of
  two distinct challenges is invertible, as the soundness argument requires.
  (The linear proof's challenge has since moved to the paper's challenge space,
  `NONZERO` coefficients of `+-1`, whose differences are invertible because
  they are short; `tau`, `rho` and `beta` keep degree below `N/2`.)
* The linear proof additionally masks the message committed alongside the
  randomness, and the verifier bounds the norm of that response. This is what
  places `sigma_i` in the set `D` that Lemma 5 requires (see below).
* Verifier norm checks are ordinary checks rather than `assert`s, which are
  soundness checks that used to vanish entirely under `NDEBUG`, and they no
  longer overflow when handed a dishonestly large response.

Four tests were added to `shuffle.c`, covering both attacks and the premise of
each:

| test | what it establishes |
| --- | --- |
| CRT-mixed list is not a permutation but passes Neff's product | the attack premise: the old identity accepts a non-permutation |
| shuffle proof rejects the CRT-mixing attack | the Lemma 5 product catches it |
| CRT-mixed sigma is outside D | the second-order attack premise |
| shuffle proof rejects the CRT-mixing attack on sigma | the norm check on `sigma_i` catches it |

The two attacks are rejected by two different mechanisms, which was confirmed by
instrumenting the verifier: the first fails the product relation, the second
fails only the norm check on `sigma_i`.

## Proving `sigma_i` in D by a norm bound

Lemma 5 requires the committed `sigma_i` to lie in a set `D` whose pairwise
differences are invertible. Protocol 1 of the same paper discharges this with a
sub-proof of `is_bin(sigma_i)` delegated to a general-purpose lattice proof
system (LaZer), which has no counterpart in this code base.

Binary coefficients are more than the lemma needs. By Lemma 1, any element of
l2-norm below `MODP^(1/2) ~ 741455` is invertible, so taking `D` to be a ball of
small norm satisfies the requirement, and shortness is exactly what
Fiat-Shamir-with-aborts proves natively. The linear proof therefore masks the
committed message with a narrow Gaussian (`SIGMA_S`, see `param.h`) and the
verifier bounds the norm of the response.

`g` maps indices to monomials rather than to binary representations because
multiplication by a monomial is a signed rotation, so `||d * x^i|| = ||d||`
exactly, which keeps `SIGMA_S` at its minimum whatever the number of messages:

| quantity | value |
| --- | --- |
| invertibility ceiling `MODP^(1/2)` | 741455 |
| extracted `||(d - d') * sigma_i|| <= 2 * (2 * sqrt(DEGREE) * SIGMA_S)` | 12416 |
| masking ratio `SIGMA_S / ||d * sigma_i||` | ~21.7 |

The extracted norm sits a factor of 59.7 below the ceiling. That margin was
4.07 at the previous `MODP = 2553802997`, where the ceiling was 50535; it
widened when `MODP` was raised to about `2^39` for binding (see
[PARAMS.md](PARAMS.md)), not by design. A binary encoding would also fit at
`N = 25`: an index of weight `w` needs `21.485 * sqrt(20) * w <= SIGMA_S <=
(MODP^(1/2) - sqrt(80) * w) / (4 * sqrt(DEGREE))`, so with up to four bits set
the window is `[384.3, 5792.3]`, a factor of about 15 (it was 1.03 at
`2^31.25`). But it would need `SIGMA_S` about four times wider, and
more as `N` grows, since `||sigma_i||_1` grows with the bit-length of the
index. The monomial encoding is independent of `N` and supports up to `DEGREE`
messages.

### Caveat: this deviates from Protocol 1

This is a change to the protocol, not just an implementation of it, and the
soundness argument changes with it. The proof is relaxed: what an extractor
obtains is not "`sigma_i` is in `D`" but "`(d - d') * sigma_i` is short" for a
challenge difference. The step of Lemma 5 that lifts `sigma_i = g(j)` from one
CRT component to the whole ring has to be redone accordingly:

> if `sigma_i = g(j) mod p_l`, then `(d - d') * (sigma_i - g(j))` is short and
> vanishes mod `p_l`, so by Lemma 1 it is zero; and since `d - d'` is
> invertible, `sigma_i = g(j)` over `R`.

The zero-knowledge argument needs updating too, since the proof now carries an
extra masked value: the simulator sets `v_sigma = <b2, z_p> + z_sigma - d * p.c2`
from responses it has already sampled.

Both arguments go through, but they restate results of the paper and **have not
been reviewed**; they are drafted in full in [SOUNDNESS.md](SOUNDNESS.md).
Anyone relying on this code should check that draft before trusting the result.
Keeping literal fidelity to Protocol 1 instead would mean an exact binary proof
via the LNP22 automorphism machinery, which this code base does not have.

## Building

Dependencies are the GMP and FLINT libraries. FLINT 3 or later is recommended;
the sources use the context-based `fmpz_mod_poly` interface, so releases older
than FLINT 2.8 will not build.

For building the code, run `make` inside the source directory. This will build the binaries for `commit`, `vericrypt` and `shuffle` to test and benchmark different modules of the code, and `voting`, which runs the voting protocol end to end.

`voting` puts the three modules together, with `MSGS` voters (25 by default):
each voter commits to a vote for one of 8 candidates and verifiably encrypts
the commitment randomness; the board verifies every ballot against its
commitment; the authority decrypts each ballot, recovers the randomness from
the relaxed witness that decryption returns (multiplying by the inverse of the
challenge difference), checks the opening, shuffles the votes and publishes
them with a proof of shuffle, as bytes; and a verifier checks that proof from
the bytes alone and tallies. It also checks that a ballot copied onto another
voter's commitment, a changed published vote, a corrupted proof and a
substituted vote are all rejected. It runs in about 25 seconds, and `make test`
includes it.

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

The proof decoder, which reads what an untrusted prover sends, has a libFuzzer
harness that needs clang. Besides the sanitizers, it checks that any input the
decoder accepts re-encodes to the same bytes, so that no proof has two
encodings:

    make fuzz_serial
    ./fuzz_serial -seeds fuzz-corpus && ./fuzz_serial -max_len=65536 fuzz-corpus

## Memory ownership

Aggregates follow the same convention as FLINT polynomials: an explicit
initializer, then any number of operations, then a release. `commit_init`,
`commit_keyinit`, `encrypt_cipher_init`, `encrypt_keyinit` and `vericrypt_init`
allocate; `commit_free`, `commit_keyfree`, `encrypt_free`, `encrypt_keyfree` and
`vericrypt_free` release. The computation routines, such as `commit_doit` and
`vericrypt_doit`, allocate nothing and may be called repeatedly on the same
object.

WARNING: This is an academic proof of concept, and in particular has not received code review. This implementation is NOT ready for any type of production use.

## Third-party code

`vcl/` vendors Agner Fog's Vector Class Library, used by the constant-time
discrete Gaussian sampler in `gaussian_ct.cpp` (by Raymond K. Zhao).
