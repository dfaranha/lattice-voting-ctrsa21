# lattice-voting-ctrsa21

Code accompannying the paper "Lattice-Based Proof of Shuffle and Applications to Electronic Voting" by Diego F. Aranha, Carsten Baum, Kristian Gjøsteen,
Tjerand Silde, and Thor Tunge accepted at CT-RSA 2021.

## The proof of shuffle in this branch is not sound

Bootle, Lyubashevsky and Merino-Gallardo, ["Efficient Verifiable Mixnets from
Lattices, Revisited"](https://eprint.iacr.org/2025/658), showed that the proof
of shuffle of the paper above is unsound, and mounted a working attack against
this implementation. The product identity the proof checks implies that the two
lists are related by a permutation only over a field, and the ring used here
never is one: it only implies that the lists are permuted within each CRT
component, possibly by two different permutations.

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
