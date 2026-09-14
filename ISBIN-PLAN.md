# Implementing `is_bin` with LNP: feasibility

## Status

This is a feasibility assessment, not a commitment to build. It estimates what
it would take to replace the relaxed norm proof currently on this branch with
the `is_bin(sigma_i)` sub-proof that Protocol 1 of ePrint 2025/658 assumes,
implemented in the style of Lyubashevsky, Nguyen and Plancon (CRYPTO 2022).

Two facts were checked against this code base rather than assumed; they are
marked **[checked]** below.

## 1. What has to be proven

Lemma 5 needs the committed `sigma_i` to lie in a set whose pairwise
differences are invertible. Protocol 1 takes that set to be the ring elements
with binary coefficients. The standard reduction is:

1. the quadratic relation `<sigma, sigma - 1> = 0` over `Z_q`, and
2. a norm bound on `sigma` strong enough that the sum does not wrap modulo `q`.

Since `s(s - 1) >= 0` for every integer `s`, with equality exactly when
`s` is 0 or 1, the two together force every coefficient to be binary.

The bound needed for (2) is mild. With `d * B * (B + 1) < q` the sum cannot
wrap, which at `DEGREE = 1024` and this modulus gives

    ||sigma||_inf <= 1952

against `||sigma||_inf = 1` for an honest binary vector: a margin of about
1952x. **The norm bound is not the hard part.** The quadratic relation is.

## 2. What LNP needs, and what is here

| component | needed for | present |
| --- | --- | --- |
| automorphisms `sigma_k : X -> X^k` | constant-coefficient extraction | no |
| constant-coefficient / inner-product proof | `<sigma, sigma - 1> = 0` | no |
| quadratic proof with garbage terms | any degree-2 relation | no |
| ABDLOP commitment (Ajtai part) | committing a *short* message | no |
| approximate range proof (JL projection) | the norm bound of (2) | no |
| linear proof, challenge space, rejection sampling | all of the above | **yes** |

The existing `commit.c` gives a BDLOP commitment whose message part carries no
norm guarantee, plus a linear proof. Everything else on that list is absent.

## 3. Two structural findings

**[checked] The automorphism swaps the CRT factors.** Constant-coefficient
extraction uses `sigma_{-1} : X -> X^{-1}`. Over `Z_q[X]/(X^n + 1)` with this
modulus, `sigma_{-1}` maps the factor `X^(n/2) + P0` onto `X^(n/2) + P1`: it
*swaps* the two CRT components rather than fixing them. Verified by applying it
to `irred[0]` and reducing: the image vanishes in component 1, not component 0.

That matters because the arithmetic here is in CRT representation. Applying an
automorphism means exchanging the two components and applying the induced map,
not acting componentwise. It is implementable, but it is exactly the kind of
detail that silently produces a proof that verifies against itself and nothing
else, so it wants a test against a coefficient-representation reference.

**[checked] None of the machinery can be reused.** Searching for automorphism,
quadratic/garbage, range-proof or Ajtai/ABDLOP code in this repository returns
nothing. The linear proof is genuinely all there is.

## 4. Two tracks

### Track A: link against LaZer

LaZer (Lyubashevsky, Seiler, Steuer, ACM CCS 2024) is the reference
implementation of this proof system, and is what Protocol 1 assumes. The work
becomes expressing the `is_bin` statement in its API and linking the
`sigma_i` commitment to it.

* Much smaller: no new proof system.
* Adds a dependency, and a representation-matching problem: LaZer has its own
  ring and commitment conventions, so either the `sigma_i` commitment is one
  LaZer already understands, or a linking proof is needed between the BDLOP
  commitment used by the shuffle and whatever LaZer commits to. That linking
  step is the real unknown and should be scoped first.

### Track B: implement natively

Roughly, in dependency order:

| stage | work | rough size |
| --- | --- | --- |
| B1 | automorphisms on `nmod_poly`, CRT-aware, with a coefficient-rep reference test | ~100 lines |
| B2 | standalone quadratic proof with garbage terms, for a toy relation | ~300 lines |
| B3 | constant-coefficient extraction, tying `<a,b>` to `sigma_{-1}` | ~300 lines |
| B4 | ABDLOP commitment: Ajtai part, opening proof that bounds it | ~400 lines |
| B5 | approximate range proof (projection to `Z^256`, rejection sampling) | ~300 lines |
| B6 | wire into the shuffle: commit `sigma_i`, run `is_bin`, connect to Lemma 5 | ~200 lines |
| B7 | parameter derivation: challenge space, rejection widths, JL soundness, modulus | research |

Order 1500-2000 lines of new and subtle cryptographic code, plus B7, which is
not a coding task. For scale, the entire optimisation, leak, API and CI effort
already applied to this repository was about 1000 inserted lines, and mostly
mechanical.

## 5. Parameter consequences

* **Modulus.** LNP's range proof and soundness analysis want more headroom than
  `q ~ 2^32` comfortably gives. A larger `q` is likely, which cascades into the
  CRT constants, the fold-based reduction, and the ciphertext sizes.
* **This would relax a constraint we already have.** The `SIGMA_S` window of the
  current relaxed proof closes at `DEGREE = 8192` because its ceiling scales as
  `sqrt(q)/sqrt(n)`. A larger `q` widens it. If the degree is going up anyway,
  the two pressures point the same way.
* **Encoding.** With `is_bin` actually proven, `g` could return to the paper's
  binary encoding instead of the monomials used here, and the `SIGMA_S`
  machinery would be replaced outright rather than supplemented.

## 6. Recommendation

Do not start at B1 with the intent of finishing B6. Run a spike first:

* **Stage 0, days not weeks:** implement B1, plus B2 against a toy statement
  such as `<s, s - 1> = 0` for a hand-committed `s`, and attempt B7 far enough
  to see whether the parameters close at a modulus you are willing to adopt.

Stage 0 answers the feasibility question at a small fraction of the cost, and
its two outputs, the automorphism layer and the parameter verdict, are the
inputs everything else depends on. If the parameters do not close at an
acceptable modulus, Track A becomes the only sensible route and no time has
been sunk.

## 7. The decision this is really about

The relaxed norm proof on this branch already closes the soundness gap, under a
restatement of Lemma 5 that is drafted in [SOUNDNESS.md](SOUNDNESS.md) and has
not been reviewed. Implementing `is_bin` buys fidelity to Protocol 1 as
published, and removes the need for anyone to check that restatement.

Whether that is worth 1500-2000 lines of new proof-system code is a judgement
about the artifact's purpose, not about the cryptography: a research artifact
accompanying a paper may be better served by the documented deviation plus a
reviewed lemma than by a second implementation of a proof system that already
exists in LaZer.
