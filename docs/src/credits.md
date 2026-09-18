# Credits

```@eval
using Quadriceps, Markdown
rs = vcat(Quadriceps.available(:gh), Quadriceps.available(:le))
k(prefix) = count(r -> startswith(r.origin, prefix), rs)
Markdown.parse("$(k("own")) of the $(length(rs)) rules were computed from scratch by the author. A further $(k("derived")) (Legendre) rules were obtained by node elimination started from Diallo and Worku's published rules. Finally, $(k("transcribed") + k("same-rule")) are rules from the literature (copied in, or found again by the author's search and recognized).")
```

The `origin` of each rule ([Stored rules](rules.md),
[`Quadriceps.ruleinfo`](@ref), `data/index.tsv`, and the `source_id` in `data/rules.bin`) is one of:

* `own` — computed by the author, with no published rule as its starting point;
* `transcribed: …` — the published rule itself, copied in;
* `same-rule: …` — the author's search converged to a rule identical to a published one
  (matched node for node; for GH up to a rotation, which the Gaussian weight permits). The rule
  belongs to the cited source;
* `derived: …` — the author's node elimination started from a published rule for the same
  ``(d, p)`` and went below its node count. The count is new; the starting point is theirs.

**If you use a rule whose origin is not `own`, cite the source named there.** The file
`NOTICE.md` at the top of the repository carries the license notice that applies to the
derived files.

## Sources of the stored rules

* R. Cools and A. Haegemans, Another step forward in searching for cubature formulae with a
  minimal number of knots for the square, *Computing* **40** (1988) 139–146; and Construction
  of symmetric cubature formulae with the number of knots (almost) equal to Möller's lower
  bound, in *Numerical Integration III*, ISNM, Birkhäuser (1988) 25–36. — GH ``d = 2``,
  ``p = 13``.
* M. Diallo and Z. A. Worku, High-order symmetric positive interior quadrature rules on two and
  three dimensional domains, arXiv:2601.14488 (2026); data at
  `github.com/mdiallo-fula/SymmetricPositiveInteriorCubatures.jl` (MIT license). — starting
  points of Le ``d = 2``, ``p = 77`` and ``d = 3``, ``p = 27, 29, …, 45``.
* M. Festa and A. Sommariva, Computing almost minimal formulas on the square, *J. Comput. Appl.
  Math.* **236** (2012) 4296–4302. — Le ``d = 2``, ``p = 9, 13, 15, 25``.
* A. Haegemans and R. Piessens, Construction of cubature formulas of degree eleven for
  symmetric planar regions, using orthogonal polynomials, *Numer. Math.* **25** (1976) 139–148.
  — GH ``d = 2``, ``p = 11``.
* A. Haegemans and R. Piessens, Construction of cubature formulas of degree seven and nine for
  symmetric planar regions, using orthogonal polynomials, *SIAM J. Numer. Anal.* **14** (1977)
  492–508. — GH ``d = 2``, ``p = 9``.
* S. I. Konyaev, Ninth-order quadrature formulas invariant with respect to the icosahedral
  group (in Russian), *Dokl. Akad. Nauk SSSR* **233** (1977) 784–787. — GH ``d = 3``,
  ``p = 9``.
* A. H. Stroud, *Approximate Calculation of Multiple Integrals*, Prentice-Hall (1971). — GH
  ``d = 2``, ``p = 5``; ``d = 3``, ``p = 5, 7``; ``d = 4``, ``p = 7``; ``d = 5``, ``p = 3, 7``.
* A. H. Stroud and D. Secrest, Approximate integration formulas for certain spherically
  symmetric regions, *Math. Comp.* **17** (1963) 105–135. — GH ``d = 5``, ``p = 5``.

## Lower bounds

* H. M. Möller, Kubaturformeln mit minimaler Knotenzahl, *Numer. Math.* **25** (1976) 185–200.
* H. M. Möller, Lower bounds for the number of nodes in cubature formulae, in *Numerische
  Integration*, ISNM 45, Birkhäuser (1979) 221–230.

## One-dimensional rules

The one-dimensional Gauss–Hermite and Gauss–Legendre rules, used for ``d = 1`` and as factors
of tensor products, are computed by
[FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl).
