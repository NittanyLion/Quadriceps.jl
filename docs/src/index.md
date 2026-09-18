# Quadriceps.jl

Positive-weight cubature rules in several dimensions, for the Gaussian weight and for the
uniform weight on the cube.

A cubature rule of degree ``p`` for a weight ``ω`` on ``\mathbb{R}^d`` is a set of nodes
``x_1, …, x_n`` and weights ``w_1, …, w_n`` with

```math
\sum_{i=1}^n w_i f(x_i) = \int f(x)\, ω(x)\, dx
\qquad \text{for every polynomial } f \text{ of total degree} ≤ p .
```

Positive weights matter in practice: the rule is then a discrete probability distribution (up
to scale), sums do not cancel, and an integrand that is nonnegative has a nonnegative
approximate integral.

The product of one-dimensional Gauss rules is such a rule with ``q^d`` nodes, where
``q = (p+1)/2``. The rules in this package need far fewer:

```@eval
using Quadriceps, Markdown
rows = ["| $(fam ≡ :gh ? "GH" : "Le"), ``d = $d``, ``q = $q`` | $(big(q)^d) | $(Quadriceps.nnodes(fam, d, q)) |"
        for (fam, d, q) in ((:gh, 3, 4), (:gh, 5, 5), (:gh, 5, 11), (:le, 2, 39), (:le, 5, 11))]
Markdown.parse("| cell | product grid | stored rule |\n|:---|---:|---:|\n" * join(rows, "\n"))
```

The package exports two functions, [`ghpos`](@ref) and [`lepos`](@ref), modeled on
`gausshermite(q)` and `gausslegendre(q)` from FastGaussQuadrature.jl: `q` is the number of nodes
of the one-dimensional Gauss rule, and the rule returned has its degree, ``p = 2q - 1``.

```julia
using Quadriceps

X, w = ghpos(3, 4)          # 27 nodes for N(0, I₃), as exact as the 4×4×4 Gauss–Hermite grid (degree 7)
X, w = lepos(2, 5)          # 17 nodes for the uniform density on [0,1]², degree 9
X, w = ghpos(3; p = 7)      # the first rule again, requested by degree
```

* [Guide](guide.md): conventions, the two keyword arguments, accuracy.
* [Stored rules](rules.md): every rule, with its node count, Möller's bound, error and origin.
* [Reference](api.md): docstrings.
* [Credits](credits.md): whose rules these are, and what to cite.
