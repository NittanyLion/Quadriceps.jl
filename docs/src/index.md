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

| cell | product grid | stored rule |
|:---|---:|---:|
| GH, ``d = 3``, ``p = 7`` | 64 | 27 |
| GH, ``d = 5``, ``p = 9`` | 3125 | 244 |
| GH, ``d = 5``, ``p = 21`` | 161051 | 13199 |
| Le, ``d = 2``, ``p = 77`` | 1521 | 1032 |
| Le, ``d = 5``, ``p = 21`` | 161051 | 10984 |

The package exports two functions, [`ghpos`](@ref) and [`lepos`](@ref), modeled on
`gausshermite` and `gausslegendre` from FastGaussQuadrature.jl.

```julia
using Quadriceps

X, w = ghpos(3, 7)          # 27 nodes for N(0, I₃), exact to degree 7
X, w = lepos(2, 9)          # 17 nodes for the uniform density on [0,1]²
```

* [Guide](guide.md): conventions, the two keyword arguments, accuracy.
* [Stored rules](rules.md): every rule, with its node count, Möller's bound, error and origin.
* [Reference](api.md): docstrings.
* [Credits](credits.md): whose rules these are, and what to cite.
