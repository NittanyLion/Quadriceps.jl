# Third-party notice

Most rules in `data/` were computed by the package author. Some are, or descend from, rules
published by others. `data/index.tsv` (column `origin`) and the header of every rule file say
which; `Quadriceps.ruleinfo` reports it at run time. The relations are:

* **transcribed** — the published rule itself, copied in; no search of the author's was involved;
* **same-rule** — the author's search converged to a rule identical to a published one (matched
  node for node, for GH up to a rotation); the rule belongs to the cited source;
* **derived** — the author's node elimination started from a published rule for the same
  `(d, p)` and went below its node count; the count is new, the starting point is theirs.

Anyone using one of these rules should cite the source named in its `origin`. Full references
are in `docs/src/credits.md`.

## Rules derived from Diallo and Worku's data (MIT license)

The Le rules `d = 2, p = 77` and `d = 3, p = 27, 29, …, 45` descend from the rules of
M. Diallo and Z. A. Worku, *High-order symmetric positive interior quadrature rules on two and
three dimensional domains*, arXiv:2601.14488 (2026), data at
`github.com/mdiallo-fula/SymmetricPositiveInteriorCubatures.jl`, which is released under the
MIT license. Those files are modified copies in the sense of that license, whose notice
follows and applies to them:

```
MIT License

Copyright (c) 2025 Moustapha Diallo

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## Rule transcribed from Festa and Sommariva

The Le rule `d = 2, p = 25` (113 nodes) is rule SMR25 of M. Festa and A. Sommariva, *Computing
almost minimal formulas on the square*, J. Comput. Appl. Math. 236 (2012) 4296–4302, taken from
J. Burkardt's `square_minimal_rule` port (distributed by its author under the MIT license) and
mapped to `[0,1]²` with weights summing to 1. The Le rules `d = 2, p = 9, 13, 15` coincide with
their SMR09, SMR13 and SMR15 but were produced by the author's own solver; nothing was copied.

## Rules transcribed from closed forms in the literature

GH `d = 2, p = 9` (Haegemans and Piessens 1977), `d = 2, p = 11` (Haegemans and Piessens 1976),
`d = 2, p = 13` (Cools and Haegemans 1988), `d = 3, p = 5` (Stroud 1971) and `d = 5, p = 5`
(Stroud and Secrest 1963) are the published rules, copied in.
