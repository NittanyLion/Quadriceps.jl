# Data format

All rules are stored in one binary file, `rules.bin` (`data/rules.bin` in Quadriceps.jl,
`src/quadriceps/data/rules.bin` in the Python package, `inst/extdata/rules.bin` in the R
package; the three files are byte-identical). Next to it, `index.tsv` is the catalog: one line
of metadata per rule. There are no per-rule text files.

The format, `QUADRICEPS1`, is the `GHBEST1` layout of the designed-quadrature project
(`bestknown_gh/NOTES.md` §7 there) with a family field added to the index. It is specified
here in enough detail to write a reader from scratch in any language. A binary float format
has two hazards, byte order and struct padding, and both are removed: byte order is fixed to
little-endian and declared in the header, and there are no structs, only flat arrays. IEEE 754
binary64 itself is portable.

```
offset              bytes            contents
-----------------------------------------------------------------------------------------
0                   H                header: ASCII, one line terminated by \n
H                   cells * 8 * 8    index: `cells` records of 8 signed 64-bit little-endian integers
H + cells * 64      …                data blocks, in index order
```

## Header

Read bytes up to the first `\n`; compute `H` by scanning for it, never hard-code it. The line
begins with the magic `QUADRICEPS1` and continues as space-separated `key=value` tokens:

| key | value | meaning |
|---|---|---|
| `fmt` | `1` | format version |
| `endian` | `little` | byte order of every integer and float |
| `cells` | e.g. `148` | number of rules |
| `index_fields` | `8` | integers per index record |
| `float` | `binary64` | IEEE 754 double precision |
| `order` | `row-major` | node by node, see below |
| `layout` | `x1..xd,w` | the values of one node |
| `families` | `0:gh,1:le` | the family codes used in the index |

A reader should check the magic, `fmt`, `endian`, `float` and `index_fields`, and take `cells`
from the header.

## Index

One record per rule, 8 × `Int64` little-endian, in this order:

| field | meaning |
|---|---|
| `family` | 0 = GH (Gaussian weight), 1 = Le (uniform weight on the cube) |
| `d` | dimension |
| `p` | degree of exactness (odd) |
| `q` | `(p+1)/2`, the argument of `ghpos(d, q)` and `lepos(d, q)` |
| `n` | number of nodes |
| `offset` | absolute byte offset of the rule's data block from the start of the file |
| `nbytes` | length of the block in bytes, always `n*(d+1)*8` |
| `source_id` | whose rule it is, table below |

Rules appear in the order family, then `d`, then `p`, all ascending. A reader should not rely on
that order; use the index.

## Data blocks

Each block is `n*(d+1)` IEEE 754 binary64 values, little-endian, **row-major**: node `i`
occupies `d+1` consecutive values, `x₁ … x_d` then `w`. No padding, no separators, no alignment
requirements, and nothing between blocks.

All rules are stored in the normalized frame: GH rules integrate against the standard normal
density `N(0, I_d)`, Le rules against the uniform density on `[0,1]^d`, and the weights of every
rule sum to 1.

## `source_id`

| id | meaning |
|---|---|
| 0 | the package author's own rule |
| 3 | derived: the author's node elimination started from Diallo and Worku's (2026) rule for the same cell |
| 10 | Stroud 1971 |
| 11 | Stroud and Secrest 1963 |
| 12 | Haegemans and Piessens 1976 |
| 13 | Haegemans and Piessens 1977 |
| 14 | Cools and Haegemans 1988 |
| 15 | Konyaev 1977 |
| 16 | Festa and Sommariva 2012 |
| 99 | another published source, named in `index.tsv` |

Ids 0 and 10–15 mean what they mean in `GHBEST1`. For ids 10 and up the rule is the published
rule, transcribed or found again by the author's search; the `origin` column of `index.tsv`
says which, and gives the details. See [Credits](credits.md).

## The catalog, `index.tsv`

Tab-separated text; lines starting with `#` are comments, the first other line names the
columns: `family`, `d`, `p`, `n`, `moller` (Möller's lower bound, `-1` where not tabulated),
`relerr` (largest relative monomial error, measured when the data were built), `minweight`,
`interior` (`yes`/`no`), `origin` (text), `source_id`, and `bankfile` and `sha256`, which name
the file of the project's rule bank that the rule was taken from. The packages read the catalog
for the metadata and `rules.bin` for the numbers, and refuse to load if the two disagree.

## Quadruple precision: `rules128.bin` and `index128.tsv`

`rules128.bin` is a second file of the same format with `float=binary128` in its header: every
number is 16 little-endian bytes, the IEEE 754 binary128 (quadruple precision) rounding of the
project's extended-precision rule, so `nbytes = n*(d+1)*16`. A C `__float128`, a Fortran
`real(16)` or a Julia `Float128` array reads a block directly. It holds the cells that have an
extended-precision file whose double-precision rounding is the rule in `rules.bin`, row for
row. `index128.tsv` is its catalog: `family`, `d`, `p`, `n`, `relerr128` (largest relative
monomial error of the binary128 rule, measured in wider arithmetic; the machine epsilon of the
format is `2^-112 ≈ 1.93e-34`), and `extendedfile` and `sha256`, which name the project's file
the numbers were rounded from. The Julia package reads both; the Python and R packages do not
ship them.

## A reader in a few lines

Python, standard library only:

```python
import struct
data = open("rules.bin", "rb").read()
h = data.split(b"\n")[0] + b"\n"                          # header; its length is H
kv = dict(t.split(b"=") for t in h.split()[1:])
assert h.startswith(b"QUADRICEPS1") and kv[b"endian"] == b"little" and kv[b"index_fields"] == b"8"
for i in range(int(kv[b"cells"])):
    fam, d, p, q, n, off, nb, sid = struct.unpack_from("<8q", data, len(h) + 64 * i)
    vals = struct.unpack_from("<%dd" % (n * (d + 1)), data, off)
    nodes = [vals[j * (d + 1): j * (d + 1) + d] for j in range(n)]
    w = [vals[j * (d + 1) + d] for j in range(n)]
```

Quadriceps.jl reads and writes the format in `src/index.jl` (`readbinindex`, `readblock`,
`writebin`), converting with `ltoh` and `htol`, which are no-ops on a little-endian machine and
byte swaps on a big-endian one.
