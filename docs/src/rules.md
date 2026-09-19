# Stored rules

Written by `build/build_data.jl`. `q`: the argument of `ghpos(d, q)` and `lepos(d, q)`;
`p = 2q - 1`: degree of exactness; `n`: number of nodes; `ρ = n^(1/d) / q`: the node count
relative to the `q^d` Gauss product grid (1.00 is that grid, smaller
is better); Möller: Möller's lower bound on `n` (**bold** `n`: bound attained, proven minimal);
rel. err.: largest relative monomial error of the stored `Float64` rule; origin: `own`, or the
published rule it is, or descends from (see [Credits](credits.md)).

## GH — Gaussian weight

### d = 2

| q | p | n | ρ | Möller | rel. err. | origin |
|---:|---:|---:|---:|---:|---:|:---|
| 1 | 1 | **1** | 1.00 | 1 | 0.0e+00 | own |
| 2 | 3 | **4** | 1.00 | 4 | 2.2e-16 | own |
| 3 | 5 | **7** | 0.88 | 7 | 2.2e-16 | same-rule: Stroud 1971, E_2^{r^2}:5-4 (identified 2026-09-11; same up to rotation, which the Gaussian weight permits) |
| 4 | 7 | **12** | 0.87 | 12 | 4.7e-16 | own |
| 5 | 9 | 18 | 0.85 | 17 | 9.5e-16 | transcribed: Haegemans and Piessens 1977 |
| 6 | 11 | 25 | 0.83 | 24 | 5.9e-16 | transcribed: Haegemans and Piessens 1976, hexagonal |
| 7 | 13 | 34 | 0.83 | 31 | 7.2e-16 | transcribed: Cools and Haegemans 1988 |
| 8 | 15 | 44 | 0.83 | 40 | 4.3e-16 | own |
| 9 | 17 | 55 | 0.82 | 49 | 9.2e-16 | own |
| 10 | 19 | 68 | 0.82 | 60 | 5.7e-16 | own |
| 11 | 21 | 82 | 0.82 | 71 | 8.6e-16 | own |
| 12 | 23 | 97 | 0.82 | 84 | 7.2e-16 | own |
| 13 | 25 | 114 | 0.82 | 97 | 7.3e-16 | own |
| 14 | 27 | 132 | 0.82 | 112 | 1.2e-15 | own |
| 15 | 29 | 153 | 0.82 | 127 | 1.1e-15 | own |
| 16 | 31 | 178 | 0.83 | 144 | 6.2e-16 | own |
| 17 | 33 | 208 | 0.85 | 161 | 1.1e-15 | own |
| 18 | 35 | 264 | 0.90 | 180 | 1.2e-15 | own |
| 19 | 37 | 357 | 0.99 | 199 | 7.8e-16 | own |
| 20 | 39 | 394 | 0.99 | 220 | 4.1e-12 | own |
| 21 | 41 | 439 | 1.00 | 241 | 7.4e-12 | own |
| 22 | 43 | 482 | 1.00 | 264 | 8.4e-16 | own |

### d = 3

| q | p | n | ρ | Möller | rel. err. | origin |
|---:|---:|---:|---:|---:|---:|:---|
| 1 | 1 | **1** | 1.00 | 1 | 0.0e+00 | own |
| 2 | 3 | **6** | 0.91 | 6 | 2.8e-17 | own |
| 3 | 5 | **13** | 0.78 | 13 | 3.0e-16 | transcribed: Stroud 1971, icosahedral (closed form) |
| 4 | 7 | 27 | 0.75 | 26 | 4.7e-16 | same-rule: Stroud 1971, E_n^{r^2}:7-1 option 1 (identified 2026-09-11 vs Burkardt en_r2_07_1) |
| 5 | 9 | 45 | 0.71 | 43 | 4.1e-16 | same-rule: Konyaev 1977, rule 1 for the weight exp(-rho^2) (identified 2026-09-16 against the closed form printed in Dokl. Akad. Nauk SSSR 233 no. 5, 784-787; same up to rotation, node displacement 9.7e-16) |
| 6 | 11 | 77 | 0.71 | 68 | 8.4e-16 | own |
| 7 | 13 | 128 | 0.72 | 99 | 4.4e-16 | own |
| 8 | 15 | 184 | 0.71 | 140 | 2.1e-15 | own |
| 9 | 17 | 264 | 0.71 | 189 | 6.9e-16 | own |
| 10 | 19 | 354 | 0.71 | 250 | 1.1e-15 | own |
| 11 | 21 | 476 | 0.71 | 321 | 4.7e-16 | own |
| 12 | 23 | 597 | 0.70 | 406 | 5.8e-16 | own |
| 13 | 25 | 776 | 0.71 | 503 | 5.5e-16 | own |
| 14 | 27 | 966 | 0.71 | 616 | 7.4e-16 | own |
| 15 | 29 | 1242 | 0.72 | 743 | 5.6e-16 | own |
| 16 | 31 | 1848 | 0.77 | 888 | 6.6e-16 | own |
| 17 | 33 | 2226 | 0.77 | 1049 | 7.8e-16 | own |
| 18 | 35 | 4749 | 0.93 | 1230 | 7.2e-16 | own |

### d = 4

| q | p | n | ρ | Möller | rel. err. | origin |
|---:|---:|---:|---:|---:|---:|:---|
| 1 | 1 | **1** | 1.00 | 1 | 0.0e+00 | own |
| 2 | 3 | **8** | 0.84 | 8 | 2.2e-16 | own |
| 3 | 5 | 22 | 0.72 | 21 | 4.4e-16 | own |
| 4 | 7 | 49 | 0.66 | 48 | 9.5e-16 | same-rule: Stroud 1971, E_n^{r^2}:7-1 option 1 (identified 2026-09-11 vs Burkardt en_r2_07_1) |
| 5 | 9 | 116 | 0.66 | 91 | 4.7e-16 | own |
| 6 | 11 | 193 | 0.62 | 160 | 1.8e-15 | own |
| 7 | 13 | 414 | 0.64 | 259 | 4.1e-16 | own |
| 8 | 15 | 577 | 0.61 | 400 | 5.0e-15 | own |
| 9 | 17 | 1056 | 0.63 | 589 | 3.9e-16 | own |
| 10 | 19 | 1505 | 0.62 | 840 | 5.7e-16 | own |
| 11 | 21 | 2318 | 0.63 | 1161 | 6.5e-16 | own |
| 12 | 23 | 3238 | 0.63 | 1568 | 7.3e-16 | own |

### d = 5

| q | p | n | ρ | Möller | rel. err. | origin |
|---:|---:|---:|---:|---:|---:|:---|
| 1 | 1 | **1** | 1.00 | 1 | 0.0e+00 | own |
| 2 | 3 | **10** | 0.79 | 10 | 4.4e-16 | same-rule: Stroud 1971, E_n^{r^2}:3-1 (identified 2026-09-11, displacement 0.0; order 2N = 10 at N=5) |
| 3 | 5 | 32 | 0.67 | 31 | 3.0e-16 | transcribed: Stroud and Secrest 1963, E5r2:5-1 |
| 4 | 7 | 83 | 0.61 | 80 | 7.1e-16 | same-rule: Stroud 1971, E_n^{r^2}:7-1 (identified 2026-09-11, compare_rules.jl over B_5) |
| 5 | 9 | 244 | 0.60 | 171 | 3.9e-15 | own |
| 6 | 11 | 485 | 0.57 | 332 | 7.9e-16 | own |
| 7 | 13 | 1135 | 0.58 | 591 | 7.0e-16 | own |
| 8 | 15 | 1767 | 0.56 | 992 | 7.9e-16 | own |
| 9 | 17 | 3986 | 0.58 | 1581 | 1.8e-15 | own |
| 10 | 19 | 7174 | 0.59 | 2422 | 4.7e-15 | own |
| 11 | 21 | 13199 | 0.61 | 3583 | 9.1e-15 | own |

## Le — uniform weight on the cube

### d = 2

| q | p | n | ρ | Möller | rel. err. | origin |
|---:|---:|---:|---:|---:|---:|:---|
| 1 | 1 | **1** | 1.00 | 1 | 0.0e+00 | own |
| 2 | 3 | **4** | 1.00 | 4 | 5.6e-17 | own |
| 3 | 5 | **7** | 0.88 | 7 | 5.6e-17 | own |
| 4 | 7 | **12** | 0.87 | 12 | 2.8e-16 | own |
| 5 | 9 | **17** | 0.82 | 17 | 5.6e-17 | same-rule: Festa and Sommariva 2012, SMR09 (count first published by Moller 1976) |
| 6 | 11 | **24** | 0.82 | 24 | 1.1e-16 | own |
| 7 | 13 | 33 | 0.82 | 31 | 1.1e-16 | same-rule: Festa and Sommariva 2012, SMR13 |
| 8 | 15 | 43 | 0.82 | 40 | 1.1e-16 | same-rule: Festa and Sommariva 2012, SMR15 |
| 9 | 17 | 54 | 0.82 | 49 | 1.1e-16 | own |
| 10 | 19 | 67 | 0.82 | 60 | 6.9e-17 | own |
| 11 | 21 | 81 | 0.82 | 71 | 2.2e-16 | own |
| 12 | 23 | 96 | 0.82 | 84 | 5.6e-17 | own |
| 13 | 25 | 113 | 0.82 | 97 | 2.8e-17 | transcribed: Festa and Sommariva 2012, SMR25 |
| 14 | 27 | 131 | 0.82 | 112 | 2.2e-16 | own |
| 15 | 29 | 150 | 0.82 | 127 | 1.1e-16 | own |
| 16 | 31 | 171 | 0.82 | 144 | 2.2e-16 | own |
| 17 | 33 | 194 | 0.82 | 161 | 2.2e-16 | own |
| 18 | 35 | 216 | 0.82 | 180 | 2.2e-16 | own |
| 19 | 37 | 242 | 0.82 | 199 | 1.1e-16 | own |
| 20 | 39 | 268 | 0.82 | 220 | 1.1e-16 | own |
| 21 | 41 | 294 | 0.82 | 241 | 5.6e-17 | own |
| 22 | 43 | 324 | 0.82 | 264 | 1.1e-16 | own |
| 23 | 45 | 356 | 0.82 | 287 | 2.2e-16 | own |
| 24 | 47 | 388 | 0.82 | 312 | 1.1e-16 | own |
| 25 | 49 | 422 | 0.82 | 337 | 5.6e-17 | own |
| 26 | 51 | 456 | 0.82 | 364 | 5.6e-17 | own |
| 27 | 53 | 492 | 0.82 | 391 | 5.6e-17 | own |
| 28 | 55 | 528 | 0.82 | 420 | 5.6e-17 | own |
| 29 | 57 | 570 | 0.82 | 449 | 5.6e-17 | own |
| 30 | 59 | 606 | 0.82 | 480 | 3.3e-16 | own |
| 31 | 61 | 642 | 0.82 | 511 | 5.6e-17 | own |
| 32 | 63 | 692 | 0.82 | 544 | 1.1e-16 | own |
| 33 | 65 | 732 | 0.82 | 577 | 1.1e-16 | own |
| 34 | 67 | 780 | 0.82 | 612 | 5.6e-17 | own |
| 35 | 69 | 826 | 0.82 | 647 | 1.1e-16 | own |
| 36 | 71 | 872 | 0.82 | 684 | 5.6e-17 | own |
| 37 | 73 | 922 | 0.82 | 721 | 1.1e-16 | own |
| 38 | 75 | 980 | 0.82 | 760 | 2.2e-16 | own |
| 39 | 77 | 1032 | 0.82 | 799 | 5.6e-17 | derived: Diallo and Worku 2026, elimination started from their rule dw_d2_p77_n1049 for the same cell |

### d = 3

| q | p | n | ρ | Möller | rel. err. | origin |
|---:|---:|---:|---:|---:|---:|:---|
| 1 | 1 | **1** | 1.00 | 1 | 0.0e+00 | own |
| 2 | 3 | **6** | 0.91 | 6 | 2.2e-16 | own |
| 3 | 5 | **13** | 0.78 | 13 | 5.6e-17 | own |
| 4 | 7 | **26** | 0.74 | 26 | 1.1e-16 | own |
| 5 | 9 | 48 | 0.73 | 43 | 2.2e-16 | own |
| 6 | 11 | 82 | 0.72 | 68 | 2.2e-16 | own |
| 7 | 13 | 128 | 0.72 | 99 | 2.2e-16 | own |
| 8 | 15 | 188 | 0.72 | 140 | 4.4e-16 | own |
| 9 | 17 | 266 | 0.71 | 189 | 5.6e-17 | own |
| 10 | 19 | 360 | 0.71 | 250 | 3.3e-16 | own |
| 11 | 21 | 476 | 0.71 | 321 | 1.1e-16 | own |
| 12 | 23 | 612 | 0.71 | 406 | 2.2e-16 | own |
| 13 | 25 | 776 | 0.71 | 503 | 4.4e-16 | own |
| 14 | 27 | 964 | 0.71 | 616 | 5.6e-17 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p27_n984 for the same cell |
| 15 | 29 | 1184 | 0.71 | 743 | 1.1e-16 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p29_n1258 for the same cell |
| 16 | 31 | 1432 | 0.70 | 888 | 2.2e-16 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p31_n1478 for the same cell |
| 17 | 33 | 1714 | 0.70 | 1049 | 2.2e-16 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p33_n1787 for the same cell |
| 18 | 35 | 2028 | 0.70 | 1230 | 3.3e-16 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p35_n2102 for the same cell |
| 19 | 37 | 2380 | 0.70 | 1429 | 4.4e-16 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p37_n2506 for the same cell |
| 20 | 39 | 2770 | 0.70 | 1650 | 2.2e-16 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p39_n2856 for the same cell |
| 21 | 41 | 3200 | 0.70 | 1891 | 2.2e-16 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p41_n3338 for the same cell |
| 22 | 43 | 3704 | 0.70 | 2156 | 1.7e-15 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p43_n3870 for the same cell |
| 23 | 45 | 4308 | 0.71 | 2443 | 1.1e-16 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p45_n4414 for the same cell |

### d = 4

| q | p | n | ρ | Möller | rel. err. | origin |
|---:|---:|---:|---:|---:|---:|:---|
| 1 | 1 | **1** | 1.00 | 1 | 0.0e+00 | own |
| 2 | 3 | **8** | 0.84 | 8 | 5.6e-17 | own |
| 3 | 5 | **21** | 0.71 | 21 | 1.1e-16 | own |
| 4 | 7 | 54 | 0.68 | 48 | 5.6e-17 | own |
| 5 | 9 | 120 | 0.66 | 91 | 5.6e-17 | own |
| 6 | 11 | 234 | 0.65 | 160 | 4.4e-16 | own |
| 7 | 13 | 416 | 0.65 | 259 | 4.4e-16 | own |
| 8 | 15 | 690 | 0.64 | 400 | 1.1e-16 | own |
| 9 | 17 | 1078 | 0.64 | 589 | 2.2e-16 | own |
| 10 | 19 | 1612 | 0.63 | 840 | 3.3e-16 | own |
| 11 | 21 | 2322 | 0.63 | 1161 | 1.1e-16 | own |
| 12 | 23 | 3244 | 0.63 | 1568 | 1.1e-16 | own |

### d = 5

| q | p | n | ρ | Möller | rel. err. | origin |
|---:|---:|---:|---:|---:|---:|:---|
| 1 | 1 | **1** | 1.00 | 1 | 0.0e+00 | own |
| 2 | 3 | **10** | 0.79 | 10 | 2.2e-16 | own |
| 3 | 5 | 32 | 0.67 | 31 | 5.6e-17 | own |
| 4 | 7 | 100 | 0.63 | 80 | 2.2e-16 | own |
| 5 | 9 | 266 | 0.61 | 171 | 5.6e-17 | own |
| 6 | 11 | 602 | 0.60 | 332 | 5.6e-16 | own |
| 7 | 13 | 1212 | 0.59 | 591 | 2.2e-16 | own |
| 8 | 15 | 2500 | 0.60 | 992 | 4.4e-16 | own |
| 9 | 17 | 3872 | 0.58 | 1581 | 4.4e-16 | own |
| 10 | 19 | 6826 | 0.58 | 2422 | 4.9e-14 | own |
| 11 | 21 | 10984 | 0.58 | 3583 | 3.3e-16 | own |

