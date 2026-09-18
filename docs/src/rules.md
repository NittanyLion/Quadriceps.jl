# Stored rules

Written by `build/build_data.jl` on 2026-09-18. `n`: number of nodes; `ρ = n^(1/d) / q` with
`q = (p+1)/2`: the node count relative to the `q^d` Gauss product grid (1.00 is that grid, smaller
is better); Möller: Möller's lower bound on `n` (**bold** `n`: bound attained, proven minimal);
rel. err.: largest relative monomial error of the stored `Float64` rule; origin: `own`, or the
published rule it is, or descends from (see [Credits](credits.md)).

## GH — Gaussian weight

### d = 2

| p | n | ρ | Möller | rel. err. | origin |
|---:|---:|---:|---:|---:|:---|
| 1 | **1** | 1.00 | 1 | 0.0e+00 | own |
| 3 | **4** | 1.00 | 4 | 2.2e-16 | own |
| 5 | **7** | 0.88 | 7 | 2.2e-16 | same-rule: Stroud 1971, E_2^{r^2}:5-4 (identified 2026-09-11; same up to rotation, which the Gaussian weight permits) |
| 7 | **12** | 0.87 | 12 | 4.7e-16 | own |
| 9 | 18 | 0.85 | 17 | 9.5e-16 | transcribed: Haegemans and Piessens 1977 |
| 11 | 25 | 0.83 | 24 | 5.9e-16 | transcribed: Haegemans and Piessens 1976, hexagonal |
| 13 | 34 | 0.83 | 31 | 7.2e-16 | transcribed: Cools and Haegemans 1988 |
| 15 | 44 | 0.83 | 40 | 4.3e-16 | own |
| 17 | 55 | 0.82 | 49 | 9.2e-16 | own |
| 19 | 68 | 0.82 | 60 | 5.7e-16 | own |
| 21 | 82 | 0.82 | 71 | 8.6e-16 | own |
| 23 | 97 | 0.82 | 84 | 7.2e-16 | own |
| 25 | 114 | 0.82 | 97 | 7.3e-16 | own |
| 27 | 132 | 0.82 | 112 | 1.2e-15 | own |
| 29 | 153 | 0.82 | 127 | 1.1e-15 | own |
| 31 | 178 | 0.83 | 144 | 6.2e-16 | own |
| 33 | 208 | 0.85 | 161 | 8.0e-13 | own |
| 35 | 264 | 0.90 | 180 | 1.2e-15 | own |
| 37 | 357 | 0.99 | 199 | 1.9e-12 | own |
| 39 | 394 | 0.99 | 220 | 4.1e-12 | own |
| 41 | 439 | 1.00 | 241 | 7.4e-12 | own |
| 43 | 482 | 1.00 | 264 | 8.5e-12 | own |

### d = 3

| p | n | ρ | Möller | rel. err. | origin |
|---:|---:|---:|---:|---:|:---|
| 1 | **1** | 1.00 | 1 | 0.0e+00 | own |
| 3 | **6** | 0.91 | 6 | 2.8e-17 | own |
| 5 | **13** | 0.78 | 13 | 3.0e-16 | transcribed: Stroud 1971, icosahedral (closed form) |
| 7 | 27 | 0.75 | 26 | 4.7e-16 | same-rule: Stroud 1971, E_n^{r^2}:7-1 option 1 (identified 2026-09-11 vs Burkardt en_r2_07_1) |
| 9 | 45 | 0.71 | 43 | 4.1e-16 | same-rule: Konyaev 1977, rule 1 for the weight exp(-rho^2) (identified 2026-09-16 against the closed form printed in Dokl. Akad. Nauk SSSR 233 no. 5, 784-787; same up to rotation, node displacement 9.7e-16) |
| 11 | 77 | 0.71 | 68 | 8.4e-16 | own |
| 13 | 128 | 0.72 | 99 | 4.4e-16 | own |
| 15 | 184 | 0.71 | 140 | 2.1e-15 | own |
| 17 | 264 | 0.71 | 189 | 6.9e-16 | own |
| 19 | 354 | 0.71 | 250 | 1.1e-15 | own |
| 21 | 476 | 0.71 | 321 | 4.7e-16 | own |
| 23 | 597 | 0.70 | 406 | 5.8e-16 | own |
| 25 | 776 | 0.71 | 503 | 5.5e-16 | own |
| 27 | 966 | 0.71 | 616 | 7.4e-16 | own |
| 29 | 1242 | 0.72 | 743 | 5.6e-16 | own |
| 31 | 1848 | 0.77 | 888 | 1.2e-12 | own |
| 33 | 2226 | 0.77 | 1049 | 3.4e-12 | own |
| 35 | 4749 | 0.93 | 1230 | 3.9e-12 | own |

### d = 4

| p | n | ρ | Möller | rel. err. | origin |
|---:|---:|---:|---:|---:|:---|
| 1 | **1** | 1.00 | 1 | 0.0e+00 | own |
| 3 | **8** | 0.84 | 8 | 2.2e-16 | own |
| 5 | 22 | 0.72 | 21 | 4.4e-16 | own |
| 7 | 49 | 0.66 | 48 | 9.5e-16 | same-rule: Stroud 1971, E_n^{r^2}:7-1 option 1 (identified 2026-09-11 vs Burkardt en_r2_07_1) |
| 9 | 116 | 0.66 | 91 | 4.7e-16 | own |
| 11 | 193 | 0.62 | 160 | 1.8e-15 | own |
| 13 | 414 | 0.64 | 259 | 4.1e-16 | own |
| 15 | 577 | 0.61 | 400 | 5.0e-15 | own |
| 17 | 1056 | 0.63 | 589 | 3.9e-16 | own |
| 19 | 1505 | 0.62 | 840 | 5.7e-16 | own |
| 21 | 2318 | 0.63 | 1161 | 6.5e-16 | own |
| 23 | 3238 | 0.63 | 1568 | 7.3e-16 | own |

### d = 5

| p | n | ρ | Möller | rel. err. | origin |
|---:|---:|---:|---:|---:|:---|
| 1 | **1** | 1.00 | 1 | 0.0e+00 | own |
| 3 | **10** | 0.79 | 10 | 4.4e-16 | same-rule: Stroud 1971, E_n^{r^2}:3-1 (identified 2026-09-11, displacement 0.0; order 2N = 10 at N=5) |
| 5 | 32 | 0.67 | 31 | 3.0e-16 | transcribed: Stroud and Secrest 1963, E5r2:5-1 |
| 7 | 83 | 0.61 | 80 | 7.1e-16 | same-rule: Stroud 1971, E_n^{r^2}:7-1 (identified 2026-09-11, compare_rules.jl over B_5) |
| 9 | 244 | 0.60 | 171 | 3.9e-15 | own |
| 11 | 485 | 0.57 | 332 | 7.9e-16 | own |
| 13 | 1135 | 0.58 | 591 | 7.0e-16 | own |
| 15 | 1767 | 0.56 | 992 | 7.9e-16 | own |
| 17 | 3986 | 0.58 | 1581 | 1.8e-15 | own |
| 19 | 7174 | 0.59 | 2422 | 4.7e-15 | own |
| 21 | 13199 | 0.61 | 3583 | 9.1e-15 | own |

## Le — uniform weight on the cube

### d = 2

| p | n | ρ | Möller | rel. err. | origin |
|---:|---:|---:|---:|---:|:---|
| 1 | **1** | 1.00 | 1 | 0.0e+00 | own |
| 3 | **4** | 1.00 | 4 | 5.6e-17 | own |
| 5 | **7** | 0.88 | 7 | 5.6e-17 | own |
| 7 | **12** | 0.87 | 12 | 2.8e-16 | own |
| 9 | **17** | 0.82 | 17 | 5.6e-17 | same-rule: Festa and Sommariva 2012, SMR09 (count first published by Moller 1976) |
| 11 | **24** | 0.82 | 24 | 1.1e-16 | own |
| 13 | 33 | 0.82 | 31 | 1.1e-16 | same-rule: Festa and Sommariva 2012, SMR13 |
| 15 | 43 | 0.82 | 40 | 1.1e-16 | same-rule: Festa and Sommariva 2012, SMR15 |
| 17 | 54 | 0.82 | 49 | 1.1e-16 | own |
| 19 | 67 | 0.82 | 60 | 6.9e-17 | own |
| 21 | 81 | 0.82 | 71 | 2.2e-16 | own |
| 23 | 96 | 0.82 | 84 | 5.6e-17 | own |
| 25 | 113 | 0.82 | 97 | 2.8e-17 | transcribed: Festa and Sommariva 2012, SMR25 |
| 27 | 131 | 0.82 | 112 | 2.2e-16 | own |
| 29 | 150 | 0.82 | 127 | 1.1e-16 | own |
| 31 | 171 | 0.82 | 144 | 2.2e-16 | own |
| 33 | 194 | 0.82 | 161 | 2.2e-16 | own |
| 35 | 216 | 0.82 | 180 | 2.2e-16 | own |
| 37 | 242 | 0.82 | 199 | 1.1e-16 | own |
| 39 | 268 | 0.82 | 220 | 1.1e-16 | own |
| 41 | 294 | 0.82 | 241 | 5.6e-17 | own |
| 43 | 324 | 0.82 | 264 | 1.1e-16 | own |
| 45 | 356 | 0.82 | 287 | 2.2e-16 | own |
| 47 | 388 | 0.82 | 312 | 1.1e-16 | own |
| 49 | 422 | 0.82 | 337 | 5.6e-17 | own |
| 51 | 456 | 0.82 | 364 | 5.6e-17 | own |
| 53 | 492 | 0.82 | 391 | 5.6e-17 | own |
| 55 | 528 | 0.82 | 420 | 5.6e-17 | own |
| 57 | 570 | 0.82 | 449 | 5.6e-17 | own |
| 59 | 606 | 0.82 | 480 | 3.3e-16 | own |
| 61 | 642 | 0.82 | 511 | 5.6e-17 | own |
| 63 | 692 | 0.82 | 544 | 1.1e-16 | own |
| 65 | 732 | 0.82 | 577 | 1.1e-16 | own |
| 67 | 780 | 0.82 | 612 | 5.6e-17 | own |
| 69 | 826 | 0.82 | 647 | 1.1e-16 | own |
| 71 | 872 | 0.82 | 684 | 5.6e-17 | own |
| 73 | 922 | 0.82 | 721 | 1.1e-16 | own |
| 75 | 980 | 0.82 | 760 | 2.2e-16 | own |
| 77 | 1032 | 0.82 | 799 | 5.6e-17 | derived: Diallo and Worku 2026, elimination started from their rule dw_d2_p77_n1049 for the same cell |

### d = 3

| p | n | ρ | Möller | rel. err. | origin |
|---:|---:|---:|---:|---:|:---|
| 1 | **1** | 1.00 | 1 | 0.0e+00 | own |
| 3 | **6** | 0.91 | 6 | 2.2e-16 | own |
| 5 | **13** | 0.78 | 13 | 5.6e-17 | own |
| 7 | **26** | 0.74 | 26 | 1.1e-16 | own |
| 9 | 48 | 0.73 | 43 | 2.2e-16 | own |
| 11 | 82 | 0.72 | 68 | 2.2e-16 | own |
| 13 | 128 | 0.72 | 99 | 2.2e-16 | own |
| 15 | 188 | 0.72 | 140 | 4.4e-16 | own |
| 17 | 266 | 0.71 | 189 | 5.6e-17 | own |
| 19 | 360 | 0.71 | 250 | 3.3e-16 | own |
| 21 | 476 | 0.71 | 321 | 1.1e-16 | own |
| 23 | 612 | 0.71 | 406 | 2.2e-16 | own |
| 25 | 776 | 0.71 | 503 | 4.4e-16 | own |
| 27 | 964 | 0.71 | 616 | 5.6e-17 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p27_n984 for the same cell |
| 29 | 1184 | 0.71 | 743 | 1.1e-16 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p29_n1258 for the same cell |
| 31 | 1432 | 0.70 | 888 | 2.2e-16 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p31_n1478 for the same cell |
| 33 | 1714 | 0.70 | 1049 | 2.2e-16 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p33_n1787 for the same cell |
| 35 | 2028 | 0.70 | 1230 | 6.1e-14 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p35_n2102 for the same cell |
| 37 | 2380 | 0.70 | 1429 | 6.7e-14 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p37_n2506 for the same cell |
| 39 | 2770 | 0.70 | 1650 | 3.3e-13 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p39_n2856 for the same cell |
| 41 | 3200 | 0.70 | 1891 | 1.4e-13 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p41_n3338 for the same cell |
| 43 | 3704 | 0.70 | 2156 | 1.7e-15 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p43_n3870 for the same cell |
| 45 | 4308 | 0.71 | 2443 | 1.0e-15 | derived: Diallo and Worku 2026, elimination started from their rule dw_d3_p45_n4414 for the same cell |

### d = 4

| p | n | ρ | Möller | rel. err. | origin |
|---:|---:|---:|---:|---:|:---|
| 1 | **1** | 1.00 | 1 | 0.0e+00 | own |
| 3 | **8** | 0.84 | 8 | 5.6e-17 | own |
| 5 | **21** | 0.71 | 21 | 1.1e-16 | own |
| 7 | 54 | 0.68 | 48 | 5.6e-17 | own |
| 9 | 120 | 0.66 | 91 | 5.6e-17 | own |
| 11 | 234 | 0.65 | 160 | 4.4e-16 | own |
| 13 | 416 | 0.65 | 259 | 4.4e-16 | own |
| 15 | 690 | 0.64 | 400 | 1.1e-16 | own |
| 17 | 1078 | 0.64 | 589 | 2.2e-16 | own |
| 19 | 1612 | 0.63 | 840 | 3.3e-16 | own |
| 21 | 2322 | 0.63 | 1161 | 6.7e-16 | own |
| 23 | 3244 | 0.63 | 1568 | 1.5e-13 | own |

### d = 5

| p | n | ρ | Möller | rel. err. | origin |
|---:|---:|---:|---:|---:|:---|
| 1 | **1** | 1.00 | 1 | 0.0e+00 | own |
| 3 | **10** | 0.79 | 10 | 2.2e-16 | own |
| 5 | 32 | 0.67 | 31 | 5.6e-17 | own |
| 7 | 100 | 0.63 | 80 | 2.2e-16 | own |
| 9 | 266 | 0.61 | 171 | 5.6e-17 | own |
| 11 | 602 | 0.60 | 332 | 5.6e-16 | own |
| 13 | 1212 | 0.59 | 591 | 2.2e-16 | own |
| 15 | 2506 | 0.60 | 992 | 3.3e-16 | own |
| 17 | 3936 | 0.58 | 1581 | 4.4e-16 | own |
| 19 | 6826 | 0.58 | 2422 | 4.9e-14 | own |
| 21 | 10984 | 0.58 | 3583 | 4.6e-13 | own |

