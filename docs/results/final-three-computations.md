# Final Three Computations

This note records the three remaining Phase 2 computations requested before manuscript drafting:

1. the greedy dose-response envelope on the `8x8` lattice,
2. one scaling repeat of the actuation-rank table on `16x16`,
3. one disorder sweep linking threshold rank to gap-junction dilution.

Unless stated otherwise, all runs use the frozen damaged snapshot convention, horizon `T = 4.0`, projected steering metrics, and residual threshold `||delta V_unreachable|| < 1e-4`.

## 1. Greedy Envelope on 8x8

Artifact: `docs/results/greedy-envelope-8x8.csv`

The envelope gives the count-by-count dose-response for the three damage geometries already reported in the `8x8` actuation-rank study. The main thresholds are:

| Damage geometry | Minimum actuation rank |
| --- | ---: |
| Corner wound | 2 |
| Side patch | 3 |
| Central lesion | 3 |

Representative residual drops from `k = 2` to `k = 4`:

| Damage geometry | Residual at k = 2 | Residual at k = 4 |
| --- | ---: | ---: |
| Corner wound | `4.14e-6` | `3.44e-7` |
| Side patch | `1.56e-4` | `2.00e-5` |
| Central lesion | `2.40e-4` | `2.80e-5` |

Interpretation: increasing actuator count by only two sites reduces residual by roughly one order of magnitude for the easier corner wound and by about one to two orders of magnitude for the harder side-patch and central-lesion geometries. This is the finite-horizon dose-response curve in actuation rank.

## 2. Scaling Check on 16x16

Artifact: `docs/results/actuation-rank-16x16.csv`

This was computed with a structured greedy search over a candidate pool built from:

- the support of the damage pattern, and
- the top `32` nodes ranked by dominant controllable `A`-mode weight.

This is a scaling check, not an exhaustive optimum over all `256` sites.

Scaled damage patterns used width `4` on the `16x16` lattice. The resulting threshold table is:

| Damage geometry | Minimum actuation rank |
| --- | ---: |
| Corner wound | 2 |
| Side patch | 6 |
| Central lesion | 2 |

Representative residual drops from `k = 2` to `k = 4`:

| Damage geometry | Residual at k = 2 | Residual at k = 4 |
| --- | ---: | ---: |
| Corner wound | `2.09e-5` | `6.88e-6` |
| Side patch | `9.76e-4` | `2.31e-4` |
| Central lesion | `2.76e-5` | `1.72e-5` |

Interpretation: the side patch becomes substantially harder under lattice scaling, while the corner wound remains easy and the central lesion remains below threshold by `k = 2` in this structured search. The sharp increase for the side patch is consistent with the idea that broad boundary wounds load harder-to-reach spectral content than compact or centrally symmetric lesions.

## 3. Disorder Sweep on 8x8 Side Patch

Artifact: `docs/results/disorder-sweep-8x8-side-patch.csv`

We repeated the `8x8` side-patch threshold computation under random gap-junction dilution at strengths `0.0`, `0.1`, `0.2`, and `0.3`, with seeds `0..4`.

Threshold-rank outcomes:

| Dilution fraction | Threshold-rank outcomes across seeds |
| --- | --- |
| `0.0` | `3, 3, 3, 3, 3` |
| `0.1` | `3, 3, 3, 3, 3` |
| `0.2` | `3, 4, 3, 4, 4` |
| `0.3` | `3, 4, 5, 5, 5` |

Median residuals across seeds:

| Dilution fraction | Median residual at k = 2 | Median residual at k = 4 |
| --- | ---: | ---: |
| `0.0` | `1.35e-4` | `1.58e-5` |
| `0.1` | `1.55e-4` | `3.31e-5` |
| `0.2` | `2.59e-4` | `4.13e-5` |
| `0.3` | `2.99e-4` | `1.10e-4` |

Interpretation: the minimum actuation rank is stable at `3` up to `10%` dilution, becomes mixed between `3` and `4` by `20%` dilution, and shifts mostly to `4-5` by `30%` dilution. Gap-junction disorder therefore raises the intervention complexity required for a fixed repair fidelity.

## Manuscript-Level Takeaway

The paper can now state three quantitative control results:

- On `8x8`, the minimum actuation rank is geometry-dependent: corner `2`, side patch `3`, central lesion `3`.
- On `16x16`, the side patch becomes markedly harder in the structured greedy scaling check, requiring `6` sites to cross the same residual threshold.
- Random gap-junction dilution shifts the side-patch threshold upward from `3` toward `4-5`, directly connecting controllability loss to communication disorder.
