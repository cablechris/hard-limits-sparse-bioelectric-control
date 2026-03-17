# Actuation Rank Study 8x8

- Horizon: `4.0`
- Modal rule candidate pair from dominant controllable A-modes: `(27, 35)`

## side_patch

- Best two-site pair: `(8, 41)`
- Best two-site projected residual / energy: `0.000137294` / `1790.89`
- Modal-rule pair projected residual / energy: `0.000260748` / `305.112`
- Greedy 4-site extension: `(8, 41, 42, 43)` with residual / energy `9.40341e-06` / `679.129`
- Greedy 8-site extension: `(8, 11, 12, 32, 40, 41, 42, 43)` with residual / energy `9.75446e-07` / `13.8887`
- Greedy 16-site extension: `(7, 8, 11, 12, 14, 16, 19, 27, 30, 32, 40, 41, 42, 43, 58, 59)` with residual / energy `0` / `0.0147135`
- Residual threshold hits from best-two seed: `1e-3 -> 2`, `5e-4 -> 2`, `1e-4 -> 3`
- Trajectory artifact: `docs/results/greedy-rank-trajectory-side_patch-8x8.csv`

## corner_wound

- Best two-site pair: `(0, 9)`
- Best two-site projected residual / energy: `2.90001e-06` / `23.77`
- Modal-rule pair projected residual / energy: `0.000111688` / `102.083`
- Greedy 4-site extension: `(0, 2, 5, 9)` with residual / energy `3.76966e-07` / `1.73461`
- Greedy 8-site extension: `(0, 2, 5, 9, 17, 23, 33, 40)` with residual / energy `2.73108e-08` / `0.135255`
- Greedy 16-site extension: `(0, 2, 3, 5, 8, 9, 11, 17, 18, 20, 23, 28, 33, 40, 53, 60)` with residual / energy `0` / `0.000344406`
- Residual threshold hits from best-two seed: `1e-3 -> 2`, `5e-4 -> 2`, `1e-4 -> 2`
- Trajectory artifact: `docs/results/greedy-rank-trajectory-corner_wound-8x8.csv`

## central_lesion

- Best two-site pair: `(11, 20)`
- Best two-site projected residual / energy: `0.000189178` / `44511.9`
- Modal-rule pair projected residual / energy: `0.000716483` / `5.14464`
- Greedy 4-site extension: `(11, 20, 28, 35)` with residual / energy `6.4411e-06` / `25.9266`
- Greedy 8-site extension: `(11, 20, 28, 29, 35, 36, 39, 56)` with residual / energy `1.37295e-07` / `9.45486`
- Greedy 16-site extension: `(4, 11, 14, 20, 24, 27, 28, 29, 34, 35, 36, 39, 42, 46, 52, 56)` with residual / energy `0` / `0.136492`
- Residual threshold hits from best-two seed: `1e-3 -> 2`, `5e-4 -> 2`, `1e-4 -> 3`
- Trajectory artifact: `docs/results/greedy-rank-trajectory-central_lesion-8x8.csv`

