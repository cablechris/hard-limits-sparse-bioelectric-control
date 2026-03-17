# Placement Study Phase 2

## 8x8 Exhaustive Sweep

- Damage pattern: low-x depolarization patch, width `2`, amplitude `1.0`
- Horizon: `4.0`
- Pairs scanned: `2016`
- Exactly reachable pairs at `1e-12`: `0`
- Exact result: no two-site pair exactly steers this wound pattern to zero on the chosen horizon.
- Best projected pair: `(8, 41)` with projected residual `0.000137294` and projected energy `1790.89`
- Median projected residual: `0.00027441`
- Median projected energy: `9939.06`
- Rule candidate from top controllable A-modes: `(27, 35)` with projected residual `0.000260748` and projected energy `305.112`
- Pair scan artifact: `docs\results\placement-scan-8x8-side-patch.csv`
- Projected-energy heatmap: `docs\results\placement-heatmap-8x8-side-patch.csv`
- Residual heatmap: `docs\results\placement-residual-8x8-side-patch.csv`

## 16x16 Structured Comparison

- default_midline: nodes `(135, 136)`, exact energy `inf`, projected residual `0.00043799`, projected energy `146.133`, reachable rank `15`
- rule_based: nodes `(120, 136)`, exact energy `inf`, projected residual `0.000438243`, projected energy `425.773`, reachable rank `16`
- corner_pair: nodes `(0, 15)`, exact energy `inf`, projected residual `0.000367894`, projected energy `3047.56`, reachable rank `12`
