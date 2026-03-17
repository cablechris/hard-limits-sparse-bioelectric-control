# Default Lattice Phase 2 Baseline

- Horizon: `4.0`
- State dimension: `256`
- Full overlap artifact: `docs\results\default-lattice-a-vs-l-overlap.csv`

## Informed Sparse Placement

- Nodes: `(135, 136)`
- Reachable rank at `1e-12`: `15`
- Reachable-subspace condition number: `2.01717e+11`
- Smallest positive Gramian eigenvalue: `1.52486e-12`
- Largest Gramian eigenvalue: `0.307589`

## Random Two-Site Sweep

- Samples: `32`
- Median reachable rank at `1e-12`: `14`
- Median condition number: `6.08812e+10`
- Best random placement: nodes `(50, 158)`, condition `1.49777e+10`, rank `12`
- Worst random placement: nodes `(207, 248)`, condition `2.78808e+11`, rank `14`

## Interpretation

The full overlap matrix records how far the frozen voltage Jacobian eigenbasis departs from a pure network-Laplacian basis.
The placement sweep records whether two-site stimulation geometry changes the conditioning of the reachable subspace enough to matter experimentally.
