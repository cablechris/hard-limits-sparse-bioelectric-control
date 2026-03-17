# Project: Minimum-Energy Bioelectric Control

## What This Is

This project asks a control question that builds directly on the adaptive gap-junction analysis completed in the parent repository. Instead of asking whether the system admits Lyapunov structure or supports an inference interpretation, this project asks how to steer the system toward a desired morphology with the smallest possible intervention.

## Core Research Question

For a nonlinear gap-junction tissue network with noise, disorder, and lesion-induced damage, what is the minimum-energy spatiotemporal control input that drives the system from a damaged initial state into the basin of a target healthy morphology, and how does that threshold scale with network size, conductance, noise, and dilution?

## Prior Work To Carry Forward

- Parent repository: `C:\Users\cable\get-physics-done`
- Key prior result 1: symmetric diffusive coupling admits a Lyapunov-like structure in the fixed-conductance sector
- Key prior result 2: that structure favors homogenization rather than target morphology
- Key prior result 3: AP pattern maintenance requires explicit non-equilibrium drive
- Key prior result 4: field-off persistence over short windows is consistent with slow relaxation, not yet proven metastability
- Key prior result 5: exact inference correspondence is dimensionally obstructed for the current coarse surrogate

## Deliverables

1. Minimum intervention threshold as a function of network size, noise, and disorder.
2. Optimal stimulation pattern and dominant steering modes.
3. Dose-response curve for recovery probability versus intervention energy.
4. Disorder robustness bound under bond dilution.
5. Manuscript-ready derivations, figures, and reproducible code.

## Constraints

- Keep the control objective explicit. Do not blur basin steering with vague "repair" language.
- Default control geometry is additive site-local voltage forcing: `dV/dt = F(V,G) + B u(t)`.
- The default actuator family is sparse site control; full-site control is retained only as a lower-bound benchmark.
- Treat quenched bond dilution as the disorder model when comparing against octanol-style blockade.
- Distinguish exact linearized control results from nonlinear numerical steering results.
- Do not claim metastable control memory unless the post-control basin structure is demonstrated.
- Carry forward the parent repo's mathematical-status labeling discipline.

## Key References

- Parent repo discussion: `C:\Users\cable\get-physics-done\research\adaptive_gap_junction\DISCUSSION.md`
- Parent Lyapunov note: `C:\Users\cable\get-physics-done\docs\derivations\lyapunov-obstruction.md`
- Parent correspondence note: `C:\Users\cable\get-physics-done\docs\derivations\inference-correspondence.md`
