# Roadmap: Minimum-Energy Bioelectric Control

## Phases

- [ ] **Phase 1: Control-Ready Model Specification** - Fix the controlled dynamics, actuator operator, intervention cost, and linearization target.
- [ ] **Phase 2: Linearized Minimum-Energy Control** - Derive the linear steering solution, controllability condition, and Gramian-based energy threshold.
- [ ] **Phase 3: Nonlinear Basin Steering** - Compute minimum-energy controls for the nonlinear network and test basin entry under finite-amplitude interventions.
- [ ] **Phase 4: Noise and Disorder Robustness** - Quantify how control thresholds and recovery probability scale with stochastic forcing and bond dilution.
- [ ] **Phase 5: Threshold Laws and Experimental Predictions** - Extract scaling laws, optimal stimulation modes, and dose-response predictions.
- [ ] **Phase 6: Manuscript and Artifact Packaging** - Prepare derivations, figures, claim-status table, and reproducibility bundle.

## Phase 1 Details

**Goal:** Turn the inherited adaptive gap-junction model into a control problem with explicit state equation, actuator operator, and cost definitions.

Plans:

- [x] 01-01: Define the controlled dynamics `dV/dt = F(V,G) + B u(t)` and intervention-energy functional.
- [x] 01-02: Define the default linearization target as a damaged pre-intervention snapshot `(V_bar, G_bar)` with frozen conductance field, and record the resulting `(A, B)` pair for Phase 2.
- [x] 01-03: Build the control-ready code scaffold and verification checks.
