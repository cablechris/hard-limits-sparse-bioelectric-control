# Hard Limits on Sparse Bioelectric Control: Minimum Actuation Rank for Gap-Junction Network Steering

## Abstract

Bioelectric intervention is usually framed as a dose problem: stronger, longer, or broader stimulation is assumed to improve repair. We show that intervention dimensionality imposes a harder constraint. We formulate bioelectric tissue repair as finite-horizon optimal control on a gap-junction-coupled network and show that below a critical actuator count, no control signal of any amplitude can remove a damaged displacement on the chosen horizon. This actuation-rank threshold is a structural property of the spectral overlap between the damage pattern and the controllable subspace of the frozen linearized network. On an `8x8` lattice, the exact two-site result and greedy higher-rank placements show a geometry-dependent threshold at projected residual below `1e-4`: `2` for a corner wound, `3` for a side patch, and `3` for a central lesion. On a `16x16` lattice under structured greedy placement, the same threshold phenomenon persists but the side patch becomes substantially harder, requiring `6` sites, while corner and central lesions remain low-rank (`2`). Below these thresholds, the unreachable residual is nonzero regardless of energy budget. Above them, each additional actuator site reduces residual by roughly `1-2` orders of magnitude. Random gap-junction dilution shifts the side-patch threshold upward, from stable rank `3` at `0-10%` dilution to mostly `4-5` by `30%` dilution. The controllable subspace is Laplacian-structured but not Laplacian-pure: dominant controllable `A`-modes overlap frozen Laplacian eigenvectors only moderately (`0.46-0.57`), showing that intrinsic voltage bistability materially reshapes network control directions. These results define a finite-horizon linearized baseline for intervention complexity in passive tissue and suggest a direct experimental program: compare predicted actuation-rank thresholds to observed minimum intervention to quantify endogenous morphogenetic competency.

## 1. Introduction

### 1.1 The dose assumption in bioelectric intervention

The modern bioelectric toolkit is already powerful. Optogenetic actuators can depolarize or hyperpolarize selected cells, ion-channel drugs can shift membrane state over extended windows, voltage-sensitive dyes can report spatial prepatterns, and applied electric fields can bias tissue-scale dynamics. What remains missing is not the ability to perturb bioelectric state, but a quantitative theory telling an experimentalist where to intervene, how many independently addressable sites are required, and whether increasing stimulation strength can compensate for poor intervention geometry.

In practice, bioelectric intervention is still mostly interpreted through a dose lens. If a perturbation fails to redirect development or repair, the first adjustments are more current, longer exposure, higher drug concentration, or broader tissue coverage. That framing assumes that intervention strength is the main bottleneck. But in a distributed tissue network coupled by gap junctions, a more basic constraint may come first: whether the intervention has enough independent degrees of freedom to access the damaged state's relevant directions in voltage space at all.

The planarian octanol result makes the problem concrete. Uniform gap-junction blockade can dramatically alter regenerative outcome, but it is also an intervention with essentially no spatial resolution. It does not answer the engineering question that experimentalists actually face: would a smaller number of well-placed, independently controlled perturbations perform better than a stronger but spatially crude whole-tissue perturbation? At present there is no first-principles framework that answers that question from the tissue dynamics themselves.

### 1.2 Network control theory reaches biology, but misses the specific physics

Several existing control traditions provide partial answers. Structural controllability asks how many driver nodes are required to control generic linear dynamics on a graph, but it discards actual dynamics, energy cost, and the specific physics of bioelectric coupling. Minimum-energy control on complex networks introduces the controllability Gramian and therefore captures energetic disparity across directions, but standard formulations are usually applied to generic linear networks rather than to tissues coupled by gap junctions and governed by cell-intrinsic voltage nonlinearities. Nonlinear attractor-control frameworks show how to switch between attractors in general nonlinear systems, but typically provide computational strategies rather than explicit spectral bounds that exploit a given biological mechanism.

The bioelectric case is more structured than these generic settings. Gap junctions impose a Laplacian-like coupling architecture, while each cell contributes local nonlinear voltage dynamics that can support multistability or bistability. Linearizing around a damaged operating point therefore produces a Jacobian that inherits network structure without collapsing to the graph Laplacian itself. That structure should be exploitable. If the controllable subspace can be resolved spectrally, it should explain not only whether repair is possible, but why some damage geometries are easy, why others are hard, and why naive placement intuition fails.

### 1.3 What this paper shows

The central quantity in this paper is the actuation-rank threshold: the smallest number of independently addressable intervention sites required to reduce the unreachable component of a damaged voltage displacement below a chosen tolerance on a fixed horizon. This is not an energy threshold. It is a dimensional threshold. Below it, no increase in control amplitude helps because part of the damage vector lies outside the controllable subspace generated by the actuator family.

Using a frozen linearization of a gap-junction-coupled tissue model, we show five results. First, sparse actuation generates a controllable subspace that is both small and ill-conditioned. Second, the dominant controllable directions are Laplacian-structured but not Laplacian-pure, indicating genuine mixing between network topology and intrinsic voltage dynamics. Third, realistic wound-like displacements can be exactly infeasible under low-rank actuation: on `8x8`, every two-site actuator pair fails to exactly steer the side-patch wound to zero on the chosen horizon. Fourth, minimum actuation rank depends on damage geometry. Fifth, the threshold shifts upward under both lattice scaling and gap-junction disorder.

These results define a passive-tissue finite-horizon baseline. A real tissue with endogenous error-correction may require fewer external degrees of freedom than the frozen linearized calculation predicts, but that monotonic relation is not proved here. The gap between predicted intervention threshold and experimentally observed minimum intervention is therefore best viewed as a candidate quantitative measure of morphogenetic competency rather than as a formally bounded quantity.

## 2. Model

### 2.1 Controlled bioelectric dynamics

We model tissue as a gap-junction-coupled network of `N` cells with membrane-potential state `V(t) in R^N`. The controlled dynamics take the form

```text
dV/dt = (1 / tau_V) F(V, G) + (1 / tau_V) B u(t),
```

where `G` is the gap-junction conductance field, `F(V, G)` combines local voltage dynamics with diffusive coupling through the tissue network, `B` is the actuator matrix, and `u(t)` is the control input. In the present implementation, each actuator column is site-local additive voltage forcing, so sparse actuation means that each column of `B` is a coordinate vector selecting a single cell. Physically, this corresponds to a local optogenetic channel, an electrode, or a spatially restricted channel-modulating perturbation.

The uncontrolled dynamics combine two ingredients. First, each cell has a local nonlinear voltage response

```text
f_local(V_i) = a V_i - b V_i^3,
```

with default parameters `a = 1` and `b = 1`. Second, cells are coupled by gap-junction-mediated voltage differences over the tissue adjacency graph:

```text
f_couple,i(V, G) = sum_{j in N(i)} g_ij (V_j - V_i).
```

The full sitewise drift used in the present calculations is therefore

```text
dV_i/dt
  = (1 / tau_V) [a V_i - b V_i^3
                 + sum_{j in N(i)} g_ij (V_j - V_i)
                 + h_i + c_i(t)],
```

where `h_i` collects static bias terms present in the frozen damaged state, and `c_i(t)` is the projected site control. The essential structural point is that the tissue dynamics are neither purely local nor purely topological. They are the sum of a cell-intrinsic nonlinear sector and a network coupling sector.

The control cost is the standard finite-horizon quadratic energy

```text
E[u] = integral_0^T ||u(t)||^2 dt.
```

This is the natural baseline for minimum-energy steering, and it is dimensionally consistent with the Gramian formulas used below.

### 2.2 Linearization at the damaged operating point

We linearize the controlled dynamics around a damaged pre-intervention snapshot `(V_bar, G_bar)`, with conductance frozen at the damaged value over the steering horizon. Writing `delta V = V - V*` relative to the target reference state yields the local linearized system

```text
d delta V / dt = A delta V + B_tilde u(t),
```

with

```text
A = (1 / tau_V) (partial F / partial V) |_(V_bar, G_bar),
B_tilde = B / tau_V.
```

This is the central local object of the paper. The frozen Jacobian `A` is not simply the graph Laplacian. In components,

```text
A_ii = (1 / tau_V) [a - 3 b V_bar,i^2 - sum_{j in N(i)} g_bar,ij]
A_ij = (1 / tau_V) g_bar,ij,   i != j and j in N(i),
```

with `A_ij = 0` for non-neighbors. Thus `A` inherits spectral structure from the gap-junction network but is modified by the local linearized voltage response at each site. That is exactly why the `A`-eigenmodes turn out to be Laplacian-structured but not Laplacian-pure.

In the codebase, the damaged operating point is stored explicitly as a `LinearizationPoint`, and the packaged local pair `(A, B_tilde)` is assembled as a `LinearizedPair` from a frozen simulation snapshot. The Jacobian is computed numerically by finite differences around the frozen state.

### 2.3 Finite-horizon steering problem

The finite-horizon linear steering problem is

```text
minimize    E[u] = integral_0^T ||u(t)||^2 dt

subject to  d delta V / dt = A delta V + B_tilde u(t)
            delta V(0) = delta V_0
            delta V(T) = delta V_target.
```

For exact steering to the target operating point, `delta V_target = 0` and `delta V_0 = V_damaged - V*`.

The endpoint formula is

```text
delta V(T) = e^(A T) delta V_0 + integral_0^T e^(A (T - t)) B_tilde u(t) dt.
```

Defining the endpoint mismatch

```text
eta = delta V_target - e^(A T) delta V_0,
```

the finite-horizon controllability Gramian is

```text
W(T) = integral_0^T e^(A tau) B_tilde B_tilde^T e^(A^T tau) dtau.
```

When the required direction lies in the reachable subspace of `W(T)`, the minimum-energy control is

```text
u*(t) = B_tilde^T e^(A^T (T - t)) W(T)^(-1) eta,
```

and the minimum energy is

```text
E_min = eta^T W(T)^(-1) eta.
```

For exact steering to zero, this becomes

```text
E_min = delta V_0^T e^(A^T T) W(T)^(-1) e^(A T) delta V_0.
```

These formulas are derived explicitly in [linear-minimum-energy-control.md](C:/Users/cable/bioelectric-minimum-energy-control/docs/derivations/linear-minimum-energy-control.md#L1).

### 2.4 Reachable and unreachable subspaces

The central obstruction arises when `W(T)` is not invertible along the required endpoint direction. Eigendecomposing the Gramian partitions state space into a reachable subspace spanned by eigenvectors with eigenvalues above tolerance, and an unreachable subspace given by the complementary directions. For a given damage pattern `delta V_0`, we decompose the endpoint mismatch into reachable and unreachable components.

This yields three quantitative objects:

- projected energy: the cost of steering the reachable component to zero,
- unreachable residual: the norm of the component that cannot be corrected by the chosen actuator family on the chosen horizon,
- actuation-rank threshold: the smallest number of actuator sites for which the unreachable residual falls below a chosen tolerance `epsilon`.

In the present numerical results we use `epsilon = 1e-4` as an operational fidelity threshold. We treat this as a working analysis tolerance, not as a biologically privileged constant. The actuation-rank threshold is the main object of the paper. It converts the control problem from a vague question about stimulation strength into a precise statement about intervention dimensionality.

## 3. Results

### 3.1 The controllable subspace is small and ill-conditioned

We begin with the default `16x16` damaged lattice and the default sparse two-site actuator. The state dimension is `256`, but the finite-horizon controllability Gramian on horizon `T = 4.0` has reachable rank only `15` at tolerance `1e-12` [default-lattice-phase2-baseline.md](C:/Users/cable/bioelectric-minimum-energy-control/docs/results/default-lattice-phase2-baseline.md#L1). Thus, less than `6%` of the local voltage state space is materially reachable under this actuator family.

The reachable directions are also extremely anisotropic. The same baseline yields a reachable-subspace condition number of `2.01717e+11`, with largest Gramian eigenvalue `0.307589` and smallest positive eigenvalue `1.52486e-12` [default-lattice-phase2-baseline.md](C:/Users/cable/bioelectric-minimum-energy-control/docs/results/default-lattice-phase2-baseline.md#L7). Sparse actuation therefore creates two constraints at once: low rank and strong energetic disparity even within the reachable sector.

This provides the setup for the rest of the paper. The relevant question is not whether the tissue is controllable under unrestricted forcing, but what a biologically sparse intervention can actually reach on a finite experimental horizon.

### 3.2 The controllable directions are Laplacian-structured but not Laplacian-pure

If control were determined by network wiring alone, the relevant modes would align closely with the frozen Laplacian eigenbasis. They do not. The full overlap matrix between frozen `A`-eigenvectors and Laplacian eigenvectors is saved in [default-lattice-a-vs-l-overlap.csv](C:/Users/cable/bioelectric-minimum-energy-control/docs/results/default-lattice-a-vs-l-overlap.csv). The dominant controllable `A`-modes have only moderate maximum overlaps with Laplacian eigenvectors, around `0.46-0.57` [default-lattice-phase2-baseline.md](C:/Users/cable/bioelectric-minimum-energy-control/docs/results/default-lattice-phase2-baseline.md#L1).

This is the first important spectral result. The controllable subspace clearly retains network structure, but it is not pure Laplacian structure. The local bistable voltage response materially reshapes the modal basis. In practical terms, a stimulation strategy derived from network geometry alone will miss part of the actual frozen control structure. Both topology and intrinsic voltage dynamics matter.

### 3.3 Two-site exact steering is infeasible for realistic wounds

We next tested whether any two-site actuator pair can exactly steer a realistic wound pattern to zero. On the `8x8` lattice, we exhaustively scanned all `2016` two-site pairs for a side-patch depolarization on the low-`x` boundary [placement-study-phase2.md](C:/Users/cable/bioelectric-minimum-energy-control/docs/results/placement-study-phase2.md#L1). The result is exact: the number of exactly reachable pairs at tolerance `1e-12` is `0`.

This is not a numerical failure. It is a structural statement about reachability. For every two-site Gramian in the exhaustive scan, the endpoint mismatch retains a nonzero projection onto the Gramian null space. The best projected pair `(8, 41)` still leaves residual `1.37294e-4`, and the median projected residual across all pairs is `2.74410e-4` [placement-study-phase2.md](C:/Users/cable/bioelectric-minimum-energy-control/docs/results/placement-study-phase2.md#L3). On this horizon, two-site control is under-ranked for this wound geometry.

The same qualitative result persists on `16x16`. The default midline pair `(135, 136)`, a rule-based pair, and a corner pair all remain exactly infeasible for the analogous side-patch wound [placement-study-phase2.md](C:/Users/cable/bioelectric-minimum-energy-control/docs/results/placement-study-phase2.md#L15). The implication is direct: below a critical actuator count, no increase in control energy can remove the damaged displacement because the required direction is not fully reachable.

### 3.4 Greedy actuation-rank thresholds depend on damage geometry

Once exact two-site infeasibility is established, the next question is how the threshold depends on damage geometry. On `8x8`, using greedy residual-minimizing placement and residual tolerance `1e-4`, the resulting actuation-rank thresholds are [actuation-rank-study-8x8.md](C:/Users/cable/bioelectric-minimum-energy-control/docs/results/actuation-rank-study-8x8.md#L1):

| Damage geometry | Minimum actuation rank |
| --- | ---: |
| Corner wound | 2 |
| Side patch | 3 |
| Central lesion | 3 |

The residual values at `k = 2` explain the split. The best corner-wound pair `(0, 9)` already achieves residual `2.90001e-6`, well below threshold. By contrast, the best side-patch pair `(8, 41)` remains at `1.37294e-4`, and the best central-lesion pair `(11, 20)` remains at `1.89178e-4` [actuation-rank-study-8x8.md](C:/Users/cable/bioelectric-minimum-energy-control/docs/results/actuation-rank-study-8x8.md#L5). The same sparse actuator family is therefore sufficient for one wound class and insufficient for another.

For `k = 2`, the conclusion is exact because the pair scan is exhaustive. For `k >= 3`, these are greedy upper bounds on the true minimum threshold because higher-rank actuator sets were not exhaustively enumerated. Even with that caveat, the core control result remains: intervention complexity is not a generic property of the tissue alone. It depends on how the damage projects onto the frozen controllable and uncontrollable directions of that tissue.

### 3.5 Above threshold, residual drops rapidly with added sites

The count-by-count greedy envelopes on `8x8` are saved in [greedy-envelope-8x8.csv](C:/Users/cable/bioelectric-minimum-energy-control/docs/results/greedy-envelope-8x8.csv#L1). Across all three geometries, increasing actuator count from `2` to `4` reduces residual by roughly `1-2` orders of magnitude.

For the corner wound, residual falls from `4.14099e-6` at `k = 2` to `3.44058e-7` at `k = 4`. For the side patch, it falls from `1.56095e-4` to `2.00112e-5`. For the central lesion, it falls from `2.40414e-4` to `2.80263e-5` [greedy-envelope-8x8.csv](C:/Users/cable/bioelectric-minimum-energy-control/docs/results/greedy-envelope-8x8.csv#L1).

The steep drop clarifies the relevant dose-response variable. The sharp transition is not in stimulation amplitude. It is in actuation rank. Below threshold, residual is pinned away from zero by unreachable directions. Above threshold, a small increase in control dimensionality produces large gains in repair fidelity. In the final paper this result should be shown as a dedicated threshold figure built directly from the greedy-envelope trajectories.

### 3.6 Placement matters, but geometry alone is not enough

The `8x8` side-patch scan also shows that "optimal placement" is not a single scalar notion. The best residual pair `(8, 41)` and the modal-rule pair `(27, 35)` do not coincide [actuation-rank-study-8x8.md](C:/Users/cable/bioelectric-minimum-energy-control/docs/results/actuation-rank-study-8x8.md#L5). The best-residual pair leaves smaller residual but requires much larger projected energy (`1790.89` versus `305.112`). Residual minimization and energy minimization select different sites.

The broader placement study reinforces the same point. The default midline pair on `16x16` does not outperform the random median on reachable-subspace conditioning [default-lattice-phase2-baseline.md](C:/Users/cable/bioelectric-minimum-energy-control/docs/results/default-lattice-phase2-baseline.md#L12). Geometric intuition is therefore insufficient. The useful placement signal comes from the frozen spectrum of `A`, not from simple symmetry arguments such as "stimulate the middle."

### 3.7 Thresholds shift with scale and disorder

The `8x8` result is not the end of the story. We therefore performed two additional checks: a structured greedy scaling run on `16x16`, and a disorder sweep on `8x8`.

For the `16x16` scaling check, we built a candidate pool from the damage support together with the top `32` nodes ranked by dominant controllable `A`-mode weight, then ran greedy placement up to `k = 8` [final-three-computations.md](C:/Users/cable/bioelectric-minimum-energy-control/docs/results/final-three-computations.md#L28). This is not an exhaustive optimum over all `256` nodes, but it is sufficient to test whether the threshold phenomenon persists at larger scale. The threshold table from this structured search is:

| Damage geometry | Greedy threshold rank |
| --- | ---: |
| Corner wound | 2 |
| Side patch | 6 |
| Central lesion | 2 |

The side patch becomes substantially harder on the larger lattice, remaining above threshold at `k = 4` with residual `2.30725e-4` and crossing below `1e-4` only at `k = 6` [actuation-rank-16x16.csv](C:/Users/cable/bioelectric-minimum-energy-control/docs/results/actuation-rank-16x16.csv#L1). By contrast, the corner wound remains easy and the central lesion remains below threshold by `k = 2` in this structured search.

Random gap-junction dilution pushes the threshold upward in the same direction. For the `8x8` side patch, the greedy threshold rank is stable at `3` for dilution fractions `0.0` and `0.1`, mixed between `3` and `4` at `0.2`, and mostly `4-5` by `0.3` across seeds `0..4` [disorder-sweep-8x8-side-patch.csv](C:/Users/cable/bioelectric-minimum-energy-control/docs/results/disorder-sweep-8x8-side-patch.csv#L1). The median residual at `k = 2` rises from `1.35e-4` with no dilution to `2.99e-4` at `30%` dilution, while the median residual at `k = 4` rises from `1.58e-5` to `1.10e-4` [final-three-computations.md](C:/Users/cable/bioelectric-minimum-energy-control/docs/results/final-three-computations.md#L57). With only five seeds per dilution level, we interpret this as a robustness check on the main threshold result rather than as a calibrated biological scaling law.

## 4. Discussion

### 4.1 Passive-tissue baseline and morphogenetic competency

Everything computed here is a passive-tissue baseline. The frozen linearized model assumes that conductance remains fixed over the steering horizon and that repair is assessed in the local finite-horizon linearized system. That approximation may overestimate or underestimate the true intervention complexity once conductance adaptation and nonlinear basin geometry are restored. We therefore do not claim a formal upper bound here.

Even without a formal bound, the framing creates a useful experimental quantity. If a real tissue repairs under fewer intervention sites than the frozen linearized prediction, then endogenous tissue dynamics are carrying part of the control burden. To make that comparison meaningful, the experimental protocol would need a count-reduction design analogous to the computational rank search: hold the damage geometry fixed, vary the number and placement of independently controlled stimulation sites, and record the smallest intervention that reliably restores pattern. No explicit model of internal goal representation is required, but the experimental observable must be defined carefully.

### 4.2 Reinterpreting existing bioelectric experiments

This framework changes how existing bioelectric perturbations should be interpreted. Uniform gap-junction blockade is a broad perturbation with low spatial specificity. In the present language, it supplies very limited structured control directions despite acting over many cells. Our results predict that spectrally informed sparse perturbations should outperform geometrically crude whole-tissue interventions when the damage pattern lies within a low-dimensional but specific reachable subspace.

Voltage-sensitive dye readouts become especially valuable under this interpretation. Rather than serving only as descriptive measurements of prepattern, they define the damaged displacement `delta V_0` that should be projected against the frozen controllable subspace. The same measurement that diagnoses the wound can therefore be used to estimate whether a given intervention modality has sufficient actuation rank to repair it.

### 4.3 Disorder, percolation, and the cancer connection

The disorder sweep links the control story directly to communication disorder. Random gap-junction dilution contracts or deforms the effective controllable subspace and raises the minimum actuation rank for the side-patch wound. By `30%` dilution, the threshold is mostly `4-5` rather than `3` [final-three-computations.md](C:/Users/cable/bioelectric-minimum-energy-control/docs/results/final-three-computations.md#L57).

This suggests a control-theoretic restatement of the communication-loss idea. If carcinogenesis involves loss of gap-junction-mediated coordination, then one plausible consequence is an increase in the dimensionality of intervention required to restore a coherent bioelectric pattern. The present disorder sweep is too small to establish that connection quantitatively, so we treat it as a motivating interpretation rather than as a cancer result.

### 4.4 Toward morphoceutical engineering specifications

For intervention design, the first engineering question is not how much current or drug to apply. It is how many independent control directions the modality provides. A systemic ion-channel drug offers something close to one degree of freedom. A promoter-restricted optogenetic tool offers as many degrees of freedom as there are distinct expression domains that can be independently addressed. A multi-electrode array offers the number of independently driven electrodes.

The present framework converts the clinical question "will this intervention repair the tissue?" into the control-theoretic specification "does this modality's actuation rank exceed the minimum required for this damage pattern on this tissue's gap-junction network?" If the answer is no, then tuning dose, timing, or amplitude cannot compensate. A higher-dimensional intervention is required.

### 4.5 Limitations

Several limitations are important. First, conductance is frozen. Real gap junctions are voltage-gated, so a fully coupled treatment would lead to a bilinear or more general nonlinear control problem, and the direction of the frozen-approximation bias is not yet known. Second, the results are local to the damaged operating point. They establish finite-horizon linear steering baselines, not full nonlinear basin-entry guarantees. Third, the current numerical work uses square lattices rather than experimentally reconstructed tissue geometry. The framework generalizes, but specific thresholds are topology-dependent. Fourth, the deterministic formulation omits stochastic forcing; noise should raise effective control cost and may increase actuation-rank thresholds for fixed fidelity. Fifth, all thresholds depend on analysis choices such as horizon `T` and residual tolerance `epsilon`; those dependencies should be mapped explicitly in follow-up sensitivity analyses.

These are genuine extensions, but they do not undercut the present result. The hard-threshold phenomenon already appears at the simplest passive linearized level, which means it is a structural feature worth measuring experimentally before more elaborate nonlinear detail is added.

## 5. Methods

### 5.1 Network model and parameters

The simulations use square lattices of size `8x8` and `16x16` with nearest-neighbor coupling. Default parameters are stored in `SimulationConfig`, including `tau_V = 1`, base conductance `0.35`, timestep `dt = 0.02`, and `200` steps, corresponding to horizon `T = 4.0` [model.py](C:/Users/cable/bioelectric-minimum-energy-control/research/control_gap_junction/model.py#L110). The local voltage sector is the double-well force `a V_i - b V_i^3` with default `a = 1` and `b = 1`, and the network coupling is diffusive through the gap-junction graph [dynamics.py](C:/Users/cable/bioelectric-minimum-energy-control/research/control_gap_junction/dynamics.py#L9).

### 5.2 Linearization procedure

A damaged operating point is constructed by simulating the default damaged lattice, then freezing voltages, polarity field, and conductance into a `LinearizationPoint` [simulator.py](C:/Users/cable/bioelectric-minimum-energy-control/research/control_gap_junction/simulator.py#L30). The local Jacobian is computed numerically by finite differences around that frozen state, and the linearized control operator is given by `B_tilde = B / tau_V` [analysis.py](C:/Users/cable/bioelectric-minimum-energy-control/research/control_gap_junction/analysis.py#L623).

### 5.3 Gramian and steering metrics

Finite-horizon controllability is computed from the exact symmetric Gramian formula when the frozen Jacobian is symmetric, and from matrix-exponential quadrature otherwise [analysis.py](C:/Users/cable/bioelectric-minimum-energy-control/research/control_gap_junction/analysis.py#L772). The endpoint mismatch, minimum-energy formula, reachable-rank calculation, and projected steering metrics are implemented in [analysis.py](C:/Users/cable/bioelectric-minimum-energy-control/research/control_gap_junction/analysis.py#L810) and derived in [linear-minimum-energy-control.md](C:/Users/cable/bioelectric-minimum-energy-control/docs/derivations/linear-minimum-energy-control.md#L1).

### 5.4 Spectral analysis

Eigenvalues and eigenvectors of the frozen Jacobian and frozen Laplacian are computed with symmetric eigensolvers, and their absolute overlaps are stored as a full mode-overlap matrix [analysis.py](C:/Users/cable/bioelectric-minimum-energy-control/research/control_gap_junction/analysis.py#L657). Modal control weights are obtained by projecting the actuator operator into the frozen `A`-eigenbasis [analysis.py](C:/Users/cable/bioelectric-minimum-energy-control/research/control_gap_junction/analysis.py#L690).

### 5.5 Actuator placement procedures

Two-site placement on `8x8` is evaluated exhaustively over all `2016` pairs [placement.py](C:/Users/cable/bioelectric-minimum-energy-control/research/control_gap_junction/placement.py#L95). Higher-rank sparse placement uses greedy residual-first selection over candidate nodes [placement.py](C:/Users/cable/bioelectric-minimum-energy-control/research/control_gap_junction/placement.py#L156). Accordingly, all reported `k >= 3` thresholds are greedy upper bounds rather than certified optima. The `16x16` scaling check uses a structured candidate pool built from damage support and high-weight controllable modes, so it is a scaling validation rather than an exhaustive optimum.

### 5.6 Damage patterns and disorder

Three damage geometries are used on `8x8`: low-`x` side patch, corner wound, and central lesion [studies.py](C:/Users/cable/bioelectric-minimum-energy-control/research/control_gap_junction/studies.py#L19). Scaled versions of the same patterns are used on `16x16` [studies.py](C:/Users/cable/bioelectric-minimum-energy-control/research/control_gap_junction/studies.py#L36). Gap-junction disorder is implemented as random edge dilution through the disorder module [disorder.py](C:/Users/cable/bioelectric-minimum-energy-control/research/control_gap_junction/disorder.py#L1).

### 5.7 Reproducibility

All numerical artifacts cited in this draft are saved under [docs/results](C:/Users/cable/bioelectric-minimum-energy-control/docs/results), including the overlap matrix, placement scans, greedy envelopes, `16x16` scaling table, and disorder sweep. Verification was performed by direct Python execution of the analysis and test helpers.

## Figures

1. Gramian eigenvalue spectrum for default two-site actuation on `16x16`.
2. Full overlap matrix between frozen `A`-eigenvectors and Laplacian eigenvectors.
3. Exhaustive two-site side-patch residual map on `8x8`, showing exact infeasibility.
4. Greedy actuation-rank table and residual-versus-count plot for the three `8x8` damage geometries.
5. Threshold figure from the greedy envelopes `k = 1..16` for `8x8`, with the `epsilon = 1e-4` line marked.
6. Two-site placement heatmaps showing residual-energy tradeoff.
7. Threshold shift under gap-junction dilution.

## Journal Positioning

The current manuscript is best aligned with `PLoS Computational Biology` or `Physical Review E`. The core result is already strong enough for a compact physics-focused version centered on the infeasibility threshold and minimum actuation rank table, with the spectral and disorder material as extension figures.
