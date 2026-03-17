# Controlled Dynamics Note

## Purpose

This note fixes the controlled state equation for the minimum-energy
bioelectric-control project. Its job is to make the control problem explicit
enough that the Phase 2 linear steering calculation is well-posed.

The key decision is the actuator operator `B`: what physical intervention is
allowed, where it acts, and how it enters the voltage dynamics.

## Base Uncontrolled System

The inherited adaptive gap-junction voltage dynamics have the form

```text
tau_V dV_i/dt = f_i(V, G),
```

where, in the current nearest-neighbor lattice model,

```text
f_i(V, G)
  = a V_i
  - b V_i^3
  + sum_{j ~ i} G_ij (V_j - V_i)
  + h_i
  + s_i
  - inj_i
  + xi_i(t).
```

For the control project, the control input must be separated from the inherited
drive and injury terms rather than mixed into them.

## Controlled State Equation

We define the controlled dynamics as

```text
tau_V dV/dt = F(V, G) + B u(t),
```

or equivalently

```text
dV/dt = (1 / tau_V) F(V, G) + (1 / tau_V) B u(t).
```

Here:

- `V in R^N` is the tissue voltage state
- `G` is the conductance state or frozen conductance field, depending on phase
- `u(t) in R^m` is the control input
- `B in R^(N x m)` maps control channels into tissue sites

The control cost is taken to be quadratic:

```text
E[u] = integral_0^T ||u(t)||^2 dt.
```

This is the default minimum-energy control objective carried into all later
phases.

## Physical Meaning Of The Control Input

Three actuator interpretations are possible in principle:

1. **Site-local injected current**
   `u_k(t)` is an externally delivered current at site or site-group `k`.
2. **Optogenetic depolarization / hyperpolarization**
   `u_k(t)` represents the effective current induced by light-gated ion
   channels on an addressed subset of cells.
3. **Gap-junction conductance modulation**
   `u_k(t)` acts on `G_ij` rather than directly on `V_i`.

Only the first two give a clean additive control term in the voltage equation.
The third is scientifically interesting but belongs to a later extension,
because it changes the control geometry from additive state forcing to bilinear
or parameter control.

## Default Phase 1 Choice

The default control model for Phases 1 and 2 is:

```text
site-local additive voltage forcing
```

with

```text
tau_V dV/dt = F(V, G) + B u(t)
```

and `B` selecting a physically addressable set of stimulated cells.

This is the best default because:

- it matches current injection and optogenetic stimulation at the model level
- it yields a standard linear control problem after linearization
- the controllability Gramian is well-defined and directly computable
- it avoids prematurely introducing bilinear control structure

## Admissible `B` Families

We will distinguish three useful actuator classes.

### 1. Full-Site Control

```text
B = I_N
```

Every cell is individually addressable.

Use:

- theoretical lower bound on minimum control energy
- control-spectrum benchmarking

Risk:

- biologically optimistic

### 2. Sparse Site Control

```text
B = [e_{i_1}, e_{i_2}, ..., e_{i_m}]
```

Only a selected subset of cells is stimulated.

Use:

- experimentally realistic actuation
- optimal targeting questions

This is likely the main actuator family for the paper.

### 3. Mode-Restricted Control

```text
B = [b_1, ..., b_m]
```

Columns of `B` are predefined spatial modes, for example low-order Laplacian
eigenvectors or region-level illumination masks.

Use:

- theoretical connection between controllability and network modes
- coarse optical actuation where individual cells are not independently targeted

## Linearization Contract

Once an operating point `(V_bar, G_bar)` is chosen, the voltage dynamics
linearize to

```text
delta V_dot = A delta V + B_tilde u(t),
```

with

```text
A = (1 / tau_V) dF/dV |_(V_bar, G_bar),
B_tilde = (1 / tau_V) B.
```

This is the system used for:

- controllability checks
- Gramian computation
- exact minimum-energy steering in the linear regime

The important implication is simple:

```text
if B is not fixed now, the Gramian problem is not fixed now.
```

So choosing `B` is not implementation detail. It is part of the mathematical
statement of the project.

## Phase 1 Decision

The project will proceed with this default:

```text
Control variable:
  u(t) = additive site-local stimulation current

Default actuator family:
  sparse site control

Reference actuator family:
  full-site control, used only as a lower-bound benchmark
```

This gives the cleanest path to the Phase 2 linearized minimum-energy control
derivation while preserving experimental relevance.

## Deferred Choices

These are intentionally not fixed in this note:

- basin-entry criterion
- steering horizon `T`
- damaged initial state family
- target-attractor neighborhood definition

Those depend on the local geometry around the chosen operating point and should
be fixed after the linearization is written down, not before.
