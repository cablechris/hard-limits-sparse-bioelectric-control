# Linear Minimum-Energy Control Note

## Purpose

This note states and solves the finite-horizon linear steering problem that
follows from the Phase 1 controlled-dynamics setup.

It is the first exact Phase 2 result: once the frozen local pair

```text
delta V_dot = A delta V + B_tilde u(t)
```

is fixed, the minimum-energy input and the controllability Gramian follow from
standard linear optimal control.

## Problem Statement

Let

```text
delta V(t) in R^N
u(t) in R^m
A in R^(N x N)
B_tilde in R^(N x m).
```

The finite-horizon linear steering problem is:

```text
minimize    E[u] = integral_0^T ||u(t)||^2 dt

subject to  delta V_dot = A delta V + B_tilde u(t)
            delta V(0) = delta V_0
            delta V(T) = delta V_target
```

For the control project, the Phase 2 local baseline is

```text
delta V_0 = V_damaged - V*
```

with `delta V_target = 0` for exact steering to the target operating point.

The nonlinear basin-entry condition

```text
delta V(T) in basin(V*)
```

is not solved exactly by the linear problem. In Phase 2, the linear result is
used in two ways:

- exact steering to a chosen linear target `delta V_target`
- local lower-bound or candidate-input construction for later nonlinear
  basin-entry tests

## Endpoint Formula

The state at time `T` is

```text
delta V(T)
  = e^(A T) delta V_0
  + integral_0^T e^(A (T - t)) B_tilde u(t) dt.
```

Define the endpoint mismatch vector

```text
eta = delta V_target - e^(A T) delta V_0.
```

Then admissible controls must satisfy

```text
integral_0^T e^(A (T - t)) B_tilde u(t) dt = eta.
```

## Controllability Gramian

The finite-horizon controllability Gramian is

```text
W(T) = integral_0^T e^(A tau) B_tilde B_tilde^T e^(A^T tau) dtau.
```

Equivalently, with the change of variables `tau = T - t`,

```text
W(T) = integral_0^T e^(A (T - t)) B_tilde B_tilde^T e^(A^T (T - t)) dt.
```

If `W(T)` is invertible, the pair `(A, B_tilde)` is controllable on the horizon
`[0, T]` and the minimum-energy problem has a unique exact solution for every
endpoint.

## Minimum-Energy Solution

Using a Lagrange multiplier `lambda in R^N` for the endpoint constraint gives

```text
L[u, lambda]
  = integral_0^T u(t)^T u(t) dt
  + 2 lambda^T (eta - integral_0^T e^(A (T - t)) B_tilde u(t) dt).
```

Stationarity with respect to `u` yields

```text
u*(t) = B_tilde^T e^(A^T (T - t)) lambda.
```

Substituting back into the endpoint constraint gives

```text
W(T) lambda = eta,
lambda = W(T)^(-1) eta.
```

Therefore the minimum-energy control is

```text
u*(t)
  = B_tilde^T e^(A^T (T - t)) W(T)^(-1)
    [delta V_target - e^(A T) delta V_0].
```

This is the standard exact finite-horizon result.

## Minimum Energy

Substituting `u*` into the cost gives

```text
E_min = eta^T W(T)^(-1) eta.
```

For exact steering to zero, `delta V_target = 0`, so

```text
eta = -e^(A T) delta V_0,
```

and the minimum energy is

```text
E_min
  = delta V_0^T e^(A^T T) W(T)^(-1) e^(A T) delta V_0.
```

Important sign and normalization note:

- this is the correct exact formula for the present state convention
- the shorter expression `delta V_0^T W(T)^(-1) delta V_0` is only a special
  case under additional conventions, for example if the dynamics are written in
  transformed coordinates where the homogeneous propagation is already removed

## Dimensional Check

Using the project conventions:

- `delta V` has voltage units
- `u` has actuator-amplitude units
- `A` has units `1 / time`
- `B_tilde` has units `voltage / (time * control)`

Then:

```text
e^(A t)                  dimensionless
B_tilde B_tilde^T        voltage^2 / (time^2 * control^2)
W(T)                     voltage^2 / (time * control^2)
W(T)^(-1)                time * control^2 / voltage^2
eta^T W(T)^(-1) eta      time * control^2
```

which matches

```text
integral_0^T ||u(t)||^2 dt.
```

So the Gramian formula is dimensionally consistent.

## Spectral Interpretation

The Gramian defines which state-space directions are cheap or expensive to
control over a finite horizon.

If

```text
W(T) q_k = w_k q_k,
```

then:

- large `w_k` means the direction `q_k` is easy to excite and cheap to steer
- small `w_k` means the direction `q_k` is weakly actuated over `[0, T]` and
  expensive to steer

For the present tissue network, `A` inherits its structure from the local
Laplacian-like coupling plus on-site terms. Therefore the expensive directions
are expected to align with the poorly actuated Laplacian modes of the frozen
network.

This statement is exact when the frozen Jacobian can be written as

```text
A = alpha I - G_bar L
```

or more generally when the on-site contribution commutes with the frozen
Laplacian term. In that case `A`, `e^(A t)`, and `W(T)` share the same modal
basis up to the actuator weighting. If the frozen on-site term is spatially
heterogeneous, the Laplacian-eigenvector interpretation remains approximate but
is still the right first comparison basis.

That is the direct route from the exact control formula to the scaling-law and
optimal-mode questions.

## Scope Limit

This note solves the exact linear endpoint problem. It does not yet prove
nonlinear basin entry. That nonlinear step remains a later phase.
