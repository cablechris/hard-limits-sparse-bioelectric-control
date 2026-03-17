# Phase 01-01 Summary

## Scope

This phase fixed the control-ready mathematical statement of the project so the
next phase can start from an explicit local linearized pair instead of vague
steering language.

## Completed

- Defined the controlled voltage equation as `dV/dt = (1 / tau_V) F(V, G) + (1 / tau_V) B u(t)`.
- Fixed the default actuator geometry as sparse site-local additive voltage forcing, with full-site control retained as a benchmark lower bound.
- Fixed the control cost as quadratic energy `E[u] = integral ||u(t)||^2 dt`.
- Fixed the default linearization target as a damaged pre-intervention snapshot with frozen conductance field.
- Added explicit code support for open-loop control traces projected through the actuator matrix.
- Added explicit code artifacts for:
  - frozen linearization points
  - local voltage-sector Jacobian `A`
  - rescaled actuator matrix `B_tilde = B / tau_V`
  - packaged local linearized pair `(A, B_tilde)`

## Main Artifacts

- Derivation note: `docs/derivations/controlled-dynamics.md`
- Control-ready model types: `research/control_gap_junction/model.py`
- Control and linearization helpers: `research/control_gap_junction/analysis.py`
- Simulation snapshot and pair packaging: `research/control_gap_junction/simulator.py`
- Verification checks: `tests/test_control_ready.py`

## Verification Status

The repository now has a consistent Phase 1 control-ready scaffold:

- the docs and code use the same additive control equation
- the actuator family is explicit
- the operating point convention is explicit
- the local linearized pair is explicit

Verification performed in this environment:

- direct Python invocation of the control-ready checks in `tests/test_control_ready.py`
- direct Python invocation of the smoke check in `tests/test_smoke.py`

`pytest` was not installed in the environment, so the test suite was not run
through `pytest`.

## Remaining Open Questions

- What basin-entry criterion is strict enough for nonlinear steering claims?
- What intervention-start convention `t_0` defines the damaged snapshot family most cleanly?
- Which comparison linearizations beyond the default damaged snapshot are worth carrying into Phase 2?

## Exit Condition

Phase 1 is complete enough for Phase 2 to begin local linear controllability and
minimum-energy steering work without ambiguity about the state equation,
actuator operator, or default operating point.
