# State

**Current Phase:** 01
**Current Phase Name:** Control-Ready Model Specification
**Status:** Phase 1 in progress
**Current Plan:** Phase 1 scaffold implemented; next work is to choose basin-entry and intervention-start conventions before Phase 2 solvers.
**Total Phases:** 6
**Total Plans in Phase:** —
**Last Activity:** 2026-03-17
**Last Activity Description:** Packaged the frozen operating point, local Jacobian, and rescaled actuator matrix into a Phase 1 linearized-pair helper and wrote the Phase 1 summary artifact.

## Open Questions

- What basin-entry criterion is strict enough to avoid counting transient near-hits as control success?
- Which intervention-start convention `t_0` gives the most experimentally meaningful damaged snapshot family?
- Which comparison linearizations beyond the default damaged snapshot are worth carrying into Phase 2: healthy-attractor or saddle-neighborhood?

**Convention Lock:**
- Metric signature: not-applicable (nonrelativistic tissue network)
- Natural units: lattice units
- Coordinate system: 2D square lattice with nearest-neighbor graph
