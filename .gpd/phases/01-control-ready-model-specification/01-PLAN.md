---
phase: 01-control-ready-model-specification
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - docs/derivations/controlled-dynamics.md
  - .gpd/PROJECT.md
  - .gpd/ROADMAP.md
  - research/control_gap_junction/model.py
  - research/control_gap_junction/analysis.py
interactive: false

conventions:
  units: "lattice units"
  metric: "not-applicable (nonrelativistic tissue network)"
  coordinates: "2D square lattice with nearest-neighbor graph"

approximations:
  - name: "additive control assumption"
    parameter: "control enters voltage sector as B u(t)"
    validity: "site-local current injection or optogenetic forcing modeled as additive voltage drive"
    breaks_when: "control acts primarily by changing conductances or other hidden variables"
    check: "verify that the control derivation and code keep B as an additive actuator operator"

contract:
  scope:
    question: "What controlled state equation, actuator operator, and cost functional define the minimum-energy steering problem?"
  claims:
    - id: "claim-controlled-model"
      statement: "The project fixes a control-ready voltage equation dV/dt = F(V,G) + B u(t) with an explicit actuator family and quadratic energy cost."
      deliverables: ["deliv-control-note", "deliv-control-ready-code"]
      acceptance_tests: ["test-controlled-equation"]
      references: ["ref-parent-lyapunov", "ref-parent-discussion"]
    - id: "claim-phase2-readiness"
      statement: "The chosen actuator operator B is specific enough that the Phase 2 linearized Gramian problem is mathematically well-posed."
      deliverables: ["deliv-control-note"]
      acceptance_tests: ["test-gramian-readiness"]
      references: ["ref-parent-lyapunov"]
  deliverables:
    - id: "deliv-control-note"
      kind: "derivation"
      path: "docs/derivations/controlled-dynamics.md"
      description: "Control-model note fixing F(V,G), B, and the quadratic intervention cost"
    - id: "deliv-control-ready-code"
      kind: "code"
      path: "research/control_gap_junction/"
      description: "Control-ready research package baseline aligned with the chosen additive actuator model"
  references:
    - id: "ref-parent-lyapunov"
      kind: "prior_artifact"
      locator: "C:/Users/cable/get-physics-done/docs/derivations/lyapunov-obstruction.md"
      role: "foundation"
      why_it_matters: "Provides the inherited voltage-sector structure and the drive-versus-smoothing baseline this control project generalizes"
      applies_to: ["claim-controlled-model", "claim-phase2-readiness"]
      must_surface: true
      required_actions: ["read", "use"]
    - id: "ref-parent-discussion"
      kind: "prior_artifact"
      locator: "C:/Users/cable/get-physics-done/research/adaptive_gap_junction/DISCUSSION.md"
      role: "scope"
      why_it_matters: "Defines the prior negative results and interpretation limits that this control project must build on without repeating"
      applies_to: ["claim-controlled-model"]
      must_surface: true
      required_actions: ["read", "use"]
  acceptance_tests:
    - id: "test-controlled-equation"
      subject: "claim-controlled-model"
      kind: "consistency"
      procedure: "Check that the derivation note, project docs, and code all use the same additive controlled voltage equation and quadratic cost definition"
      pass_condition: "The controlled equation, actuator family, and cost are explicit and identical across the note and code-facing project artifacts"
      evidence_required: ["deliv-control-note", "deliv-control-ready-code"]
    - id: "test-gramian-readiness"
      subject: "claim-phase2-readiness"
      kind: "consistency"
      procedure: "Verify that B is fixed as a site-local additive actuator family so the linearized pair (A,B) is well-defined for the next phase"
      pass_condition: "The project does not leave B ambiguous between additive voltage forcing and conductance modulation"
      evidence_required: ["deliv-control-note"]
  forbidden_proxies:
    - id: "fp-vague-control"
      subject: "claim-controlled-model"
      proxy: "Talking about intervention or stimulation without defining the actual actuator operator B"
      reason: "That would leave the minimum-energy control problem ill-posed"
    - id: "fp-bilinear-drift"
      subject: "claim-phase2-readiness"
      proxy: "Quietly switching between additive voltage forcing and conductance control depending on convenience"
      reason: "The Gramian and controllability analysis depend directly on which control geometry is chosen"
  links:
    - id: "link-controlled-model-note"
      source: "claim-controlled-model"
      target: "deliv-control-note"
      relation: "supports"
      verified_by: ["test-controlled-equation"]
    - id: "link-phase2-readiness-note"
      source: "claim-phase2-readiness"
      target: "deliv-control-note"
      relation: "supports"
      verified_by: ["test-gramian-readiness"]
  uncertainty_markers:
    weakest_anchors: ["The eventual experimental realism of the sparse site-control actuator family"]
    disconfirming_observations: ["If the relevant intervention mechanism acts mainly through conductance modulation, the additive B u(t) baseline will not be the right control geometry"]
---

<objective>
Fix the control-ready mathematical statement of the project so the next phase can compute linearized controllability and minimum-energy steering without ambiguity.

Purpose: turn the inherited adaptive gap-junction dynamics into an explicit optimal-control problem.
Output: a pinned controlled state equation, actuator operator family, and cost functional that make the Phase 2 Gramian problem well-posed.
</objective>

<tasks>

<task type="auto">
  <name>Task 1: Lock the controlled state equation and actuator family</name>
  <files>docs/derivations/controlled-dynamics.md, .gpd/PROJECT.md</files>
  <action>State the controlled voltage equation, identify the physical meaning of u(t), and fix the default actuator family as sparse site-local additive control with full-site control retained only as a benchmark lower bound.</action>
  <verify>Check that the note and project description use the same control geometry and do not mix additive voltage forcing with conductance modulation.</verify>
  <done>The repo contains one unambiguous control equation and one default B family for downstream derivations.</done>
</task>

<task type="auto">
  <name>Task 2: Make the Phase 2 linearization target explicit</name>
  <files>docs/derivations/controlled-dynamics.md, .gpd/ROADMAP.md</files>
  <action>Record the linearized pair (A, B_tilde) that will define the next phase and ensure the roadmap names the Gramian-based minimum-energy control problem in those terms.</action>
  <verify>Check that A and B_tilde are defined from the same controlled equation and that the roadmap references linearized minimum-energy control rather than vague steering language.</verify>
  <done>The next phase can start directly from a well-posed linearized control problem.</done>
</task>

<task type="auto">
  <name>Task 3: Keep the code scaffold aligned with the control formulation</name>
  <files>research/control_gap_junction/model.py, research/control_gap_junction/analysis.py</files>
  <action>Verify that the copied research package remains compatible with an additive-control extension and note any code surfaces that will need B-aware helpers in the next phase.</action>
  <verify>Check that no copied module hard-codes a conflicting control interpretation.</verify>
  <done>The research package is ready for linear-control utilities without a model rewrite.</done>
</task>

</tasks>

<verification>
The control equation, actuator operator, and cost functional must be explicit, consistent across docs, and sufficient to define the linearized controllability pair (A, B).
</verification>

<success_criteria>
The project has one clear control-ready model statement, one default actuator family, and one well-posed linearized target for Phase 2. No ambiguity remains about whether control acts additively on voltages or indirectly through conductances.
</success_criteria>

<output>
After completion, create `.gpd/phases/01-control-ready-model-specification/01-01-SUMMARY.md`.
</output>
