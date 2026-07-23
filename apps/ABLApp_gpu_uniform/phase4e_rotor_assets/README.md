<!-- Reason for edit start: document the expanded Phase 4e rotor-only STL assets so both the original one-blade prototype and the newer three-blade rotor set can be reused without re-deriving the orientation. -->
# Phase 4e Rotor Assets

This folder contains application-ready STL copies derived from:

- `blade_hub.stl`
- `matt-bladerev5.stl`

## Files

- `blade_hub_phase4e_rotor_aligned.stl`
  - Hub recentered at the origin and rotated so the rotor axis is along `+X`.
- `matt_bladerev5_phase4e_rotor_aligned.stl`
  - One blade shifted into the same local rotor frame and attached to the hub region.
- `hub_plus_1_blade_phase4e_preview.stl`
  - Combined preview mesh for quick inspection outside the solver.
- `matt_bladerev5_phase4e_rotor3_blade_a.stl`
  - Blade A for the rotor-only `hub + 3 blades` test.
- `matt_bladerev5_phase4e_rotor3_blade_b.stl`
  - Blade B pre-rotated around `+X` from the same shared local rotor frame.
- `matt_bladerev5_phase4e_rotor3_blade_c.stl`
  - Blade C pre-rotated around `+X` from the same shared local rotor frame.
- `hub_plus_3_blades_phase4e_preview.stl`
  - Combined preview mesh for quick inspection of the `hub + 3 blades` assembly.

## Intended Phase 4e Use

Use the files as one of these two convex-polyhedron assemblies under the same moving-body runtime:

- `hub + 1 blade`
- `hub + 3 blades`

The important assumption is:

- both files already share one local rotor coordinate frame,
- the rotor spin axis is `X`,
- no extra per-piece rotation is expected inside the current application path.

## Practical Meaning

These assets now support two practical prototypes:

- `1 hub`
- `1 blade`
- `1 hub`
- `3 blades`

The current `3-blade` set uses a small common phase shift so the assembly is more balanced for the existing Phase 4e centering logic, even though the later rotation-center code fix is still the proper long-term solution.
<!-- Reason for edit end: document the expanded Phase 4e rotor-only STL assets so both the original one-blade prototype and the newer three-blade rotor set can be reused without re-deriving the orientation. -->
