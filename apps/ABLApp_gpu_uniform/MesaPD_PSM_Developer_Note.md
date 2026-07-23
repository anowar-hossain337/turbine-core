<!-- Reason for edit start: document the full Phase 1 to Phase 4 MesaPD/PSM integration in one developer-facing note so the app structure, file edits, and runtime flow are easy to understand later. -->
# MesaPD + PSM Developer Note

## 1. Purpose

This note documents the MesaPD and PSM integration that was added to `ABLApp_gpu_uniform` from the initial scaffold up to the current Phase 4 moving-body setup.

The goal of the work was:

- keep the original ABL-generated numerics available from day one
- preserve the existing static urban obstacle path
- add one optional moving object using MesaPD + PSM
- make the moving-object feature easy to switch on and off from `input.prm`

## 2. Short Answer

`PSM` stands for `Partially Saturated Cells Method`.

## 3. High-Level Design

The implementation uses one executable with two runtime branches:

- the original ABL-only branch
- the optional MesaPD + PSM branch

The two branches are compiled together, and `input.prm` selects which branch runs.

This means:

- no separate application is required
- switching `PSM.enabled` on or off does not require a rebuild
- static urban obstacles and the moving body can coexist in the same run

## 4. Phase Summary

### Phase 1

Purpose:

- prepare the build system, code generation, and runtime config parsing
- keep real runtime coupling disabled

Outcome:

- the app could build both the plain ABL path and the future PSM path
- `input.prm` gained the new MesaPD / PSM / MovingBody control blocks
- enabling the feature at runtime still aborted intentionally

### Phase 2

Purpose:

- activate real MesaPD + PSM runtime coupling
- start with one stationary sphere only
- keep the existing static box obstacle active

Outcome:

- when all three feature toggles were enabled together, the app created one MesaPD sphere
- the sphere was mapped into the flow using PSM
- the app replaced the normal ABL stream-collide step with the PSM sweep sequence

### Phase 3

Purpose:

- extend the stationary sphere into a kinematic moving sphere
- drive the sphere from trajectory keyframes in `input.prm`

Outcome:

- the sphere state is updated every timestep
- the static urban box still comes from the normal boundary path
- the moving sphere comes from MesaPD + PSM

### Phase 4

Purpose:

- extend the kinematic moving-body path from a sphere-only case to a sphere-or-box case
- add prescribed rotation through a constant angular-velocity input
- keep the existing static urban box obstacle active

Outcome:

- the moving object can now be a sphere, box, or cube alias
- the body orientation is updated every timestep together with the translation state
- the original ABL branch and the Phase 3 sphere path remain available through prm settings

## 5. Files Edited And Why

### `ABLApp_gpu_uniform/CMakeLists.txt`

Purpose of edit:

- add a second Python code-generation target for the local PSM kernels
- link the executable against MesaPD and PSM dependencies

What changed:

- kept the original `LBMGeneration.py` target
- added `ABLApp_gpu_uniform_PSM_Generation` from `ABLPSMGeneration.py`
- linked:
  - `walberla::mesa_pd`
  - `walberla::lbm_mesapd_coupling`
  - `ABLApp_gpu_uniform_PSM_Generation`

Why it matters:

- this is what allows one executable to contain both generated code paths

### `ABLApp_gpu_uniform/ABLCodegenCommon.py`

Purpose of edit:

- centralize the shared ABL numerical setup

What changed:

- moved common stencil, layout, force model, SGS model, field typedefs, kernel-info metadata, and codegen settings into one helper module

Why it matters:

- both the plain ABL generator and the PSM generator now use the same numerical setup
- this prevents the PSM path from drifting away from the established ABL numerics

### `ABLApp_gpu_uniform/LBMGeneration.py`

Purpose of edit:

- refactor the original ABL generator to use `ABLCodegenCommon.py`

What changed:

- replaced duplicated local setup with `create_abl_common_setup(...)`
- reused common field typedef builders and kernel-info builders

Why it matters:

- the legacy ABL path stayed behavior-equivalent
- the generator became reusable beside the new PSM generator

### `ABLApp_gpu_uniform/ABLPSMGeneration.py`

Purpose of edit:

- create a local PSM generator that matches the ABL application numerics

What changed:

- generated:
  - `waLBerlaABLPSM_Sweep`
  - `waLBerlaABLPSM_InitializeDomainForPSM`
  - `waLBerlaABLPSM_KernelInfo`
- used the same stencil, layout, force model, and optimization metadata as the ABL path

Why it matters:

- this made PSM compatible with the current ABL setup instead of forcing a switch to a stock tutorial-style PSM solver

### `ABLApp_gpu_uniform/MovingBodyConfig.h`

Purpose of edit:

- add a small shared runtime parser for moving-body related prm blocks

What changed:

- added parsing for:
  - `MesaPD`
  - `PSM`
  - `MovingBody`
  - `Trajectory`
- in Phase 3, extended parsing to store sorted trajectory keyframes as real data
- in Phase 4, added `boxEdgeLength`, `initialRotation`, and `angularVelocity`

Why it matters:

- `main.cu` receives one structured runtime config instead of manually parsing many blocks in many places

### `ABLApp_gpu_uniform/main.cu`

Purpose of edit:

- integrate the optional MesaPD + PSM runtime path into the ABL application

What changed:

- added early validation of MesaPD / PSM / MovingBody config
- added optional MesaPD particle storage, shape storage, accessor, and PSM helper objects
- added one moving-body creation path with sphere and box support
- added PSM initialization sweeps before timestep 0
- added a timeloop branch:
  - legacy ABL stream-collide if PSM is off
  - PSM sweep sequence if PSM is on
- in Phase 3:
  - added moving-sphere state evaluation
  - added trajectory interpolation
  - added per-timestep sphere position and velocity update before particle mapping
- in Phase 4:
  - generalized the sphere state into a moving-body state
  - added box and cube-alias support through MesaPD `Box`
  - added per-timestep body rotation and angular-velocity updates before particle mapping

Why it matters:

- this file is the runtime dispatcher and the real integration point between the ABL solver, the static obstacle path, and the moving MesaPD body

### `ABLApp_gpu_uniform/input.prm`

Purpose of edit:

- expose the new feature to the user through prm settings

What changed:

- added:
  - `MesaPD { ... }`
  - `PSM { ... }`
  - `MovingBody { ... }`
  - `Trajectory { ... }`
- kept the existing static obstacle definition under `Boundaries -> Body`
- in Phase 3, added multiple trajectory keyframes and adjusted VTK output cadence for visible motion
- in Phase 4, switched the example to a kinematic rotating cube case

Why it matters:

- this file is now the user-facing control surface for enabling, disabling, and configuring the moving object

## 6. Generated Files Produced By The Integration

These are not hand-edited directly, but they are created from the Python generators and are part of the final application build:

- from `LBMGeneration.py`
  - `waLBerlaABL_*`
- from `ABLPSMGeneration.py`
  - `waLBerlaABLPSM_Sweep.*`
  - `waLBerlaABLPSM_InitializeDomainForPSM.*`
  - `waLBerlaABLPSM_KernelInfo.h`

These generated files live in the build tree and are compiled into the executable by CMake.

## 7. How The Files Are Connected

The connection chain is:

1. `ABLCodegenCommon.py`
   - defines the shared ABL numerical contract

2. `LBMGeneration.py`
   - generates the normal ABL kernels and boundary handling

3. `ABLPSMGeneration.py`
   - generates the PSM kernels using the same numerical contract

4. `CMakeLists.txt`
   - builds both generated code paths into one executable

5. `input.prm`
   - tells the executable whether the moving-body branch should stay off or be activated

6. `MovingBodyConfig.h`
   - parses the user prm blocks into a runtime config object

7. `main.cu`
   - decides which runtime branch to execute
   - allocates and updates the MesaPD moving body
   - runs the PSM sweep sequence when enabled

## 8. Runtime Application Flow

### Build Time

At build time:

- `LBMGeneration.py` generates the plain ABL kernels
- `ABLPSMGeneration.py` generates the PSM kernels
- CMake compiles both into one application

### Startup Time

At startup:

- `main.cu` reads `input.prm`
- `MovingBodyConfig.h` parses the optional MesaPD / PSM / MovingBody / Trajectory blocks
- `MovingBodyConfig.h` also reads the optional box-size and rotational controls for Phase 4
- the app checks whether the moving-body branch is:
  - fully disabled
  - partially enabled by mistake
  - fully enabled

### Legacy ABL Branch

If the moving-body feature is off:

- the app behaves like the normal ABL GPU uniform solver
- static obstacles still come from the generated boundary path

### PSM Branch

If the moving-body feature is on:

- `main.cu` creates one MesaPD moving body
- PSM helper fields are allocated
- the body is mapped into the fluid domain
- PDFs in intersecting cells are initialized consistently with the moving-body state

Then the timeloop uses:

1. particle mapping
2. particle velocity setup
3. PSM sweep
4. particle force reduction

instead of the single normal ABL stream-collide step

### Static And Moving Obstacles Together

The static urban object and the moving body are handled by different mechanisms:

- static urban box:
  - comes from `Boundaries -> Body`
  - handled through the standard generated boundary / flag-field path

- moving body:
  - comes from MesaPD particle storage
  - coupled through PSM

This separation is why they can coexist in the same simulation.

## 9. Phase 4 Motion Logic

In the current Phase 4 version:

- the body is kinematic
- `representation sphere`, `representation box`, and `representation cube` are supported
- the translation trajectory is read from `Trajectory -> point` blocks
- points are sorted by timestep
- the body center is interpolated linearly between keyframes
- the linear velocity is computed from the active segment
- the orientation starts from `MovingBody.initialRotation`
- the orientation is advanced from `MovingBody.angularVelocity`
- the body state is applied before the PSM particle-mapping sweep each timestep

## 10. Important User Controls In `input.prm`

### Feature toggles

- `MesaPD.enabled`
- `PSM.enabled`
- `PSM.phase1ScaffoldOnly`
- `MovingBody.enabled`

These must be consistent. The real coupling path expects all three main feature toggles to be enabled together.

### Moving body definition

- `MovingBody.mode`
- `MovingBody.representation`
- `MovingBody.radius`
- `MovingBody.boxEdgeLength`
- `MovingBody.initialPosition`
- `MovingBody.initialVelocity`
- `MovingBody.initialRotation`
- `MovingBody.angularVelocity`

### Trajectory definition

- `Trajectory.type`
- repeated `Trajectory.point` blocks containing:
  - `step`
  - `center`

### Static urban obstacle

- remains under:
  - `Boundaries -> Body`

## 11. Current Known Scope

The current implementation is intentionally narrow:

- one moving object
- sphere or box representation
- kinematic motion
- trajectory from prm keyframes
- prescribed constant angular velocity from prm

It is not yet a general rigid-body framework for arbitrary aircraft geometry.

## 12. Practical Interpretation

The final design should be understood like this:

- the original ABL application remains the base solver
- MesaPD + PSM was added as an optional overlay for one moving body
- the moving-body path was introduced gradually so the original numerics stayed available from the beginning
- the static urban obstacle path was preserved instead of being replaced

That was the main design intention throughout all four phases.
<!-- Reason for edit end: document the full Phase 1 to Phase 4 MesaPD/PSM integration in one developer-facing note so the app structure, file edits, and runtime flow are easy to understand later. -->
