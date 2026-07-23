<!-- Reason for edit start: document the Phase 4a to Phase 4d moving-body work and the OpenMesh/CMake compatibility fix in one app-local summary so the file-level responsibilities and phase-by-phase capabilities are easy to revisit later. -->
# Phase 4a to 4d Modification Summary

## Overview

This note summarizes the moving-body work completed for `ABLApp_gpu_uniform` from Phase `4a` through `4d`, plus the separate OpenMesh/CMake compatibility fix that was needed to remove the convex-polyhedron build/runtime blocker.

The same `input.prm` file was reused as the rolling test case while the work advanced from box placement to rotating box, then to prototype convex polyhedron, and finally to STL-driven convex polyhedron. That means the current checked-in `input.prm` reflects the latest `4d` state, not the exact older `4a`, `4b`, or `4c` test values.

The common runtime chain across these phases is:

`input.prm` -> `MovingBodyConfig.h` -> `main.cu` -> `PSMSweepCollection.h` -> `ParticleAndVolumeFractionMappingSweepsGPU.h`

In simple words:

- `input.prm` describes what body should exist and how it should move.
- `MovingBodyConfig.h` parses those settings into a shared runtime config object.
- `main.cu` validates the config, creates the MesaPD particle and shape, and chooses which mapping path to use.
- `PSMSweepCollection.h` dispatches to the right mapping sweep while keeping the shared velocity and force-reduction sweeps unchanged.
- `ParticleAndVolumeFractionMappingSweepsGPU.h` contains the actual GPU overlap/mapping implementation for sphere, box, and convex-polyhedron occupancy.

## Phase 4a

### 1. Goal

Add a `box` moving-body representation without removing the existing `sphere` path, so the application can place and run a box obstacle through the same MesaPD + PSM runtime pipeline.

### 2. Which files were modified

- `apps/ABLApp_gpu_uniform/main.cu`
- `apps/ABLApp_gpu_uniform/MovingBodyConfig.h`
- `apps/ABLApp_gpu_uniform/input.prm`
- `src/lbm_mesapd_coupling/partially_saturated_cells_method/codegen/PSMSweepCollection.h`
- `src/lbm_mesapd_coupling/partially_saturated_cells_method/codegen/ParticleAndVolumeFractionMappingSweepsGPU.h`

### 3. How they are connected

- `MovingBodyConfig.h` parses the moving-body representation, box size, position, velocity, and trajectory data from `input.prm`.
- `main.cu` reads that config, creates the MesaPD box particle, computes a safe interaction radius from the box half-diagonal, and instantiates the shared PSM collection with the `box` representation token.
- `PSMSweepCollection.h` keeps one collection interface but selects either the sphere mapper or the box mapper at construction time.
- `ParticleAndVolumeFractionMappingSweepsGPU.h` implements the actual GPU-side `BoxFractionMappingSweep`, fills `mappingUIDs`, uploads particle positions, and drives the shared overlap fields used by the downstream PSM sweep.
- `input.prm` was used as the runtime test driver for the box case.

### 4. Why the modification was needed

The stock `PSMSweepCollection` path was sphere-oriented. To place a box in the application without duplicating the entire downstream PSM workflow, the mapper selection had to become representation-aware while still reusing the same particle-velocity and force-reduction sweeps.

### 5. What capabilities the application has after Phase 4a

- It can run either a `sphere` or a `box` through the same PSM moving-body pipeline.
- The box is mapped into the `B`, `Bs`, and overlap/index fields needed by the generated PSM fluid update.
- The original sphere path remains available instead of being replaced.

## Phase 4b

### 1. Goal

Add proper box rotation so the box is no longer treated as axis-aligned in world coordinates.

### 2. Which files were modified

- `apps/ABLApp_gpu_uniform/main.cu`
- `apps/ABLApp_gpu_uniform/input.prm`
- `src/lbm_mesapd_coupling/partially_saturated_cells_method/codegen/ParticleAndVolumeFractionMappingSweepsGPU.h`

### 3. How they are connected

- `input.prm` provides `initialRotation` and `angularVelocity`.
- `main.cu` evaluates the moving-body state each timestep, including accumulated angular rotation, and writes the updated MesaPD particle orientation before mapping.
- `ParticleAndVolumeFractionMappingSweepsGPU.h` sends the box center, half-edge lengths, and inverse rotation matrix to the GPU box kernel so occupancy is tested in box body coordinates instead of with a world-axis shortcut.

### 4. Why the modification was needed

Phase `4a` was enough for box placement, but not for a rotating object. Without forwarding orientation to the GPU mapping kernel, a rotating box would still behave like an axis-aligned box in the overlap test, which would be physically and geometrically wrong.

### 5. What capabilities the application has after Phase 4b

- The box can stay fixed and rotate.
- The box can translate and rotate through the same kinematic state evaluator.
- The mapped obstacle shape now follows the actual MesaPD orientation instead of an axis-aligned approximation.

## Phase 4c

### 1. Goal

Add a first `convex_polyhedron` moving-body path using one hard-coded prototype mesh, while keeping the sphere and box paths intact.

### 2. Which files were modified

- `apps/ABLApp_gpu_uniform/main.cu`
- `apps/ABLApp_gpu_uniform/input.prm`
- `src/lbm_mesapd_coupling/partially_saturated_cells_method/codegen/PSMSweepCollection.h`
- `src/lbm_mesapd_coupling/partially_saturated_cells_method/codegen/ParticleAndVolumeFractionMappingSweepsGPU.h`

### 3. How they are connected

- `main.cu` accepts `convex_polyhedron` as a valid representation token, checks that ConvexPolyhedron support exists, and builds one prototype mesh from the current `boxEdgeLength`.
- `main.cu` creates a MesaPD `ConvexPolyhedron` shape and passes the new representation token into the shared PSM collection.
- `PSMSweepCollection.h` adds a third mapping branch, `convexPolyhedron`, beside `sphere` and `box`.
- `ParticleAndVolumeFractionMappingSweepsGPU.h` adds `ConvexPolyhedronFractionMappingSweep`, which extracts face normals and face points from the OpenMesh-backed convex mesh, uploads them to the GPU, and runs the convex-polyhedron occupancy kernel.
- `input.prm` was used to switch the test case from box to prototype convex polyhedron.

### 4. Why the modification was needed

This was the bridge from simple analytic shapes to mesh-based shapes. Before loading a real turbine mesh, the application first needed to prove that one convex-polyhedron body could pass through the same MesaPD + PSM runtime path.

### 5. What capabilities the application has after Phase 4c

- It can run `sphere`, `box`, and `convex_polyhedron` from one shared PSM collection.
- The convex-polyhedron path can map one prototype OpenMesh-backed convex body.
- The app fails early if convex-polyhedron is requested without OpenMesh-backed ConvexPolyhedron support.

## Phase 4d

### 1. Goal

Load a user-provided mesh file, convert it into one convex hull, and run that hull through the existing single-convex-polyhedron path.

### 2. Which files were modified

- `apps/ABLApp_gpu_uniform/MovingBodyConfig.h`
- `apps/ABLApp_gpu_uniform/main.cu`
- `apps/ABLApp_gpu_uniform/input.prm`

### 3. How they are connected

- `input.prm` now provides `meshFile`, `representation convex_polyhedron`, and the target `boxEdgeLength` fit envelope.
- `MovingBodyConfig.h` adds the new `meshFile` setting to the shared runtime config and normalizes surrounding quotes so standard prm syntax like `"windTurbineRFP.stl"` works correctly.
- `main.cu` resolves the mesh path relative to the run directory or executable directory, validates that the file exists, reads it with OpenMesh, builds one convex hull with `QHull`, scales that hull to fit inside `boxEdgeLength`, recenters it, and creates the MesaPD `ConvexPolyhedron` from the result.
- `main.cu` also logs whether the convex-polyhedron body came from the old prototype path or from the new mesh-file path.

### 4. Why the modification was needed

Phase `4c` proved the convex-polyhedron runtime path, but it still used only a hard-coded prototype shape. Phase `4d` was needed to connect a real external geometry file to that path without redesigning the downstream mapper or fluid update.

### 5. What capabilities the application has after Phase 4d

- The application can load an external STL/OBJ/OFF file through `meshFile`.
- The loaded mesh is converted into one single convex hull and scaled to fit the configured `boxEdgeLength`.
- The resulting body can use the same prescribed translation and rotation controls as the earlier sphere and box cases.
- The `meshFile` parser now tolerates the normal quoted prm form.

### Important current limitation of Phase 4d

The loaded turbine mesh is not yet used as its exact non-convex surface. The current implementation converts the whole file into one single convex hull. That is why this phase is a bridge toward more complex future geometry work, not the final multi-piece turbine representation.

## OpenMesh / CMake compatibility fix

### 1. Goal

Remove the build/runtime blocker that prevented the convex-polyhedron path from being recognized even when OpenMesh was actually available.

### 2. Which file was modified

- `walberla-my-walberla-work/CMakeLists.txt`

### 3. How the modification is connected to Phases 4c and 4d

- The convex-polyhedron runtime path in `main.cu` and `PSMSweepCollection.h` depends on `WALBERLA_MESAPD_CONVEX_POLYHEDRON_AVAILABLE`.
- That feature flag depends on OpenMesh being detected correctly during the waLBerla build.
- The root `CMakeLists.txt` OpenMesh section was updated so both the system-package path and the `FetchContent` path export the legacy `OPENMESH_*` variables expected by older downstream checks.

### 4. Why the modification was needed

Modern OpenMesh discovery often exposes imported CMake targets like `OpenMeshCore` and `OpenMeshTools`, but older waLBerla logic still looked for legacy variables such as `OPENMESH_FOUND`, `OPENMESH_LIBRARIES`, and `OPENMESH_INCLUDE_DIRS`. Without the compatibility bridge, convex-polyhedron support stayed disabled even though OpenMesh itself was present.

### 5. What capabilities the build now has

- waLBerla can recognize OpenMesh in both the system-package path and the fetched-source path.
- The MesaPD ConvexPolyhedron code path can be enabled reliably.
- Phases `4c` and `4d` can build and run without the earlier OpenMesh recognition mismatch.

## Current capability summary after Phase 4d

The application now has these moving-body capabilities inside the one-way-coupled MesaPD + PSM runtime path:

- `sphere` representation
- `box` representation
- rotating `box` representation
- prototype `convex_polyhedron` representation
- file-driven single-convex-hull `convex_polyhedron` representation from STL/OBJ/OFF

The main missing step for a more realistic turbine body is no longer “can the app read a mesh file?” The remaining gap is “can the app use multiple convex pieces or another richer mesh representation instead of collapsing the entire turbine into one convex hull?”
<!-- Reason for edit end: document the Phase 4a to Phase 4d moving-body work and the OpenMesh/CMake compatibility fix in one app-local summary so the file-level responsibilities and phase-by-phase capabilities are easy to revisit later. -->
