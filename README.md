# waLBerla-wind

**waLBerla-wind** is a [waLBerla](https://www.walberla.net/) extension for
wind flow simulation. It provides a C++ module set that models atmospheric
boundary-layer wind fields using the Lattice Boltzmann Method (LBM) and is
designed to couple with turbine-performance models.

## Features

- **Logarithmic wind profile** – initialises the domain with a physically
  motivated log-law velocity profile.
- **LBM coupling** – `WindLBMModel` initialises PDF fields from the wind
  profile and applies a body-force term during the collision step.
- **Inflow boundary condition** – `WindInflowBoundary` prescribes equilibrium
  populations on inlet faces.
- **Block-structured parallelism** – built on waLBerla's block-forest, so it
  runs on shared-memory, distributed-memory (MPI), and hybrid architectures.

## Repository Layout

```
.
├── CMakeLists.txt                  Top-level CMake (requires waLBerla)
├── src/
│   └── wind/
│       ├── WindParameters.h/.cpp   Log-law profile & simulation parameters
│       ├── WindField.h/.cpp        Block-based 3-D wind velocity field
│       ├── lbm/
│       │   └── WindLBMModel.h      LBM ↔ wind-field coupling (header-only)
│       └── boundary/
│           └── WindBoundary.h      Inflow boundary condition (header-only)
├── apps/
│   └── WindSimulation/
│       └── WindSimulation.cpp      Example D3Q19 SRT wind simulation
└── tests/
    ├── WindParametersTest.cpp      Unit tests for WindParameters
    └── WindFieldTest.cpp           Unit tests for WindField
```

## Requirements

| Dependency | Version |
|------------|---------|
| waLBerla   | ≥ 6.0   |
| CMake      | ≥ 3.15  |
| C++        | ≥ 17    |
| MPI        | optional (inherited from waLBerla build) |

## Build

```bash
# Clone this repository next to your waLBerla source tree, then:
mkdir build && cd build
cmake .. -DWALBERLA_DIR=/path/to/walberla
make -j$(nproc)

# Run tests
ctest --output-on-failure
```

## Quick Start

```bash
# After building, run the bundled example:
./apps/WindSimulation/WindSimulation
```

The example simulates incompressible wind flow on a 64 × 32 × 32 lattice
using a D3Q19 SRT model.  Inflow velocity is set by the log-law profile with
a reference speed of 0.05 (lattice units) at height 10.

## Wind Profile

The logarithmic wind profile is:

```
u(z) = u_ref * ln(z / z0) / ln(z_ref / z0)
```

| Symbol   | Parameter            | Default |
|----------|----------------------|---------|
| `u_ref`  | `inflowSpeed`        | 0.05    |
| `z_ref`  | `referenceHeight`    | 10.0    |
| `z0`     | `roughnessLength`    | 0.1     |

## License

waLBerla-wind is distributed under the **GNU General Public License v3**.
See the individual source file headers for the full license text.