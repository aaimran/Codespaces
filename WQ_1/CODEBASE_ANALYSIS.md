# WaveQLab3D Codebase Analysis Report

**Generated**: 2026-06-10  
**Source**: `/workspaces/Codespaces/WQ/src`

---

## 1. Project Overview

**WaveQLab3D** is a Fortran-based computational code for **3D seismic wave propagation and earthquake rupture dynamics**. It solves the elastic wave equation in curvilinear coordinates (supporting complex geometries) with support for:

- Nonplanar frictional fault interfaces
- Off-fault viscoplasticity (Drucker-Prager plasticity)
- Spatially variable elastic properties
- Multiple friction laws (rate-and-state, slip-weakening, linear friction)
- Advanced attenuation modeling (anelastic-Q, frequency-dependent Q, constant-Q)
- Perfectly Matched Layers (PML) for absorbing boundaries
- Moment tensor point sources

**Authors**: Kenneth Duru, Sam Bydlon, Eric Dunham, Kyle Withers (parallelization by Hari Radhakrishnan)  
**License**: MIT  
**Validation suite**: SCEC Spontaneous Rupture Code Validation (SRCV) benchmarks (TPV26-34 series)

---

## 2. Tech Stack

| Component | Details |
|-----------|---------|
| **Language** | Fortran 2008+ (free-form) |
| **Build System** | CMake 3.5+ |
| **Compilers** | gfortran, ifort, ftn (Cray) |
| **Parallelization** | MPI (OpenMPI / MPICH / Intel MPI / Cray MPI) |
| **Precision** | Double (64-bit) working precision; single-precision output option |
| **Optional I/O** | HDF5 (currently disabled in build) |
| **Docs** | Doxygen |

---

## 3. Directory Structure

```
WQ/
├── src/                          # Main source (~45,866 lines, 40+ modules)
│   ├── main.f90                 # Entry point
│   ├── domain.f90               # Domain coordinator
│   ├── block.f90                # Computational block operations
│   ├── datatypes.f90            # Central type definitions
│   ├── common.f90               # Global constants & precision
│   ├── elastic.f90              # Wave equation discretization
│   ├── plastic.f90              # Drucker-Prager plasticity
│   ├── RHS_Interior.f90         # Spatial derivative RHS (3478 lines)
│   ├── JU_xJU_yJU_z6.f90        # High-order FD stencils (23,392 lines)
│   ├── material.f90             # Material properties & Q models
│   ├── grid.f90                 # Grid generation (2740 lines)
│   ├── metrics.f90              # Jacobian/metric derivatives
│   ├── Interface_Condition.f90  # Fault slip laws (2889 lines)
│   ├── iface.f90                # Interface data structures
│   ├── BoundaryConditions*.f90  # Boundary conditions (3 variants)
│   ├── time_step.f90            # Runge-Kutta time integration
│   ├── pml.f90                  # Perfectly Matched Layer
│   ├── mpi3dbasic.f90           # MPI initialization
│   ├── mpi3dcomm.f90            # MPI communication patterns
│   ├── mpi3dio.f90              # Parallel I/O
│   ├── mpi3d_interface.f90      # Interface-specific MPI
│   ├── seismogram.f90           # Station output
│   ├── fault_output.f90         # Fault surface output
│   ├── slice_output.f90         # 2D slice output
│   ├── plane_output.f90         # Plane-based output
│   ├── mms.f90                  # Method of Manufactured Solutions
│   ├── moment_tensor.f90        # Point moment tensor sources
│   └── ...
├── build/                        # CMake build artifacts
├── bin/                          # Executables (waveqlab3d, pre_wql3d)
├── cmake/                        # CMake scripts
├── conf/                         # Doxygen config
├── inputfile/                    # 50+ test problem input files (TPV series)
├── simulation/                   # SLURM batch scripts & test runners
├── python/                       # Analysis & visualization (Python)
├── auxilary/                     # Post-processing tools
├── test_problems/                # Reference truth solutions
├── build.sh                      # Build automation
└── check_requirements.sh         # Dependency verification
```

---

## 4. Key Files

| File | Purpose | Lines |
|------|---------|-------|
| `main.f90` | Entry point — init domain, run time steps | 137 |
| `domain.f90` | Domain coordinator (blocks + interfaces) | 808 |
| `datatypes.f90` | All type definitions | 345 |
| `block.f90` | Block init, time step, MMS evaluation | 526 |
| `time_step.f90` | RK integration (1st–4th order) | 200+ |
| `RHS_Interior.f90` | Wave equation RHS computation | 3,478 |
| `JU_xJU_yJU_z6.f90` | High-order FD operators (pre-generated) | 23,392 |
| `material.f90` | Material & attenuation initialization | 1,820 |
| `grid.f90` | Grid generation (Cartesian, curvilinear, file) | 2,740 |
| `Interface_Condition.f90` | Fault friction laws (slip-weakening, RSF) | 2,889 |
| `metrics.f90` | Metric tensor & Jacobian derivatives | 1,284 |
| `mpi3dcomm.f90` | MPI neighbor exchange | 781 |
| `seismogram.f90` | Seismometer station I/O | 634 |

---

## 5. Architecture

### Domain Decomposition Pattern

```
domain_type
  ├── block_type[1..N]  (each with grid, fields, material, boundaries, PML)
  └── iface_type[1..M]  (each with slip S, velocity V, traction T, state W)
```

Each block can be decomposed across MPI ranks. Interfaces (faults) couple neighbouring blocks.

### Module Dependency Hierarchy

```
common.f90
    └── datatypes.f90
            ├── grid.f90
            ├── material.f90
            ├── fields.f90
            ├── elastic.f90 / plastic.f90
            ├── iface.f90
            ├── boundary.f90
            └── block.f90
                    └── domain.f90
                            └── time_step.f90
                                    └── main.f90

MPI layer (orthogonal):
    mpi3dbasic.f90 → mpi3dcomm.f90 → mpi3dio.f90 → mpi3d_interface.f90
```

### Pluggable Physics (via `response=` config key)

| Physics | Modules |
|---------|---------|
| `elastic` | `elastic.f90`, `RHS_Interior.f90`, `JU_xJU_yJU_z6.f90` |
| `plastic` | + `plastic.f90`, `plastic_material.f90` |
| `anelastic`, `anelastic-Q`, `anelastic-Q8`, etc. | + `material.f90` Q models + memory variables |

### Fault Friction Laws (`Interface_Condition.f90`)

- `locked` — no slip
- `linear` — Coulomb friction
- `slip-weakening_friction` — linear-weakening
- `rate-and-state_friction` — aging law with state variable θ

---

## 6. Computation Flow

```
Initialize Domain
  ├─ Read namelists (&problem_list, &block_list, &output_list)
  ├─ Create blocks (grids, materials, fields, MPI decomposition)
  ├─ Create interfaces (fault parameters, friction law)
  └─ Set up output structures (seismograms, fault, plane slices)

Time Loop (n = 1 → nt):
  └─ RK stage loop (typically 5 stages, 4th-order):
      ├─ MPI ghost exchange (6 neighbors per block)
      ├─ Boundary conditions (free surface, characteristic)
      ├─ Interior RHS:
      │    ├─ Spatial FD derivatives (JU_xJU_yJU_z6)
      │    ├─ Elastic stress-strain
      │    └─ PML damping
      ├─ Interface coupling:
      │    ├─ Transform to fault-aligned coordinates
      │    ├─ Compute hat variables
      │    └─ Enforce friction law → update slip rates
      ├─ Plasticity update (Drucker-Prager, if enabled)
      ├─ Attenuation update (memory variables, if enabled)
      └─ RK field update: F += dt * B_rk * DF

  Output (stride-based):
      ├─ Seismograms (binary)
      ├─ Fault surface (slip, velocity, traction, state)
      └─ Plane/slice data
```

---

## 7. Physics & Algorithms

### Elastic Wave Equation (first-order system)
```
∂v/∂t  = ∇·σ/ρ + f/ρ
∂σ/∂t  = λ(∇·v)I + 2μ ε̇
```

### Attenuation (GSLS)
- Multiple relaxation mechanisms (4 or 8) with memory variables η_k
- `dη_k/dt = ω_k(σ − σ_R)` — frequency-dependent effective moduli

### Slip-Weakening Friction
- `f(d) = f_s + (f_0 − f_s)·exp(−d/L_w)` until slip d ≥ L_w, then f = f_s

### Rate-and-State Friction
- `f = f_0 + a·ln(V/V_0) + b·ln(θV_0/L)`
- State evolution: `dθ/dt = 1 − θ·V/L` (aging law)

### Drucker-Prager Plasticity
- Yield: `F = √J₂ + μ·p − σ_y`
- Viscous regularization with viscosity η

---

## 8. Spatial Discretization

- Finite Difference on structured curvilinear grids
- Order configurable (2–6), via `order=` config key
- FD type: `traditional` (centered) or `upwind` or `DRP` (Dispersion Relation Preserving)
- Jacobian-scaled operators handle grid curvature (metrics.f90)
- Summation-by-Parts (SBP) operators for stability at boundaries

---

## 9. Input/Output

### Input Format (Fortran namelists)
```fortran
&problem_list
  name='data/output_prefix', problem='TPV26',
  response='elastic', nblocks=2, nt=10000, CFL=0.5
/
&block_list
  btp(1)%nqrs = 201, 201, 401
  btp(1)%rho_s_p = 2.67, 3.464, 6.0
/
```

### Output Formats
- Binary Fortran unformatted — primary large data (fault surfaces, planes)
- Parallel MPI-IO — distributed writes across ranks
- ASCII text — diagnostics, station seismograms

### Optional HDF5 (disabled by default)
- Module `hdf5_output.f90` present but not linked in CMake

---

## 10. Build & Run

```bash
# Build
./build.sh
# → bin/waveqlab3d   (main solver)
# → bin/pre_wql3d   (preprocessor)

# Run
mpirun -n 4 ./bin/waveqlab3d inputfile/tpv26.in
```

**Note**: MPI rank count must be even (minimum 2); serial mode uses 1.

---

## 11. Potential Issues & Concerns

| Issue | Details |
|-------|---------|
| **Code duplication** | Three `BoundaryConditions*.f90` variants; backup `original_RHS_interior.f90` still present |
| **Module name conflicts** | Multiple BC modules compiled — potential symbol clash risk |
| **Memory management** | Manual allocation/deallocation; no cleanup in error paths |
| **Large generated file** | `JU_xJU_yJU_z6.f90` at 23K lines — hard to verify or maintain by hand |
| **Documentation gaps** | Limited inline docs; Doxygen template may be stale |
| **No CI** | Tests defined in CMake but no automation pipeline visible |
| **Static load balancing** | No runtime MPI load redistribution between blocks |
| **MPI rank constraint** | Even-number constraint is strict and not documented prominently |

---

## 12. Test Coverage

- **50+ pre-configured test problems** in `inputfile/` (TPV26–34, elasticPML, constant-Q, rate-and-state, etc.)
- **Method of Manufactured Solutions** (MMS) module for numerical verification
- **Reference truth solutions** in `test_problems/truth/`
- **SCEC SRCV benchmarks** as primary validation standard
- CMake test targets: `mpi_test1`, `mpi_test2`, `premesh_test1`, etc.

---

*Report generated by Claude Code analysis of `/workspaces/Codespaces/WQ/src`*
