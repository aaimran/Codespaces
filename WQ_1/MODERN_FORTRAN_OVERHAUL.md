# Modernizing WaveQLab3D with the Latest Fortran Ecosystem — Analysis & Migration Plan

*Date: 2026-06-11. Companion to `IMPROVEMENT_ROADMAP.md` (which identified GPU execution and engineering health as the top parity gaps) and `CODEBASE_ANALYSIS.md` (current architecture).*

**Scope of the proposed overhaul:** coarrays / PGAS "shared global arrays", native CPU/GPU parallelization without MPI (`do concurrent`), Fortran Package Manager (fpm), GASNet-EX, Caffeine, unit testing, CI, and Julienne.

---

## 0. Executive summary — the one fact that decides the architecture

The proposed stack is coherent — it is essentially the Berkeley Lab / LLVM-Flang modern-Fortran ecosystem (Caffeine, PRIF, Julienne, Assert from Rouson's group; GASNet-EX underneath; fpm and `do concurrent` from the broader community). But **its two pillars are mature on *disjoint* compilers today**:

| Capability | Mature on | Not available / immature on |
|---|---|---|
| `do concurrent` **GPU offload** | **nvfortran** (`-stdpar=gpu`, unified memory — most mature); HPE Cray `ftn`; ifx (Intel GPUs only, manual data movement, reported slow) | gfortran (CPU only), LLVM flang (under development) |
| **Coarrays** (native parallel Fortran) | gfortran+OpenCoarrays, ifx, HPE Cray `ftn`; **LLVM flang via PRIF/Caffeine (first demonstrated SC'25 workshop — very new)** | **nvfortran (no coarray support at all)** |

So "coarrays for communication + `do concurrent` for GPU" — the obvious reading of the proposal — **cannot be compiled by any mainstream compiler today** except HPE Cray's `ftn` (vendor-locked to Cray systems). The flang+PRIF+Caffeine convergence point will eventually fix this, but its multi-image lowering was first demonstrated at an SC'25 workshop — it is research-grade, not production-grade, in mid-2026.

**Therefore the recommendation is tiered:**

- **Adopt now (low risk, high value):** fpm, Julienne + Assert, CI, `do concurrent` for all interior kernels (CPU parallel everywhere, GPU offload via nvfortran). This directly closes Gap 1 (GPU) and Gap 5 (engineering health) of the roadmap.
- **Adopt behind an abstraction (hedge):** introduce a thin communication-interface module so that `mpi3dcomm`'s halo exchange can be swapped between MPI (production today) and coarrays/Caffeine (experimental branch) without touching physics code.
- **Do not do yet:** delete MPI wholesale and bet production on coarrays+Caffeine+GASNet. Revisit when flang ships both PRIF-based multi-image support *and* `do concurrent` offload in a release-quality compiler (plausibly 1–2 years out).

**Bottom line (see §6 for the full case):** the recommended stack is **MPI (kept) + ISO `do concurrent` kernels + a thin layer of OpenMP `target data` directives for GPU residency** — same modernization, different transport. Coarrays/Caffeine graduate from the experimental lane only when proven.

---

## 1. Baseline — what the overhaul is replacing

Current stack (from `CODEBASE_ANALYSIS.md`):

- Fortran 2008+, CMake, **MPI required** (block-structured domain decomposition; halo exchange in `mpi3dcomm`; even-rank constraint; static load balancing).
- No threading, no GPU path. Serial triple loops in `RHS_Interior.f90` with **per-point physics dispatch** (`if (M%anelastic_*)` chains inside the loops) and the 23K-line generated stencil file `JU_xJU_yJU_z6.f90`.
- Tests exist (MMS, 50+ TPV inputs, truth solutions) but assertions are ad-hoc and **no CI runs them**.
- HDF5 output optional/disabled; binary output via MPI-IO patterns.

---

## 2. Technology-by-technology analysis

### 2.1 Coarrays / PGAS — the "shared global array" model

**What it is.** Fortran 2008/2018/2023 native parallelism: the program runs as *images* (SPMD, like MPI ranks); arrays declared with codimensions (`real :: u(nx,ny,nz)[*]`) are directly addressable across images (`u(:, :, 1)[neighbor]`). Fortran 2018 adds teams, events, and collectives (`co_sum`, `co_max`, `co_broadcast`).

**Important conceptual correction:** coarrays are **PGAS (partitioned global address space), not shared memory**. There is no single "shared global array" — each image owns its partition, and remote access is communication (one-sided put/get). Treating a coarray as a free shared array (fine-grained remote indexing inside hot loops) produces catastrophic performance. The correct usage for WaveQLab3D is exactly its current pattern: bulk halo puts/gets at step boundaries — i.e., coarrays change the *syntax and runtime* of halo exchange, not its structure.

**What it buys WQ:**
- `mpi3dcomm`'s send/recv + derived-datatype machinery collapses to array-syntax one-sided puts: `u_halo(:,:,:)[east] = u(nx-nh+1:nx,:,:)` plus `sync images` — dramatically less code, compiler-checked.
- Reductions (`co_sum` for time-step minima, seismogram gathers) become one-liners.
- Removes the undocumented even-rank constraint only if the decomposition logic is also reworked (the constraint is algorithmic, not MPI's fault).

**How (if/when adopted):**
1. First isolate all communication behind an interface module (`comm_halo_exchange`, `comm_reduce_min`, `comm_gather_fault`) — pure refactor, MPI implementation unchanged.
2. Add a coarray implementation of that interface as a build option.
3. Validate with gfortran+OpenCoarrays and ifx on CPU clusters; benchmark halo-exchange latency vs. MPI at production scale (this is where PGAS implementations historically disappoint on some networks).

**Pros**
- Native language feature — no external API; compiler checks remote accesses; far less boilerplate than MPI derived datatypes.
- One-sided semantics map naturally to halo exchange and to modern RDMA networks (especially via GASNet-EX/Cray).
- Future-proof: this is where the Fortran standard and LLVM investment are going (PRIF).

**Cons**
- **No nvfortran support — directly conflicts with the only mature `do concurrent` GPU path.** This is the blocking con.
- Performance varies wildly by runtime (OpenCoarrays-over-MPI can be *slower* than direct MPI; Cray native is excellent; Caffeine/GASNet promising but young).
- WQ's halo buffers live in derived types with allocatable components — coarrays of such types are exactly the corner where compiler bugs concentrate.
- Parallel HDF5 and most I/O middleware are MPI-based; pure-coarray builds lose parallel I/O unless GASNet's MPI-interoperability mode is used (it can coexist, but that re-introduces MPI as a dependency — undermining "without MPI").
- Ecosystem debugging/profiling tooling (Score-P, Vampir, Nsight) understands MPI far better than coarrays.

### 2.2 Caffeine + GASNet-EX + PRIF

**What it is.** [Caffeine](https://github.com/BerkeleyLab/caffeine) (Berkeley Lab) is a parallel runtime implementing Fortran's multi-image features atop **GASNet-EX**, the exascale one-sided communication middleware (supports Slingshot, InfiniBand, EFA, Ethernet, shared memory). It is the first implementation of **PRIF** (Parallel Runtime Interface for Fortran), a compiler-agnostic spec letting any compiler target any runtime. Currently supports **LLVM flang ≥ 19 and gfortran ≥ 13**; flang lowering of multi-image features through PRIF/Caffeine was demonstrated in an SC'25 workshop paper.

**What it buys WQ:** an MPI-free communication substrate for the coarray version of the comm layer, with potentially better one-sided performance than OpenCoarrays-over-MPI, and alignment with where LLVM Fortran is headed. GASNet-EX is proven exascale technology (it underlies UPC++, Legion, Chapel).

**How:** only as the runtime behind the experimental coarray build (§2.1 step 3): flang + Caffeine + GASNet-EX on a test cluster; compare halo-exchange microbenchmarks and a TPV regression against the MPI build.

**Pros**
- GASNet-EX is mature, exascale-proven middleware with broad network support.
- PRIF decouples compiler from runtime — strategic insurance against vendor lock-in.
- Active, funded development (Berkeley Lab / CASS); responsive maintainers.

**Cons**
- The *Fortran-facing* layer is very new: flang multi-image lowering is a 2025 workshop demonstration, not a hardened production path.
- Adds two build dependencies (Caffeine, GASNet) with their own configure complexity — friction against the fpm simplicity goal.
- No GPU-aware communication story yet comparable to CUDA-aware MPI (relevant once fields live in GPU memory).
- Team/event coverage and performance at thousands of images are still being established.

### 2.3 `do concurrent` — native CPU/GPU parallelization without MPI directives

**What it is.** Standard Fortran loop-parallelism (`do concurrent` with Fortran 2023 `reduce`/locality specifiers). Compilers map it to multicore (all majors), or GPUs: **nvfortran `-stdpar=gpu`** (most mature — implicit unified-memory data movement, tunable with `-gpu=nomanaged`), **HPE Cray ftn**, **ifx** (Intel GPUs; currently requires manual OpenMP `target data` regions to avoid per-loop transfers and has known performance/correctness rough edges), **flang** (offload under development). The arXiv portability study ("Portability of Fortran's do concurrent on GPUs") confirms nvfortran as the reference implementation.

**What it buys WQ:** this is **the implementation vehicle for Roadmap Gap 1 (GPU)** — single-source Fortran, no CUDA/OpenACC/OpenMP dialect, runs parallel on CPU with every compiler and on NVIDIA GPUs today.

**How (concrete, in order):**
1. **Hoist the physics dispatch** out of the spatial loops in `RHS_Interior.f90` (the `if (M%anelastic_*)` chains around lines 255–262 and the per-point PML routing). Specialize loop nests per response type at block level. Required for *any* parallelization; benefits CPU immediately.
2. **Regenerate `JU_xJU_yJU_z6.f90`** with `do concurrent (k=..., j=..., i=...)` tight nests, no per-point calls into assumed-shape helpers; pure functions or inlined expressions only (offload compilers reject or deoptimize opaque calls).
3. Convert RK stage updates, memory-variable (eta4–eta9) updates, and PML auxiliary updates — all embarrassingly point-local.
4. Keep fault-interface and boundary SAT code on CPU initially (2D, cheap, branchy); move later if profiling justifies.
5. Manage residency: with nvfortran start with unified memory (`-stdpar=gpu`), then pin the nine field arrays + memory variables device-resident across RK stages and exchange only halos.
6. CI builds: gfortran (CPU, correctness) + nvfortran (GPU, performance) on every commit.

**Pros**
- Pure ISO Fortran — no directive dialects, no second source; CPU-parallel for free on all compilers.
- nvfortran path is genuinely production-grade (NVIDIA's flagship stdpar showcase) and unified memory makes incremental porting safe.
- Composes with MPI (current) *and* with future coarrays — orthogonal to the communication decision.

**Cons**
- Performance portability is not automatic: ifx needs manual data-region management; AMD GPUs need Cray `ftn` or future flang — today's practical GPU target is NVIDIA.
- Less control than CUDA/OpenACC (shared-memory staging, kernel fusion) — typically reaches a good fraction of hand-tuned performance, not 100%; stencil codes generally fare well.
- Requires the dispatch-hoisting + stencil-regeneration refactor first (significant but valuable work regardless).
- Hides data movement — careless adoption silently thrashes PCIe; requires profiling discipline.

### 2.4 Fortran Package Manager (fpm)

**What it is.** The community-standard build tool and package manager (`fpm build`, `fpm test`, `fpm.toml` manifest, automatic module-dependency resolution, dev-dependency support).

**What it buys WQ:** trivially correct module-order builds (currently hand-maintained in `CMakeLists.txt`), one-command test running (`fpm test` — integrates Julienne directly), frictionless dependency on Julienne/Assert (declared as dev-dependencies, fetched automatically), and drastically lower barrier for new contributors.

**How:**
1. Add `fpm.toml`; restructure to fpm conventions (`src/`, `app/` for the two executables `waveqlab3d` and `pre_wql3d`, `test/`).
2. Express HDF5/MPI via fpm metapackages (`dependencies: mpi`) and link flags.
3. **Keep CMake in parallel** for HPC-site installs (CMake is what supercomputing centers script against); fpm becomes the developer/CI workflow. Dual-build is low-cost since fpm needs almost no configuration.

**Pros:** near-zero config; standard in the modern Fortran ecosystem; makes Julienne adoption one manifest line; great CI ergonomics.
**Cons:** less flexible than CMake for exotic site toolchains/cross-compilation; preprocessor/generated-file workflows (the stencil generator) need custom build steps that fpm handles less gracefully; two build systems must be kept in sync (mitigate: CI builds both).

### 2.5 Julienne + Assert — unit testing and correctness checking

**What it is.** [Julienne](https://github.com/BerkeleyLab/julienne) (Berkeley Lab; Rouson, Bonachea, Rasmussen) is a Fortran-2023 correctness-checking and unit-testing framework with natural-language-style test idioms, **zero dependencies**, deliberate compiler portability, and — uniquely among Fortran test frameworks — support for testing *parallel* (multi-image) code. Its companion **Assert** library provides runtime-checkable design-by-contract assertions that compile away in production builds.

**What it buys WQ:** WQ has excellent *integration* assets (MMS, TPV truth solutions) but **no unit tests** — the Q-dispatch bugs fixed on 2026-06-11 (flag aliasing, tau denominator, weight lookup) were exactly the kind of module-level logic errors unit tests catch and integration seismogram comparisons smear out.

**How (highest-value targets first):**
1. **Attenuation coefficients:** for each `response=` variant, unit-test relaxation weights, tau spacing, and unrelaxed-modulus correction against analytically known values (would have caught the tau/16 bug and the 8M/4M flag aliasing directly).
2. **Friction laws:** slip-weakening and rate-and-state point solvers against closed-form/reference solutions.
3. **SBP operators:** accuracy order and the SBP property (∑ wᵢ uᵢ Dvᵢ + ∑ wᵢ vᵢ Duᵢ = boundary terms) on random fields — catches stencil-regeneration regressions (critical for the §2.3 rewrite).
4. **PML damping profiles** and `pml_damping_and_q` routing.
5. Sprinkle Assert preconditions on the namelist-input paths (positive moduli, Q bounds, npml vs. block size — turning silent misconfigurations into errors).
6. Wire into `fpm test` and CI; keep MMS/TPV as the integration tier above the unit tier.

**Pros:** zero dependencies; portable across compilers (important given the multi-compiler strategy above); parallel-aware testing matches the coarray experiments; assertions document physics invariants.
**Cons:** young v1.x project — API may move; small community vs. pFUnit (the alternative — pFUnit is more featureful but heavier, NASA-maintained, preprocessor-based); Fortran-2023 idioms require recent compilers (fine for CI, may constrain old site compilers).

### 2.6 CI (continuous integration)

**What/why:** already Roadmap Phase-0. All of the above only pays off if enforced automatically.

**How:** GitHub Actions matrix — {gfortran-13/-14, flang, nvfortran (self-hosted or container runner for GPU)} × {fpm, CMake}; stages: build → Julienne unit tests → MMS convergence check → 2–3 TPV regressions vs. truth solutions (tolerance-based seismogram comparison). Cache GASNet/Caffeine builds for the experimental lane. **Precondition: fix the WQ/ gitlink** — CI is meaningless while git tracks nothing inside `WQ/`.

---

## 3. The decision matrix

| Component | Verdict | Why |
|---|---|---|
| `do concurrent` (CPU + NVIDIA GPU) | **Adopt now** | Closes the #1 parity gap; pure ISO Fortran; mature on nvfortran; CPU-parallel everywhere |
| fpm (dual with CMake) | **Adopt now** | Near-zero cost; unlocks Julienne; developer/CI ergonomics |
| Julienne + Assert | **Adopt now** | Directly targets the demonstrated bug class; zero deps; parallel-aware |
| CI | **Adopt now** (after git fix) | Phase-0 roadmap item; everything else depends on it |
| Communication abstraction layer | **Adopt now** | Cheap refactor that makes the MPI-vs-coarray question reversible |
| Coarrays (replacing `mpi3dcomm`) | **Experimental branch only** | Blocked by nvfortran's lack of coarrays; runtime performance unproven for this code |
| Caffeine + GASNet-EX | **Experimental branch only** | flang/PRIF lowering is a 2025 workshop result; not production-hardened |
| Dropping MPI entirely | **Not yet** | Loses parallel HDF5, tooling, and the only mature GPU compiler; revisit when flang ships offload + PRIF together |

## 4. Overall pros and cons of the overhaul

**Pros**
- Single-source ISO-standard code: no OpenMP/OpenACC/CUDA dialects, no MPI boilerplate in physics code — substantially smaller, more readable, more teachable codebase (fits WQ's academic-niche strategy).
- Rides the funded LLVM/Berkeley trajectory (flang, PRIF, Caffeine, Julienne are all actively invested in) instead of WQ-specific infrastructure.
- The testing/CI/fpm layer directly remedies the engineering-credibility gap and the exact bug class found in the Q dispatch.
- GPU capability lands without forking the code per vendor.

**Cons / risks**
- **The compiler-matrix hole**: coarrays and GPU offload don't meet in one production compiler (except Cray). Mitigated by the tiered plan + comm abstraction.
- Bleeding-edge dependencies (Caffeine, flang offload, Julienne v1) can stall the project on upstream bugs — keep MPI/CMake escape hatches until two independent compilers pass the full suite per feature.
- nvfortran-centric GPU path is NVIDIA-locked near-term (AMD/Intel need Cray ftn or future flang/ifx maturity).
- Refactor risk: stencil regeneration + dispatch hoisting touch the numerical core — must land *after* the SBP-property and attenuation unit tests exist (sequence the safety net first).
- Two build systems and a comm abstraction add maintenance surface; justified only because they keep strategic options open.

## 5. Sequenced plan — overview (detailed phase-by-phase plan in §7)

1. **Week 0:** fix WQ/ git tracking. *(Blocks everything.)*
2. **Weeks 1–3 (Phase 0):** fpm.toml dual build; Julienne+Assert with attenuation/friction/SBP unit tests; GitHub Actions CI (gfortran+flang, fpm+CMake, unit+MMS+TPV tiers); delete dead code.
3. **Weeks 3–6:** communication abstraction module around `mpi3dcomm` (pure refactor, CI green throughout).
4. **Months 2–5 (Phase 1):** dispatch hoisting → stencil regeneration with `do concurrent` (validated by the new SBP unit tests) → nvfortran `-stdpar=gpu` bring-up with unified memory → device-resident fields with halo-only transfers. Deliverable: single-node GPU speedup on a TPV benchmark.
5. **Parallel experimental lane (time-boxed):** coarray implementation of the comm interface; build with gfortran/OpenCoarrays and flang/Caffeine/GASNet-EX; microbenchmark halo exchange vs. MPI; report. **Go/no-go criterion:** adopt for production only when (a) one compiler supports both PRIF-based images and `do concurrent` offload at release quality, and (b) halo-exchange performance ≥ MPI at target scale.
6. **Reassess at each flang release** — flang is the designated convergence point for the whole stack.

---

## 6. Recommended alternative: MPI + `do concurrent` + OpenMP `target data` residency

Having weighed the full proposed stack (§0–§5), the stack I recommend is deliberately more boring than the Berkeley-ecosystem endpoint:

> **Keep MPI as the communication transport. Write all kernels in ISO `do concurrent`. Manage GPU memory residency with a thin layer of OpenMP `target enter/exit data` directives. Skip coarrays/Caffeine in production entirely for now.**

### Why this beats the proposed stack

**1. Vendor-portable today, not in two years.** The same `do concurrent` kernels offload via nvfortran (NVIDIA), ifx (Intel GPUs), HPE Cray `ftn` (AMD GPUs on Frontier-class machines), and AMD's next-generation flang — *provided* data movement is managed explicitly. It is the pure `-stdpar` unified-memory approach that locks the code to NVIDIA; ifx in particular requires manual data regions anyway (§2.3). Adding OpenMP `target enter/exit data` directives around the time loop fixes residency on every vendor at once. The directives:

- are standard-conforming OpenMP that every relevant compiler accepts;
- compile away to nothing in CPU-only builds (no second code path);
- touch only ~5 places — field/memory-variable allocation, RK time-loop entry/exit, and the halo pack/unpack buffers — never the physics kernels themselves.

**2. Keeps the parts of MPI that cannot be replaced.** Parallel HDF5 (MPI-IO), CUDA-/GPU-aware halo exchange, profiling and debugging tooling (Nsight, Score-P, Vampir), and every HPC site's launch ecosystem assume MPI. Coarrays+GASNet would replace the *one* MPI feature WQ uses — halo exchange, the part that already works — while forfeiting the rest.

**3. Removes all bleeding-edge dependencies from the critical path.** Every component (MPI, OpenMP 5.x data directives, `do concurrent`) has two or more independent production implementations. Nothing load-bearing is a workshop paper.

### Alternatives considered and rejected

| Alternative | Verdict | Reason |
|---|---|---|
| **OpenACC + MPI** | Fallback only | More kernel-tuning control on NVIDIA, but a directive dialect with one real vendor; `do concurrent` reaches comparable performance on stencil codes with zero dialect. Reserve for individual hot kernels if profiling demands it. |
| **C++ rewrite on Kokkos/RAJA** (the SW4/EQSIM route) | Reject | Maximum portability and control, but a multi-year rewrite that abandons the team's Fortran expertise and the verified numerical core. Wildly disproportionate cost for a niche-strategy code. |
| **Full Berkeley stack now** (coarrays + Caffeine + GASNet, no MPI) | Reject for production; keep as experimental lane (§5 step 5) | Elegant endpoint, wrong timing: you become the production guinea pig for flang's PRIF lowering while losing nvfortran — the one compiler that makes the GPU port easy. |

### What survives from the original proposal unchanged

- **fpm, Julienne + Assert, and CI** — unambiguously right, independent of the parallelism decision (§2.4–§2.6).
- **The communication-abstraction module** (§5 step 3) — cheap insurance. If flang+Caffeine matures into a release-quality compiler shipping both PRIF multi-image support and `do concurrent` offload, swap the halo layer then: coarrays become a Phase-3 *simplification adopted on evidence*, not a Phase-1 gamble.

**One-line summary:** same modernization, different transport — modern ISO Fortran kernels and modern tooling now; MPI as the communication workhorse until the coarray ecosystem earns its place; directives (not unified memory) for cross-vendor GPU residency.

---

## 7. Detailed implementation plan (executes the §6 recommendation)

Assumes a small team (1–3 developers). Phases 0–1 are pure CPU work and carry no GPU-hardware dependency; Phase 4 runs as a time-boxed lane parallel to Phases 2–3. Every phase has explicit exit criteria — do not start the next phase until they hold, because each phase is the safety net for the one after it.

### Phase 0 — Foundation & safety net (weeks 1–3)

*Goal: make every later change verifiable and reversible. No numerics are touched.*

| # | Task | Detail | Exit criterion |
|---|---|---|---|
| 0.1 | **Fix git tracking** | Remove the stale `WQ/` gitlink, `git add` the full tree, commit (includes the 2026-06-11 Q-dispatch fixes, currently unversioned) | `git ls-files WQ/src` lists all sources; clean status |
| 0.2 | **Repo hygiene** | Delete `original_RHS_interior.f90`, `RHS_Interior.f90.bak`; consolidate the three `BoundaryConditions*.f90` into one module; document the even-rank constraint in the README | Single BC module compiles; no `.bak`/backup files in `src/` |
| 0.3 | **fpm dual build** | Add `fpm.toml` (apps: `waveqlab3d`, `pre_wql3d`; `dependencies: mpi`); restructure to fpm layout (`src/`, `app/`, `test/`); keep CMake in sync | `fpm build` and CMake both produce working binaries |
| 0.4 | **Julienne + Assert** | Dev-dependencies in `fpm.toml`. Unit tests, in value order: (a) attenuation — relaxation weights, tau spacing, unrelaxed-modulus correction per `response=` variant against analytic values; (b) friction — slip-weakening and rate-and-state point solvers vs. reference solutions; (c) **SBP operators** — accuracy order + discrete SBP property on random fields (this is the safety net for Phase 1.3); (d) PML damping profiles / `pml_damping_and_q`. Assert preconditions on namelist input (positive moduli, Q bounds, `npml` vs. block size) | `fpm test` green; SBP and attenuation suites exist **before** any kernel refactor |
| 0.5 | **CI** | GitHub Actions: matrix {gfortran-13, gfortran-14} × {fpm, CMake}; tiers: build → unit tests → MMS convergence-order check → 2–3 TPV regressions vs. truth solutions (tolerance-based seismogram comparison) | CI green on every push; failures block merge |
| 0.6 | **Performance baseline** | Profile one production-size TPV run on CPU: % time in RHS interior, halo exchange, fault interface, I/O; record wall-clock per simulated second | Baseline document committed — all later speedups are measured against this |

### Phase 1 — Kernel restructuring on CPU (months 1–3)

*Goal: make the code kernel-shaped. All work validated against Phase-0 tests; results must match the baseline within tolerance on the full TPV suite at every step. No GPU yet.*

| # | Task | Detail |
|---|---|---|
| 1.1 | **Communication abstraction** | New `comm` module wrapping `mpi3dcomm`: `comm_halo_exchange`, `comm_reduce_min`, `comm_gather_fault`, `comm_io_*`. Physics code never references MPI directly. Pure refactor — CI stays green throughout |
| 1.2 | **Dispatch hoisting** | Replace the per-point `if (M%anelastic_*)` chains in `RHS_Interior.f90` (lines ~255–262 and the curvilinear/PML variants) with a block-level `select case (response_type)` choosing specialized loop nests. Hoist the point-in-PML routing to region-level loops (interior vs. PML zones iterate separately — the zone bounds logic at lines ~206–213 already exists). Expect a measurable CPU speedup from branch elimination alone |
| 1.3 | **Stencil regeneration** | Recover or rewrite the generator for `JU_xJU_yJU_z6.f90` (fallback: regenerate from the SBP coefficient tables — the operators are standard 6th-order-interior/3rd-order-boundary SBP). Emit tight `do concurrent (k,j,i)` nests; no per-point calls to assumed-shape helpers; pure inlined expressions. **Validated by the 0.4(c) SBP unit tests and MMS convergence orders** |
| 1.4 | **Convert remaining loops** | RK stage updates, memory-variable (eta4–eta9) updates, Drucker–Prager return map, PML auxiliary updates → `do concurrent` with Fortran-2023 `reduce`/locality specifiers where needed. All point-local, low risk |
| 1.5 | **CPU-parallel validation** | Build with nvfortran `-stdpar=multicore` and ifx; verify thread-count scaling on one node; confirm no MPI calls outside the `comm` module (grep-enforced in CI) |

**Exit criteria:** full TPV suite matches baseline within tolerance on gfortran, ifx, nvfortran-multicore; measured single-node CPU speedup vs. 0.6 baseline; `RHS_Interior.f90` contains no per-point physics branching.

### Phase 2 — GPU bring-up (months 3–6)

*Goal: production single- and multi-GPU capability on NVIDIA, with the residency layer written portably.*

| # | Task | Detail |
|---|---|---|
| 2.1 | **Unified-memory bring-up** | nvfortran `-stdpar=gpu` (implicit data movement). Correctness first: MMS + TPV within tolerance. Profile with Nsight Systems — expect migration thrash; this step is for correctness, not speed |
| 2.2 | **Residency layer** | OpenMP `target enter data` at startup for the 9 field arrays, material arrays, memory variables, and PML auxiliaries; `target exit data` at shutdown; `target update` restricted to halo pack/unpack buffers and seismogram extraction. **~5 code sites, all in allocation/time-loop/comm code — physics kernels untouched.** This is the §6 portability move: the same directives serve ifx and Cray `ftn` |
| 2.3 | **Halo path on device** | Pack/unpack kernels as `do concurrent` on device into contiguous (pinned) buffers; use GPU-aware MPI where the site supports it, host-staged transfers otherwise — selected inside the `comm` module |
| 2.4 | **Fault & boundary kernels** | Profile first. The 2D fault-interface SAT solve is cheap but sits on the critical path each RK stage; if host↔device transfer of fault planes dominates, port the interface solve to device; otherwise keep on CPU and transfer only the fault-plane slices |
| 2.5 | **Cross-vendor check + CI lane** | Build (not necessarily tune) with ifx on Intel GPU and Cray `ftn` if accessible; add an nvfortran GPU lane to CI (self-hosted or container runner; nightly if per-push is impractical) |

**Exit criteria:** full TPV suite passes on GPU within tolerances; measured single-GPU speedup vs. one CPU node (report the number — do not promise one in advance); multi-GPU run on ≥8 GPUs with correct results; no unified-memory dependence remaining (builds run with `-gpu=nomanaged`).

### Phase 3 — Scale-out & production hardening (months 6–9)

| # | Task | Detail |
|---|---|---|
| 3.1 | **Comm/compute overlap** | Interior `do concurrent` kernels overlap with halo exchange (boundary-strip kernels wait; interior doesn't). Multi-node scaling study; weak/strong scaling plots |
| 3.2 | **Parallel I/O + restart** | Enable the dormant HDF5 path with parallel HDF5; checkpoint/restart for device-resident runs (fields + memory variables + fault state + RK phase) |
| 3.3 | **Performance regression CI** | Track per-kernel times in CI; alert on >10% regressions |
| 3.4 | **Release** | Modernized v1.0: documentation, versioned release, Zenodo DOI, updated citation. This is the adoption milestone the niche strategy (`IMPROVEMENT_ROADMAP.md` edges) builds on |

### Phase 4 — Experimental coarray lane (time-boxed, parallel to Phases 2–3)

*Goal: keep the Berkeley-stack option alive on evidence, not hope. Budget: ≤15% of team time.*

| # | Task | Detail |
|---|---|---|
| 4.1 | **Coarray `comm` backend** | Second implementation of the Phase-1.1 interface using coarray halo puts + `sync images` + `co_*` collectives; gfortran+OpenCoarrays first (most stable today) |
| 4.2 | **Caffeine/GASNet build** | flang + Caffeine + GASNet-EX on one test cluster; halo-exchange microbenchmarks vs. the MPI backend at 8–64 images; one TPV regression |
| 4.3 | **Go/no-go review at each flang release** | Adopt for production only when **all** hold: (a) one release-quality compiler ships both PRIF-based multi-image support and `do concurrent` offload; (b) halo performance ≥ MPI at target scale; (c) a parallel-I/O story exists (GASNet/MPI interop or HDF5 alternative). Otherwise: re-review in 6 months, lane stays frozen |

### Risk register

| Risk | Mitigation |
|---|---|
| Stencil generator is lost / unrecoverable | Regenerate from published SBP coefficient tables (operators are standard); 0.4(c) unit tests verify equivalence to the legacy file before it is deleted |
| nvfortran limits on derived-type components in `-stdpar` kernels | Flatten hot fields to plain arrays passed as arguments to kernel routines (also helps every other compiler) |
| Numerical differences from loop-order/reduction reordering | Tolerance-based (not bitwise) regression comparisons from day one; document tolerances per benchmark |
| GPU CI runner availability | Nightly self-hosted lane instead of per-push; per-push stays CPU-only |
| Two build systems drift | CI builds both on every push (0.5) |
| ifx/Cray tuning stalls (vendor-specific bugs) | NVIDIA lane is the supported product; other vendors are best-effort until a driving user appears |

### Effort summary

| Phase | Duration | Dominant skill |
|---|---|---|
| 0 | 2–3 weeks | Build/CI engineering |
| 1 | 2–3 months | Fortran refactoring + numerical verification |
| 2 | 2–3 months | GPU performance engineering |
| 3 | 2–3 months | HPC scaling + I/O |
| 4 | ≤15% time, parallel | Parallel-runtime experimentation |

Total: a modernized, GPU-capable, CI-protected WaveQLab3D in **roughly 9 months** for a small team — with the coarray future kept open at near-zero cost via the `comm` abstraction.

---

## Sources

- [Caffeine repository (BerkeleyLab)](https://github.com/BerkeleyLab/caffeine) · [Caffeine: parallel runtime for modern Fortran compilers (paper)](https://www.researchgate.net/publication/390305336_Caffeine_A_parallel_runtime_library_for_supporting_modern_Fortran_compilers) · [original Caffeine paper (IEEE/PAW-ATM)](https://ieeexplore.ieee.org/abstract/document/10027066) · [CASS software page](https://cass.community/software/caffeine.html)
- [Lowering and runtime support for Fortran's multi-image features via LLVM Flang, PRIF, and Caffeine (SC'25 workshops)](https://dl.acm.org/doi/10.1145/3731599.3767480)
- [Portability of Fortran's `do concurrent` on GPUs (arXiv)](https://arxiv.org/html/2408.07843)
- [Accelerating Fortran DO CONCURRENT with GPUs (NVIDIA blog)](https://developer.nvidia.com/blog/accelerating-fortran-do-concurrent-with-gpus-and-the-nvidia-hpc-sdk/)
- [Intel Fortran (ifx) 2025 release notes](https://www.intel.com/content/www/us/en/developer/articles/release-notes/fortran-compiler/2025.html) · [ifx DO CONCURRENT data-management discussion](https://community.intel.com/t5/Intel-Fortran-Compiler/Does-ifx-support-data-transfer-management-for-quot-DO-CONCURRENT/td-p/1679419)
- [Julienne repository (BerkeleyLab)](https://github.com/BerkeleyLab/julienne) · [Julienne + Assert == Correctness-Checking for Functional Fortran](https://escholarship.org/uc/item/4m0270sj) · [LBL software highlight](https://cs.lbl.gov/news-and-events/news/2025/software-highlight-julienne-and-assert-strengthen-fortran-code-reliability/)
- [GPU offloading in Fortran (Fortran-lang Discourse)](https://fortran-lang.discourse.group/t/gpu-offloading-in-fortran/9120)
