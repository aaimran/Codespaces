# WaveQLab3D Improvement Roadmap — Reaching SOTA Parity and Building a Competitive Edge

*Date: 2026-06-11. Companion to `SOTA_COMPARISON.md` (the landscape survey), `CODEBASE_ANALYSIS.md` (architecture), `RESPONSE_ANALYSIS.md` (attenuation variants), and `SRC_BUG_REVIEW.md` (bug hunt & fixes).*

---

## Strategic framing

WaveQLab3D **cannot and should not chase SeisSol on its own turf** — unstructured intersecting fault networks at exascale represent roughly a decade of dedicated team effort (mesh infrastructure, code generation, LTS clustering, multi-physics). The realistic strategy is:

1. **Close a small number of *disqualifying* gaps** — things that currently prevent anyone from choosing WaveQLab3D regardless of its strengths.
2. **Double down on the defensible niche** where it is already ahead of every leader: provably stable SBP-SAT rupture dynamics with energy-stable PML and rich viscoelastic attenuation on curvilinear nonplanar faults.

The goal is not "trailing generalist with fewer features than SeisSol" but "clear leader in a niche the leaders have structurally abandoned."

---

## Part 1 — Parity gaps (ordered by how disqualifying they are)

### Gap 1: GPU execution — the non-negotiable one

**Problem.** Every current leader is GPU-resident (SeisSol via CUDA/HIP code generation; SW4/EQSIM via RAJA; AWP-ODC via CUDA; even the small FD3D_TSN). Pure-MPI Fortran caps WaveQLab3D 1–2 orders of magnitude below them in effective throughput on current machines. This single gap removes the code from consideration for any allocation-funded campaign.

**Why it's tractable.** Structured-grid FD is the *easiest* method class to port — regular memory access, no indirection, no gather/scatter. AWP-ODC and SW4 both made the transition; FD3D_TSN did it with a tiny team. WaveQLab3D does not need a CUDA rewrite: modern Fortran offload via **OpenMP `target` directives or `do concurrent`** is now well supported by NVIDIA (nvfortran) and AMD (flang/AOMP) compilers and preserves the single-source Fortran codebase.

**Prerequisite refactor — the code is not kernel-shaped today:**

- The RHS performs **per-point physics dispatch**: chains of `if (M%anelastic_*)` branches *inside* the triple spatial loop in `src/RHS_Interior.f90` (e.g., around lines 255–262 where five dispatch flags are tested per grid point before the PML call). These flags are constant for the entire run; they must be **hoisted out of the loops** into specialized loop variants (or a select-case at block level). This is worth doing for CPU performance alone (branch elimination, vectorization).
- The 23K-line generated stencil file `src/JU_xJU_yJU_z6.f90` must be **regenerated in a loop-fused, offload-friendly form** (collapse-able tight loops, no per-point subroutine calls with assumed-shape arrays). Recovering or rewriting the generator script is part of this task.
- Memory-variable updates (6 components × N mechanisms) should be restructured for coalesced access — the current `(x,y,z,component,mechanism)` access pattern needs auditing.
- MPI halo exchange must overlap with interior computation (communication/computation overlap) to scale on GPU nodes where compute is fast relative to network.

**Effort:** the dominant Phase-1 item; months for a small team, but each sub-step (dispatch hoisting → stencil regeneration → offload directives → halo overlap) delivers standalone value.

### Gap 2: Multi-fault geometry

**Problem.** The single conforming fault interface between two blocks is the biggest *scientific* limiter. Post-Kahramanmaraş (2023 Mw 7.8/7.7 doublet), multi-segment and branching rupture is the headline problem in the field, and it is the reason SeisSol dominates real-event studies.

**What's out of reach.** Fully unstructured intersecting fault networks are fundamentally incompatible with block-structured curvilinear grids. Do not attempt.

**Tractable intermediate.** Generalize the two-block interface machinery (`Interface_Condition.f90`, `Impose_Interface_Condition`) to **N blocks with multiple frictional interfaces**:

- Parallel and en-echelon fault segments (each a block-block interface) — covers stepovers, the most common multi-segment configuration.
- Branch geometries where faults meet at block corners/edges (requires careful SAT treatment at the junction lines — a genuine research contribution in the SBP-SAT setting, publishable in itself).
- The block-structured domain decomposition in `domain.f90`/`block.f90` already supports multiple blocks; the work is in the interface bookkeeping, the friction-coupling SATs at junctions, and proving stability.

This covers a large fraction of real events (stepovers, parallel strands, simple branches) without abandoning structured-grid efficiency.

**Effort:** Phase-2; the junction stability analysis is research-grade, the plumbing is engineering.

### Gap 3: Thermal pressurization (missing fault physics)

**Problem.** Thermal pressurization (TP) of pore fluids is now standard in SeisSol-class dynamic rupture studies and increasingly expected by reviewers for large-event modeling.

**Why it's cheap.** TP is a **1D diffusion solve per fault node** (temperature and pore pressure profiles normal to the fault), coupled to the friction law through effective normal stress. It touches only the fault interface code (`Interface_Condition.f90` and the rate-and-state friction path) — **no change to the wave solver, the attenuation machinery, or the PML**. The fault-local state-variable infrastructure (rate-and-state already carries per-node state) provides the pattern to follow.

**Effort:** medium-small; highest value per unit effort of any physics addition. Validate against SCEC TPV benchmark problems that include TP (TPV105 series).

### Gap 4: Attenuation generality — independent, file-loadable Q

**Problem.** Two hard-coded limits (documented in `RESPONSE_ANALYSIS.md`):

1. Q is always **derived from Vs** (`Q_S = c/Vs`) rather than independently specified.
2. `Q_P = 2·Q_S` is **hard-coded** — no independent P-wave attenuation.

Community velocity models (SCEC CVM, USGS models) ship independent Q_S and Q_P volumes; AWP-ODC and SW4 both consume them. Parity requires reading Q_S and Q_P as separate volumetric input fields.

**Why it's small.** The GSLS memory-variable machinery already supports spatially heterogeneous coefficients (the unrelaxed-modulus correction is computed per point). The work is: namelist/input plumbing for two more 3D fields, separating the P and S relaxation-weight computations, and regression-testing against the existing constant-Q cases. Also fix the misleadingly named `lookup` weight method (ignores fmin/fmax) or document it honestly.

**Effort:** small (days-to-weeks). Do early — it strengthens the attenuation edge (Part 2, Edge 2).

### Gap 5: Engineering credibility

**Problem.** No external group adopts a code in the current state, regardless of its numerics:

- **No CI** — tests are defined in CMake (`mpi_test1`, `mpi_test2`, …) and a strong validation suite exists (50+ TPV inputs, MMS, SCEC truth solutions) but nothing runs automatically.
- **Dead/duplicated code** — three `BoundaryConditions*.f90` variants compiled simultaneously (symbol-clash risk), `original_RHS_interior.f90` and `RHS_Interior.f90.bak` still in `src/`.
- **Flag aliasing hacks** — `constant-Q-8M` rides on the `anelastic_const_Q_4M` flag, `frequency-Q-8M` on `anelastic_Qf` (partially cleaned up in the 2026-06-11 bug-fix pass, but the pattern remains fragile).
- **Version control is broken** — `WQ/` is a stale gitlink: git tracks *nothing* inside it, so the recent critical Q-dispatch bug fixes exist only in the working tree. **This is the single most urgent item: one disk failure loses the fixes.**
- Undocumented constraints (even-MPI-rank requirement), manual memory management with no error-path cleanup, stale Doxygen.

**Fixes, in order:**
1. Repair git tracking (remove the gitlink, commit the tree) — hours.
2. CI (GitHub Actions: build matrix with gfortran + one vendor compiler, run the MMS + 2–3 TPV regression cases against truth solutions) — days. The test assets already exist; this is pure plumbing and it transforms the adoption story.
3. Delete dead code, consolidate BC modules, document the rank constraint — days.
4. Parallel HDF5 output + checkpoint/restart (the HDF5 path exists but is disabled) — needed for any long campaign on shared machines.

**Effort:** Phase-0; weeks total; zero research content; disproportionate credibility payoff.

### Explicitly *not* parity priorities

- **Local time stepping** — far less valuable on near-uniform structured grids than on graded tetrahedral meshes; SeisSol needs it because of mesh-induced timestep disparity that WaveQLab3D largely doesn't have.
- **Poroelasticity, tsunami/seismo-acoustic coupling** — SeisSol differentiators serving problem classes WaveQLab3D should not compete for.
- **Unstructured meshing** — abandoning block-structured grids forfeits the code's efficiency identity.

---

## Part 2 — Competitive edges (where to lead, not follow)

### Edge 1: Own the PML niche — and prove it

WaveQLab3D is nearly alone in having **energy-stable PML coupled to viscoelasticity *and* dynamic rupture**: `src/pml.f90` wires PML blocks into all six faces, and dedicated `apply_anelastic_point_pml` routines in `RHS_Interior.f90` update GSLS memory variables *inside* the PML. Every leader uses diffuse absorbers instead (SW4: supergrid; AWP: sponges; SeisSol: absorbing fluxes + domain extension) — all of which force substantially padded domains.

**The move:** a published head-to-head benchmark. Same rupture problem (e.g., TPV-class with attenuation), WaveQLab3D-with-PML versus sponge/supergrid-padded equivalents, demonstrating:
- equal seismogram accuracy at a **meaningfully smaller domain** (quantify the cost reduction — domains for rupture runs are commonly padded 2–3× per dimension);
- **long-time stability** in regimes where naive elastic PMLs are known to blow up (grazing incidence, surface waves, heterogeneity), leveraging the SBP-SAT energy estimates.

A concrete "X% cheaper for the same answer, provably stable" number is a marketable result no competitor can currently match, and it advertises the code's mathematical foundations.

### Edge 2: Best-in-class frequency-dependent attenuation for high-frequency ground motion

With the Q-dispatch bugs fixed (2026-06-11), the attenuation menu — 4- and 8-mechanism constant-Q (`anelastic-Q`, `anelastic-Q8`, `constant-Q-4M/8M`) plus power-law frequency-dependent Q(f) (`anelastic-Qf`, `frequency-Q-8M`) — is already **richer than SW4's and on par with AWP-ODC's**, and unlike AWP it sits in a spontaneous-rupture code.

The frontier in ground-motion simulation is pushing past 5–10 Hz (EQSIM's stated regime), exactly where Q(f) dominates the spectrum and where the Withers-lineage physics in this code matters. Adding independent, file-loadable Q_P/Q_S (Gap 4) makes WaveQLab3D **the most flexible attenuation engine available in any dynamic-rupture code** — a citable, defensible claim.

Target demonstration: a >10 Hz near-field ground-motion study from a dynamic source with Q(f), compared against constant-Q, showing the bias constant-Q introduces.

### Edge 3: Ensembles and inference, not hero runs

The field's direction (e.g., fused-ensemble Bayesian dynamic-rupture inference) is **hundreds-to-thousands of moderate forward runs**, not single hero simulations. For smooth single-fault problems, high-order structured FD costs far less per unit accuracy than DG on tetrahedra: no mesh generation per realization, no LTS clustering overhead, trivially predictable performance.

**Positioning:** "the fast, provably stable forward engine for rupture-parameter inference." In this market the geometry limitation (Gap 2) largely disappears — inference studies use simplified geometry by design. This pairs directly with the GPU port: the target is *1000 cheap runs on a handful of nodes*, not one exascale run. Deliverables:
- scripted parameter-sweep / ensemble driver (input-file templating, batched submission, automatic extraction of rupture metrics);
- a demonstration Bayesian inversion (e.g., recover stress-drop / friction parameters of a TPV-style event from synthetic seismograms).

### Edge 4: The rigor brand — dynamic solver for earthquake cycles (SEAS)

The SEAS community (Tandem, Thrase, Erickson-group SBP codes; the SCEC SEAS benchmarks) is built on **exactly WaveQLab3D's mathematical foundations** — SBP-SAT, provable stability of nonlinear friction coupling — and currently lacks mature 3D *fully dynamic* solvers with attenuation. Quasi-dynamic cycle codes approximate inertia; the known gap is coseismic fidelity.

**The move:** a hybrid workflow where a quasi-dynamic cycle code hands off interseismically-evolved fault state (stress, state variable) to WaveQLab3D for the fully dynamic coseismic phase, then returns the post-event state. Long-time stability guarantees — including the stable PML, since cycle-relevant domains must not accumulate boundary artifacts over many events — are precisely what this application demands. There is a ready-made benchmark community and publication venue (SEAS benchmark papers), and the rate-and-state aging-law implementation already matches SEAS conventions.

### Edge 5 (opportunistic): adjoint/differentiable capability

SBP discretizations are adjoint-consistent by construction; energy-stable PML has a well-defined adjoint. A discrete-adjoint capability for the velocity-stress + GSLS system would enable gradient-based dynamic source inversion and attenuation (Q) tomography from rupture data — territory none of the production rupture codes occupy. Higher risk/effort than Edges 1–4; pursue only if there is a driving science user.

---

## Part 3 — Sequencing

| Phase | Timescale | Items | Outcome |
|---|---|---|---|
| **0 — Credibility** | weeks | Fix git tracking (urgent — fixes are unversioned); CI on existing MMS+TPV suite; delete dead code / consolidate BC modules; independent Q_P/Q_S from file; document rank constraint | Adoptable, trustworthy code; attenuation edge sharpened |
| **1 — Performance & first edge** | months | Dispatch-hoisting refactor in `RHS_Interior.f90`; regenerate `JU_xJU_yJU_z6.f90` kernel-friendly; OpenMP target / `do concurrent` GPU offload; halo-overlap; thermal pressurization (+ TPV105 validation); **PML-vs-sponge benchmark paper** | Competitive throughput on GPU nodes; headline differentiator published |
| **2 — Market position** | 6–18 mo | Multi-interface faults (stepovers, simple branches) with junction stability analysis; ensemble/inference driver + demonstration inversion; SEAS hybrid coseismic handoff | Clear leadership in the defensible niche |

**Honest calibration:** Phases 0–1 are achievable for a small group (1–3 people). They do **not** reach SeisSol parity on geometry or multi-physics — and shouldn't try. They make WaveQLab3D the clear leader in: stable PML + viscoelastic dynamic rupture, frequency-dependent attenuation in a rupture code, cheap rigorous forward modeling for inference, and dynamic coseismic solves for the SEAS community.

---

## Cross-references

- `SOTA_COMPARISON.md` — the landscape this roadmap responds to
- `CODEBASE_ANALYSIS.md` §11 — issues/concerns underlying Gap 5
- `RESPONSE_ANALYSIS.md` — attenuation details underlying Gap 4 / Edge 2
- `SRC_BUG_REVIEW.md` — the 2026-06-11 Q-dispatch fixes assumed throughout
