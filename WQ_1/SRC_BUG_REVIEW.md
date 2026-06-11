# WaveQLab3D `src/` — Bug Hunt & Upgrade Review

**Date:** 2026-06-10
**Scope:** Full review of `WQ/src/*.f90` + `CMakeLists.txt`, focused on correctness bugs and upgrade opportunities.
**Companions:** `CODEBASE_ANALYSIS.md` (architecture), `RESPONSE_ANALYSIS.md` (Q-model physics), `FQ_PAPER_VS_CODE.md` (Withers-2015 comparison).

Two bugs from `FQ_PAPER_VS_CODE.md` have since been **fixed** in `material.f90`: the `f_trans` shift is now applied (`M%tau_Qf = M%tau_Qf / M%f_trans_Qf`, lines 511, 1794) and the Withers `w_k = N·λ_k` normalization is corrected (`M%weight_Qf / real(N_MECH, wp)`, line 517). The bugs below are new findings.

---

## 1. Critical bugs

> **Fix status (2026-06-11):** all four critical bugs are fixed in the working tree.
> C1 — `select case` in `block.f90`. C2 — `apply_anelastic_Qf_point_any` /
> `apply_const_Q_point_any` (module `anelastic_point`) called at all 32 interior
> dispatch sites in `JU_xJU_yJU_z6.f90`; `anelastic_point.f90` added to
> `src/CMakeLists.txt`. C3 — duplicate dispatch line deleted. C4 — `in_block_comm`
> guards in `domain.f90` `set_rates` and per-block guarded MMS print in `main.f90`.
> Build verified (gfortran, RELEASE). The serial regression tests mismatch the
> stored truth files by ~2-6 absolute, but this predates these fixes (verified by
> rebuilding without them: identical diff 1.9621861) — likely stale truth files /
> different compiler. C2 still needs physical validation against an analytic Q
> solution (see `FQ_tests.md`).

### C1. `anelastic-Qf` response is silently disabled — and leaves wrong moduli behind
**File:** `block.f90:166-194`

The response checks are a chain of independent `if/else` blocks. Three responses share the `anelastic_Qf` flag, and the *later* checks reset what the earlier ones set:

```fortran
if (trim(response) == 'anelastic-Qf') then
   call init_anelastic_Qf_properties(B%M, B%G, infile)   ! sets anelastic_Qf = .true.
else
   B%M%anelastic_Qf = .false.
end if
...
if (trim(response) == 'frequency-Q-4M') then
   call init_anelastic_Qf_properties(B%M, B%G, infile)
else
   B%M%anelastic_Qf = .false.        ! <-- wipes the flag set above for 'anelastic-Qf'
end if
```

For `response = 'anelastic-Qf'`: the flag is set true at line 166, then reset to `.false.` at line 188. Every dispatch site (`RHS_Interior.f90`, `fields.f90`) is gated on `M%anelastic_Qf`, so **no attenuation is ever applied**.

It is worse than "runs elastic": `init_anelastic_Qf_properties` has already **stiffened λ and μ to their unrelaxed values** (material.f90:519-546) on the assumption that the memory variables will relax them back. The simulation therefore runs purely elastically with wave speeds that are systematically too fast.

`frequency-Q-4M` and `frequency-Q-8M` happen to survive only because their init calls come after the resets.

**Fix:** replace the if/else chain with a single `select case(trim(response))` that initializes exactly one model, and initialize all flags `= .false.` once before it (they already default to `.false.` in `datatypes.f90`, so the `else` branches can simply be deleted).

### C2. Newer Q models are never applied at interior points
**Files:** `JU_xJU_yJU_z6.f90` (all `JJU_*_interior*` routines), `RHS_Interior.f90`

The generated interior-stencil routines dispatch only the three legacy models:

```fortran
if (M%anelastic)    ...
if (M%anelastic_Q)  ...
if (M%anelastic_Q8) ...
```

`grep` confirms **zero** occurrences of `anelastic_Qf`, `anelastic_const_Q_4M/8M`, or `Qf8` in the entire 1.3 MB file. Those dispatches exist only in `RHS_near_boundaries` (`RHS_Interior.f90`), which `cycle`s past interior points.

Consequence for `constant-Q-4M`, `constant-Q-8M`, `frequency-Q-4M`, `frequency-Q-8M` (and `anelastic-Qf` once C1 is fixed): memory-variable attenuation is applied **only in the thin shell of points near block faces** (≈6 points for the traditional scheme). The interior bulk evolves elastically — again with unrelaxed (stiffened) moduli from the init correction. Waves effectively see almost no attenuation regardless of the requested Q.

**Fix:** add the four missing dispatch calls next to the three legacy ones at every interior dispatch site in `JU_xJU_yJU_z6.f90` (the file is generated — fix the generator if it still exists, otherwise a scripted edit; the pattern is identical at each site). A verification run measuring amplitude decay against the analytic Q solution should be part of the fix (see `FQ_tests.md`).

### C3. Constant-Q applied twice at near-boundary points (traditional scheme)
**File:** `RHS_Interior.f90:257-258`

```fortran
if (M%anelastic_const_Q_4M) call apply_const_Q_4M_point_dispatch(F, M, G, x, y, z, Ux, Uy, Uz, DU)
if (M%anelastic_const_Q_4M) call apply_const_Q_4M_point_dispatch(F, M, G, x, y, z, Ux, Uy, Uz, DU)
```

Duplicated line in the `'traditional'` branch of `RHS_near_boundaries` only (the `upwind` branches call it once). The dispatch both increments `Deta*` and subtracts `sum(eta*)` from `DU`, so near-boundary points receive **double the attenuation contribution** relative to interior points for `constant-Q-4M` and `constant-Q-8M` (8M routes through the same dispatch). Combined with C2 this means: interior = no attenuation, boundary shell = 2× attenuation.

**Fix:** delete line 258.

### C4. Unguarded `D%B(2)` access — out-of-bounds in 1-block mode, undefined on block-1 ranks
**Files:** `domain.f90:738-767` (`set_rates`), `main.f90:107-108`

`set_rates` tests `D%B(1)%MT%use_moment_tensor` and `D%B(2)%MT%use_moment_tensor` with no guard:

- With `nblocks = 1`, `D%B` is allocated with size 1 → `D%B(2)` is an **out-of-bounds array access**.
- With 2 blocks and >1 rank, `init_block` is skipped for the block a rank doesn't own (`if (.not.in_block_comm(i)) cycle`), so `D%B(2)%MT%use_moment_tensor` is an **uninitialized logical** on block-1 ranks (and vice versa). Whether the body executes is luck of the memory.

Similarly `main.f90` prints `D%B(1)%sum ... D%B(2)%sum` unconditionally when MMS is on — out-of-bounds in 1-block mode, and undefined for ranks not in that block's communicator.

**Fix:** guard with `if (in_block_comm(i))` / loop `i = 1, D%nblocks` (matching how `enforce_bound_iface_conditions` and `update_fields` already do it), and reduce the MMS norm across ranks before printing on the master.

---

## 2. High-severity bugs

### H1. `type_of_mesh` has no default value
**File:** `domain.f90:75,101-121`

Every other `problem_list` variable gets a default before the namelist read; `type_of_mesh` does not. If the user's input file omits it, an **uninitialized string** is compared in `init_block` (`'cartesian'` vs `'curvilinear'`) — with `mesh_source='compute'` neither branch may run, leaving the grid unallocated and crashing later with no useful message. Set `type_of_mesh = 'cartesian'` (and validate it like `response` is validated).

### H2. Stress-component output files are mislabeled
**File:** `seismogram.f90:334-353`

Field storage order (from `RHS_Interior` rates) is `4=sxx, 5=syy, 6=szz, 7=sxy, 8=sxz, 9=syz`. The output naming uses `4=sxx, 5=sxy, 6=sxz, 7=syy, 8=syz, 9=szz`. Five of the six stress files contain a different component than their filename says (e.g. `*_sxy.dat` actually holds σyy). Anyone post-processing these files gets silently wrong science.

### H3. Block-2 field output files are opened but never written
**File:** `seismogram.f90:361-365` vs `write_seismogram` (398-405)

`init_seismogram` opens `file_unit_block2(1:9)` when `output_fields_block2` is set, but `write_seismogram` only has an `output_fields_block1 .and. block_num==1` branch. Block-2 files are created empty and closed at the end. Add the symmetric branch (or generalize to one per-block array).

### H4. Fault `Svel` output actually contains `Uhat` components
**File:** `domain.f90:611-614` (`write_hat_output`)

```fortran
call write_hats(D%fault%Uhat_pluspres(...), D%fault%Vhat_pluspres(...), &
                D%fault%Uhat_pluspres(mr1:pr1,ms1:ps1,1:3), &   ! <-- passed as "Svel"
                D%fault%time_rup(...), D%fault)
```

The third argument is written to `*_interface.Svel`, but it is `Uhat(:,:,1:3)`, not the slip-velocity array (`fault%Svel` / `I%Svel` exists and is maintained elsewhere). If intentional (hat-velocity output), the file name is misleading; if not, slip-velocity output is wrong. Either way it needs resolving.

### H5. Legacy `anelastic` response: unsupported `weight_exp` leaves τ = 0 → division by zero
**File:** `material.f90:86-110`

`tau` and `weight` are populated only for `weight_exp ≈ 0` and `weight_exp ≈ 0.6`. Any other value (e.g. 0.3) sails through with `tau = 0`, and the memory-variable update divides by `M%tau(i)` → immediate Inf/NaN with no diagnostic. Validate and `stop` with a clear message.

### H6. `nt` input parameter is read and silently ignored
**File:** `domain.f90:88,288`

`nt` is in the `problem_list` namelist but `D%nt = floor(D%t_final/dtmin)` unconditionally. Users who set `nt` get a different number of steps than requested. Remove it from the namelist or honor it (`if (nt > 0) D%nt = nt`).

### H7. Missing namelists are silently ignored everywhere
All namelist reads use `if (stat > 0) stop ...`. A *missing* namelist (or a typo'd group name) returns `stat < 0` (EOF) and silently leaves defaults in place. For physics-critical groups (`anelastic_Qf_list`, `moment_list`, …) this hides input-file mistakes. At minimum warn on `stat < 0`; ideally track which groups are required per response.

### H8. `find_neighbor3d` falls back to an arbitrary rank
**File:** `mpi3d_interface.f90:376-379`

If no opposite-side process with matching interface coordinates is found, the code assigns `rank_neighbor = ranks(min(comm_size, 2*n(1)*n(2)))` — a random member of the communicator. Subsequent `MPI_SendRecv` will exchange data with the wrong rank, producing corrupt interface values or deadlock rather than a clean error. This can trigger whenever the two blocks' decompositions don't align (which is never validated — the `@todo` in `iface.f90:30`). Replace the fallback with `call error(...)`.

---

## 3. Medium / latent issues

| # | Location | Issue |
|---|----------|-------|
| M1 | `domain.f90:191` | Time-step loop reads `btp(i)` for `i = 1..max(1,nblocks)` **before** `nblocks` is validated (line 208). `nblocks = 3` in the input reads `btp(3)` out of bounds. Validate first. |
| M2 | `domain.f90:240-241` | `in_block_comm(1) = rank < (nprocs+1)/2`, `in_block_comm(2) = rank >= nprocs/2`: for odd `nprocs` the middle rank joins **both** blocks, and the loop-carried `coord`/`cart_size` used for `new_interface` then belong to block 2 while the rank also represents block 1. Currently latent only because `start_mpi` rejects odd `nprocs` (mpi3dbasic.f90:38) — but that guard contradicts 1-block mode, where odd counts would be fine. Make the split explicit and remove the loop-carried-variable dependence. |
| M3 | `time_step.f90:26` | RK order hard-coded `order = 4` with comment "read from input file". Williamson(3,3) and the other schemes are unreachable. Plumb through the input file. |
| M4 | `main.f90:84-89` | `cpu_time` used for per-step timing in an MPI program (measures CPU, not wall, and only on the master). Use `MPI_Wtime` (already wrapped as `time_elapsed()`). The local `dt = 0.1_wp` parameter is dead. |
| M5 | `moment_tensor.f90:894-929` | `set_moment_tensor` loops over **every grid point** for every tensor at every RK stage, calling `source_time` (constant per tensor/stage) and `singular_source` (compact support, a few points wide) everywhere. Hoist `source_time` and restrict the spatial loop to the support of the discrete delta — this is a large constant cost on big grids. |
| M6 | `moment_tensor.f90:36-68,111-125` | Tensor-list parsing `read(infile,'(a)') temp` without `iostat`: missing `!---begin:tensor_list*---` markers crash with a raw runtime error. `n_mom` is undefined if `block_id` ∉ {1,2}. |
| M7 | `seismogram.f90:400-402` | Leftover debug `print *, S%file_unit_block1(n)` in the per-step output path; also full-field ASCII dump `write(...,*) F(:,:,:,n)` includes **ghost nodes** and is enormous — should slice `mq:pq,...` and be binary. |
| M8 | `seismogram.f90` | `stride_fields` is read but never used (fields written every `w_stride` instead). `i_phys/j_phys/k_phys` never deallocated. Station write format `4f15.10` silently prints `***` for values ≥ 10⁴; use `es` format. |
| M9 | `plastic.f90:184` | `Pf = 9.8d0*y_ijk` (pore pressure) computed in `plastic_flow2` but **never used** — the yield stress uses total stress `sigma`, not effective stress `sigma + Pf`. Either the SCEC effective-stress term is missing, or the line is dead; verify against the TPV26/27 spec. Material constants (`c`, `nu`, `Tv`) are hard-coded per problem string. |
| M10 | `Interface_Condition.f90:1216-1237` | Bisection solver: sign tests written as `fv*fl > tol` (product vs 1e-12 — wrong near roots); non-convergence only `print`s and continues with the last iterate. Should at least flag/error. |
| M11 | `mpi3dcomm.f90:236-248` | `block3d_qr` built with `MPI_Type_struct(C%nb, ...)` while the metadata arrays are sized `C%lnq` — works only because `nb ≤ lnq`; confusing and fragile. Whole file uses deprecated MPI-1 `MPI_Type_struct` with default-integer byte displacements (overflow for > 2 GiB local arrays). Replace with `MPI_Type_create_subarray`. |
| M12 | `fault_output.f90` + `domain.f90` | The same 8 fault files are opened by **two different communicators** (side-1 and side-2), with handles 1,3,4,5,6,7,8 written only by side 1 and handle 2 only by side 2. Works, but any future change that makes writes collective across the "other" comm deadlocks. Worth restructuring around a single fault communicator. |
| M13 | `common.f90` | `sp = selected_real_kind(9,49)` — 9 significant digits already forces a **double-precision** kind; the "short precision" name is misleading (real single is `selected_real_kind(6,37)`). |
| M14 | `CMakeLists.txt` | `fault_output.f90` listed twice in `WQL3D_SRC`; `-O5` is not a real gfortran level (clamped to `-O3`); `cmake_minimum_required(VERSION 3.5)` is below modern CMake's deprecation floor. |
| M15 | `src/` hygiene | Dead files that still pattern-match builds/greps: `RHS_Interior.f90.bak`, `original_RHS_interior.f90`, `BoundaryConditions2.f90`, `BoundaryConditions3.f90` (the latter two define the **same module name** `BoundaryConditions` — adding either to the build breaks it). Move to an `attic/` or delete; git history preserves them. |
| M16 | `grid.f90` | 12 leftover debug `print *` statements in production paths. |
| M17 | memory | All `eta*`/`Deta*` arrays are allocated **with ghost nodes** (`ghost_nodes=.true.`) but are never halo-exchanged. For `*-8M` that is 12 arrays × 8 mechanisms of wasted ghost storage; allocate interior-only. |
| M18 | `fields.f90` / `iface.f90` | `DF`, `DS`, `DW` initialized to `1.0e40` sentinels; correctness relies on every RK scheme having `A(1) = 0` so stage 1 zeroes them. True for all current schemes, but a fragile invariant — worth an explicit zeroing at step start or a comment in `init_RK`. |

---

## 4. Suggested upgrades

**Correctness / maintainability**

1. **Single attenuation framework.** The 8 Q-variants duplicate ~12 arrays + 6 scalars + a flag each in `block_material`, with 6 near-identical code blocks in `fields.f90` and per-variant `apply_*_point` routines in `RHS_Interior.f90`. Replace with one set of arrays `eta(:,:,:,1:n_mech,1:6)`, one `n_mech`, one `tau(:)/weight(:)`, and one enum for the weight/τ recipe. This single change would have made bugs C1–C3 impossible.
2. **Select-case response dispatch** in `block.f90` (fixes C1 structurally) and validation of every response-specific namelist.
3. **Stop generating/committing the 1.3 MB `JU_xJU_yJU_z6.f90`** as opaque text: either commit the generator script and run it at build time, or factor the per-point physics (rate assembly + anelastic dispatch) into one inlined routine so it exists in exactly one place (fixes C2 structurally).
4. **Error handling:** use the existing `error()` (MPI_Abort) wrapper instead of bare `stop` in MPI-collective contexts — a lone `stop` on one rank hangs the rest of the job on most MPI implementations.
5. **Input validation pass at startup** (master rank): required namelists present, `nblocks` range, `type_of_mesh`/`fd_type`/`order` combinations supported, station/tensor list markers present. Most current failure modes are crashes deep in init.

**Performance**

6. Compact-support moment-tensor application (M5) — likely the cheapest large win for source-driven runs.
7. Use the already-computed work ratios `ratio1/ratio2` for the block rank split (`domain.f90:233-241` computes them and then hard-codes 50/50); blocks of unequal size currently idle half the machine.
8. Non-blocking halo exchange (`MPI_Isend/Irecv` + overlap with `RHS_Center`) — the structure (separate center/boundary RHS) is already perfect for it.
9. Hybrid parallelism (OpenMP on the triple loops in `RHS_near_boundaries` / interior stencils) for modern many-core nodes.
10. Binary (MPI-IO or HDF5 — scaffolding exists, currently commented out) field/plane output instead of per-rank ASCII.

**Modernization**

11. `use mpi_f08` instead of `use mpi`; replace MPI-1 derived-type constructors with `MPI_Type_create_subarray` (also fixes M11).
12. Read RK order from input (M3); document supported `fd_type`/`order` pairs (the `66`/`679` magic orders deserve names).
13. CI pipeline running the existing `ctest` suite (`mpi_test1/2`, `premesh`, `serial`) plus an MMS convergence check and a Q-decay check; none of the bugs above survive a Q-amplitude regression test.
14. `cmake_minimum_required(VERSION 3.20)`, dedupe source list, prefer `-O3` and add an explicit `-fcheck=bounds` debug preset.

---

## 5. Priority fix order

1. **C1 + C2 + C3** — anything using the new Q models today is producing wrong physics (no/double attenuation, stiffened moduli). One PR: select-case dispatch in `block.f90`, add 4 dispatches to interior stencils, delete duplicated line; verify with Q-decay test.
2. **C4, H1** — crash/UB bugs triggered by 1-block mode or omitted input keys.
3. **H2, H3, H4** — output correctness (mislabeled stresses is the nastiest: silently wrong post-processing).
4. **H5–H8** — input robustness and MPI neighbor fallback.
5. Medium items opportunistically; upgrades 1–3 as the structural follow-up.
