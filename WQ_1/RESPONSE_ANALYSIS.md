# WaveQLab3D — Response Types: Detailed Analysis

**Source files**: [elastic.f90](src/elastic.f90), [material.f90](src/material.f90), [RHS_Interior.f90](src/RHS_Interior.f90), [datatypes.f90](src/datatypes.f90), [withers_tables.f90](src/withers_tables.f90)  
**Generated**: 2026-06-10

---

## 1. Overview of All Valid Response Values

Set via `response=` in `&problem_list`. Validated at runtime by `domain.f90:130`:

| `response=` value | Physics | Mechanisms | Key feature |
|---|---|---|---|
| `elastic` | Pure elastic wave equation | — | No attenuation |
| `plastic` | Elastic + Drucker-Prager plasticity | — | Off-fault yield |
| `anelastic` / `low-pass` | GSLS attenuation (legacy) | N = 4 | Hand-fitted weights, has tau formula bug |
| `anelastic-Q` | GSLS attenuation (corrected) | N = 4 | NNLS-fitted weights, fixed tau formula |
| `anelastic-Q8` | GSLS attenuation (broadband) | N = 8 | Better constant-Q accuracy |
| `anelastic-Qf` / `frequency-Q-4M` | Frequency-dependent Q | N = 4 | Power-law Q(f), Withers tables |
| `frequency-Q-8M` | Frequency-dependent Q | N = 8 | Same, 8 mechanisms |
| `constant-Q-4M` | Configurable constant-Q | N = 4 | User fmin/fmax, NNLS solver |
| `constant-Q-8M` | Configurable constant-Q | N = 8 | Same, 8 mechanisms |

---

## 2. Elastic Response (`response='elastic'`)

### Physics

Solves the **first-order hyperbolic velocity-stress system** on a block-structured grid:

```
∂v/∂t  = (1/ρ) ∇·σ + f/ρ          (velocity update)
∂σ/∂t  = λ(∇·v) I + 2μ ε̇          (stress update)
```

Where:
- **v** = (v_x, v_y, v_z) — particle velocity
- **σ** = (σ_xx, σ_yy, σ_zz, σ_xy, σ_xz, σ_yz) — Cauchy stress tensor
- **λ, μ** — Lamé parameters, **ρ** — density

### Material Array Layout

`M%M(i,j,k,1:3)`:

| Index | Quantity | Formula from input |
|---|---|---|
| 1 | λ (Lamé first) | `ρ*(Vp² − 2*Vs²)` |
| 2 | μ (shear modulus) | `ρ*Vs²` |
| 3 | ρ (density) | direct |

### Material Profiles (Hardcoded by Problem)

Initialized in `material.f90:init_material` under `case('elastic')`:

| Problem name | Profile type |
|---|---|
| `default` | Homogeneous — uniform λ, μ, ρ from `rho_s_p` input |
| `LOH1` | 3-layer model at x=0: Layer of Heterogeneity benchmark |
| `LOH1_Harmonic` | LOH1 with harmonic-mean interface averaging |
| `TPV31` | Depth-gradient 7-layer profile (Vp/Vs/ρ all increase with depth) |
| `TPV32` | 9-layer strongly heterogeneous profile (near-surface soft sediments) |
| `TPV33` | Low-velocity zone flanked by standard rock (fault zone proxy) |
| `OKLAHOMA` | 12-layer crustal profile from Oklahoma seismicity region |

Material can also be loaded from file via `material_source='from_file'` using `init_material_from_file` (parallel MPI-IO binary format).

### Spatial Discretization

Entry point: `elastic.f90:set_rates_elastic` → `RHS_Interior.f90`:

**Two-stage computation**:
1. `RHS_Center` — interior points, uses optimized repeated-stencil loops
2. `RHS_Near_Boundaries` — boundary-adjacent points, uses general stencil per-point

**FD schemes** (`fd_type=` in input):

| `fd_type` | Orders available | Stencil type |
|---|---|---|
| `traditional` | 6 only | Centered SBP (Summation-by-Parts) |
| `upwind` | 2, 3, 4, 5, 6, 7, 8, 9 | Upwind-biased |
| `upwind_drp` | 3, 4, 5, 6, 7, 66, 679 | Upwind + Dispersion Relation Preserving |

Special orders `66` and `679` are multi-stencil DRP composites (mixed accuracy).

All schemes are implemented in `JU_xJU_yJU_z6.f90` as Jacobian-scaled operators that handle curvilinear coordinate transformations automatically.

### Memory Usage (Elastic Only)

- Fields **F**: 9 arrays (3 velocity + 6 stress), each of shape `(nq, nr, ns)`
- Material **M%M**: 3 arrays (λ, μ, ρ)
- No memory variables — purely instantaneous constitutive law

---

## 3. Anelastic Responses — Common Infrastructure

All anelastic variants share the same conceptual framework: **Generalized Standard Linear Solid (GSLS)** with N Maxwell elements in parallel.

### Physical Model

The GSLS stress-strain relation introduces **memory (history) effects**:

```
σ_effective(t) = σ_elastic(t) + Σ_k η_k(t)
```

Each memory variable `η_k` evolves by:

```
dη_k/dt = −η_k/τ_k + w_k · σ̇_elastic
```

Where:
- `τ_k` — relaxation time for mechanism k (seconds)
- `w_k` — dimensionless weight for mechanism k
- `σ̇_elastic` — rate of elastic stress

This gives **frequency-dependent effective moduli** — the medium is stiffer (faster) at high frequency and softer at low frequency, mimicking seismic attenuation.

### Q Definition

The quality factor Q measures energy loss per cycle:

```
Q(ω) ≈ 1 / (2 · Σ_k w_k / (1 + (ω·τ_k)²))
```

This function is what the different variants are fitting/targeting.

### Universal Q Ratio

In **all** variants, Q_P is hard-coded to:

```fortran
Q_P⁻¹ = 0.5 × Q_S⁻¹      (Q_P = 2 × Q_S)
```

This is a seismologically common approximation (P-wave attenuation half of S-wave).

### Q from Velocity

Q_S is always derived from the shear velocity field:

```fortran
Q_S⁻¹(i,j,k) = 1 / (c × sqrt(μ/ρ)) = 1 / (c × Vs(i,j,k))
```

The parameter `c` scales Q spatially (higher c = higher Q = less attenuation).

### Memory Variable Layout

All 6 stress components have their own set of memory variables:

| Variable | Stress component |
|---|---|
| `eta4(i,j,k,1:N)` | σ_xx |
| `eta5(i,j,k,1:N)` | σ_yy |
| `eta6(i,j,k,1:N)` | σ_zz |
| `eta7(i,j,k,1:N)` | σ_xy |
| `eta8(i,j,k,1:N)` | σ_xz |
| `eta9(i,j,k,1:N)` | σ_yz |

Rate arrays (`Deta4` through `Deta9`) are updated each RK stage alongside velocity and stress.

### Unrelaxed Moduli Correction

After weights and taus are set, all variants correct the elastic moduli to ensure the **relaxed (low-frequency) velocity matches the user-input velocity**:

```fortran
val_S = Σ_k w_k / ((ω_ref²·τ_k² + 1) × Q_S)
μ_unrelaxed = ρ × Vs² / (1 − val_S)
```

The stored moduli are then overwritten with the unrelaxed values. This is a critical step — without it, the simulated velocity would be lower than intended at the reference frequency.

---

## 4. `anelastic` / `low-pass` (Legacy, N=4)

**File**: [material.f90:11](src/material.f90#L11) — `init_anelastic_properties`  
**Namelist**: `&anelastic_list { c, weight_exp, fref }`

### Tau Formula (has a known bug)

```fortran
tau(k) = exp( ln(taumin) + (2k−1)/16 × ln(taumax/taumin) )
```

**Bug**: The denominator should be `2*N = 8` for 4 mechanisms, but is hard-coded as `16`. This spaces the relaxation times twice as densely as intended, reducing the effective Q bandwidth. `anelastic-Q` fixes this.

### Two Pre-fitted Weight Sets (selected by `weight_exp`)

| `weight_exp` | Weights (k=1..4) | Intended use |
|---|---|---|
| `0.0` (default) | 1.6126, 0.6255, 0.6382, 1.5969 | Standard attenuation |
| `0.6` | 0.0336, 0.6873, 0.8767, 1.5202 | Frequency-weighted attenuation |

Frequency band implied: τ_min ≈ 1/(2π×15) ≈ 0.0106 s, τ_max ≈ 200/(2π×0.08) ≈ 398 s — very wide but the buggy tau spacing means fewer independent mechanisms.

### Summary

- Kept for backward compatibility with old input files
- Use `anelastic-Q` instead for new simulations

---

## 5. `anelastic-Q` (Corrected, N=4)

**File**: [material.f90:145](src/material.f90#L145) — `init_anelastic_Q_properties`  
**Namelist**: `&anelastic_Q_list { c, fref }`

### Fixes over `anelastic`

1. **Correct tau formula**: denominator = `2*N = 8`
2. **NNLS-fitted weights**: optimized via Non-Negative Least Squares over [0.08, 15] Hz

### Relaxation Times

```
τ_min = 1/(2π×20) ≈ 0.007958 s    (covers up to 20 Hz)
τ_max = 1/(2π×0.05) ≈ 3.183 s     (covers down to 0.05 Hz)
τ_k = exp(ln(τ_min) + (2k−1)/(2N) × ln(τ_max/τ_min))
```

Band: **[0.05, 20] Hz** — suitable for most near-field strong-motion simulations.

### Pre-fitted Weights

```
w = [1.549360, 0.804277, 0.887718, 1.464160]
```

Fitted for constant-Q (γ=0), Q₀=1 (scaled by 1/Q_S at runtime).  
Accuracy: **max error < 18%, mean error < 0.3%** over [0.08, 15] Hz.

### Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `c` | 1.0 | Q_S = c × Vs (spatial scaling) |
| `fref` | 1.0 Hz | Reference freq for modulus correction |

---

## 6. `anelastic-Q8` (Broadband, N=8)

**File**: [material.f90:276](src/material.f90#L276) — `init_anelastic_Q8_properties`  
**Namelist**: `&anelastic_Q8_list { c, fref }`

Identical architecture to `anelastic-Q` but with **N=8 mechanisms**, providing better constant-Q approximation:

### Weights

```
w = [1.685770, 0.682533, 0.769700, 0.850033,
     0.916467, 0.971533, 1.067600, 1.528133]
```

Accuracy: **max error < 5%, mean error < 0.1%** over [0.05, 20] Hz.

### Cost vs Benefit

| | N=4 (`anelastic-Q`) | N=8 (`anelastic-Q8`) |
|---|---|---|
| Memory per grid point | 6 × 4 = 24 memory vars | 6 × 8 = 48 memory vars |
| Extra flops/step | 24 ODE updates | 48 ODE updates |
| Max Q error | < 18% | < 5% |
| Mean Q error | < 0.3% | < 0.1% |

Double the memory and compute cost for 3.6× better worst-case accuracy.

---

## 7. `anelastic-Qf` / `frequency-Q-4M` (Frequency-Dependent Q, N=4)

**File**: [material.f90:409](src/material.f90#L409) — `init_anelastic_Qf_properties`  
**Namelist**: `&anelastic_Qf_list { c, gamma, f_trans, fref }`

### Physics

Models **power-law Q**:

```
Q(f) = Q₀                      for f ≤ f_trans
Q(f) = Q₀ × (f / f_trans)^γ   for f > f_trans
```

Where γ ∈ [0.0, 0.9]. γ=0 reduces to constant Q. This empirically matches observations of frequency-dependent attenuation in the shallow crust.

### Withers Tables

Relaxation times and weights are looked up from pre-computed tables (**Withers et al., 2015**, BSSA):

```fortran
call get_relaxation_times_Qf(gamma, tau_Qf)
call get_withers_weights_Qf(gamma, Q_ref, weight_Qf)
```

Tables cover γ ∈ {0.0, 0.1, 0.2, ..., 0.9} with **N=4 mechanisms**.  
Band: approximately 0.4–100 Hz (τ_min ≈ 0.0032 s, τ_max ≈ 15.9 s for γ < 0.7).

### Parameters

| Parameter | Default | Range | Meaning |
|---|---|---|---|
| `c` | 1.0 | > 0 | Q_S scale factor |
| `gamma` | 0.0 | [0.0, 0.9] | Power-law exponent (clamped) |
| `f_trans` | 1.0 Hz | > 0 | Transition frequency |
| `fref` | 1.0 Hz | > 0 | Modulus correction reference |

---

## 8. `frequency-Q-8M` (Frequency-Dependent Q, N=8)

**File**: [material.f90:1715](src/material.f90#L1715) — `init_anelastic_Qf8_properties`  
**Namelist**: `&anelastic_Qf8_list { c, gamma, f_trans, fref }`

Same power-law Q model as `anelastic-Qf` but with **N=8 mechanisms** using the full 8-mechanism Withers tables:

```fortran
call get_relaxation_times(gamma, tau_Qf8)    ! 8-mechanism table
call get_withers_weights(gamma, 1.0, weight_Qf8)
```

Note: this routine also sets `M%anelastic_Qf = .true.` — a **flag aliasing bug**. The `anelastic-Qf8` state rides on the `anelastic_Qf` code path for RHS updates. This is intentional but implicit.

---

## 9. `constant-Q-4M` (Configurable Constant-Q, N=4)

**File**: [material.f90:1157](src/material.f90#L1157) — `init_const_Q_4M_properties`  
**Namelist**: `&constant_Q_4M_list { c, fmin, fmax, target_Q, weight_method, manual_weights, fref }`

The most flexible constant-Q variant. User specifies the target frequency band and Q value explicitly.

### Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `c` | 1.0 | Q_S = c / Vs (spatial scaling) |
| `fmin` | 0.05 Hz | Lower frequency bound of target band |
| `fmax` | 20.0 Hz | Upper frequency bound of target band |
| `target_Q` | 50.0 | Desired constant Q value |
| `weight_method` | `'nnls'` | Weight computation method |
| `manual_weights` | 0,0,0,0 | If all > 0: bypass auto-computation |
| `fref` | 1.0 Hz | Modulus correction reference frequency |

### Three Weight Computation Methods

#### `'nnls'` (default, recommended)
Solves a Non-Negative Least Squares problem at runtime, minimizing Q(f) error over `[fmin, fmax]`:

```
Minimize: Σ_f (Q_target − Q_approx(f, w))²
subject to: w_k ≥ 0 ∀k
```

Implementation: simplified active-set method, max 100 iterations, 30 log-spaced test frequencies.  
Uses custom Gaussian elimination solver (`solve_normal_equations`) with partial pivoting.

#### `'withers'`
Returns the first N weights from the 8-mechanism Withers table (γ=0), rescaled to `[fmin, fmax]`. Fast but less accurate than NNLS for narrow bands.

#### `'lookup'`
Same as `'withers'` but scales weights by `target_Q / 50.0`. Essentially a proportional rescaling of the Withers table — only appropriate as a rough approximation.

### Manual Override

If `manual_weights` are all positive, all auto-computation is bypassed. Useful for reproducing specific published results or expert tuning.

---

## 10. `constant-Q-8M` (Configurable Constant-Q, N=8)

**File**: [material.f90:1580](src/material.f90#L1580) — `init_const_Q_8M_properties`  
**Namelist**: `&constant_Q_8M_list { c, fmin, fmax, target_Q, weight_method, manual_weights, fref }`

Identical to `constant-Q-4M` but with N=8 mechanisms and 8 `manual_weights`.

**Note**: also sets `M%anelastic_const_Q_4M = .true.` — same flag aliasing pattern as `frequency-Q-8M`.

---

## 11. Withers Tables (`withers_tables.f90`)

**Reference**: Withers et al. (2015), *Bull. Seismol. Soc. Am.*, 105(6), 3129–3142.

Contains two pre-computed lookup tables for γ ∈ {0.0, …, 0.9}:

| Table | Usage | Condition |
|---|---|---|
| `W_HIGH_Q` | weights divided by Q | Q > 200 (high Q) |
| `A_COEF` / `B_COEF` | `w_k = a_k/Q² + b_k/Q` | 15 < Q < 200 |

Relaxation time bounds vary with γ:

| γ | τ_min (s) | τ_max (s) | Approx band |
|---|---|---|---|
| 0.0–0.6 | 0.0032 | 15.9155 | 0.01–50 Hz |
| 0.7–0.8 | 0.0066 | 3.979 | 0.04–24 Hz (narrower) |
| 0.9 | 0.0085 | 3.979 | 0.04–19 Hz |

---

## 12. Choosing the Right Response

```
Is attenuation needed?
├── No  → response='elastic'
│         Fastest, no overhead
│
└── Yes → Is Q frequency-dependent?
          ├── No (constant Q)
          │   ├── Quick test / legacy input?  → 'anelastic' (but has tau bug)
          │   ├── Standard accuracy needed?   → 'anelastic-Q' (N=4)
          │   ├── Better accuracy needed?     → 'anelastic-Q8' (N=8)
          │   └── Custom fmin/fmax/Q target?  → 'constant-Q-4M' or 'constant-Q-8M'
          │
          └── Yes (Q increases with frequency, power law)
              ├── gamma ∈ [0.0, 0.9], N=4?   → 'anelastic-Qf' / 'frequency-Q-4M'
              └── gamma ∈ [0.0, 0.9], N=8?   → 'frequency-Q-8M'
```

For **off-fault energy dissipation** in earthquake rupture dynamics, add:
- `response='plastic'` on top of any elastic configuration

---

## 13. Key Limitations & Issues

| Issue | Detail |
|---|---|
| **anelastic tau bug** | `init_anelastic_properties` uses denominator 16 instead of 2N=8 for 4 mechanisms, compressing tau spacing. Documented in `anelastic-Q` comments. |
| **Hard-coded Q_P/Q_S ratio** | Q_P = 2×Q_S in every variant. No way to set independent P-wave attenuation. |
| **Flag aliasing** | `constant-Q-8M` sets `anelastic_const_Q_4M=.true.`; `frequency-Q-8M` sets `anelastic_Qf=.true.`. The 8M variants piggyback on 4M flags — correct but fragile. |
| **NNLS implementation** | Active-set solver is a simplified approximation; convergence in 100 iterations not guaranteed for difficult Q targets. |
| **`lookup` method** | Just rescales Withers weights by `target_Q/50.0` — ignores the user's `fmin`/`fmax` entirely. Label is misleading. |
| **Memory overhead** | 8-mechanism variants use 6×8 = 48 memory variables per grid point. On a 200×200×400 grid this is ~1.9 GB of additional state. |
| **No spatially varying Q from file** | Q is always derived from velocity via `c/Vs`. No independent Q file input. |

---

## 14. Memory Budget (Per Block)

| Response | Extra arrays | Extra memory variables | Note |
|---|---|---|---|
| `elastic` | 0 | 0 | Baseline |
| `anelastic` | 2 (Qp_inv, Qs_inv) | 6×4×2 = 48 | eta + Deta |
| `anelastic-Q` | 2 | 6×4×2 = 48 | Same as anelastic |
| `anelastic-Q8` | 2 | 6×8×2 = 96 | 2× anelastic |
| `anelastic-Qf` | 2 | 6×4×2 = 48 | Same as anelastic |
| `frequency-Q-8M` | 4 (Qf + Qf8) | 6×8×2 + 6×8×2 = 192 | Both Qf and Qf8 sets |
| `constant-Q-4M` | 2 | 6×4×2 = 48 | Same as anelastic |
| `constant-Q-8M` | 4 (4M + 8M) | 96 + 96 = 192 | Both 4M and 8M sets |

---

*Analysis of source code at `/workspaces/Codespaces/WQ/src`. See also [CODEBASE_ANALYSIS.md](CODEBASE_ANALYSIS.md) for the broader architecture report.*
