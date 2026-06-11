# FQ.pdf vs Implemented Responses — Detailed Comparison

**Paper**: Withers, K. B., Olsen, K. B., & Day, S. M. (2015).  
*Memory-Efficient Simulation of Frequency-Dependent Q.*  
Bull. Seismol. Soc. Am., 105(6), 3129–3142. doi: 10.1785/0120150020

**Source files examined**: [withers_tables.f90](src/withers_tables.f90), [material.f90](src/material.f90), [JU_xJU_yJU_z6.f90](src/JU_xJU_yJU_z6.f90)

---

## 1. The Paper in Brief

### Physical Model

The paper targets a **piecewise power-law Q(f)**:

```
Q(f) = Q₀                       for  f ≤ f_T
Q(f) = Q₀ × (f / f_T)^γ        for  f > f_T
```

where γ ∈ [0.0, 0.9] is the power-law exponent (γ=0 = constant Q).

This is approximated by N=8 memory variables obeying:

```
τ_k · dξ_k/dt = −ξ_k + λ_k · ε(t)          (eq. 4)
σ(t) = M_u · ε(t) − Σ_k ξ_k(t)              (eq. 3)
```

The Q function from this model (low-loss approximation, eq. 8):

```
Q⁻¹(ω) ≈ Σ_k  λ_k · ωτ_k / (1 + (ωτ_k)²)
```

### Two-Regime Fitting

| Q range | Method | Formula |
|---|---|---|
| Q > 200 (high Q) | Linear NNLS (low-loss approximation) | `λ_k = w_k / Q₀` (scale from Q₀=1 fit) |
| 15 ≤ Q ≤ 200 (low Q) | Nonlinear conjugate-gradient (effective Q, eq. 12) | `λ_k = a_k/Q² + b_k/Q` (eq. 16) |

### Relaxation Times (Paper eq. 15)

```
ln(τ_k) = ln(τ_m) + (2k − 1)/16 × ln(τ_M / τ_m)
```

where τ_m, τ_M are the band limits (both γ-dependent, see Table 1).

---

## 2. Table-by-Table Verification

### Table 1 (High-Q weights) — Code: `W_HIGH_Q`

Cross-checking all 10 γ values against the paper. **Exact match.**

Sample (γ=0.0):

| k | Paper w_k | Code W_HIGH_Q | Match |
|---|---|---|---|
| 1 | 0.8867 | 0.8867 | ✓ |
| 2 | 0.8323 | 0.8323 | ✓ |
| 3 | 0.5615 | 0.5615 | ✓ |
| 4 | 0.8110 | 0.8110 | ✓ |
| 5 | 0.4641 | 0.4641 | ✓ |
| 6 | 1.0440 | 1.0440 | ✓ |
| 7 | 0.0423 | 0.0423 | ✓ |
| 8 | 1.7275 | 1.7275 | ✓ |

Relaxation time bounds (τ_m, τ_M) are also reproduced exactly from the table:

| γ | Paper τ_m | Code TAU_MIN | Paper τ_M | Code TAU_MAX |
|---|---|---|---|---|
| 0.0–0.6 | 0.0032 | 0.0032 | 15.9155 | 15.9155 | 
| 0.7–0.8 | 0.0066 | 0.0066 | 3.9789 | 3.9789 |
| 0.9 | 0.0085 | 0.0085 | 3.9789 | 3.9789 |

### Table 2 (Low-Q coefficients) — Code: `A_COEF`, `B_COEF`

**Exact match** for all 8 mechanisms × 10 γ values.

Sample (γ=0.0):

| k | Paper a_k | Code A_COEF | Paper b_k | Code B_COEF |
|---|---|---|---|---|
| 1 | −27.5 | −27.50 | 7.41 | 7.410 | ✓ |
| 5 | 14.6 | 14.60 | 3.88 | 3.880 | ✓ |
| 8 | −82.8 | −82.80 | 13.19 | 13.190 | ✓ |

---

## 3. Tau Formula

**Paper eq. 15:**
```
ln(τ_k) = ln(τ_m) + (2k−1)/16 × ln(τ_M/τ_m)
```

**Code (`withers_tables.f90:167`):**
```fortran
tau(k) = exp(log(taumin) + real(2*k-1, wp) / 16.0_wp * (log(taumax) - log(taumin)))
```

**Status: Exact match** ✓

The denominator **16 = 2×N** where N=8 is correct for the Withers 8-mechanism table. *(Note: the `anelastic` legacy variant incorrectly uses 16 for N=4 mechanisms — see [RESPONSE_ANALYSIS.md](RESPONSE_ANALYSIS.md).)*

---

## 4. Memory Variable ODE — Code vs Paper

### Paper (eq. 17, from Day & Bradley 2001)

For 3D stress tensor with separate shear (Q_s) and bulk (Q_p):

```
τ_k · dξ_ij/dt = −ξ_ij + λ_k [ 2μ_u Q_s⁻¹ ε̇_ij + (κ_u Q_p⁻¹ − ⅔ μ_u Q_s⁻¹) δ_ij ε̇_kk ]
```

### Code (`JU_xJU_yJU_z6.f90:16326-16349` for `anelastic`)

For the diagonal components (σ_xx):
```fortran
Deta4(x,y,z,i) += ( weight(i)*2μ*Qs_inv*ε̇_xx
                   + weight(i)*((λ+2μ)*Qp_inv − 2μ*Qs_inv)*tr
                   − eta4(x,y,z,i) ) / tau(i)
```

For the off-diagonal shear (σ_xy):
```fortran
Deta7(x,y,z,i) += (weight(i)*μ*Qs_inv*(ε̇_yx + ε̇_xy) − eta7(x,y,z,i)) / tau(i)
```

where `tr = ε̇_xx + ε̇_yy + ε̇_zz`.

**Mapping to paper notation:**

| Code term | Paper term |
|---|---|
| `weight(i) * 2μ * Qs_inv` | `λ_k * 2μ_u * Q_s⁻¹` (deviatoric shear forcing) |
| `weight(i) * ((λ+2μ)*Qp_inv − 2μ*Qs_inv)` | `λ_k * (κ_u Q_p⁻¹ − ⅔ μ_u Q_s⁻¹)` (volumetric forcing) |
| `eta4 / tau(i)` | `ξ_ij / τ_k` (relaxation) |

**Status: Structurally correct** ✓ — the code correctly separates P and S attenuation channels. The combination `weight(i) * Qs_inv` plays the role of `λ_k`.

### Stress update

```fortran
DF(4) -= (eta4(x,y,z,1) + eta4(x,y,z,2) + eta4(x,y,z,3) + eta4(x,y,z,4))
```

This implements `∂σ_xx/∂t = ... − Σ_k ξ_k`, matching paper eq. 18. ✓

---

## 5. Q Normalization — Conventional vs Coarse-Grained

This is the most important structural difference between the paper and the code.

### Paper's approach: Coarse-grained

The Withers code distributes memory variables across space:
- One relaxation time per stress node (cycling period = 2)
- For N=8 mechanisms in 3D, each coarse-grained cell has 23 stress nodes
- The coarse-grained scaling: **w_k = N_cell × λ_k** (paper eq. 43 from Day 1998)
- Paper Table 1 caption: *"wk = Nλk, in which N = 8"*

The factor N=8 in the caption refers to the number of mechanisms, not the cell volume. The table weights already absorb the normalization needed for the coarse-grained approach.

### Code's approach: Conventional (all N at every node)

The code stores and updates **all N memory variables at every grid node**. In the conventional approach (non-coarse-grained), the Q formula is:

```
Q⁻¹(ω) = Σ_k  λ_k · f(ωτ_k)
```

where λ_k are the bare relaxation weights from the paper's fitting.

The paper fits λ_k such that this sum ≈ 1/Q₀ (with Q₀=1). The Table 1 values (labeled w_k) are explicitly stated to be w_k = N·λ_k = 8·λ_k.

### Implication for the code

The code uses `W_HIGH_Q(k, gamma) / Q_S` directly as its effective λ_k in the RHS. This means:

```
Effective λ_k in code = W_HIGH_Q(k, gamma) / Q_S = w_k / Q_S = 8·λ_k_paper / Q_S
```

For the code's Q to match the intended Q_S:

```
Q⁻¹_code(ω) = Σ_k (8·λ_k_paper / Q_S) · f(ωτ_k)
             = (8 / Q_S) · Σ_k λ_k_paper · f(ωτ_k)
             ≈ 8 / Q_S   [since Σ_k λ_k · f ≈ 1/Q₀=1]
```

This would give an effective Q ≈ Q_S/8 — **8× more attenuation than the Q_S specified by `c`**.

**However**, this is only the case for variants using `get_withers_weights` directly (`frequency-Q-8M` and `anelastic-Qf8`). The variants with independently fitted weights (`anelastic-Q`, `anelastic-Q8`, `constant-Q-4M`, `constant-Q-8M`) use NNLS optimization that finds weights directly fitting Q, without the wk = N·λk ambiguity.

| Response | Weight source | N factor concern |
|---|---|---|
| `anelastic` | Hand-fitted, Q0-scale unclear | Likely affected |
| `anelastic-Q` | NNLS over [0.05, 20 Hz] — self-consistent | No issue |
| `anelastic-Q8` | NNLS over [0.05, 20 Hz] — self-consistent | No issue |
| `anelastic-Qf` / `frequency-Q-4M` | First 4 of Withers Table 1 (w_k = N·λ_k) | **Potentially 8× error** |
| `frequency-Q-8M` | Full Withers Table 1 (w_k = N·λ_k) | **Potentially 8× error** |
| `constant-Q-4M` | NNLS or Withers method | NNLS: No issue; Withers/lookup: potentially affected |
| `constant-Q-8M` | NNLS or Withers method | Same as above |

---

## 6. `f_trans` — Not Applied to Relaxation Times

### Paper

The transition frequency f_T is implemented by **rescaling τ_k**:

> "This frequency shift is achieved by dividing the relaxation times by the desired factor; the weights remain unchanged."

So for a target f_T ≠ 1 Hz (the base fitting frequency):
```
τ_k(f_T) = τ_k(1 Hz) / f_T
```

### Code

In `init_anelastic_Qf_properties` ([material.f90:510](src/material.f90#L510)):

```fortran
M%f_trans_Qf = f_trans          ! stored ...
call get_relaxation_times_Qf(M%gamma_Qf, M%tau_Qf)  ! ... but NOT passed here
call get_withers_weights_Qf(M%gamma_Qf, 1.0_wp, M%weight_Qf)
```

`f_trans` is read from the namelist and stored in `M%f_trans_Qf`, but **never used** in the tau or weight computation. The same issue exists in `init_anelastic_Qf8_properties`.

**Result: The transition frequency is always f_T = 1 Hz regardless of the user's `f_trans` setting. This is a silent bug.**

The fix would be to divide the computed relaxation times by `f_trans`:

```fortran
call get_relaxation_times_Qf(M%gamma_Qf, M%tau_Qf)
M%tau_Qf = M%tau_Qf / M%f_trans_Qf     ! Apply frequency shift
```

---

## 7. 4-Mechanism Extraction from 8-Mechanism Table

### Paper

The paper exclusively uses N=8 mechanisms. No N=4 variant is discussed.

### Code

`get_withers_weights_Qf` and `get_relaxation_times_Qf` extract the **first 4** of the 8 mechanisms from the Withers tables. This is used by `anelastic-Qf` / `frequency-Q-4M`.

```fortran
subroutine get_relaxation_times_Qf(gamma, tau_Qf)
   call get_relaxation_times(gamma, tau_full)
   tau_Qf(1:4) = tau_full(1:4)   ! first 4 only
end subroutine
```

This is an approximation with no theoretical justification from the paper. The first 4 mechanisms span the low-frequency half of the relaxation band. For γ=0.0 (constant Q), the first 4 Withers weights are:

```
w₁=0.8867, w₂=0.8323, w₃=0.5615, w₄=0.8110
```

The last 4 (which are dropped):
```
w₅=0.4641, w₆=1.0440, w₇=0.0423, w₈=1.7275
```

Mechanism 8 (weight 1.7275) is the **largest**, covering the highest frequencies. Dropping it means the 4M variant is particularly inaccurate at high frequencies (near f_max of the band).

---

## 8. Smooth Transition Region

### Paper

To help the NNLS fitting near f_T, the paper applies a smooth Q(f) transition between 0.8 f_T and 1.2 f_T (power γ/2 in the transition zone), rather than a sharp step.

### Code

The code uses a **sharp step** at f_trans. The parameter `f_trans` is stored but (as noted above) not applied at all. There is no smoothing.

---

## 9. Independent Q_P and Q_S

### Paper

The paper treats Q_S and Q_P as independent. For the validation tests it uses Q_S = Q_P = 50, but the framework supports different Q values for each. The memory variable ODE eq. 17 carries both Q_s⁻¹ and Q_p⁻¹ terms explicitly.

### Code

**Q_P is hard-coded as 2×Q_S everywhere:**

```fortran
M%Qp_inv_Qf(i,j,k) = 0.5_wp * M%Qs_inv_Qf(i,j,k)    ! Q_P = 2*Q_S always
```

Users cannot set an independent P-wave attenuation. The Q_P = 2Q_S approximation is physically reasonable for the crust but is a constraint the paper does not impose.

---

## 10. Gamma Interpolation (Code Extension Beyond Paper)

### Paper

The paper provides discrete tables for γ ∈ {0.0, 0.1, 0.2, ..., 0.9}.

### Code

The code adds **linear interpolation in γ** between tabulated values:

```fortran
subroutine find_gamma_indices(gamma, idx_low, idx_high, alpha)
   ! Linear interpolation weight
   alpha = (gamma - GAMMA_VALUES(idx_low)) / (GAMMA_VALUES(idx_high) - GAMMA_VALUES(idx_low))
```

This allows any γ value, not just multiples of 0.1. The linear interpolation is reasonable since the weights vary smoothly with γ, but it is an extrapolation beyond the paper.

---

## 11. FD Order and Time Integration

### Paper

Implements the method in a **4th-order staggered-grid FD** code (Cui et al., 2010), 2nd-order accurate in time.

### Code

Supports FD orders 2–9 (upwind) or 6 (traditional SBP), with Runge-Kutta time integration (1st–4th order RK, configurable). This is more flexible than the paper.

The memory variable ODE is integrated as part of the same RK stage as the velocity-stress update, which is consistent with the paper's approach (eq. 5).

---

## 12. Summary Table

| Aspect | Paper | Code | Status |
|---|---|---|---|
| Q(f) model | Q₀·(f/f_T)^γ, γ∈[0,0.9] | Same | ✓ |
| N mechanisms | 8 only | 4 or 8 | Extended |
| Table 1 weights | Exact | Exact match | ✓ |
| Table 2 (a_k, b_k) | Exact | Exact match | ✓ |
| Tau formula eq. 15 | Denominator 2N=16 | Same | ✓ |
| High-Q scaling (÷Q) | Yes | Yes | ✓ |
| Low-Q formula (a/Q²+b/Q) | Yes (eq. 16) | Yes | ✓ |
| Gamma interpolation | Discrete (step of 0.1) | Linear between table values | Extension |
| Memory variable ODE | eq. 17 (Day & Bradley 2001) | Matches | ✓ |
| P/S forcing separation | Independent Q_P, Q_S | Q_P = 2×Q_S hard-coded | Simplified |
| Implementation style | Coarse-grained | Conventional (all N/node) | Different |
| wk = N·λk normalization | Explicit (Table 1 caption) | Used directly — **potential 8× error in Qf/Qf8 variants** | ⚠ Bug risk |
| f_trans application | τ_k divided by f_T | f_trans stored but **never applied** | ✗ Bug |
| Smooth transition near f_T | Yes (γ/2 power law) | No (sharp step) | Simplified |
| 4M extraction from 8M | Not described | First 4 only (drops highest-freq mechanism) | ⚠ Degraded |
| FD order | 4th order | 2–9th order configurable | Extended |
| Time integration | 2nd-order (RK2) | 1st–4th order RK configurable | Extended |

---

## 13. Recommended Actions

### Critical (likely incorrect physics)

1. **`f_trans` is never applied** — `anelastic-Qf`, `anelastic-Qf8`, `frequency-Q-4M`, `frequency-Q-8M` all simulate with f_T = 1 Hz regardless of the `f_trans` input parameter. Fix: divide tau arrays by `f_trans` after calling `get_relaxation_times`.

2. **wk normalization in Qf variants** — The Withers Table 1 weights w_k = N·λ_k (N=8) are used directly as λ_k in the conventional approach. If this interpretation is correct, the actual Q is approximately Q_specified/8. Verify against a known analytical or f-k solution with explicit Q.

### Moderate (approximations worth documenting)

3. **4M extraction from 8M table** — Dropping mechanisms 5–8 (especially w₈ = 1.73, the largest) introduces error at high frequencies. Should be documented or replaced with independently fitted 4M weights.

4. **Q_P = 2×Q_S** — Hard-coded ratio prevents independent bulk attenuation. Expose as a configurable parameter.

5. **Sharp vs smooth transition at f_T** — The current sharp step is harder for memory variables to approximate. Could implement the 0.8f_T–1.2f_T smooth ramp from the paper.

### Minor (informational only)

6. **Gamma interpolation** — Linear interpolation between table values is reasonable but undocumented. Could note accuracy bounds (likely < 1% error relative to discrete table).

---

*See also: [RESPONSE_ANALYSIS.md](RESPONSE_ANALYSIS.md) for the full response-type reference, [CODEBASE_ANALYSIS.md](CODEBASE_ANALYSIS.md) for architecture overview.*
