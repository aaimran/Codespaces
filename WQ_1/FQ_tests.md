# FQ_tests.md — Numerical Test Specifications from Withers et al. (2015)

Reference: Withers, K. B., Olsen, K. B., & Day, S. M. (2015).
Memory-Efficient Simulation of Frequency-Dependent Q.
*Bull. Seismol. Soc. Am.*, 105(6), 3129–3142. doi: 10.1785/0120150020

---

## Common Setup (Tests 1–5)

### Grid and Numerics

| Parameter | Value |
|---|---|
| Grid spacing (Δx) | 40 m |
| Time step (Δt) | 0.002 s |
| FD accuracy | 4th-order in space, 2nd-order in time |
| Memory variable approach | Coarse-grained (one mechanism per stress node, cycling with period 2 nodes) |
| Number of mechanisms (N) | 8 |
| Band-pass filter | 0.2–10 Hz, 4th-order zero-phase Butterworth |
| Reference frequency (f₀) | 1.0 Hz (for unrelaxed modulus computation) |
| Domain size | Large enough to exclude boundary reflections at receiver during simulation |

### Point Source

| Parameter | Value |
|---|---|
| Type | Double-couple buried point source |
| Fault geometry | Right-lateral strike-slip: strike 90°, dip 90°, rake 0° |
| Depth | 1.8 km |
| Source time function | Cosine-bell moment rate |
| Duration (T) | 0.2 s |
| Formula | Ṁ(t) = M₀/2 · (1 − cos(2πt/T)) for 0 ≤ t ≤ T, else 0 |
| Seismic moment (M₀) | 10¹⁶ N·m |
| Moment magnitude (Mw) | 4.6 |
| Grid averaging | Source averaged across two depth points; horizontal receiver components averaged to compensate for staggered-grid half-point offset |

### Receiver Station (Tests 1–5)

| Parameter | Value |
|---|---|
| Location | Free surface |
| Horizontal distance from source | 15 km |
| Azimuth from north | 53.13° |
| Components | Radial, Transverse, Vertical |

### Reference Solution

- **Method**: Frequency-wavenumber (f-k) integration (Zhu & Rivera 2002), modified to include power-law Q dispersion
- **Misfit metrics**: Envelope misfit (EM) and phase misfit (PM) from Kristekova et al. (2009)
- **Interpretation**: Misfit values reported for the first 5 s of signal (body wave and early surface wave arrivals)

---

## Test 1 — Elastic Half-Space (Baseline Accuracy)

**Purpose**: Establish FD solver accuracy for prescribed grid parameters before adding attenuation. This is the reference for all subsequent tests.

### Model

| Parameter | Value |
|---|---|
| Geometry | Homogeneous half-space |
| Vp | 6000 m/s |
| Vs | 3464 m/s |
| Density (ρ) | 2700 kg/m³ |
| Q | ∞ (elastic; no memory variables) |

### Expected Output

- Time-series for radial, transverse, and vertical velocity components (cm/s vs. s)
- Fourier amplitude spectra for each component (log scale, 10⁻⁵ to 10⁰ cm·s⁻¹·Hz⁻¹, 10⁰ to 10¹ Hz)
- Rayleigh wave arrival: 4.5–5.5 s on radial and vertical components, free of spurious oscillations

### Acceptance Criteria

| Metric | Threshold |
|---|---|
| Phase misfit (PM) | < 1% for all major phases |
| Amplitude / Envelope misfit (EM) | < 2% peak amplitude relative error |
| Fourier spectral amplitudes | Visually indistinguishable from f-k solution |

### WaveQLab3D Configuration

```fortran
&problem_list
  response = 'elastic'
  physics  = 'elastic'
/
! Vp=6000, Vs=3464, rho=2700 set via material profile
```

---

## Test 2 — Constant-Q Half-Space (γ = 0.0)

**Purpose**: Verify that constant-Q (frequency-independent) memory-variable implementation produces accurate attenuation at low Q (Q = 50).

### Model

| Parameter | Value |
|---|---|
| Geometry | Homogeneous half-space |
| Vp | 6000 m/s |
| Vs | 3464 m/s |
| Density (ρ) | 2700 kg/m³ |
| Qs₀ | 50 |
| Qp₀ | 50 |
| Power-law exponent (γ) | 0.0 (constant Q — frequency independent) |
| Transition frequency (f_T) | 1.0 Hz (irrelevant for γ = 0) |
| Q regime | Low-Q (Table 2, Q = 50 < 200) |

### Q Fitting Parameters

- Weight formula (eq. 16): λ_k = a_k/Q² + b_k/Q
- Coefficients from Table 2, column γ = 0.0
- Relaxation times (eq. 15): τ_k = exp(ln(τ_m) + (2k−1)/16 · ln(τ_M/τ_m)), with τ_m = 0.0032 s, τ_M = 15.9155 s
- Unrelaxed modulus correction applied (eq. 14) using effective Q formula and f_ref = 1 Hz

### Expected Output

- More attenuated waveforms than Test 1; amplitudes reduced (especially at high frequencies)
- Fourier spectra fall steeply above f₀ = 1 Hz

### Acceptance Criteria

| Metric | Threshold |
|---|---|
| Envelope misfit (EM) | < 6% for all components |
| Phase misfit (PM) | < 2% for all components |

Note: Slightly lower accuracy than Test 3 because the constant-Q spectrum is harder to fit at low Q than the smoother power-law case.

### WaveQLab3D Configuration

```fortran
&problem_list
  response = 'frequency-Q-8M'
/

&anelastic_Qf8_list
  c       = 50.0 / 3464.0   ! Qs0 = c * Vs → set c = Q_target / Vs
  gamma   = 0.0
  f_trans = 1.0
  fref    = 1.0
/
```

**Note**: WaveQLab3D hard-codes Qp = 2 × Qs. The paper uses Qp₀ = Qs₀ = 50. This discrepancy cannot be resolved with the current code; Qp will be 100 instead of 50.

---

## Test 3 — Power-Law Q Half-Space (γ = 0.6)

**Purpose**: Verify frequency-dependent Q (Q increases with frequency above f_T). Demonstrates that less energy is attenuated at f > f_T compared to constant Q.

### Model

| Parameter | Value |
|---|---|
| Geometry | Homogeneous half-space |
| Vp | 6000 m/s |
| Vs | 3464 m/s |
| Density (ρ) | 2700 kg/m³ |
| Qs₀ | 50 |
| Qp₀ | 50 |
| Power-law exponent (γ) | 0.6 |
| Transition frequency (f_T) | 1.0 Hz |
| Q law | Q(f) = 50 for f ≤ 1 Hz; Q(f) = 50·(f/1)^0.6 for f > 1 Hz |
| Transition smoothing | Smooth ramp 0.8 f_T to 1.2 f_T |
| Q regime | Low-Q (Table 2, Q = 50) |

### Expected Qualitative Behavior

- More high-frequency energy than Test 2 (constant Q) — less attenuation above f_T
- Fourier spectra show smaller roll-off above 1 Hz compared to Test 2

### Acceptance Criteria

| Metric | Threshold |
|---|---|
| Envelope misfit (EM) | < 2% for all components |
| Phase misfit (PM) | < 1% for all components |

Note: Better accuracy than Test 2 because the power-law Q(f) is smoother and better approximated by the 8-mechanism basis.

### WaveQLab3D Configuration

```fortran
&problem_list
  response = 'frequency-Q-8M'
/

&anelastic_Qf8_list
  c       = 50.0 / 3464.0
  gamma   = 0.6
  f_trans = 1.0
  fref    = 1.0
/
```

---

## Test 4 — Layered Model, Elastic (Baseline for Layered Case)

**Purpose**: Establish FD accuracy for a layered medium that generates reflected and converted phases. Baseline for Test 5.

### Model (Table 3 from paper)

| Layer | Vp (m/s) | Vs (m/s) | ρ (kg/m³) | Thickness (m) |
|---|---|---|---|---|
| 1 (shallow) | 5196 | 3000 | 2550 | 1000 |
| 2 (half-space) | 6000 | 3464 | 2700 | ∞ |

No Q (elastic, Q = ∞ in both layers).

### Expected Output

- Body wave phases including converted P-to-S and S-to-P arrivals from the velocity discontinuity at 1 km depth
- Love waves on transverse component
- Rayleigh waves on radial and vertical components

### Acceptance Criteria

| Metric | Threshold |
|---|---|
| Phase misfit (PM) | < 1% for all major phases (reference from paper: max EM 3%, PM 1%) |
| Amplitude / Envelope misfit (EM) | < 3% |

---

## Test 5 — Layered Model with Frequency-Dependent Q (γ = 0.6)

**Purpose**: Verify that Q(f) works correctly in a layered model with large Q contrast across layers and verify Love wave accuracy.

### Model (Table 3 with Q)

| Layer | Vp (m/s) | Vs (m/s) | ρ (kg/m³) | Qs₀ | Qp₀ | γ | Thickness (m) |
|---|---|---|---|---|---|---|---|
| 1 (shallow) | 5196 | 3000 | 2550 | 20 | 20 | 0.6 | 1000 |
| 2 (half-space) | 6000 | 3464 | 2700 | 210 | 210 | 0.6 | ∞ |

Q contrast: 10× (Q = 20 in shallow layer, 210 in half-space).
Q averaging across coarse-grained cell: harmonic average in shallow layer, arithmetic average in half-space.

### Expected Output

- Waveforms show Q = 20 (low-Q shallow layer) substantially attenuates early arrivals
- Love wave on transverse component (4.5+ s) clearly present but attenuated
- Good agreement between FD and f-k for all phases including post-critical reflections

### Acceptance Criteria

| Phase window | EM | PM |
|---|---|---|
| First 5 s (body waves and early surface waves) | < 5% | < 2% |
| Later surface wave arrivals (5–7 s) | < 10% (8–10% per paper, classified "excellent" by Kristekova scale > 9) | < 2% |

Note: The coarse-grained scheme is known to lose accuracy near sharp Q discontinuities at very low Q. The paper finds the accuracy acceptable for Q ≥ 20 and notes that heterogeneous Q within the coarse-grain cell would reduce the later-arrival misfits.

---

## Test 6 — Chino Hills Earthquake Application (Mw 5.4, 2008)

**Purpose**: Demonstrate the practical importance of frequency-dependent Q for ground-motion prediction. Compare γ = 0.0 vs γ = 0.8 models against recorded strong-motion data.

### Event and Source

| Parameter | Value |
|---|---|
| Event | 2008 Mw 5.4 Chino Hills, California |
| Source model | Finite-fault model adapted from Shao et al. (2012) |
| Source constraints | Constrained by data up to 2.5 Hz; significant energy above that threshold |
| Target frequency | Up to 4 Hz (6.25 points per minimum wavelength) |

### 3D Model

| Parameter | Value |
|---|---|
| Velocity model | SCEC Community Velocity Model v.4 (Magistrale et al. 2000; Kohler et al. 2003) |
| Min Vs imposed | 200 m/s |
| Vp when Vs = 200 m/s | 600 m/s |
| Q₀ relation | Qs₀ = Vs (km/s) × 10; Qp₀ = 2 × Qs₀ |
| Min Q₀ | > 20 (imposed by Vs lower limit) |
| Boundary conditions | Cerjan sponge zones (no free-surface reflection at edges) |
| Grid spacing | 8 m |
| Domain | 56 km (E–W) × 40 km (N–S) × 24 km depth |
| Memory variables | Coarse-grained, period-2 cycling, N = 8 mechanisms |

### Q Models Compared

| Label | γ | Description |
|---|---|---|
| Constant-Q | 0.0 | Frequency-independent Q throughout |
| Power-law Q | 0.8 | Frequency-dependent Q (upper bound of Song & Jordan 2013 estimate for southern California) |

### Observation Dataset

- 110 strong-motion stations from the Center for Engineering Strong Motion Data
- East–west component used for comparison
- Seismograms bandpassed 0.1–4 Hz
- Energy plots bandpassed 1.0–4.0 Hz (to emphasize frequency-dependent Q regime)

### Comparison Metrics

| Metric | Description |
|---|---|
| Cumulative energy | Geometric mean of energy (integral of v²) over horizontal components per station vs. Rrup |
| Fourier amplitude | Narrow-band (±0.05 Hz) average spectral amplitude at 0.25, 2.5, 3.5 Hz vs. Rrup |
| Spectral acceleration | GMRotD50 at periods 3.0 s, 1.0 s, 0.4 s vs. Rrup; compared to GMPE range (Boore & Atkinson 2008, Campbell & Bozorgnia 2008, Chiou & Youngs 2008) |

### Expected Outcomes

| Observation | Expected result |
|---|---|
| Waveform differences | Diverge above 1 Hz; γ = 0.8 has more energy in coda |
| Cumulative energy (1–4 Hz) vs Rrup | Constant-Q too attenuating at Rrup > 25 km; γ = 0.8 better matches data |
| SA at 3 s (< 1 Hz regime) | Both models nearly identical (Q(f) negligible at low f) |
| SA at 0.4 s (> 1 Hz regime) | γ = 0.8 clearly better than γ = 0.0 at larger distances |
| Fourier amplitude at 3.5 Hz | Constant-Q deficient by factor 3–5 relative to data and γ = 0.8 at Rrup > 25 km |

---

## Q Fitting Accuracy (Paper Figures 1 and 2)

### High-Q Regime (Q > 200, low-loss approximation, Figure 1)

Weights λ_k fit by least squares to Q_approx(f) / Q_target ≈ 1.

| γ | Max fit error over [0.1–10 Hz] |
|---|---|
| 0.0 | < 5% |
| 0.2 | < 5% |
| 0.5 | < 5% |
| 0.8 | < 5% (ratio Q_approx/Q_target oscillates about 1.0 across Debye peaks) |

Tabulated as w_k = N·λ_k in Table 1 (N = 8). **For the conventional approach (all N mechanisms at every node), the code must use λ_k = w_k / N.**

### Low-Q Regime (15 ≤ Q ≤ 200, Table 2 interpolation, Figure 2)

Effective weights λ_k computed from quadratic formula (eq. 16):

```
λ_k = a_k / Q² + b_k / Q
```

Coefficients a_k, b_k tabulated in Table 2 for k = 1…8 and γ ∈ {0.0, 0.1, …, 0.9}.

| γ | Max fit error over [0.1–10 Hz] at Q = 20 |
|---|---|
| 0.0 | < 5% |
| 0.2 | < 5% |
| 0.5 | < 5% |
| 0.8 | Largest error (sharp Q transition + low Q is hardest case) |

---

## Relaxation Time Formula (Paper Equation 15)

```
τ_k = exp( ln(τ_m) + (2k−1)/16 · ln(τ_M/τ_m) )    k = 1…8
```

Denominator 16 = 2 × N with N = 8.

Band limits per γ (Table 1):

| γ | τ_m (s) | τ_M (s) |
|---|---|---|
| 0.0 | 0.0032 | 15.9155 |
| 0.1 | 0.0032 | 15.9155 |
| 0.2 | 0.0020 | 15.9155 |
| 0.3 | 0.0016 | 15.9155 |
| 0.4 | 0.0013 | 15.9155 |
| 0.5 | 0.0010 | 15.9155 |
| 0.6 | 0.0008 | 15.9155 |
| 0.7 | 0.0005 | 15.9155 |
| 0.8 | 0.0004 | 15.9155 |
| 0.9 | 0.0002 | 15.9155 |

Transition frequency shift: divide τ_k by f_T to center absorption band at the desired frequency.

---

## Code Implementation Notes for WaveQLab3D

### Bug Fixes Required (applied in material.f90)

Two bugs affecting `frequency-Q-4M` / `frequency-Q-8M` / `anelastic-Qf` have been fixed:

**Fix 1 — f_trans not applied (silent bug)**

The paper states τ_k must be divided by f_T to shift the band. Before the fix, f_T was read but never used.

```fortran
! In init_anelastic_Qf8_properties (and init_anelastic_Qf_properties):
call get_relaxation_times(M%gamma_Qf8, M%tau_Qf8)
M%tau_Qf8 = M%tau_Qf8 / M%f_trans_Qf8   ! <-- fix: shift to f_T
```

**Fix 2 — w_k = N·λ_k normalization error (~8× Q error)**

Table 1 stores w_k = N·λ_k (N = 8). The conventional approach requires λ_k.

```fortran
call get_withers_weights(M%gamma_Qf8, 1.0_wp, M%weight_Qf8)
M%weight_Qf8 = M%weight_Qf8 / real(N_MECH, wp)   ! <-- fix: w_k → λ_k
```

### Known Code–Paper Discrepancies

| Issue | Paper | WaveQLab3D | Impact |
|---|---|---|---|
| Q_P independence | Qp₀ and Qs₀ set independently | Qp = 2 × Qs hard-coded | Cannot reproduce Qp₀ = Qs₀ = 50 from Tests 2–5 |
| f_T transition shape | Smooth ramp 0.8 f_T – 1.2 f_T | Sharp step at f_T | Minor spectral artifact near f_T |
| 4M variant | Uses first 4 of 8 Withers mechanisms | Drops mechanism 8 (largest weight, w_k = 1.73) | Reduced high-frequency fit accuracy |
| Memory variable approach | Coarse-grained (1 mechanism per stress node) | Conventional (all N at every node) | Different memory cost; both physically valid |

### Suggested Input File (Test 3 — power-law γ = 0.6)

```fortran
&problem_list
  response = 'frequency-Q-8M'
  physics  = 'elastic'
/

&anelastic_Qf8_list
  c       = 50.0 / 3464.0   ! Qs0 = c * Vs → Qs0 = 50 for Vs = 3464 m/s
  gamma   = 0.6
  f_trans = 1.0             ! transition frequency (Hz) — now applied after fix
  fref    = 1.0             ! reference frequency for modulus correction
/
```

For Test 2 (constant Q), change `gamma = 0.0`. For Test 5 (layered model), specify two material regions with the parameters from Table 3 above.

---

*Derived from Withers, Olsen & Day (2015), BSSA 105(6). See also FQ_PAPER_VS_CODE.md for the full paper-vs-code comparison and RESPONSE_ANALYSIS.md for a description of all response types.*
