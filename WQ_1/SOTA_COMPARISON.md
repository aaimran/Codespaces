# SOTA of Dynamic Rupture & Seismic Wave Propagation Simulation — and Where WaveQLab3D Stands

*Survey date: 2026-06-11. Based on a literature scan plus direct inspection of the WaveQLab3D source in this repository.*

## Scope

State of the art in **dynamic earthquake rupture** and **seismic wave propagation** simulation, with emphasis on two specific capabilities:

1. **Anelastic attenuation** (viscoelasticity, typically via Generalized Standard Linear Solid / memory variables)
2. **PML** (perfectly matched layer) absorbing boundaries

---

## The SOTA landscape (mid-2026)

The field has split into two overlapping fronts: **dynamic rupture** (spontaneous, friction-governed sources) and **large-scale ground-motion propagation** (usually kinematic sources).

### SeisSol (LMU/TUM) — the overall dynamic-rupture leader

- ADER-DG on unstructured tetrahedra; arbitrary high order in space and time; clustered local time stepping.
- Code-generation backend (YATeTo / ChainForge fused-GEMM) giving portable CUDA/HIP performance on LUMI- and Frontier-class machines (>33 SP-TFLOP/s demonstrated on 5M-element meshes, 2024).
- Attenuation: GSLS viscoelasticity integrated into the high-order scheme; also anisotropy, poroelasticity, off-fault Drucker–Prager plasticity, thermal pressurization.
- Defining capability is **geometric**: intersecting fault networks and topography on unstructured meshes — e.g., the rapid 3D dynamic-rupture models of the 2023 Kahramanmaraş Mw 7.8/7.7 doublet, and multiscale fracture-network cascade studies.
- Pushing into ensemble/Bayesian dynamic-rupture inference via fused simulations.
- Notably does **not** use PML — absorbing fluxes on extended domains, since stable PML for high-order DG remains awkward.

### EQSIM / SW4 (LBNL/LLNL) — the exascale ground-motion leader

- 4th-order SBP finite differences, curvilinear mesh refinement, topography; GSLS attenuation.
- GPU-resident via RAJA; DOE ECP EQSIM workflow runs regional fault-to-structure simulations at 5–10 Hz on Frontier/Aurora-class machines.
- **Kinematic sources only** — no spontaneous rupture.
- Deliberately avoids PML in favor of Petersson & Sjögreen's "supergrid" layers, because of PML's long-time stability problems.

### AWP-ODC (SCEC/SDSU) — the throughput workhorse

- Staggered-grid FD, GPU-optimized; used for CyberShake and the largest ground-motion ensembles.
- Frequency-dependent Q via coarse-grained memory variables (Withers/Olsen lineage — the same family as WaveQLab3D's Qf models).
- Has a dynamic-rupture mode (SGSN) but is mostly used kinematically.

### SPECFEM3D

- Spectral elements, GPU (CUDA + HIP), viscoelastic attenuation.
- One of the few production codes with **C-PML actually deployed**.
- Includes a dynamic-rupture solver (Galvez/Kame lineage), though that is a smaller part of its use.

### The SBP/SEAS research tier

- **Tandem** (DG, earthquake sequences), **Thrase** and Erickson-group SBP codes, **FD3D_TSN** (lean GPU rupture code).
- Dual-pairing and upwind SBP methods (Duru's newer ANU work) — the direct intellectual successors of WaveQLab3D.
- These lead on *numerical rigor* (provable stability of friction coupling) rather than scale.

### PML state of the art specifically

- PML remains the most accurate absorbing boundary per grid point; the unsolved problem is **long-time stability** in elastic media (grazing incidence, surface waves, anisotropy).
- SOTA formulations: complex-frequency-shifted C-PML, M-PML, and **energy-stable PML within SBP-SAT frameworks** — largely Duru's research program.
- Tellingly, most production rupture codes (SeisSol, SW4, AWP) *avoid* PML. The combination of dynamic rupture + viscoelasticity + PML in one code is rare.

---

## Where WaveQLab3D stands

Based on the in-depth analyses in `WQ/CODEBASE_ANALYSIS.md` and `WQ/RESPONSE_ANALYSIS.md`, plus direct source verification.

### Genuine differentiator — exactly the asked-about combination

- PML module (`src/pml.f90`) wired into the RHS on **all six block faces**.
- Dedicated `apply_anelastic_point_pml` variants in `src/RHS_Interior.f90` so GSLS memory-variable updates are computed *inside* the PML region.
- I.e., **viscoelastic–PML coupling for dynamic rupture in a provably stable SBP-SAT discretization on curvilinear grids with nonplanar rate-and-state faults**. Almost no other code offers that specific stack; the leaders sidestep PML entirely.
- Attenuation menu (4/8-mechanism constant-Q and frequency-dependent Q(f), Withers-style) is competitive with AWP-ODC's and more flexible than SW4's.

### Where it trails the leaders

| Dimension | Gap |
|---|---|
| Geometry | Block-structured curvilinear grids, single fault interface. No intersecting fault networks, branching, or topography-conforming unstructured meshes (SeisSol's domain). |
| Hardware | Pure MPI Fortran, no GPU port, static load balancing, no local time stepping. Petascale CPU scaling demonstrated ~2015–2017; leaders are now exascale-GPU. Roughly a 1–2 order-of-magnitude effective-throughput gap on current machines. |
| Physics breadth | No thermal pressurization, poroelasticity, or multi-physics coupling (tsunami, seismo-acoustic). |
| Engineering health | No CI; duplicated BC modules; 23K-line generated stencil file; Q-dispatch bugs found and fixed locally 2026-06-11 (working tree only — the WQ/ gitlink is stale and git tracks nothing inside it). |

### Bottom line

- **SeisSol** is the current SOTA for dynamic rupture broadly; **EQSIM/SW4** and **AWP-ODC** for exascale ground motion.
- WaveQLab3D is no longer competitive on scale or geometric generality, but it **remains state of the art in a narrow, mathematically meaningful niche**: high-order provably stable SBP-SAT rupture dynamics with stable PML and viscoelastic attenuation on curvilinear nonplanar faults.
- Best viewed today as a rigorous research code whose ideas live on in Duru's upwind/dual-pairing SBP successors.

---

## Sources

- [SeisSol documentation](https://seissol.readthedocs.io/) · [SeisSol citation/feature list](https://seissol.org/about/howtocite/)
- [Fused ensembles of dynamic-rupture simulations for Bayesian inference](https://pmc.ncbi.nlm.nih.gov/articles/PMC13065561/)
- [Rapid 3D dynamic rupture modeling of the 2023 Kahramanmaraş doublet (arXiv)](https://arxiv.org/pdf/2308.02144)
- [Fused GEMMs for GPU ADER-DG in SeisSol (2024)](https://onlinelibrary.wiley.com/doi/full/10.1002/cpe.8037)
- [Petascale high-order dynamic rupture on heterogeneous supercomputers (IEEE)](https://ieeexplore.ieee.org/document/7012188/)
- [PML review and recent developments (arXiv)](https://arxiv.org/pdf/2104.09854)
- [C-PML for the viscoelastic wave equation (GJI)](https://academic.oup.com/gji/article/179/1/333/738594)
- [OpenSWPC: viscoelastic FD with PML](https://link.springer.com/article/10.1186/s40623-017-0687-2)
- [SCEC/USGS dynamic rupture code verification exercise](https://pubs.usgs.gov/publication/70035802) · [benchmark suite](https://authors.library.caltech.edu/records/mz6nk-pxh25)
- [WaveQLab3D: nonplanar faults in complex elastic solids (Duru & Dunham, JCP 2016)](https://www.sciencedirect.com/science/article/abs/pii/S0021999115006853)
- [Upwind SBP methods for 3D elastic waves (arXiv)](https://arxiv.org/pdf/2011.02600) · [Dual-pairing SBP rupture methods (arXiv)](https://arxiv.org/pdf/2207.07891)
- [FD3D_TSN: GPU dynamic rupture code](https://geo.mff.cuni.cz/~gallovic/abst/Premus.etal.SRL.2020.pdf)
