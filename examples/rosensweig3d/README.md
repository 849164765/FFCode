# rosensweig3d — 3D Rosensweig instability in a ferrofluid layer

3D extension of **Section III.D** of *Z. Guo et al., "Phase-field lattice Boltzmann
model with adaptive mesh refinement for ferrofluid interfacial dynamics,"
Phys. Fluids 37, 022148 (2025)*.

## Physics setup

A rectangular cavity (X, Y periodic, Z walls) holds an immiscible two-phase system:
the **bottom 1/3 is ferrofluid** (order parameter φ=1, heavy) and the **top 2/3 is
organic solvent** (φ=0, light). A uniform **vertical magnetic field H0** is applied
along +Z (Neumann wall BC ∂ψ/∂z = −H0). Gravity acts downward and stabilises the
interface; the magnetic normal-field stress destabilises it. Above the critical
field Hc the flat interface breaks into a **hexagonal array of peaks** (the
Rosensweig / normal-field instability).

The lattice parameters follow the paper (lattice units):
ρ_l=1.0, ρ_h=1.975, η_l=0.0073, η_h=0.292, σ=0.0071, W=4, M=0.01,
μ_l=1.0 (organic), μ_h=2.2 (ferrofluid relative permeability), χ_h=μ_h−1=1.2.

### Derived quantities

Gravity is set from the critical wavelength `LambdaC` (Cowley–Rosensweig,
paper Eq. 71):
```
λc = 2π√(σ/(gΔρ))   →   g = σ/Δρ·(2π/λc)²
Hc = √(2(μ0/μ+1)/(μ0/μ−1)²)·(σgΔρ)^{1/4}   (μ0=1)
```
For the shipped `.ini` (`LambdaC=25`): g ≈ 4.60e−4, **Hc ≈ 0.132**.
H0/Hc controls the peak growth:
- H0=0.15 (1.14·Hc): weak ripples
- H0=0.20 (1.51·Hc): clear hexagonal peaks (production default)
- H0=0.30 (2.3·Hc): strong peaks, grows fast

### Initial condition

The interface is a slightly perturbed plane at `Interface_Z` (bottom 1/3 of the
height). `Perturb_Mode` seeds the instability:
- `hex` (default): three plane waves at 60° → hexagonal peak array at λc
- `xwave`: parallel ridges (2D-like)
- `noise`: deterministic low-mode mix

The NS field is initialised with the **hydrostatic pressure** profile so the
interface starts near equilibrium and the instability develops from a quiet base
(no settling transient).

## Numerics

Same framework as `bubbleMag3d`: Allen–Cahn phase field (D3Q19 MRT),
velocity-based Navier–Stokes (D3Q19 MRT, Guo force), and magnetic scalar
potential ψ (D3Q7 MRT diffusion, sub-iterated). ψ is pinned to −H0·z at the Z
walls and the MF periodic ghost planes in X/Y are synchronised across MPI ranks
(see the `SyncMFPeriodicGhosts` lambda, copied from bubbleMag3d).

Per output step the run prints `z_if max/min/amp` (global φ=0.5 crossing height),
`|u|max`, `|H|max/H0`, and `<φ>` (ferrofluid volume fraction).

## Run

On the cluster (targets 96 ranks, one 32³ block per rank for the shipped mesh):

```bash
sbatch dlf.sh
```

Locally (any number of ranks; `BlockCellLen` sets the block size):

```bash
make            # needs mpic++
mpirun -n N ./rosensweig3d.exe
```

VTK output goes to `vtkoutput/vtidata/` (PHI, HMAG, PSI, Velocity).

## Analyse

```bash
python3 analyze_rosensweig.py vtkoutput/vtidata --mesh 128,128,192 --iface 64
```

Prints per-step interface statistics and, with `--plot`, saves a 3D surface plot
of the peak array plus a 2D peak map. The peak-to-peak spacing can be compared to
the critical wavelength `LambdaC`.

## Tuning notes

- **Run length**: peak amplitude grows roughly as exp(0.0023·t) at H0=0.17
  (measured on a 48×48×96 proxy). On the production mesh the hexagonal array
  organises over ~10⁴ steps. Raise H0 to grow peaks faster.
- **Block layout** is computed from `BlockCellLen` (each block ≈ BlockCellLen³),
  so the same binary works for any mesh and rank count.
- The magnetic force is zeroed within `MFWallBand` rows of the Z walls (solver
  artifact region); the interface is far from the walls so this is safe.
