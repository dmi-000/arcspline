# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build commands

All test and diagnostic executables are single-file C++17 with no external dependencies:

```bash
# Standard compile prefix (macOS, use clang++ or g++):
CXXFLAGS="-std=c++17 -O2 -I."

# Primary regression tests
g++ $CXXFLAGS -o test_clifford_nd   test_clifford_nd.cpp   && ./test_clifford_nd
g++ $CXXFLAGS -o test_cylinder_edge test_cylinder_edge.cpp && ./test_cylinder_edge

# Stress-test diagnostic for CliffordWindow (30 parametric cases, reports pass/fail)
g++ $CXXFLAGS -o diag_clifford diag_clifford.cpp && ./diag_clifford

# Helix/cylinder diagnostics
g++ $CXXFLAGS -o diag_cylinder  diag_cylinder.cpp  && ./diag_cylinder
g++ $CXXFLAGS -o diag_cylinder2 diag_cylinder2.cpp && ./diag_cylinder2

# Integration demos
g++ $CXXFLAGS -o demo      demo.cpp      && ./demo
g++ $CXXFLAGS -o demo_nd   demo_nd.cpp   && ./demo_nd
g++ $CXXFLAGS -o demo_conic demo_conic.cpp && ./demo_conic
```

**Regression baseline**: `diag_clifford` (target: 26–27/30), `test_clifford_nd` (21/21), `test_cylinder_edge` (29/29).

## Architecture: four header layers

The library is header-only C++20, all in namespace `fc`, all `template<int Dim>`:
```

The solver fits 5 points in ℝᴺ to a Clifford torus `S¹(r₁) × S¹(r₂)`:

1. **Affine hull** (`Hull4D`): project 5 points into their 4D span via Gram-Schmidt. If rank < 4 → fall to lower-dim window.
2. **Sphere fit** (`sphere_fit_4d`): fit circumscribed 4-sphere to hull coords.
3. **Moment matrix**: compute `M = Σ (p-c)(p-c)ᵀ`, eigendecompose into two planes.
4. **ESPRIT propagator**: `A₀ = Q₁ Q₀⁻¹`; symmetrized `S = (A+Aᵀ)/2` → eigenvalues give `cos(ω₁·dt)`, `cos(ω₂·dt)`.
5. **Frame/radii extraction**: angles `φᵢ` from eigendecomposition → LS fit `uv(t) = a·t + b` → extract frames `f1,f2,f3,f4` and `r1,r2`.
6. **Grid search refinement** (`extract_and_probe`): for borderline cases, grid over `(ω₁,ω₂)` using GJ direct formula to minimize residual.
7. **Accept** if `best_err < 1e-6 · spread`.

**Direct formula path** (`use_direct_`): stores `o1d_,o2d_,Av_,Bv_,Av2_,Bv2_` and evaluates `p(t) = center + Av·cos(ω₁t) + Bv·sin(ω₁t) + Av2·cos(ω₂t) + Bv2·sin(ω₂t)` bypassing angle extraction. Used when GJ formula gives a lower residual than the angle-parameterization path.

## Known diag_clifford failure cases (correct behavior)

These 4 cases fail by design — the geometry is genuinely lower-dimensional or aliased:
- **om2=2π−1** (hull_rank=2, circle): algebraically correct clf=NO
- **om1=π** (hull_rank=3, sin(πk)=0 for integer k blocks identification)
- **om1=2, om2=π** (hull_rank=3): same aliasing
- **om1·dt≈2π** (near-repeat knot): clf=NO is correct

## Development workflow

The key metric is `diag_clifford` (currently 26/30 — cases 1–3 and 5 above are unfixable; case 4 `dt=0.01` is the active work item (passes as of 2026-07-07)).

The `.bak` files (e.g. `arcspline.hpp.bak`) represent the last confirmed-good state. Do not overwrite them.
