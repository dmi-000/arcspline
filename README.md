# arcspline

C^N-continuous parametric curve interpolation using local geometric-arc windows.
Header-only C++17, no dependencies, works in any dimension ≥ 2.

---

## What it does

Given a sequence of control points and their parameter values, `arcspline` builds a
smooth curve that passes **exactly** through every control point.  Smoothness is C^N —
position and the first N derivatives are continuous everywhere, not just at knots.

Each local segment is built from a best-fit conic arc (ellipse, parabola, or hyperbola)
through 5 consecutive control points.  Adjacent windows overlap and are blended with a
smoothstep weight, giving exact C^N continuity without solving any global system.

**Key properties:**
- **Exact interpolation** — curve passes through all interior control points
- **C^N continuity** — exact, algebraic, no tolerance (default N=2, quintic blend)
- **Affine invariant** — stretching or shearing the control points gives the same curve,
  stretched (~2e-13 error); similarity-invariant for circle windows
- **Exact conic reproduction** — if input lies on a conic, output is machine-epsilon exact
- **nD support** — works in ℝ² through ℝⁿ via `template<int Dim>`
- **Header-only C++17** — single `#include`, no build step, no dependencies
- **Local support** — changing one control point affects at most 4 consecutive segments

---

## Quick start

```cpp
#include "arcspline_nd.hpp"   // all windows; includes arcspline.hpp + cylinder
#include <vector>
#include <cmath>

int main()
{
    // 3D helix: P(t) = (cos t, sin t, 0.3 t)
    int n = 20;
    std::vector<fc::Vec3>   ctrl(n);
    std::vector<double>     times(n);
    for (int i = 0; i < n; ++i) {
        double t = 4.0 * M_PI * i / (n - 1);
        ctrl[i]  = fc::Vec3(std::cos(t), std::sin(t), 0.3 * t);
        times[i] = t;
    }

    // 3D: automatic best window (Clifford/cylinder/conic/FH fallback)
    auto r3 = fc::blend_curve<3>(ctrl, times, fc::nd_tag{}, 80, 2);

    // 3D: explicit cylinder windows (best for helices and space curves)
    auto r_cyl = fc::blend_curve(ctrl, times, fc::cylinder_tag{}, 80, 2);

    // N-dimensional: works for any Dim >= 2
    using V4 = fc::VecN<4>;
    std::vector<V4>     ctrl4(n);
    std::vector<double> times4(n);
    // ... fill ctrl4 ...
    auto r4 = fc::blend_curve<4>(ctrl4, times4, fc::nd_tag{}, 80, 2);

    // r.pts   — std::vector<VecN<Dim>>, dense output points
    // r.times — corresponding parameter values
}
```

Compile:
```bash
clang++ -std=c++17 -O2 -I. -o my_program my_program.cpp
```

---

## Four headers

| Header | Window type | Points | Min n | Exact for | Invariance | Dim |
|---|---|---|---|---|---|---|
| `arcspline_circle.hpp` | Circle arc | 3 | 4 | Circles | Similarity | any ≥ 2 |
| `arcspline.hpp` | Conic arc | 5 | 6 | All conics | **Affine** | any ≥ 2 |
| `arcspline_cylinder.hpp` | Cylinder geodesic | 5 | 6 | Helices, geodesics | Euclidean | 3 only |
| `arcspline_nd.hpp` | Clifford torus S¹×S¹ | 5 | 6 | Torus helices | Euclidean | any ≥ 4 |

Each header includes the ones above it. `arcspline_nd.hpp` is the recommended single include.

### Tag dispatch

```cpp
#include "arcspline_nd.hpp"

auto r = fc::blend_curve(ctrl, times);                     // circle (3D, no tag)
auto r = fc::blend_curve(ctrl, times, fc::conic_tag{});    // conic (affine-invariant)
auto r = fc::blend_curve(ctrl, times, fc::cylinder_tag{}); // cylinder (3D, best for helices)
auto r = fc::blend_curve<N>(ctrl, times, fc::nd_tag{});    // automatic (N=2→conic, 3→cylinder, ≥4→Clifford)
```

### Fallback template parameter

Every `blend_curve` overload accepts an optional `Fallback` template parameter
(default: `LagrangeWindow`) for windows where the geometric fit fails:

```cpp
// Use Floater-Hormann Depth=3 rational fallback instead of degree-4 polynomial:
auto r = fc::blend_curve<4, fc::FHWindow3>(ctrl, times, fc::nd_tag{}, 80, 2);

// Explicit Depth:
template<int Dim> using FH1 = fc::FHWindow<Dim, 1>;
auto r = fc::blend_curve<3, FH1>(ctrl, times, fc::cylinder_tag{}, 80, 2);
```

`FHWindow<Dim, Depth>` is a pole-free barycentric rational interpolant (Floater-Hormann
2007). Unlike per-coordinate rationals, its weights depend only on the node times —
not coordinate values — making it rotationally invariant. `Depth` is the degree of the
local polynomial sub-interpolants being blended (not the degree of the resulting rational).
`Depth=4` reduces exactly to `LagrangeWindow`; `Depth=3` (default via `FHWindow3`)
matches Lagrange's O(h⁵) accuracy.

---

## Build and test

```bash
CXX="g++ -std=c++17 -O2 -I."

# Primary regression suites
$CXX -o test_clifford_nd   test_clifford_nd.cpp   && ./test_clifford_nd    # 21/21
$CXX -o test_cylinder_edge test_cylinder_edge.cpp && ./test_cylinder_edge  # 29/29

# Clifford torus stress test (30 parametric cases)
$CXX -o diag_clifford diag_clifford.cpp && ./diag_clifford  # 27/30

# Integration demos
$CXX -o demo       demo.cpp       && ./demo
$CXX -o demo_conic demo_conic.cpp && ./demo_conic
$CXX -o demo_cylinder demo_cylinder.cpp && ./demo_cylinder

# Or via CMake
cmake -B build && cmake --build build && ctest --test-dir build
```

Test coverage: exact interpolation, C^N continuity, similarity/affine invariance,
collinearity guards, all fallback levels, FHWindow identity (Depth=4=Lagrange),
nD (Dim=2…6), Clifford torus helices, cylinder helices and twisted cubics.

---

## Documentation

Full documentation — architecture, mathematical properties, invariance proofs, API
reference, comparison with NURBS — is in [`arcspline_doc.md`](arcspline_doc.md).

---

## License

MIT — see [LICENSE](LICENSE).
