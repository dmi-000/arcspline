// test_clifford_nd.cpp — 5 tests for conicblend_nd.hpp
//
// Build:
//   g++ -std=c++17 -O2 -I. \
//       -isystem /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1 \
//       -o test_clifford_nd test_clifford_nd.cpp && ./test_clifford_nd

#include "conicblend_nd.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

// ── Helper ─────────────────────────────────────────────────────────────────

static int n_pass = 0, n_fail = 0;

static void check(bool ok, const char* label, const char* detail = "")
{
    if (ok) {
        std::printf("  PASS  %s  %s\n", label, detail);
        ++n_pass;
    } else {
        std::printf("  FAIL  %s  %s\n", label, detail);
        ++n_fail;
    }
}

// ── Test 1: Helix on Clifford torus in ℝ⁴ ─────────────────────────────────
//
// Curve: p(t) = (r1·cos(ω1·t), r1·sin(ω1·t), r2·cos(ω2·t), r2·sin(ω2·t))
// with r1=2, r2=1, ω1=1, ω2=√2.
// 5 consecutive knots → CliffordWindow<4> → knot error < 1e-12.

static void test1()
{
    std::printf("\n-- T1: Helix on Clifford torus in R^4 --\n");
    using V4 = fc::VecN<4>;

    const double r1 = 2.0, r2 = 1.0;
    const double om1 = 1.0, om2 = std::sqrt(2.0);

    auto torus_pt = [&](double t) -> V4 {
        return V4{r1*std::cos(om1*t), r1*std::sin(om1*t),
                  r2*std::cos(om2*t), r2*std::sin(om2*t)};
    };

    const double dt = 1.0;
    double ts[5];
    V4     pts[5];
    for (int k = 0; k < 5; ++k) {
        ts[k]  = k * dt;
        pts[k] = torus_pt(ts[k]);
    }

    fc::CliffordWindow<4> win(pts[0],pts[1],pts[2],pts[3],pts[4],
                               ts[0], ts[1], ts[2], ts[3], ts[4]);

    check(win.valid(),         "T1a: valid()");
    check(win.used_clifford(), "T1b: used_clifford()");
    check(win.hull_rank() == 4, "T1c: hull_rank==4");

    // Knot error
    double max_err = 0.0;
    for (int k = 0; k < 5; ++k)
        max_err = std::max(max_err, (win(ts[k]) - pts[k]).norm());

    char buf[64]; std::snprintf(buf, sizeof(buf), "(max_err=%.2e)", max_err);
    check(max_err < 1e-12, "T1d: knot error < 1e-12", buf);
}


// ── Test 2: Helix in 3D sub-hull of ℝ⁵ ────────────────────────────────────
//
// Curve: p(t) = (r·cos(t), r·sin(t), h·t, 0, 0) embedded in ℝ⁵.
// Hull rank = 3 (two circular dims + axial dim, no 4th independent direction).
// CliffordWindow<5> should use level 2 (CylinderWindow<3>).

static void test2()
{
    std::printf("\n-- T2: Helix in 3D sub-hull of R^5 --\n");
    using V5 = fc::VecN<5>;

    const double r = 1.5, h = 0.4;

    auto helix_pt = [&](double t) -> V5 {
        return V5{r*std::cos(t), r*std::sin(t), h*t, 0.0, 0.0};
    };

    const double dt = 1.0;
    double ts[5];
    V5     pts[5];
    for (int k = 0; k < 5; ++k) {
        ts[k]  = k * dt;
        pts[k] = helix_pt(ts[k]);
    }

    fc::CliffordWindow<5> win(pts[0],pts[1],pts[2],pts[3],pts[4],
                               ts[0], ts[1], ts[2], ts[3], ts[4]);

    check(win.valid(),          "T2a: valid()");
    check(win.used_cylinder(),  "T2b: used_cylinder()");
    check(win.hull_rank() == 3, "T2c: hull_rank==3");

    double max_err = 0.0;
    for (int k = 0; k < 5; ++k)
        max_err = std::max(max_err, (win(ts[k]) - pts[k]).norm());

    char buf[64]; std::snprintf(buf, sizeof(buf), "(max_err=%.2e)", max_err);
    check(max_err < 1e-12, "T2d: knot error < 1e-12", buf);
}


// ── Test 3: Parabola in 2D sub-hull of ℝ⁵ ─────────────────────────────────
//
// Curve: p(t) = (t, t², 0, 0, 0) in ℝ⁵.  Hull rank = 2.
// CliffordWindow<5> should use level 3 (ConicWindow<2>).

static void test3()
{
    std::printf("\n-- T3: Parabola in 2D sub-hull of R^5 --\n");
    using V5 = fc::VecN<5>;

    auto parab_pt = [](double t) -> V5 {
        return V5{t, t*t, 0.0, 0.0, 0.0};
    };

    // Use one-sided arm to ensure monotone φ in ConicWindow<2>
    double ts[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
    V5 pts[5];
    for (int k = 0; k < 5; ++k) pts[k] = parab_pt(ts[k]);

    fc::CliffordWindow<5> win(pts[0],pts[1],pts[2],pts[3],pts[4],
                               ts[0], ts[1], ts[2], ts[3], ts[4]);

    check(win.valid(),          "T3a: valid()");
    check(win.used_planar(),    "T3b: used_planar()");
    check(win.hull_rank() == 2, "T3c: hull_rank==2");

    double max_err = 0.0;
    for (int k = 0; k < 5; ++k)
        max_err = std::max(max_err, (win(ts[k]) - pts[k]).norm());

    char buf[64]; std::snprintf(buf, sizeof(buf), "(max_err=%.2e)", max_err);
    check(max_err < 1e-12, "T3d: knot error < 1e-12", buf);
}


// ── Test 4: 4D torus helix embedded in ℝ⁶ ────────────────────────────────
//
// Same torus helix as T1, embedded in ℝ⁶ as (p1,p2,p3,p4,0,0).
// affine_hull_4d should reduce to rank 4; Clifford torus should succeed.

static void test4()
{
    std::printf("\n-- T4: 4D torus helix embedded in R^6 --\n");
    using V6 = fc::VecN<6>;

    const double r1 = 2.0, r2 = 1.0;
    const double om1 = 1.0, om2 = std::sqrt(2.0);

    auto torus6_pt = [&](double t) -> V6 {
        return V6{r1*std::cos(om1*t), r1*std::sin(om1*t),
                  r2*std::cos(om2*t), r2*std::sin(om2*t), 0.0, 0.0};
    };

    const double dt = 1.0;
    double ts[5];
    V6     pts[5];
    for (int k = 0; k < 5; ++k) {
        ts[k]  = k * dt;
        pts[k] = torus6_pt(ts[k]);
    }

    fc::CliffordWindow<6> win(pts[0],pts[1],pts[2],pts[3],pts[4],
                               ts[0], ts[1], ts[2], ts[3], ts[4]);

    check(win.valid(),          "T4a: valid()");
    check(win.used_clifford(),  "T4b: used_clifford()");
    check(win.hull_rank() == 4, "T4c: hull_rank==4");

    double max_err = 0.0;
    for (int k = 0; k < 5; ++k)
        max_err = std::max(max_err, (win(ts[k]) - pts[k]).norm());

    char buf[64]; std::snprintf(buf, sizeof(buf), "(max_err=%.2e)", max_err);
    check(max_err < 1e-12, "T4d: knot error < 1e-12", buf);
}


// ── Test 5: Full blend regression, 4D helix n=8 ───────────────────────────
//
// 8 uniformly-spaced control points on the torus helix in ℝ⁴.
// blend_curve<4> with nd_tag; verify:
//   (a) used_clifford = true
//   (b) max_dev from true torus curve < 1e-13

static void test5()
{
    std::printf("\n-- T5: Full blend regression, 4D helix n=8 --\n");
    using V4 = fc::VecN<4>;

    const double r1 = 2.0, r2 = 1.0;
    const double om1 = 1.0, om2 = std::sqrt(2.0);

    auto torus_pt = [&](double t) -> V4 {
        return V4{r1*std::cos(om1*t), r1*std::sin(om1*t),
                  r2*std::cos(om2*t), r2*std::sin(om2*t)};
    };

    const int n = 8;
    const double dt = 1.0;
    std::vector<V4>    ctrl(n);
    std::vector<double> times(n);
    for (int k = 0; k < n; ++k) {
        times[k] = k * dt;
        ctrl[k]  = torus_pt(times[k]);
    }

    bool used_clifford = false;
    auto result = fc::blend_curve<4>(ctrl, times, fc::nd_tag{},
                                     120, 2, &used_clifford);

    check(used_clifford, "T5a: used_clifford=true");

    // Measure maximum deviation from true torus curve
    double max_dev = 0.0;
    for (size_t i = 0; i < result.pts.size(); ++i) {
        double t   = result.times[i];
        V4 true_pt = torus_pt(t);
        double dev = (result.pts[i] - true_pt).norm();
        if (dev > max_dev) max_dev = dev;
    }

    char buf[64]; std::snprintf(buf, sizeof(buf), "(max_dev=%.2e)", max_dev);
    check(max_dev < 1e-13, "T5b: max_dev < 1e-13", buf);
}


// ── Test 6: FHWindow<2,3> exact knot interpolation on parabola ────────────
//
// FHWindow<2,3> (Floater-Hormann D=3) must interpolate exactly at all 5
// knots for any input.  Also verify that FHWindow<2,4> gives the same
// result as LagrangeWindow<2> (D=4 is a special case of FH = Lagrange).

static void test6()
{
    std::printf("\n-- T6: FHWindow barycentric rational fallback --\n");
    using V2 = fc::VecN<2>;

    // Parabola: p(t) = (t, t^2)
    auto par = [](double t) -> V2 { return V2{t, t*t}; };
    double ts[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
    V2 pts[5];
    for (int k = 0; k < 5; ++k) pts[k] = par(ts[k]);

    // T6a: FHWindow<2,3> interpolates exactly at all 5 knots
    fc::FHWindow<2,3> fh3(pts[0],pts[1],pts[2],pts[3],pts[4],
                           ts[0],ts[1],ts[2],ts[3],ts[4]);
    double max_err3 = 0.0;
    for (int k = 0; k < 5; ++k)
        max_err3 = std::max(max_err3, (fh3(ts[k]) - pts[k]).norm());
    char buf[64]; std::snprintf(buf, sizeof(buf), "(max_err=%.2e)", max_err3);
    check(max_err3 < 1e-12, "T6a: FHWindow<2,3> knots exact", buf);

    // T6b: FHWindow<2,4> matches LagrangeWindow<2> on 100 sample points
    fc::FHWindow<2,4>     fh4(pts[0],pts[1],pts[2],pts[3],pts[4],
                               ts[0],ts[1],ts[2],ts[3],ts[4]);
    fc::LagrangeWindow<2> lag(pts[0],pts[1],pts[2],pts[3],pts[4],
                               ts[0],ts[1],ts[2],ts[3],ts[4]);
    double max_diff = 0.0;
    for (int k = 0; k <= 100; ++k) {
        double t = ts[0] + k * (ts[4]-ts[0]) / 100.0;
        max_diff = std::max(max_diff, (fh4(t) - lag(t)).norm());
    }
    char buf2[64]; std::snprintf(buf2, sizeof(buf2), "(max_diff=%.2e)", max_diff);
    check(max_diff < 1e-12, "T6b: FHWindow<2,4> == LagrangeWindow<2>", buf2);

    // T6c: blend_curve<2, FHWindow3> runs without errors (smoke test)
    std::vector<V2>     ctrl6(8);
    std::vector<double> times6(8);
    for (int k = 0; k < 8; ++k) { times6[k] = k; ctrl6[k] = par(k); }
    auto r = fc::blend_curve<2, fc::FHWindow3>(ctrl6, times6, fc::conic_tag{}, 10, 2);
    check(!r.pts.empty(), "T6c: blend_curve<2,FHWindow3> produces output");
}


// ── Main ───────────────────────────────────────────────────────────────────

int main()
{
    std::printf("=== test_clifford_nd ===\n");

    test1();
    test2();
    test3();
    test4();
    test5();
    test6();

    std::printf("\n%d/%d passed\n", n_pass, n_pass + n_fail);
    return (n_fail == 0) ? 0 : 1;
}
