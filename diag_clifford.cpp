// diag_clifford.cpp — stress-test CliffordWindow on known Clifford tori
//
// Generates 5-point windows on exact Clifford tori with varied parameters
// and reports whether the solver finds them (used_clifford=true, knot error < 1e-12).
//
// Build:
//   g++ -std=c++17 -O2 -I. \
//       -isystem /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1 \
//       -o diag_clifford diag_clifford.cpp && ./diag_clifford

#include "conicblend_nd.hpp"
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>

using V4 = fc::VecN<4>;

static int n_pass = 0, n_fail = 0;

static void probe(const char* label,
                  double r1, double r2, double om1, double om2,
                  double ts[5])
{
    // Generate exact Clifford torus points
    V4 pts[5];
    for (int k = 0; k < 5; ++k) {
        double t = ts[k];
        pts[k] = V4{r1*std::cos(om1*t), r1*std::sin(om1*t),
                    r2*std::cos(om2*t), r2*std::sin(om2*t)};
    }

    fc::CliffordWindow<4> win(pts[0],pts[1],pts[2],pts[3],pts[4],
                               ts[0], ts[1], ts[2], ts[3], ts[4]);

    bool clf  = win.used_clifford();
    double max_err = 0.0;
    if (win.valid()) {
        for (int k = 0; k < 5; ++k)
            max_err = std::max(max_err, (win(ts[k]) - pts[k]).norm());
    }

    bool ok = clf && (max_err < 1e-10);
    if (ok) ++n_pass; else ++n_fail;

    std::printf("  %-46s  clf=%-3s  knot_err=%.2e  %s\n",
                label, clf ? "yes" : "NO", max_err, ok ? "PASS" : "FAIL");
}

int main()
{
    std::printf("=== diag_clifford: Clifford torus solver stress test ===\n\n");

    // ── 1. Baseline: uniform spacing ─────────────────────────────────────────
    std::printf("-- Baseline (uniform dt=1.0) --\n");
    {
        double ts[5] = {0,1,2,3,4};
        probe("r1=2 r2=1 om1=1 om2=sqrt2",       2,1, 1, std::sqrt(2.0), ts);
        probe("r1=1 r2=1 om1=1 om2=sqrt2",       1,1, 1, std::sqrt(2.0), ts);
        probe("r1=3 r2=0.5 om1=1 om2=sqrt2",     3,0.5, 1, std::sqrt(2.0), ts);
        probe("r1=0.5 r2=3 om1=1 om2=sqrt2",     0.5,3, 1, std::sqrt(2.0), ts);
        probe("r1=2 r2=1 om1=0.5 om2=0.3",       2,1, 0.5, 0.3, ts);
        probe("r1=2 r2=1 om1=2.0 om2=sqrt2",     2,1, 2.0, std::sqrt(2.0), ts);
        probe("r1=2 r2=1 om1=3.0 om2=sqrt2",     2,1, 3.0, std::sqrt(2.0), ts);
    }

    // ── 2. Equal radii (r1=r2): moment matrix has 4-fold degenerate eigenvalue ─
    std::printf("\n-- Equal radii r1=r2 --\n");
    {
        double ts[5] = {0,1,2,3,4};
        probe("r1=r2=1.0 om1=1 om2=sqrt2",   1,1, 1, std::sqrt(2.0), ts);
        probe("r1=r2=2.0 om1=1 om2=2.5",     2,2, 1, 2.5, ts);
        probe("r1=r2=1.5 om1=0.7 om2=1.3",   1.5,1.5, 0.7, 1.3, ts);
    }

    // ── 3. cos(om1*dt) ≈ cos(om2*dt): ESPRIT eigenvalues nearly degenerate ───
    std::printf("\n-- Eigenvalue near-degeneracy: cos(om1) ≈ cos(om2) --\n");
    {
        double ts[5] = {0,1,2,3,4};
        // cos(1.0) ≈ 0.5403;  cos(2π-1) = cos(5.283) ≈ 0.5403  → degenerate!
        probe("om1=1.0 om2=2pi-1=5.283 (exact degen)", 2,1, 1.0, 2*M_PI-1.0, ts);
        // near-degenerate
        probe("om1=1.0 om2=5.0 (cos≈0.28 vs 0.54)",    2,1, 1.0, 5.0, ts);
        probe("om1=1.0 om2=5.1 (near 2pi-1)",           2,1, 1.0, 5.1, ts);
        probe("om1=1.0 om2=5.2 (closer to 2pi-1)",      2,1, 1.0, 5.2, ts);
    }

    // ── 4. Non-uniform time spacing ──────────────────────────────────────────
    std::printf("\n-- Non-uniform time spacing --\n");
    {
        double ts_mild[5]   = {0, 0.8, 1.7, 2.9, 4.0};   // mild non-uniform
        double ts_heavy[5]  = {0, 0.1, 2.0, 2.3, 4.0};   // heavy clustering
        double ts_jitter[5] = {0, 1.1, 1.9, 3.2, 3.8};   // jittered
        double ts_log[5]    = {0, 0.2, 0.6, 1.4, 4.0};   // logarithmic-like
        probe("mild nonuniform  {0,.8,1.7,2.9,4}",   2,1, 1, std::sqrt(2.0), ts_mild);
        probe("heavy clustering {0,.1,2,2.3,4}",     2,1, 1, std::sqrt(2.0), ts_heavy);
        probe("jittered         {0,1.1,1.9,3.2,3.8}",2,1, 1, std::sqrt(2.0), ts_jitter);
        probe("log-spaced       {0,.2,.6,1.4,4}",    2,1, 1, std::sqrt(2.0), ts_log);
    }

    // ── 5. Large angular steps (wrapping) ────────────────────────────────────
    std::printf("\n-- Large angular steps (possible wrap) --\n");
    {
        double ts[5] = {0,1,2,3,4};
        probe("om1=pi-0.1 ≈ 3.04 (near half-turn)", 2,1, M_PI-0.1, 1.0, ts);
        probe("om1=pi     (exact half-turn)",         2,1, M_PI,     1.0, ts);
        probe("om1=pi+0.1 (just over half-turn)",    2,1, M_PI+0.1, 1.0, ts);
        probe("om2=pi-0.1",                           2,1, 1.0, M_PI-0.1, ts);
        probe("om1=2 om2=pi",                         2,1, 2.0, M_PI,    ts);
    }

    // ── 6. Very small / very large radii ratio ───────────────────────────────
    std::printf("\n-- Extreme radii ratios --\n");
    {
        double ts[5] = {0,1,2,3,4};
        probe("r1=10 r2=0.1",   10, 0.1, 1, std::sqrt(2.0), ts);
        probe("r1=0.1 r2=10",   0.1, 10, 1, std::sqrt(2.0), ts);
        probe("r1=100 r2=1",    100, 1,  1, std::sqrt(2.0), ts);
        probe("r1=1 r2=100",    1, 100,  1, std::sqrt(2.0), ts);
    }

    // ── 7. Small dt (nearly collinear points) ────────────────────────────────
    std::printf("\n-- Very small dt (near-flat window) --\n");
    {
        double ts_tiny[5]   = {0, 0.01, 0.02, 0.03, 0.04};
        double ts_medium[5] = {0, 0.1,  0.2,  0.3,  0.4};
        probe("dt=0.01",  2,1, 1, std::sqrt(2.0), ts_tiny);
        probe("dt=0.1",   2,1, 1, std::sqrt(2.0), ts_medium);
    }

    // ── 8. Large dt (nearly full-circle steps) ───────────────────────────────
    std::printf("\n-- Large dt (near full-period) --\n");
    {
        // om1*dt ≈ 2pi → points nearly repeat
        double ts_full[5] = {0, 2*M_PI/1.0-0.1, 2*2*M_PI/1.0-0.2, 3*2*M_PI/1.0, 4*2*M_PI/1.0};
        probe("om1*dt≈2pi (near repeat)", 2,1, 1.0, std::sqrt(2.0), ts_full);
    }

    std::printf("\n%d/%d passed\n", n_pass, n_pass+n_fail);
    return (n_fail > 0) ? 1 : 0;
}
