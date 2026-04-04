// demo_nd_plot.cpp — output CSV data for demo_nd_clifford.py
//
// Sections printed to stdout:
//   CTRL    — n=16 control points on torus helix (dt=1)
//   ND      — blend_curve<4>(nd_tag)   dense output
//   CONIC   — blend_curve<4>(conic_tag) dense output
//   TRUE    — true torus curve (dense, same t range)
//   ERRORS  — max knot error vs dt for nd_tag and conic_tag
//
// Build:
//   g++ -std=c++17 -O2 -I. \
//       -isystem /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1 \
//       -o demo_nd_plot demo_nd_plot.cpp

#include "conicblend_nd.hpp"
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>

using V4 = fc::VecN<4>;

static const double R1 = 2.0, R2 = 1.0;
static const double OM1 = 1.0, OM2 = 1.4142135623730951; // sqrt(2)

static V4 torus_pt(double t) {
    return V4{R1*std::cos(OM1*t), R1*std::sin(OM1*t),
              R2*std::cos(OM2*t), R2*std::sin(OM2*t)};
}

// Measure max knot error for a set of windows over n control points
static double knot_error_nd(std::vector<V4> const& c, std::vector<double> const& t)
{
    int n = (int)c.size();
    double err = 0;
    for (int i = 0; i < n-4; ++i) {
        fc::CliffordWindow<4> win(c[i],c[i+1],c[i+2],c[i+3],c[i+4],
                                   t[i],t[i+1],t[i+2],t[i+3],t[i+4]);
        if (!win.valid()) continue;
        for (int k = 0; k < 5; ++k)
            err = std::max(err, (win(t[i+k]) - c[i+k]).norm());
    }
    return err;
}

static double knot_error_conic(std::vector<V4> const& c, std::vector<double> const& t)
{
    int n = (int)c.size();
    double err = 0;
    for (int i = 0; i < n-4; ++i) {
        fc::ConicWindow<4> win(c[i],c[i+1],c[i+2],c[i+3],c[i+4],
                                t[i],t[i+1],t[i+2],t[i+3],t[i+4]);
        if (!win.valid()) continue;
        for (int k = 0; k < 5; ++k)
            err = std::max(err, (win(t[i+k]) - c[i+k]).norm());
    }
    return err;
}

int main()
{
    // ── Control points ────────────────────────────────────────────────────
    const int    N  = 16;
    const double DT = 1.0;
    std::vector<V4>     ctrl(N);
    std::vector<double> times(N);
    for (int i = 0; i < N; ++i) {
        times[i] = i * DT;
        ctrl[i]  = torus_pt(times[i]);
    }

    std::printf("CTRL\nt,x,y,z,w\n");
    for (int i = 0; i < N; ++i)
        std::printf("%.6f,%.8f,%.8f,%.8f,%.8f\n",
                    times[i], ctrl[i][0], ctrl[i][1], ctrl[i][2], ctrl[i][3]);

    // ── nd_tag blend ──────────────────────────────────────────────────────
    auto r_nd = fc::blend_curve<4>(ctrl, times, fc::nd_tag{}, 120, 2);
    std::printf("ND\nt,x,y,z,w\n");
    for (size_t i = 0; i < r_nd.pts.size(); ++i)
        std::printf("%.8f,%.8f,%.8f,%.8f,%.8f\n",
                    r_nd.times[i],
                    r_nd.pts[i][0], r_nd.pts[i][1],
                    r_nd.pts[i][2], r_nd.pts[i][3]);

    // ── conic_tag blend ───────────────────────────────────────────────────
    auto r_cn = fc::blend_curve<4>(ctrl, times, fc::conic_tag{}, 120, 2);
    std::printf("CONIC\nt,x,y,z,w\n");
    for (size_t i = 0; i < r_cn.pts.size(); ++i)
        std::printf("%.8f,%.8f,%.8f,%.8f,%.8f\n",
                    r_cn.times[i],
                    r_cn.pts[i][0], r_cn.pts[i][1],
                    r_cn.pts[i][2], r_cn.pts[i][3]);

    // ── True torus curve (dense) ──────────────────────────────────────────
    std::printf("TRUE\nt,x,y,z,w\n");
    double t_lo = times[2], t_hi = times[N-3];
    int n_true = 600;
    for (int i = 0; i <= n_true; ++i) {
        double t = t_lo + (t_hi - t_lo) * i / n_true;
        V4 p = torus_pt(t);
        std::printf("%.8f,%.8f,%.8f,%.8f,%.8f\n", t, p[0], p[1], p[2], p[3]);
    }

    // ── Knot error vs dt ──────────────────────────────────────────────────
    std::printf("ERRORS\ndt,nd_err,conic_err\n");
    double dt_vals[] = {0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.3, 1.6, 2.0};
    for (double d : dt_vals) {
        std::vector<V4>     c(N);
        std::vector<double> t(N);
        for (int i = 0; i < N; ++i) { t[i] = i*d; c[i] = torus_pt(t[i]); }
        double nd_e  = knot_error_nd(c, t);
        double cn_e  = knot_error_conic(c, t);
        std::printf("%.4f,%.6e,%.6e\n", d, nd_e, cn_e);
    }

    return 0;
}
