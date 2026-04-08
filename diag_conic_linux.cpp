// diag_conic_linux.cpp — diagnose conic window failure on Linux g++
// Build: g++ -std=c++17 -O2 -I. -o diag_conic_linux diag_conic_linux.cpp
#include "conicblend.hpp"
#include <cstdio>
#include <cmath>

using V2 = fc::VecN<2>;

static const double PI = 3.14159265358979323846;

// Test: single ConicWindow on tilted ellipse window i=0 (pts 0..4)
int main()
{
    double ca = std::cos(PI/6.0), sa = std::sin(PI/6.0);
    int n = 16;
    std::vector<V2> ctrl(n);
    std::vector<double> times(n);
    for (int i = 0; i < n; ++i) {
        double t = 2.0*PI * i / (n-1);
        double x = 3.0*std::cos(t), y = 1.0*std::sin(t);
        ctrl[i]  = V2(ca*x - sa*y, sa*x + ca*y);
        times[i] = t;
    }

    // Print all pts for reference
    std::printf("=== tilted ellipse pts ===\n");
    for (int i = 0; i < n; ++i)
        std::printf("  ctrl[%2d] t=%.6f  (%.8f, %.8f)\n",
                    i, times[i], ctrl[i][0], ctrl[i][1]);

    // Test each window
    std::printf("\n=== window validity ===\n");
    for (int w = 0; w < n-4; ++w) {
        fc::ConicWindow<2> win(ctrl[w],ctrl[w+1],ctrl[w+2],ctrl[w+3],ctrl[w+4],
                               times[w],times[w+1],times[w+2],times[w+3],times[w+4]);
        std::printf("  win[%2d]: valid=%d  fit_err=%.3e  c_infty=%d\n",
                    w, (int)win.valid(), win.fit_error(), (int)win.uses_c_infty());
    }

    // Try with -ffast-math effect: same under ffast?
    std::printf("\nDone.\n");
}
