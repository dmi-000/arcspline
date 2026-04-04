// conicblend_nd.hpp — C^N n-dimensional curve interpolation via Clifford torus
// Header-only C++17 library.  Requires conicblend_cylinder.hpp.
//
// Generalises the conic-on-cylinder architecture to any N≥4:
//   5 points in ℝᴺ always span at most a 4D affine subspace.
//   In that 4D hull the Clifford torus S¹(r₁)×S¹(r₂) has exactly 10 DOF = 5×2
//   constraints, matching the conic-on-cylinder DOF count in ℝ³.
//   Solver: project to 4D hull → linear sphere fit → moment-matrix eigendecomposition
//           → Grassmannian gradient refinement → unroll (θ₁,θ₂) → ConicWindow<2>.
//
// Fallback chain per window (CliffordWindow<N>, N≥4):
//   1. Clifford torus in 4D affine hull        (hull rank = 4)
//   2. CylinderWindow<3> in 3D sub-hull         (hull rank = 3)
//   3. ConicWindow<2>    in 2D sub-hull          (hull rank = 2)
//   4. LagrangeWindow<N>                         (CliffordWindow::valid() == false)
//
// Usage:
//   auto r = fc::blend_curve<4>(ctrl, times, fc::nd_tag{});
//
// Compile:
//   g++ -std=c++17 -O2 -I.
//       -isystem /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1

#pragma once
#include "conicblend_cylinder.hpp"

namespace fc {

struct nd_tag {};

namespace detail {

// ── Affine hull of 5 points in ℝᴺ ────────────────────────────────────────
//
// Returns an orthonormal basis for the affine hull of pts[0..4].
// rank = number of independent directions found (0–4).
// coords[k][j] = coordinate of pts[k] along basis[j]  (j < rank; 0 otherwise).
// origin = pts[0] (all coordinates are relative to pts[0]).
//
// Algorithm: Gram-Schmidt on the 4 difference vectors d[k]=pts[k+1]-pts[0].
// A vector is accepted if its residual after projection exceeds 1e-10·|d[k]|.

template <int N>
struct Hull4D {
    VecN<N> origin;       // pts[0]
    VecN<N> basis[4];     // orthonormal; first `rank` valid
    int     rank;         // 0–4
    double  coords[5][4]; // coords of pts[k] in hull basis
};

template <int N>
inline Hull4D<N> affine_hull_4d(VecN<N> const pts[5])
{
    Hull4D<N> h;
    h.origin = pts[0];
    h.rank = 0;
    for (int k = 0; k < 5; ++k)
        for (int j = 0; j < 4; ++j)
            h.coords[k][j] = 0.0;

    for (int di = 0; di < 4; ++di) {
        VecN<N> v = pts[di + 1] - pts[0];
        double orig_norm = v.norm();
        if (orig_norm < 1e-14) continue;          // coincident points

        // Orthogonalize v against current basis
        for (int j = 0; j < h.rank; ++j)
            v = v - h.basis[j] * h.basis[j].dot(v);

        double vn = v.norm();
        if (vn < 1e-10 * orig_norm) continue;     // nearly in existing subspace

        h.basis[h.rank] = v * (1.0 / vn);
        h.rank++;
    }

    // Project all 5 points onto hull basis
    for (int k = 0; k < 5; ++k) {
        VecN<N> q = pts[k] - pts[0];
        for (int j = 0; j < h.rank; ++j)
            h.coords[k][j] = h.basis[j].dot(q);
        for (int j = h.rank; j < 4; ++j)
            h.coords[k][j] = 0.0;
    }

    return h;
}


// ── Circumscribed sphere fit through 5 points in ℝ⁴ (linear) ────────────
//
// All 5 points on a Clifford torus lie on a sphere of radius √(r₁²+r₂²).
// System: (p_k − p_0)·c = (|p_k|² − |p_0|²)/2  for k=1..4.
// 4×4 Gaussian elimination with partial pivoting.
// Returns false if the system is (near-)degenerate.

inline bool sphere_fit_4d(VecN<4> const pts[5], VecN<4>& center, double& R)
{
    // Augmented matrix [A | b], 4×5
    double A[4][5];
    for (int k = 1; k <= 4; ++k) {
        for (int j = 0; j < 4; ++j)
            A[k-1][j] = pts[k][j] - pts[0][j];
        double pk2 = pts[k].dot(pts[k]);
        double p02 = pts[0].dot(pts[0]);
        A[k-1][4] = (pk2 - p02) * 0.5;
    }

    // Gaussian elimination with partial pivoting
    for (int col = 0; col < 4; ++col) {
        int piv = col;
        for (int row = col + 1; row < 4; ++row)
            if (std::abs(A[row][col]) > std::abs(A[piv][col])) piv = row;
        if (std::abs(A[piv][col]) < 1e-14) return false;
        if (piv != col)
            for (int k = 0; k <= 4; ++k) std::swap(A[col][k], A[piv][k]);
        for (int row = col + 1; row < 4; ++row) {
            double f = A[row][col] / A[col][col];
            for (int k = col; k <= 4; ++k) A[row][k] -= f * A[col][k];
        }
    }

    // Back substitution → c
    double c[4];
    for (int row = 3; row >= 0; --row) {
        c[row] = A[row][4];
        for (int k = row + 1; k < 4; ++k) c[row] -= A[row][k] * c[k];
        c[row] /= A[row][row];
    }
    for (int j = 0; j < 4; ++j) center[j] = c[j];

    // R = RMS distance (should equal √(r₁²+r₂²) for exact Clifford torus data)
    R = 0.0;
    for (int k = 0; k < 5; ++k) R += (pts[k] - center).norm2();
    R = std::sqrt(R / 5.0);
    return R > 1e-12;
}


// Orthonormal complement of span(f1,f2) in ℝ⁴ → f3,f4.  Returns false if
// f1 or f2 are degenerate.  Used both inside clifford_best_fit_4d and in
// try_clifford_ when exploring all 6 eigenspace splits.
inline bool make_complement_4d(VecN<4> const& f1, VecN<4> const& f2,
                                VecN<4>& f3, VecN<4>& f4)
{
    const VecN<4> stdbasis[4] = {VecN<4>{1,0,0,0}, VecN<4>{0,1,0,0},
                                  VecN<4>{0,0,1,0}, VecN<4>{0,0,0,1}};
    int found = 0;
    f3 = VecN<4>{}; f4 = VecN<4>{};
    for (auto const& b : stdbasis) {
        if (found == 2) break;
        VecN<4> v = b;
        v = v - f1 * f1.dot(v) - f2 * f2.dot(v);
        if (found == 1) v = v - f3 * f3.dot(v);
        double vn = v.norm();
        if (vn < 1e-10) continue;
        if (found == 0) f3 = v * (1.0 / vn);
        else            f4 = v * (1.0 / vn);
        ++found;
    }
    return found == 2;
}

// ── Clifford torus best-fit solver ────────────────────────────────────────
//
// Given 5 points in ℝ⁴, fits the Clifford torus S¹(r₁)×S¹(r₂) of least
// variance, returning the frame (e1,e2,e3,e4), radii, and quality metric.
//
// Algorithm:
//   1. sphere_fit_4d → center c, R; q_k = pts[k] − c
//   2. Moment matrix M = Σ q_k q_kᵀ; SymEig<4>::compute
//   3. Try all 3 eigenspace 2+2 splits; pick best max_surf_res candidate
//   4. Gradient descent on Gr(2,4) to minimise Var(|Π₁q_k|²)
//   5. Compute e3,e4 as orthogonal complement of span(e1,e2)
//
// max_surf_res = max_k ||Π₁q_k|² − r₁²| / R²  (dimensionless, 1e30 on failure)

struct CliffordBFSol4D {
    VecN<4> e1, e2;          // P₁ = span(e1,e2)
    VecN<4> e3, e4;          // P₂ = P₁⊥ = span(e3,e4)
    VecN<4> center;
    double  r1, r2;
    double  max_surf_res;    // 1e30 if fit failed
    VecN<4> ev[4];           // eigenvectors of moment matrix (ascending eigenvalue order)
};

inline CliffordBFSol4D clifford_best_fit_4d(VecN<4> const pts[5])
{
    CliffordBFSol4D sol;
    sol.max_surf_res = 1e30;
    sol.r1 = 0.0; sol.r2 = 0.0;

    // Step 1: circumscribed sphere
    VecN<4> center;
    double R;
    if (!sphere_fit_4d(pts, center, R)) return sol;
    double R2 = R * R;

    VecN<4> q[5];
    for (int k = 0; k < 5; ++k) q[k] = pts[k] - center;

    // Step 2: moment matrix M = Σ q_k q_kᵀ (copy before passing to SymEig)
    double Mcopy[4][4] = {};
    for (int k = 0; k < 5; ++k)
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                Mcopy[i][j] += q[k][i] * q[k][j];

    SymEig<4> eig;
    eig.compute(Mcopy);  // ascending order; Mcopy modified in place

    // Extract eigenvectors
    VecN<4> ev[4];
    for (int j = 0; j < 4; ++j)
        for (int i = 0; i < 4; ++i)
            ev[j][i] = eig.vec[j][i];

    // Helper: compute r1² and max_surf_res for a split P1 = span(f1,f2)
    auto eval_split = [&](VecN<4> const& f1, VecN<4> const& f2,
                          double& r1sq_out, double& msr_out)
    {
        r1sq_out = 0.0;
        for (int k = 0; k < 5; ++k) {
            double a1 = f1.dot(q[k]), a2 = f2.dot(q[k]);
            r1sq_out += a1*a1 + a2*a2;
        }
        r1sq_out /= 5.0;
        msr_out = 0.0;
        for (int k = 0; k < 5; ++k) {
            double a1 = f1.dot(q[k]), a2 = f2.dot(q[k]);
            double d = std::abs(a1*a1 + a2*a2 - r1sq_out) / (R2 + 1e-30);
            if (d > msr_out) msr_out = d;
        }
    };

    // Steps 3–4: multi-start gradient descent on Gr(2,4)
    //
    // Starting points: all 3 eigenspace 2+2 splits, plus 20 random 2-frames.
    // For each start, run Riemannian gradient descent on the cost
    //   F = Σ_k (|Π₁q_k|² − r₁²)²
    // (r₁² recomputed each step as the mean).  Keep the frame with the lowest
    // final cost.  This multi-start strategy reliably finds the global minimum
    // (F = 0) for exact Clifford torus data even when the moment matrix is
    // poorly conditioned (large off-diagonal cross-terms).

    // F (circle residual): Σ(|Π_P1 q_k|² − r1²)²
    auto eval_cost = [&](VecN<4> const& f1, VecN<4> const& f2, double r1sq_in) {
        double cost = 0.0;
        for (int k = 0; k < 5; ++k) {
            double a1 = f1.dot(q[k]), a2 = f2.dot(q[k]);
            double diff = a1*a1 + a2*a2 - r1sq_in;
            cost += diff * diff;
        }
        return cost;
    };

    // Equal-chord regularisation: G'' = Var(|Π_P1 dq_k|²), dq_k = q[k+1]−q[k].
    // G'' = 0 for the correct Clifford torus (constant-angular-velocity ↔ equal
    // chord lengths for equal time steps).  G'' > 0 for spurious solutions.
    // Weight eps_G is small enough not to disturb convergence for non-Clifford
    // curves, and large enough to break the F=0 degeneracy.
    const double eps_G = 1e-6;
    VecN<4> dq[4];
    for (int k = 0; k < 4; ++k) dq[k] = q[k+1] - q[k];

    auto eval_G_chord = [&](VecN<4> const& f1, VecN<4> const& f2) {
        double ck[4]; double cmean = 0.0;
        for (int k = 0; k < 4; ++k) {
            double c1 = f1.dot(dq[k]), c2 = f2.dot(dq[k]);
            ck[k] = c1*c1 + c2*c2; cmean += ck[k];
        }
        cmean /= 4.0;
        double G = 0.0;
        for (int k = 0; k < 4; ++k) { double d = ck[k]-cmean; G += d*d; }
        return G;
    };

    // Combined objective H = F + eps_G * G''.  H = 0 uniquely at the correct
    // Clifford torus (both circle and equal-chord conditions satisfied).
    auto eval_H = [&](VecN<4> const& f1, VecN<4> const& f2, double r1sq_in) {
        return eval_cost(f1, f2, r1sq_in) + eps_G * eval_G_chord(f1, f2);
    };

    // Gradient descent from a given 2-frame on H; modifies (f1,f2) in place.
    auto grad_descent = [&](VecN<4>& f1, VecN<4>& f2) {
        double r1sq = 0.0;
        for (int k = 0; k < 5; ++k) {
            double a1 = f1.dot(q[k]), a2 = f2.dot(q[k]);
            r1sq += a1*a1 + a2*a2;
        }
        r1sq /= 5.0;
        for (int step = 0; step < 60; ++step) {
            VecN<4> g1{}, g2{};
            // ∂F/∂f: 4·Σ(|Πq|²−r1²)·(f·q)·q
            for (int k = 0; k < 5; ++k) {
                double a1 = f1.dot(q[k]), a2 = f2.dot(q[k]);
                double fk = a1*a1 + a2*a2 - r1sq;
                g1 = g1 + q[k] * (4.0 * fk * a1);
                g2 = g2 + q[k] * (4.0 * fk * a2);
            }
            // ∂G''/∂f: 4·eps_G·Σ(ck−cmean)·(f·dq_k)·dq_k
            {
                double ck[4]; double cmean = 0.0;
                for (int k = 0; k < 4; ++k) {
                    double c1 = f1.dot(dq[k]), c2 = f2.dot(dq[k]);
                    ck[k] = c1*c1 + c2*c2; cmean += ck[k];
                }
                cmean /= 4.0;
                for (int k = 0; k < 4; ++k) {
                    double diff = ck[k] - cmean;
                    double c1 = f1.dot(dq[k]), c2 = f2.dot(dq[k]);
                    g1 = g1 + dq[k] * (4.0 * eps_G * diff * c1);
                    g2 = g2 + dq[k] * (4.0 * eps_G * diff * c2);
                }
            }
            // Tangent projection
            VecN<4> d1 = g1 - f1*(f1.dot(g1)) - f2*(f2.dot(g1));
            VecN<4> d2 = g2 - f1*(f1.dot(g2)) - f2*(f2.dot(g2));
            double gnorm = std::sqrt(d1.norm2() + d2.norm2());
            if (gnorm < 1e-15) break;
            double alpha = 0.1 / gnorm;
            double H0 = eval_H(f1, f2, r1sq);
            bool improved = false;
            for (int ls = 0; ls < 16; ++ls) {
                VecN<4> ne1 = f1 - d1 * alpha;
                VecN<4> ne2 = f2 - d2 * alpha;
                double n1 = ne1.norm();
                if (n1 < 1e-12) { alpha *= 0.5; continue; }
                ne1 = ne1 * (1.0 / n1);
                ne2 = ne2 - ne1 * ne1.dot(ne2);
                double n2 = ne2.norm();
                if (n2 < 1e-12) { alpha *= 0.5; continue; }
                ne2 = ne2 * (1.0 / n2);
                double r1sq_new = 0.0;
                for (int k = 0; k < 5; ++k) {
                    double a1 = ne1.dot(q[k]), a2 = ne2.dot(q[k]);
                    r1sq_new += a1*a1 + a2*a2;
                }
                r1sq_new /= 5.0;
                double H1 = eval_H(ne1, ne2, r1sq_new);
                if (H1 < H0) {
                    f1 = ne1; f2 = ne2; r1sq = r1sq_new;
                    improved = true; break;
                }
                alpha *= 0.5;
            }
            if (!improved) break;
        }
    };

    VecN<4> e1{}, e2{};
    double  r1sq = 0.0;
    double  best_cost = 1e30;

    // Step 3: try all 6 eigenspace 2-frame starting points.
    // The 3 distinct 2+2 partitions of {0,1,2,3} give 6 starting 2-frames
    // (both halves of each partition): {01},{23},{02},{13},{03},{12}.
    const int splits[6][2] = {{0,1},{2,3},{0,2},{1,3},{0,3},{1,2}};
    for (int s = 0; s < 6; ++s) {
        VecN<4> f1 = ev[splits[s][0]], f2 = ev[splits[s][1]];
        grad_descent(f1, f2);
        double rs = 0.0;
        for (int k = 0; k < 5; ++k) {
            double a1 = f1.dot(q[k]), a2 = f2.dot(q[k]);
            rs += a1*a1 + a2*a2;
        }
        rs /= 5.0;
        double cost = eval_H(f1, f2, rs);
        if (cost < best_cost) { best_cost = cost; e1 = f1; e2 = f2; r1sq = rs; }
    }

    // Step 4: random starts on Gr(2,4) — ensures global minimum is found
    // Deterministic LCG (no <random> header needed).
    uint64_t rng = 0xdeadbeef42ULL;
    auto rng_step = [&]() {
        rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
        return rng;
    };
    auto rng_f = [&]() {   // uniform in (-1,1)
        return (static_cast<double>(rng_step() & 0xFFFF) / 32767.5) - 1.0;
    };

    const int NUM_RAND = 24;
    for (int rs_idx = 0; rs_idx < NUM_RAND; ++rs_idx) {
        VecN<4> f1, f2;
        for (int i = 0; i < 4; ++i) f1[i] = rng_f();
        for (int i = 0; i < 4; ++i) f2[i] = rng_f();
        double n1 = f1.norm(); if (n1 < 1e-12) continue;
        f1 = f1 * (1.0 / n1);
        f2 = f2 - f1 * f1.dot(f2);
        double n2 = f2.norm(); if (n2 < 1e-12) continue;
        f2 = f2 * (1.0 / n2);
        grad_descent(f1, f2);
        double rs = 0.0;
        for (int k = 0; k < 5; ++k) {
            double a1 = f1.dot(q[k]), a2 = f2.dot(q[k]);
            rs += a1*a1 + a2*a2;
        }
        rs /= 5.0;
        double cost = eval_H(f1, f2, rs);
        if (cost < best_cost) { best_cost = cost; e1 = f1; e2 = f2; r1sq = rs; }
        if (best_cost < 1e-15) break;   // exact fit found — no need for more starts
    }

    // Helper: compute orthonormal complement e3,e4 of span(e1,e2)
    auto make_complement = [](VecN<4> const& f1, VecN<4> const& f2,
                               VecN<4>& f3, VecN<4>& f4) -> bool {
        VecN<4> stdbasis[4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
        int found = 0;
        f3 = {}; f4 = {};
        for (auto& b : stdbasis) {
            if (found == 2) break;
            VecN<4> v = b;
            v = v - f1 * f1.dot(v) - f2 * f2.dot(v);
            if (found == 1) v = v - f3 * f3.dot(v);
            double vn = v.norm();
            if (vn < 1e-10) continue;
            if (found == 0) f3 = v * (1.0 / vn);
            else            f4 = v * (1.0 / vn);
            ++found;
        }
        return found == 2;
    };

    // Orthogonal complement of span(e1,e2)
    VecN<4> e3{}, e4{};
    if (!make_complement(e1, e2, e3, e4)) return sol;

    // Newton refinement: solve d_k^T P s_k = 0 exactly.
    //
    // For exact Clifford torus data, d_k^T P s_k = 0 iff P = e1e1^T+e2e2^T
    // (the true frame) or P = I-P (complementary frame).  Gradient descent
    // converges to a nearby near-minimum; Newton steps bring it to machine
    // precision in 3–5 iterations (quadratic convergence).
    //
    // Parameterise tangent perturbation as:
    //   δe1 = β₁ e3 + β₂ e4   (preserves e1·e1 to first order)
    //   δe2 = γ₁ e3 + γ₂ e4
    // Residuals:  ε_k = (e1·d_k)(e1·s_k) + (e2·d_k)(e2·s_k)
    //   where d_k = q[k+1]−q[0],  s_k = q[k+1]+q[0],  k=0..3.
    // Jacobian J[k][{β₁,β₂,γ₁,γ₂}]:
    //   ∂ε_k/∂β₁ = (e1·s_k)(e3·d_k) + (e1·d_k)(e3·s_k)   (analogously for rest)
    {
        VecN<4> dv[4], sv[4];
        for (int k = 0; k < 4; ++k) { dv[k] = q[k+1]-q[0]; sv[k] = q[k+1]+q[0]; }

        for (int iter = 0; iter < 8; ++iter) {
            double u1[4], u2[4], v1_[4], v2_[4];
            double d3[4], d4_[4], s3[4], s4_[4];
            for (int k = 0; k < 4; ++k) {
                u1[k]  = e1.dot(dv[k]); u2[k]  = e2.dot(dv[k]);
                v1_[k] = e1.dot(sv[k]); v2_[k] = e2.dot(sv[k]);
                d3[k]  = e3.dot(dv[k]); d4_[k] = e4.dot(dv[k]);
                s3[k]  = e3.dot(sv[k]); s4_[k] = e4.dot(sv[k]);
            }
            double eps_[4]; double max_eps = 0;
            for (int k = 0; k < 4; ++k) {
                eps_[k] = u1[k]*v1_[k] + u2[k]*v2_[k];
                double ae = std::abs(eps_[k]);
                if (ae > max_eps) max_eps = ae;
            }
            if (max_eps < 1e-14) break;

            // Jacobian J (4×4): columns ∂/∂β1, ∂/∂β2, ∂/∂γ1, ∂/∂γ2
            double Ja[4][5];
            for (int k = 0; k < 4; ++k) {
                Ja[k][0] = u1[k]*s3[k]  + v1_[k]*d3[k];
                Ja[k][1] = u1[k]*s4_[k] + v1_[k]*d4_[k];
                Ja[k][2] = u2[k]*s3[k]  + v2_[k]*d3[k];
                Ja[k][3] = u2[k]*s4_[k] + v2_[k]*d4_[k];
                Ja[k][4] = -eps_[k];
            }
            // Gaussian elimination with partial pivoting → solve Ja δ = rhs
            bool ok = true;
            for (int col = 0; col < 4; ++col) {
                int piv = col;
                for (int row = col+1; row < 4; ++row)
                    if (std::abs(Ja[row][col]) > std::abs(Ja[piv][col])) piv = row;
                if (std::abs(Ja[piv][col]) < 1e-14) { ok = false; break; }
                if (piv != col) for (int c2 = 0; c2 <= 4; ++c2) std::swap(Ja[col][c2], Ja[piv][c2]);
                for (int row = col+1; row < 4; ++row) {
                    double f = Ja[row][col] / Ja[col][col];
                    for (int c2 = col; c2 <= 4; ++c2) Ja[row][c2] -= f * Ja[col][c2];
                }
            }
            if (!ok) break;
            double delta[4];
            for (int row = 3; row >= 0; --row) {
                delta[row] = Ja[row][4];
                for (int c2 = row+1; c2 < 4; ++c2) delta[row] -= Ja[row][c2] * delta[c2];
                delta[row] /= Ja[row][row];
            }

            // Update frame
            VecN<4> ne1 = e1 + e3*delta[0] + e4*delta[1];
            VecN<4> ne2 = e2 + e3*delta[2] + e4*delta[3];
            double n1 = ne1.norm(); if (n1 < 1e-12) break;
            ne1 = ne1 * (1.0/n1);
            ne2 = ne2 - ne1 * ne1.dot(ne2);
            double n2 = ne2.norm(); if (n2 < 1e-12) break;
            ne2 = ne2 * (1.0/n2);
            VecN<4> ne3{}, ne4{};
            if (!make_complement(ne1, ne2, ne3, ne4)) break;
            e1 = ne1; e2 = ne2; e3 = ne3; e4 = ne4;
        }
    }

    // Final r1sq, r2sq, max_surf_res
    r1sq = 0.0;
    for (int k = 0; k < 5; ++k) {
        double a1 = e1.dot(q[k]), a2 = e2.dot(q[k]);
        r1sq += a1*a1 + a2*a2;
    }
    r1sq /= 5.0;

    double msr_final = 0.0;
    for (int k = 0; k < 5; ++k) {
        double a1 = e1.dot(q[k]), a2 = e2.dot(q[k]);
        double d = std::abs(a1*a1 + a2*a2 - r1sq) / (R2 + 1e-30);
        if (d > msr_final) msr_final = d;
    }

    // r2sq = mean |Π₂ q_k|²
    double r2sq = 0.0;
    for (int k = 0; k < 5; ++k) {
        double a3 = e3.dot(q[k]), a4 = e4.dot(q[k]);
        r2sq += a3*a3 + a4*a4;
    }
    r2sq /= 5.0;

    sol.e1 = e1;  sol.e2 = e2;
    sol.e3 = e3;  sol.e4 = e4;
    sol.center = center;
    sol.r1 = std::sqrt(std::max(0.0, r1sq));
    sol.r2 = std::sqrt(std::max(0.0, r2sq));
    sol.max_surf_res = msr_final;
    for (int j = 0; j < 4; ++j) sol.ev[j] = ev[j];
    return sol;
}

} // namespace detail


// ── CliffordWindow<N>: 5-point Clifford-torus window (N≥4) ───────────────
//
// Constructor tries four fallback levels in order:
//   1. Clifford torus S¹(r₁)×S¹(r₂) in the 4D affine hull (hull rank ≥ 4)
//   2. CylinderWindow<3> on the 3D sub-hull projection          (hull rank ≥ 3)
//   3. ConicWindow<2>    on the 2D sub-hull projection          (hull rank ≥ 2)
//   4. valid() == false  → blend_curve uses LagrangeWindow<N>
//
// For levels 2 and 3, the sub-hull coordinates are computed in the orthonormal
// hull basis, so Euclidean geometry is preserved exactly.

template <int N>
class CliffordWindow {
    static_assert(N >= 4, "CliffordWindow requires N >= 4");

    bool valid_        = false;
    bool use_cylinder_ = false;  // true → level 2 (CylinderWindow<3>) was used
    bool use_planar_   = false;  // true → level 3 (ConicWindow<2>) was used
                                 // neither → level 1 (Clifford torus) was used

    // Affine hull (common to all levels)
    VecN<N> origin_{};
    VecN<N> basis_[4]{};
    int     hull_rank_ = 0;

    // Level 1: Clifford torus (coordinates in 4D hull basis)
    VecN<4> e1_{}, e2_{}, e3_{}, e4_{}, center4_{};
    double  r1_ = 0.0, r2_ = 0.0;
    VecN<2> inner_a_lin_{};   // linear LS slope:     uv(t) = a*t + b
    VecN<2> inner_b_lin_{};   // linear LS intercept

    // Direct formula (set when extract_and_probe wins, avoids angle-extraction
    // round-trip that is O(1/dt²)-conditioned for small angular steps).
    // p(t) = center4_ + Av_*cos(o1d_*t) + Bv_*sin(o1d_*t)
    //                 + Av2_*cos(o2d_*t) + Bv2_*sin(o2d_*t)
    bool    use_direct_ = false;
    VecN<4> Av_{}, Bv_{}, Av2_{}, Bv2_{};
    double  o1d_ = 0.0, o2d_ = 0.0;

    // Ambient direct formula: same GJ system solved in ℝᴺ using original pts,
    // bypassing lift() to eliminate O(1e-10) Gram-Schmidt reconstruction error
    // that occurs when basis_[3] is ill-conditioned (e.g. dt=0.01).
    bool    use_direct_amb_ = false;
    VecN<N> center_amb_{}, Av_amb_{}, Bv_amb_{}, Av2_amb_{}, Bv2_amb_{};

    // Level 2: CylinderWindow<3> on 3D hull coords
    CylinderWindow<3> cyl_win3d_;

    // Level 3: ConicWindow<2> on 2D hull coords
    ConicWindow<2> plane_win_;

    // Lift from (4D) hull coordinates to ambient ℝᴺ
    VecN<N> lift(double c0, double c1, double c2, double c3) const
    {
        VecN<N> r = origin_;
        if (hull_rank_ > 0) r = r + basis_[0] * c0;
        if (hull_rank_ > 1) r = r + basis_[1] * c1;
        if (hull_rank_ > 2) r = r + basis_[2] * c2;
        if (hull_rank_ > 3) r = r + basis_[3] * c3;
        return r;
    }

    // try_clifford_: attempt Clifford torus window.
    // pts4[5] are the hull-coordinate representations of the 5 control points.
    //
    // Strategy: the optimizer in clifford_best_fit_4d may converge to a wrong
    // Clifford frame (different F=0 minimum) when the moment matrix is
    // ill-conditioned.  We therefore try ALL 7 candidate frames:
    //   – the gradient-descent result (bf.e1, bf.e2)
    //   – all 6 eigenspace 2+2 splits from bf.ev[0..3]
    // For each candidate we fit a linear LS model uv(t) = a*t + b to the
    // unrolled angles.  For the CORRECT frame angles are exactly linear in t
    // (constant angular velocity), so the roundtrip error is ≈0.  For wrong
    // frames the angles are non-linear, the LS residuals are large, and the
    // roundtrip fails.  We accept the frame with the smallest roundtrip error.
    bool try_clifford_(VecN<4> const pts4[5], VecN<N> const pts_amb[5],
                       detail::CliffordBFSol4D const& bf,
                       double t0, double t1, double t2, double t3, double t4)
    {
        // Window spread for roundtrip gate
        double spread = 0.0;
        for (int i = 0; i < 5; ++i)
            for (int j = i+1; j < 5; ++j)
                spread = std::max(spread, (pts4[i]-pts4[j]).norm());
        if (spread < 1e-12) return false;

        const double ts5[5] = {t0, t1, t2, t3, t4};

        // Precompute q_k = pts4[k] - center (sphere center, same for all frames)
        VecN<4> q[5];
        for (int k = 0; k < 5; ++k) q[k] = pts4[k] - bf.center;

        // LS coefficients for time vector (same denominant for all frames)
        double sum_t = 0.0, sum_t2 = 0.0;
        for (int k = 0; k < 5; ++k) { sum_t += ts5[k]; sum_t2 += ts5[k]*ts5[k]; }
        double ls_det = sum_t2 * 5.0 - sum_t * sum_t;
        if (std::abs(ls_det) < 1e-14) return false;

        // try_one: given frame (f1,f2), returns max roundtrip error (1e30 on failure)
        // and writes a_lin, b_lin if successful.
        auto try_one = [&](VecN<4> const& f1, VecN<4> const& f2,
                           VecN<4>& o_f3, VecN<4>& o_f4,
                           double& o_r1, double& o_r2,
                           VecN<2>& o_a, VecN<2>& o_b) -> double
        {
            VecN<4> f3{}, f4{};
            if (!detail::make_complement_4d(f1, f2, f3, f4)) return 1e30;

            double r1sq = 0.0, r2sq = 0.0;
            for (int k = 0; k < 5; ++k) {
                double a1 = f1.dot(q[k]), a2 = f2.dot(q[k]);
                r1sq += a1*a1 + a2*a2;
                double a3 = f3.dot(q[k]), a4 = f4.dot(q[k]);
                r2sq += a3*a3 + a4*a4;
            }
            r1sq /= 5.0; r2sq /= 5.0;
            if (r1sq < 1e-24 || r2sq < 1e-24) return 1e30;
            double r1 = std::sqrt(r1sq), r2 = std::sqrt(r2sq);

            // Unroll angles
            double th1[5], th2[5];
            for (int k = 0; k < 5; ++k) {
                th1[k] = std::atan2(f2.dot(q[k]), f1.dot(q[k]));
                th2[k] = std::atan2(f4.dot(q[k]), f3.dot(q[k]));
            }
            for (int k = 1; k < 5; ++k) {
                double d = th1[k]-th1[k-1];
                while (d >  detail::PI) { d -= 2*detail::PI; th1[k] -= 2*detail::PI; }
                while (d < -detail::PI) { d += 2*detail::PI; th1[k] += 2*detail::PI; }
            }
            for (int k = 1; k < 5; ++k) {
                double d = th2[k]-th2[k-1];
                while (d >  detail::PI) { d -= 2*detail::PI; th2[k] -= 2*detail::PI; }
                while (d < -detail::PI) { d += 2*detail::PI; th2[k] += 2*detail::PI; }
            }

            // upts
            VecN<2> upts[5];
            for (int k = 0; k < 5; ++k)
                upts[k] = VecN<2>{r1 * th1[k], r2 * th2[k]};

            // Linear LS: uv(t) = a * t + b
            VecN<2> sum_u{}, sum_tu{};
            for (int k = 0; k < 5; ++k) {
                sum_u  = sum_u  + upts[k];
                sum_tu = sum_tu + upts[k] * ts5[k];
            }
            VecN<2> a = (sum_tu * 5.0 - sum_u * sum_t) * (1.0 / ls_det);
            VecN<2> b = (sum_u * sum_t2 - sum_tu * sum_t) * (1.0 / ls_det);

            // Max roundtrip error at knots
            double max_err = 0.0;
            for (int k = 0; k < 5; ++k) {
                VecN<2> uv = a * ts5[k] + b;
                double ph1 = uv[0] / r1, ph2 = uv[1] / r2;
                VecN<4> rec = bf.center
                            + f1 * (r1 * std::cos(ph1))
                            + f2 * (r1 * std::sin(ph1))
                            + f3 * (r2 * std::cos(ph2))
                            + f4 * (r2 * std::sin(ph2));
                double e = (rec - pts4[k]).norm();
                if (e > max_err) max_err = e;
            }

            o_f3 = f3; o_f4 = f4;
            o_r1 = r1; o_r2 = r2;
            o_a  = a;  o_b  = b;
            return max_err;
        };

        // Candidates: gradient-descent frame + 6 eigenspace splits
        const int splits[6][2] = {{0,1},{2,3},{0,2},{1,3},{0,3},{1,2}};

        double best_err = 1e30;
        VecN<4> best_f1{}, best_f2{}, best_f3{}, best_f4{};
        double  best_r1 = 0.0, best_r2 = 0.0;
        VecN<2> best_a{}, best_b{};
        VecN<4> best_center = bf.center;
        bool    best_use_direct = false;
        VecN<4> best_Av{}, best_Bv{}, best_Av2{}, best_Bv2{};
        double  best_o1d = 0.0, best_o2d = 0.0;

        auto probe = [&](VecN<4> const& f1, VecN<4> const& f2) {
            VecN<4> f3{}, f4{};
            double r1 = 0, r2 = 0;
            VecN<2> a{}, b{};
            double err = try_one(f1, f2, f3, f4, r1, r2, a, b);
            if (err < best_err) {
                best_err = err;
                best_f1 = f1; best_f2 = f2; best_f3 = f3; best_f4 = f4;
                best_r1 = r1; best_r2 = r2;
                best_a  = a;  best_b  = b;
                best_center = bf.center; // try_one uses bf.center
            }
        };

        probe(bf.e1, bf.e2);
        for (auto& s : splits) probe(bf.ev[s[0]], bf.ev[s[1]]);

        // ESPRIT propagator: A = Q1 * inv(Q0), S = (A + A^T)/2.
        // For uniformly-spaced t, S has eigenspaces = correct Clifford planes.
        {
            // Q0 columns = q[0..3], Q1 columns = q[1..4]
            // Build augmented [Q0 | I] for Gauss-Jordan inversion
            double aug[4][8];
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) aug[i][j] = q[j][i];
                for (int j = 0; j < 4; ++j) aug[i][j+4] = (i == j) ? 1.0 : 0.0;
            }
            bool ok = true;
            for (int col = 0; col < 4 && ok; ++col) {
                int piv = col;
                for (int row = col+1; row < 4; ++row)
                    if (std::abs(aug[row][col]) > std::abs(aug[piv][col])) piv = row;
                if (std::abs(aug[piv][col]) < 1e-12 * spread) { ok = false; break; }
                if (piv != col)
                    for (int k = 0; k < 8; ++k) std::swap(aug[col][k], aug[piv][k]);
                double s = 1.0 / aug[col][col];
                for (int k = 0; k < 8; ++k) aug[col][k] *= s;
                for (int row = 0; row < 4; ++row) {
                    if (row == col) continue;
                    double f = aug[row][col];
                    for (int k = 0; k < 8; ++k) aug[row][k] -= f * aug[col][k];
                }
            }
            if (ok) {
                // A = Q1 * Q0inv
                double A[4][4] = {};
                for (int i = 0; i < 4; ++i)
                    for (int j = 0; j < 4; ++j)
                        for (int kk = 0; kk < 4; ++kk)
                            A[i][j] += q[kk+1][i] * aug[kk][j+4];
                // S = (A + A^T) / 2
                double S[4][4];
                for (int i = 0; i < 4; ++i)
                    for (int j = 0; j < 4; ++j)
                        S[i][j] = 0.5 * (A[i][j] + A[j][i]);
                detail::SymEig<4> seig;
                seig.compute(S);
                VecN<4> sev[4];
                for (int j = 0; j < 4; ++j)
                    for (int i = 0; i < 4; ++i)
                        sev[j][i] = seig.vec[j][i];
                probe(sev[0], sev[1]);
                probe(sev[2], sev[3]);
            }
        }

        // ── Direct-formula refinement using best-found frequencies ───────────
        // After the initial angle-parameterized probes (try_one + ESPRIT), extract
        // ω₁ = |best_a[0]| / best_r1 and ω₂ = |best_a[1]| / best_r2, then solve
        // the 5×5 GJ system and apply mean-residual center refinement.
        //
        // This bypasses the angle-extraction round-trip (atan2 → LS → cos/sin),
        // which is O(1/dt²)-conditioned for small angular steps (dt≤0.01) and
        // causes O(1e-10) knot error even when the frame vectors are correct.
        // After center refinement the error drops to O(dt² · ε_machine · cond) ≈ 1e-12.
        //
        // Applied unconditionally so that best_use_direct is set even when the
        // angle-parameterized probe already passes the 1e-6·spread gate.
        if (best_r1 > 1e-12 && best_r2 > 1e-12) {
            double om1_try = std::abs(best_a[0]) / best_r1;
            double om2_try = std::abs(best_a[1]) / best_r2;
            if (om1_try > 1e-6 && om2_try > 1e-6 && std::abs(om1_try - om2_try) > 1e-8) {
                // Try both (om1,om2) and (om2,om1) orderings.
                for (int swap_ab = 0; swap_ab < 2; ++swap_ab) {
                    double o1 = swap_ab ? om2_try : om1_try;
                    double o2 = swap_ab ? om1_try : om2_try;
                    double aug2[5][9];
                    for (int k = 0; k < 5; ++k) {
                        aug2[k][0] = 1.0;
                        aug2[k][1] = std::cos(o1*ts5[k]); aug2[k][2] = std::sin(o1*ts5[k]);
                        aug2[k][3] = std::cos(o2*ts5[k]); aug2[k][4] = std::sin(o2*ts5[k]);
                        for (int jj = 0; jj < 4; ++jj) aug2[k][5+jj] = pts4[k][jj];
                    }
                    bool ok2 = true;
                    for (int col = 0; col < 5 && ok2; ++col) {
                        int piv = col;
                        for (int row = col+1; row < 5; ++row)
                            if (std::abs(aug2[row][col]) > std::abs(aug2[piv][col])) piv = row;
                        if (std::abs(aug2[piv][col]) < 1e-12) { ok2 = false; break; }
                        if (piv != col)
                            for (int kk = 0; kk < 9; ++kk) std::swap(aug2[col][kk], aug2[piv][kk]);
                        double inv2 = 1.0 / aug2[col][col];
                        for (int kk = 0; kk < 9; ++kk) aug2[col][kk] *= inv2;
                        for (int row = 0; row < 5; ++row) {
                            if (row == col) continue;
                            double ff = aug2[row][col];
                            for (int kk = 0; kk < 9; ++kk) aug2[row][kk] -= ff * aug2[col][kk];
                        }
                    }
                    if (!ok2) continue;
                    VecN<4> e_gj{}, c1{}, s1{}, c2{}, s2{};
                    for (int jj = 0; jj < 4; ++jj) {
                        e_gj[jj] = aug2[0][5+jj]; c1[jj] = aug2[1][5+jj]; s1[jj] = aug2[2][5+jj];
                        c2[jj]   = aug2[3][5+jj]; s2[jj] = aug2[4][5+jj];
                    }
                    VecN<4> mean_res{};
                    for (int k = 0; k < 5; ++k) {
                        VecN<4> rec = e_gj
                            + c1*std::cos(o1*ts5[k]) + s1*std::sin(o1*ts5[k])
                            + c2*std::cos(o2*ts5[k]) + s2*std::sin(o2*ts5[k]);
                        mean_res = mean_res + (rec - pts4[k]);
                    }
                    VecN<4> e_ref = e_gj - mean_res * 0.2;
                    double max_err_d = 0.0;
                    for (int k = 0; k < 5; ++k) {
                        VecN<4> rec = e_ref
                            + c1*std::cos(o1*ts5[k]) + s1*std::sin(o1*ts5[k])
                            + c2*std::cos(o2*ts5[k]) + s2*std::sin(o2*ts5[k]);
                        double e = (rec - pts4[k]).norm();
                        if (e > max_err_d) max_err_d = e;
                    }
                    if (max_err_d < best_err) {
                        best_err = max_err_d;
                        best_center = e_ref;
                        best_Av = c1; best_Bv = s1;
                        best_Av2 = c2; best_Bv2 = s2;
                        best_o1d = o1; best_o2d = o2;
                        best_use_direct = true;
                    }
                }
            }
        }

        // ── Orthogonality grid search (non-uniform spacing fallback) ──────────
        // When ESPRIT and moment-matrix splits all fail, search (ω₁, ω₂) ∈ [0.1,5]²
        // for the Clifford orthogonality condition: span(A₁,B₁) ⊥ span(A₂,B₂).
        //
        // For fixed (ω₁,ω₂) the model p(tₖ)=e+A₁cosω₁tₖ+B₁sinω₁tₖ+… is linear.
        // Solve the 5×5 system B(ω₁,ω₂)·X = P_data.  The correct (ω₁,ω₂) satisfy
        // Score = ||[A₁,B₁]ᵀ[A₂,B₂]||² = 0 exactly.
        // Non-uniform sampling avoids aliasing: the grid can cover any ω range.
        if (best_err > 1e-6 * spread) {
            // Data matrix: raw hull-coord points (NOT sphere-centered)
            double P[5][4];
            for (int k = 0; k < 5; ++k)
                for (int j = 0; j < 4; ++j)
                    P[k][j] = pts4[k][j];

            // score_fn: 5×9 Gauss-Jordan on [B(ω1,ω2) | P], return orthogonality score.
            // Returns 1e30 if B is singular.
            auto score_fn = [&](double om1, double om2) -> double {
                double aug[5][9];
                for (int k = 0; k < 5; ++k) {
                    aug[k][0] = 1.0;
                    aug[k][1] = std::cos(om1*ts5[k]);
                    aug[k][2] = std::sin(om1*ts5[k]);
                    aug[k][3] = std::cos(om2*ts5[k]);
                    aug[k][4] = std::sin(om2*ts5[k]);
                    for (int j = 0; j < 4; ++j) aug[k][5+j] = P[k][j];
                }
                for (int col = 0; col < 5; ++col) {
                    int piv = col;
                    for (int row = col+1; row < 5; ++row)
                        if (std::abs(aug[row][col]) > std::abs(aug[piv][col])) piv = row;
                    if (std::abs(aug[piv][col]) < 1e-12) return 1e30;
                    if (piv != col)
                        for (int kk = 0; kk < 9; ++kk) std::swap(aug[col][kk], aug[piv][kk]);
                    double inv = 1.0 / aug[col][col];
                    for (int kk = 0; kk < 9; ++kk) aug[col][kk] *= inv;
                    for (int row = 0; row < 5; ++row) {
                        if (row == col) continue;
                        double f = aug[row][col];
                        for (int kk = 0; kk < 9; ++kk) aug[row][kk] -= f * aug[col][kk];
                    }
                }
                // X rows: 0=e, 1=A₁, 2=B₁, 3=A₂, 4=B₂
                double score = 0.0;
                for (int ia = 1; ia <= 2; ++ia)
                    for (int ib = 3; ib <= 4; ++ib) {
                        double dot = 0.0;
                        for (int j = 0; j < 4; ++j) dot += aug[ia][5+j] * aug[ib][5+j];
                        score += dot * dot;
                    }
                return score;
            };

            // Evaluate the full coarse grid (50×50 upper triangle, ω₂ > ω₁).
            // Then find grid-local-minima: score strictly less than all 4-directional
            // neighbours that exist in the upper triangle.  Run Newton's method from
            // each (sorted by score, top-5), probing both planes per candidate.
            // This handles multiple basins — the global grid minimum may be spurious.
            const int NOM = 50;
            const double om_step = 0.1;
            double gscore[NOM][NOM]; // [i][j]: o1=(i+1)*step, o2=(j+1)*step
            for (int i = 0; i < NOM; ++i)
                for (int j = 0; j < NOM; ++j)
                    gscore[i][j] = (j > i) ? score_fn((i+1)*om_step, (j+1)*om_step) : 1e30;

            struct Cand { double sc, o1, o2; };
            Cand locals[30]; int n_locals = 0;
            for (int i = 0; i < NOM && n_locals < 30; ++i) {
                for (int j = i+1; j < NOM && n_locals < 30; ++j) {
                    double s = gscore[i][j];
                    if (s >= 1e20) continue;
                    // Check 4-directional neighbours within the upper triangle
                    if (i   >   0 && gscore[i-1][j  ] <= s) continue;
                    if (i+1 <   j && gscore[i+1][j  ] <= s) continue;
                    if (j   > i+1 && gscore[i  ][j-1] <= s) continue;
                    if (j+1 < NOM && gscore[i  ][j+1] <= s) continue;
                    locals[n_locals++] = {s, (i+1)*om_step, (j+1)*om_step};
                }
            }
            std::sort(locals, locals + n_locals,
                      [](Cand const& a, Cand const& b){ return a.sc < b.sc; });
            int n_try = std::min(n_locals, 5);

            // Refine (o1,o2) and probe the frame it defines.
            // Uses frequency-aware angle unwrapping (known o1,o2 from Newton) so
            // large inter-sample gaps (ω·Δt > π) are handled correctly.
            auto extract_and_probe = [&](double o1, double o2) {
                double aug2[5][9];
                for (int k = 0; k < 5; ++k) {
                    aug2[k][0] = 1.0;
                    aug2[k][1] = std::cos(o1*ts5[k]); aug2[k][2] = std::sin(o1*ts5[k]);
                    aug2[k][3] = std::cos(o2*ts5[k]); aug2[k][4] = std::sin(o2*ts5[k]);
                    for (int j = 0; j < 4; ++j) aug2[k][5+j] = P[k][j];
                }
                bool ok2 = true;
                for (int col = 0; col < 5 && ok2; ++col) {
                    int piv = col;
                    for (int row = col+1; row < 5; ++row)
                        if (std::abs(aug2[row][col]) > std::abs(aug2[piv][col])) piv = row;
                    if (std::abs(aug2[piv][col]) < 1e-12) { ok2 = false; break; }
                    if (piv != col)
                        for (int kk = 0; kk < 9; ++kk) std::swap(aug2[col][kk], aug2[piv][kk]);
                    double inv2 = 1.0 / aug2[col][col];
                    for (int kk = 0; kk < 9; ++kk) aug2[col][kk] *= inv2;
                    for (int row = 0; row < 5; ++row) {
                        if (row == col) continue;
                        double f2 = aug2[row][col];
                        for (int kk = 0; kk < 9; ++kk) aug2[row][kk] -= f2 * aug2[col][kk];
                    }
                }
                if (!ok2) return;
                // Use GJ constant term as center (self-consistent with frame extraction).
                // aug2[0][5..8] = e_gj: the Clifford center in the GJ decomposition of
                // pts4.  Using e_gj avoids the mismatch between bf.center (sphere fit)
                // and the GJ constant term that causes O(1e-10) errors for ill-conditioned
                // windows (e.g. very small dt).
                VecN<4> e_gj{};
                for (int j = 0; j < 4; ++j) e_gj[j] = aug2[0][5+j];
                VecN<4> q_gj[5];
                for (int k = 0; k < 5; ++k)
                    q_gj[k] = VecN<4>{P[k][0]-e_gj[0], P[k][1]-e_gj[1],
                                      P[k][2]-e_gj[2], P[k][3]-e_gj[3]};

                // Direct formula: the GJ solution satisfies
                //   e_gj + c1*cos(o1*t[k]) + s1*sin(o1*t[k])
                //        + c2*cos(o2*t[k]) + s2*sin(o2*t[k]) = pts4[k]   (exact)
                // Using it directly bypasses the angle-extraction round-trip
                // (atan2 → LS fit → cos/sin), which is O(1/dt²)-conditioned
                // for small angular steps and causes O(1e-10) knot error.
                //
                // Center refinement: the GJ errors have a large constant-mode
                // component (delta_e ≈ GJ center error) plus small O(dt²)
                // oscillatory residuals.  Subtracting mean(rec[k]-pts4[k]) from
                // the center eliminates the dominant constant error and leaves
                // only the O(dt²) term ≈ O(1e-12), well below 1e-10.
                {
                    VecN<4> c1{}, s1{}, c2{}, s2{};
                    for (int j = 0; j < 4; ++j) {
                        c1[j] = aug2[1][5+j]; s1[j] = aug2[2][5+j];
                        c2[j] = aug2[3][5+j]; s2[j] = aug2[4][5+j];
                    }
                    // Compute residuals and their mean
                    VecN<4> mean_res{};
                    for (int k = 0; k < 5; ++k) {
                        VecN<4> rec = e_gj
                            + c1 * std::cos(o1 * ts5[k]) + s1 * std::sin(o1 * ts5[k])
                            + c2 * std::cos(o2 * ts5[k]) + s2 * std::sin(o2 * ts5[k]);
                        mean_res = mean_res + (rec - pts4[k]);
                    }
                    VecN<4> e_ref = e_gj - mean_res * 0.2;
                    // Measure error with refined center
                    double max_err_d = 0.0;
                    for (int k = 0; k < 5; ++k) {
                        VecN<4> rec = e_ref
                            + c1 * std::cos(o1 * ts5[k]) + s1 * std::sin(o1 * ts5[k])
                            + c2 * std::cos(o2 * ts5[k]) + s2 * std::sin(o2 * ts5[k]);
                        double e = (rec - pts4[k]).norm();
                        if (e > max_err_d) max_err_d = e;
                    }
                    if (max_err_d < best_err) {
                        best_err = max_err_d;
                        best_center = e_ref;
                        best_Av  = c1; best_Bv  = s1;
                        best_Av2 = c2; best_Bv2 = s2;
                        best_o1d = o1; best_o2d = o2;
                        best_use_direct = true;
                    }
                }

                // Columns 1,2 carry ω₁=o1; columns 3,4 carry ω₂=o2.
                // plane=0 → (f1e,f2e) plane has frequency o1, complement has o2.
                // plane=1 → (f1e,f2e) plane has frequency o2, complement has o1.
                for (int plane = 0; plane < 2; ++plane) {
                    VecN<4> Av{}, Bv{};
                    for (int j = 0; j < 4; ++j) {
                        Av[j] = aug2[1+2*plane][5+j];
                        Bv[j] = aug2[2+2*plane][5+j];
                    }
                    double nA = Av.norm();
                    if (nA < 1e-12 * spread) continue;
                    VecN<4> f1e = Av * (1.0/nA);
                    Bv = Bv - f1e * f1e.dot(Bv);
                    double nB = Bv.norm();
                    if (nB < 1e-12 * spread) continue;
                    VecN<4> f2e = Bv * (1.0/nB);

                    VecN<4> f3e{}, f4e{};
                    if (!detail::make_complement_4d(f1e, f2e, f3e, f4e)) continue;

                    // Radii — projected onto q_gj (centered at e_gj for consistency)
                    double r1sq = 0.0, r2sq = 0.0;
                    for (int k = 0; k < 5; ++k) {
                        double a1 = f1e.dot(q_gj[k]), a2 = f2e.dot(q_gj[k]);
                        r1sq += a1*a1 + a2*a2;
                        double a3 = f3e.dot(q_gj[k]), a4 = f4e.dot(q_gj[k]);
                        r2sq += a3*a3 + a4*a4;
                    }
                    r1sq /= 5.0; r2sq /= 5.0;
                    if (r1sq < 1e-24 || r2sq < 1e-24) continue;
                    double r1e = std::sqrt(r1sq), r2e = std::sqrt(r2sq);

                    // Frequency-aware angle unwrapping.
                    // Standard unwrapping (consecutive-diff < π) fails when
                    // ω·Δt > π (large gaps).  Using the known frequency o_p/o_q
                    // allows exact disambiguation for any Δt.
                    double om_p = (plane == 0) ? o1 : o2; // freq of (f1e,f2e) plane
                    double om_q = (plane == 0) ? o2 : o1; // freq of complement plane
                    double th1e[5], th2e[5];
                    for (int k = 0; k < 5; ++k) {
                        th1e[k] = std::atan2(f2e.dot(q_gj[k]), f1e.dot(q_gj[k]));
                        th2e[k] = std::atan2(f4e.dot(q_gj[k]), f3e.dot(q_gj[k]));
                    }
                    double phi1 = th1e[0] - om_p * ts5[0];
                    for (int k = 1; k < 5; ++k) {
                        double diff = th1e[k] - (om_p * ts5[k] + phi1);
                        th1e[k] -= std::round(diff / (2*detail::PI)) * (2*detail::PI);
                    }
                    double phi2 = th2e[0] - om_q * ts5[0];
                    for (int k = 1; k < 5; ++k) {
                        double diff = th2e[k] - (om_q * ts5[k] + phi2);
                        th2e[k] -= std::round(diff / (2*detail::PI)) * (2*detail::PI);
                    }

                    // LS fit: uv(t) = a_e·t + b_e
                    VecN<2> upts[5];
                    for (int k = 0; k < 5; ++k)
                        upts[k] = VecN<2>{r1e * th1e[k], r2e * th2e[k]};
                    VecN<2> sum_u{}, sum_tu{};
                    for (int k = 0; k < 5; ++k) {
                        sum_u  = sum_u  + upts[k];
                        sum_tu = sum_tu + upts[k] * ts5[k];
                    }
                    VecN<2> a_e = (sum_tu * 5.0 - sum_u * sum_t) * (1.0 / ls_det);
                    VecN<2> b_e = (sum_u * sum_t2 - sum_tu * sum_t) * (1.0 / ls_det);

                    // Roundtrip error — reconstruct using e_gj (consistent with frames)
                    double max_err_e = 0.0;
                    for (int k = 0; k < 5; ++k) {
                        VecN<2> uv = a_e * ts5[k] + b_e;
                        double ph1 = uv[0] / r1e, ph2 = uv[1] / r2e;
                        VecN<4> rec = e_gj
                                    + f1e * (r1e * std::cos(ph1))
                                    + f2e * (r1e * std::sin(ph1))
                                    + f3e * (r2e * std::cos(ph2))
                                    + f4e * (r2e * std::sin(ph2));
                        double e_err = (rec - pts4[k]).norm();
                        if (e_err > max_err_e) max_err_e = e_err;
                    }
                    if (max_err_e < best_err) {
                        best_err = max_err_e;
                        best_f1 = f1e; best_f2 = f2e;
                        best_f3 = f3e; best_f4 = f4e;
                        best_r1 = r1e; best_r2 = r2e;
                        best_a  = a_e; best_b  = b_e;
                        best_center = e_gj; // GJ center for this frame
                    }
                }
            };

            // Newton's method (numerical 2×2 Hessian) from each candidate.
            // Quadratic convergence drives Score to machine epsilon in ~20 iterations.
            // hfd adapts to current residual: h ≈ 0.01·‖Δω‖ keeps FD accurate even
            // when Δω → 0, preventing the fixed-h gradient cancellation that stalls
            // convergence once score ≲ 1e-12.
            for (int ci = 0; ci < n_try; ++ci) {
                double o1 = locals[ci].o1, o2 = locals[ci].o2;
                for (int iter = 0; iter < 60; ++iter) {
                    double s0  = score_fn(o1,      o2);
                    const double hfd = std::max(1e-10, std::sqrt(s0) * 0.01);
                    if (s0 < 1e-28) break;
                    double sp1 = score_fn(o1+hfd,  o2);
                    double sm1 = score_fn(o1-hfd,  o2);
                    double sp2 = score_fn(o1,      o2+hfd);
                    double sm2 = score_fn(o1,      o2-hfd);
                    if (sp1 >= 1e20 || sm1 >= 1e20 || sp2 >= 1e20 || sm2 >= 1e20) break;
                    double g1  = (sp1 - sm1) / (2*hfd);
                    double g2  = (sp2 - sm2) / (2*hfd);
                    double H11 = (sp1 - 2*s0 + sm1) / (hfd*hfd);
                    double H22 = (sp2 - 2*s0 + sm2) / (hfd*hfd);
                    double spp = score_fn(o1+hfd, o2+hfd);
                    double spm = score_fn(o1+hfd, o2-hfd);
                    double smp = score_fn(o1-hfd, o2+hfd);
                    double smm = score_fn(o1-hfd, o2-hfd);
                    if (spp >= 1e20 || spm >= 1e20 || smp >= 1e20 || smm >= 1e20) break;
                    double H12 = (spp - spm - smp + smm) / (4*hfd*hfd);
                    double det = H11*H22 - H12*H12;
                    double d1, d2;
                    if (std::abs(det) > 1e-20 * (H11*H11 + H22*H22 + 1e-60)) {
                        d1 = -(H22*g1 - H12*g2) / det; // Newton step
                        d2 = -(-H12*g1 + H11*g2) / det;
                    } else {
                        double gnorm = std::sqrt(g1*g1 + g2*g2);
                        if (gnorm < 1e-30) break;
                        d1 = -0.01*g1/gnorm; d2 = -0.01*g2/gnorm; // gradient fallback
                    }
                    double alpha = 1.0;
                    bool improved = false;
                    for (int ls = 0; ls < 30; ++ls) {
                        double n1 = o1 + alpha*d1, n2 = o2 + alpha*d2;
                        if (n1 > 0 && n2 > 0 && score_fn(n1, n2) < s0) {
                            o1 = n1; o2 = n2; improved = true; break;
                        }
                        alpha *= 0.5;
                    }
                    if (!improved) break;
                }
                extract_and_probe(o1, o2);
                if (best_err <= 1e-6 * spread) break;
            }
        }

        if (best_err > 1e-6 * spread) return false;

        // Accept
        e1_ = best_f1; e2_ = best_f2;
        e3_ = best_f3; e4_ = best_f4;
        center4_ = best_center;
        r1_ = best_r1; r2_ = best_r2;
        inner_a_lin_ = best_a;
        inner_b_lin_ = best_b;
        if (best_use_direct) {
            use_direct_ = true;
            Av_ = best_Av; Bv_ = best_Bv;
            Av2_ = best_Av2; Bv2_ = best_Bv2;
            o1d_ = best_o1d; o2d_ = best_o2d;

            // Ambient GJ: solve same interpolation system in ℝᴺ using original pts.
            // Bypasses lift() which accumulates O(1e-10) error when basis_[3] is
            // ill-conditioned (happens for small dt where the 4th span direction
            // is a dt³-order Taylor residual).
            double aug_amb[5][5 + N];
            for (int k = 0; k < 5; ++k) {
                aug_amb[k][0] = 1.0;
                aug_amb[k][1] = std::cos(best_o1d * ts5[k]);
                aug_amb[k][2] = std::sin(best_o1d * ts5[k]);
                aug_amb[k][3] = std::cos(best_o2d * ts5[k]);
                aug_amb[k][4] = std::sin(best_o2d * ts5[k]);
                for (int jj = 0; jj < N; ++jj) aug_amb[k][5 + jj] = pts_amb[k][jj];
            }
            bool ok_amb = true;
            for (int col = 0; col < 5 && ok_amb; ++col) {
                int piv = col;
                for (int row = col+1; row < 5; ++row)
                    if (std::abs(aug_amb[row][col]) > std::abs(aug_amb[piv][col])) piv = row;
                if (std::abs(aug_amb[piv][col]) < 1e-12) { ok_amb = false; break; }
                if (piv != col)
                    for (int kk = 0; kk < 5+N; ++kk) std::swap(aug_amb[col][kk], aug_amb[piv][kk]);
                double inv_a = 1.0 / aug_amb[col][col];
                for (int kk = 0; kk < 5+N; ++kk) aug_amb[col][kk] *= inv_a;
                for (int row = 0; row < 5; ++row) {
                    if (row == col) continue;
                    double ff = aug_amb[row][col];
                    for (int kk = 0; kk < 5+N; ++kk) aug_amb[row][kk] -= ff * aug_amb[col][kk];
                }
            }
            if (ok_amb) {
                for (int jj = 0; jj < N; ++jj) {
                    center_amb_[jj] = aug_amb[0][5+jj];
                    Av_amb_[jj]     = aug_amb[1][5+jj];
                    Bv_amb_[jj]     = aug_amb[2][5+jj];
                    Av2_amb_[jj]    = aug_amb[3][5+jj];
                    Bv2_amb_[jj]    = aug_amb[4][5+jj];
                }
                use_direct_amb_ = true;
            }
        }
        return true;
    }

public:
    CliffordWindow() = default;

    CliffordWindow(
        VecN<N> const& p0, VecN<N> const& p1, VecN<N> const& p2,
        VecN<N> const& p3, VecN<N> const& p4,
        double t0, double t1, double t2, double t3, double t4)
    {
        VecN<N> pts5[5] = {p0, p1, p2, p3, p4};

        // Compute affine hull
        auto hull = detail::affine_hull_4d<N>(pts5);
        origin_    = hull.origin;
        for (int j = 0; j < 4; ++j) basis_[j] = hull.basis[j];
        hull_rank_ = hull.rank;

        // ── Level 1: Clifford torus in 4D hull ───────────────────────────
        if (hull_rank_ >= 4) {
            VecN<4> p4[5];
            for (int k = 0; k < 5; ++k)
                p4[k] = VecN<4>{hull.coords[k][0], hull.coords[k][1],
                                hull.coords[k][2], hull.coords[k][3]};
            auto bf = detail::clifford_best_fit_4d(p4);
            if (try_clifford_(p4, pts5, bf, t0,t1,t2,t3,t4)) {
                valid_ = true;
                return;
            }
        }

        // ── Level 2: CylinderWindow<3> in 3D sub-hull ────────────────────
        // The 3D hull coordinates are exact (orthonormal basis) so Euclidean
        // geometry is preserved: CylinderWindow<3> operates on true distances.
        if (hull_rank_ >= 3) {
            VecN<3> p3[5];
            for (int k = 0; k < 5; ++k)
                p3[k] = VecN<3>{hull.coords[k][0], hull.coords[k][1],
                                hull.coords[k][2]};
            CylinderWindow<3> cw(p3[0],p3[1],p3[2],p3[3],p3[4],
                                 t0,t1,t2,t3,t4);
            if (cw.valid()) {
                cyl_win3d_   = std::move(cw);
                use_cylinder_ = true;
                valid_        = true;
                return;
            }
        }

        // ── Level 3: ConicWindow<2> in 2D sub-hull ───────────────────────
        if (hull_rank_ >= 2) {
            VecN<2> p2[5];
            for (int k = 0; k < 5; ++k)
                p2[k] = VecN<2>{hull.coords[k][0], hull.coords[k][1]};
            ConicWindow<2> pw(p2[0],p2[1],p2[2],p2[3],p2[4],
                              t0,t1,t2,t3,t4);
            if (pw.valid()) {
                plane_win_ = std::move(pw);
                use_planar_ = true;
                valid_      = true;
                return;
            }
        }

        // valid_ stays false → blend_curve uses LagrangeWindow<N>
    }

    bool valid()         const { return valid_; }
    bool used_clifford() const { return valid_ && !use_cylinder_ && !use_planar_; }
    bool used_cylinder() const { return use_cylinder_; }
    bool used_planar()   const { return use_planar_; }
    int  hull_rank()     const { return hull_rank_; }

    VecN<N> operator()(double t) const
    {
        if (use_planar_) {
            VecN<2> uv = plane_win_(t);
            return lift(uv[0], uv[1], 0.0, 0.0);
        }
        if (use_cylinder_) {
            VecN<3> p3 = cyl_win3d_(t);
            return lift(p3[0], p3[1], p3[2], 0.0);
        }
        // Level 1: Clifford torus
        // Ambient path: GJ solved directly in ℝᴺ — no lift() needed.
        if (use_direct_amb_) {
            return center_amb_
                + Av_amb_  * std::cos(o1d_ * t) + Bv_amb_  * std::sin(o1d_ * t)
                + Av2_amb_ * std::cos(o2d_ * t) + Bv2_amb_ * std::sin(o2d_ * t);
        }
        VecN<4> p4;
        if (use_direct_) {
            // Direct GJ formula in hull coords (fallback when ambient solve fails).
            p4 = center4_
               + Av_  * std::cos(o1d_ * t) + Bv_  * std::sin(o1d_ * t)
               + Av2_ * std::cos(o2d_ * t) + Bv2_ * std::sin(o2d_ * t);
        } else {
            VecN<2> uv  = inner_a_lin_ * t + inner_b_lin_;
            double  ph1 = uv[0] / r1_;
            double  ph2 = uv[1] / r2_;
            p4 = center4_
               + e1_ * (r1_ * std::cos(ph1))
               + e2_ * (r1_ * std::sin(ph1))
               + e3_ * (r2_ * std::cos(ph2))
               + e4_ * (r2_ * std::sin(ph2));
        }
        return lift(p4[0], p4[1], p4[2], p4[3]);
    }
};


// ── blend_curve(..., nd_tag{}) ────────────────────────────────────────────
//
// Dispatches by dimension:
//   N == 2 → conic_tag    (ConicWindow<2>)
//   N == 3 → cylinder_tag (CylinderWindow<3>)
//   N ≥  4 → CliffordWindow<N>
//
// used_clifford: for N≥4, set to true if at least one window used the
//   Clifford-torus path (not CylinderWindow or ConicWindow fallback).
//   For N=2/3, set the same way as the underlying used_conic/used_cylinder.

template <int N>
inline BlendResultND<N> blend_curve(
#if __cplusplus >= 202002L
    std::span<VecN<N> const> ctrl,
    std::span<double const>  times,
#else
    std::vector<VecN<N>> const& ctrl,
    std::vector<double>  const& times,
#endif
    nd_tag,
    int   pts_per_seg   = 60,
    int   smooth_N      = 2,
    bool* used_clifford = nullptr)
{
    if constexpr (N == 2) {
        return blend_curve<2>(ctrl, times, conic_tag{},
                              pts_per_seg, smooth_N, used_clifford);
    } else if constexpr (N == 3) {
        return blend_curve<3>(ctrl, times, cylinder_tag{},
                              pts_per_seg, smooth_N, used_clifford);
    } else {
        static_assert(N >= 4);

        int n = static_cast<int>(ctrl.size());
        if (n != static_cast<int>(times.size()))
            FC_THROW("fc::blend_curve(nd): ctrl and times must have equal length");
        if (n < 6)
            FC_THROW("fc::blend_curve(nd): need at least 6 control points");
        if (pts_per_seg < 2)
            FC_THROW("fc::blend_curve(nd): pts_per_seg must be >= 2");

        using Win = std::variant<CliffordWindow<N>, LagrangeWindow<N>>;
        std::vector<Win> wins;
        wins.reserve(n - 4);
        int n_clifford = 0;

        for (int i = 0; i < n-4; ++i) {
            CliffordWindow<N> cw(ctrl[i],ctrl[i+1],ctrl[i+2],ctrl[i+3],ctrl[i+4],
                                 times[i],times[i+1],times[i+2],times[i+3],times[i+4]);
            if (cw.valid()) {
                if (cw.used_clifford()) ++n_clifford;
                wins.emplace_back(std::move(cw));
            } else {
                wins.emplace_back(LagrangeWindow<N>(
                    ctrl[i],ctrl[i+1],ctrl[i+2],ctrl[i+3],ctrl[i+4],
                    times[i],times[i+1],times[i+2],times[i+3],times[i+4]));
            }
        }
        if (used_clifford) *used_clifford = (n_clifford > 0);

        auto eval_win = [](Win const& w, double t) -> VecN<N> {
            return std::visit([t](auto const& ww) -> VecN<N> { return ww(t); }, w);
        };

        int n_segs = n - 5;
        BlendResultND<N> out;
        out.pts.reserve(n_segs * pts_per_seg + 1);
        out.times.reserve(n_segs * pts_per_seg + 1);

        for (int j = 2; j <= n-4; ++j) {
            Win const& wA = wins[j-2];
            Win const& wB = wins[j-1];
            double t_lo = times[j], t_hi = times[j+1];
            int k_max = (j < n-4) ? pts_per_seg - 1 : pts_per_seg;
            for (int k = 0; k <= k_max; ++k) {
                double s = static_cast<double>(k) / pts_per_seg;
                double t = t_lo + s * (t_hi - t_lo);
                double w = smoothstep(s, smooth_N);
                VecN<N> A = eval_win(wA, t), B = eval_win(wB, t);
                out.pts.push_back(A*(1.0-w) + B*w);
                out.times.push_back(t);
            }
        }
        return out;
    }
}

} // namespace fc
