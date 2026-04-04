"""
demo_nd_clifford.py — Visualise 4D torus-helix blending with nd_tag.

Compiles and runs demo_nd_plot.cpp, then plots three panels:
  Left:   x-y projection (first S¹ factor)
  Centre: z-w projection (second S¹ factor)
  Right:  max knot error vs dt  (nd_tag vs conic_tag, log-log)

Usage:
    python3 demo_nd_clifford.py
"""

import subprocess, sys, os, io
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ── 1. Compile C++ data source ────────────────────────────────────────────────
HERE = os.path.dirname(os.path.abspath(__file__))
SRC  = os.path.join(HERE, "demo_nd_plot.cpp")
EXE  = os.path.join(HERE, "demo_nd_plot")

ret = subprocess.run(
    ["g++", "-std=c++17", "-O2", f"-I{HERE}",
     "-isystem", "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1",
     "-o", EXE, SRC],
    capture_output=True, text=True)
if ret.returncode != 0:
    print("Compile error:\n", ret.stderr); sys.exit(1)

raw = subprocess.run([EXE], capture_output=True, text=True).stdout

# ── 2. Parse sections ─────────────────────────────────────────────────────────
sections = {}
current  = None
for line in raw.splitlines():
    if line in ("CTRL", "ND", "CONIC", "TRUE", "ERRORS"):
        current = line; sections[current] = []; continue
    if current and "," in line and not line.startswith("t,") and not line.startswith("dt,"):
        sections[current].append(line)

def to_array(key, cols=5):
    return np.array([[float(v) for v in r.split(",")] for r in sections[key]])

ctrl  = to_array("CTRL")          # t x y z w
nd    = to_array("ND")
conic = to_array("CONIC")
true_ = to_array("TRUE")
errs  = to_array("ERRORS", 3)     # dt  nd_err  conic_err

# Replace 0.0 conic errors (window failed → not interpolating) with NaN
conic_err = errs[:, 2].copy()
conic_err[conic_err == 0] = np.nan

# ── 3. Figure layout ──────────────────────────────────────────────────────────
fig = plt.figure(figsize=(15, 5))
fig.suptitle(
    r"4D torus helix  $p(t)=(2\cos t,\;2\sin t,\;\cos\sqrt{2}\,t,\;\sin\sqrt{2}\,t)$"
    "\nblended with nd_tag (CliffordWindow) vs conic_tag (ConicWindow)",
    fontsize=11)

gs   = gridspec.GridSpec(1, 3, figure=fig, wspace=0.35, left=0.06, right=0.97)
ax_xy = fig.add_subplot(gs[0])
ax_zw = fig.add_subplot(gs[1])
ax_err = fig.add_subplot(gs[2])

# Projection helpers
def proj_xy(arr): return arr[:, 1], arr[:, 2]
def proj_zw(arr): return arr[:, 3], arr[:, 4]

CTRL_RANGE = slice(2, -3)  # interior control points (exact interpolation)

for ax, proj, xlabel, ylabel, title in [
    (ax_xy, proj_xy, "x", "y", r"$xy$-projection  (first $S^1$ factor)"),
    (ax_zw, proj_zw, "z", "w", r"$zw$-projection  (second $S^1$ factor)"),
]:
    # True torus curve (faint grey)
    ax.plot(*proj(true_), color="0.80", lw=1.2, zorder=1, label="true torus")
    # nd_tag blend
    ax.plot(*proj(nd),    color="#1f77b4", lw=1.5, zorder=3, label="nd_tag")
    # conic_tag blend
    ax.plot(*proj(conic), color="#ff7f0e", lw=1.5, ls="--", zorder=2, label="conic_tag")
    # Control points
    cx, cy = proj(ctrl)
    ax.scatter(cx[CTRL_RANGE], cy[CTRL_RANGE],
               s=28, color="#2ca02c", zorder=5, label="control pts")
    # Mark first and last interior knot
    ax.scatter([cx[2], cx[-3]], [cy[2], cy[-3]],
               s=55, color="#d62728", marker="*", zorder=6)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.set_aspect("equal")
    ax.legend(fontsize=8, loc="lower right")

# ── Error vs dt panel ─────────────────────────────────────────────────────────
dt_v  = errs[:, 0]
nd_e  = errs[:, 1]

ax_err.loglog(dt_v, nd_e,    "o-",  color="#1f77b4", lw=1.8, ms=5,
              label="nd_tag (CliffordWindow)")
ax_err.loglog(dt_v, conic_err, "s--", color="#ff7f0e", lw=1.8, ms=5,
              label="conic_tag (ConicWindow)\n[NaN = window failed]")

# Reference slopes
dt_ref = np.array([dt_v[0], dt_v[-1]])
ax_err.loglog(dt_ref, 1e-15 * np.ones(2), ":", color="0.6", lw=1)
ax_err.text(dt_ref[0]*1.1, 3e-15, "machine ε", fontsize=7, color="0.5")

ax_err.set_xlabel("dt  (control-point spacing)")
ax_err.set_ylabel("max knot error")
ax_err.set_title("Knot interpolation error vs dt", fontsize=10)
ax_err.legend(fontsize=8)
ax_err.grid(True, which="both", ls=":", alpha=0.4)
ax_err.set_ylim(1e-16, 10)

# Annotation box
ax_err.text(0.03, 0.08,
    "nd_tag: ≤ 5e-14 (machine ε)\nconic_tag: O(dt²) error",
    transform=ax_err.transAxes, fontsize=8,
    bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="0.7", alpha=0.9))

out = os.path.join(HERE, "demo_nd_clifford.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved {out}")
plt.show()
