#!/usr/bin/env python3
"""Figure: the Pisot boundary for a_p = 1/(2^(1/p) - 1).  Output: fig/p343_boundary.png, fig/p343_boundary.pdf"""
import os

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42  # TrueType, not Type 3 (arXiv flags Type 3 fonts)
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

os.makedirs("fig", exist_ok=True)
P = np.arange(2, 13)


def conj(p):
    t = 2 ** (1 / p)
    return np.array([1 / (t * np.exp(2j * np.pi * s / p) - 1) for s in range(1, p)])


mx = np.array([np.abs(conj(p)).max() for p in P])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1.15, 1]})
ax1.axhspan(0, 1, color="#e8f1fb")
ax1.axhline(1, color="k", lw=1)
ax1.plot(P, mx, color="#555", lw=1, zorder=1)
ax1.scatter(P, mx, c=["#1f77b4" if v < 1 else "#d62728" for v in mx], s=46, zorder=2)
ax1.annotate("p = 5: Ramanujan (7.4.1)\na alone is Pisot", xy=(5, mx[3]), xytext=(2.3, 1.55),
             arrowprops=dict(arrowstyle="->", color="#1f77b4"), color="#1f77b4", fontsize=9)
ax1.annotate("p = 7: Ramanujan (7.4.2)\nneeds extra units b, c", xy=(7, mx[5]), xytext=(7.6, 0.45),
             arrowprops=dict(arrowstyle="->", color="#d62728"), color="#d62728", fontsize=9)
ax1.text(10.2, 0.08, "Pisot region", color="#1f77b4", fontsize=9)
ax1.set_xticks(P)
ax1.set_xlabel("p")
ax1.set_ylabel(r"$\max_{s\neq 0}\,|\sigma_s(a_p)|$")
ax1.set_title(r"$a_p = 1/(2^{1/p}-1)$ is Pisot $\Leftrightarrow 2^{1/p} > 2\cos(2\pi/p) \Leftrightarrow p \leq 6$", fontsize=10)
ax1.set_ylim(0, 2.0)

u = np.linspace(0, 2 * np.pi, 400)
ax2.plot(np.cos(u), np.sin(u), color="k", lw=1)
c5, c7 = conj(5), conj(7)
ax2.scatter(c5.real, c5.imag, s=46, color="#1f77b4", label=r"$p=5$: all $|\sigma_s| < 1$", zorder=3)
ax2.scatter(c7.real, c7.imag, s=46, marker="s", facecolors="none", edgecolors="#d62728",
            label=r"$p=7$: one pair $|\sigma_s| = 1.0896$", zorder=3)
ax2.axhline(0, color="#bbb", lw=0.6)
ax2.axvline(0, color="#bbb", lw=0.6)
ax2.set_aspect("equal")
ax2.set_xlim(-1.35, 1.35)
ax2.set_ylim(-1.35, 1.35)
ax2.set_title("non-real conjugates vs. the unit circle", fontsize=10)
ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=1, fontsize=8.5, frameon=False)

fig.suptitle("Page 343 of Ramanujan's Lost Notebook: why 5 and 7", fontsize=12)
fig.tight_layout()
fig.savefig("fig/p343_boundary.png", dpi=200)
fig.savefig("fig/p343_boundary.pdf")
print("max |sigma_s(a_p)|, p = 2..12:", np.round(mx, 7).tolist())
print("p=7 moduli:", np.round(np.sort(np.abs(c7)), 4).tolist())
