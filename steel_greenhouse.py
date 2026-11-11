#!/usr/bin/env python3
"""The steel greenhouse, one line at a time.

Reproduces the three runs of J. E. Postma's "Steel Greenhouse" simulation with
his own parameters (sphere r = 1.00 m, P = 2953 W = 235 W/m^2, hollow steel
shell 1.50-1.52 m, black surfaces, vacuum, space at 0 K, shell starting at 0 K)
and a radial heat-conduction model of the shell wall:

  original      sphere held at 253.7 K; shell inner face absorbs sigma*T_sp^4*(r_sp/r_in)^2
                and does not emit inward                      (his "original model")
  backradiation same shell boundary, but the sphere also receives A_sp*sigma*T_in^4
                from the shell, which the shell never loses   (his "back-radiation model")
  physics       the shell's inner face emits as well, so what the sphere receives
                the shell loses: net inner flux A_sp*sigma*(T_sp^4 - T_in^4)
                                                              (his "energy-conserving variant")

The only difference between "backradiation" and "physics" is one term in the
shell's inner boundary condition. Energy ledger for every run:
  created = energy to space + energy stored - energy produced by the core

Usage:
  python3 steel_greenhouse.py              # runs, prints the numbers, writes figures to results/
  python3 steel_greenhouse.py --hours 500  # longer runs
"""
import argparse
import os

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SIGMA = 5.670374419e-8                     # W m^-2 K^-4

# ---- Postma's parameters (from the captions of his figures) -----------------
R_SP = 1.00                                # sphere radius, m
P_INT = 2953.0                             # core power, W (235 W/m^2 x 4 pi m^2)
C_SP = 0.95e6                              # sphere heat capacity, J/K
R_IN, R_OUT = 1.50, 1.52                   # shell radii, m
RHO, C_STEEL, K_STEEL = 7850.0, 490.0, 45.0
T_SP0 = 253.7                              # K, sphere alone: (235 / sigma)^(1/4)
N_CELLS = 10                               # finite-volume cells across the steel wall

A_SP = 4 * np.pi * R_SP**2
A_OUT = 4 * np.pi * R_OUT**2

# finite-volume grid across the wall
r_faces = np.linspace(R_IN, R_OUT, N_CELLS + 1)
r_mid = 0.5 * (r_faces[1:] + r_faces[:-1])
vol = 4 / 3 * np.pi * (r_faces[1:] ** 3 - r_faces[:-1] ** 3)
cap = RHO * C_STEEL * vol                  # J/K per cell
C_SH = cap.sum()
A_int = 4 * np.pi * r_faces[1:-1] ** 2     # interior faces
G_int = K_STEEL * A_int / np.diff(r_mid)   # conductance between neighbouring cells, W/K


def inner_heat(mode, T_sp, T_in):
    """Heat (W) entering the shell through its inner face."""
    if mode == "physics":
        return A_SP * SIGMA * (T_sp**4 - T_in**4)      # absorbs the sphere, emits inward too
    return A_SP * SIGMA * T_sp**4                      # Postma: absorbs, never emits inward


def back_radiation(mode, T_in):
    """Radiation (W) the sphere receives from the shell."""
    return 0.0 if mode == "original" else A_SP * SIGMA * T_in**4


def rhs(t, y, mode):
    T_sp, T = y[0], y[1:]
    if mode == "original":
        dT_sp = 0.0                                    # held constant
    else:
        dT_sp = (P_INT + back_radiation(mode, T[0]) - A_SP * SIGMA * T_sp**4) / C_SP
    q = np.zeros(N_CELLS)
    flow = G_int * (T[:-1] - T[1:])                    # conduction outward between cells
    q[:-1] -= flow
    q[1:] += flow
    q[0] += inner_heat(mode, T_sp, T[0])
    q[-1] -= A_OUT * SIGMA * T[-1] ** 4                # outer face to space at 0 K
    return np.concatenate(([dT_sp], q / cap))


def run(mode, hours):
    t_end = hours * 3600.0
    t_eval = np.linspace(0, t_end, 2000)
    y0 = np.concatenate(([T_SP0], np.zeros(N_CELLS)))
    sol = solve_ivp(rhs, (0, t_end), y0, args=(mode,), method="BDF",
                    t_eval=t_eval, rtol=1e-8, atol=1e-6)
    t, T_sp, T = sol.t, sol.y[0], sol.y[1:]
    p_out = A_OUT * SIGMA * T[-1] ** 4
    p_back = np.array([back_radiation(mode, x) for x in T[0]])
    produced = P_INT * t if mode != "original" else A_SP * SIGMA * T_SP0**4 * t
    to_space = np.concatenate(([0.0], np.cumsum(0.5 * (p_out[1:] + p_out[:-1]) * np.diff(t))))
    stored = C_SP * (T_sp - T_sp[0]) + (cap[:, None] * (T - T[:, :1])).sum(axis=0)
    created = to_space + stored - produced
    return dict(t=t / 3600, T_sp=T_sp, T_in=T[0], T_out=T[-1], p_out=p_out, p_back=p_back,
                created=created)


def steady_state(mode):
    """Long-run limit, solved directly (thin wall: conduction drop included)."""
    from scipy.optimize import brentq
    if mode == "original":
        Ts = T_SP0 * np.sqrt(R_SP / R_OUT)
        return T_SP0, Ts
    if mode == "physics":
        # outer face emits P; spherical conduction drop across the wall; sphere balance
        T_out = (P_INT / (A_OUT * SIGMA)) ** 0.25
        dT = P_INT / (4 * np.pi * K_STEEL) * (1 / R_IN - 1 / R_OUT)
        T_in = T_out + dT
        T_sp = ((P_INT / A_SP + SIGMA * T_in**4) / SIGMA) ** 0.25
        return T_sp, T_out
    # backradiation: shell passes A_sp*sigma*T_sp^4 outward; sphere: sigma T_sp^4 = P/A_sp + sigma T_in^4
    def f(T_sp):
        P_sh = A_SP * SIGMA * T_sp**4
        T_out = (P_sh / (A_OUT * SIGMA)) ** 0.25
        T_in = T_out + P_sh / (4 * np.pi * K_STEEL) * (1 / R_IN - 1 / R_OUT)
        return SIGMA * T_sp**4 - P_INT / A_SP - SIGMA * T_in**4
    T_sp = brentq(f, 200, 600)
    return T_sp, (A_SP * SIGMA * T_sp**4 / (A_OUT * SIGMA)) ** 0.25


def at(r, hours):
    i = int(np.argmin(np.abs(r["t"] - hours)))
    return {k: (v[i] if hasattr(v, "__len__") else v) for k, v in r.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hours", type=float, default=73.1, help="run length in hours (default: his 73.1 h)")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "results"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    print(f"Sphere: r = {R_SP} m, P = {P_INT:.0f} W ({P_INT / A_SP:.1f} W/m^2), C = {C_SP / 1e6:.2f} MJ/K")
    print(f"Shell:  r = {R_IN}-{R_OUT} m steel, C = {C_SH / 1e6:.2f} MJ/K, {N_CELLS} cells across the wall\n")

    runs = {m: run(m, a.hours) for m in ("original", "backradiation", "physics")}
    lines = []
    for m, r in runs.items():
        e = at(r, a.hours)
        Tsp_inf, Tsh_inf = steady_state(m)
        lines.append(f"{m:14s} t = {a.hours:.1f} h: T_sp = {e['T_sp']:.1f} K, T_sh = {e['T_out']:.1f} K, "
                     f"P_out = {e['p_out']:.0f} W, P_back = {e['p_back']:.0f} W, "
                     f"created = {e['created'] / 1e6:+.1f} MJ")
        lines.append(f"{'':14s} long-run limit: T_sp = {Tsp_inf:.1f} K, T_sh = {Tsh_inf:.1f} K")
    print("\n".join(lines))
    with open(os.path.join(a.out, "results.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")

    # ---- figure 1: same code, one term -------------------------------------
    BG, FG, DIM = "#0f1419", "#e7e9ea", "#8b98a5"
    plt.rcParams.update({"text.color": FG, "axes.labelcolor": FG, "xtick.color": DIM, "ytick.color": DIM,
                         "axes.edgecolor": DIM, "figure.facecolor": BG, "axes.facecolor": BG,
                         "font.family": "DejaVu Sans"})
    B, Ph = runs["backradiation"], runs["physics"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8.2), dpi=160, sharex=True)
    ax1.axhline(T_SP0, color=DIM, ls=":", lw=1)
    ax1.text(1, T_SP0 + 3, "sphere alone 253.7 K", color=DIM, fontsize=9)
    ax1.plot(B["t"], B["T_sp"], color="#ff6b6b", lw=2.2, label="sphere, Postma's back-radiation model")
    ax1.plot(B["t"], B["T_out"], color="#ff6b6b", lw=1.4, ls="--", label="shell, Postma's back-radiation model")
    ax1.plot(Ph["t"], Ph["T_sp"], color="#4fd1c5", lw=2.2, label="sphere, shell emits from both faces")
    ax1.plot(Ph["t"], Ph["T_out"], color="#4fd1c5", lw=1.4, ls="--", label="shell, shell emits from both faces")
    ax1.set_ylabel("Temperature (K)")
    ax1.legend(loc="lower right", frameon=False, fontsize=9)
    ax1.grid(alpha=0.15)
    ax1.set_title("The steel greenhouse: same code, same parameters, one term in the shell's inner boundary",
                  loc="left", fontsize=12.5, pad=10)

    ax2.axhline(P_INT, color=DIM, ls=":", lw=1)
    ax2.text(1, P_INT + 60, "core power 2953 W", color=DIM, fontsize=9)
    ax2.plot(B["t"], B["p_out"], color="#ff6b6b", lw=2, label="power to space, Postma's model")
    ax2.plot(B["t"], B["p_back"], color="#f6ad55", lw=1.6, ls="--", label="back-radiation the shell never pays for")
    ax2.plot(Ph["t"], Ph["p_out"], color="#4fd1c5", lw=2, label="power to space, shell emits from both faces")
    ax2.set_ylabel("Power (W)")
    ax2.set_xlabel("time (hours)")
    ax2.legend(loc="upper left", frameon=False, fontsize=9)
    ax2.grid(alpha=0.15)
    eB, eP = at(B, a.hours), at(Ph, a.hours)
    sB, sP = steady_state("backradiation"), steady_state("physics")
    ax2.text(0.99, 0.04,
             f"Postma's model: energy created {eB['created'] / 1e6:+.0f} MJ, long run {sB[0]:.1f} K, "
             f"{A_SP * SIGMA * sB[0] ** 4:.0f} W to space\n"
             f"Both faces emit: energy created {eP['created'] / 1e6:+.1f} MJ, long run {sP[0]:.1f} K, "
             f"{P_INT:.0f} W to space",
             transform=ax2.transAxes, ha="right", va="bottom", fontsize=9.5,
             bbox=dict(boxstyle="round,pad=0.4", fc="#1c232b", ec=DIM, lw=0.6))
    fig.text(0.01, 0.012, "Sphere r = 1.00 m, P = 2953 W, C = 0.95 MJ/K · steel shell 1.50-1.52 m, k = 45 W/(m K), "
             "C = 2.20 MJ/K · vacuum, black surfaces, space at 0 K   ·   Made with Python by @splugged",
             fontsize=8, color=DIM)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(os.path.join(a.out, "one_term.png"), facecolor=BG)

    # ---- figure 2: any finite shell warms the sphere ----------------------
    k = np.linspace(1.0, 6.0, 200)                    # shell radius / sphere radius
    T_phys = ((P_INT / A_SP) * (1 + 1 / k**2) / SIGMA) ** 0.25
    fig, ax = plt.subplots(figsize=(10, 4.6), dpi=160)
    ax.axhline(T_SP0, color=DIM, ls=":", lw=1)
    ax.text(5.95, T_SP0 - 1.0, "no shell: 253.7 K", color=DIM, fontsize=9, ha="right", va="top")
    ax.plot(k, T_phys, color="#4fd1c5", lw=2.4, label="steady state, black thin shell, both faces emit")
    ax.plot([R_OUT / R_SP], [((P_INT / A_SP) * (1 + (R_SP / R_OUT) ** 2) / SIGMA) ** 0.25], "o", color=FG)
    ax.annotate("Postma's geometry (1.52 m)", (R_OUT, 277.6), xytext=(2.0, 290), color=FG, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=DIM))
    ax.set_ylim(T_SP0 - 6, 305)
    ax.set_xlabel("shell radius / sphere radius")
    ax.set_ylabel("sphere temperature (K)")
    ax.set_title("Any finite shell warms the sphere:  sigma T^4 = 235 x (1 + (r_sp/r_sh)^2)", loc="left", fontsize=12.5)
    ax.grid(alpha=0.15)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    fig.text(0.01, 0.01, "Made with Python by @splugged", fontsize=8, color=DIM)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(os.path.join(a.out, "any_radius.png"), facecolor=BG)
    print(f"\nfigures written to {a.out}")


if __name__ == "__main__":
    main()
