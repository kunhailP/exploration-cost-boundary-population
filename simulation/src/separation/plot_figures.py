"""Draw the manuscript figures from simulation/results/separation/.

  python plot_figures.py

Main text
  fig1_mechanism.pdf           Theorem 3 construction: exact functions, no data
  fig2_separation_map.pdf      (a) exact HT risk, (b) first-order RMSE approximation of the boundary
                               estimator, with exploration-loss contours; from phase.csv
  fig3_exploration_costs.pdf   closed-form design costs at V_sp = 0.1^2; from costs.csv
  fig4_temperature_sweep.pdf   Wald coverage and exploration at fixed n; from curve.csv
Appendix
  figA1_exploration_costs_small_delta.pdf   as Figure 3 at V_sp = 0.03^2
  figA2_temperature_sweep_error.pdf         error diagnostics for the Figure 4 runs

Figures are sized for a 6.3-inch text width, so fonts print at roughly their nominal size.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, "..", "..", "results", "separation"))
FIGS = os.path.normpath(os.path.join(HERE, "..", "..", "..", "paper", "tex", "figs"))
DELTA = 0.1
LOSS_LEVELS = [1, 10, 100, 1000]
TEXTWIDTH = 6.3

plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
                     "legend.fontsize": 8, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
                     "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
C_POP, C_AIPW, C_BND, C_LOSS, C_OFF = "#b2182b", "#ef8a62", "#2166ac", "#4d4d4d", "#7b3294"
BLUE, ORANGE, GREEN = "#17628B", "#C56523", "#437F5F"
BOX = dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85)


def save(fig, name, tight=True):
    os.makedirs(FIGS, exist_ok=True)
    fig.savefig(os.path.join(FIGS, name), dpi=300, bbox_inches="tight" if tight else None)
    plt.close(fig)


# ------------------------------------------------------------------ Figure 1
def figure_mechanism():
    h = np.linspace(-1, 1, 2401)
    u1, u2, t = 0.35, 0.85, 1.0
    b = np.where((h >= u1) & (h <= u2), np.sin(np.pi * (h - u1) / (u2 - u1)) ** 2, 0.0)
    hh = np.linspace(u1, u2, 2_000_001)
    theta1 = t * np.trapz(np.sin(np.pi * (hh - u1) / (u2 - u1)) ** 2, hh) / 2   # E[b(H)], H ~ U(-1,1)

    fig, ax = plt.subplots(1, 2, figsize=(TEXTWIDTH, 2.7), constrained_layout=True)
    for a in ax:
        a.axvspan(u1, u2, color=ORANGE, alpha=0.13, lw=0)
        a.axvline(0, color="0.55", ls=":", lw=1)
        a.set_xlabel(r"context score $h$")
        a.set_xlim(-1, 1)
    ax[0].plot(h, np.zeros_like(h), color=BLUE, lw=1.8, ls="--", label=r"Model I: $c^{\mathrm{I}}(h)=0$")
    ax[0].plot(h, t * b, color=ORANGE, lw=2, label=r"Model II: $c^{\mathrm{II}}(h)=b(h)$")
    ax[0].scatter([0], [0], color="black", s=18, zorder=5)
    ax[0].set(ylim=(-0.12, 1.45), ylabel=r"causal effect $c(h)$")
    ax[0].set_title("(a) effects in two outcome models", loc="left")
    ax[0].text(-0.96, 0.98, "boundary effect: 0 in both\n"
               + rf"population effect: 0 (I), {theta1:.3f} (II)", fontsize=8, va="top")
    ax[0].text((u1 + u2) / 2, 1.42, "models differ", fontsize=8, ha="center", va="top", color=ORANGE)
    ax[0].legend(loc="center left", bbox_to_anchor=(0, 0.38), frameon=False)
    for tau, color in ((0.2, BLUE), (0.1, GREEN), (0.05, ORANGE)):
        ax[1].semilogy(h, 1 / (1 + np.exp(np.abs(h) / tau)), color=color, lw=1.8, label=rf"$\tau={tau:g}$")
    ax[1].set(ylim=(1e-9, 1.5), ylabel=r"off-greedy probability $q_\tau(h)$")
    ax[1].set_title("(b) probability of the unpreferred action", loc="left")
    ax[1].text((u1 + u2) / 2, 3e-9, "models\ndiffer", fontsize=8, ha="center", va="bottom", color=ORANGE)
    ax[1].legend(loc="upper left", frameon=False, ncol=1)
    save(fig, "fig1_mechanism.pdf")
    return dict(u1=u1, u2=u2, t=t, theta_model_I=0.0, theta_model_II=float(theta1), beta0=0.0)


# ------------------------------------------------------------------ Figure 2
def figure_map():
    ph = pd.read_csv(os.path.join(RES, "phase.csv")).sort_values("inv_tau")
    inv_tau = ph.inv_tau.to_numpy()
    log_n = np.linspace(2, 9, 281)
    N = 10 ** log_n[:, None]
    rmse_pop = np.sqrt(ph.v_ht.to_numpy()[None, :] / N)
    rmse_bnd = np.sqrt(ph.bias_b.to_numpy()[None, :] ** 2 + ph.v_b.to_numpy()[None, :] / N)
    loss = N * ph.loss_1.to_numpy()[None, :]

    fig, axes = plt.subplots(1, 2, figsize=(TEXTWIDTH, 3.0), sharey=True, constrained_layout=True)
    panels = [(axes[0], rmse_pop, r"(a) population effect $\theta$ (exact, HT)", (1.35, 6.8)),
              (axes[1], rmse_bnd, r"(b) boundary effect $\beta_0$ (first-order)", (17, 6.6))]
    out = {}
    for ax, rmse, title, label_xy in panels:
        z = np.clip(np.log10(rmse), -2.5, 1.0)
        pc = ax.pcolormesh(inv_tau, log_n, z, cmap="viridis_r", shading="auto", vmin=-2.5, vmax=1.0,
                           rasterized=True)
        ax.contour(inv_tau, log_n, rmse, levels=[DELTA], colors="white", linewidths=2.2)
        ax.text(*label_xy, rf"RMSE $\leq {DELTA:g}$", fontsize=8.5, color="black", bbox=BOX)
        ax.contour(inv_tau, log_n, loss, levels=LOSS_LEVELS, colors="#f0f0f0", linewidths=0.9,
                   linestyles="--")
        for lev in LOSS_LEVELS:
            y = np.log10(lev / ph.loss_1.to_numpy()[-1])
            if log_n[0] < y < log_n[-1]:
                ax.text(inv_tau[-1] * 0.9, y + 0.1, rf"$R_n={lev:g}$", ha="right", va="bottom",
                        fontsize=7.5, bbox=BOX)
        feas = rmse <= DELTA
        i, j = np.unravel_index(np.argmin(np.where(feas, loss, np.inf)), loss.shape)
        ax.plot(inv_tau[j], log_n[i], marker="*", ms=12, color="#ffd92f", mec="black", mew=0.7, zorder=5)
        right = j > len(inv_tau) // 2
        ax.annotate(rf"smallest loss $\approx{loss[i, j]:.3g}$", (inv_tau[j], log_n[i]),
                    xytext=(-8 if right else 8, -16), textcoords="offset points", fontsize=8,
                    ha="right" if right else "left", bbox=BOX)
        ax.set_xscale("log")
        ax.set_xlabel(r"policy sharpness $1/\tau$")
        ax.set_title(title, loc="left")
        ax.set_ylim(log_n[0], log_n[-1])
        out["population" if ax is axes[0] else "boundary"] = dict(
            min_loss=float(loss[i, j]), log10_n=float(log_n[i]), inv_tau=float(inv_tau[j]))
    y3 = np.log10(3 * inv_tau)
    show = y3 > log_n[0]
    axes[1].plot(inv_tau[show], y3[show], color="black", lw=1.1, ls=":")
    axes[1].fill_between(inv_tau[show], log_n[0], y3[show], color="white", alpha=0.35, lw=0)
    axes[1].text(inv_tau[-1] * 0.9, log_n[0] + 0.15, r"$n\tau<3$", fontsize=8, ha="right", va="bottom",
                 bbox=BOX)
    axes[0].set_ylabel(r"$\log_{10} n$")
    cb = fig.colorbar(pc, ax=axes, shrink=0.92, pad=0.01)
    cb.set_label(r"$\log_{10}$ RMSE")
    save(fig, "fig2_separation_map.pdf")
    out.update(delta=DELTA, grid_log10_n=[float(log_n[0]), float(log_n[-1])],
               grid_inv_tau=[float(inv_tau[0]), float(inv_tau[-1])])
    return out


# ------------------------------------------------------------------ Figure 3 and A1
def figure_costs(delta, name, width=4.3):
    rows = pd.read_csv(os.path.join(RES, "costs.csv"))
    sel = rows[np.isclose(rows.delta, delta)].sort_values("n")
    fig, ax = plt.subplots(figsize=(width, 3.0), constrained_layout=True)
    for col, lab, color, marker in (("R_common_temperature", "one common temperature", ORANGE, "o"),
                                    ("R_uniform", "uniform mixing", BLUE, "s"),
                                    ("R_optimal", "gap-based allocation", GREEN, "^")):
        ax.loglog(sel.n, sel[col], label=lab, color=color, marker=marker, markersize=4.5, lw=1.8)
    ax.set_title(rf"fixed sparse-exploration criterion $V_{{\rm sp}}={delta:g}^2$", loc="left")
    ax.set_xlabel(r"deployment size $n$")
    ax.set_ylabel(r"cumulative exploration loss $R_n$")
    ax.grid(which="major", color=".9", lw=.6)
    ax.legend(loc="upper left", frameon=False)
    save(fig, name)


# ------------------------------------------------------------------ Figure 4 and A2
def figure_sweep():
    cu = pd.read_csv(os.path.join(RES, "curve.csv")).sort_values("inv_tau")
    x = cu.inv_tau.to_numpy()
    n = int(cu.n.iloc[0])
    reps = int(cu.reps.iloc[0])

    fig, axes = plt.subplots(1, 2, figsize=(TEXTWIDTH, 2.9), constrained_layout=True)
    ax = axes[0]
    for key, col, mk, lab in (("ht", C_POP, "o", r"HT for $\theta$"),
                              ("aipw", C_AIPW, "s", r"fitted AIPW for $\theta$"),
                              ("bnd", C_BND, "^", r"boundary estimator for $\beta_0$")):
        ax.errorbar(x, cu[f"{key}_coverage"], yerr=1.96 * cu[f"{key}_coverage_mcse"], fmt=mk + "-",
                    color=col, ms=4, mfc="none", lw=1.0, capsize=1.5, label=lab)
    ax.axhline(0.95, color="grey", lw=0.8, ls=":")
    ax.set_ylim(-0.02, 1.04)
    ax.set_ylabel("Wald coverage")
    ax.set_title("(a) coverage of 95% intervals", loc="left")
    ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.45), frameon=False)
    top = ax.secondary_xaxis("top", functions=(lambda v: n / np.maximum(v, 1e-12),
                                               lambda v: n / np.maximum(v, 1e-12)))
    top.set_xlabel(r"$n\tau$", fontsize=8.5)
    ax = axes[1]
    ax.plot(x, cu.loss_exact, color=C_LOSS, lw=1.6, label=r"exploration loss $R_n$")
    ax.plot(x, cu.loss_mc, "D", color=C_LOSS, ms=3.5, mfc="none")
    ax.plot(x, cu.offgreedy_exact, color=C_OFF, lw=1.6, ls="--", label="off-greedy actions")
    ax.plot(x, cu.offgreedy_mc, "v", color=C_OFF, ms=3.5, mfc="none")
    ax.set_yscale("log")
    ax.set_title("(b) exploration (lines exact, markers MC)", loc="left")
    ax.legend(loc="lower left", frameon=False)
    for a in axes:
        a.set_xscale("log")
        a.set_xlabel(r"policy sharpness $1/\tau$")
    save(fig, "fig4_temperature_sweep.pdf")

    fig, ax = plt.subplots(figsize=(TEXTWIDTH * 0.8, 3.2), constrained_layout=True)
    ax.plot(x, cu.ht_sd_exact, color=C_POP, lw=1.4, label="HT: exact SD")
    ax.plot(x, cu.aipw_oracle_sd_exact, color=C_AIPW, lw=1.4, ls="-.", label="oracle AIPW (true regression): exact SD")
    ax.plot(x, cu.bnd_rmse_first_order, color=C_BND, lw=1.4, label="boundary estimator: first-order RMSE")
    ax.plot(x, cu.ht_rmse_mc, "o", color=C_POP, ms=4, mfc="none", label="HT: Monte Carlo RMSE")
    ax.plot(x, cu.aipw_rmse_mc, "s", color=C_AIPW, ms=4, mfc="none", label="fitted AIPW: Monte Carlo RMSE")
    ax.plot(x, cu.bnd_rmse_mc, "^", color=C_BND, ms=4, mfc="none", label="boundary estimator: Monte Carlo RMSE")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(5e-3, 1e4)
    ax.set_xlabel(r"policy sharpness $1/\tau$")
    ax.set_ylabel("RMSE or SD")
    ax.set_title(rf"error diagnostics, $n={n:,}$, {reps:,} replications per $\tau$", loc="left")
    ax.legend(loc="upper left", frameon=False, fontsize=7.8)
    save(fig, "figA2_temperature_sweep_error.pdf")


if __name__ == "__main__":
    mech = figure_mechanism()
    with open(os.path.join(RES, "fig_mechanism_numbers.json"), "w") as fh:
        json.dump(mech, fh, indent=2)
    fmap = figure_map()
    with open(os.path.join(RES, "fig1_numbers.json"), "w") as fh:     # file name kept for the table script
        json.dump(fmap, fh, indent=2)
    figure_costs(0.1, "fig3_exploration_costs.pdf")
    figure_costs(0.03, "figA1_exploration_costs_small_delta.pdf")
    if os.path.exists(os.path.join(RES, "curve.csv")):
        figure_sweep()
    print(mech)
    print(fmap)
    print("figures written to", FIGS)
