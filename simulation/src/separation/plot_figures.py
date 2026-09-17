"""Draw the four manuscript figures.

  python plot_figures.py

  fig1_mechanism.pdf        exact functions of the Theorem 3 construction (no data)
  fig2_separation_map.pdf   exact risks and loss contours from results/separation/phase.csv
  fig3_exploration_costs.pdf closed-form design costs from results/separation/costs.csv
  fig4_temperature_sweep.pdf fixed-n Monte Carlo and exact values from results/separation/curve.csv

Figures 1 and 3 follow a review revision (source bundle based on commit df5007f); Figure 1's
axis label bug (a tab character in "q_tau") is fixed here and the bump integral is computed with
numpy instead of scipy.
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

plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
                     "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
                     "pdf.fonttype": 42})
C_POP, C_AIPW, C_BND, C_LOSS = "#b2182b", "#ef8a62", "#2166ac", "#4d4d4d"


def figure_mechanism():
    """Theorem 3 construction: same boundary effect, different population effect."""
    h = np.linspace(-1, 1, 2401)
    u1, u2, t = 0.35, 0.85, 1.0
    b = np.where((h >= u1) & (h <= u2), np.sin(np.pi * (h - u1) / (u2 - u1)) ** 2, 0.0)
    hh = np.linspace(u1, u2, 2_000_001)
    theta1 = t * np.trapz(np.sin(np.pi * (hh - u1) / (u2 - u1)) ** 2, hh) / 2   # E[b(H)], H ~ U(-1,1)
    blue, orange, green = "#17628B", "#C56523", "#437F5F"
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)
    for a in ax:
        a.axvspan(u1, u2, color=orange, alpha=0.12, lw=0)
        a.axvline(0, color="0.55", ls=":", lw=1)
        a.set_xlabel(r"score / context $h$")
        a.spines[["top", "right"]].set_visible(False)
    ax[0].plot(h, t * b, color=orange, lw=2, label=r"$c^1(h)=b(h)$")
    ax[0].plot(h, np.zeros_like(h), color=blue, lw=1.7, ls="--", label=r"$c^0(h)=0$")
    ax[0].scatter([0], [0], color="black", s=20, zorder=5)
    ax[0].set(xlim=(-1, 1), ylim=(-0.12, 1.28), ylabel=r"causal effect $c(h)$")
    ax[0].set_title("(a) same boundary, different population", loc="left")
    ax[0].text(-0.95, 0.86, r"$\beta_0^0=\beta_0^1=0$" + "\n" + rf"$\theta^0=0,\quad\theta^1={theta1:.3f}$",
               fontsize=9)
    ax[0].legend(loc="upper left", frameon=False, bbox_to_anchor=(0, 0.60))
    for tau, color in ((0.2, blue), (0.1, green), (0.05, orange)):
        ax[1].semilogy(h, 1 / (1 + np.exp(np.abs(h) / tau)), color=color, lw=1.8, label=rf"$\tau={tau:g}$")
    ax[1].set(xlim=(-1, 1), ylim=(1e-9, 1), ylabel=r"off-greedy probability $q_\tau(h)$")
    ax[1].set_title("(b) off-greedy probability", loc="left")
    ax[1].legend(loc="lower center", frameon=False, ncol=3, columnspacing=0.8)
    fig.savefig(os.path.join(FIGS, "fig1_mechanism.pdf"), dpi=300)
    plt.close(fig)
    return dict(u1=u1, u2=u2, t=t, theta0=0.0, theta1=float(theta1), beta0=0.0)


def figure_costs():
    """Closed-form exploration loss vs deployment size at a fixed sparse-exploration criterion."""
    rows = pd.read_csv(os.path.join(RES, "costs.csv"))
    blue, orange, green = "#17628B", "#C56523", "#437F5F"
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.2), constrained_layout=True)
    for ax, delta, label in zip(axes, (0.1, 0.03), ("(a)", "(b)")):
        sel = rows[np.isclose(rows.delta, delta)].sort_values("n")
        for col, lab, color, marker in (("R_common_temperature", "one common temperature", orange, "o"),
                                        ("R_uniform", "uniform mixing", blue, "s"),
                                        ("R_optimal", "gap-based allocation", green, "^")):
            ax.loglog(sel.n, sel[col], label=lab, color=color, marker=marker, markersize=4, lw=1.8)
        ax.set_title(label + rf" fixed $V_{{\rm sp}}={delta:g}^2$", loc="left")
        ax.set_xlabel(r"deployment size $n$")
        ax.set_ylabel(r"expected cumulative loss $R_n$")
        ax.set_ylim(25 if delta == 0.1 else 300, 1e6)
        ax.grid(which="major", color=".9", lw=.6)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(loc="upper left", frameon=False, fontsize=7.4)
    fig.savefig(os.path.join(FIGS, "fig3_exploration_costs.pdf"), dpi=300)
    plt.close(fig)


def figure1():
    ph = pd.read_csv(os.path.join(RES, "phase.csv")).sort_values("inv_tau")
    inv_tau = ph.inv_tau.to_numpy()
    log_n = np.linspace(2, 9, 281)
    N = 10 ** log_n[:, None]
    rmse_pop = np.sqrt(ph.v_ht.to_numpy()[None, :] / N)
    rmse_bnd = np.sqrt(ph.bias_b.to_numpy()[None, :] ** 2 + ph.v_b.to_numpy()[None, :] / N)
    loss = N * ph.loss_1.to_numpy()[None, :]

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.3), sharey=True, constrained_layout=True)
    panels = [(axes[0], rmse_pop, r"(a) population effect $\theta$ (HT)"),
              (axes[1], rmse_bnd, r"(b) boundary effect $\beta_0$ (Theorem 1)")]
    for ax, rmse, title in panels:
        z = np.clip(np.log10(rmse), -2.5, 1.0)                 # dark = large error
        pc = ax.pcolormesh(inv_tau, log_n, z, cmap="viridis_r", shading="auto", vmin=-2.5, vmax=1.0,
                           rasterized=True)
        ax.contourf(inv_tau, log_n, (rmse <= DELTA).astype(float), levels=[0.5, 1.5],
                    hatches=["////"], colors="none")
        ax.contour(inv_tau, log_n, rmse, levels=[DELTA], colors="white", linewidths=2.0)
        ax.contour(inv_tau, log_n, loss, levels=LOSS_LEVELS, colors="#f0f0f0", linewidths=0.9,
                   linestyles="--")
        for lev in LOSS_LEVELS:                      # label each loss contour at the right edge
            y = np.log10(lev / ph.loss_1.to_numpy()[-1])
            if log_n[0] < y < log_n[-1]:
                ax.text(inv_tau[-1] * 0.93, y + 0.12, rf"$R_n={lev:g}$", ha="right", va="bottom",
                        fontsize=6.5, color="black",
                        bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.8))
        feas = rmse <= DELTA
        if feas.any():                               # smallest exploration loss meeting the precision
            i, j = np.unravel_index(np.argmin(np.where(feas, loss, np.inf)), loss.shape)
            ax.plot(inv_tau[j], log_n[i], marker="*", ms=11, color="#ffd92f", mec="black", mew=0.7,
                    zorder=5)
            right = j > len(inv_tau) // 2
            ax.annotate(rf"min $R_n\approx{loss[i, j]:.3g}$", (inv_tau[j], log_n[i]),
                        xytext=(-8 if right else 6, -14), textcoords="offset points", fontsize=7,
                        ha="right" if right else "left",
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
        ax.set_xscale("log")
        ax.set_xlabel(r"sharpness $1/\tau$ (one common temperature)")
        ax.set_title(title)
    axes[0].set_ylabel(r"$\log_{10} n$")
    cb = fig.colorbar(pc, ax=axes, shrink=0.9, pad=0.01)
    cb.set_label(r"$\log_{10}$ RMSE (exact)")
    os.makedirs(FIGS, exist_ok=True)
    fig.savefig(os.path.join(FIGS, "fig2_separation_map.pdf"), dpi=300)
    plt.close(fig)

    # numbers quoted with the figure
    feas_pop, feas_bnd = rmse_pop <= DELTA, rmse_bnd <= DELTA
    def argmin_point(feas):
        i, j = np.unravel_index(np.argmin(np.where(feas, loss, np.inf)), loss.shape)
        return dict(min_loss=float(loss[i, j]), log10_n=float(log_n[i]), inv_tau=float(inv_tau[j]))
    out = dict(delta=DELTA, grid_log10_n=[float(log_n[0]), float(log_n[-1])],
               grid_inv_tau=[float(inv_tau[0]), float(inv_tau[-1])],
               population=argmin_point(feas_pop), boundary=argmin_point(feas_bnd))
    with open(os.path.join(RES, "fig1_numbers.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    return out


def figure2():
    cu = pd.read_csv(os.path.join(RES, "curve.csv")).sort_values("inv_tau")
    x = cu.inv_tau.to_numpy()
    n = int(cu.n.iloc[0])
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.7), constrained_layout=True)

    ax = axes[0]
    ax.plot(x, cu.ht_sd_exact, color=C_POP, lw=1.3, label=r"HT, exact SD")
    ax.plot(x, cu.aipw_sd_exact, color=C_AIPW, lw=1.3, ls="-.", label=r"AIPW, exact SD")
    ax.plot(x, cu.bnd_rmse_exact, color=C_BND, lw=1.3, label=r"boundary, exact RMSE")
    ax.plot(x, cu.ht_rmse_mc, "o", color=C_POP, ms=3.5, mfc="none")
    ax.plot(x, cu.aipw_rmse_mc, "s", color=C_AIPW, ms=3.2, mfc="none")
    ax.plot(x, cu.bnd_rmse_mc, "^", color=C_BND, ms=3.5, mfc="none")
    ax.set_yscale("log")
    ax.set_ylim(5e-3, 1e4)
    ax.set_ylabel("RMSE / SD")
    ax.set_title("(a) error (lines exact, markers MC)")
    ax.legend(loc="upper left", frameon=False)

    ax = axes[1]
    for key, col, mk, lab in (("ht", C_POP, "o", r"HT for $\theta$"), ("aipw", C_AIPW, "s", r"AIPW for $\theta$"),
                              ("bnd", C_BND, "^", r"boundary for $\beta_0$")):
        ax.errorbar(x, cu[f"{key}_coverage"], yerr=1.96 * cu[f"{key}_coverage_mcse"], fmt=mk, color=col,
                    ms=3.5, mfc="none", lw=0.8, capsize=1.5, label=lab)
    ax.axhline(0.95, color="grey", lw=0.8, ls=":")
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Wald coverage (MC)")
    ax.set_title("(b) coverage")
    ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.42), frameon=False)

    ax = axes[2]
    ax.plot(x, cu.loss_exact, color=C_LOSS, lw=1.3, label=r"loss $R_n$, exact")
    ax.plot(x, cu.loss_mc, "D", color=C_LOSS, ms=3, mfc="none")
    ax.plot(x, cu.offgreedy_exact, color="#7b3294", lw=1.3, ls="--", label="off-greedy actions, exact")
    ax.plot(x, cu.offgreedy_mc, "v", color="#7b3294", ms=3, mfc="none")
    ax.set_yscale("log")
    ax.set_title("(c) exploration")
    ax.legend(loc="lower left", frameon=False)

    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel(r"$1/\tau$")
    fig.suptitle(rf"Common temperature, $n={n:,}$, {int(cu.reps.iloc[0]):,} replications per $\tau$",
                 fontsize=9)
    fig.savefig(os.path.join(FIGS, "fig4_temperature_sweep.pdf"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs(FIGS, exist_ok=True)
    mech = figure_mechanism()
    with open(os.path.join(RES, "fig_mechanism_numbers.json"), "w") as fh:
        json.dump(mech, fh, indent=2)
    print(mech)
    figure_costs()
    print(figure1())
    if os.path.exists(os.path.join(RES, "curve.csv")):
        figure2()
    print("figures written to", FIGS)
