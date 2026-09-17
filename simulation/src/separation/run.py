"""Numerical checks for the separation results (theory/separation/).

Every subcommand writes <out>/<name>.csv plus <out>/<name>.json with the exact
arguments, seed, git commit, numpy version and wall time, so each table in the
manuscript can be traced to one invocation.

  python run.py costs        # Thm A' constants and Thm B continuous rate (deterministic)
  python run.py finite       # (A.1) tightness in finite contexts; Thm B finite slope (deterministic)
  python run.py kantorovich  # Lemma R bound attained (deterministic)
  python run.py designs      # uniform / optimal / context-temperature / common-temperature, MC
  python run.py boundary     # Thm D: RMSE, coverage of beta_0, exploration loss vs n (tau_n = n^-zeta), MC
  python run.py coverage     # Thm D coverage diagnostics (bias, SE/SD, normality), MC
  python run.py audit        # counterexamples from the 2026-09-17 audit (quadrature)
  python run.py phase        # Figure 1: exact risks and loss over temperatures (quadrature)
  python run.py curve        # Figure 4: estimators vs temperature at fixed n (MC + exact/first-order)
  python run.py firstorder   # accuracy of the first-order RMSE of the boundary estimator (MC)
  python run.py rdcompare    # boundary effect: Theorem 1 estimator vs fuzzy local-linear RD (MC)

Naming: `costs` and `finite` evaluate closed forms / quadrature under the sparse-exploration
design-variance criterion E[sigma^2/p]/n = delta^2; they do not run estimators. `designs`,
`boundary` and `coverage` run estimators by Monte Carlo.

DGP for MC subcommands: H ~ U(-1,1), score Delta = H, greedy a* = 1{H > 0},
operational gap g = |H|, mu_0(h) = h, effect c(h) = 1 + |h| (ATE 1.5, boundary 1),
Y = mu_A(H) + N(0, 1).
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUT = os.path.normpath(os.path.join(HERE, "..", "..", "results", "separation"))
ATE, BETA0 = 1.5, 1.0


def sig(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -700, 700)))


def bisect_decreasing(fun, lo, hi, target, iters=200):
    """Geometric bisection for fun decreasing on [lo, hi]; returns x with fun(x) ~ target."""
    for _ in range(iters):
        mid = np.sqrt(lo * hi)
        if fun(mid) > target:
            lo = mid
        else:
            hi = mid
    return hi


def grid(m=400_000):
    h = np.linspace(-1, 1, m + 1)
    h = h[h != 0.0]
    return h


def git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip()
    except Exception:
        return "unknown"


def save(df, name, args, t0):
    os.makedirs(args.out, exist_ok=True)
    df.to_csv(os.path.join(args.out, f"{name}.csv"), index=False)
    meta = {k: v for k, v in vars(args).items() if k != "func"}
    meta.update(subcommand=name, git_commit=git_commit(), numpy=np.__version__,
                python=sys.version.split()[0], wall_seconds=round(time.time() - t0, 1),
                finished=time.strftime("%Y-%m-%dT%H:%M:%S"))
    with open(os.path.join(args.out, f"{name}.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print(df.to_string(index=False))
    print(f"-> {args.out}/{name}.csv (+ .json)")


# ---------------------------------------------------------------- design helpers
def designs_for(n, delta, s2=1.0, hg=None):
    """Design functions p(h) (off-greedy probability) with E[s2/p]/n = delta^2 (design variance)."""
    hg = grid() if hg is None else hg
    g = np.abs(hg)
    T = n * delta ** 2 / s2
    c = bisect_decreasing(lambda c: np.mean(1 / np.minimum(0.5, c / np.sqrt(g))), 1e-12, 1e3, T)
    t = bisect_decreasing(lambda t: np.mean(1 / sig(-g / t)), 1e-4, 1e4, T)

    def p_opt(h):
        return np.minimum(0.5, c / np.sqrt(np.abs(h) + 1e-300))

    def p_ctx(h):  # Proposition C: context temperature implementing p_opt
        p = p_opt(h)
        with np.errstate(divide="ignore", invalid="ignore"):
            tau_h = np.abs(h) / np.log((1 - p) / p)
        out = sig(-np.abs(h) / np.where(p < 0.5, tau_h, 1.0))
        return np.where(p < 0.5, out, 0.5)

    return {"uniform": lambda h: np.full_like(h, min(0.5, 1 / T)),
            "optimal": p_opt,
            "ctx_temperature": p_ctx,
            "common_temperature": lambda h: sig(-np.abs(h) / t)}, t


# ---------------------------------------------------------------- deterministic
def cmd_costs(args):
    t0 = time.time()
    hg = grid()
    g = np.abs(hg)
    rows = []
    for delta in args.deltas:
        for n in args.ns:
            if n * delta ** 2 < 3:
                continue
            ds, t = designs_for(n, delta, hg=hg)
            cost = {k: n * np.mean(g * f(hg)) for k, f in ds.items()}
            T = n * delta ** 2
            rows.append(dict(delta=delta, n=n, common_tau=t,
                             **{f"R_{k}": v for k, v in cost.items()},
                             R_optimal_closed=np.mean(np.sqrt(g)) ** 2 / delta ** 2,
                             R_uniform_closed=np.mean(g) / delta ** 2,
                             ratio_common_over_optimal=cost["common_temperature"] / cost["optimal"],
                             nd2_over_log2=T / np.log(T) ** 2))
    save(pd.DataFrame(rows), "costs", args, t0)


def cmd_finite(args):
    t0 = time.time()
    f = np.array(args.f)
    gap = np.array(args.gap)
    s2, B, delta = 1.0, args.B, args.delta
    rows = [dict(kind="A1_bound", n=np.nan,
                 value=s2 * (f @ np.sqrt(gap)) ** 2 / delta ** 2 - s2 * np.pi ** 2 / B ** 2 * gap.sum()),
            dict(kind="asymptotic_constant_closed_form", n=np.nan,
                 value=s2 * (f @ np.sqrt(gap)) ** 2 / delta ** 2)]
    prev = None
    for k in range(3, 10):
        n = 10 ** k
        t = bisect_decreasing(lambda t: (f * s2 / (n * sig(-gap / t))).sum(), 1e-4, 1e3, delta ** 2)
        R = n * (f * gap * sig(-gap / t)).sum()
        rows.append(dict(kind="common_temperature_R", n=n, value=R, tau=t,
                         local_slope=np.nan if prev is None else np.log10(R / prev)))
        prev = R
    rows.append(dict(kind="predicted_slope", n=np.nan, value=1 - gap.min() / gap.max()))
    save(pd.DataFrame(rows), "finite", args, t0)


def cmd_kantorovich(args):
    t0 = time.time()
    g = np.abs(grid())
    opt = np.mean(np.sqrt(g)) ** 2
    rows = []
    for a, b in ((0.5, 2), (0.25, 4), (0.1, 10), (0.01, 100)):
        worst = 0.0
        for cut in np.linspace(0, 1, 401):
            for low_first in (True, False):
                r = np.where(g < cut, a, b) if low_first else np.where(g < cut, b, a)
                gh = r * g
                worst = max(worst, np.mean(g / np.sqrt(gh)) * np.mean(np.sqrt(gh)) / opt)
        rho = b / a
        rows.append(dict(a=a, b=b, rho=rho, worst_ratio=worst,
                         K_rho=(1 + np.sqrt(rho)) ** 2 / (4 * np.sqrt(rho)),
                         uniform_over_optimal=np.mean(g) / opt))
    save(pd.DataFrame(rows), "kantorovich", args, t0)


# ---------------------------------------------------------------- Monte Carlo
def cov_rate(covered):
    return float(np.mean(covered))


def cov_mcse(covered):
    """Monte Carlo SE of an observed coverage proportion (uses the observed rate, not 0.95)."""
    c = float(np.mean(covered))
    return float(np.sqrt(c * (1 - c) / len(covered)))


def aipw_once(h, u_sel, noise, p):
    """One replication given shared draws (common random numbers across designs)."""
    greedy = h > 0
    off = u_sel < p
    a = np.where(off, ~greedy, greedy)
    e1 = np.where(greedy, 1 - p, p)
    y = h + (1 + np.abs(h)) * a + noise
    X = np.column_stack([np.ones_like(h), h, np.abs(h)])
    b1 = np.linalg.lstsq(X[a], y[a], rcond=None)[0] if a.sum() > 3 else np.zeros(3)
    b0 = np.linalg.lstsq(X[~a], y[~a], rcond=None)[0] if (~a).sum() > 3 else np.zeros(3)
    m1, m0 = X @ b1, X @ b0
    psi = m1 - m0 + a / e1 * (y - m1) - (~a) / (1 - e1) * (y - m0)
    return psi.mean(), psi.std(ddof=1) / np.sqrt(len(h)), (np.abs(h) * off).sum(), off.sum()


def cmd_designs(args):
    """Designs compared on common random numbers: in each replication all designs share the
    same H, the same selection uniforms and the same outcome noise. `optimal` and
    `ctx_temperature` are the same probability design (Prop. C); the second is an
    implementation check, and max |p_opt - p_ctx| is recorded."""
    t0 = time.time()
    rng = np.random.default_rng(args.seed)
    hg = grid(200_000)
    rows = []
    for n in args.ns:
        deltas = [n ** -0.25] if args.rate else args.deltas
        for delta in deltas:
            ds, t = designs_for(n, delta, hg=hg)
            impl_gap = float(np.max(np.abs(ds["optimal"](hg) - ds["ctx_temperature"](hg))))
            res = {name: [] for name in ds}
            for _ in range(args.reps):
                h = rng.uniform(-1, 1, n)
                u_sel = rng.uniform(size=n)
                noise = rng.normal(size=n)
                for name, pf in ds.items():
                    res[name].append(aipw_once(h, u_sel, noise, pf(h)))
            for name, r in res.items():
                r = np.array(r)
                est, se = r[:, 0], r[:, 1]
                rows.append(dict(n=n, delta=delta, design=name, reps=args.reps, common_tau=t,
                                 max_abs_p_opt_minus_p_ctx=impl_gap,
                                 rmse=np.sqrt(np.mean((est - ATE) ** 2)), bias=est.mean() - ATE,
                                 coverage=cov_rate(np.abs(est - ATE) <= 1.96 * se),
                                 coverage_mcse=cov_mcse(np.abs(est - ATE) <= 1.96 * se),
                                 explored=r[:, 3].mean(), cum_loss=r[:, 2].mean()))
    save(pd.DataFrame(rows), "designs_rate" if args.rate else "designs", args, t0)


def boundary_rep(rng, n, tau):
    h = rng.uniform(-1, 1, n)
    p = sig(h / tau)
    a = rng.uniform(size=n) < p
    y = h + (1 + np.abs(h)) * a + rng.normal(size=n)
    w1, w0 = a * (1 - p), (~a) * p
    if w1.sum() == 0 or w0.sum() == 0:
        return np.nan, np.nan, np.nan, np.nan
    m1, m0 = (w1 * y).sum() / w1.sum(), (w0 * y).sum() / w0.sum()
    psi = w1 * (y - m1) / w1.mean() - w0 * (y - m0) / w0.mean()
    neff = min(w1.sum() ** 2 / (w1 ** 2).sum(), w0.sum() ** 2 / (w0 ** 2).sum())
    loss = (np.abs(h) * (a != (h > 0))).sum()
    return m1 - m0, psi.std(ddof=1) / np.sqrt(n), loss, neff


def cmd_boundary(args):
    t0 = time.time()
    rng = np.random.default_rng(args.seed)
    rows = []
    for zeta in args.zetas:
        for n, reps in zip(args.ns, args.reps_list):
            tau = n ** -zeta
            out = np.array([boundary_rep(rng, n, tau) for _ in range(reps)])
            est, se, loss = out[:, 0], out[:, 1], out[:, 2]
            ok = ~np.isnan(est)
            covered = np.zeros(reps, dtype=bool)
            covered[ok] = np.abs(est[ok] - BETA0) <= 1.96 * se[ok]
            rows.append(dict(zeta=zeta, n=n, n_tau=n * tau, tau=tau, reps=reps, failures=int((~ok).sum()),
                             rmse=np.sqrt(np.mean((est[ok] - BETA0) ** 2)),
                             sd_pred=np.sqrt(2 / (n * tau)),
                             coverage_beta0=cov_rate(covered), coverage_mcse=cov_mcse(covered),
                             cum_loss=np.nanmean(loss), cum_loss_pred=n * tau ** 2 * np.pi ** 2 / 12))
    save(pd.DataFrame(rows), "boundary", args, t0)


def cmd_coverage(args):
    t0 = time.time()
    rng = np.random.default_rng(args.seed)
    hg = np.linspace(-1, 1, 4_000_001)[1:-1]
    rows = []
    for n, reps in zip(args.ns, args.reps_list):
        tau = n ** -args.zeta
        w = sig(hg / tau) * (1 - sig(hg / tau))
        beta_ov = (w * (1 + np.abs(hg))).sum() / w.sum()
        out = np.array([boundary_rep(rng, n, tau) for _ in range(reps)])
        est, se, neff = out[:, 0], out[:, 1], out[:, 3]
        sd = np.nanstd(est, ddof=1)
        z = (est - np.nanmean(est)) / sd
        cov_hat = np.abs(est - BETA0) <= 1.96 * se
        rows.append(dict(zeta=args.zeta, n=n, n_tau=n * tau, tau=tau, reps=reps,
                         kish_neff_per_arm=np.nanmean(neff),
                         bias_beta0=np.nanmean(est) - BETA0, bias_beta_ov=np.nanmean(est) - beta_ov,
                         sd_mc=sd, sd_pred=np.sqrt(2 / (n * tau)), se_over_sd=np.nanmean(se) / sd,
                         coverage_se_hat=cov_rate(cov_hat),
                         coverage_mc_sd=np.nanmean(np.abs(est - BETA0) <= 1.96 * sd),
                         coverage_mcse=cov_mcse(cov_hat),
                         skew=np.nanmean(z ** 3), excess_kurtosis=np.nanmean(z ** 4) - 3))
    save(pd.DataFrame(rows), "coverage", args, t0)


def cmd_audit(args):
    """Counterexamples raised in the 2026-09-17 audit, by quadrature."""
    t0 = time.time()
    rows = []
    # (1) Theorem B scope: mixture of two finite temperatures, 1 and 1/n, keeps cost O(1/delta^2).
    u = np.linspace(0, 1, 2_000_001)[1:]
    q = sig(-u)
    for delta in (0.1,):
        for n in (10 ** 4, 10 ** 5, 10 ** 6):
            w = np.mean(1 / q) / (n * delta ** 2)
            p = w * q + (1 - w) * sig(-n * u)
            rows.append(dict(check="B_mixture_temperatures", n=n, delta=delta, weight_hot=w,
                             design_variance=np.mean(1 / p) / n, cum_loss=n * np.mean(u * p),
                             bound=np.mean(1 / q) * np.mean(u * q) / delta ** 2))
    # (2) Theorem D remark: c(h) = 1 + |h|^{3/2} is differentiable at 0 but has tau^{3/2} bias.
    h = np.linspace(-1, 1, 40_000_001)
    for tau in (0.04, 0.02, 0.01, 0.005):
        om = sig(h / tau) * (1 - sig(h / tau))
        b = (om * np.abs(h) ** 1.5).sum() / om.sum()
        rows.append(dict(check="D_smoothness_remark", tau=tau, bias=b,
                         bias_over_tau_1p5=b / tau ** 1.5, bias_over_tau_2=b / tau ** 2))
    # (3) Lemma R: with atoms the Kantorovich bound need not be attained, so "iff" fails.
    f = np.array([0.99, 0.01]); g = np.array([1.0, 100.0]); rho = 36.0
    w = f * np.sqrt(g); w = w / w.sum()
    s_lo, s_hi = 1.0, np.sqrt(rho)
    worst = max((w @ s) * (w @ (1 / s)) for s in
                (np.array([s_lo, s_hi]), np.array([s_hi, s_lo]), np.array([s_lo, s_lo])))
    rows.append(dict(check="R_iff_with_atoms", uniform_over_optimal=(f @ g) / (f @ np.sqrt(g)) ** 2,
                     worst_misspecified_over_optimal=worst,
                     K_rho=(1 + np.sqrt(rho)) ** 2 / (4 * np.sqrt(rho))))
    save(pd.DataFrame(rows), "audit", args, t0)


# ---------------------------------------------------------------- figures (manuscript Figures 1-2)
def exact_moments(tau, m=4_000_000):
    """Exact (quadrature) moments for the simulation DGP under common temperature tau.

    H ~ U(-1,1), Delta = h, e = sig(h/tau) = P(A=1|H), q = sig(-|h|/tau), g = |h|,
    mu0 = h, mu1 = h + 1 + |h|, sigma = 1.  Returns per-unit quantities; divide by n where noted.
      v_ht     : n * Var(HT estimator of theta)                (exact: HT is unbiased)
      v_aipw   : n * Var(oracle AIPW with the TRUE regression)  (exact; = Var(c) + E[1/e + 1/(1-e)])
      bias_b   : beta_ov,tau - beta_0                           (exact)
      v_b      : n * FIRST-ORDER (linearization) variance of the ratio-type boundary estimator.
                 sqrt(bias_b**2 + v_b/n) is a first-order RMSE approximation, not the exact
                 finite-sample RMSE; see `run.py firstorder` for its accuracy.
      loss_1   : E[g q], so R_n = n * loss_1                     (exact)
      offgreedy: E[q], expected fraction of off-greedy actions  (exact)
    """
    h = -1 + (np.arange(m) + 0.5) * (2.0 / m)          # midpoint rule on (-1, 1), density 1/2
    e = sig(h / tau)
    q = sig(-np.abs(h) / tau)
    mu0, mu1 = h, h + 1 + np.abs(h)
    c = mu1 - mu0
    theta = c.mean()
    inv_e = 1.0 + np.exp(np.clip(-h / tau, -700, 700))    # 1/e without underflow
    inv_1me = 1.0 + np.exp(np.clip(h / tau, -700, 700))   # 1/(1-e)
    v_ht = np.mean((mu1 ** 2 + 1) * inv_e + (mu0 ** 2 + 1) * inv_1me) - theta ** 2
    v_aipw = np.var(c) + np.mean(inv_e + inv_1me)
    w = e * (1 - e)
    mt = w.mean()
    bar1, bar0 = (w * mu1).mean() / mt, (w * mu0).mean() / mt
    bias_b = (w * c).mean() / mt - 1.0
    v_b = np.mean(w * ((1 - e) * (1 + (mu1 - bar1) ** 2) + e * (1 + (mu0 - bar0) ** 2))) / mt ** 2
    return dict(tau=tau, inv_tau=1 / tau, theta=theta, v_ht=v_ht, v_aipw=v_aipw, bias_b=bias_b,
                v_b=v_b, loss_1=np.mean(np.abs(h) * q), offgreedy=q.mean())


def cmd_phase(args):
    """Figure 1 input: exact per-unit moments on a grid of temperatures (no Monte Carlo)."""
    t0 = time.time()
    inv_taus = np.unique(np.round(np.geomspace(args.inv_tau_min, args.inv_tau_max, args.points), 6))
    rows = [exact_moments(1 / it, m=args.quad) for it in inv_taus]
    save(pd.DataFrame(rows), "phase", args, t0)


def _curve_rep(job):
    n, tau, seed = job
    rng = np.random.default_rng(seed)
    h = rng.uniform(-1, 1, n)
    e = sig(h / tau)
    a = rng.uniform(size=n) < e
    y = h + (1 + np.abs(h)) * a + rng.normal(size=n)
    out = {}
    # population effect: Horvitz-Thompson with known e. Inverse probabilities are computed as
    # 1 + exp(-+h/tau) and used only for the arm actually taken, so no 0 * inf arises when e
    # underflows (such arms are never observed).
    with np.errstate(over="ignore"):
        inv_e = 1.0 + np.exp(np.clip(-h / tau, -700, 700))
        inv_1me = 1.0 + np.exp(np.clip(h / tau, -700, 700))
    w_1 = np.where(a, inv_e, 0.0)
    w_0 = np.where(a, 0.0, inv_1me)
    psi = w_1 * y - w_0 * y
    out["ht"], out["ht_se"] = psi.mean(), psi.std(ddof=1) / np.sqrt(n)
    # population effect: AIPW with a linear working model in (1, h, |h|) per arm
    X = np.column_stack([np.ones(n), h, np.abs(h)])
    if a.sum() > 3 and (~a).sum() > 3:
        b1 = np.linalg.lstsq(X[a], y[a], rcond=None)[0]
        b0 = np.linalg.lstsq(X[~a], y[~a], rcond=None)[0]
        m1, m0 = X @ b1, X @ b0
        psi = m1 - m0 + w_1 * (y - m1) - w_0 * (y - m0)
        out["aipw"], out["aipw_se"] = psi.mean(), psi.std(ddof=1) / np.sqrt(n)
    else:
        out["aipw"], out["aipw_se"] = np.nan, np.nan
    # boundary effect: overlap-weighted Hajek estimator of Theorem 1
    w1, w0 = a * (1 - e), (~a) * e
    if w1.sum() > 0 and w0.sum() > 0:
        mm1, mm0 = (w1 * y).sum() / w1.sum(), (w0 * y).sum() / w0.sum()
        psi = w1 * (y - mm1) / w1.mean() - w0 * (y - mm0) / w0.mean()
        out["bnd"], out["bnd_se"] = mm1 - mm0, psi.std(ddof=1) / np.sqrt(n)
    else:
        out["bnd"], out["bnd_se"] = np.nan, np.nan
    out["loss"] = (np.abs(h) * (a != (h > 0))).sum()
    out["offgreedy"] = (a != (h > 0)).sum()
    return out


def cmd_curve(args):
    """Figure 2 input: Monte Carlo at fixed n over temperatures, plus exact moments."""
    from multiprocessing import Pool
    t0 = time.time()
    ss = np.random.SeedSequence(args.seed)
    rows = []
    with Pool(args.procs) as pool:
        for inv_tau in args.inv_taus:
            tau = 1.0 / inv_tau
            seeds = [int(s.generate_state(1)[0]) for s in ss.spawn(args.reps)]
            reps = pool.map(_curve_rep, [(args.n, tau, sd) for sd in seeds], chunksize=4)
            df = pd.DataFrame(reps)
            ex = exact_moments(tau)
            row = dict(n=args.n, inv_tau=inv_tau, tau=tau, reps=args.reps, n_tau=args.n * tau)
            for key, target in (("ht", ATE), ("aipw", ATE), ("bnd", BETA0)):
                est, se = df[key].to_numpy(), df[key + "_se"].to_numpy()
                ok = np.isfinite(est) & np.isfinite(se)
                covered = np.zeros(len(est), dtype=bool)
                covered[ok] = np.abs(est[ok] - target) <= 1.96 * se[ok]
                row[f"{key}_rmse_mc"] = float(np.sqrt(np.mean((est[ok] - target) ** 2))) if ok.any() else np.nan
                row[f"{key}_sd_mc"] = float(np.std(est[ok], ddof=1)) if ok.sum() > 1 else np.nan
                row[f"{key}_coverage"] = cov_rate(covered)
                row[f"{key}_coverage_mcse"] = cov_mcse(covered)
                row[f"{key}_failures"] = int((~ok).sum())
            row.update(ht_sd_exact=np.sqrt(ex["v_ht"] / args.n),
                       aipw_oracle_sd_exact=np.sqrt(ex["v_aipw"] / args.n),
                       bnd_bias_exact=ex["bias_b"], bnd_sd_first_order=np.sqrt(ex["v_b"] / args.n),
                       bnd_rmse_first_order=np.sqrt(ex["bias_b"] ** 2 + ex["v_b"] / args.n),
                       loss_mc=float(df["loss"].mean()), loss_exact=args.n * ex["loss_1"],
                       offgreedy_mc=float(df["offgreedy"].mean()), offgreedy_exact=args.n * ex["offgreedy"])
            rows.append(row)
            print(f"1/tau={inv_tau}: done", flush=True)
    save(pd.DataFrame(rows), "curve", args, t0)


def _firstorder_rep(job):
    n, tau, seed = job
    rng = np.random.default_rng(seed)
    h = rng.uniform(-1, 1, n)
    e = sig(h / tau)
    a = rng.uniform(size=n) < e
    y = h + (1 + np.abs(h)) * a + rng.normal(size=n)
    w1, w0 = a * (1 - e), (~a) * e
    if w1.sum() == 0 or w0.sum() == 0:
        return np.nan              # undefined; the manuscript's convention sets the estimate to 0
    return (w1 * y).sum() / w1.sum() - (w0 * y).sum() / w0.sum()


def cmd_firstorder(args):
    """Accuracy of the first-order RMSE approximation for the boundary estimator (Figure 2b)."""
    from multiprocessing import Pool
    t0 = time.time()
    ss = np.random.SeedSequence(args.seed)
    rows = []
    with Pool(args.procs) as pool:
        for n in args.ns:
            for inv_tau in args.inv_taus:
                tau = 1.0 / inv_tau
                ex = exact_moments(tau, m=2_000_000)
                fo = float(np.sqrt(ex["bias_b"] ** 2 + ex["v_b"] / n))
                seeds = [int(s.generate_state(1)[0]) for s in ss.spawn(args.reps)]
                est = np.array(pool.map(_firstorder_rep, [(n, tau, sd) for sd in seeds], chunksize=8))
                ok = ~np.isnan(est)
                conv = np.where(ok, est, 0.0)       # manuscript convention: estimate 0 when undefined
                mc_conv = float(np.sqrt(np.mean((conv - BETA0) ** 2)))
                mc_def = float(np.sqrt(np.mean((est[ok] - BETA0) ** 2)))
                rows.append(dict(n=n, inv_tau=inv_tau, n_tau=n * tau, reps=args.reps,
                                 rmse_first_order=fo,
                                 rmse_mc=mc_conv, ratio_mc_over_first_order=mc_conv / fo,
                                 rmse_mc_defined_only=mc_def, ratio_defined_only=mc_def / fo,
                                 undefined_rate=float(1 - ok.mean())))
    save(pd.DataFrame(rows), "firstorder", args, t0)



def _rd_rep(job):
    n, tau, bw, seed = job
    rng = np.random.default_rng(seed)
    h = rng.uniform(-1, 1, n)
    e = sig(h / tau)
    a = rng.uniform(size=n) < e
    y = h + (1 + np.abs(h)) * a + rng.normal(size=n)
    res = {}
    # Theorem 1 estimator (known propensities)
    w1, w0 = a * (1 - e), (~a) * e
    if w1.sum() > 0 and w0.sum() > 0:
        m1, m0 = (w1 * y).sum() / w1.sum(), (w0 * y).sum() / w0.sum()
        psi = w1 * (y - m1) / w1.mean() - w0 * (y - m0) / w0.mean()
        res["ow"], res["ow_se"] = m1 - m0, psi.std(ddof=1) / np.sqrt(n)
    else:
        res["ow"], res["ow_se"] = 0.0, np.inf
    # fuzzy local-linear RD at 0 (ignores the known propensities); triangular kernel, bandwidth bw
    fits = []
    for side in (h >= 0, h < 0):
        xs, ys, as_ = h[side], y[side], a[side].astype(float)
        k = np.clip(1 - np.abs(xs) / bw, 0, None)
        keep = k > 0
        xs, ys, as_, k = xs[keep], ys[keep], as_[keep], k[keep]
        X = np.column_stack([np.ones_like(xs), xs])
        M = np.linalg.inv((X.T * k) @ X) @ (X.T * k)    # coefficients = M @ v
        fits.append((X, M, ys, as_))
    jy = fits[0][1][0] @ fits[0][2] - fits[1][1][0] @ fits[1][2]
    ja = fits[0][1][0] @ fits[0][3] - fits[1][1][0] @ fits[1][3]
    if abs(ja) < 1e-8:
        res["rd"], res["rd_se"] = 0.0, np.inf
        return res
    est = jy / ja
    var = 0.0
    for X, M, ys, as_ in fits:                          # delta method: Z = Y - est * A
        z = ys - est * as_
        r = z - X @ (M @ z)
        var += np.sum((M[0] * r) ** 2)                   # HC0 variance of the intercept
    res["rd"], res["rd_se"] = est, np.sqrt(var) / abs(ja)
    return res


def cmd_rdcompare(args):
    """Boundary effect: Theorem 1 estimator (known propensities) vs fuzzy local-linear RD."""
    from multiprocessing import Pool
    t0 = time.time()
    ss = np.random.SeedSequence(args.seed)
    rows = []
    with Pool(args.procs) as pool:
        for inv_tau in args.inv_taus:
            tau = 1.0 / inv_tau
            for bw in args.bandwidths:
                seeds = [int(s.generate_state(1)[0]) for s in ss.spawn(args.reps)]
                df = pd.DataFrame(pool.map(_rd_rep, [(args.n, tau, bw, sd) for sd in seeds], chunksize=4))
                row = dict(n=args.n, inv_tau=inv_tau, tau=tau, n_tau=args.n * tau, bandwidth=bw, reps=args.reps)
                for key in ("ow", "rd"):
                    est, se = df[key].to_numpy(), df[key + "_se"].to_numpy()
                    covered = np.abs(est - BETA0) <= 1.96 * se
                    row[f"{key}_bias"] = float(np.mean(est) - BETA0)
                    row[f"{key}_rmse"] = float(np.sqrt(np.mean((est - BETA0) ** 2)))
                    row[f"{key}_coverage"] = cov_rate(covered)
                    row[f"{key}_coverage_mcse"] = cov_mcse(covered)
                rows.append(row)
                print(f"1/tau={inv_tau} bw={bw}: done", flush=True)
    save(pd.DataFrame(rows), "rdcompare", args, t0)


def cmd_rdlimit(args):
    """Probability limits of both boundary estimators in the rdcompare design (quadrature, no Monte Carlo).

    With finite tau the assignment probability is continuous at 0, so the fuzzy RD estimator has no
    discontinuity to exploit; its limit is the ratio of the population local-linear intercept gaps
    of E[Y|H] = h + e_tau(h)(1 + |h|) and E[A|H] = e_tau(h), which differs from beta_0 = 1."""
    t0 = time.time()
    rows = []
    for inv_tau in args.inv_taus:
        e = lambda h: sig(h * inv_tau)
        my = lambda h: h + e(h) * (1 + np.abs(h))
        hh = np.linspace(-1, 1, 400001)
        om = e(hh) * (1 - e(hh))
        ow_limit = float(np.sum(om * (1 + np.abs(hh))) / np.sum(om))
        for bw in args.bandwidths:
            gaps = {}
            for name, fun in (("y", my), ("a", e)):
                icpt = []
                for sgn in (1, -1):
                    h = np.linspace(0, bw, 200001)[1:] * sgn
                    k = 1 - np.abs(h) / bw
                    X = np.column_stack([np.ones_like(h), h])
                    icpt.append(np.linalg.solve((X.T * k) @ X, (X.T * k) @ fun(h))[0])
                gaps[name] = icpt[0] - icpt[1]
            rows.append(dict(inv_tau=inv_tau, bandwidth=bw, first_stage_gap=gaps["a"],
                             rd_limit=gaps["y"] / gaps["a"], ow_limit=ow_limit, beta0=BETA0))
    save(pd.DataFrame(rows), "rdlimit", args, t0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=DEFAULT_OUT)
    sub = ap.add_subparsers(required=True)

    s = sub.add_parser("costs"); s.set_defaults(func=cmd_costs)
    s.add_argument("--ns", type=int, nargs="+", default=[10 ** k for k in range(4, 9)])
    s.add_argument("--deltas", type=float, nargs="+", default=[0.1, 0.03])

    s = sub.add_parser("finite"); s.set_defaults(func=cmd_finite)
    s.add_argument("--f", type=float, nargs="+", default=[0.5, 0.5])
    s.add_argument("--gap", type=float, nargs="+", default=[0.2, 1.0])
    s.add_argument("--B", type=float, default=5.0)
    s.add_argument("--delta", type=float, default=0.1)

    s = sub.add_parser("kantorovich"); s.set_defaults(func=cmd_kantorovich)

    s = sub.add_parser("designs"); s.set_defaults(func=cmd_designs)
    s.add_argument("--ns", type=int, nargs="+", default=[100_000])
    s.add_argument("--deltas", type=float, nargs="+", default=[0.1, 0.05])
    s.add_argument("--rate", action="store_true", help="use delta_n = n^-1/4 instead of --deltas")
    s.add_argument("--reps", type=int, default=400)
    s.add_argument("--seed", type=int, default=20260917)

    s = sub.add_parser("boundary"); s.set_defaults(func=cmd_boundary)
    s.add_argument("--zetas", type=float, nargs="+", default=[0.6, 0.75],
                   help="temperature exponents: tau_n = n^-zeta")
    s.add_argument("--ns", type=int, nargs="+", default=[10 ** 4, 10 ** 5, 10 ** 6])
    s.add_argument("--reps-list", type=int, nargs="+", default=[1000, 500, 200])
    s.add_argument("--seed", type=int, default=20260918)

    s = sub.add_parser("audit"); s.set_defaults(func=cmd_audit)

    s = sub.add_parser("phase"); s.set_defaults(func=cmd_phase)
    s.add_argument("--inv-tau-min", type=float, default=1.0)
    s.add_argument("--inv-tau-max", type=float, default=300.0)
    s.add_argument("--points", type=int, default=90)
    s.add_argument("--quad", type=int, default=4_000_000, help="midpoint-rule points on (-1, 1)")

    s = sub.add_parser("firstorder"); s.set_defaults(func=cmd_firstorder)
    s.add_argument("--ns", type=int, nargs="+", default=[100, 1000, 10000])
    s.add_argument("--inv-taus", type=float, nargs="+", default=[1, 3, 10, 30, 100, 300])
    s.add_argument("--reps", type=int, default=4000)
    s.add_argument("--seed", type=int, default=20260922)
    s.add_argument("--procs", type=int, default=min(64, os.cpu_count() or 8))

    s = sub.add_parser("rdlimit"); s.set_defaults(func=cmd_rdlimit)
    s.add_argument("--inv-taus", type=float, nargs="+", default=[10, 30, 120, 1000])
    s.add_argument("--bandwidths", type=float, nargs="+", default=[0.05, 0.2])

    s = sub.add_parser("rdcompare"); s.set_defaults(func=cmd_rdcompare)
    s.add_argument("--n", type=int, default=100_000)
    s.add_argument("--inv-taus", type=float, nargs="+", default=[10, 30, 120, 1000])
    s.add_argument("--bandwidths", type=float, nargs="+", default=[0.05, 0.2])
    s.add_argument("--reps", type=int, default=1000)
    s.add_argument("--seed", type=int, default=20260923)
    s.add_argument("--procs", type=int, default=min(64, os.cpu_count() or 8))

    s = sub.add_parser("curve"); s.set_defaults(func=cmd_curve)
    s.add_argument("--n", type=int, default=100_000)
    s.add_argument("--inv-taus", type=float, nargs="+",
                   default=[1, 2, 3, 5, 7, 10, 14, 20, 30, 50, 80, 120])
    s.add_argument("--reps", type=int, default=1000)
    s.add_argument("--seed", type=int, default=20260921)
    s.add_argument("--procs", type=int, default=min(64, os.cpu_count() or 8))

    s = sub.add_parser("coverage"); s.set_defaults(func=cmd_coverage)
    s.add_argument("--zeta", type=float, default=0.6, help="tau_n = n^-zeta")
    s.add_argument("--ns", type=int, nargs="+", default=[10 ** 4, 10 ** 5, 10 ** 6])
    s.add_argument("--reps-list", type=int, nargs="+", default=[3000, 3000, 1000])
    s.add_argument("--seed", type=int, default=20260919)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
