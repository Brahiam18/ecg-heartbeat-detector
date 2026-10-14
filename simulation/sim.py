"""Behavioral simulation of the analog ECG heartbeat detector.

Models the signal chain stage by stage -- gain stage (TP1), precision
full-wave rectifier (TP2/TP3), peak-envelope detector with R9*C1 discharge
(TP4) and Schmitt-trigger comparator (TP5) -- and evaluates it with
rectangular pulses, sinusoids and a synthetic PWL QRS complex. It also
computes theoretical envelope discharge times / maximum heart rate and a
rectifier resistor-tolerance analysis (worst case + Monte Carlo).

Usage:
    pip install numpy matplotlib
    python sim.py

Figures (f1_rectificador.png ... f7_qrs_pwlin.png) and results.json are
written to the fig/ folder next to this script; a summary is printed.
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, json, os

OUT = os.path.dirname(os.path.abspath(__file__)) + "/fig"
os.makedirs(OUT, exist_ok=True)

# ---------- design parameters ----------
G = 1 + 39e3/10e3            # 4.9
VSAT = 10.5
R7, R9, C1 = 100.0, 100e3, 1e-6
TAU = R9*C1                  # 0.1 s
TAUC = R7*C1                 # 0.1 ms
R10, R11 = 10e3, 430e3
VREF = 12*1e3/(15e3+1e3)     # 0.75
k = R10/R11
VTH = VREF*(1+k) + VSAT*k
VTL = VREF*(1+k) - VSAT*k
FS = 20000.0
DT = 1/FS

BLUE, ORANGE, AQUA, YEL, MAG, VIO = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e6e5e0"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.8, "figure.facecolor": "white", "axes.facecolor": "white",
                     "legend.frameon": False, "lines.linewidth": 1.8})

def chain(vin):
    vg = np.clip(G*vin, -VSAT, VSAT)
    vhw = np.where(vg > 0, -vg, 0.0)
    vr = np.clip(-(vg + 2*vhw), -VSAT, VSAT)
    env = np.zeros_like(vr); vo = np.zeros_like(vr)
    e = 0.0; s = -VSAT
    dec = np.exp(-DT/TAU); chg = 1-np.exp(-DT/TAUC)
    for i, x in enumerate(vr):
        e *= dec
        if x > e: e += (x-e)*chg
        env[i] = e
        if s < 0 and e > VTH: s = VSAT
        elif s > 0 and e < VTL: s = -VSAT
        vo[i] = s
    return vg, vhw, vr, env, vo

def pulse(t, f, A, w=0.075, pol=1):
    return pol*A*(((t % (1/f)) < w).astype(float))

def stats(t, vo, nbeats_expected, f):
    hi = vo > 0
    rises = np.where(np.diff(hi.astype(int)) == 1)[0]
    falls = np.where(np.diff(hi.astype(int)) == -1)[0]
    widths = []
    for r in rises:
        fl = falls[falls > r]
        if len(fl): widths.append((fl[0]-r)*DT)
    stuck = hi[-1] and len(falls) < len(rises)-1
    return len(rises), (np.median(widths)*1e3 if widths else float("nan")), stuck

# ---------- pulse test-case table ----------
rows = []
T_SIM = 10.0
t = np.arange(0, T_SIM, DT)
for bpm in [40, 72, 120, 150, 180]:
    f = bpm/60
    for A in [0.5, 1.0, 2.0]:
        for pol in [1, -1]:
            vin = pulse(t+0.2, f, A, pol=pol)   # +0.2 s: first pulse not at t=0
            vg, vhw, vr, env, vo = chain(vin)
            n_exp = int(np.sum(np.diff((((t+0.2) % (1/f)) < 0.075).astype(int)) == 1))
            n, wms, stuck = stats(t, vo, n_exp, f)
            # envelope minimum just before each beat (last complete beat)
            starts = np.where(np.diff((((t+0.2) % (1/f)) < 0.075).astype(int)) == 1)[0]
            pre = [env[s0-1] for s0 in starts[1:]]
            rows.append(dict(bpm=bpm, Vpp=A, pol=pol, expected=n_exp, detected=n,
                             width_ms=round(wms, 0), peak=round(env.max(), 2),
                             env_min_before=round(min(pre), 2) if pre else None))

# ---------- theoretical vs simulated discharge time ----------
theo = []
for A in [0.5, 1.0, 2.0]:
    Vp = min(G*A, VSAT)
    td = TAU*np.log(Vp/VTL)
    theo.append(dict(Vpp=A, Vp=round(Vp, 2), t_desc_ms=round(td*1e3, 0),
                     width_theo_ms=round((0.075+td)*1e3, 0),
                     fmax_bpm=round(60/(0.075+td), 0)))

# ---------- sinusoid ----------
sine_rows = []
for fs_ in [0.5, 1.0, 1.5]:
    for A in [0.5, 1.0, 2.0]:
        ts = np.arange(0, 8, DT)
        vin = (A/2)*np.sin(2*np.pi*fs_*ts)
        *_, env, vo = chain(vin)
        n, wms, stuck = stats(ts, vo, 0, fs_)
        sine_rows.append(dict(f=fs_, Vpp=A, cycles=int(8*fs_), pulses=n,
                              pulses_per_cycle=round(n/(8*fs_), 2),
                              env_min=round(env[int(2/DT):].min(), 2)))

# ---------- rectifier tolerances (Monte Carlo + worst case) ----------
def rect_gains(R2, R3, R4, R5, R6):
    gpos = (R6/R5)*(R3/R2) - R6/R4
    gneg = R6/R4
    return gpos, gneg
worst = {}
for tol in [0.01, 0.05]:
    best = 0
    for s2 in [-1, 1]:
        for s3 in [-1, 1]:
            for s4 in [-1, 1]:
                for s5 in [-1, 1]:
                    for s6 in [-1, 1]:
                        gp, gn = rect_gains(10e3*(1+s2*tol), 10e3*(1+s3*tol), 20e3*(1+s4*tol),
                                            10e3*(1+s5*tol), 20e3*(1+s6*tol))
                        best = max(best, abs(gp-gn)/gn)
    rng = np.random.default_rng(1)
    sd = tol/3
    R = [10e3, 10e3, 20e3, 10e3, 20e3]
    m = np.array([R]*20000)*(1+rng.normal(0, sd, (20000, 5)))
    gp, gn = rect_gains(*m.T)
    worst[tol] = dict(worst=round(best*100, 1), p95=round(np.percentile(np.abs(gp-gn)/gn, 95)*100, 2))

res = dict(G=G, TAU=TAU, VTH=VTH, VTL=VTL, VREF=VREF, k=k, theo=theo, rows=rows,
           sine=sine_rows, tol=worst)
json.dump(res, open(OUT + "/results.json", "w"), indent=1)

def thr(ax):
    ax.axhline(VTH, color=INK2, ls="--", lw=1)
    ax.axhline(VTL, color=INK2, ls=":", lw=1.2)
    ax.text(ax.get_xlim()[1], VTH, f" VTH {VTH:.2f} V", va="center", fontsize=8.5, color=INK2)
    ax.text(ax.get_xlim()[1], VTL, f" VTL {VTL:.2f} V", va="center", fontsize=8.5, color=INK2)

def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

# F1 rectifier transfer characteristic
vgx = np.linspace(-10, 10, 401)
_, vhw, vr, *_ = chain(vgx/G) if False else (None, np.where(vgx > 0, -vgx, 0), np.abs(vgx))
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
ax[0].plot(vgx, np.where(vgx > 0, -vgx, 0), color=ORANGE, label="Vhw (TP2)")
ax[0].plot(vgx, np.abs(vgx), color=BLUE, label="Vr (TP3)")
ax[0].set_xlabel("Vg (V)"); ax[0].set_ylabel("Output (V)")
ax[0].set_title("Transfer curve", fontsize=10.5, color=INK, loc="left")
ax[0].legend(loc="lower right")
ts = np.arange(0, 2, DT); vin = 0.5*np.sin(2*np.pi*1*ts)
vg, vhw2, vr2, *_ = chain(vin)
ax[1].plot(ts, vg, color=MUTED, lw=1.2, label="Vg (TP1)")
ax[1].plot(ts, vhw2, color=ORANGE, label="Vhw (TP2)")
ax[1].plot(ts, vr2, color=BLUE, label="Vr (TP3)")
ax[1].set_xlabel("Time (s)"); ax[1].set_ylabel("Voltage (V)")
ax[1].set_title("Sinusoid 1 Vpp, 1 Hz", fontsize=10.5, color=INK, loc="left")
ax[1].set_ylim(-3.7, 2.9); ax[1].legend(loc="lower center", ncol=3, fontsize=8.5)
fig.tight_layout(); save(fig, "f1_rectificador")

# F2 waveforms, 1 Vpp pulse at 72 bpm
ts = np.arange(0, 3.4, DT); f = 1.2
vin = pulse(ts+0.6, f, 1.0)
vg, vhw, vr, env, vo = chain(vin)
fig, ax = plt.subplots(5, 1, figsize=(9, 8.6), sharex=True)
for a, y, c, lab in zip(ax, [vin, vg, vr, env, vo], [INK2, BLUE, BLUE, ORANGE, AQUA],
                        ["TP0 input (V)", "TP1 Vg (V)", "TP3 Vr (V)", "TP4 Venv (V)", "TP5 Vo (V)"]):
    a.plot(ts, y, color=c); a.set_ylabel(lab)
ax[2].plot(ts, vhw, color=ORANGE, lw=1.2, label="TP2 Vhw"); ax[2].legend(loc="lower right", fontsize=8.5)
ax[3].set_xlim(0, 3.4); thr(ax[3])
ax[-1].set_xlabel("Time (s)")
ax[0].set_title("1 Vpp pulse at 72 bpm (1.2 Hz), 75 ms width", fontsize=11, color=INK, loc="left")
fig.tight_layout(); save(fig, "f2_pulso_72lpm")

# F3 envelope discharge
fig, ax = plt.subplots(figsize=(9, 3.6))
tt = np.linspace(0, 0.6, 600)
for A, c in zip([0.5, 1.0, 2.0], [AQUA, BLUE, VIO]):
    Vp = G*A; v = Vp*np.exp(-tt/TAU); td = TAU*np.log(Vp/VTL)
    ax.plot(tt*1e3, v, color=c, label=f"{A} Vpp (peak {Vp:.2f} V)")
    ax.plot([td*1e3], [VTL], "o", color=c, ms=6, mec="white", mew=1.5)
    ax.annotate(f"{td*1e3:.0f} ms", (td*1e3, VTL), textcoords="offset points", xytext=(4, 8),
                fontsize=8.5, color=INK2)
ax.set_xlim(0, 600); ax.set_ylim(0, 10.5); thr(ax)
ax.set_xlabel("Time since end of pulse (ms)"); ax.set_ylabel("Venv (V)")
ax.set_title("C1 discharge with τ = 100 ms down to VTL", fontsize=11, color=INK, loc="left")
ax.legend(loc="upper right"); fig.tight_layout(); save(fig, "f3_descarga")

# F4 maximum rate vs amplitude and vs tau
fig, ax = plt.subplots(figsize=(9, 3.6))
A = np.linspace(0.25, 2.1, 200)
for tau_, c, lab in [(0.068, AQUA, "τ = 68 ms"), (0.1, BLUE, "τ = 100 ms (design)"), (0.15, ORANGE, "τ = 150 ms")]:
    Vp = np.minimum(G*A, VSAT)
    fmax = 60/(0.075 + tau_*np.log(np.maximum(Vp/VTL, 1.0001)))
    ax.plot(A, fmax, color=c, label=lab)
ax.axhline(150, color=INK2, ls="--", lw=1); ax.text(2.1, 150, " 150 bpm", va="center", fontsize=8.5, color=INK2)
ax.axvline(VTH/G, color=MUTED, ls=":", lw=1.2)
ax.text(VTH/G, 30, " detection\n threshold", fontsize=8.5, color=INK2)
ax.set_ylim(0, 420); ax.set_xlim(0.2, 2.1)
ax.set_xlabel("Input amplitude (Vpp)"); ax.set_ylabel("Maximum rate (bpm)")
ax.set_title("Maximum heart rate with timely release: 60 / (w + τ·ln(Vp/VTL))", fontsize=10.5, color=INK, loc="left")
ax.legend(loc="upper right"); fig.tight_layout(); save(fig, "f4_fmax")

# F5 hysteresis loop
fig, ax = plt.subplots(figsize=(6, 3.6))
x = np.concatenate([np.linspace(0, 2, 400), np.linspace(2, 0, 400)])
s = -VSAT; y = []
for xi in x:
    if s < 0 and xi > VTH: s = VSAT
    elif s > 0 and xi < VTL: s = -VSAT
    y.append(s)
ax.plot(x, y, color=BLUE)
for xv, lab, ha, dx in [(VTL, f"VTL = {VTL:.2f} V", "right", -0.08), (VTH, f"VTH = {VTH:.2f} V", "left", 0.08)]:
    ax.axvline(xv, color=MUTED, ls=":", lw=1.2); ax.text(xv+dx, 0.6, lab, fontsize=8.5, color=INK2, ha=ha)
ax.annotate("", xy=(VTH+0.06, 6), xytext=(VTH+0.06, -6), arrowprops=dict(arrowstyle="->", color=INK2, lw=1.4))
ax.annotate("", xy=(VTL-0.06, -6), xytext=(VTL-0.06, 6), arrowprops=dict(arrowstyle="->", color=INK2, lw=1.4)); ax.text(1.25, VSAT-2.2, "Venv rising: trips at VTH", fontsize=8.5, color=INK2); ax.text(0.03, -VSAT+1.4, "Venv falling: releases at VTL", fontsize=8.5, color=INK2)
ax.set_xlabel("Venv (TP4) (V)"); ax.set_ylabel("Vo (TP5) (V)")
ax.set_title("Schmitt trigger (U6) hysteresis loop", fontsize=11, color=INK, loc="left")
fig.tight_layout(); save(fig, "f5_histeresis")

# F6 sinusoid 0.5 Hz and 1.5 Hz, 2 Vpp
fig, ax = plt.subplots(3, 2, figsize=(9, 6), sharex="col")
for j, fs_ in enumerate([0.5, 1.5]):
    ts = np.arange(0, 4, DT); vin = 1.0*np.sin(2*np.pi*fs_*ts)
    vg, vhw, vr, env, vo = chain(vin)
    ax[0, j].plot(ts, vg, color=BLUE); ax[0, j].set_title(f"Sinusoid 2 Vpp, {fs_} Hz", fontsize=10.5, color=INK, loc="left")
    ax[1, j].plot(ts, vr, color=MUTED, lw=1.2, label="Vr"); ax[1, j].plot(ts, env, color=ORANGE, label="Venv")
    ax[1, j].set_xlim(0, 4); thr(ax[1, j]) if j == 1 else (ax[1, j].axhline(VTH, color=INK2, ls="--", lw=1), ax[1, j].axhline(VTL, color=INK2, ls=":", lw=1.2))
    ax[2, j].plot(ts, vo, color=AQUA); ax[2, j].set_xlabel("Time (s)")
ax[0, 0].set_ylabel("TP1 Vg (V)"); ax[1, 0].set_ylabel("TP3/TP4 (V)"); ax[2, 0].set_ylabel("TP5 Vo (V)")
ax[1, 0].legend(loc="upper right", fontsize=8.5)
fig.tight_layout(); save(fig, "f6_sinusoide")

# F7 synthetic PWLIN heartbeat, both polarities
pts = np.array([[0, 0], [0.300, 0], [0.320, -0.10], [0.345, 1.00], [0.370, -0.25], [0.390, 0],
                [0.500, 0], [0.580, 0.20], [0.660, 0], [0.833, 0]])
def pwl(t):
    return np.interp(t % 0.8333, pts[:, 0], pts[:, 1])
fig, ax = plt.subplots(3, 2, figsize=(9, 6), sharex=True)
for j, pol in enumerate([1, -1]):
    ts = np.arange(0, 2.5, DT); vin = pol*pwl(ts)
    vg, vhw, vr, env, vo = chain(vin)
    ax[0, j].plot(ts, vin, color=INK2)
    ax[0, j].set_title("Positive QRS" if pol > 0 else "Inverted QRS", fontsize=10.5, color=INK, loc="left")
    ax[1, j].plot(ts, vr, color=MUTED, lw=1.2, label="Vr"); ax[1, j].plot(ts, env, color=ORANGE, label="Venv")
    ax[1, j].axhline(VTH, color=INK2, ls="--", lw=1); ax[1, j].axhline(VTL, color=INK2, ls=":", lw=1.2)
    ax[2, j].plot(ts, vo, color=AQUA); ax[2, j].set_xlabel("Time (s)")
ax[0, 0].set_ylabel("TP0 (V)"); ax[1, 0].set_ylabel("TP3/TP4 (V)"); ax[2, 0].set_ylabel("TP5 Vo (V)")
ax[1, 0].legend(loc="upper right", fontsize=8.5)
fig.tight_layout(); save(fig, "f7_qrs_pwlin")

print(json.dumps(dict(G=G, VTH=VTH, VTL=VTL, theo=theo, tol=worst), indent=1))
for r in rows: print(r)
for r in sine_rows: print(r)
