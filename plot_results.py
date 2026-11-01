"""Create figures used by nozzle_fvm_report.tex.

Run from the project root with:
    pvpython plot_results.py
"""

import csv
import os

import matplotlib.pyplot as plt


OUT = "results/visualizations"
os.makedirs(OUT, exist_ok=True)


def read_csv(name):
    with open(name, newline="") as stream:
        rows = list(csv.DictReader(stream))
    return rows


def values(name, column):
    return [float(row[column]) for row in read_csv(name)]


plt.rcParams.update({"font.size": 10, "figure.dpi": 150})

# Grid-refinement profiles.
fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
for prefix, label in [("grid_100", "N = 100"), ("grid_200", "N = 200"),
                      ("grid_400", "N = 400")]:
    x = values(prefix + ".csv", "x")
    p = values(prefix + ".csv", "pressure")
    m = values(prefix + ".csv", "mach")
    axes[0].plot(x, [v / 101325.0 for v in p], label=label)
    axes[1].plot(x, m, label=label)
axes[0].set(xlabel="x [m]", ylabel=r"$p/p_0$", title="Pressure profile")
axes[1].set(xlabel="x [m]", ylabel="Mach number", title="Mach-number profile")
for ax in axes:
    ax.grid(True, alpha=0.3)
    ax.legend()
fig.savefig(os.path.join(OUT, "grid_profiles.png"), bbox_inches="tight")
plt.close(fig)

# Grid residual histories.
fig, ax = plt.subplots(figsize=(6.5, 4))
for prefix, label in [("grid_100", "N = 100"), ("grid_200", "N = 200"),
                      ("grid_400", "N = 400")]:
    step = values(prefix + "_residual.csv", "step")
    residual = values(prefix + "_residual.csv", "residual")
    ax.semilogy(step, residual, label=label)
ax.set(xlabel="Pseudo-time iteration", ylabel="Residual",
       title="Grid-refinement residual histories")
ax.grid(True, which="both", alpha=0.3)
ax.legend()
fig.savefig(os.path.join(OUT, "grid_residuals.png"), bbox_inches="tight")
plt.close(fig)

# CFL residual histories. The CFL=1 case is intentionally shown as a
# non-converged stability case.
fig, ax = plt.subplots(figsize=(6.5, 4))
for prefix, label in [("cfl_02", "CFL = 0.2"), ("cfl_05", "CFL = 0.5"),
                      ("cfl_08", "CFL = 0.8"), ("cfl_10", "CFL = 1.0")]:
    step = values(prefix + "_residual.csv", "step")
    residual = values(prefix + "_residual.csv", "residual")
    ax.semilogy(step, residual, label=label)
ax.set(xlabel="Pseudo-time iteration", ylabel="Residual",
       title="CFL-sensitivity residual histories")
ax.grid(True, which="both", alpha=0.3)
ax.legend()
fig.savefig(os.path.join(OUT, "cfl_residuals.png"), bbox_inches="tight")
plt.close(fig)

# Final profiles for converged CFL cases. CFL=1 is excluded because it did not
# produce a converged steady solution.
fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
for prefix, label in [("cfl_02", "CFL = 0.2"), ("cfl_05", "CFL = 0.5"),
                      ("cfl_08", "CFL = 0.8")]:
    x = values(prefix + ".csv", "x")
    p = values(prefix + ".csv", "pressure")
    m = values(prefix + ".csv", "mach")
    axes[0].plot(x, [v / 101325.0 for v in p], label=label)
    axes[1].plot(x, m, label=label)
axes[0].set(xlabel="x [m]", ylabel=r"$p/p_0$", title="Converged pressure profiles")
axes[1].set(xlabel="x [m]", ylabel="Mach number", title="Converged Mach profiles")
for ax in axes:
    ax.grid(True, alpha=0.3)
    ax.legend()
fig.savefig(os.path.join(OUT, "cfl_profiles.png"), bbox_inches="tight")
plt.close(fig)

# Baseline overview: geometry, pressure, and Mach number.
fig, ax1 = plt.subplots(figsize=(7, 4))
x = values("results_baseline.csv", "x")
area = values("results_baseline.csv", "area")
p = values("results_baseline.csv", "pressure")
mach = values("results_baseline.csv", "mach")
ax1.plot(x, area, color="black", label="Area")
ax1.set(xlabel="x [m]", ylabel=r"Area [m$^2$]")
ax2 = ax1.twinx()
ax2.plot(x, [v / 101325.0 for v in p], color="tab:blue", label=r"$p/p_0$")
ax2.plot(x, mach, color="tab:red", label="Mach")
ax2.set_ylabel(r"$p/p_0$ and Mach")
ax1.axvline(0.4, color="gray", linestyle="--", linewidth=0.8)
ax1.set_title("Baseline nozzle solution")
ax1.grid(True, alpha=0.3)
lines = ax1.get_lines() + ax2.get_lines()
ax1.legend(lines, [line.get_label() for line in lines], loc="best")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "baseline_overview.png"), bbox_inches="tight")
plt.close(fig)
