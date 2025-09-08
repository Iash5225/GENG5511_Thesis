# pip install seaborn pandas matplotlib

from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ---- Paste your ST dict ----
ST = {
    "Vm":       {"v00": 0.872, "g0": 0.162, "q0": 0.105, "a1": 0.095, "cc": 0.040, "aa": 0.039},
    "KappaT":   {"cc": 0.840, "Th0": 0.617, "q0": 0.556, "a1": 0.470, "g0": 0.432, "bb": 0.259},
    "KappaS":   {"v00": 1.408, "Th0": 1.038, "a1": 0.920, "cc": 0.710, "g0": 0.640, "q0": 0.520},
    "Alpha":    {"Th0": 721.583, "q0": 473.591, "a1": 8.140, "bb": 2.935, "a3": 2.756, "g0": 1.001},
    "cp":       {"Th0": 1249575.008, "q0": 16968.244, "bb": 2977.267, "a1": 1667.928, "a3": 1.797, "cc": 1.658},
    "cv":       {"cc": 16.115, "Th0": 2.248, "v00": 1.029, "g0": 1.006, "a1": 1.004, "aa": 1.001},
    "Gruneisen": {"q0": 3.884, "aa": 1.142, "Th0": 1.126, "cc": 0.876, "g0": 0.858, "bb": 0.708},
    "U":        {"cc": 16.645, "Th0": 2.255, "v00": 1.028, "g0": 1.006, "a1": 1.004, "aa": 1.001},
    "S":        {"cc": 15.838, "Th0": 2.253, "v00": 1.020, "g0": 1.006, "a1": 1.004, "aa": 1.001},
    "A":        {"cc": 15.074, "Th0": 2.250, "v00": 1.013, "g0": 1.005, "a1": 1.003, "aa": 1.001},
    "H":        {"cc": 16.645, "Th0": 2.255, "v00": 1.028, "g0": 1.006, "a1": 1.004, "aa": 1.001},
    "G":        {"cc": 15.074, "Th0": 2.250, "v00": 1.013, "g0": 1.005, "a1": 1.003, "aa": 1.001},
}

# ---- Tidy DataFrame ----
rows = []
for prop, d in ST.items():
    for param, val in d.items():
        rows.append({"Property": prop, "Parameter": param, "ST": float(val)})
df = pd.DataFrame(rows)

# Groups
elastic = ["v00", "a1", "a2", "a3"]
vibrational = ["Th0", "g0", "q0"]
anharmonic = ["aa", "bb", "cc"]


def group_param(p):
    if p in elastic:
        return "Elastic"
    if p in vibrational:
        return "Vibrational"
    if p in anharmonic:
        return "Anharmonic"
    return "Other"


df["Group"] = df["Parameter"].map(group_param)

# Row-normalised sizes
df["ST_norm"] = df.groupby("Property")["ST"].transform(
    lambda x: x / (x.max() if x.max() > 0 else 1.0)
)

# Order columns by groups
param_order = [p for p in (elastic + vibrational +
                           anharmonic) if p in df["Parameter"].unique()]
prop_order = list(ST.keys())

# LaTeX labels with subscripts
latex_map = {
    "v00": r"$v_{00}$", "a1": r"$a_{1}$", "a2": r"$a_{2}$", "a3": r"$a_{3}$",
    "Th0": r"$\Theta_{0}$", "g0": r"$g_{0}$", "q0": r"$q_{0}$",
    "aa": r"$a_{a}$", "bb": r"$b_{b}$", "cc": r"$c_{c}$",
}
param_labels = [latex_map.get(p, p) for p in param_order]

# Palette
palette = {"Elastic": "#4C78A8", "Vibrational": "#59A14F",
           "Anharmonic": "#E15759", "Other": "#B07AA1"}

# --- Figure with legend column ---
sns.set_theme(style="whitegrid", context="talk")
fig = plt.figure(figsize=(14, 7))
gs = GridSpec(1, 2, width_ratios=[4, 1.1], wspace=0.05)

ax = fig.add_subplot(gs[0, 0])
ax_leg = fig.add_subplot(gs[0, 1])
ax_leg.axis('off')  # this axis is only for legends

# Scatter on main axis
sns.scatterplot(
    data=df,
    x=pd.Categorical(df["Parameter"], categories=param_order, ordered=True),
    y=pd.Categorical(df["Property"], categories=prop_order, ordered=True),
    hue="Group",
    size="ST_norm", sizes=(70, 1100),
    palette=palette,
    marker="o",
    edgecolor="black", linewidth=0.5,
    alpha=0.9,
    ax=ax
)

ax.set_xlabel("Equation of State Parameters", labelpad=10)
ax.set_ylabel("Thermodynamic Properties", labelpad=10)
ax.set_xticks(range(len(param_order)))
ax.set_xticklabels(param_labels, rotation=45, ha="right")
ax.set_yticks(range(len(prop_order)))
ax.set_yticklabels(prop_order)
ax.grid(axis="y", linestyle=":", alpha=0.35)
ax.grid(axis="x", visible=False)

# Remove auto legend from main axis
if ax.legend_:
    ax.legend_.remove()

# --- Build legends in the right panel with the same width ---
# 1) Size legend (horizontal row)
size_vals = [0.3, 0.6, 1.0]
size_handles = [ax.scatter([], [], s=70+(1100-70)*v,
                           facecolor="#888888", edgecolor="black", linewidth=0.5)
                for v in size_vals]
size_labels = [f"{v:.1f}" for v in size_vals]
leg_size = ax_leg.legend(size_handles, size_labels, title="Relative Sensitivity",
                         loc="upper left", bbox_to_anchor=(0.0, 1.0),
                         frameon=True, ncol=3, columnspacing=1.2, handlelength=1.5, handletextpad=0.6)

# 2) Group legend below
group_handles = [Line2D([0], [0], marker='o', color='none',
                        markerfacecolor=palette[g], markeredgecolor="black",
                        markersize=10, label=g)
                 for g in ["Elastic", "Vibrational", "Anharmonic"] if g in df["Group"].unique()]
leg_group = ax_leg.legend(handles=group_handles, title="Parameter Group",
                          loc="upper left", bbox_to_anchor=(0.0, 0.70),
                          frameon=True, handletextpad=0.6)

# Keep both legends
ax_leg.add_artist(leg_size)
ax_leg.add_artist(leg_group)

plt.tight_layout()

# out_path = Path("/mnt/data/sensitivity_bubble_matrix_thesis_legends_panel.png")
# plt.savefig(out_path, dpi=300, bbox_inches="tight")
plt.show()

#
