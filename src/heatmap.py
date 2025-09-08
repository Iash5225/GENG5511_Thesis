# pip install seaborn pandas matplotlib

from matplotlib.lines import Line2D
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

# Order columns: by group (Elastic → Vibrational → Anharmonic), then mean ST_norm desc within group


def order_within(names):
    sub = df[df["Parameter"].isin(names)].groupby(
        "Parameter")["ST_norm"].mean().sort_values(ascending=False).index.tolist()
    # keep only those that actually appear
    return [p for p in names if p in sub] + [p for p in sub if p not in names]


elastic_cols = [p for p in elastic if p in df["Parameter"].unique()]
vibrational_cols = [p for p in vibrational if p in df["Parameter"].unique()]
anharmonic_cols = [p for p in anharmonic if p in df["Parameter"].unique()]

elastic_cols = order_within(elastic_cols)
vibrational_cols = order_within(vibrational_cols)
anharmonic_cols = order_within(anharmonic_cols)

param_order = elastic_cols + vibrational_cols + anharmonic_cols
prop_order = list(ST.keys())  # keep user order

# Palette (muted thesis-style)
palette = {
    "Elastic":     "#4C78A8",  # muted blue
    "Vibrational": "#59A14F",  # muted green
    "Anharmonic":  "#E15759",  # muted red
    "Other":       "#B07AA1",  # muted purple (unused here)
}

# Figure + style
sns.set_theme(style="whitegrid", context="talk")  # larger fonts
fig, ax = plt.subplots(
    figsize=(1.2 + 0.55*len(param_order), 0.8 + 0.55*len(prop_order)))

# Subtle vertical band shading per group
x_ticks = np.arange(len(param_order))


def span_for(cols):
    if not cols:
        return None
    start = param_order.index(cols[0]) - 0.5
    end = param_order.index(cols[-1]) + 0.5
    return start, end


for cols, color in [(elastic_cols, "#4C78A810"), (vibrational_cols, "#59A14F10"), (anharmonic_cols, "#E1575910")]:
    rng = span_for(cols)
    if rng:
        ax.axvspan(rng[0], rng[1], color=color, zorder=0)

# Scatter (circles, black edges)
sns.scatterplot(
    data=df,
    x=pd.Categorical(df["Parameter"], categories=param_order, ordered=True),
    y=pd.Categorical(df["Property"], categories=prop_order, ordered=True),
    hue="Group",
    size="ST_norm", sizes=(80, 1400),
    palette=palette,
    marker="o",
    edgecolor="black", linewidth=0.5,
    alpha=0.9,
    ax=ax
)

# Axes cleanup
ax.set_xlabel("Equation of State Parameters", labelpad=10)
ax.set_ylabel("Thermodynamic Properties", labelpad=10)
ax.set_title(
    "Global Sensitivity (Sobol ST) — row-normalised per property", pad=14)
ax.set_xticklabels(param_order, rotation=45, ha="right")
# Keep horizontal grid, lighten it; remove vertical grid lines
ax.grid(axis="y", linestyle=":", alpha=0.4)
ax.grid(axis="x", visible=False)

# Limits & ticks
ax.set_xlim(-0.5, len(param_order)-0.5)
ax.set_ylim(-0.5, len(prop_order)-0.5)
ax.set_xticks(range(len(param_order)))
ax.set_yticks(range(len(prop_order)))
ax.set_yticklabels(prop_order)

# Annotate top bubble per row with raw ST value
for i, prop in enumerate(prop_order):
    sub = df[df["Property"] == prop]
    if sub.empty:
        continue
    # get param with max ST (raw, not normalized)
    idx = sub["ST"].idxmax()
    if pd.isna(idx):
        continue
    row = sub.loc[idx]
    x = param_order.index(row["Parameter"])
    y = i
    # format raw ST: 3 sig figs or scientific for big numbers
    val = row["ST"]
    if val >= 1e4:
        label = f"{val:.2e}"
    else:
        label = f"{val:.3g}"
    ax.text(x, y, label, ha="center", va="center",
            fontsize=10, color="white", weight="bold")

# Legends: size + group
# Remove seaborn's auto legend and rebuild cleaner ones
if ax.legend_:
    ax.legend_.remove()

# Size legend (relative sensitivity)
size_vals = [0.3, 0.6, 1.0]
size_handles = [plt.scatter([], [], s=80 + (1400-80)*v, edgecolor="black",
                            linewidth=0.5, facecolor="#777777", alpha=0.8) for v in size_vals]
size_labels = [f"{v:.1f}" for v in size_vals]
leg1 = ax.legend(size_handles, size_labels,
                 title="Relative sensitivity\n(ST, row-normalised)", loc="upper left", frameon=True)
ax.add_artist(leg1)

# Group legend
group_handles = [Line2D([0], [0], marker='o', color='none', label=lab,
                        markerfacecolor=palette[lab], markeredgecolor="black",
                        markersize=10) for lab in ["Elastic", "Vibrational", "Anharmonic"] if lab in df["Group"].unique()]
leg2 = ax.legend(handles=group_handles, title="Parameter group",
                 loc="lower right", frameon=True)

plt.tight_layout()

# Save figure
#out_path = Path("C:\Users\iashb\OneDrive - The University of Western Australia\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\img\outputsensitivity_bubble_matrix_thesis.png")
#plt.savefig(out_path, dpi=300, bbox_inches="tight")
plt.show()

#
