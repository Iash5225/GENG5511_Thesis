from enum import IntEnum
import numpy as np
from dataclasses import dataclass, field
import matplotlib.pyplot as plt
import math
# from math import ceil
class Metric(IntEnum):
    VM_SUB = 0
    VM_MELT = 1
    VM_HIGHP = 2
    CP_SUB = 3
    ALPHA_SUB = 4
    BETAT_SUB = 5
    BETAS_SUB = 6
    H_SOLID_SUB = 7
    H_SOLID_MELT = 8
    P_SUB = 9
    P_MELT = 10
    GAMMA_T = 11
    TOTAL = 12


_METRIC_LABELS = {
    Metric.VM_SUB:       "Vm_sub",
    Metric.VM_MELT:      "Vm_melt",
    Metric.VM_HIGHP:     "Vm_highp",
    Metric.CP_SUB:       "cp_sub",
    Metric.ALPHA_SUB:    "alpha_sub",
    Metric.BETAT_SUB:    "kappaT_sub",
    Metric.BETAS_SUB:    "kappaS_sub",
    Metric.H_SOLID_SUB:  "H_solid_sub",
    Metric.H_SOLID_MELT: "H_solid_melt",
    Metric.P_SUB:        "P_sub",
    Metric.P_MELT:       "P_melt",
    Metric.GAMMA_T:      "Gamma_T",
    Metric.TOTAL:        "Total cost",
}
@dataclass
class DeviationRecorder:
    metrics: np.ndarray = field(default_factory=lambda: np.zeros(len(Metric)))
    counts: np.ndarray = field(default_factory=lambda: np.zeros(len(Metric), dtype=int))
    history: list = field(default_factory=lambda: [
                          [] for _ in range(len(Metric))])

    def record(self, metric: Metric, deviation: float):
        self.metrics[metric] += float(deviation)
        self.counts[metric] += 1
        self.history[metric].append(float(deviation))

    def total_deviation(self) -> float:
        return np.sum(self.metrics)
    def reset(self):
        self.metrics.fill(0)
        self.counts.fill(0)
        self.history = [
            [] for _ in range(len(Metric))
        ]

    def plot_history(self, ncols: int = 8, col_w: float = 3.2, row_h: float = 2.6):
        """
        Make a horizontally wide subplot grid for deviation histories.

        Parameters
        ----------
        ncols : int
            Number of subplot columns (increase for wider layout).
        col_w : float
            Width per column in inches.
        row_h : float
            Height per row in inches.
        """
        n_metrics = len(Metric)
        nrows = int(math.ceil(n_metrics / ncols))

        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(col_w * ncols, row_h * nrows),
            sharex=False, sharey=False
        )

        # Normalize axes to 2D array
        axes = np.atleast_2d(axes)

        for i, metric in enumerate(Metric):
            r, c = divmod(i, ncols)
            ax = axes[r, c]
            y = np.asarray(self.history[metric], dtype=float)
            x = np.arange(1, y.size + 1)
            if y.size:
                ax.plot(x, y, marker='o', ms=3, lw=1.5)
            else:
                ax.text(0.5, 0.5, "no data", ha="center", va="center", fontsize=9)
            ax.set_title(_METRIC_LABELS[metric], fontsize=10)
            ax.set_xlabel("Evaluation", fontsize=9)
            ax.set_ylabel("Deviation", fontsize=9)
            ax.grid(True, alpha=0.25)
            ax.tick_params(labelsize=8)

        # Hide any unused axes
        for j in range(n_metrics, nrows * ncols):
            r, c = divmod(j, ncols)
            axes[r, c].axis("off")

        fig.tight_layout()
        plt.show()
        return fig, axes
    
    def _ensure_store(self):
        if not hasattr(self, "_pointwise_rows"):
            self._pointwise_rows = []  # list of dict rows

    
