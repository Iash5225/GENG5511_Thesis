from enum import IntEnum
import numpy as np
from dataclasses import dataclass, field
import matplotlib.pyplot as plt

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

    def plot_history(self):
        plt.figure(figsize=(12, 8))
        for i, metric in enumerate(Metric):
            plt.subplot(4, 4, i + 1)
            plt.plot(self.history[metric], marker='o', linestyle='-')
            plt.title(_METRIC_LABELS[metric])
            plt.xlabel('Iteration')
            plt.ylabel('Deviation')
            plt.grid(True)
        plt.tight_layout()
        plt.show()
