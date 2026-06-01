"""
visualization.py
----------------
Plotting utilities for SAR change detection outputs.

Author: Axel Franke
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.patches import Patch
from pathlib import Path
from typing import Optional, Union


COLORS = {
    "accumulation": "#2ecc71",
    "drawdown": "#e74c3c",
    "stable": "#95a5a6",
    "backscatter": "#2c3e50",
    "threshold": "#e67e22",
    "ground_truth": "#8e44ad",
}


def plot_backscatter_timeseries(
    df: pd.DataFrame,
    title: str = "SAR Backscatter Time Series — AOI",
    output_path: Optional[Union[str, Path]] = None,
    figsize: tuple = (14, 5),
) -> plt.Figure:
    """Plot mean AOI backscatter over time with significant change events highlighted.

    Parameters
    ----------
    df : pd.DataFrame
        Output from change_detection.build_time_series().
        Required columns: date, backscatter_db, significant, direction.
    title : str
        Plot title.
    output_path : str or Path, optional
        If provided, save figure to this path.

    Returns
    -------
    plt.Figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    ax.plot(df["date"], df["backscatter_db"], color=COLORS["backscatter"],
            linewidth=1.8, marker="o", markersize=5, zorder=3, label="Mean σ⁰ (dB)")

    # Highlight significant change events
    sig = df[df["significant"] == True]
    for _, row in sig.iterrows():
        color = COLORS.get(row["direction"], COLORS["stable"])
        ax.axvline(row["date"], color=color, alpha=0.4, linewidth=6, zorder=1)

    ax.set_xlabel("Acquisition Date", fontsize=11)
    ax.set_ylabel("Mean Backscatter σ⁰ (dB)", fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    plt.xticks(rotation=30, ha="right")

    legend_elements = [
        Patch(facecolor=COLORS["accumulation"], alpha=0.5, label="Significant accumulation"),
        Patch(facecolor=COLORS["drawdown"], alpha=0.5, label="Significant drawdown"),
    ]
    ax.legend(handles=[ax.lines[0]] + legend_elements, fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")

    return fig


def plot_change_events(
    df: pd.DataFrame,
    ground_truth_df: Optional[pd.DataFrame] = None,
    title: str = "Log-Ratio Change Detection Results",
    output_path: Optional[Union[str, Path]] = None,
    figsize: tuple = (14, 6),
) -> plt.Figure:
    """Bar chart of log-ratio change per image pair with optional ground truth overlay.

    Parameters
    ----------
    df : pd.DataFrame
        Output from build_time_series(). Required: date, log_ratio, significant.
    ground_truth_df : pd.DataFrame, optional
        Public operational data for validation. Required columns: date, value.
        Will be plotted on a secondary axis.
    title : str
        Plot title.
    output_path : str or Path, optional
        Save path.

    Returns
    -------
    plt.Figure
    """
    fig, ax1 = plt.subplots(figsize=figsize)

    plot_df = df.dropna(subset=["log_ratio"])
    colors = [
        COLORS["accumulation"] if (r > 0 and s) else
        COLORS["drawdown"] if (r < 0 and s) else
        COLORS["stable"]
        for r, s in zip(plot_df["log_ratio"], plot_df["significant"])
    ]

    ax1.bar(plot_df["date"], plot_df["log_ratio"], color=colors,
            width=12, alpha=0.85, zorder=2)

    threshold = plot_df["threshold_used_db"].iloc[0] if "threshold_used_db" in plot_df.columns else 2.0
    ax1.axhline(threshold, color=COLORS["threshold"], linestyle="--",
                linewidth=1.2, label=f"+{threshold:.1f} dB threshold", zorder=3)
    ax1.axhline(-threshold, color=COLORS["threshold"], linestyle="--",
                linewidth=1.2, label=f"-{threshold:.1f} dB threshold", zorder=3)
    ax1.axhline(0, color="black", linewidth=0.8, alpha=0.5)

    ax1.set_xlabel("Image Pair Date", fontsize=11)
    ax1.set_ylabel("Log-Ratio Change (dB)", fontsize=11, color=COLORS["backscatter"])
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    plt.xticks(rotation=30, ha="right")

    if ground_truth_df is not None:
        ax2 = ax1.twinx()
        ax2.plot(ground_truth_df["date"], ground_truth_df["value"],
                 color=COLORS["ground_truth"], linewidth=2, marker="D",
                 markersize=6, label="BHP Quarterly Production (kt)", zorder=4)
        ax2.set_ylabel("Reported Production (kt)", fontsize=11, color=COLORS["ground_truth"])
        ax2.tick_params(axis="y", labelcolor=COLORS["ground_truth"])
        lines2, labels2 = ax2.get_legend_handles_labels()
    else:
        lines2, labels2 = [], []

    legend_elements = [
        Patch(facecolor=COLORS["accumulation"], alpha=0.85, label="Significant accumulation"),
        Patch(facecolor=COLORS["drawdown"], alpha=0.85, label="Significant drawdown"),
        Patch(facecolor=COLORS["stable"], alpha=0.6, label="No significant change"),
    ]
    ax1.legend(handles=legend_elements + lines2,
               labels=[e.get_label() for e in legend_elements] + labels2,
               fontsize=9, loc="upper left")

    ax1.set_title(title, fontsize=13, fontweight="bold")
    ax1.grid(True, alpha=0.2, axis="y")
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")

    return fig


def plot_volumetric_summary(
    results: list,
    output_path: Optional[Union[str, Path]] = None,
    figsize: tuple = (12, 5),
) -> plt.Figure:
    """Plot estimated mass change per image pair with error bounds.

    Parameters
    ----------
    results : list of dict
        List of outputs from change_detection.backscatter_to_volume().
        Each dict should also contain a 'date' key.
    output_path : str or Path, optional
        Save path.

    Returns
    -------
    plt.Figure
    """
    dates = [r["date"] for r in results]
    masses = [r["mass_change_t"] for r in results]
    errors = [r["error_bounds_t"] for r in results]

    fig, ax = plt.subplots(figsize=figsize)

    colors = [COLORS["accumulation"] if m >= 0 else COLORS["drawdown"] for m in masses]
    ax.bar(dates, masses, color=colors, width=12, alpha=0.8, zorder=2)
    ax.errorbar(dates, masses, yerr=errors, fmt="none",
                color="black", capsize=4, linewidth=1.2, zorder=3)
    ax.axhline(0, color="black", linewidth=0.8)

    ax.set_xlabel("Image Pair Date", fontsize=11)
    ax.set_ylabel("Estimated Mass Change (metric tonnes)", fontsize=11)
    ax.set_title("Volumetric Change Estimates with Uncertainty Bounds", fontsize=13, fontweight="bold")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    plt.xticks(rotation=30, ha="right")

    note = "Note: Error bars reflect combined uncertainty from spatial resolution,\nbulk density assumption (±10%), and sensitivity parameter (±30%)"
    ax.text(0.01, 0.02, note, transform=ax.transAxes, fontsize=8,
            color="gray", verticalalignment="bottom")

    ax.grid(True, alpha=0.2, axis="y")
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")

    return fig
