#!/usr/bin/env python3
from __future__ import annotations

import math
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from matplotlib.backends.backend_pdf import PdfPages

OUT = Path("dd_50um_pb_figures")
OUT.mkdir(exist_ok=True)

SOURCES = {
    "DD, 50 $\\mu$m, no Pb": Path("build/DD-50um-noPb.root"),
    "DD, 50 $\\mu$m, with Pb": Path("build/DD-50um-Pb.root"),
}

CATEGORY_TOKENS = {
    0: "e-(gamma)",
    1: "e-(ion)",
    2: "e+",
    3: "proton",
    4: "lightIon",
    5: "alpha",
    6: "Li",
    7: "C",
    8: "Be/B",
    9: "otherHeavyIon",
    10: "neutralLocal",
    11: "other",
}

CATEGORY_NAMES = {
    0: "e$^-$ (gamma-mediated)",
    1: "e$^-$ (ion-mediated)",
    2: "e$^+$",
    3: "proton",
    4: "light ion",
    5: "alpha",
    6: "lithium",
    7: "carbon",
    8: "Be/B",
    9: "other heavy ion",
    10: "neutral local",
    11: "other",
}

# Fixed category colors are used identically in both conditions so that the
# two particle-resolved depth figures can be placed side by side.
CATEGORY_COLORS = {
    0: "#1f77b4",
    1: "#ff7f0e",
    3: "#2ca02c",
    7: "#d62728",
    99: "#9467bd",
}

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "legend.fontsize": 8.1,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.55,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
    }
)

ROOTS = {label: uproot.open(path) for label, path in SOURCES.items()}


def find_object(root: uproot.ReadOnlyDirectory, title: str):
    for _, obj in root.items(recursive=True):
        if getattr(obj, "title", "") == title:
            return obj
    raise KeyError(f"Histogram not found: {title}")


def find_objects_starting(root: uproot.ReadOnlyDirectory, prefix: str):
    found = []
    for _, obj in root.items(recursive=True):
        if getattr(obj, "title", "").startswith(prefix):
            found.append(obj)
    return found


def hist_data(obj):
    return (
        np.asarray(obj.values(flow=False), dtype=float),
        np.asarray(obj.axis().edges(), dtype=float),
    )


def normalized(values: np.ndarray) -> np.ndarray:
    total = float(np.nansum(values))
    return values / total if total > 0 else np.zeros_like(values)


def finish_axes(ax: plt.Axes, *, logx=False, logy=False, legend=True) -> None:
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    ax.grid(True, which="major", alpha=0.24, linewidth=0.55)
    ax.grid(True, which="minor", alpha=0.10, linewidth=0.4)
    ax.tick_params(direction="in", top=True, right=True)
    if legend:
        ax.legend(frameon=False, handlelength=2.5)


def save_figure(fig: plt.Figure, stem: str, pdf_pages: PdfPages) -> None:
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.svg")
    fig.savefig(OUT / f"{stem}.png", dpi=400)
    pdf_pages.savefig(fig)
    plt.close(fig)


def aggregate_let(root: uproot.ReadOnlyDirectory) -> tuple[np.ndarray, np.ndarray]:
    objs = find_objects_starting(root, "LETcalc Eion-weighted spectrum, absorber 2,")
    total = None
    edges = None
    for obj in objs:
        values, obj_edges = hist_data(obj)
        if total is None:
            total = np.zeros_like(values)
            edges = obj_edges
        total += values
    if total is None or edges is None:
        raise RuntimeError("No particle-resolved LETcalc histograms were found")
    return total, edges


def particle_partition(root: uproot.ReadOnlyDirectory) -> dict[int, float]:
    tree = root["TrackLET"]
    a = tree.arrays(
        ["absorberID", "depositCategory", "eion_keV", "trackWeight"],
        library="np",
    )
    mask = a["absorberID"] == 2
    cats = a["depositCategory"][mask].astype(int)
    eion = a["eion_keV"][mask].astype(float) * a["trackWeight"][mask].astype(float)
    total = float(np.sum(eion))
    return {
        int(cat): float(np.sum(eion[cats == cat]) / total)
        for cat in np.unique(cats)
        if total > 0
    }


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(mask):
        return math.nan
    values = values[mask]
    weights = weights[mask]
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    return float(values[np.searchsorted(np.cumsum(weights), 0.5 * np.sum(weights))])


def event_summary(label: str, root: uproot.ReadOnlyDirectory) -> dict[str, float | str]:
    tree = root["EventLET"]
    cols = [
        "absorberID",
        "trackWeight",
        "eion_keV",
        "niel_keV",
        "letDcalc_keV_per_um",
        "maxLetCalc_keV_per_um",
        "fractionEionAbove10",
        "fractionEionAbove100",
        "fractionEionAbove1000",
        "meanDepth_um",
        "depthSigma_um",
    ]
    a = tree.arrays(cols, library="np")
    mask = a["absorberID"] == 2
    w = a["trackWeight"][mask].astype(float)
    eion = a["eion_keV"][mask].astype(float)
    niel = a["niel_keV"][mask].astype(float)
    hit = eion > 0
    eweight = w * eion
    total_eion = float(np.sum(eweight))
    total_weight = float(np.sum(w))
    total_dep = float(np.sum(w * (eion + niel)))

    def ew(name: str) -> float:
        x = a[name][mask].astype(float)
        return float(np.sum(eweight * x) / total_eion) if total_eion > 0 else math.nan

    return {
        "source": label.replace("$\\mu$", "u"),
        "incident_events_weighted": total_weight,
        "hit_events_weighted": float(np.sum(w[hit])),
        "hit_fraction": float(np.sum(w[hit]) / total_weight) if total_weight > 0 else math.nan,
        "mean_Eion_per_primary_keV": total_eion / total_weight if total_weight > 0 else math.nan,
        "median_Eion_per_hit_keV": weighted_median(eion[hit], w[hit]),
        "energy_weighted_median_LETcalc_keV_per_um": weighted_median(
            a["letDcalc_keV_per_um"][mask][hit].astype(float), eweight[hit]
        ),
        "mean_event_LETcalc_Eion_weighted_keV_per_um": ew("letDcalc_keV_per_um"),
        "mean_max_LETcalc_Eion_weighted_keV_per_um": ew("maxLetCalc_keV_per_um"),
        "Eion_fraction_above_10_keV_per_um": ew("fractionEionAbove10"),
        "Eion_fraction_above_100_keV_per_um": ew("fractionEionAbove100"),
        "Eion_fraction_above_1000_keV_per_um": ew("fractionEionAbove1000"),
        "Eion_weighted_mean_depth_um": ew("meanDepth_um"),
        "Eion_weighted_mean_depth_sigma_um": ew("depthSigma_um"),
        "NIEL_fraction_of_Edep": float(np.sum(w * niel) / total_dep) if total_dep > 0 else math.nan,
    }


def normalize_branch(name: str) -> str:
    return "".join(ch.lower() for ch in name if ch.isalnum())


def branch_by_tokens(tree, required_tokens: tuple[str, ...]) -> str | None:
    for name in tree.keys():
        norm = normalize_branch(str(name))
        if all(token in norm for token in required_tokens):
            return str(name)
    return None


def step_level_depth_data(root: uproot.ReadOnlyDirectory):
    if "StepLET" not in root:
        return None
    tree = root["StepLET"]
    absorber = branch_by_tokens(tree, ("absorber",))
    category = branch_by_tokens(tree, ("deposit", "category"))
    eion = branch_by_tokens(tree, ("eion", "kev"))
    weight = branch_by_tokens(tree, ("weight",))
    norm_depth = branch_by_tokens(tree, ("norm", "depth"))
    if not all([absorber, category, eion, weight, norm_depth]):
        return None
    a = tree.arrays([absorber, category, eion, weight, norm_depth], library="np")
    mask = np.asarray(a[absorber]) == 2
    return (
        np.asarray(a[norm_depth][mask], dtype=float) * 50.0,
        np.asarray(a[category][mask], dtype=int),
        np.asarray(a[eion][mask], dtype=float) * np.asarray(a[weight][mask], dtype=float),
        "step-level local deposition",
    )


def track_level_depth_data(root: uproot.ReadOnlyDirectory):
    tree = root["TrackLET"]
    a = tree.arrays(
        ["absorberID", "depositCategory", "eion_keV", "trackWeight", "meanDepth_um"],
        library="np",
    )
    mask = a["absorberID"] == 2
    return (
        np.asarray(a["meanDepth_um"][mask], dtype=float),
        np.asarray(a["depositCategory"][mask], dtype=int),
        np.asarray(a["eion_keV"][mask], dtype=float)
        * np.asarray(a["trackWeight"][mask], dtype=float),
        "Eion-weighted track-mean depth",
    )


def grouped_particle_depth(root: uproot.ReadOnlyDirectory):
    raw = step_level_depth_data(root)
    if raw is None:
        raw = track_level_depth_data(root)
    depth, category, energy, method = raw
    finite = np.isfinite(depth) & np.isfinite(energy) & (energy > 0) & (depth >= 0) & (depth <= 50)
    depth = depth[finite]
    category = category[finite]
    energy = energy[finite]
    edges = np.linspace(0.0, 50.0, 101)
    widths = np.diff(edges)
    total_energy = float(np.sum(energy))

    groups = {
        0: [0],
        1: [1],
        3: [3],
        7: [7],
        99: [2, 4, 5, 6, 8, 9, 10, 11],
    }
    names = {
        0: "e$^-$ (gamma-mediated)",
        1: "e$^-$ (ion-mediated)",
        3: "proton",
        7: "carbon",
        99: "other particles",
    }
    profiles = {}
    fractions = {}
    for group, cats in groups.items():
        select = np.isin(category, cats)
        hist, _ = np.histogram(depth[select], bins=edges, weights=energy[select])
        profiles[group] = hist / total_energy / widths if total_energy > 0 else np.zeros_like(widths)
        fractions[group] = float(np.sum(energy[select]) / total_energy) if total_energy > 0 else 0.0
    return edges, profiles, fractions, names, method


def plot_aggregate_let(pdf_pages: PdfPages) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.7))
    for label, root in ROOTS.items():
        values, edges = aggregate_let(root)
        ax.stairs(normalized(values), edges, label=label)
    ax.set_xlabel(r"Electronic stopping power, $LET_{calc}$ (keV $\mu$m$^{-1}$)")
    ax.set_ylabel("Fraction of ionizing energy per logarithmic bin")
    ax.set_xlim(1e-2, 1e4)
    finish_axes(ax, logx=True, logy=True)
    fig.tight_layout()
    save_figure(fig, "05_DD_50um_with_without_Pb_Eion_weighted_LET", pdf_pages)


def plot_particle_let(pdf_pages: PdfPages, partitions: dict[str, dict[int, float]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.5), sharex=True, sharey=True)
    categories = [0, 1, 3, 7, 9]
    for ax, (label, root) in zip(axes, ROOTS.items()):
        for cat in categories:
            title = f"LETcalc Eion-weighted spectrum, absorber 2, {CATEGORY_TOKENS[cat]}"
            obj = find_object(root, title)
            values, edges = hist_data(obj)
            if np.sum(values) <= 0:
                continue
            frac = partitions[label].get(cat, 0.0)
            ax.stairs(normalized(values), edges, label=f"{CATEGORY_NAMES[cat]} ({frac:.1%})")
        ax.set_title(label)
        ax.set_xlabel(r"$LET_{calc}$ (keV $\mu$m$^{-1}$)")
        ax.set_xlim(1e-2, 1e4)
        finish_axes(ax, logx=True, logy=True, legend=True)
    axes[0].set_ylabel("Fraction of category ionizing energy per logarithmic bin")
    fig.tight_layout()
    save_figure(fig, "07_DD_50um_particle_resolved_Eion_weighted_LET", pdf_pages)


def plot_depth_overlay(pdf_pages: PdfPages) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.7))
    for label, root in ROOTS.items():
        values, edges = hist_data(find_object(root, "dEion/dz absolute depth profile, absorber 2"))
        density = normalized(values) / np.diff(edges)
        ax.stairs(density, edges, label=label)
    ax.set_xlabel(r"Depth in active polymer ($\mu$m)")
    ax.set_ylabel(r"Normalized ionizing-energy density ($\mu$m$^{-1}$)")
    ax.set_xlim(0, 50)
    finish_axes(ax)
    fig.tight_layout()
    save_figure(fig, "09_DD_50um_with_without_Pb_depth_profile", pdf_pages)


def plot_mean_eion(pdf_pages: PdfPages, summary: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(5.8, 4.7))
    labels = ["No Pb", "With Pb"]
    vals = summary["mean_Eion_per_primary_keV"].to_numpy()
    bars = ax.bar(np.arange(2), vals, width=0.58)
    ax.set_xticks(np.arange(2), labels)
    ax.set_ylabel(r"Mean $E_{ion}$ per incident primary (keV)")
    ax.grid(True, axis="y", alpha=0.24, linewidth=0.55)
    ax.tick_params(direction="in", top=True, right=True)
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, val, f"{val:.3g}", ha="center", va="bottom", fontsize=8.5)
    fig.tight_layout()
    save_figure(fig, "11_DD_50um_mean_Eion_per_primary", pdf_pages)


def plot_high_let(pdf_pages: PdfPages, summary: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.7))
    x = np.arange(2)
    width = 0.24
    columns = [
        ("Eion_fraction_above_10_keV_per_um", r"$LET_{calc}\geq10$ keV $\mu$m$^{-1}$"),
        ("Eion_fraction_above_100_keV_per_um", r"$LET_{calc}\geq100$ keV $\mu$m$^{-1}$"),
        ("Eion_fraction_above_1000_keV_per_um", r"$LET_{calc}\geq1000$ keV $\mu$m$^{-1}$"),
    ]
    for offset, (col, legend) in zip([-width, 0, width], columns):
        ax.bar(x + offset, summary[col], width, label=legend)
    ax.set_xticks(x, ["No Pb", "With Pb"])
    ax.set_ylabel("Fraction of deposited ionizing energy")
    ax.set_ylim(0, 1.05)
    finish_axes(ax, legend=True)
    fig.tight_layout()
    save_figure(fig, "12_DD_50um_high_LET_energy_fractions", pdf_pages)


def plot_hit_eion(pdf_pages: PdfPages) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.7))
    for label, root in ROOTS.items():
        values, edges = hist_data(find_object(root, "Eion per event (hit-only), absorber 2"))
        ax.stairs(normalized(values), edges, label=label)
    ax.set_xlabel("Ionizing energy deposited per hit (keV)")
    ax.set_ylabel("Fraction of hit events per logarithmic bin")
    ax.set_xlim(1e-2, 2e4)
    finish_axes(ax, logx=True, logy=True)
    fig.tight_layout()
    save_figure(fig, "13_DD_50um_hit_event_Eion", pdf_pages)


def plot_particle_depth_figures(pdf_pages: PdfPages, depth_results) -> pd.DataFrame:
    ymax = 0.0
    for _, (_, profiles, _, _, _) in depth_results.items():
        for profile in profiles.values():
            ymax = max(ymax, float(np.nanmax(profile)))
    ymax = ymax * 1.08 if ymax > 0 else 1.0

    rows = []
    stems = {
        "DD, 50 $\\mu$m, no Pb": "14_DD_50um_noPb_particle_depth_profile",
        "DD, 50 $\\mu$m, with Pb": "15_DD_50um_withPb_particle_depth_profile",
    }
    for label, (edges, profiles, fractions, names, method) in depth_results.items():
        fig, ax = plt.subplots(figsize=(7.2, 4.7))
        for group in [0, 1, 3, 7, 99]:
            if fractions[group] <= 0:
                continue
            ax.stairs(
                profiles[group],
                edges,
                label=f"{names[group]} ({fractions[group]:.1%})",
                color=CATEGORY_COLORS[group],
            )
            rows.append(
                {
                    "source": label.replace("$\\mu$", "u"),
                    "particle_group": names[group].replace("$", ""),
                    "Eion_fraction": fractions[group],
                    "depth_method": method,
                }
            )
        ax.set_title(label)
        ax.set_xlabel(r"Depth in active polymer ($\mu$m)")
        ax.set_ylabel(r"Fraction of total $E_{ion}$ per unit depth ($\mu$m$^{-1}$)")
        ax.set_xlim(0, 50)
        ax.set_ylim(0, ymax)
        finish_axes(ax, legend=True)
        ax.text(
            0.99,
            0.02,
            method,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=7.5,
            alpha=0.72,
        )
        fig.tight_layout()
        save_figure(fig, stems[label], pdf_pages)
    return pd.DataFrame(rows)


def write_readme(summary: pd.DataFrame, depth_methods: dict[str, str]) -> None:
    no_pb = summary.iloc[0]
    with_pb = summary.iloc[1]

    def pct_change(column: str) -> float:
        base = float(no_pb[column])
        return 100.0 * (float(with_pb[column]) / base - 1.0) if base != 0 else math.nan

    text = f"""DD 50 um with/without Pb device-physics figure package
========================================================

Inputs
------
- build/DD-50um-noPb.root
- build/DD-50um-Pb.root

Generated comparisons
---------------------
05: Aggregate Eion-weighted LETcalc distribution.
07: Particle-resolved Eion-weighted LETcalc, matched two-panel layout.
09: Total ionizing-energy deposition depth profile.
11: Mean Eion deposited per incident primary.
12: Fractions of Eion deposited above 10, 100, and 1000 keV/um.
13: Eion deposited per hit event.
14: Particle-resolved depth profile, no Pb.
15: Particle-resolved depth profile, with Pb.

Particle depth normalization
----------------------------
Each particle-depth curve is divided by total Eion from all particle classes and by bin width.
Therefore, the integral of a curve equals that particle group's fraction of total ionizing energy.
The two separate depth figures use identical axes, binning, particle colors, and y-axis limits.

Depth scoring method
--------------------
- No Pb: {depth_methods[list(SOURCES.keys())[0]]}
- With Pb: {depth_methods[list(SOURCES.keys())[1]]}

If StepLET was present with the required branches, local step-level Eion deposition was used.
Otherwise, the TrackLET meanDepth_um value was weighted by track Eion; this is a track-mean
approximation and is labeled as such on the figures.

Selected numerical changes with Pb
----------------------------------
- Hit fraction change: {pct_change('hit_fraction'):+.2f}%
- Mean Eion per primary change: {pct_change('mean_Eion_per_primary_keV'):+.2f}%
- Median Eion per hit change: {pct_change('median_Eion_per_hit_keV'):+.2f}%
- Energy-weighted median LETcalc change: {pct_change('energy_weighted_median_LETcalc_keV_per_um'):+.2f}%
- Fraction above 10 keV/um change: {100*(with_pb['Eion_fraction_above_10_keV_per_um']-no_pb['Eion_fraction_above_10_keV_per_um']):+.2f} percentage points
- Fraction above 100 keV/um change: {100*(with_pb['Eion_fraction_above_100_keV_per_um']-no_pb['Eion_fraction_above_100_keV_per_um']):+.2f} percentage points
- Fraction above 1000 keV/um change: {100*(with_pb['Eion_fraction_above_1000_keV_per_um']-no_pb['Eion_fraction_above_1000_keV_per_um']):+.2f} percentage points

Formats
-------
Each figure is supplied as vector PDF, editable SVG, and 400 dpi PNG.
The combined PDF places the figures in numerical order.
"""
    (OUT / "README.txt").write_text(text, encoding="utf-8")


def main() -> None:
    summaries = [event_summary(label, root) for label, root in ROOTS.items()]
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT / "DD_50um_with_without_Pb_event_summary.csv", index=False)

    partitions = {label: particle_partition(root) for label, root in ROOTS.items()}
    partition_rows = []
    for label, fractions in partitions.items():
        for cat, fraction in fractions.items():
            partition_rows.append(
                {
                    "source": label.replace("$\\mu$", "u"),
                    "category_id": cat,
                    "category": CATEGORY_NAMES.get(cat, str(cat)),
                    "Eion_fraction": fraction,
                }
            )
    pd.DataFrame(partition_rows).to_csv(
        OUT / "DD_50um_with_without_Pb_particle_partition.csv", index=False
    )

    depth_results = {label: grouped_particle_depth(root) for label, root in ROOTS.items()}
    depth_methods = {label: result[4] for label, result in depth_results.items()}

    combined = OUT / "DD_50um_with_without_Pb_figures_combined.pdf"
    with PdfPages(combined) as pdf_pages:
        plot_aggregate_let(pdf_pages)
        plot_particle_let(pdf_pages, partitions)
        plot_depth_overlay(pdf_pages)
        plot_mean_eion(pdf_pages, summary)
        plot_high_let(pdf_pages, summary)
        plot_hit_eion(pdf_pages)
        depth_table = plot_particle_depth_figures(pdf_pages, depth_results)

    depth_table.to_csv(OUT / "DD_50um_particle_depth_summary.csv", index=False)
    write_readme(summary, depth_methods)

    package = Path("DD_50um_with_without_Pb_figure_package")
    if package.exists():
        shutil.rmtree(package)
    shutil.copytree(OUT, package)
    shutil.make_archive(str(package), "zip", root_dir=package)


if __name__ == "__main__":
    main()
