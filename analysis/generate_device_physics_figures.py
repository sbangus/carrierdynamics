#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import shutil
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from matplotlib.backends.backend_pdf import PdfPages

OUT = Path("device_physics_figures")
OUT.mkdir(exist_ok=True)

SOURCES = {
    "Gamma, 662 keV": Path("build/bare-detector-662keV-gamma.root"),
    "Gamma, 2.2 MeV": Path("build/bare-detector-2_2MeV-gamma.root"),
    "DD neutrons, 2.45 MeV": Path("build/bare-detector-DD-neutron.root"),
    "DT neutrons, 14.1 MeV": Path("build/bare-detector-DT-neutron.root"),
    "Thermal + fast, no B$_4$C": Path("build/bare-detector-thermal-fast-noB4C.root"),
    "Thermal + fast, with B$_4$C": Path("build/bare-detector-thermal-fast-withB4C.root"),
    "Thermal only, with B$_4$C": Path("build/bare-detector-thermal-only.root"),
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

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "legend.fontsize": 8.3,
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


def find_objects_starting(root: uproot.ReadOnlyDirectory, title_prefix: str):
    found = []
    for _, obj in root.items(recursive=True):
        title = getattr(obj, "title", "")
        if title.startswith(title_prefix):
            found.append(obj)
    return found


def hist_data(obj):
    values = np.asarray(obj.values(flow=False), dtype=float)
    edges = np.asarray(obj.axis().edges(), dtype=float)
    return values, edges


def bin_centers(edges: np.ndarray) -> np.ndarray:
    positive = np.all(edges > 0)
    if positive:
        return np.sqrt(edges[:-1] * edges[1:])
    return 0.5 * (edges[:-1] + edges[1:])


def normalized(values: np.ndarray) -> np.ndarray:
    total = float(np.nansum(values))
    return values / total if total > 0 else np.zeros_like(values)


def save_figure(fig: plt.Figure, stem: str, pdf_pages: PdfPages) -> None:
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.svg")
    fig.savefig(OUT / f"{stem}.png", dpi=400)
    pdf_pages.savefig(fig)
    plt.close(fig)


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


def plot_hist_comparison(
    pdf_pages: PdfPages,
    stem: str,
    title: str,
    source_labels: Iterable[str],
    xlabel: str,
    ylabel: str,
    *,
    normalize: bool = True,
    logx: bool = True,
    logy: bool = True,
    xlim: tuple[float, float] | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.7))
    plotted = 0
    for source_label in source_labels:
        obj = find_object(ROOTS[source_label], title)
        values, edges = hist_data(obj)
        if np.nansum(values) <= 0:
            continue
        y = normalized(values) if normalize else values
        ax.stairs(y, edges, label=source_label)
        plotted += 1
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xlim is not None:
        ax.set_xlim(*xlim)
    finish_axes(ax, logx=logx, logy=logy, legend=plotted > 1)
    fig.tight_layout()
    save_figure(fig, stem, pdf_pages)


def aggregate_let(root: uproot.ReadOnlyDirectory) -> tuple[np.ndarray, np.ndarray]:
    objects = find_objects_starting(root, "LETcalc Eion-weighted spectrum, absorber 2,")
    total = None
    edges = None
    for obj in objects:
        values, obj_edges = hist_data(obj)
        if total is None:
            total = np.zeros_like(values)
            edges = obj_edges
        total += values
    if total is None or edges is None:
        raise RuntimeError("No LETcalc Eion-weighted spectra found")
    return total, edges


def plot_aggregate_let(pdf_pages: PdfPages, stem: str, labels: list[str]) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.7))
    for label in labels:
        values, edges = aggregate_let(ROOTS[label])
        ax.stairs(normalized(values), edges, label=label)
    ax.set_xlabel(r"Electronic stopping power, $LET_{calc}$ (keV $\mu$m$^{-1}$)")
    ax.set_ylabel("Fraction of ionizing energy per logarithmic bin")
    ax.set_xlim(1e-2, 1e4)
    finish_axes(ax, logx=True, logy=True)
    fig.tight_layout()
    save_figure(fig, stem, pdf_pages)


def plot_particle_resolved_let(pdf_pages: PdfPages, stem: str, source: str, categories: list[str]) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.7))
    for category in categories:
        title = f"LETcalc Eion-weighted spectrum, absorber 2, {category}"
        obj = find_object(ROOTS[source], title)
        values, edges = hist_data(obj)
        if np.nansum(values) <= 0:
            continue
        ax.stairs(normalized(values), edges, label=CATEGORY_NAMES.get({
            "e-(gamma)": 0, "e-(ion)": 1, "e+": 2, "proton": 3,
            "lightIon": 4, "alpha": 5, "Li": 6, "C": 7,
            "Be/B": 8, "otherHeavyIon": 9, "neutralLocal": 10, "other": 11,
        }[category], category))
    ax.set_xlabel(r"Electronic stopping power, $LET_{calc}$ (keV $\mu$m$^{-1}$)")
    ax.set_ylabel("Fraction of category ionizing energy per logarithmic bin")
    ax.set_xlim(1e-2, 1e4)
    finish_axes(ax, logx=True, logy=True)
    fig.tight_layout()
    save_figure(fig, stem, pdf_pages)


def plot_depth(pdf_pages: PdfPages, stem: str, labels: list[str]) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.7))
    title = "dEion/dz absolute depth profile, absorber 2"
    for label in labels:
        values, edges = hist_data(find_object(ROOTS[label], title))
        if np.nansum(values) <= 0:
            continue
        widths = np.diff(edges)
        profile = normalized(values) / widths
        ax.stairs(profile, edges, label=label)
    ax.set_xlabel(r"Depth in active polymer ($\mu$m)")
    ax.set_ylabel(r"Normalized ionizing-energy density ($\mu$m$^{-1}$)")
    ax.set_xlim(0, 60)
    finish_axes(ax, logx=False, logy=False)
    fig.tight_layout()
    save_figure(fig, stem, pdf_pages)


def plot_event_eion(pdf_pages: PdfPages, stem: str, labels: list[str]) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.7))
    title = "Eion per event (hit-only), absorber 2"
    for label in labels:
        values, edges = hist_data(find_object(ROOTS[label], title))
        if np.nansum(values) <= 0:
            continue
        ax.stairs(normalized(values), edges, label=label)
    ax.set_xlabel(r"Ionizing energy deposited per hit (keV)")
    ax.set_ylabel("Fraction of hit events per logarithmic bin")
    ax.set_xlim(1e-2, 2e4)
    finish_axes(ax, logx=True, logy=True)
    fig.tight_layout()
    save_figure(fig, stem, pdf_pages)


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(mask):
        return math.nan
    v = values[mask]
    w = weights[mask]
    order = np.argsort(v)
    v = v[order]
    w = w[order]
    return float(v[np.searchsorted(np.cumsum(w), 0.5 * np.sum(w))])


def event_summary(label: str, root: uproot.ReadOnlyDirectory) -> dict[str, float | str]:
    tree = root["EventLET"]
    cols = [
        "absorberID", "trackWeight", "eion_keV", "niel_keV",
        "letDcalc_keV_per_um", "maxLetCalc_keV_per_um",
        "fractionEionAbove10", "fractionEionAbove100", "fractionEionAbove1000",
        "meanDepth_um", "depthSigma_um", "nChargedTracks",
    ]
    a = tree.arrays(cols, library="np")
    mask = a["absorberID"] == 2
    w = a["trackWeight"][mask].astype(float)
    eion = a["eion_keV"][mask].astype(float)
    niel = a["niel_keV"][mask].astype(float)
    hit = eion > 0
    eweight = w * eion
    total_eion = np.sum(eweight)
    total_dep_proxy = np.sum(w * (eion + niel))

    def eweighted(name: str) -> float:
        x = a[name][mask].astype(float)
        return float(np.sum(eweight * x) / total_eion) if total_eion > 0 else math.nan

    return {
        "source": label.replace("$_4$", "4"),
        "incident_events": int(np.sum(w)),
        "hit_events": int(np.sum(w[hit])),
        "hit_fraction": float(np.sum(w[hit]) / np.sum(w)) if np.sum(w) else math.nan,
        "mean_Eion_per_primary_keV": float(total_eion / np.sum(w)) if np.sum(w) else math.nan,
        "median_Eion_per_hit_keV": weighted_median(eion[hit], w[hit]),
        "energy_weighted_median_LETcalc_keV_per_um": weighted_median(
            a["letDcalc_keV_per_um"][mask][hit].astype(float), eweight[hit]
        ),
        "mean_event_LETcalc_Eion_weighted_keV_per_um": eweighted("letDcalc_keV_per_um"),
        "mean_max_LETcalc_Eion_weighted_keV_per_um": eweighted("maxLetCalc_keV_per_um"),
        "Eion_fraction_above_10_keV_per_um": eweighted("fractionEionAbove10"),
        "Eion_fraction_above_100_keV_per_um": eweighted("fractionEionAbove100"),
        "Eion_fraction_above_1000_keV_per_um": eweighted("fractionEionAbove1000"),
        "Eion_weighted_mean_depth_um": eweighted("meanDepth_um"),
        "Eion_weighted_mean_depth_sigma_um": eweighted("depthSigma_um"),
        "NIEL_fraction_of_Edep": float(np.sum(w * niel) / total_dep_proxy) if total_dep_proxy > 0 else math.nan,
    }


def particle_energy_partition(label: str, root: uproot.ReadOnlyDirectory) -> list[dict[str, float | str]]:
    tree = root["TrackLET"]
    a = tree.arrays(["absorberID", "depositCategory", "eion_keV", "trackWeight"], library="np")
    mask = a["absorberID"] == 2
    cat = a["depositCategory"][mask]
    e = a["eion_keV"][mask] * a["trackWeight"][mask]
    total = float(np.sum(e))
    rows = []
    for category in sorted(np.unique(cat)):
        amount = float(np.sum(e[cat == category]))
        rows.append(
            {
                "source": label.replace("$_4$", "4"),
                "category_id": int(category),
                "category": CATEGORY_NAMES.get(int(category), str(category)),
                "Eion_keV": amount,
                "Eion_fraction": amount / total if total > 0 else math.nan,
            }
        )
    return rows


def plot_summary_bars(pdf_pages: PdfPages, summary: pd.DataFrame) -> None:
    source_order = [label.replace("$_4$", "4") for label in SOURCES]
    s = summary.set_index("source").loc[source_order]

    fig, ax = plt.subplots(figsize=(7.2, 4.7))
    x = np.arange(len(s))
    ax.bar(x, s["mean_Eion_per_primary_keV"])
    ax.set_xticks(x, [
        "662 keV $\gamma$", "2.2 MeV $\gamma$", "DD", "DT",
        "thermal+fast\nno B$_4$C", "thermal+fast\nwith B$_4$C", "thermal only\nwith B$_4$C",
    ], rotation=25, ha="right")
    ax.set_ylabel(r"Mean $E_{ion}$ per incident primary (keV)")
    ax.set_yscale("log")
    ax.grid(True, axis="y", which="both", alpha=0.24, linewidth=0.55)
    ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout()
    save_figure(fig, "11_mean_Eion_per_primary", pdf_pages)

    fig, ax = plt.subplots(figsize=(7.2, 4.7))
    width = 0.25
    ax.bar(x - width, s["Eion_fraction_above_10_keV_per_um"], width, label=r"$LET_{calc}\geq10$ keV $\mu$m$^{-1}$")
    ax.bar(x, s["Eion_fraction_above_100_keV_per_um"], width, label=r"$LET_{calc}\geq100$ keV $\mu$m$^{-1}$")
    ax.bar(x + width, s["Eion_fraction_above_1000_keV_per_um"], width, label=r"$LET_{calc}\geq1000$ keV $\mu$m$^{-1}$")
    ax.set_xticks(x, [
        "662 keV $\gamma$", "2.2 MeV $\gamma$", "DD", "DT",
        "thermal+fast\nno B$_4$C", "thermal+fast\nwith B$_4$C", "thermal only\nwith B$_4$C",
    ], rotation=25, ha="right")
    ax.set_ylabel("Fraction of deposited ionizing energy")
    ax.set_ylim(0, 1.03)
    ax.grid(True, axis="y", alpha=0.24, linewidth=0.55)
    ax.tick_params(direction="in", top=True, right=True)
    ax.legend(frameon=False)
    fig.tight_layout()
    save_figure(fig, "12_high_LET_energy_fractions", pdf_pages)


def write_captions() -> None:
    text = """# Device-physics simulation figures

All spectra are taken from absorber 2, the charge-active polymer layer. Shape comparisons are normalized to unit integral so that different Monte Carlo histories can be compared without implying equal physical source fluence. The event-level summary CSV retains absolute per-primary normalization.

## 01_gamma_first_generation_electron_KE
Kinetic-energy distributions at creation for first-generation secondary electrons generated by 662 keV and 2.2 MeV photons. The figure isolates the initial electron spectrum that carries gamma energy into the polymer and establishes the low-LET electron reference for comparison with neutron-generated recoils.

## 02_DD_DT_first_generation_proton_KE
Kinetic-energy distributions at creation for first-generation recoil protons under 2.45 MeV DD and 14.1 MeV DT neutron irradiation. The higher-energy DT field extends the recoil-proton spectrum to larger energies, changing both track range and the spatial distribution of ionization.

## 03_DD_DT_first_generation_carbon_KE
Kinetic-energy distributions at creation for first-generation carbon recoils under DD and DT irradiation. Carbon recoils are shorter-ranged and higher-LET than recoil protons; this comparison identifies the heavy-recoil component that can produce locally dense ionization.

## 04_thermal_capture_product_KE
Kinetic-energy spectra of first-generation alpha and lithium reaction products for the B4C-containing thermal-neutron configurations. These short-range ions are the principal high-LET products of boron neutron capture and deposit energy very differently from gamma-generated electrons and fast-neutron recoil protons.

## 05_gamma_fast_neutron_Eion_weighted_LET
Total ionizing-energy-weighted electronic stopping-power distributions for the two gamma fields and the DD and DT neutron fields. Each curve answers what fraction of the ionizing energy in the polymer was deposited at a given LET, rather than how many Geant4 transport steps occurred there.

## 06_thermal_configuration_Eion_weighted_LET
Total ionizing-energy-weighted LET distributions for mixed thermal/fast irradiation with and without B4C and for thermal-only irradiation with B4C. The comparison isolates how the boron converter and neutron-energy composition shift deposited energy toward high-LET alpha and lithium tracks.

## 07_DD_particle_resolved_Eion_weighted_LET
Particle-resolved ionizing-energy-weighted LET distributions for DD neutrons. Recoil protons, carbon recoils, and ion-mediated electrons occupy distinct LET regimes and together determine the initial charge-density distribution.

## 08_DT_particle_resolved_Eion_weighted_LET
Particle-resolved ionizing-energy-weighted LET distributions for DT neutrons. Relative to DD irradiation, the higher neutron energy changes the recoil-energy spectrum and the balance between long proton tracks, heavy recoils, and secondary electrons.

## 09_gamma_fast_neutron_depth_profile
Normalized depth distribution of ionizing-energy deposition through the 60 micrometer polymer layer for gamma, DD, and DT irradiation. This figure separates LET effects from deposition-depth effects, which can independently influence trapping, contact sensitivity, and charge-collection distance.

## 10_thermal_configuration_depth_profile
Normalized depth distribution of ionizing-energy deposition for the thermal-neutron configurations. The profile shows whether capture-product energy is deposited locally near converter-containing regions or distributed through the active layer.

## 11_mean_Eion_per_primary
Mean ionizing energy deposited in the active polymer per incident simulated primary. This plot retains absolute Monte Carlo normalization and therefore reports interaction probability and deposited-energy yield together. It should not be interpreted as collected charge without a transport and recombination model.

## 12_high_LET_energy_fractions
Fractions of total deposited ionizing energy occurring above 10, 100, and 1000 keV per micrometer. These compact event-level metrics quantify the high-LET content of each radiation field and provide direct independent variables for comparison with experimental overshoot, bias response, and relaxation behavior.

## Supplementary event-energy figures
The `supplementary_*_hit_event_Eion` figures show ionizing energy deposited per hit event. They distinguish the energy scale of individual interactions from the probability that an incident primary interacts in the polymer.
"""
    (OUT / "figure_descriptions.md").write_text(text, encoding="utf-8")


def main() -> None:
    pdf_path = OUT / "device_physics_figures_combined.pdf"
    with PdfPages(pdf_path) as pdf_pages:
        plot_hist_comparison(
            pdf_pages,
            "01_gamma_first_generation_electron_KE",
            "KE at generation of first-generation secondary e-",
            ["Gamma, 662 keV", "Gamma, 2.2 MeV"],
            "Electron kinetic energy at creation (MeV)",
            "Fraction of generated electrons per logarithmic bin",
            xlim=(1e-4, 3),
        )
        plot_hist_comparison(
            pdf_pages,
            "02_DD_DT_first_generation_proton_KE",
            "KE at generation of first-generation secondary proton",
            ["DD neutrons, 2.45 MeV", "DT neutrons, 14.1 MeV"],
            "Recoil-proton kinetic energy at creation (MeV)",
            "Fraction of recoil protons per logarithmic bin",
            xlim=(1e-3, 20),
        )
        plot_hist_comparison(
            pdf_pages,
            "03_DD_DT_first_generation_carbon_KE",
            "KE at generation of first-generation secondary carbon",
            ["DD neutrons, 2.45 MeV", "DT neutrons, 14.1 MeV"],
            "Carbon-recoil kinetic energy at creation (MeV)",
            "Fraction of carbon recoils per logarithmic bin",
            xlim=(1e-4, 20),
        )

        fig, ax = plt.subplots(figsize=(7.2, 4.7))
        thermal_sources = ["Thermal + fast, with B$_4$C", "Thermal only, with B$_4$C"]
        for source in thermal_sources:
            for particle, readable, linestyle in [
                ("alpha", r"$\alpha$", "-"),
                ("Li7", r"$^7$Li", "--"),
            ]:
                obj = find_object(ROOTS[source], f"KE at generation of first-generation secondary {particle}")
                values, edges = hist_data(obj)
                if np.nansum(values) <= 0:
                    continue
                ax.stairs(normalized(values), edges, label=f"{source}: {readable}", linestyle=linestyle)
        ax.set_xlabel("Capture-product kinetic energy at creation (MeV)")
        ax.set_ylabel("Fraction of reaction products per logarithmic bin")
        ax.set_xlim(0.1, 3)
        finish_axes(ax, logx=True, logy=True)
        fig.tight_layout()
        save_figure(fig, "04_thermal_capture_product_KE", pdf_pages)

        plot_aggregate_let(
            pdf_pages,
            "05_gamma_fast_neutron_Eion_weighted_LET",
            ["Gamma, 662 keV", "Gamma, 2.2 MeV", "DD neutrons, 2.45 MeV", "DT neutrons, 14.1 MeV"],
        )
        plot_aggregate_let(
            pdf_pages,
            "06_thermal_configuration_Eion_weighted_LET",
            ["Thermal + fast, no B$_4$C", "Thermal + fast, with B$_4$C", "Thermal only, with B$_4$C"],
        )
        plot_particle_resolved_let(
            pdf_pages,
            "07_DD_particle_resolved_Eion_weighted_LET",
            "DD neutrons, 2.45 MeV",
            ["e-(gamma)", "e-(ion)", "proton", "alpha", "C", "Be/B", "otherHeavyIon"],
        )
        plot_particle_resolved_let(
            pdf_pages,
            "08_DT_particle_resolved_Eion_weighted_LET",
            "DT neutrons, 14.1 MeV",
            ["e-(gamma)", "e-(ion)", "proton", "lightIon", "alpha", "C", "Be/B", "otherHeavyIon"],
        )
        plot_depth(
            pdf_pages,
            "09_gamma_fast_neutron_depth_profile",
            ["Gamma, 662 keV", "Gamma, 2.2 MeV", "DD neutrons, 2.45 MeV", "DT neutrons, 14.1 MeV"],
        )
        plot_depth(
            pdf_pages,
            "10_thermal_configuration_depth_profile",
            ["Thermal + fast, no B$_4$C", "Thermal + fast, with B$_4$C", "Thermal only, with B$_4$C"],
        )

        summary_rows = [event_summary(label, root) for label, root in ROOTS.items()]
        summary = pd.DataFrame(summary_rows)
        summary.to_csv(OUT / "event_level_summary.csv", index=False)
        plot_summary_bars(pdf_pages, summary)

        plot_event_eion(
            pdf_pages,
            "supplementary_gamma_fast_hit_event_Eion",
            ["Gamma, 662 keV", "Gamma, 2.2 MeV", "DD neutrons, 2.45 MeV", "DT neutrons, 14.1 MeV"],
        )
        plot_event_eion(
            pdf_pages,
            "supplementary_thermal_hit_event_Eion",
            ["Thermal + fast, no B$_4$C", "Thermal + fast, with B$_4$C", "Thermal only, with B$_4$C"],
        )

    partition_rows = []
    for label, root in ROOTS.items():
        partition_rows.extend(particle_energy_partition(label, root))
    pd.DataFrame(partition_rows).to_csv(OUT / "particle_Eion_partition.csv", index=False)
    write_captions()

    # Include analysis provenance and a compact manifest.
    manifest = sorted(p.name for p in OUT.iterdir())
    (OUT / "manifest.txt").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    shutil.make_archive("device_physics_thesis_figure_package", "zip", OUT)


if __name__ == "__main__":
    main()
