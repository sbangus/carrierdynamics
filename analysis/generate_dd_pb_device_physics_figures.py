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

OUT = Path("dd_pb_device_physics_figures")
OUT.mkdir(exist_ok=True)

SOURCES = {
    "DD without Pb": Path("build/DD-without-Pb.root"),
    "DD with Pb": Path("build/DD-with-Pb.root"),
}
ROOTS = {k: uproot.open(v) for k, v in SOURCES.items()}

CATEGORY_TITLES = [
    ("e-(gamma)", "e$^-$ (gamma-mediated)"),
    ("e-(ion)", "e$^-$ (ion-mediated)"),
    ("proton", "proton"),
    ("lightIon", "light ion"),
    ("alpha", "alpha"),
    ("Li", "lithium"),
    ("C", "carbon"),
    ("Be/B", "Be/B"),
    ("otherHeavyIon", "other heavy ion"),
]

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.5,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "legend.fontsize": 8.2,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.55,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.bbox": "tight",
})


def find_object(root, title: str):
    for _, obj in root.items(recursive=True):
        if getattr(obj, "title", "") == title:
            return obj
    raise KeyError(title)


def find_prefix(root, prefix: str):
    return [obj for _, obj in root.items(recursive=True)
            if getattr(obj, "title", "").startswith(prefix)]


def hist(obj):
    return np.asarray(obj.values(flow=False), float), np.asarray(obj.axis().edges(), float)


def norm(v):
    s = np.nansum(v)
    return v / s if s > 0 else np.zeros_like(v)


def style(ax, logx=False, logy=False, legend=True):
    if logx: ax.set_xscale("log")
    if logy: ax.set_yscale("log")
    ax.grid(True, which="major", alpha=.24, linewidth=.55)
    ax.grid(True, which="minor", alpha=.10, linewidth=.4)
    ax.tick_params(direction="in", top=True, right=True)
    if legend: ax.legend(frameon=False, handlelength=2.5)


def save(fig, stem, pdf):
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.svg")
    fig.savefig(OUT / f"{stem}.png", dpi=400)
    pdf.savefig(fig)
    plt.close(fig)


def aggregate_let(root):
    objs = find_prefix(root, "LETcalc Eion-weighted spectrum, absorber 2,")
    total = None; edges = None
    for obj in objs:
        v, e = hist(obj)
        if total is None:
            total = np.zeros_like(v); edges = e
        total += v
    if total is None: raise RuntimeError("No LET spectra")
    return total, edges


def weighted_median(x, w):
    m = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not np.any(m): return math.nan
    x, w = x[m], w[m]
    o = np.argsort(x); x, w = x[o], w[o]
    return float(x[np.searchsorted(np.cumsum(w), .5*np.sum(w))])


def summary(label, root):
    t = root["EventLET"]
    cols = ["absorberID", "trackWeight", "eion_keV", "niel_keV",
            "letDcalc_keV_per_um", "maxLetCalc_keV_per_um",
            "fractionEionAbove10", "fractionEionAbove100", "fractionEionAbove1000",
            "meanDepth_um", "depthSigma_um"]
    a = t.arrays(cols, library="np")
    m = a["absorberID"] == 2
    w = a["trackWeight"][m].astype(float)
    e = a["eion_keV"][m].astype(float)
    ew = w*e; hit=e>0; tot=ew.sum()
    def eavg(name):
        x=a[name][m].astype(float)
        return float(np.sum(ew*x)/tot) if tot>0 else math.nan
    return {
        "condition": label,
        "incident_events": float(w.sum()),
        "hit_fraction": float(w[hit].sum()/w.sum()),
        "mean_Eion_per_primary_keV": float(tot/w.sum()),
        "median_Eion_per_hit_keV": weighted_median(e[hit], w[hit]),
        "energy_weighted_median_LETcalc_keV_per_um": weighted_median(a["letDcalc_keV_per_um"][m][hit], ew[hit]),
        "mean_event_LETcalc_Eion_weighted_keV_per_um": eavg("letDcalc_keV_per_um"),
        "mean_max_LETcalc_Eion_weighted_keV_per_um": eavg("maxLetCalc_keV_per_um"),
        "Eion_fraction_above_10_keV_per_um": eavg("fractionEionAbove10"),
        "Eion_fraction_above_100_keV_per_um": eavg("fractionEionAbove100"),
        "Eion_fraction_above_1000_keV_per_um": eavg("fractionEionAbove1000"),
        "Eion_weighted_mean_depth_um": eavg("meanDepth_um"),
        "Eion_weighted_mean_depth_sigma_um": eavg("depthSigma_um"),
    }


def particle_partition(label, root):
    t=root["TrackLET"]
    a=t.arrays(["absorberID","depositCategory","eion_keV","trackWeight"], library="np")
    m=a["absorberID"]==2
    c=a["depositCategory"][m]; e=a["eion_keV"][m]*a["trackWeight"][m]
    total=e.sum(); rows=[]
    names={0:"e- (gamma-mediated)",1:"e- (ion-mediated)",2:"e+",3:"proton",4:"light ion",5:"alpha",6:"lithium",7:"carbon",8:"Be/B",9:"other heavy ion",10:"neutral local",11:"other"}
    for k in sorted(np.unique(c)):
        amt=e[c==k].sum()
        rows.append({"condition":label,"category_id":int(k),"category":names.get(int(k),str(k)),"Eion_keV":float(amt),"Eion_fraction":float(amt/total)})
    return rows


summ = pd.DataFrame([summary(k,r) for k,r in ROOTS.items()])
summ.to_csv(OUT/"dd_pb_event_level_summary.csv", index=False)
pd.DataFrame([row for k,r in ROOTS.items() for row in particle_partition(k,r)]).to_csv(OUT/"dd_pb_particle_Eion_partition.csv", index=False)

with PdfPages(OUT/"dd_pb_figures_combined.pdf") as pdf:
    # Figure 5 replacement: aggregate Eion-weighted LET
    fig,ax=plt.subplots(figsize=(7.2,4.7))
    for label,root in ROOTS.items():
        v,e=aggregate_let(root); ax.stairs(norm(v),e,label=label)
    ax.set_xlabel(r"Electronic stopping power, $LET_{calc}$ (keV $\mu$m$^{-1}$)")
    ax.set_ylabel("Fraction of ionizing energy per logarithmic bin")
    ax.set_xlim(1e-2,1e4); style(ax,True,True); fig.tight_layout(); save(fig,"05_DD_Pb_Eion_weighted_LET",pdf)

    # Figure 7 replacement: particle-resolved LET, matched panels
    fig,axs=plt.subplots(1,2,figsize=(10.6,4.5),sharex=True,sharey=True)
    for ax,(label,root) in zip(axs,ROOTS.items()):
        for cat,display in CATEGORY_TITLES:
            try: v,e=hist(find_object(root,f"LETcalc Eion-weighted spectrum, absorber 2, {cat}"))
            except KeyError: continue
            if v.sum()>0: ax.stairs(norm(v),e,label=display)
        ax.set_title(label); ax.set_xlim(1e-2,1e4); style(ax,True,True,False)
        ax.set_xlabel(r"$LET_{calc}$ (keV $\mu$m$^{-1}$)")
    axs[0].set_ylabel("Fraction of category ionizing energy per logarithmic bin")
    handles,labels=axs[1].get_legend_handles_labels(); fig.legend(handles,labels,loc="center right",frameon=False,bbox_to_anchor=(1.0,.5))
    fig.tight_layout(rect=(0,0,.84,1)); save(fig,"07_DD_Pb_particle_resolved_Eion_weighted_LET",pdf)

    # Figure 9 replacement: depth profile
    fig,ax=plt.subplots(figsize=(7.2,4.7))
    for label,root in ROOTS.items():
        v,e=hist(find_object(root,"dEion/dz absolute depth profile, absorber 2")); p=norm(v)/np.diff(e)
        ax.stairs(p,e,label=label)
    ax.set_xlabel(r"Depth in active polymer ($\mu$m)"); ax.set_ylabel(r"Normalized ionizing-energy density ($\mu$m$^{-1}$)")
    ax.set_xlim(0,60); style(ax); fig.tight_layout(); save(fig,"09_DD_Pb_depth_profile",pdf)

    # Figure 11 replacement: mean Eion per primary
    fig,ax=plt.subplots(figsize=(6.2,4.6)); x=np.arange(2)
    ax.bar(x,summ["mean_Eion_per_primary_keV"]); ax.set_xticks(x,summ["condition"])
    ax.set_ylabel(r"Mean $E_{ion}$ per incident primary (keV)"); ax.set_yscale("log"); style(ax,False,False,False)
    fig.tight_layout(); save(fig,"11_DD_Pb_mean_Eion_per_primary",pdf)

    # Figure 12 replacement: high LET fractions
    fig,ax=plt.subplots(figsize=(6.8,4.7)); width=.24
    ax.bar(x-width,summ["Eion_fraction_above_10_keV_per_um"],width,label=r"$LET_{calc}\geq10$ keV $\mu$m$^{-1}$")
    ax.bar(x,summ["Eion_fraction_above_100_keV_per_um"],width,label=r"$LET_{calc}\geq100$ keV $\mu$m$^{-1}$")
    ax.bar(x+width,summ["Eion_fraction_above_1000_keV_per_um"],width,label=r"$LET_{calc}\geq1000$ keV $\mu$m$^{-1}$")
    ax.set_xticks(x,summ["condition"]); ax.set_ylim(0,1.05); ax.set_ylabel("Fraction of deposited ionizing energy"); style(ax)
    fig.tight_layout(); save(fig,"12_DD_Pb_high_LET_energy_fractions",pdf)

    # Figure 13 replacement: Eion per hit
    fig,ax=plt.subplots(figsize=(7.2,4.7))
    for label,root in ROOTS.items():
        v,e=hist(find_object(root,"Eion per event (hit-only), absorber 2")); ax.stairs(norm(v),e,label=label)
    ax.set_xlabel("Ionizing energy deposited per hit (keV)"); ax.set_ylabel("Fraction of hit events per logarithmic bin")
    ax.set_xlim(1e-2,2e4); style(ax,True,True); fig.tight_layout(); save(fig,"13_DD_Pb_hit_event_Eion",pdf)

# Add a compact text summary to package.
with open(OUT/"README.txt","w",encoding="utf-8") as f:
    f.write("DD with/without Pb device-physics comparison figures.\n")
    f.write("Figures correspond to replacements for original Figures 5, 7, 9, 11, 12, and 13.\n")
    f.write("Each is supplied as PDF, SVG, and 400-dpi PNG.\n\n")
    f.write(summ.to_string(index=False))

shutil.make_archive("dd_pb_device_physics_figure_package","zip",root_dir=OUT)
print(summ.to_string(index=False))
