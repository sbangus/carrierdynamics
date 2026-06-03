#!/usr/bin/env python3
# =============================================================================
# Extract a 1D ROOT histogram as (energy, weight) pairs for Geant4 GPS-style
# tabulated spectra (e.g. /gps/ene/type Arb with a bias-format ASCII file).
#
# Dependencies: python3 -m pip install uproot numpy   (or: pip3 install uproot numpy)
#
# Typical usage from the build directory (where hadr05.root is produced):
#   python3 ../export_hist_for_gps.py -r hadr05.root -k 321 -o neutron_frontface_gps.txt
#   python3 ../export_hist_for_gps.py -r hadr05.root -k 321 --gps-dat -o neutron_frontface_gps.dat
#
# Histogram 321 axis is in eV (Run.cc). /gps/hist/file expects MeV; --gps-dat applies
# energy-scale 1e-6 (eV -> MeV) unless you override --energy-scale.
# =============================================================================

from __future__ import annotations

import argparse
import sys


def main() -> int:
    try:
        import numpy as np
        import uproot
    except ImportError:
        print(
            "Requires uproot and numpy:  python3 -m pip install uproot numpy",
            file=sys.stderr,
        )
        return 1

    parser = argparse.ArgumentParser(
        description="Export ROOT TH1 bin centers and weights for GPS tabulated energy."
    )
    parser.add_argument(
        "-r",
        "--root",
        default="hadr05.root",
        help="Input ROOT file (default: hadr05.root)",
    )
    parser.add_argument(
        "-k",
        "--hist",
        default="321",
        help='Histogram object name (default: "321")',
    )
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help='Output path, or "-" for stdout (default: -)',
    )
    parser.add_argument(
        "--center",
        choices=("lin", "log"),
        default="log",
        help="Bin centre: arithmetic mid (lin) or sqrt(low*high) for log bins (default: log)",
    )
    parser.add_argument(
        "--raw-counts",
        action="store_true",
        help="Export histogram counts; default is probabilities summing to 1",
    )
    parser.add_argument(
        "--keep-empty",
        action="store_true",
        help="Include bins with zero weight (default: omit them)",
    )
    parser.add_argument(
        "--energy-scale",
        type=float,
        default=1.0,
        help="Multiply energies by this factor (e.g. 1e-6 if axis is eV but GPS needs MeV) (default: 1)",
    )
    parser.add_argument(
        "--gps-dat",
        action="store_true",
        help="Write only Energy Weight lines (no # comments) for /gps/hist/file",
    )

    args = parser.parse_args()

    if args.gps_dat and args.energy_scale == 1.0:
        args.energy_scale = 1e-6

    path = args.root
    key = args.hist

    try:
        with uproot.open(path) as rf:
            if key not in rf:
                print(f"Key '{key}' not in file. Top-level keys:", file=sys.stderr)
                print(list(rf.keys()), file=sys.stderr)
                return 1
            obj = rf[key]
    except FileNotFoundError:
        print(f"File not found: {path}", file=sys.stderr)
        return 1

    class_name = obj.classname
    if class_name not in ("TH1D", "TH1F", "TH1I"):
        print(f"Object '{key}' is {class_name}, expected TH1D/TH1F/TH1I.", file=sys.stderr)
        return 1

    counts, edges = obj.to_numpy(flow=False)
    counts = np.asarray(counts, dtype=np.float64)
    edges = np.asarray(edges, dtype=np.float64)

    if args.center == "lin":
        energy = 0.5 * (edges[:-1] + edges[1:])
    else:
        energy = np.sqrt(np.maximum(edges[:-1] * edges[1:], np.nextafter(0.0, 1.0)))

    energy *= args.energy_scale

    if args.raw_counts:
        weights = counts.copy()
    else:
        total = counts.sum()
        if total <= 0.0:
            print("Histogram has zero integral; nothing to export.", file=sys.stderr)
            return 1
        weights = counts / total

    lines_out = []
    if not args.gps_dat:
        lines_out.append("# Energy Weight_or_probability")
        try:
            title = obj.member("fTitle")
            if isinstance(title, bytes):
                title = title.decode("utf-8", errors="replace")
            title_str = str(title).strip()
        except Exception:
            title_str = ""
        lines_out.append("# ROOT file: {}  histogram: {} {}".format(path, key, title_str))
        lines_out.append(
            "# Units: energy multiplied by {}; weights are {}; {} bins.".format(
                args.energy_scale,
                "raw bin counts" if args.raw_counts else "normalized (sum=1 over full hist)",
                "all" if args.keep_empty else "non-zero only",
            )
        )

    for e, w in zip(energy, weights):
        if not args.keep_empty and w <= 0.0:
            continue
        lines_out.append(f"{e:.17g} {w:.17g}")

    text = "\n".join(lines_out) + "\n"

    if args.output == "-":
        sys.stdout.write(text)
    else:
        with open(args.output, "w", encoding="utf-8") as out:
            out.write(text)

    return 0


if __name__ == "__main__":
    sys.exit(main())
