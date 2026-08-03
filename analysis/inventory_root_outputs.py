#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
import uproot

FILES = [
    "build/bare-detector-2_2MeV-gamma.root",
    "build/bare-detector-662keV-gamma.root",
    "build/bare-detector-DD-neutron.root",
    "build/bare-detector-DT-neutron.root",
    "build/bare-detector-thermal-fast-noB4C.root",
    "build/bare-detector-thermal-fast-withB4C.root",
    "build/bare-detector-thermal-only.root",
]


def clean(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def axis_summary(obj: Any) -> dict[str, Any]:
    out: dict[str, Any] = {}
    try:
        axes = obj.axes
        out["ndim"] = len(axes)
        for i, axis in enumerate(axes):
            edges = np.asarray(axis.edges())
            out[f"axis{i}_bins"] = len(edges) - 1
            out[f"axis{i}_min"] = float(edges[0])
            out[f"axis{i}_max"] = float(edges[-1])
            out[f"axis{i}_label"] = clean(getattr(axis, "label", ""))
    except Exception:
        pass
    return out


def inspect_file(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with uproot.open(path) as root:
        for key, obj in root.items(recursive=True):
            classname = clean(getattr(obj, "classname", type(obj).__name__))
            row: dict[str, Any] = {
                "file": path.name,
                "key": key,
                "classname": classname,
                "title": clean(getattr(obj, "title", "")),
            }
            row.update(axis_summary(obj))
            try:
                values = np.asarray(obj.values(flow=False))
                row["sum"] = float(np.nansum(values))
                row["nonzero_bins"] = int(np.count_nonzero(values))
                row["max_bin"] = float(np.nanmax(values)) if values.size else 0.0
            except Exception:
                row["sum"] = ""
                row["nonzero_bins"] = ""
                row["max_bin"] = ""
            try:
                row["num_entries"] = int(obj.num_entries)
                row["branches"] = ";".join(obj.keys())
            except Exception:
                row["num_entries"] = ""
                row["branches"] = ""
            rows.append(row)
    return rows


def main() -> None:
    outdir = Path("root_inventory")
    outdir.mkdir(exist_ok=True)
    rows: list[dict[str, Any]] = []
    missing: list[str] = []
    for name in FILES:
        path = Path(name)
        if not path.exists():
            missing.append(name)
            continue
        rows.extend(inspect_file(path))

    fields = sorted({key for row in rows for key in row})
    with (outdir / "root_inventory.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "files_requested": FILES,
        "files_missing": missing,
        "object_count": len(rows),
        "objects_by_file": {
            filename: sum(1 for row in rows if row["file"] == Path(filename).name)
            for filename in FILES
        },
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # Human-readable filtered inventory for rapid review.
    terms = ("LET", "Eion", "secondary", "Secondary", "depth", "Depth", "track", "Track", "fraction")
    with (outdir / "relevant_objects.txt").open("w", encoding="utf-8") as f:
        for row in rows:
            haystack = f"{row.get('key', '')} {row.get('title', '')}"
            if any(term in haystack for term in terms):
                f.write(
                    f"{row['file']}\t{row['key']}\t{row['classname']}\t"
                    f"{row.get('title', '')}\tsum={row.get('sum', '')}\n"
                )


if __name__ == "__main__":
    main()
