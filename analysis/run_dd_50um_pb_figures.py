#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

import generate_dd_50um_pb_figures as figures

_original_step_level_depth_data = figures.step_level_depth_data


def step_level_depth_data_with_empty_check(root):
    """Use StepLET only when it contains nonzero in-layer ionizing deposition."""
    raw = _original_step_level_depth_data(root)
    if raw is None:
        return None
    depth, category, energy, method = raw
    valid = (
        np.isfinite(depth)
        & np.isfinite(energy)
        & (energy > 0)
        & (depth >= 0)
        & (depth <= 50)
    )
    if not np.any(valid) or float(np.sum(energy[valid])) <= 0:
        return None
    return raw


figures.step_level_depth_data = step_level_depth_data_with_empty_check
figures.main()
