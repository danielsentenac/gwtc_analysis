"""Offline test of the skymap lightening done before pesummary reads a PE file."""
from __future__ import annotations

import h5py
import numpy as np
import pytest

from gwtc_analysis.parameters_estimation import lighten_skymaps


def test_lighten_skymaps_sums_nested_blocks(tmp_path):
    rng = np.random.default_rng(0)
    big = rng.random(12 * 16 ** 2); big /= big.sum()
    src = tmp_path / "pe.h5"
    with h5py.File(src, "w") as f:
        f.attrs["version"] = "x"
        f["C00:A/posterior_samples"] = np.arange(5.0)
        f["C00:A/skymap/data"] = big
        f["C00:A/skymap/meta_data/nest"] = [b"True"]
        f["C00:B/skymap/data"] = np.full(12 * 4 ** 2, 1 / 192)
    lite = lighten_skymaps(src, max_nside=4)
    assert lite != src and lighten_skymaps(src, max_nside=16) == src
    with h5py.File(lite) as f:
        d = f["C00:A/skymap/data"][()]
        assert d.shape == (192,) and d.sum() == pytest.approx(1.0)
        assert d[0] == pytest.approx(big[:16].sum())
        assert f["C00:A/posterior_samples"][()].tolist() == [0, 1, 2, 3, 4] and f.attrs["version"] == "x"
        assert f["C00:A/skymap/meta_data/nest"][0] == b"True" and f["C00:B/skymap/data"].shape == (192,)
    assert lighten_skymaps(src, max_nside=4) == lite          # reused
