"""Offline tests for supplementary PSDs and label-to-approximant parsing.

The supplementary release is replaced by a pre-filled cache or a stub, so no
network access is needed.
"""
from __future__ import annotations

import json

import numpy as np
import pytest
from pesummary.gw.file.psd import PSDDict

from gwtc_analysis import pe_supplements as ps
from gwtc_analysis.gwpe_utils import approximant_from_label_rhs, label_has_psd


class _PEData:
    """Minimal stand-in for a pesummary result with empty PSD groups."""

    def __init__(self, labels, with_psd=()):
        self.labels = list(labels)
        self.psd = {lab: PSDDict({}) for lab in labels}
        for lab in with_psd:
            self.psd[lab] = PSDDict({"L1": np.column_stack([np.arange(1.0, 5.0), np.ones(4)])})


PSDS = {"L1": np.column_stack([np.arange(20.0, 30.0), np.full(10, 1e-46)])}


@pytest.mark.parametrize("name, event", [
    ("GW230529_181500", "GW230529_181500"),
    ("GW230529", "GW230529_181500"),
    ("GW190425_081805-v3", "GW190425_081805"),
    ("GW200105_162426", "GW200105_162426"),
    ("GW150914", None),
])
def test_registry_lookup(name, event):
    """Events are found by full or short name; unregistered events return None."""
    spec = ps.get_psd_supplement(name)
    assert (spec.event if spec else None) == event


def test_registry_sources_are_public():
    """Every supplement points to a public DCC or Zenodo file."""
    for spec in ps.PSD_SUPPLEMENTS.values():
        assert spec.url.startswith(("https://dcc.ligo.org/public/", "https://zenodo.org/records/"))


def test_fill_missing_psds_fills_only_empty_labels(monkeypatch):
    """Labels without a PSD get the supplementary PSDs; others are left alone."""
    monkeypatch.setattr(ps, "load_supplementary_psds", lambda spec, **kw: PSDS)
    pe = _PEData(["C00:A", "C00:B", "C00:Mixed"], with_psd=["C00:B"])
    original_b = pe.psd["C00:B"]
    logs = []
    filled = ps.fill_missing_psds(pe, "GW230529_181500", log_cb=logs.append)
    assert filled == ["C00:A", "C00:Mixed"]
    assert all(label_has_psd(pe, lab) for lab in pe.labels)
    assert pe.psd["C00:B"] is original_b
    assert np.array_equal(pe.psd["C00:A"]["L1"], PSDS["L1"])
    assert any("Zenodo 10845779" in m for m in logs)


def test_fill_missing_psds_warns_without_supplement():
    """An unregistered event with empty PSDs is reported, not modified."""
    pe = _PEData(["C00:A"])
    logs = []
    assert ps.fill_missing_psds(pe, "GW150914_095045", log_cb=logs.append) == []
    assert not label_has_psd(pe, "C00:A")
    assert any("no supplementary" in m for m in logs)


def test_fill_missing_psds_noop_when_complete(monkeypatch):
    """Files that have their PSDs never trigger a download."""
    def boom(spec, **kw):
        raise AssertionError("should not fetch")

    monkeypatch.setattr(ps, "load_supplementary_psds", boom)
    pe = _PEData(["C00:A"], with_psd=["C00:A"])
    assert ps.fill_missing_psds(pe, "GW230529_181500") == []


def test_supplementary_psds_come_from_cache(tmp_path, monkeypatch):
    """A cached extraction for the same source URL is reused without network."""
    spec = ps.PSD_SUPPLEMENTS["GW190425_081805"]
    np.savez(tmp_path / f"{spec.event}_psds.npz", **PSDS)
    (tmp_path / f"{spec.event}_psds.json").write_text(json.dumps({"url": spec.url}))

    import fsspec

    def no_network(*a, **k):
        raise AssertionError("should not fetch")

    monkeypatch.setattr(fsspec, "filesystem", no_network)
    out = ps.load_supplementary_psds(spec, cache_dir=tmp_path)
    assert np.array_equal(out["L1"], PSDS["L1"])


@pytest.mark.parametrize("rhs, approximant", [
    ("IMRPhenomPv2_NRTidal:HighSpin", "IMRPhenomPv2_NRTidal"),
    ("IMRPhenomXPHM:LowSpinSecondary", "IMRPhenomXPHM"),
    ("IMRPhenomPv2-NRTidalv2", "IMRPhenomPv2_NRTidalv2"),
    ("IMRPhenomXPHM-SpinTaylor", "IMRPhenomXPHM"),
    ("SEOBNRv5PHM:HighSpin", "SEOBNRv5PHM"),
    ("SEOBNRv4_ROM_NRTidalv2_NSBH", "SEOBNRv4_ROM_NRTidalv2_NSBH"),
])
def test_approximant_from_label_rhs(rhs, approximant):
    """PE label suffixes are stripped to a waveform approximant name."""
    assert approximant_from_label_rhs(rhs) == approximant


def test_clean_psd_drops_repeated_copy_and_restores_grid():
    """A table stored twice with 6-digit frequencies becomes one strictly increasing grid."""
    df = 1.0 / 128
    f = np.arange(0, 2048, df)
    table = np.column_stack([np.array([float(f"{x:.6g}") for x in f]), np.linspace(1, 2, len(f))])
    doubled = np.vstack([table, table])
    out = ps._clean_psd(doubled)
    assert out.shape == table.shape
    assert np.all(np.diff(out[:, 0]) > 0)
    assert np.array_equal(out[:, 0], f)
    assert np.array_equal(out[:, 1], table[:, 1])
