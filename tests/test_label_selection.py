"""Default PE label of parameters_estimation (select_label): Mixed for the samples, and for the strain a
PSD-capable label, the default engine's when the Mixed label has no PSD (GWTC-4.0 and later)."""
from __future__ import annotations

import pytest

from gwtc_analysis import gwpe_utils as gu


class FakePE:
    def __init__(self, labels, no_psd=()):
        self.labels = list(labels)
        self.samples_dict = {lab: {} for lab in labels}
        self.approximant = [lab.split(":", 1)[1] for lab in labels]
        self.psd = {lab: (None if lab in no_psd else {"H1": [[20.0, 1e-46], [21.0, 1e-46]]}) for lab in labels}


@pytest.fixture(autouse=True)
def _no_remembered_engine(monkeypatch):
    monkeypatch.setattr(gu, "_LAST_REQUESTED_WAVEFORM_ENGINE", None)


GWTC4 = ["C00:IMRPhenomTPHM", "C00:IMRPhenomXO4a", "C00:IMRPhenomXPHM-SpinTaylor", "C00:Mixed", "C00:Mixed+XO4a",
         "C00:NRSur7dq4", "C00:SEOBNRv5PHM"]


def _sel(pe, **kw):
    return gu.select_label(pe, show_labels=False, **kw)


def test_gwtc4_layout_without_psd_in_mixed():
    """GW231123: samples from C00:Mixed (not the first label, IMRPhenomTPHM); strain from IMRPhenomXPHM."""
    pe = FakePE(GWTC4, no_psd=("C00:Mixed", "C00:Mixed+XO4a"))
    assert _sel(pe, require_psd=False) == "C00:Mixed"
    assert _sel(pe, require_psd=True) == "C00:IMRPhenomXPHM-SpinTaylor"


def test_gwtc21_layout_with_psd_in_mixed():
    pe = FakePE(["C01:IMRPhenomXPHM", "C01:Mixed", "C01:SEOBNRv4PHM"])
    assert _sel(pe, require_psd=False) == "C01:Mixed" and _sel(pe, require_psd=True) == "C01:Mixed"


def test_plain_mixed_preferred_over_variants():
    """GW200115 has C01:Mixed and C01:Mixed:NSBH:*: the plain one."""
    pe = FakePE(["C01:IMRPhenomNSBH:HighSpin", "C01:Mixed:NSBH:HighSpin", "C01:Mixed", "C01:Mixed:NSBH:LowSpin"])
    assert _sel(pe, require_psd=False) == "C01:Mixed"


def test_explicit_choices_unchanged():
    pe = FakePE(GWTC4, no_psd=("C00:Mixed", "C00:Mixed+XO4a"))
    assert _sel(pe, pe_label="C00:NRSur7dq4") == "C00:NRSur7dq4"
    assert _sel(pe, waveform_engine="SEOBNRv5PHM", require_psd=True) == "C00:SEOBNRv5PHM"


def test_no_mixed_label():
    """GW170817 bundle: no Mixed label, the first one as before; no PSD at all: the first label."""
    pe = FakePE(["C02:IMRPhenomPv2_NRTidal-HighSpin", "C02:IMRPhenomPv2_NRTidal-LowSpin"])
    assert _sel(pe) == "C02:IMRPhenomPv2_NRTidal-HighSpin"
    assert _sel(FakePE(["C00:A", "C00:B"], no_psd=("C00:A", "C00:B")), require_psd=True) == "C00:A"
