"""Offline tests for matching catalog event ids to skymap event names."""
from __future__ import annotations

import pytest

from gwtc_analysis.gw_stat import skymap_event_key

KNOWN = {"GW150914_095045", "GW170817_124104", "GW190412_053044", "GW190521_030229", "GW190521_074359"}


@pytest.mark.parametrize("ev_id, key", [
    ("GW190412_053044-v3", "GW190412_053044"),   # full name: unchanged behaviour
    ("GW150914-v4", "GW150914_095045"),          # GWTC-1 short name in GWOSC
    ("GW150914", "GW150914_095045"),
    ("GW190521-v3", None),                       # ambiguous short name: two events that day
    ("GW151226-v3", None),                       # no skymap for it
    ("", None),
    (None, None),
])
def test_skymap_event_key(ev_id, key):
    """Full names pass through; a short name maps to its unique full-name skymap."""
    assert skymap_event_key(ev_id, KNOWN) == key


def test_short_name_skymap_kept_as_is():
    """A skymap indexed under a short name is matched directly."""
    assert skymap_event_key("GW170817-v3", {"GW170817"}) == "GW170817"
