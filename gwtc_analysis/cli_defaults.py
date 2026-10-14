"""Constants of the command line, without numpy, pandas or the analysis modules.

`gwtc_analysis -h` builds the whole parser; importing the modes for their defaults would load numpy, pandas,
healpy, ligo.skymap and the waveform packages first, which takes seconds on a cold start. The modes import
their defaults from here, so each is defined once; the analysis modules are imported only when a mode runs.
"""
from __future__ import annotations

from . import catalog_registry as _reg

# rates: the LVK cumulative search-sensitivity releases (real-injection mixture with Cartesian spins)
RATES_SENSITIVITY_RELEASES = {k: (r.record, r.label) for k, r in _reg.SENSITIVITY_RELEASES.items()}
RATES_DEFAULT_RELEASE = _reg.DEFAULT_RATES_RELEASE

# hubble_constant: the releases with semi-analytic O1+O2 injections, so that the whole catalog since O1 can be used
H0_SENSITIVITY_RELEASES = {
    k: dict(record=r.record, label=r.semi_label, file_re=r.semi_file_re, runs=r.runs,
            catalogs=_reg.gwosc_lists(_reg.release_catalogs(k)), published=r.published or None)
    for k, r in sorted(_reg.SENSITIVITY_RELEASES.items())
}
H0_DEFAULT_RELEASE = _reg.DEFAULT_H0_RELEASE
# GWTC-4.0 cosmology: GW231123 (PE systematics) and GW200105 (NSBH) are not used
H0_DEFAULT_EXCLUDE = ("GW231123_135430", "GW200105_162426")
H0_STAGES = ("prepare", "sample", "combine", "reweight", "report")
H0_MASS_MODELS = {"plp": "Power Law + Peak", "mltp": "Multi Peak", "fullpop4": "FullPop-4.0"}
H0_DEFAULT_MASS_MODEL = "plp"

# counterpart and hubble_constant --method bright: the events of counterpart.COUNTERPARTS
COUNTERPART_EVENTS = ("GW170817", "GW190521")
# neutron_star_eos: the events of ns_eos.EVENTS
NS_EOS_EVENTS = ("GW170817", "GW190425")

# event_selection: source classes
NS_MAX_MASS = 3.0                # as catalog_statistics --ns-threshold and the binary types
MASS_GAP = (3.0, 5.0)            # lower mass gap (M☉)
PISN_GAP_MIN = 50.0              # lower edge of the pair-instability gap (M☉); uncertain, ~45-65 in the literature
SOURCE_PRESETS = ("neutron-stars", "mass-gap", "hierarchical")

# lens_pairs: the O4a events of the LVK GWTC-4.0 lensing search
LENS_DEFAULT_CATALOGS = ("GWTC-4",)

ALLOWED_CATALOGS = _reg.allowed_catalogs()
