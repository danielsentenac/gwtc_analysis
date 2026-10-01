"""Catalog-level source classes and remnants, from the GWOSC medians and 90% intervals.

- Remnants and energetics: the radiated energy E_rad = (M_total - M_final) c^2 and the radiated fraction,
  from the GWOSC medians, and an estimate of the final spin from the mass ratio and chi_eff with the
  aligned-spin fit of Rezzolla et al. 2008 (arXiv:0712.3541). The final spin of the PE files (and the peak
  luminosity, which GWOSC does not list) are read with `parameters_estimation`.
- Event classes for `event_selection --preset`:
    neutron-stars  a component below the maximum neutron-star mass (median);
    mass-gap       a component in the lower mass gap between neutron stars and black holes (median);
    hierarchical   a primary in the pair-instability mass gap, or a confidently negative chi_eff: the
                   signatures of black holes formed in earlier mergers, in dense environments.

The GWOSC values are medians of the marginal posteriors: the difference of two medians is not the median of
the difference, and a median in a range is not a probability of being in it. The 90% intervals are given
alongside, and the flags ``*_confident`` require the whole interval in the range.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

MSUN_C2_ERG = 1.7877e54          # M☉ c² in erg
NS_MAX_MASS = 3.0                # as catalog_statistics --ns-threshold and the binary types
MASS_GAP = (3.0, 5.0)            # lower mass gap (M☉)
PISN_GAP_MIN = 50.0              # lower edge of the pair-instability gap (M☉); uncertain, ~45-65 in the literature
PRESETS = ("neutron-stars", "mass-gap", "hierarchical")


# ---------------------------------------------------------------------------
# remnants
# ---------------------------------------------------------------------------
def final_spin_aligned(eta, chi):
    """Final spin of a black-hole merger with aligned spins (Rezzolla et al. 2008, arXiv:0712.3541):

    a_f = a + s4 a² η + s5 a η² + t0 a η + 2√3 η + t2 η² + t3 η³,

    with a the spin of both black holes (here chi_eff). 0.686 for equal masses without spin; a for η -> 0.
    """
    eta, a = np.asarray(eta, float), np.asarray(chi, float)
    s4, s5, t0, t2, t3 = -0.1229, 0.4537, -2.8904, -3.5171, 2.5763
    af = a + s4 * a ** 2 * eta + s5 * a * eta ** 2 + t0 * a * eta + 2 * np.sqrt(3) * eta + t2 * eta ** 2 + t3 * eta ** 3
    return np.clip(af, -1.0, 1.0)


def add_remnant_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Radiated energy (M☉ c² and erg), radiated fraction, and final-spin estimate of BBH-like events."""
    out = df.copy()
    # the GWOSC total mass (median of m1 + m2), else the sum of the medians (GWTC-1 lists no total mass)
    tot = pd.to_numeric(out["total_mass_source"], errors="coerce") if "total_mass_source" in out else \
        pd.Series(np.nan, index=out.index)
    if "total_mass_source_gwosc" in out:
        tot = pd.to_numeric(out["total_mass_source_gwosc"], errors="coerce").fillna(tot)
    fin = pd.to_numeric(out.get("final_mass_source"), errors="coerce")
    e = tot - fin
    out["radiated_energy_msun"] = e.where(e > 0)
    out["radiated_energy_erg"] = out["radiated_energy_msun"] * MSUN_C2_ERG
    out["radiated_fraction"] = out["radiated_energy_msun"] / tot
    m1 = pd.to_numeric(out["mass_1_source"], errors="coerce")
    m2 = pd.to_numeric(out["mass_2_source"], errors="coerce")
    eta = m1 * m2 / (m1 + m2) ** 2
    chi = pd.to_numeric(out.get("chi_eff"), errors="coerce").fillna(0.0)
    af = pd.Series(final_spin_aligned(eta, chi), index=out.index)
    # neutron stars do not form a black hole of this spin (tides, ejecta): BBH only
    bbh = (m2 >= NS_MAX_MASS) if "binary_type" not in out else out["binary_type"].eq("BBH")
    out["final_spin_estimate"] = af.where(bbh & eta.notna())
    return out


# ---------------------------------------------------------------------------
# presets
# ---------------------------------------------------------------------------
def _num(df: pd.DataFrame, col: str) -> pd.Series:
    return pd.to_numeric(df[col], errors="coerce") if col in df else pd.Series(np.nan, index=df.index)


def _in(df, col, lo, hi):
    return (_num(df, col) >= lo) & (_num(df, col) <= hi)


def _interval_in(df, col, lo, hi):
    return (_num(df, f"{col}_lo90") >= lo) & (_num(df, f"{col}_hi90") <= hi)


def apply_preset(df: pd.DataFrame, preset: str, *, ns_max_mass: float = NS_MAX_MASS,
                 mass_gap: tuple[float, float] = MASS_GAP, pisn_gap_min: float = PISN_GAP_MIN
                 ) -> tuple[pd.Series, pd.DataFrame, str]:
    """(mask, extra columns, description) of a preset on a GWOSC event table with *_lo90 / *_hi90 bounds."""
    m1, m2 = _num(df, "mass_1_source"), _num(df, "mass_2_source")
    extra = pd.DataFrame(index=df.index)
    if preset == "neutron-stars":
        mask = m2 < ns_max_mass
        extra["class"] = np.where(m1 < ns_max_mass, "BNS", "NSBH")
        extra["mass_2_source_hi90"] = _num(df, "mass_2_source_hi90")
        extra["m2_below_ns_max_90"] = _num(df, "mass_2_source_hi90") < ns_max_mass
        return mask, extra, (f"a component below {ns_max_mass:g} M☉ (median); m2_below_ns_max_90: the 90% upper bound "
                             "of m2 below it too")
    if preset == "mass-gap":
        lo, hi = mass_gap
        g1, g2 = _in(df, "mass_1_source", lo, hi), _in(df, "mass_2_source", lo, hi)
        extra["gap_component"] = np.select([g1 & g2, g1, g2], ["both", "primary", "secondary"], "")
        extra["gap_confident"] = (g1 & _interval_in(df, "mass_1_source", lo, hi)) | \
                                 (g2 & _interval_in(df, "mass_2_source", lo, hi))
        for c in ("mass_1_source", "mass_2_source"):
            extra[f"{c}_lo90"], extra[f"{c}_hi90"] = _num(df, f"{c}_lo90"), _num(df, f"{c}_hi90")
        return g1 | g2, extra, (f"a component in the lower mass gap {lo:g}–{hi:g} M☉ (median); gap_confident: its "
                                "whole 90% interval in the gap")
    if preset == "hierarchical":
        gap = m1 >= pisn_gap_min
        neg = _num(df, "chi_eff_hi90") < 0
        extra["pisn_gap"] = gap
        extra["pisn_gap_confident"] = _num(df, "mass_1_source_lo90") >= pisn_gap_min
        extra["negative_chi_eff"] = neg
        for c in ("mass_1_source_lo90", "chi_eff", "chi_eff_lo90", "chi_eff_hi90"):
            extra[c] = _num(df, c)
        return gap | neg, extra, (f"primary in the pair-instability gap, m1 ≥ {pisn_gap_min:g} M☉ (median; "
                                  "pisn_gap_confident: the 90% lower bound too), or chi_eff < 0 at 90% "
                                  "(spins anti-aligned: dynamical assembly)")
    raise ValueError(f"unknown preset {preset!r}; choose from {', '.join(PRESETS)}")
