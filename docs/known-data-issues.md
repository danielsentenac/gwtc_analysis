# Known issues in the public releases

Working with the releases exposed a few events whose files lack products that the release
descriptions announce. `gwtc_analysis` works around them, and this page records what is missing,
where it can be found, and how the tool handles it.

## Missing PSDs and calibration envelopes

Each PESummary `PEDataRelease` file is described as including the noise power spectral densities
(PSDs) and the calibration uncertainty envelopes of every analysis. For three low-mass events, the
`psds` and `calibration_envelope` groups exist but are empty in all analyses:

| Event | Release(s) affected | Public source of the missing data | Caveat |
|---|---|---|---|
| GW230529_181500 | GWTC-4.0 v2 and v3, GWTC-4.1 | Discovery release, [Zenodo 10845779](https://zenodo.org/records/10845779) ([LIGO-P2300352](https://dcc.ligo.org/LIGO-P2300352/public)): L1 PSD in all 15 runs, calibration envelopes in the single-waveform runs | The L1 PSD is identical in all 15 runs |
| GW190425_081805 | GWTC-2.1 v2 (the only version with this event) | Discovery release [LIGO-P2000026](https://dcc.ligo.org/LIGO-P2000026/public), also in the GWTC-2 release [LIGO-P2000223](https://dcc.ligo.org/LIGO-P2000223/public) | PSDs of the earlier LALInference analysis, not those of the GWTC-2.1 bilby runs, which are not public |
| GW200105_162426 | GWTC-3 v1 and v2 | **Fixed in GWTC-3 v3** ([Zenodo 22685054](https://zenodo.org/records/22685054)); also the discovery release [LIGO-P2100143](https://dcc.ligo.org/LIGO-P2100143/public) | The discovery PSDs agree with the v3 ones to ~2% (median), differing mostly on lines |

A scan of every PE file of GWTC-4.0 (86 events), GWTC-4.1 (88) and GWTC-5.0 (104) found
GW230529 to be the only event with empty PSDs in those releases. In the affected files the
configuration metadata still point to the original input files on the analysts' cluster accounts: the
inputs were apparently not embedded when the results were packaged.

**How the tool handles it.** The whitened strain overlay of `parameters_estimation` needs a PSD. If a
PE label has none, the event is looked up in a registry of public supplementary releases
(`gwtc_analysis/pe_supplements.py`), usually the data release of its discovery paper. The missing
PSD (or skymap) is taken from there and attached to every label that lacks one. Only the PSD group is
read, over HTTP range requests where the server allows it, and the result is cached in
`~/.gwcache/psd_supplements`. The log names the source and its caveat. Events with missing PSDs and
no registered supplement are reported, and the overlay then whitens with a PSD estimated from the
strain.

## GW200105: only in the marginal list

The NSBH event GW200105_162426 is listed by GWOSC only in `GWTC-3-marginal` (p_astro = 0.36,
FAR = 0.2 per year), not in `GWTC-3-confident` nor in the cumulative list of confident events. Yet
the GWTC-3 PE release presents it together with the confident events ("plus GW200105_162426, which is
a clear outlier from the noise background"), and it is one of the two NSBH detections of its
discovery paper [\[44\]](references.md#ref-44).

A tool that builds its event list from the confident lists therefore misses GW200105, while its
sibling GW200115 is there. The `rates` and `hubble_constant` modes read the marginal lists too and
select events by false-alarm rate, as the LVK population analyses do [\[12\]](references.md#ref-12) [\[13\]](references.md#ref-13). (The `hubble_constant` mode then
leaves GW200105 out, as the GWTC-4.0 cosmology analysis [\[29\]](references.md#ref-29) does, because it is an NSBH.)

## Rounded false-alarm rates

The GWOSC event lists give, for each event, its lowest false-alarm rate over the search pipelines (checked
against the per-pipeline values of the GWOSC v2 API for 206 of 208 events of GWTC-2.1, GWTC-3 and
GWTC-4.0), **rounded to two decimals**. The LVK analyses cut at full precision, so an event listed at
exactly a threshold may be inside it: GW191127_050227 is listed at FAR 0.25 per year (its GstLAL value)
and is part of the GWTC-4.0 and GWTC-5.0 cosmology samples (FAR < 0.25), as is GW240824_205609 for
GWTC-5.0. `rates` and `hubble_constant` therefore compare the published values inclusively
(FAR ≤ threshold), which reproduces the event selections of those papers exactly (137 and 231 BBHs);
the injections, whose FARs are at full precision, keep the strict cut.

## Version changes without a changelog

The GWTC-3 PE release v3 (September 2026) fixed the GW200105 PSDs, but neither the Zenodo record nor
the GWOSC GWTC-3 page says what changed. Because the tool uses the latest version of each release by
default and can select any older one ([Zenodo release versions](data-sources.md#zenodo-release-versions)),
results can always be traced to a given version.
