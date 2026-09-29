# parameters_estimation

Everything about one event: its posterior distributions, the whitened strain of each detector with
the best-fit waveform overlaid, a q-transform (time–frequency map), and the matched-filter SNR.

```bash
python -m gwtc_analysis.cli parameters_estimation --src-name GW231223_032836
python -m gwtc_analysis.cli parameters_estimation --src-name GW150914_095045 \
    --pe-vars chi_eff chi_p --pe-pairs mass_1_source:mass_2_source
```

## Posterior samples

The event's PESummary `PEDataRelease` file is read from `--data-repo`. It holds several analyses,
called **labels** (`C00:Mixed`, `C00:IMRPhenomXPHM-SpinTaylor`, `C00:SEOBNRv5PHM`, …), each with its
posterior samples, priors, PSDs, calibration envelopes and configuration. `--pe-vars` adds 1-D
posteriors and `--pe-pairs` 2-D ones (`x:y`).

### Choosing `--pe-label` and `--waveform-engine`

The mode distinguishes **which label is read** (posteriors and metadata) from **which waveform engine
synthesizes** the signal of the strain overlay.

1. **`--pe-label` given:** that label is used for the posteriors, and as the source of the PSDs and
   of the maximum-likelihood parameters of the overlay.
2. **Only `--waveform-engine` given:** the label whose name best matches the engine (substring
   match, no hard-coded mapping) is used for everything. `--waveform-engine IMRPhenomXPHM` selects
   `C00:IMRPhenomXPHM-SpinTaylor` if present.
3. **Neither:** the `Mixed` label.

If the requested engine cannot be instantiated (for instance outside its parameter range), the mode
logs a warning, falls back to another engine when possible, and reports both the requested and the
used engine in the logs and plot titles.

## Strain overlay and q-transform

The strain around the event is downloaded from GWOSC, whitened with the PSD of the label, band-passed,
and the maximum-likelihood waveform projected on each detector is overlaid. The time windows and
frequency bands have shared defaults and per-product overrides:

| Applies to | Options |
|---|---|
| both products | `--start`, `--stop` (seconds before and after the merger), `--fmin`, `--fmax` |
| overlay only | `--overlay-start`, `--overlay-stop`, `--overlay-fmin`, `--overlay-fmax` |
| q-transform only | `--q-start`, `--q-stop`, `--q-fmin`, `--q-fmax`, `--q-fscale {linear,log}` |

When the posterior is BNS-like (median chirp mass below 5 M☉), the windows not set explicitly switch
to a BNS profile: longer windows and a wider frequency range.

The waveform is aligned to the data in time and phase with the matched filter below, so that the
overlay stays coherent even when the stored maximum-likelihood extrinsic parameters are approximate.

## Matched-filter SNR

For each detector, the maximum-likelihood projected waveform is matched-filtered against the strain
(PyCBC; Allen et al. 2012 [\[51\]](../references.md#ref-51),
Usman et al. 2016 [\[52\]](../references.md#ref-52)). The resulting |ρ(t)| should peak at the
coalescence time, at about the detector's recovered SNR.

- The strain is conditioned the standard PyCBC way (high-pass at 15 Hz, resampled to 2048 Hz, edges
  cropped of filter transients), so that the off-source |ρ(t)| has unit-scale RMS (~0.7). A guard
  warns if it strays from that range.
- Loud glitches are gated out before filtering, and the SNR is only reported where the matched
  filter is valid.
- **Short (BBH-like) signals only.** A single template cannot coherently recover a long BNS inspiral:
  over thousands of cycles, small parameter and phase differences accumulate and the SNR is lost
  (that requires a template bank). When the template is longer than the conditioned data, the
  matched-filter SNR is skipped with a warning; the other products are still made.

## Missing PSDs

A few official PE files have no PSDs. They are then taken from public supplementary releases: see
[Known issues in the public releases](../known-data-issues.md#missing-psds-and-calibration-envelopes).

## GW170817

GW170817 has no catalog PE file. The mode transparently uses the bundle rebuilt from public GWTC-1
products: see [build_unofficial_pe](unofficial-pe.md).

All options: [CLI reference](../cli-reference.md#parameters_estimation).
