# bright_siren

The Hubble constant from **GW170817 and its host galaxy NGC 4993** (bright siren), alone or combined
with the [spectral siren](hubble-constant.md). It takes seconds: no sampling, no injections.

```bash
gwtc_analysis bright_siren
gwtc_analysis bright_siren --spectral-posterior hubble_constant_run     # combined with the spectral siren
```

## Method

The GW signal gives the luminosity distance \(d\) with no distance ladder; the host galaxy gives the
recession velocity \(v_r\). Removing the galaxy's peculiar velocity \(\langle v_p \rangle\) leaves the
Hubble-flow velocity, and at \(z \approx 0.01\) the Hubble law is linear:

\[
v_H = v_r - \langle v_p \rangle = H_0\, d .
\]

The likelihood follows LVK 2017 [\[27\]](../references.md#ref-27). With Gaussian velocity
measurements, the marginalization over the true peculiar velocity gives one Gaussian of width
\(\sigma = \sqrt{\sigma_r^2 + \sigma_p^2}\). For sources uniform in volume and a detection limited by
the GW signal, the selection term and the volume prior on the distance cancel
(Chen, Fishbach & Holz 2018 [\[72\]](../references.md#ref-72); Mandel, Farr & Gair 2019
[\[35\]](../references.md#ref-35)). The PE samples \(d_i\) were drawn with a \(d^2\) prior, so with a
flat H₀ prior

\[
p(H_0 \mid \text{data}) \propto \sum_i \mathcal{N}\!\left(v_H;\; H_0 d_i,\; \sigma\right).
\]

The distance must be the one **at the counterpart's sky position**: the distance and the sky
position are correlated through the antenna pattern. The GWTC-1 GW170817 samples were produced with
the sky fixed to AT2017gfo, so they are used as they are. Samples that are not are restricted to
those within `--sky-radius` of the counterpart.

| Input | GW170817 default | Source |
|---|---|---|
| PE samples | `C02:IMRPhenomPv2_NRTidal-LowSpin` and `-HighSpin` of the [GW170817 bundle](unofficial-pe.md) | GWTC-1 public release |
| \(v_r\) | 3327 ± 72 km/s (NGC 4993 group, CMB frame) | [\[27\]](../references.md#ref-27) |
| \(\langle v_p \rangle\) | 310 ± 150 km/s | [\[27\]](../references.md#ref-27) |
| \(v_H\) | 3017 ± 166 km/s | |
| H₀ prior | flat, 10–200 km/s/Mpc | that of the spectral siren |

## Result

| Analysis | Maximum a posteriori, 68% interval | Median, 90% interval | \(d_L\) median (90%) |
|---|---|---|---|
| LowSpin | **69.2** (+23.4 / −8.1) | 79.3, 62.5–134.3 | 40.0 Mpc (24.9–47.3) |
| HighSpin | **68.8** (+14.8 / −7.3) | 74.2, 61.8–111.4 | 41.7 Mpc (29.1–47.5) |
| LVK 2017 [\[27\]](../references.md#ref-27) | **70.0** (+12.0 / −8.0) | | 43.8 Mpc (+2.9 / −6.9, 68%) |

![H0 from GW170817 and NGC 4993](../img/modes/bright_siren_GW170817.png)

The maximum a posteriori and the lower side agree with the published value. The upper side is wider
because the public GWTC-1 samples come from a later reanalysis (GWTC-1, IMRPhenomPv2_NRTidal),
whose distance has a longer low-distance tail than the 2017 analysis. That tail comes from the
**distance–inclination degeneracy**: an inclined binary is fainter than a face-on one at the same
distance, so a nearby inclined source and a farther face-on one look alike. Constraining the
inclination with the radio jet of GW170817 narrows H₀ considerably (Hotokezaka et al. 2019
[\[73\]](../references.md#ref-73)); this mode does not use EM inclination constraints.

## Combination with the spectral siren

`--spectral-posterior` takes a `hubble_constant` work directory (its `posterior_reweighted.tsv`, else
`posterior.tsv`) or any TSV with an `H0` column. The two measurements use different events (the
spectral siren uses binary black holes only) and the same flat prior, so their posteriors multiply.
The spectral-siren density is a Gaussian kernel estimate, reflected at the prior bounds.

![Bright siren combined with the spectral siren](../img/modes/bright_siren_combined.png)

With the *Power Law + Peak* spectral siren (GWTC-4.0 setup): H₀ = 71.1 (+22.5 / −7.9) km/s/Mpc
(maximum a posteriori, 68%). GW170817 dominates; the broad spectral siren shifts the posterior slightly upwards.

## Options

| Option | Default | Meaning |
|---|---|---|
| `--src-name` | `GW170817` | event with a registered host galaxy (`COUNTERPARTS` in `bright_siren.py`) |
| `--pe-label` | all the labels, LowSpin first | PE labels; the first one is the main result and the one combined |
| `--pe-file` | the event's bundle | another PE file (PESummary layout) |
| `--v-recession V SIGMA` | 3327 72 | recession velocity of the host, km/s |
| `--v-peculiar V SIGMA` | 310 150 | peculiar velocity of the host, km/s |
| `--sky-radius` | 3° | for samples not fixed to the counterpart |
| `--spectral-posterior` | none | spectral-siren posterior to combine with |
| `--h0-range MIN MAX` | 10 200 | flat prior; keep the spectral siren's for a combination |

The published comparison is shown only with the default velocities.

## Outputs

- `--out-report` (HTML): the result, the inputs and the plots;
- `--out-summary` (TSV): maximum a posteriori, 68.3% highest-density interval, median and 90%
  interval of each analysis, with the distance of each label;
- `<summary>.posterior.tsv`: the posterior densities on an H₀ grid;
- `--plots-dir`: `h0_bright_siren_GW170817.png`, and `h0_combined.png` with `--spectral-posterior`.
