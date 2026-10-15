# continuous_waves

Continuous waves (CW) are long-lived, nearly monochromatic signals: from spinning neutron stars with a small
deformation (a "mountain", of ellipticity ε), or from clouds of ultralight bosons around spinning black holes, which
annihilate into gravitons at twice the boson's Compton frequency. **None has been detected.** The mode gathers the
upper limits of the LVK data releases and shows them as the papers do, with the boson mass found by
[scalar_cloud](scalar-cloud.md) for GW190728 drawn on the boson-cloud limits.

```bash
gwtc_analysis continuous_waves                                  # downloads 3 MB once
gwtc_analysis continuous_waves --boson-mass 2e-13 3e-13 4e-13 --boson-label "my band"
```

## Data

| Search | Paper | Data release (LIGO DCC, public) | What the mode reads |
|---|---|---|---|
| 45 known pulsars, O4a | [\[105\]](../references.md#ref-105) | [LIGO-G2402400](https://dcc.ligo.org/LIGO-G2402400/public) | the 95% upper limits on h₀ (Bayesian pipeline); the spin-down limits, distances and ellipticity limits from the paper's table (its arXiv LaTeX source); the narrowband upper limits of 16 of them (and 2 post-glitch analyses) |
| all-sky, isolated neutron stars, O4a | [\[106\]](../references.md#ref-106) | [LIGO-P2500416](https://dcc.ligo.org/LIGO-P2500416/public) | h₀ upper limits against frequency, 20–2000 Hz: FrequencyHough, PowerFlux (population average), SOAP |
| all-sky, boson clouds, O3 | [\[107\]](../references.md#ref-107) | [LIGO-P2100343](https://dcc.ligo.org/LIGO-P2100343/public) | h₀ upper limits per 1 Hz band (20–610 Hz); the excluded black-hole and boson masses (a distance of 1 or 15 kpc, an initial spin of 0.5 or 0.9, an age of the cloud); the largest excluded distance against the boson mass |

## Results

| Result | Value | Published |
|---|---|---|
| known pulsars with h₀ below the spin-down limit | 26 with the Bayesian pipeline; 29 with any of the three | 29 |
| narrowband search (the frequency allowed to wander around twice the rotation frequency) below the spin-down limit | 12 of 16; the Crab at 3% of it | the same (Crab 0.03) |
| lowest h₀ upper limit | 6.4 × 10⁻²⁷, J0537−6910 (LMC, 62 Hz rotation) | the same |
| lowest ellipticity upper limit | 8.8 × 10⁻⁹, J0437−4715 (millisecond pulsar at 160 pc) | the same |
| all-sky O4a, lowest population-averaged h₀ | 9.5 × 10⁻²⁶ at 401 Hz (single 1/16 Hz band); 9.8 × 10⁻²⁶ near 298 Hz | 9.7 × 10⁻²⁶ near 290 Hz: the curve is flat within a few per cent from about 250 to 400 Hz |

![Known pulsars and the all-sky sensitivity](../img/modes/cw_known_pulsars.png)

*Each known pulsar's 95% upper limit (blue) against its spin-down limit (grey triangle: the strain if all the
rotational energy it loses went into gravitational waves); the stars beat it, the Crab by a factor of 130. The open diamonds are the narrowband search, less sensitive than the
targeted one but robust to a mismatch between the GW and the radio or X-ray timing (the Crab at 3% of its spin-down
limit, about 0.1% of the spin-down energy in gravitational waves; the later O4 narrowband search of 34 pulsars brings
it to 2%, ApJ 1005, 221 (2026)). The
curves are the all-sky O4a sensitivity to unknown neutron stars, which must search every direction and spin-down.*

## Boson clouds and GW190728

A boson of mass \(m\) grows a cloud by superradiance around a black hole whose size matches its Compton
wavelength, \(\alpha = G M m / (\hbar c) \sim 0.1\), and the cloud then annihilates into a continuous wave at

\[
f = \frac{2 m c^2}{h} = 483.6\ \text{Hz}\ \frac{m}{10^{-12}\ \text{eV}} .
\]

The scalar cloud fitted to GW190728 ([scalar_cloud](scalar-cloud.md#gw190728_064510)) has
\(m_\phi = 0.89\) (0.70–1.17) × 10⁻¹² eV with the superradiance prior (0.73 (0.56–1.14) × 10⁻¹² eV agnostic), so
the same particle would make clouds around Galactic black holes radiate at **340–570 Hz**, inside the band of the
O3 boson-cloud search, whose h₀ upper limits there are 1.2–1.9 × 10⁻²⁵.

![Excluded black-hole and boson masses](../img/modes/cw_boson_exclusion.png)

*O3 exclusions for a black hole of initial spin 0.9 at 1 and 15 kpc, for clouds of different ages, with the
GW190728 band. At these boson masses only light black holes are concerned, 8–28 M☉ at 1 kpc and 12–28 M☉ at
15 kpc for a young cloud: heavier ones are not superradiant for so heavy a boson. That is the mass range of
GW190728's own components (12 and 8 M☉).*

![Largest excluded distance](../img/modes/cw_boson_distance.png)

*In the GW190728 band, clouds younger than 10³ yr are excluded out to 16–32 kpc, 10⁴ yr 13–15 kpc,
10⁵ yr 5–8 kpc, 10⁶ yr 2–3 kpc, 10⁷ yr about 1 kpc, for black holes up to 50 or 100 M☉.*

**What it means.** The continuous-wave limits do not exclude the GW190728 boson: they say that no young cloud
around a light, fast-spinning black hole sits within a few kpc of the Earth. A boson of that mass is compatible with
them as long as such systems are rare or old; turning the limits into a constraint on the boson needs a model of the
Galactic black holes (their number, masses, spins and ages). The two probes are complementary: the merger tests a
cloud carried by a binary, the continuous wave a cloud around an isolated black hole. Isolated Galactic black holes
found by microlensing (one so far, many more expected with Rubin) come with a mass and a distance, the inputs of a directed
search for their clouds.

## Outputs

- `--out-report` (HTML): the summary, the known pulsars, the all-sky limits and the boson clouds, with the plots;
- `--out-summary` (TSV): the summary; `<name>.pulsars.tsv`, `<name>.narrowband.tsv` and `<name>.boson_reach.tsv`:
  the pulsar upper limits, the narrowband limits and the distance reach for boson clouds;
- `--plots-dir`: `cw_known_pulsars.png`, `cw_boson_exclusion.png`, `cw_boson_distance.png`.

## Options

| Option | Default | Meaning |
|---|---|---|
| `--boson-mass LOW MID HIGH` | GW190728, superradiance prior: 0.70e-12 0.89e-12 1.17e-12 | the band drawn on the boson plots, eV |
| `--boson-label` | `GW190728 (scalar_cloud)`, or `given band` with `--boson-mass` | its label |
| `--refresh` | off | download the data releases again (they are cached in `~/.cache_gwtc_analysis/continuous_waves`, or `$GWTC_CW_CACHE`) |
| `--out-report`, `--out-summary`, `--plots-dir` | `continuous_waves.html`, `continuous_waves.tsv`, `continuous_waves_plots` | outputs; the summary has three companions, `<name>.pulsars.tsv` (the 45 pulsars), `<name>.narrowband.tsv` and `<name>.boson_reach.tsv` (the largest excluded distance in the band) |
| `--data-repo {zenodo,s3}` | `zenodo` | where the LVK upper-limit files are downloaded from: their public sources, or `s3`, the MinIO mirror (faster), when it holds the same file ([S3 mirror of the public files](../data-sources.md#s3-mirror-of-the-public-files)) |
