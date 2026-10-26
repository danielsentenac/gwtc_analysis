# hubble_constant: bright siren

The Hubble constant from an event with an **identified host galaxy** (bright siren), with
`hubble_constant --method bright`. It takes seconds: no sampling. The result goes to a work directory that
`--method joint` combines with the [spectral or dark siren](hubble-constant.md#joint-posterior); where the
counterpart sits in the GW posterior (sky, distance, viewing angle) is the [counterpart](counterpart.md) mode.

```bash
gwtc_analysis hubble_constant --method bright                          # GW170817 and NGC 4993
gwtc_analysis hubble_constant --method bright --viewing-angle 20 3     # with a viewing-angle constraint
gwtc_analysis hubble_constant --method bright --event GW190521         # candidate AGN flare, z = 0.438
gwtc_analysis hubble_constant --method joint \
    --inputs hubble_constant_bright h0_dark_plp                        # combined with the dark siren
```

## Method

The GW signal gives the luminosity distance \(d_L\) with no distance ladder; the host gives the
Hubble-flow redshift \(z\). For a nearby galaxy, \(z\) comes from its recession velocity minus its
peculiar velocity, \(c z = v_r - \langle v_p \rangle\). In flat ΛCDM, with \(\Omega_m = 0.3065\) as in
the spectral siren,

\[
d_L = \frac{c}{H_0}\, D(z; \Omega_m).
\]

For sources uniform in comoving volume and source-frame time,
\(p_\text{pop}(z) \propto \frac{dV_c}{dz} \frac{1}{1+z}\), and with a flat H₀ prior
(Chen, Fishbach & Holz 2018 [\[72\]](../references.md#ref-72);
Mandel, Farr & Gair 2019 [\[35\]](../references.md#ref-35)):

\[
p(H_0 \mid \text{data}) \propto \frac{1}{\beta(H_0)}
\int \mathcal{N}(z_\text{obs};\, z,\, \sigma_z)\; p_\text{pop}(z)\; \mathcal{L}_\text{GW}\big(d_L(z, H_0)\big)\, dz .
\]

\(\mathcal{L}_\text{GW}\) is the posterior of \(d_L\) divided by the PE distance prior, read from the
file (\(d_L^2\) up to O3, uniform in comoving volume and source-frame time in O4), or set by `--pe-distance-prior`. Over the PE samples \(d_i\), whose
redshift at a trial H₀ is \(z_i\), the integral is a kernel sum:

\[
\sum_i w_i\, \mathcal{N}(z_\text{obs};\, z_i,\, \sigma)\,
\frac{p_\text{pop}(z_i)}{\pi_\text{PE}(d_i)\; \partial d_L/\partial z\,(z_i)},
\qquad \sigma = \max(\sigma_z, h),
\]

where \(h\) is the kernel width of the \(z_i\) (Silverman's rule). It only matters for a distant host
with a spectroscopic redshift, much narrower than the spread of the samples. The weights \(w_i\) are 1,
or an independent constraint on the viewing angle ([below](#viewing-angle-constraint)).

**The selection term** \(\beta(H_0)\) is the fraction of the population that would be detected. At a
given redshift a larger H₀ means a smaller distance, so a louder signal. Leaving β out biases H₀:
see [Validation](#validation-mock-bright-sirens). `--selection`:

| Value | β(H₀) | Use |
|---|---|---|
| `euclidean` | ∝ H₀³: GW-limited detection of nearby sources. It cancels the volume factor of the numerator, and with a \(d_L^2\) prior the posterior becomes \(\sum_i \mathcal{N}(v_H;\, H_0 d_i,\, \sigma_v)\), the LVK 2017 formula [\[27\]](../references.md#ref-27) | z < 0.05 |
| `injections` | the found LVK sensitivity injections of the event's observing run, carried to the source frame at 77 values of H₀ (interpolated between them) and reweighted to the population of the event's class: Power Law + Peak for black holes (as in `rates`), uniform 1–2.5 M☉ for neutron stars; the injected isotropic spins | distant events |
| `auto` (default) | `euclidean` below z = 0.05, `injections` above; `injections` with `--population fullpop4` | |

**The population** (`--population`). By default the sources are uniform in comoving volume and source-frame time.
The masses enter only with the `injections` selection, where the mass model of the class also weights each PE sample
(with the detector-frame masses carried to the source frame at each H₀), so that event and selection see the same
population; with `euclidean` they do not enter. `fullpop4` is the population of the GWTC-4.0 reanalysis of GW170817
[\[29\]](../references.md#ref-29) (Appendix E): the FullPop-4.0 mass distribution, which covers BNS, NSBH and BBH in
one (a broken power law with a dip between neutron stars and black holes, two Gaussian peaks, and a pairing function
in the mass ratio), and the Madau–Dickinson rate, fixed to the medians of the paper's spectral siren. Each sample then
also carries \(p_m\big(m_1/(1+z_i), m_2/(1+z_i)\big)\, \psi(z_i) / (1+z_i)^2\), with the detector-frame masses of
the samples (the PE mass prior is uniform in them), and the injections are reweighted to the same population. The
density is icarogw's `m1m2_paired_massratio_bplmulti_dip`, reimplemented with numpy (it agrees with icarogw to
10⁻¹³ in ln p), so that the bright siren does not need the icarogw environment.

**Sky position.** The distance must be the one at the counterpart's position: distance and sky
position are correlated through the antenna pattern. The GWTC-1 GW170817 samples were produced with
the sky fixed to AT2017gfo and are used as they are. Other samples are restricted to those within
`--sky-radius` of the counterpart.

### What the mode does

**The mode computes the posterior above itself, on a grid of H₀, with numpy: no sampler and no icarogw.** It reads
the released PE samples, as they are, and fixes the population; it does not fit the population with H₀.

1. **PE samples.** GW170817 has no PE release file: the mode reads a bundle rebuilt from the public GWTC-1 files of
   the LIGO DCC ([data sources](../data-sources.md#special-case-gw170817)), with the labels
   `C02:IMRPhenomPv2_NRTidal-LowSpin` and `-HighSpin`. Other events read their Zenodo PE file, the one the spectral
   siren reads. The labels are `--pe-label`, else those of the list, else all the labels of the file, LowSpin ones
   first. From each label: \(d_L\), the sky position, the viewing angle (`theta_jn`, else `iota`), the detector-frame
   masses and the recorded distance prior.
2. **Host redshift.** `--redshift`, else \(z = (v_r - \langle v_p \rangle)/c\) with the two velocity errors added in
   quadrature, else the redshift of the list. An event outside the list takes its run from the date in its name and
   its class from the catalog masses: BNS when both are below 3 M☉, else BBH.
3. **Sky.** When the median separation of the samples from the counterpart is below 0.1°, the sky is taken as fixed
   and all the samples are kept; otherwise those within `--sky-radius`, and the mode stops if fewer than 200 remain.
4. **Grid.** H₀ every 0.05 km/s/Mpc over `--h0-range` (3801 points by default); flat ΛCDM tabulated up to
   \(z = 10\). For each label, the kernel sum above at every grid point, with the kernel width
   \(1.06\, \mathrm{sd}(z_i)\, n^{-1/5}\) of the \(n\) samples.
5. **Selection.** `euclidean`: \(3 \ln H_0\). `injections`: the found injections of the event's run, prepared as
   for the spectral siren ([Method](hubble-constant.md#method)), weighted at 77 values of H₀ by the population
   density in the detector frame and interpolated between them; \(\beta = \sum_j w_j / N_\text{gen}\). The report
   gives the effective number of injections and warns below 50. The fixed mass models: for BBH the GWTC-3 Power
   Law + Peak of `rates` (\(\alpha = 3.4\), \(m_\text{min} = 5.1\), \(m_\text{max} = 87\) M☉, a peak at 34 M☉);
   for BNS uniform in 1–2.5 M☉; with `fullpop4`, the GWTC-4.0 medians written in `h0_bright.py`
   (Madau–Dickinson \(\gamma = 3.57\), \(\kappa = 2.97\), \(z_p = 2.59\)).
6. **Viewing angle.** With `--viewing-angle`, each sample is weighted by the Gaussian at its viewing angle, folded
   to 0–90°.
7. **Results.** One posterior per label, normalised on the grid: maximum a posteriori, 68.3% highest-density
   interval, median and 90% equal-tailed interval. The first label is the result written as `p` in
   `posterior_grid.tsv` and used by `--method joint`.

## Registered events

The bright siren uses the built-in list of known counterparts (the `COUNTERPARTS` dictionary of
`gwtc_analysis/counterpart.py`, described on the [counterpart](counterpart.md#the-list-of-known-counterparts)
page): for each event, the counterpart's position, the host redshift, the observing run, the PE file and the PE
labels used by default. In the options below, "from the list" refers to it.

| Event | Host | z | Status | Selection (auto) |
|---|---|---|---|---|
| GW170817 | NGC 4993, AT2017gfo | 3017 ± 166 km/s (\(v_r\) 3327 ± 72, \(\langle v_p \rangle\) 310 ± 150) [\[27\]](../references.md#ref-27) | confirmed (kilonova) | euclidean |
| GW190521 | AGN J124942.3+344929, ZTF19abanrhr | 0.438 ± 0.0015 [\[74\]](../references.md#ref-74) | **candidate** | injections (O3a, BBH) |

**GW190521 is a candidate.** Graham et al. proposed the ZTF flare in an AGN disk as the counterpart
[\[74\]](../references.md#ref-74). The odds of a common source are only between 1 and 12, depending on the
waveform model (Ashton et al. 2021 [\[75\]](../references.md#ref-75)). Its result is H₀ *if* the flare is
the counterpart, and the report says so. The [counterpart](counterpart.md) mode shows where the flare sits
in the GW posterior.

**Any other event** of the GWTC PE releases is run by giving its counterpart on the command line, for example the
afterglow of a coincident gamma-ray burst and the redshift of its host (found with
[search_skymaps --counterparts](search-skymaps.md#multi-messenger-coincidences)):

```bash
gwtc_analysis hubble_constant --method bright --event GW190425_081805 --ra 251.2 --dec -24.2 \
    --redshift 0.035 0.003 --host "host name" --viewing-angle 10 5      # a GRB: a jet seen nearly on-axis
```

The event takes its full name (that of its PE release, whose file is read from Zenodo); the observing run comes
from the name, and the class of the binary (BNS or BBH, the mass model of the selection) from the catalog masses,
or `--counterpart-class`. A detected GRB means a jet pointed within its opening angle of the line of sight:
`--viewing-angle` turns that into a constraint, which narrows the distance as the radio jet does for GW170817
([below](#viewing-angle-constraint)). The example is a pipeline check: GW190425 had no counterpart. A permanent
entry of the list (`COUNTERPARTS` in `gwtc_analysis/counterpart.py`) remains the way to fix a published result.

## Results

| Event, PE label | Maximum a posteriori, 68% interval | Median, 90% interval |
|---|---|---|
| GW170817, LowSpin | **69.8** (+23.6 / −8.2) | 79.9, 63.0–135.5 |
| GW170817, HighSpin | **69.4** (+14.8 / −7.4) | 74.8, 62.2–112.2 |
| GW170817, LVK 2017 [\[27\]](../references.md#ref-27) | **70.0** (+12.0 / −8.0) | |
| GW170817, LowSpin, viewing angle 20 ± 3° | **68.2** (+4.7 / −4.4) | 68.4, 61.1–76.1 |
| GW170817 and the radio jet, Hotokezaka et al. 2019 [\[73\]](../references.md#ref-73) | 70.3 (+5.3 / −5.0) | |
| GW190521 + ZTF19abanrhr, IMRPhenomXPHM | **32.9** (+115.6 / −4.6) | 94.9, 31.9–184.2 |

![H0 from GW170817 and NGC 4993](../img/modes/h0_bright_GW170817.png)

**GW170817.** The maximum a posteriori and the lower side agree with the published value. The upper
side of the LowSpin row is wider because of the **spin prior**, not of the waveform or the reanalysis: the 2017
measurement used the high-spin samples, and quoted its interval with a prior uniform in log H₀ (∝ 1/H₀), which
lowers the upper tail. The HighSpin row with that prior gives **68.8** (+12.4 / −7.0), against 70.0 (+12.0 /
−8.0). The low-spin samples have a longer low-distance tail (90% lower bound 24.9 Mpc, against 29.1 Mpc with high
spins), which comes from the **distance–inclination degeneracy**: an inclined binary is fainter than a face-on one
at the same distance (see the [counterpart](counterpart.md#the-distanceinclination-degeneracy) mode). The full ΛCDM \(d_L\) raises H₀ by 0.8% (0.6 km/s/Mpc) compared with the linear Hubble law of the
2017 paper.

### Against the GWTC-4.0 reanalysis

The GWTC-4.0 cosmology paper [\[29\]](../references.md#ref-29) reanalyses GW170817 with the same LowSpin samples
(H₀ = 77.8 (+25.7 / −11.6), median and 68%). Its setup differs from the default here in the population and the
selection term, which `--population fullpop4` reproduces:

| GW170817, LowSpin | Median, 68% | 90% |
|---|---|---|
| default (uniform in comoving volume, `euclidean`) | 79.9 (+27.6 / −12.4) | +55.6 / −16.9 |
| `--selection injections` | 80.0 | |
| `--population fullpop4 --selection euclidean` | 79.9 | |
| **`--population fullpop4`** (selection from injections) | **78.5 (+25.6 / −11.5)** | +51.9 / −15.9 |
| `--pe-distance-prior comoving` or `source-frame` | 79.8, 79.7 | |
| GWTC-4.0 paper | 77.8 (+25.7 / −11.6) | +52.2 / −16.2 |

The shift comes from the selection term with the full population, not from the mass term of the numerator: with
FullPop-4.0 most of the detectable binaries are black holes, seen out to z ≈ 0.3, and β(H₀) no longer goes as H₀³.
With it the widths are the paper's to 0.3 km/s/Mpc and the median is 0.7 above; the paper's result is an equal
mixture of icarogw and gwcosmo. The distance prior does not matter at 40 Mpc, where a uniform-in-volume prior is
\(d_L^2\) to 1%.

### Viewing-angle constraint

`--viewing-angle MEAN SIGMA` multiplies the GW likelihood by an independent Gaussian constraint on the viewing
angle (degrees, folded to 0–90°): the PE samples are weighted by it. For GW170817 the superluminal motion of the
radio jet gives \(0.25 < \theta_\text{obs}\,(d_L / 41\ \text{Mpc}) < 0.45\) rad, i.e. 14–26° at 41 Mpc
(Hotokezaka et al. 2019 [\[73\]](../references.md#ref-73)). Approximated as 20 ± 3°, it removes the inclined
orbits, hence the upper tail:

![H0 from GW170817 with the viewing angle constrained to 20 ± 3 degrees](../img/modes/h0_bright_viewing_angle_GW170817.png)

H₀ = 68.2 (+4.7 / −4.4) km/s/Mpc, against 70.3 (+5.3 / −5.0) in the paper. The Gaussian ignores the distance
dependence of the jet constraint (the angle scales as 1/d_L) and the jet modelling behind it; the constraint
must not come from the GW data.

![H0 from GW190521 and ZTF19abanrhr](../img/modes/h0_bright_GW190521.png)

**GW190521.** The posterior is bimodal because the distance posterior is, in the direction of the flare.
At z = 0.438, H₀ ≈ 30 corresponds to \(d_L\) ≈ 5.6 Gpc and H₀ ≈ 75 to 2.25 Gpc; the samples within 3°
of the flare have a median of 4.6 Gpc with 16% below 2.2 Gpc. With injections, the source-frame mass model
of the BBH class (Power Law + Peak, m_max = 87 M☉) weights both the selection and the event: GW190521's
primary, about 95 M☉ in the source frame at the flare's redshift, lies beyond it, so only a few samples
carry weight (hence the jagged curve) and the population, not the flare, shapes the upper side. Until
October 2026 the event carried no mass term (30.4, median 63.2), inconsistent with the selection. The result
depends strongly on the sky selection (median 55.5 within 5°, against 94.9 within 3°) and on the waveform
model; only IMRPhenomXPHM has enough samples near the flare in GWTC-2.1. Read it as an illustration of the
method, not as a measurement.

The published values from GW190521 and the flare alone span the same range: 43–62 km/s/Mpc depending on
the waveform (Mukherjee et al. [\[112\]](../references.md#ref-112)), 48 (+23 / −10) (Chen et al.
[\[113\]](../references.md#ref-113)), both with the 2020 PE release, and 102 (+27 / −25) in a later
reanalysis (Morton et al. [\[114\]](../references.md#ref-114)). The 90% interval here
contains all of them. The closeness of the median to the last one is partly a coincidence: here the shift
upwards comes from the mass cutoff (a higher H₀ puts the source at a higher redshift and lowers the
source-frame mass below 87 M☉), not from the distance posterior.
The selection term is computed from 106,376 found O3a injections, with
13,000–16,000 effective injections over the H₀ range.

## Validation: mock bright sirens

The analysis is tested on simulated detections with a known H₀ (`tests/test_h0_bright.py`).

- **Sources:** uniform in comoving volume up to z = 4, masses uniform in 20–40 M☉, random orientations.
- **Detection:** each source has a network SNR \(\rho = A(\mathcal{M}_\text{det})\, w(\iota) / d_L\)
  plus unit Gaussian noise, and is detected above 12. The detected events have a median z ≈ 0.5.
- **Observations:** each detected event gets the posterior of \(d_L\) given its observed SNR. It uses a
  \(d_L^2\) prior and is marginalized over the inclination, which reproduces the real
  distance–inclination degeneracy. The host redshift has an error of 0.001.
- **Injections:** they are made at H₀ = 70 with broader masses than the population (5–120 M☉), as the
  LVK ones, so that they cover the population at every trial H₀.

| Test | Result |
|---|---|
| β from the reweighted injections against the detected fraction simulated directly at each H₀ (40–120) | equal within 3% (β varies by more than a factor of 10) |
| 100 events, H₀ = 70: with the selection term | 70.9 ± 1.1 |
| the same without the selection term | 80.7 ± 2.7: biased by more than 3σ |
| P–P test, 100 events with H₀ drawn from the prior | uniform with the selection term (KS p = 0.16), not without (p = 10⁻⁴) |
| nearby sources (horizon z ≈ 0.02) | β from injections ∝ H₀^3.0, the `euclidean` selection |
| distance correlated with the viewing angle, constraint 20 ± 3° | H₀ moves to the distance of the constrained angles and narrows |

![Mock bright sirens: catalog posterior and P–P test](../img/modes/h0_bright_mock_validation.png)

The injections must cover the population at every trial H₀. In a first version of the mock, with
injection masses drawn from the population's own 20–40 M☉, β was wrong by up to a factor of 2.5 at the
edges of the H₀ range, and H₀ was biased by 3σ. The LVK injections are drawn much more broadly than any
population; the report gives the effective number of injections over the H₀ range.

## Outputs

- The HTML report ends with the **PE files of this run**: for each event the file read, its type (`nocosmo`, `cosmo`, `combined` or the GW170817 bundle) and where it was taken from (Zenodo, the S3 mirror, Galaxy, a PE cache), then the options of the run; see [Which PE file is read](../data-sources.md#which-pe-file-is-read);
- `--out-report` (HTML): the result, the caveat of a candidate counterpart, the inputs, the selection
  term and the plot;
- `--out-summary` (TSV): maximum a posteriori, 68.3% highest-density interval, median and 90%
  interval of each PE label, with its distance;
- `<workdir>/posterior_grid.tsv`: the posterior densities on an H₀ grid (`p`: the first label);
- `<workdir>/bright.json`: the event, host, redshift, labels, selection, prior range and constraint, read by
  `--method joint`;
- `<workdir>/plots/h0_posterior.png`.

## Options

Options of `--method bright` (the others of `hubble_constant` are refused):

| Option | Default | Meaning |
|---|---|---|
| `--event` | `GW170817` | an event of the [list of known counterparts](counterpart.md#the-list-of-known-counterparts) (`GW170817`, `GW190521`), or any event of the GWTC PE releases (full name) with `--ra`, `--dec` and `--redshift` |
| `--ra`, `--dec` | the position from the list | position of the counterpart, degrees (required for an event outside the list) |
| `--host` | from the list | name of the host galaxy, for the report |
| `--counterpart-class` | from the list, else from the catalog masses | `BNS` or `BBH`: the mass model of the selection (an NSBH takes the BBH one) |
| `--pe-label` | the labels from the list, else all the labels of the file | which analyses of the PE file to use ([PE labels](counterpart.md#pe-labels)): one H₀ per label; the first one is the result, the one `--method joint` uses |
| `--pe-file` | the event's bundle or Zenodo file | another PE file (PESummary layout) |
| `--cache-dir` | `.cache_gwosc` | cache of the unofficial GW170817 bundle |
| `--pe-cache` | that of the spectral siren | where Zenodo PE files are downloaded |
| `--v-recession V SIGMA`, `--v-peculiar V SIGMA` | from the list (GW170817) | velocities of a nearby host, km/s: measured recession velocity, and the peculiar velocity subtracted from it |
| `--redshift Z SIGMA` | from the list (GW190521) | Hubble-flow redshift of the host, instead of the velocities |
| `--viewing-angle MEAN SIGMA` | none | independent Gaussian constraint on the viewing angle, degrees |
| `--pe-distance-prior` | the file's, else \(d_L^2\) up to O3 and source-frame (uniform in comoving volume and source-frame time) from O4 | distance prior of the PE samples: `dl2`, `comoving` or `source-frame` (Planck15_LAL) |
| `--selection` | `auto` | `euclidean` (β ∝ H₀³), `injections` (the found injections of the event's run, reweighted to the population) or `auto` (`euclidean` below z = 0.05, `injections` above or with `fullpop4`); [Method](#method) |
| `--population` | `volume` | `volume` (uniform in comoving volume and source-frame time; with `injections`, the mass model of the class weights both the PE samples and the selection) or `fullpop4` (FullPop-4.0 and Madau–Dickinson at the GWTC-4.0 medians, [Method](#method)) |
| `--sensitivity-release`, `--sensitivity-file` | `gwtc4` | injections of the selection term |
| `--far-threshold`, `--snr-threshold` | 0.25 /yr, 10 | found injections (real searches, semi-analytic O1+O2) |
| `--sky-radius` | 3° | for samples not fixed to the counterpart |
| `--h0-range MIN MAX` | 10 200 | flat prior; `--method joint` needs the default |
| `--workdir` | `hubble_constant_bright` | the work directory |
| `--settings FILE` | none | JSON or YAML file of option values (keys: the long option names); the command line overrides it |
| `--out-report`, `--out-summary` | `hubble_constant_bright.html`, `.tsv` | the report and the summary |
| `--data-repo {zenodo,s3,galaxy}` | `zenodo` | where the PE file and the search-sensitivity injections are downloaded from: their public sources, or a mirror holding the same files, `s3` (the MinIO bucket, faster) or `galaxy` (the public usegalaxy.org history), else the public source ([S3 mirror of the public files](../data-sources.md#s3-mirror-of-the-public-files)) |

The report ends with all these options, the values used and their defaults.

The published comparison is shown only with the velocities or redshift from the list, no viewing-angle constraint,
the default population and the file's distance prior.
