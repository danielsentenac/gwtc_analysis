# rates

Merger rates of binary neutron stars (BNS), neutron star–black hole binaries (NSBH) and binary black
holes (BBH), per Gpc³ per year, as R = N / ⟨VT⟩. The method is explained in
[Merger rates](../science/merger-rates.md).

```bash
gwtc_analysis rates                               # GWTC-5.0 injections (O3 + O4a + O4b)
gwtc_analysis rates --sensitivity-release gwtc4   # GWTC-4.0 injections (O3 + O4a)
```

## What is counted

- **N:** the GWOSC candidates of GWTC-2.1, GWTC-3, GWTC-4.0 and GWTC-5.0, confident and marginal
  lists, inside the observing periods covered by the injections and with a false-alarm rate below
  `--far-threshold` (1 per year by default). They are classified by their median source-frame masses,
  neutron stars being below `--ns-max-mass` (2.5 M☉). Candidates without masses are left out.
- **⟨VT⟩:** the sensitive volume-time of each population, from the LVK injections, downloaded
  automatically (`--sensitivity-release gwtc5` by default, or `gwtc4`; `--sensitivity-file` for a
  local file).

The observing periods are those of the injections (their times split at gaps longer than a week), so
candidates of the engineering run ER15, just before O4a, are not counted.

## Selecting catalogs

`--catalogs` restricts the rates to the observing runs of some catalogs, for the events **and** the
injections, so that the counts and the sensitive volume-time describe the same observing time:

| Catalog | Runs |
|---|---|
| `GWTC-1` | O1, O2 |
| `GWTC-2.1` | O3a |
| `GWTC-3` | O3b |
| `GWTC-4` | O4a |
| `GWTC-5` | O4b |

```bash
gwtc_analysis rates --catalogs GWTC-5          # O4b only
gwtc_analysis rates --catalogs GWTC-4 GWTC-5   # O4a + O4b
gwtc_analysis rates --catalogs ALL             # O1 to O4b
```

Without `--catalogs`, the rates cover the runs of the real-injection mixture of the release (O3 onward).
Selecting GWTC-1 uses the release's mixture with semi-analytic O1+O2 injections, found above
`--snr-threshold` (10). The restriction keeps the importance sums of the selected runs: with the
mixture weights, they are the ⟨VT⟩ of those runs, and the ⟨VT⟩ of separate catalogs add up to that of
all of them. An event listed in two catalogs (the O1–O2 events of GWTC-1 and GWTC-2.1) is counted
once.

BBH rates of each catalog (GWTC-5.0 injections, R ∝ (1+z)^2.9, at z = 0.2, Gpc⁻³ yr⁻¹, median [90%]):

| Catalog | BBH candidates | BBH rate |
|---|---|---|
| GWTC-2.1 (O3a) | 37 | 31.1 [23.5, 40.3] |
| GWTC-3 (O3b) | 23 | 20.9 [14.5, 28.8] |
| GWTC-4 (O4a) | 84 | 25.2 [21.0, 30.0] |
| GWTC-5 (O4b) | 104 | 24.8 [21.0, 29.0] |
| O3a–O4b | 248 | 25.2 [22.7, 27.9] |

The catalogs agree with one another. O4b has no BNS or NSBH candidate below 1 per year: its BNS and
NSBH rates are upper limits (90%: 113 and 40 Gpc⁻³ yr⁻¹).

## Population models

| Population | Mass model | Spins | Redshift |
|---|---|---|---|
| BNS | both masses uniform in [1, 2.5] M☉ | isotropic, \|χ\| < 0.4 | constant rate |
| NSBH | black hole ∝ m^−2.35 on [2.5, 40] M☉, neutron star uniform in [1, 2.5] M☉ | BH \|χ\| < 0.99, NS \|χ\| < 0.4 | constant rate |
| BBH | GWTC-3 Power Law + Peak [\[12\]](../references.md#ref-12) | isotropic, \|χ\| < 0.99 | R ∝ (1+z)^κ, reported at z = 0.2 (`--bbh-kappa`, `--bbh-z-ref`), and constant |

## Results

| Release | Candidates | BNS | NSBH | BBH at z = 0.2 | BBH at z = 0 |
|---|---|---|---|---|---|
| `gwtc5` (O3–O4b, 2.59 yr) | 259 | 26 [4, 87] | 33 [13, 67] | 25 [23, 28] | 15 [13, 17] |
| `gwtc4` (O3–O4a, 1.7 yr) | 155 | 43 [6, 143] | 54 [21, 109] | 26 [23, 30] | 15 [13, 18] |

Rates in Gpc⁻³ yr⁻¹, median [90%]. Compare with the LVK population papers:
GWTC-5.0 [\[14\]](../references.md#ref-14) (BBH 27.5–49.4 at z = 0.2 for masses 2.5–200 M☉),
GWTC-4.0 [\[13\]](../references.md#ref-13) (z = 0: BNS 7.6–250, NSBH 9.1–84, BBH 14–26) and
GWTC-3 [\[12\]](../references.md#ref-12). The BNS and NSBH rates fall inside the LVK intervals. The
BBH rate with the fixed GWTC-3 Power Law + Peak sits at the low edge of the GWTC-5.0 interval, with a
much narrower interval, and so does its value at z = 0 against the GWTC-4.0 interval: see the next
section. The BBH rate at z = 0 is the same evolving rate, R ∝ (1+z)^κ, taken at z = 0; it is given for
comparison with local rates. The constant-rate row is a different assumption: without evolution, the
rate is in effect averaged over the redshifts of the detected binaries (z ≈ 0.3–0.5), hence its
higher value (41 [37, 46] with `gwtc5`).

## Population uncertainty and the Multi Peak model

With fixed population shapes the intervals are Poisson only, and the BBH rate depends on the mass
model through ⟨VT⟩: light binaries are hard to detect, so a model with more of them (the peak near
10 M☉ found since GWTC-4.0) has a smaller ⟨VT⟩ and gives a higher rate. `--population-posterior`
takes a posterior of the BBH mass model and of the rate evolution (γ, κ, z_p), for instance a
`hubble_constant` work directory, and computes ⟨VT⟩ for `--n-draws` posterior samples (200 by
default), mixing the Poisson posterior of each: the interval then includes the population uncertainty.
`--mass-model` chooses Power Law + Peak (`plp`, default) or Multi Peak (`mltp`: the power law with two
Gaussian peaks near 10 and 35 M☉, the model of the LVK cosmology papers); `mltp` needs a posterior.

```bash
gwtc_analysis rates --mass-model mltp --population-posterior h0_mltp_workdir            # a hubble_constant work directory
gwtc_analysis rates --mass-model mltp --population-posterior posterior_reweighted.tsv   # or its posterior file
```

### The population posterior

`--population-posterior` takes posterior samples of the BBH population: the mass-model parameters and
the rate-evolution parameters. A `hubble_constant` run of the same `--mass-model` produces them.

**What to give.** Either the `hubble_constant` work directory, in which case `posterior_reweighted.tsv` is
read (or `posterior.tsv` if the run was not reweighted), or a TSV file directly.

**Columns the file must contain**, one row per sample:

| `--mass-model` | Mass-model columns | Rate-evolution columns |
|---|---|---|
| `plp` | `alpha beta mmin mmax delta_m mu_g sigma_g lambda_peak` | `gamma kappa zp` |
| `mltp` | `alpha beta mmin mmax delta_m mu_g_low sigma_g_low mu_g_high sigma_g_high lambda_g lambda_g_low` | `gamma kappa zp` |

These are the names `hubble_constant` writes. Other columns, such as `H0`, are ignored. A missing
column stops the run with its name.

| BBH mass model | ⟨VT⟩ (Gpc³ yr) | BBH at z = 0.2 | BBH at z = 0 |
|---|---|---|---|
| Power Law + Peak, GWTC-3 values, R ∝ (1+z)^2.9 | 16.7 | 25 [23, 28] | 15 [13, 17] |
| Multi Peak, 200 posterior draws, Madau–Dickinson evolution | 14.4 | 31 [21, 51] | 17 [9, 28] |
| LVK, GWTC-5.0 / GWTC-4.0 | | 27.5–49.4 | 14–26 |

*GWTC-5.0 injections, O3 to O4b, 248 BBH candidates. The Multi Peak posterior is that of a
spectral-siren run on the GWTC-4.0 events, with H₀ free; a posterior of the same events, with the
cosmology fixed, would be the self-consistent choice. With the population uncertainty the BBH rate
agrees with the LVK intervals.*

## Options

| Option | Default | Meaning |
|---|---|---|
| `--sensitivity-release` | `gwtc5` | LVK injections downloaded from Zenodo: `gwtc5` = O3 + O4a + O4b (~900 MB), `gwtc4` = O3 + O4a (~400 MB) |
| `--sensitivity-file` | none | a local LVK search-sensitivity injection file (same format as the release files), read instead of downloading the `--sensitivity-release` file from Zenodo, which is then ignored; with `--catalogs GWTC-1` it must be the mixture that includes the semi-analytic O1+O2 injections |
| `--catalogs` | the runs of the real-injection mixture (O3 onward) | catalog keys (GWTC-1 … GWTC-5, GWTC-4.1, or ALL): events and injections restricted to their observing runs (see [Selecting catalogs](#selecting-catalogs)) |
| `--far-threshold` | 1 per year | events and injections below this false-alarm rate |
| `--snr-threshold` | 10 | semi-analytic O1+O2 injections above this network SNR (with GWTC-1) |
| `--ns-max-mass` | 2.5 M☉ | neutron stars below it, black holes above, for the classification and the BNS and NSBH models |
| `--mass-model` | `plp` | BBH mass model: `plp` (Power Law + Peak, GWTC-3 values unless `--population-posterior`) or `mltp` (Multi Peak, two peaks near 10 and 35 M☉; needs `--population-posterior`) |
| `--population-posterior` | none | samples of the BBH population (a `hubble_constant` run, or a TSV file): the BBH rate is then computed over them, see [The population posterior](#the-population-posterior) |
| `--n-draws` | 200 | posterior draws with `--population-posterior` |
| `--bbh-kappa` | 2.9 | BBH rate evolution R ∝ (1+z)^κ, for the fixed Power Law + Peak (not used with `--population-posterior`, whose evolution is fitted) |
| `--bbh-z-ref` | 0.2 | redshift at which the evolving BBH rate is reported (also at z = 0) |
| `--out-rates`, `--out-events`, `--out-report` | `merger_rates.tsv`, `merger_rates_events.tsv`, `merger_rates.html` | output files (see below) |
| `--plots-dir` | `rates_plots` | directory of the plots |

## Outputs

- `--out-rates` (TSV): one row per population model, with N, ⟨VT⟩, the effective number of
  injections and the rate quantiles; the evolving BBH rate appears twice, at `--bbh-z-ref` and at
  z = 0 (the `stochastic` mode uses the first, at `--bbh-z-ref`);
- `--out-events` (TSV): the events counted, with their class;
- `--out-report` (HTML): the tables, and the observed and selection-corrected primary-mass
  distributions.

## Example

```bash
gwtc_analysis rates
```

The HTML report shows, with the rates table, the observed and the selection-corrected primary-mass
distributions:

![Observed and selection-corrected primary-mass distributions](../img/modes/rates_mass_distribution.png){ width="640" }

*O3 to O4b, GWTC-5.0 injections. Top: detected candidates per bin of primary mass. Bottom: merger
rate per logarithmic mass interval, the counts divided by the sensitive volume-time of each bin. The
peaks near 10 and 35 M☉ stand out once the heavy binaries' larger detection volume is removed.*

All options: [CLI reference](../cli-reference.md#rates).
