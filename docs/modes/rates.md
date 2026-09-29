# rates

Merger rates of binary neutron stars (BNS), neutron star–black hole binaries (NSBH) and binary black
holes (BBH), per Gpc³ per year, as R = N / ⟨VT⟩. The method is explained in
[Merger rates](../science/merger-rates.md).

```bash
python -m gwtc_analysis.cli rates                               # GWTC-5.0 injections (O3 + O4a + O4b)
python -m gwtc_analysis.cli rates --sensitivity-release gwtc4   # GWTC-4.0 injections (O3 + O4a)
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

## Population models

| Population | Mass model | Spins | Redshift |
|---|---|---|---|
| BNS | both masses uniform in [1, 2.5] M☉ | isotropic, \|χ\| < 0.4 | constant rate |
| NSBH | black hole ∝ m^−2.35 on [2.5, 40] M☉, neutron star uniform in [1, 2.5] M☉ | BH \|χ\| < 0.99, NS \|χ\| < 0.4 | constant rate |
| BBH | GWTC-3 Power Law + Peak | isotropic, \|χ\| < 0.99 | R ∝ (1+z)^κ, reported at z = 0.2 (`--bbh-kappa`, `--bbh-z-ref`), and constant |

## Results

| Release | Candidates | BNS | NSBH | BBH at z = 0.2 |
|---|---|---|---|---|
| `gwtc5` (O3–O4b, 2.59 yr) | 258 | 26 [4, 87] | 33 [13, 67] | 25 [23, 28] |
| `gwtc4` (O3–O4a, 1.7 yr) | 154 | 43 [6, 143] | 54 [21, 109] | 26 [23, 30] |

Rates in Gpc⁻³ yr⁻¹, median [90%]. They are consistent with the LVK population papers:
[GWTC-5.0](https://arxiv.org/abs/2605.27226) (BBH 27.5–49.4 at z = 0.2 for masses 2.5–200 M☉),
[GWTC-4.0](https://arxiv.org/abs/2508.18083) (z = 0: BNS 7.6–250, NSBH 9.1–84, BBH 14–26) and
[GWTC-3](https://arxiv.org/abs/2111.03634).

## Outputs

- `--out-rates` (TSV): one row per population model, with N, ⟨VT⟩, the effective number of
  injections and the rate quantiles;
- `--out-events` (TSV): the events counted, with their class;
- `--out-report` (HTML): the tables, and the observed and selection-corrected primary-mass
  distributions.

All options: [CLI reference](../cli-reference.md#rates).
