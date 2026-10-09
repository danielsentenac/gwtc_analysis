# Changelog

## 0.8.0 (unreleased)

### Hubble constant: one mode, four methods

`hubble_constant` now holds every way of measuring H₀, chosen with `--method`:

- `spectral` (default): the BBH mass spectrum, as before.
- `dark`: the mass spectrum and a galaxy catalog (`--galaxy-catalog`, built by the new `galaxy_catalog` mode).
  Later stages take the method from the work directory, without `--galaxy-catalog` again.
- `bright`: an event with an identified host (GW170817 and NGC 4993, the candidate GW190521 flare), in seconds.
  New `--viewing-angle MEAN SIGMA`: an independent constraint on the viewing angle (e.g. the radio jet of
  GW170817), applied as weights on the PE samples. The result is a work directory (`posterior_grid.tsv`,
  `bright.json`).
- `joint`: the product of independent results (`--inputs DIR DIR ...`: work directories of the other methods, or
  posterior TSVs), with the table of each input and of the joint posterior. The inputs are checked to be
  independent: at most one spectral or dark siren, no event in two inputs.

The options of another method are refused. The default work directory and outputs are named after the method:
`hubble_constant_<method>`, `hubble_constant_<method>.html`, `hubble_constant_<method>.tsv` (they were
`hubble_constant_run`, `hubble_constant.html`, `hubble_constant.tsv`).

### New mode: `counterpart`

Where an electromagnetic counterpart sits in the GW posterior of an event: the searched probability of its position
(from the LVK sky map, else the PE samples), the distance along its line of sight against the distance of the host
redshift for Planck and SH0ES H₀, the viewing angle and the distance–inclination degeneracy, and an optional
viewing-angle constraint. `--ra`/`--dec` test another position.

### Removed: `bright_siren`

Replaced by `hubble_constant --method bright` (H₀), `hubble_constant --method joint` (combination with the spectral
or dark siren, formerly `--spectral-posterior`) and `counterpart` (degeneracy plot and table). `--src-name` is now
`--event`; `--plots-dir` is replaced by `<workdir>/plots`. The module `bright_siren.py` is split into `h0_bright.py`,
`h0_joint.py` and `counterpart.py`.

### Dark sirens

- `galaxy_catalog` mode: an icarogw line-of-sight galaxy catalog (GLADE+ K band by default, or any catalog read in
  chunks), local or as Slurm jobs; GLADE+ selection options (`--glade-types`, `--glade-redshift`, `--glade-sigmaz`,
  `--sigmaz`, `--where`, `--zmin`, `--ptype`).
- `hubble_constant` dark siren with icarogw's `CBC_catalog_vanilla_rate`; Slurm executor (`--executor slurm`,
  chains of batch jobs); likelihood thresholds `--neff-pe`, `--neff-inj`; settings files (`--settings`); the
  resolved options recorded in `<workdir>/options_<mode>.json`.
- Validation: the GWTC-4.0 Power Law + Peak dark siren with GLADE+, 114.6 (+41.5 / −33.8) km/s/Mpc with 5000 PE
  samples per event, against 115.4 (+40.1 / −33.8) published.
- `--catalogs` outside the observing runs of the injection release is refused (it was silently restricted).

### Other changes

- `rates`: Multi Peak mass model, BBH rate over the draws of a `hubble_constant` posterior
  (`--population-posterior`), evolving BBH rate also at z = 0.
- `stochastic`: Multi Peak mass model, read from the `hubble_constant` posterior.
- `parameters_estimation`: 3D sky map from the release FITS archive, with GLADE+ host-galaxy candidates.
- `search_skymaps`, `skymap3d`: one sky map per event, the requested waveform label or a shared fallback.
- New modes `neutron_star_eos` (joint equation of state from GW170817 and GW190425) and `spin_population` (BBH
  effective-spin population).
- Documentation: Waveform models page, the stochastic background of compact binaries, dark sirens, the bright siren
  and the counterpart pages.
