# build_unofficial_pe: GW170817

GW170817, the binary neutron star seen with its electromagnetic counterpart
(Abbott et al. 2017 [\[38\]](../references.md#ref-38)), has no PESummary file in the catalog
releases. `build_unofficial_pe` builds a PESummary-compatible `PEDataRelease` bundle for it, from
public GWTC-1 products only, so that `parameters_estimation` can treat it like any other event.

```bash
python -m gwtc_analysis.cli build_unofficial_pe --src-name GW170817          # build, or reuse the cache
python -m gwtc_analysis.cli build_unofficial_pe --src-name GW170817 --force  # rebuild
```

## Sources

The files are downloaded into `~/.gwcache` the first time they are needed:

| Product | DCC release | File |
|---|---|---|
| Posterior samples (IMRPhenomPv2_NRTidal, low and high spin) | [LIGO-P1800370](https://dcc.ligo.org/LIGO-P1800370/public) | `GW170817_GWTC-1.hdf5` |
| PSDs (H1, L1, V1) | [LIGO-P1900011](https://dcc.ligo.org/LIGO-P1900011/public) | `GWTC1_GW170817_PSDs.dat` |
| Calibration uncertainty envelopes | [LIGO-P1900040](https://dcc.ligo.org/LIGO-P1900040/public) | `GWTC1_GW170817_CalEnv/GWTC1_GW170817_{H,L,V}_CalEnv.txt` |
| Skymap | [LIGO-P1800381](https://dcc.ligo.org/LIGO-P1800381/public) | `GW170817_skymap.fits.gz` |

## Reconstructing the missing parameters

The public samples carry no polarization angle, phase, coalescence time or likelihood. The builder:

1. draws the polarization `psi` and the phase from their priors, sets `geocent_time` to the trigger
   time, and derives a synthetic ranking `log_likelihood`;
2. for the sample ranked first (the maximum-likelihood sample used by the strain overlay), fits
   `geocent_time`, `phase` and `psi` to the GWOSC strain by maximizing the coherent network
   likelihood with the public PSDs, keeping the intrinsic parameters and the sky position. The build
   log reports the recovered matched-filter SNRs: about 18 in H1, 25 in L1 and 31 for the network;
3. derives the detector arrival times (`H1_time`, `L1_time`, `V1_time`) from `geocent_time`, `ra`
   and `dec` with the LAL detector delays.

All other samples keep their prior-drawn extrinsic values. If the strain fit fails (it needs network
access to GWOSC), the bundle is still built and the log warns that the overlay will not be coherent.

## Cache

The bundle is cached with a fingerprint of its recipe (`<bundle>.recipe.json`) and rebuilt when the
recipe or a source file changes. A missing source that cannot be downloaded stops the build, and the
overlays of this special case are then not produced.

All options: [CLI reference](../cli-reference.md#build_unofficial_pe).
