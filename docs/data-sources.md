# Data sources and versions

`gwtc_analysis` uses only public data products of the LVK collaboration.

## Event lists: GWOSC

Event metadata (names, GPS times, false-alarm rates, median masses and distances, network SNR) come
from the event API of the [Gravitational Wave Open Science Center](https://gwosc.org)
(`https://gwosc.org/eventapi/jsonfull/<list>/`). The catalog keys map to these lists:

| Key | GWOSC list | Zenodo PE and skymap release |
|---|---|---|
| `GWTC-1` | `GWTC-1-confident` | (no separate skymaps; GWTC-2.1 covers O1–O2) |
| `GWTC-2.1` | `GWTC-2.1-confident` | [Zenodo 6513631](https://zenodo.org/records/6513631) |
| `GWTC-3` | `GWTC-3-confident` | [Zenodo 22685054](https://zenodo.org/records/22685054) |
| `GWTC-4` | `GWTC-4.0` | [Zenodo 17602505](https://zenodo.org/records/17602505) |
| `GWTC-5` | `GWTC-5.0` | [Zenodo 20348005](https://zenodo.org/records/20348005) (part 1, with the skymaps) and [20348006](https://zenodo.org/records/20348006) (part 2) |

The `rates` and `hubble_constant` modes also read the **marginal** lists (`GWTC-2.1-marginal`,
`GWTC-3-marginal`): some events used by the LVK population and cosmology analyses are only there
(see [GW200105](known-data-issues.md#gw200105-only-in-the-marginal-list)).

## Parameter estimation and skymaps: three repositories

`--data-repo` chooses where PE files and skymaps are read from:

`zenodo` (default)
:   The official catalog releases on Zenodo (table above).

`s3`
:   The `gwtc` bucket of an S3/MinIO mirror at `https://minio-dev.odahub.fr`.

`galaxy`
:   Collections staged by Galaxy (`./galaxy_inputs/<CATALOG>-PE`, `-SKYMAPS`); when none is staged,
    PE files are downloaded over HTTP from the public usegalaxy.org "GWTC" published history
    (anonymous, no API key).

## Zenodo release versions

The Zenodo releases are versioned: GWTC-3, for instance, has v1, v2 and v3. With
`--data-repo zenodo`, each catalog uses its **latest** version by default. The version listing is
fetched from the Zenodo API (`/api/records/<id>/versions`) and cached for one day in
`~/.cache_gwtc_analysis/zenodo`, so a new release is picked up automatically.

An older version is selected with `--zenodo-version CATALOG=VERSION`, versions being numbered from the
oldest (`v1`); `latest` is also accepted:

```bash
python -m gwtc_analysis.cli zenodo_releases --catalogs GWTC-3 GWTC-4     # list the versions
python -m gwtc_analysis.cli search_skymaps --catalogs GWTC-3 --ra-deg 40 --dec-deg -30 --zenodo-version GWTC-3=v2
python -m gwtc_analysis.cli parameters_estimation --src-name GW200105_162426 --zenodo-version GWTC-3=v2
```

Skymap tarballs are cached per Zenodo record, and the PE index is rebuilt when the selected records
change. If zenodo.org is unreachable, the cached listing is used; with no cache at all, the latest
version falls back to the record of the table above.

## Search-sensitivity injections

The `rates` and `hubble_constant` modes use the LVK **search-sensitivity estimates**: large sets of
simulated signals injected into the detector data and searched by the LVK pipelines (see
[Selection effects and injections](science/selection-effects.md)). They are downloaded automatically
from the latest version of their Zenodo record and cached:

| Release key | Record | File used by `rates` | File used by `hubble_constant` |
|---|---|---|---|
| `gwtc4` | [Zenodo 16740128](https://zenodo.org/records/16740128) | `mixture-real_o3_o4a-cartesian_spins_*.hdf` (~400 MB) | `mixture-semi_o1_o2-real_o3_o4a-cartesian_spins_*.hdf` |
| `gwtc5` | [Zenodo 19500052](https://zenodo.org/records/19500052) | `mixture-real_o3_o4a_o4b-cartesian_spins_*-clipped.hdf` (~900 MB) | `mixture-semi_o1_o2-real_o3_o4a_o4b-cartesian_spins_*-clipped.hdf` |

The *real* mixtures cover O3 onward, with signals injected in the real data. The *semi* mixtures add
semi-analytic O1+O2 injections (detection from the SNR computed on the noise spectra of the time),
so that the whole catalog since 2015 can be used. The release is described in
Essick et al. 2025 [\[34\]](references.md#ref-34).

## Special case: GW170817

GW170817 has no PESummary file in the catalog releases. `build_unofficial_pe` rebuilds one from the
public GWTC-1 products of the LIGO DCC: see [build_unofficial_pe](modes/unofficial-pe.md).
