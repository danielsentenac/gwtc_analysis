# Data sources and versions

`gwtc_analysis` uses only public data products of the LVK collaboration.

## Event lists: GWOSC

Event metadata (names, GPS times, false-alarm rates, median masses and distances, network SNR) come
from the event API of the [Gravitational Wave Open Science Center](https://gwosc.org)
(`https://gwosc.org/eventapi/jsonfull/<list>/`). The catalog keys map to these lists:

| Key | GWOSC list | Observing runs of its events | Events | Zenodo PE and skymap release |
|---|---|---|---|---|
| `GWTC-1` | `GWTC-1-confident` | O1 (3), O2 (8) | 11 | none of its own: the GWTC-2.1 release, [Zenodo 6513631](https://zenodo.org/records/6513631), which re-analysed O1–O2 |
| `GWTC-2.1` | `GWTC-2.1-confident` | O3a (44), plus O1 (3) and O2 (7) re-analysed | 54 | [Zenodo 6513631](https://zenodo.org/records/6513631) |
| `GWTC-3` | `GWTC-3-confident` | O3b | 35 | [Zenodo 22685054](https://zenodo.org/records/22685054) |
| `GWTC-4` | `GWTC-4.0` | O4a (128), plus GW230518 from the engineering run ER15 | 129 | [Zenodo 17602505](https://zenodo.org/records/17602505) |
| `GWTC-4.1` | `GWTC-4.1` | O4a (138), plus GW230517 and GW230518 from ER15: GWTC-4.0 re-analysed, with 11 new events | 140 | [Zenodo 20275769](https://zenodo.org/records/20275769) |
| `GWTC-5` | `GWTC-5.0` | O4b (156), plus 5 events of 6–8 April 2024, just before the start of O4b | 161 | [Zenodo 20348005](https://zenodo.org/records/20348005) (part 1, with the skymaps) and [20348006](https://zenodo.org/records/20348006) (part 2) |

`GWTC-4.1` is an **update** of GWTC-4.0: the same O4a data re-analysed, with the 129 events of
GWTC-4.0 and 11 new ones. It is used only when named (`--catalogs GWTC-4.1`), in place of `GWTC-4` for
the O4a events: `ALL` and the defaults of every mode keep GWTC-4.0, the catalog of the published
analyses. Its PE files are read only on request (`--zenodo-version GWTC-4.1=latest`).

All the keys are **confident** catalogs: every event of these lists has a probability of
astrophysical origin p_astro ≥ 0.5 (the re-analysed O1–O2 events of GWTC-2.1 carry no p_astro value).
Only the `rates`, `hubble_constant` and `spin_population` modes select events from the marginal lists (`counterpart` and `search_skymaps --counterparts` also look events up in them) (below).

Event counts of the GWOSC lists in September 2026. A catalog key holds the events of *its* release
only: GWTC-4.0, for instance, does not repeat the O1–O3 events. The exception is GWTC-2.1, whose list
also has 10 O1–O2 events re-analysed with the GWTC-2.1 methods; the same events are in GWTC-1. Events outside the observing runs
(engineering runs) are left out by the `rates` and `hubble_constant` modes, which select events by
the periods covered by the injections or by the run dates.

## Observing runs

| Run | Start (UTC) | End (UTC) |
|---|---|---|
| O1 | 2015-09-12 | 2016-01-19 |
| O2 | 2016-11-30 | 2017-08-25 |
| O3a | 2019-04-01 | 2019-10-01 |
| O3b | 2019-11-01 | 2020-03-27 |
| O4a | 2023-05-24 | 2024-01-16 |
| O4b | 2024-04-10 | 2025-01-28 |

These are the GWOSC run boundaries used by the `hubble_constant` event selection.

The `rates`, `hubble_constant` and `spin_population` modes also read the **marginal** lists (`GWTC-2.1-marginal`,
`GWTC-3-marginal`): some events used by the LVK population and cosmology analyses are only there
(see [GW200105](known-data-issues.md#gw200105-only-in-the-marginal-list)).

## Parameter estimation and skymaps: three repositories

`--data-repo` chooses where PE files and skymaps are downloaded from. It does not choose **which** PE
file is read: that depends only on the mode and on the release version (below,
[Which PE file is read](#which-pe-file-is-read)).

`zenodo` (default)
:   The official catalog releases on Zenodo (table above).

`s3`
:   The `gwtc` bucket of an S3/MinIO mirror at `https://minio-dev.odahub.fr`, much faster than Zenodo (the
    other public files are mirrored too: [S3 mirror of the public files](#s3-mirror-of-the-public-files)). Zenodo
    remains the reference and the fallback: a PE file is taken from the mirror when it holds the file of that
    name, and the skymaps of `catalog_statistics` and `search_skymaps` when it holds those of the selected
    release version (their file names carry it, e.g. `IGWN-GWTC4p0-38214bd95_724-…`); a file or a catalog the
    mirror does not have, a mirror that cannot be reached and a failed transfer all fall back to Zenodo, with
    a log line. Older versions kept on the mirror are thus never read by mistake, and `--zenodo-version` applies
    with `s3` as with `zenodo`.

`galaxy`
:   Collections staged by Galaxy (`./galaxy_inputs/<CATALOG>-PE`, `-SKYMAPS`): the staged collection is the
    input of the Galaxy tool, and the rule of the mode picks the file among its files of the event. When none
    is staged, the file chosen in the Zenodo listing is downloaded over HTTP from the public usegalaxy.org
    "GWTC" published history (anonymous, no API key) when the history has it, else from Zenodo.

## Zenodo release versions

The Zenodo releases are versioned: GWTC-3, for instance, has v1, v2 and v3. With
`--data-repo zenodo`, each catalog uses its **latest** version by default. The version listing is
fetched from the Zenodo API (`/api/records/<id>/versions`) and cached for one day in
`~/.cache_gwtc_analysis/zenodo`, so a new release is picked up automatically.

An older version is selected with `--zenodo-version CATALOG=VERSION`, versions being numbered from the
oldest (`v1`); `latest` is also accepted. It is taken by `catalog_statistics` and `search_skymaps` (with
`--data-repo zenodo` or `s3`) and by `parameters_estimation` (with any `--data-repo`), and `zenodo_releases` lists the
versions. The other modes that read PE files take no `--zenodo-version`: they read the latest versions,
whatever `--data-repo`. `zenodo_releases` has one option, `--catalogs` (default `ALL`):

```bash
gwtc_analysis zenodo_releases --catalogs GWTC-3 GWTC-4     # list the versions
gwtc_analysis search_skymaps --catalogs GWTC-3 --ra-deg 40 --dec-deg -30 --zenodo-version GWTC-3=v2
gwtc_analysis parameters_estimation --src-name GW200105_162426 --zenodo-version GWTC-3=v2
```

Skymap tarballs are cached per Zenodo record, and the PE index is rebuilt when the selected records
change. If zenodo.org is unreachable, the cached listing is used; with no cache at all, the latest
version falls back to the record of the table above.

## Which PE file is read

The PE release of an event holds one file per event from GWTC-4.0 on (`*-combined_PEDataRelease.hdf5`),
but two for GWTC-2.1 and GWTC-3, with the same analyses and different distance priors:

| File | Samples | Distance prior recorded in the file |
|---|---|---|
| `*_PEDataRelease_mixed_nocosmo.h5` | as sampled, with a prior ∝ D_L² (uniform in Euclidean volume) | `PowerLaw(alpha=2)`: the prior of the samples |
| `*_PEDataRelease_mixed_cosmo.h5` | reweighted to a source population uniform in comoving volume and source-frame time (Planck15): the posteriors quoted in the GWTC-2.1 and GWTC-3 papers; fewer samples after the reweighting | still `PowerLaw(alpha=2)`: **not** the prior of the samples |

A hierarchical analysis divides each event's samples by the prior they were drawn with, so it needs the
`nocosmo` file, whose recorded prior is the right one. The file read is therefore chosen in the Zenodo
listing of the release (its latest version, or that of `--zenodo-version`), by a rule that depends only on
the mode:

| Modes | Rule | GWTC-2.1, GWTC-3 | GWTC-4.0 on |
|---|---|---|---|
| `hubble_constant` (all methods), `counterpart`, `neutron_star_eos`, `scalar_cloud`, `spin_population`, `lens_pairs` | samples under their recorded prior | `*_mixed_nocosmo.h5` | the combined file |
| `parameters_estimation` | the catalog's posteriors | `*_mixed_cosmo.h5` | the combined file |

`--data-repo` then only says where that file is downloaded from (above), and a result does not depend on
what the S3 mirror or the Galaxy history hold: whenever the mirror lacks a file, or cannot be reached, or a
transfer from it fails, the file is downloaded from Zenodo. Every HTML report of these modes ends with the **PE files of
this run**: for each event the file read, its type and where this run took it from. GW170817 has no release file: the modes read the bundle
rebuilt by [`build_unofficial_pe`](modes/unofficial-pe.md).

Two safeguards keep it so:

- **Caches.** A PE file is reused from the PE cache (`files/`) only under its exact name. The samples
  extracted from it (`samples/`, `samples_spin/`, `samples_lens/`) record the name of their source file, and
  an extract made from another file (an earlier release version, or the other file of the pair) is made
  again. Extracts made before this rule, which may come from a `cosmo` file of the mirror, are therefore
  redone on the next run.
- **A `cosmo` file given explicitly** (for instance `--pe-file` of the bright siren) is read with the prior
  its samples follow, uniform in comoving volume and source-frame time (Planck15), not with the recorded
  one. Reweighting the `nocosmo` samples of GW151226 by this prior over D_L² gives the distance quantiles of
  its `cosmo` file (5, 50 and 95%: 263, 461 and 621 Mpc against 264, 460 and 622).

## S3 mirror of the public files

The other public files the modes download are mirrored in the same bucket, under a key derived from their URL:

| URL | Key in the `gwtc` bucket |
|---|---|
| `https://zenodo.org/records/<record>/files/<name>` | `public/zenodo/<record>/<name>` |
| `https://dcc.ligo.org/public/<nnnn>/<id>/<version>/<name>` | `public/dcc/<id>/<version>/<name>` |
| `https://arxiv.org/e-print/<id>` | `public/arxiv/<id>.tar.gz` |

The file is chosen as without the mirror (the same URL; for the injections, the latest version of the Zenodo
record, so a new version is a new key), and `--data-repo s3` downloads it from the bucket when the bucket holds
it, else from its public source. The modes that read these files all take `--data-repo {zenodo,s3}`:

| Files | Modes | Source |
|---|---|---|
| Search-sensitivity injections, `real` mixtures (`gwtc4` ~400 MB, `gwtc5` ~890 MB) | `rates` | Zenodo [16740128](https://zenodo.org/records/16740128), [19500052](https://zenodo.org/records/19500052) |
| Search-sensitivity injections, `semi` mixtures (~620 MB, ~1.1 GB) | `hubble_constant` (all methods), `spin_population` | the same records |
| GW250114 data release, and `posterior_samples_NRSur7dq4.h5` (`--with-imr`) | `area_law` | Zenodo [16877102](https://zenodo.org/records/16877102) |
| LVK continuous-wave upper limits (known pulsars, narrowband, all-sky O4a, O3 boson clouds), and the arXiv source of 2501.01495 | `continuous_waves` | LIGO DCC G2402400, P2500416, P2100343; arXiv |
| PSD and skymap supplements of GW230529, GW190425, GW200105 | `parameters_estimation` | Zenodo [10845779](https://zenodo.org/records/10845779); DCC P2000026, P2100143 |
| GWTC-1 samples, PSDs, skymap and calibration envelopes of GW170817 | `build_unofficial_pe` (and `parameters_estimation`, `counterpart`, `hubble_constant --method bright`, `neutron_star_eos` through it) | LIGO DCC P1800370, P1900011, P1800381, P1900040 |

The reports list these files under **Public data files of this run**, with where each one was taken from.
Queries are not mirrored: the GWOSC event lists and strain data, the GLADE+ galaxies (VizieR), and the GRB and
neutrino catalogs (HEASARC, GCN).

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
Essick et al. 2025 [\[38\]](references.md#ref-38).

## Special case: GW170817

GW170817 has no PESummary file in the catalog releases. `build_unofficial_pe` rebuilds one from the
public GWTC-1 products of the LIGO DCC: see [build_unofficial_pe](modes/unofficial-pe.md).

## Adding a new catalog

All the catalogs, observing runs and sensitivity releases are described in one module,
`gwtc_analysis/catalog_registry.py`; every other table of the package (catalog keys, GWOSC list names,
Zenodo records, S3 prefixes, run mappings, the `ALL` expansions, the releases of `rates` and
`hubble_constant`, the help texts) is derived from it. When the LVK publishes a new catalog (for
instance the one of the O4c run):

0. **`check_catalogs`** compares what GWOSC and Zenodo publish with the registry:

    ```bash
    gwtc_analysis check_catalogs --out-json check.json
    ```

    It lists the registry catalogs with their GWOSC event counts and Zenodo records (latest version,
    concept ID, number of versions, publication date), then reports the GWTC event lists and observing
    runs that the registry does not describe, and the
    registry records that have a newer Zenodo version. For each new list it drafts the registry entry:
    the observing runs of its events, its Zenodo records (found from the PE links of its events in the
    GWOSC v2 API) with their concept IDs and skymap tarball, and `update_of` when the list covers the runs
    of an existing catalog (a re-analysis, as GWTC-4.1 of GWTC-4.0). The draft is a starting point: check
    it against the release notes. On the current registry it reports nothing new.

    | Option | Default | Description |
    |---|---|---|
    | `--out-json` | none | JSON file with the full report (registry catalogs, new lists and runs, newer Zenodo versions, draft entries) |
    | `--sample-events` | 3 | events of each new list whose PE links are followed to find its Zenodo records |

1. **`OBSERVING_RUNS`**: add the run with its official start and end of observing (the GWOSC run
   limits may include the engineering run before it).
2. **`CATALOGS`**: add an entry with the catalog key, its GWOSC event list, its run(s), its Zenodo
   record(s) (one per part of a split release, the first one with the skymap tarball name) and its S3
   prefix. The record IDs are only starting points: the latest version of each is resolved at run time.
   An update of an existing catalog (a re-analysis of the same runs) takes `update_of`: it is then
   used only when named, and `ALL` and the defaults are unchanged.
3. **`SENSITIVITY_RELEASES`**, if the catalog comes with injections: its Zenodo record, the covered runs
   and the patterns of its real and semi-analytic mixture files; set `DEFAULT_RATES_RELEASE` to it. The
   spectral-siren values of the matching cosmology paper, when published, go in `published`.
4. Check what the registry cannot know in advance: the PE labels of the new files (`hubble_constant`
   prefers `IMRPhenomXPHM-SpinTaylor`), any new PE distance-prior class (unknown ones stop the
   `prepare` stage with an error), the priors of a new cosmology paper (`h0_icarogw.PRIOR_SETS`), and
   events with missing products (`pe_supplements.py`).
5. Upload the new release to the S3 bucket and the Galaxy history if these repositories are used, then
   update the catalog tables of these pages, set `REGISTRY_CHECKED` to the date of the check, and
   regenerate the README and docs blocks (`python gwtc_analysis/gen_readme_cli_tables.py`), which name
   the catalogs of the version.

`tests/test_catalog_registry.py` checks the consistency of the registry, and that a new entry reaches
every derived table.
