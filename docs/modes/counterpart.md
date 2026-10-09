# counterpart

Where an **electromagnetic counterpart** and its host sit in the GW posterior of an event: on the sky, along the
line of sight, and against the orientation of the binary. It takes seconds. The Hubble constant of the event is
[`hubble_constant --method bright`](hubble-constant-bright.md); this mode is about the association and the
quantities the H₀ depends on.

```bash
gwtc_analysis counterpart                                   # GW170817 and AT2017gfo
gwtc_analysis counterpart --viewing-angle 20 3              # with the viewing angle of the radio jet
gwtc_analysis counterpart --event GW190521                  # candidate flare ZTF19abanrhr
gwtc_analysis counterpart --event GW190521 --ra 190 --dec 30   # another position, same host redshift
```

## What it computes

For each PE label of the event (the counterparts are those of the bright siren:
[registered events](hubble-constant-bright.md#registered-events)):

| Quantity | How |
|---|---|
| **Searched probability** of the counterpart's position | the credible level of the smallest sky region that contains it. From the LVK sky map of the event for this label (Zenodo archive of its catalog, ligo.skymap crossmatch); without a map (`--pe-file`, `--sky-map none`), a kernel density estimate on the PE samples' sky positions |
| **Distance along its line of sight** | the PE samples within `--sky-radius` of the position (all of them when the sky was fixed to the counterpart): median and 90% interval |
| **Distance of the host redshift** | \(d_L(z; H_0)\) for Planck (67.4) and SH0ES (73.0), with its percentile in the distance posterior: a host redshift far in the tails of the GW distance makes the association, or the cosmology, doubtful |
| **Viewing angle** | the angle between the line of sight and the total angular momentum, θ_JN folded to 0–90° (0° face-on) |
| **Distance–inclination degeneracy** | the distance samples against the viewing angle, each colored by the H₀ it implies at the host redshift, and a table by angle |
| **With `--viewing-angle MEAN SIGMA`** | an independent Gaussian constraint on the viewing angle (degrees), as weights on the samples: the distance and the angle with it, and the effective number of samples |

The searched probability from the samples converges slowly for a sky with fine structure: for GW190521 and
ZTF19abanrhr, 0.55 from 20,000 samples and 0.62 from 50,000, against 0.64 from the LVK sky map. The map is used
whenever the event's PE file comes from Zenodo.

## GW170817 and AT2017gfo

The GWTC-1 samples have the sky fixed to AT2017gfo: there is no searched probability, and all the samples are on
its line of sight.

| PE label | \(d_L\) (Mpc), 90% | Host at Planck H₀: 45.1 Mpc | Host at SH0ES H₀: 41.6 Mpc | Viewing angle, 90% |
|---|---|---|---|---|
| LowSpin | 40.0 (24.9–47.3) | percentile 82 | percentile 60 | 33° (9–62°) |
| HighSpin | 41.7 (29.1–47.5) | percentile 79 | percentile 49 | 28° (8–55°) |
| LowSpin, viewing angle 20 ± 3° | 44.5 (41.9–47.3) | | | 20° (15–25°) |

![Distance of GW170817 toward AT2017gfo](../img/modes/counterpart_distance_GW170817.png)

The host redshift is well inside the GW distance posterior for both values of H₀. The viewing-angle constraint
(the radio jet, 14–26° at 41 Mpc, Hotokezaka et al. 2019 [\[73\]](../references.md#ref-73), here 20 ± 3°) keeps
1,600 effective samples of 8,078 and removes the low-distance tail.

### The distance–inclination degeneracy

![Distance against viewing angle for GW170817](../img/modes/counterpart_degeneracy_GW170817.png)

| Viewing angle | Samples | \(d_L\) median | H₀ implied |
|---|---|---|---|
| 0–30° | 43% | 44.6 Mpc | 68 |
| 30–60° | 51% | 35.9 Mpc | 85 |
| 60–90° | 6% | 23 Mpc | 132 |

*The samples form one band: the distance falls as the orbit is seen more inclined, because an inclined binary is
fainter than a face-on one at the same distance. The inclination is measured, from all three detectors (θ_JN =
147°, 90%: 118–171°), but broadly: Virgo's response to GW170817 was small (the wave came close to one of its blind
directions), which constrained the sky position but carried little polarization information, and near face-on the
two polarizations are almost equal, (1 + cos²ι)/2 against cos ι. The upper tail of the bright-siren H₀ comes from
the inclined orbits; the shaded band is the jet constraint, which removes them.*

## GW190521 and ZTF19abanrhr

| Quantity | Value |
|---|---|
| Searched probability of ZTF19abanrhr (LVK sky map, IMRPhenomXPHM) | 0.64 |
| Samples within 3° of the flare | 3,436 of 148,527 |
| \(d_L\) along that line of sight | 4.64 Gpc (90%: 1.75–7.12) |
| Host at z = 0.438: Planck H₀, SH0ES H₀ | 2.51 Gpc (percentile 21), 2.31 Gpc (percentile 18) |
| Viewing angle | 36° (90%: 10–84°) |

![Distance of GW190521 toward ZTF19abanrhr](../img/modes/counterpart_distance_GW190521.png)

The flare is inside the 64% credible region of the sky, and the distance of its host for the usual values of H₀
falls on the lower of the two modes of the GW distance in that direction. Neither confirms nor excludes the
association, whose odds are 1 to 12 depending on the waveform model (Ashton et al. 2021
[\[75\]](../references.md#ref-75)); the bimodal distance is why the bright-siren H₀ of GW190521 is bimodal.

## Options

| Option | Default | Meaning |
|---|---|---|
| `--event` | `GW170817` | registered event (`GW170817`, `GW190521`) |
| `--ra`, `--dec` | the counterpart's | another position to test (degrees), with the same host redshift |
| `--pe-label` | the registry's, else all the labels (LowSpin first) | PE labels; the plots use the first one |
| `--pe-file` | the event's bundle or Zenodo file | another PE file (PESummary layout) |
| `--cache-dir`, `--pe-cache` | `.cache_gwosc`, that of `hubble_constant` | caches of the GW170817 bundle and of the Zenodo PE files |
| `--v-recession V SIGMA`, `--v-peculiar V SIGMA`, `--redshift Z SIGMA` | the registry's | the host redshift |
| `--sky-radius` | 3° | line of sight: samples within this angle of the position |
| `--sky-map` | the event's LVK map | FITS sky map for the searched probability; `none`: from the PE samples |
| `--viewing-angle MEAN SIGMA` | none | independent Gaussian constraint on the viewing angle, degrees |
| `--out-report`, `--out-summary`, `--plots-dir` | `counterpart.html`, `counterpart.tsv`, `counterpart_plots` | outputs |

## Outputs

- `--out-report` (HTML): the counterpart and host, the caveat of a candidate, the searched probability, the
  distance along the line of sight against the host redshift, the degeneracy and the constraint;
- `--out-summary` (TSV): one row per PE label with all the quantities above;
- `--plots-dir`: `distance_<event>.png` and `distance_inclination_<event>.png`.
