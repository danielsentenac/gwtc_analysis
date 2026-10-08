# How the redshift is obtained

A GW signal alone does **not** measure the redshift. It measures the luminosity distance \(D_L\), from
its amplitude, and the detector-frame masses \(m_\text{det} = (1+z)\, m_\text{src}\), from its
frequency evolution ([What a GW signal measures](gw-signals.md)). Cosmic expansion stretches the
waveform, which looks exactly like a heavier binary: \(z\) and the source-frame mass are perfectly
degenerate. The redshift always comes from an additional ingredient, and there are five ways to get
it.

## 1. Assume a cosmology (what the catalogs do)

In a flat ΛCDM cosmology,

\[
D_L(z) = (1+z)\,\frac{c}{H_0} \int_0^z \frac{dz'}{\sqrt{\Omega_m (1+z')^3 + 1 - \Omega_m}} .
\]

With \(H_0\) and \(\Omega_m\) fixed (the catalogs use the Planck 2015 values), each PE sample of
\(D_L\) is inverted into a redshift, and the source masses follow as
\(m_\text{src} = m_\text{det}/(1+z)\). The `redshift` and `mass_1_source` values of GWOSC and of the PE
files are obtained this way: they are **not measurements of the redshift**, and they depend on the
assumed \(H_0\).

### The Planck 2015 cosmology

It is the flat ΛCDM cosmology whose parameters were measured by the Planck satellite from the cosmic
microwave background, in its 2015 release
(Planck Collaboration 2016 [\[32\]](../references.md#ref-32)). The GW software uses it as a fixed
reference, in two versions:

| Name | H₀ (km/s/Mpc) | Ω_m | Definition | Where it appears |
|---|---|---|---|---|
| `Planck15` (astropy) | 67.74 | 0.3075, including massive neutrinos (Σm_ν = 0.06 eV) | Paper XIII, Table 4, TT,TE,EE+lowP+lensing+ext; with radiation (T_CMB = 2.7255 K) | redshifts and source-frame quantities computed with astropy |
| `Planck15_LAL` (LAL, bilby) | 67.90 | 0.3065 | Planck 2015 values as defined in LAL; no radiation (T_CMB = 0) | the PE distance prior of O4a (`UniformSourceFrame(cosmology='Planck15_LAL')`); the Ω_m fixed by the GWTC-4.0 cosmology analysis (LVK 2025 [\[29\]](../references.md#ref-29)) and by `hubble_constant` |

The 2015 values became the convention of the LVK software during O1–O2 and were kept so that the
catalogs stay consistent with one another. The newer Planck 2018 values
(Planck Collaboration 2020 [\[33\]](../references.md#ref-33): H₀ = 67.4, Ω_m = 0.315) change the
redshift inferred from a given distance by about 0.5%, far below the typical 20–40% distance
uncertainty of an event. When H₀ itself is measured, the reference cosmology enters only through the
PE prior, which is divided out exactly.

## 2. Identify the host galaxy (bright siren)

If an electromagnetic counterpart is observed, as for GW170817 and its kilonova in the galaxy NGC 4993
(LVK et al. 2017 [\[41\]](../references.md#ref-41)), the redshift comes from the spectrum of the host
galaxy (z ≈ 0.01), corrected for the galaxy's own peculiar velocity. Combined with the GW distance,
this measures \(H_0\) directly: 70 (+12 / −8) km/s/Mpc
(LVK et al. 2017 [\[27\]](../references.md#ref-27)).

## 3. Use a galaxy catalog statistically (dark siren)

Without a counterpart, every galaxy inside the three-dimensional localization volume of the event is a
possible host. Each contributes its redshift, weighted by its probability of being the host (its
position in the skymap, its luminosity, the completeness of the catalog), and many events together
constrain \(H_0\) (Gray et al. 2023 [\[26\]](../references.md#ref-26)).

Galaxies fainter than the catalog's limit are missing, so the redshift prior of each line of sight adds a
completeness correction: the galaxies the catalog misses, uniform in comoving volume. Where the catalog is
incomplete, the prior is nearly uniform and the event carries almost no redshift information from it. Since GWTC-3
the LVK fits the mass distribution at the same time, so a "dark siren" result combines the galaxy catalog with the
spectral siren below. In `gwtc_analysis`, the [galaxy_catalog](../modes/galaxy-catalog.md) mode builds the catalog
and `hubble_constant --galaxy-catalog` uses it; with GLADE+ the GWTC-4.0 result is reproduced
([validation](../modes/hubble-constant.md#validation-gwtc-40-power-law-peak-glade-k-band)), and the catalog narrows
the interval by only a few percent, because GLADE+ is nearly empty at the gigaparsec distances of the black holes.

## 4. Use features of the mass distribution (spectral siren)

If the source-frame masses have a feature at a fixed mass, such as the peak near 30 M☉ of the
binary black holes, its observed position \(m_\text{det} = (1+z)\, m_\text{peak}\) gives the redshift
statistically, over the whole population
(Taylor, Gair & Mandel 2012 [\[20\]](../references.md#ref-20);
Ezquiaga & Holz 2022 [\[22\]](../references.md#ref-22)). In practice the redshift of an individual
event is never computed: for each trial \(H_0\), the analysis converts every \(D_L\) into a \(z\) and
checks that the source-frame masses form one consistent population. This is the method of the
`hubble_constant` mode ([Hubble constant](spectral-siren.md)).

## 5. Break the degeneracy with matter effects (future, neutron stars)

Tidal deformation of neutron stars depends on their source-frame masses, through the equation of state
of dense matter, and not on \(1+z\). A well-measured tidal signal could therefore give the redshift of a
binary neutron star from the GW signal alone
(Messenger & Read 2012 [\[19\]](../references.md#ref-19)). This requires a known equation of state and
the sensitivity of the next generation of detectors.

## In gwtc_analysis

| Mode | Redshift |
|---|---|
| `rates` | approach 1: the redshifts of the injections and of the events are given in the Planck 2015 cosmology |
| `hubble_constant` | approach 4: \(z\) is a function of the trial \(H_0\), never fixed |
| `event_selection`, `catalog_statistics` | approach 1: the catalog values, for the Planck 2015 cosmology |
