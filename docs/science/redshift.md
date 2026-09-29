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

## 2. Identify the host galaxy (bright siren)

If an electromagnetic counterpart is observed, as for GW170817 and its kilonova in the galaxy NGC 4993
(LVK et al. 2017 [\[38\]](../references.md#ref-38)), the redshift comes from the spectrum of the host
galaxy (z ≈ 0.01), corrected for the galaxy's own peculiar velocity. Combined with the GW distance,
this measures \(H_0\) directly: 70 (+12 / −8) km/s/Mpc
(LVK et al. 2017 [\[25\]](../references.md#ref-25)).

## 3. Use a galaxy catalog statistically (dark siren)

Without a counterpart, every galaxy inside the three-dimensional localization volume of the event is a
possible host. Each contributes its redshift, weighted by its probability of being the host (its
position in the skymap, its luminosity, the completeness of the catalog), and many events together
constrain \(H_0\) (Gray et al. 2023 [\[24\]](../references.md#ref-24)).

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
