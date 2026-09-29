# What a GW signal measures

This page gives the physical background that the [merger rates](merger-rates.md) and the
[Hubble constant](spectral-siren.md) rely on.

## Amplitude and frequency

A compact binary emits a "chirp": a wave whose frequency and amplitude grow until the merger. To
leading (Newtonian) order, its evolution is set by the **chirp mass**

\[
\mathcal{M} = \frac{(m_1 m_2)^{3/5}}{(m_1 + m_2)^{1/5}},
\qquad
\dot f = \frac{96}{5}\,\pi^{8/3} \left(\frac{G\mathcal{M}}{c^3}\right)^{5/3} f^{11/3},
\]

and its strain amplitude, for a source at luminosity distance \(D_L\), is

\[
h \propto \frac{(G\mathcal{M})^{5/3} (\pi f)^{2/3}}{c^4\, D_L} \times (\text{orientation and antenna factors}).
\]

- **The frequency evolution** gives the masses: the chirp mass from the inspiral, the mass ratio and
  spins from higher-order effects and from the merger and ringdown.
- **The amplitude** falls as \(1/D_L\). Once the masses are known from the phase evolution, the
  amplitude gives the distance directly, with no cosmological model and no distance ladder: GW
  sources are **standard sirens** (Schutz 1986 [\[17\]](../references.md#ref-17),
  Holz & Hughes 2005 [\[18\]](../references.md#ref-18)). The inclination of the orbit is
  partly degenerate with the distance, which is why distances are often uncertain by tens of percent.

## Redshift and the detector frame

The expansion of the Universe stretches the signal on its way: every time scale is multiplied by
\(1+z\), where \(z\) is the redshift, and every frequency divided by \(1+z\). A binary of masses
\(m_\text{src}\) at redshift \(z\) produces exactly the signal of a binary of masses

\[
m_\text{det} = (1+z)\, m_\text{src}
\]

at rest. The detector measures the **detector-frame** masses \(m_\text{det}\) (also called redshifted
masses); the physical **source-frame** masses require the redshift.

The signal alone cannot give \(z\): a heavy nearby binary and a lighter distant one can produce the
same signal. The distance is measured, but converting it to a redshift requires a cosmology.

## Distance and cosmology

In a flat ΛCDM cosmology, the luminosity distance of a source at redshift \(z\) is

\[
D_L(z) = (1+z)\,\frac{c}{H_0} \int_0^z \frac{dz'}{\sqrt{\Omega_m (1+z')^3 + 1 - \Omega_m}}
\;\approx\; \frac{c\,z}{H_0} \quad (z \ll 1),
\]

where \(H_0\) is the Hubble constant and \(\Omega_m\) the matter density. So:

- a measured \(D_L\) gives \(z\) **only for an assumed \(H_0\)** (and \(\Omega_m\));
- source-frame masses, and every population property in the source frame, depend on the assumed
  cosmology. The catalogs quote them for the [Planck 2015 cosmology](redshift.md#the-planck-2015-cosmology).

Conversely, if the redshift of GW sources can be found by another route, the \(D_L\)–\(z\) relation
measures \(H_0\):

| Method | Where the redshift comes from | Example |
|---|---|---|
| Bright siren | an electromagnetic counterpart and its host galaxy | GW170817 (Abbott et al. 2017 [\[25\]](../references.md#ref-25)) |
| Dark siren | a statistical association with the galaxies of a catalog | Gray et al. 2023 [\[24\]](../references.md#ref-24) |
| Spectral siren | features in the source-frame mass distribution of the population | [Hubble constant](spectral-siren.md) |

## Comoving volume and the rate of mergers

The number of mergers observed per unit time from redshifts between \(z\) and \(z+dz\) is

\[
\frac{dN}{dt_\text{det}\,dz} = \frac{R(z)}{1+z}\,\frac{dV_c}{dz},
\]

where \(R(z)\) is the merger rate per unit comoving volume and per unit source-frame time, \(V_c\) the
comoving volume, and the factor \(1/(1+z)\) the time dilation between the source and the detector.
This is the redshift distribution used by the [rates](merger-rates.md) and
[spectral-siren](spectral-siren.md) analyses.
