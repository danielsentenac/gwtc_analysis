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
| Bright siren | an electromagnetic counterpart and its host galaxy | GW170817 (Abbott et al. 2017 [\[27\]](../references.md#ref-27)) |
| Dark siren | a statistical association with the galaxies of a catalog | Gray et al. 2023 [\[26\]](../references.md#ref-26) |
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

## The stochastic background

The catalogs list the mergers loud enough to be detected one by one. Most mergers in the universe are too
distant and too faint: summed over the whole universe, a binary black hole merges every few minutes and a binary
neutron star every few tens of seconds, far more often than the detectors resolve. Their signals add up, with
random phases and from every direction, into a **stochastic background**: not a signal with a shape, but an extra,
persistent noise-like strain common to all detectors. Other, more speculative sources would add to it:
gravitational waves from inflation, phase transitions or cosmic strings in the early universe. The background of
compact binaries is the one that must exist, given the merger rates the catalogs measure.

**What is measured.** The strength of a background is its energy density per logarithmic frequency interval, as
a fraction of the critical density that makes the universe flat, Ω_GW(f): a pure number, usually quoted at 25 Hz,
where the LIGO–Virgo network is most sensitive to it.

**Why its spectrum rises.** During the inspiral a binary radiates an energy per unit frequency dE/df ∝ f^(−1/3);
weighted by f, as Ω_GW is, this gives a background growing as **f^(2/3)** up to the frequencies where the
redshifted binaries merge. The searches therefore look for a power law of index 2/3 for compact binaries, and of
index 0 (a flat spectrum) for many cosmological models.

**Continuous or "popcorn".** Binary-neutron-star signals last minutes in band and overlap, so their background is
nearly continuous and Gaussian. Binary-black-hole signals last seconds and rarely overlap: their background is a
sequence of faint, separate events, a **popcorn** background, which searches designed for Gaussian noise treat
less efficiently.

**How it is searched for: cross-correlation.** A background cannot be told apart from noise in one detector. With
two detectors, the noise of each is independent, but the background is common to both: multiplying their strains
and averaging over time, the noise products average away while the background adds up. The signal-to-noise ratio
of the correlation grows as the square root of the observing time, so the sensitivity improves with every run.
Two detectors see the same background only if they are close and similarly oriented compared with the wavelength;
the **overlap reduction function** measures this, and for the LIGO pair it falls above a few tens of hertz, which
is why the searches are most sensitive around 20–30 Hz.

**The same idea at other frequencies.** Pulsar timing arrays use the arrival times of radio pulses from dozens of
millisecond pulsars as a galaxy-sized detector, sensitive at nanohertz frequencies. In 2023 they reported evidence
for a background there, most likely from supermassive black-hole binaries; their cross-correlation between pairs
of pulsars plays the role of the overlap reduction function.

The background of compact binaries has not been detected yet in the LIGO–Virgo–KAGRA band. Its prediction from
the catalog, the formula behind it and the comparison with the current upper limits are in the
[stochastic](../modes/stochastic.md) mode.
