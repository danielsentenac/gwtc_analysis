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

## Masses: primary, secondary, and what the signal measures

The two objects of a binary are labelled by mass, not by history: the **primary** is the heavier one, of mass
\(m_1\), the **secondary** the lighter one, \(m_2\), so that \(m_1 \ge m_2\) always. The catalogs and the PE
files give them with the combinations below:

| Quantity | Definition | PE column | GW150914 (source frame) | GW170817 (source frame) |
|---|---|---|---|---|
| Primary mass | the heavier object, \(m_1\) | `mass_1_source` | 35.6 (+4.7 / −3.1) M☉ | 1.46 (+0.12 / −0.10) M☉ |
| Secondary mass | the lighter object, \(m_2 \le m_1\) | `mass_2_source` | 30.6 (+3.0 / −4.4) M☉ | 1.27 (+0.09 / −0.09) M☉ |
| Total mass | \(M = m_1 + m_2\) | `total_mass_source` | ≈ 66 M☉ | ≈ 2.7 M☉ |
| Mass ratio | \(q = m_2 / m_1\), between 0 and 1 (some files also give \(1/q \ge 1\)) | `mass_ratio` | ≈ 0.86 | ≈ 0.87 |
| Chirp mass | \(\mathcal{M} = (m_1 m_2)^{3/5} / M^{1/5}\) | `chirp_mass_source` | 28.6 (+1.7 / −1.5) M☉ | 1.186 ± 0.001 M☉ |
| Remnant (final) mass | the black hole left after the merger, \(M_f = M - E_\text{rad}/c^2\) | `final_mass_source` | 63.1 (+3.4 / −3.0) M☉ | — (not a black-hole binary fit) |

*Medians and 90% intervals, GWTC-1 (GWOSC event metadata).*

- **The signal measures combinations, not the two masses separately.** The inspiral phase is set by the chirp
  mass, which is therefore measured best: to 0.1% for GW170817, against ±7% for each component. The mass ratio
  enters only at higher order, and is partly degenerate with the spins aligned with the orbit (an unequal-mass
  binary with aligned spins can mimic a more equal one), so \(m_1\) and \(m_2\) are strongly correlated and
  each is uncertain. For heavy black holes the detectors see few inspiral cycles and mostly the merger and
  ringdown, which constrain the total mass rather than the chirp mass.
- **Detector frame or source frame.** The signal gives the redshifted masses \((1+z)\,m\)
  ([below](#redshift-and-the-detector-frame)). PE files give both: `mass_1`, `chirp_mass`… in the detector
  frame, `mass_1_source`, `chirp_mass_source`… in the source frame, computed with the distance and a reference
  cosmology. Population analyses use the source frame.
- **The remnant is lighter than the sum.** About 5% of the total mass is radiated as gravitational waves
  (GW150914: 66 M☉ in, 63.1 M☉ out, 3.1 M☉c² radiated in a fraction of a second). The remnant mass and spin
  come from fits to numerical relativity ([IMR](glossary.md)); the next section is about the spin.
- **Why population models are written for the primary.** The mass models of the catalogs (Power Law + Peak,
  Multi Peak, FullPop-4.0) describe the distribution of \(m_1\), and the secondary through the mass ratio,
  \(p(m_2 \mid m_1) \propto q^{\beta}\), or through a pairing function in FullPop-4.0. The primary is the better
  measured component and the one that carries the features: the peaks of black-hole masses near 10 and 35 M☉,
  the upper end of the distribution. These source-frame features are what the [spectral siren](spectral-siren.md)
  uses as a ruler. The secondary is often lighter than the features, and its distribution is broad.
- **Where the objects are.** Below about 2.5–3 M☉ an object is taken to be a neutron star (GW170817: both
  components); between about 3 and 5 M☉ lies the *lower mass gap*, where few compact objects are known
  (GW230529: 3.7 M☉ primary); above about 50 M☉ the *pair-instability gap*, where stars should leave no black hole
  (GW190521 and GW231123 have components in or near it, which points to black holes born in earlier mergers). The class of a
  binary (BNS, NSBH, BBH) is a statement about \(m_1\) and \(m_2\) against these limits.

## Primordial black holes and sub-solar masses

**What they are.** A primordial black hole (PBH) would form in the first fraction of a second after the Big Bang,
when an exceptionally dense region of the hot, radiation-dominated Universe collapses under its own gravity,
long before the first stars. Its mass is roughly the mass inside the horizon at that moment,

\[
M \sim \frac{c^3 t}{G} \approx 2 \times 10^5\, M_\odot \left(\frac{t}{1\ \text{s}}\right),
\]

so the formation time sets the mass and **any mass is possible**: about 1 M☉ at \(t \approx 10^{-5}\) s (the
epoch of the quark–hadron transition), asteroid masses at \(10^{-18}\) s, 10⁵ M☉ at 1 s.

**Why they can be dark matter.** The evidence for dark matter is only gravitational: something cold, dark and
collisionless, about five times as abundant as ordinary matter. Ordinary (baryonic) matter is fixed at ~5% of
the cosmic budget by the light elements made in the first minutes and by the microwave background, so no object
made from baryons afterwards (faint stars, planets, the black holes of the GW catalogs) can be the dark matter.
A PBH escapes that argument: it forms *before* the light elements, from radiation, and is not counted as
baryons. Once formed it is massive, dark, collisionless and stable (if heavier than ~10¹² kg, below which it
has evaporated by Hawking radiation). Their share is written \(f_\text{PBH}\), the fraction of the dark matter
they make up.

**Why "sub-solar" comes up: two different reasons.**

1. **It is the mass range where a merger cannot come from a star.** Stars leave white dwarfs (below 1.4 M☉, but
   too large to merge in the LIGO–Virgo band), neutron stars (no lighter than ~1 M☉ from core collapse) and
   black holes (above the neutron stars, ≳ 3 M☉). A compact object *below about 1 M☉* merging at
   tens to thousands of hertz is therefore not an ordinary stellar remnant, and is the cleanest PBH signature.
   Above 1 M☉ a PBH binary looks like any other BBH: it can only be told apart statistically (spins near zero,
   a mass distribution or a redshift evolution unlike that of stars). Since 2024 the sub-solar test is no
   longer airtight: the disk of a collapsing massive star (a collapsar) may fragment into neutron stars down to
   ~0.1 M☉ that merge (Metzger, Hui & Cantiello 2024 [\[103\]](../references.md#ref-103)); unlike PBHs, that channel
   comes with light (a long gamma-ray burst and a kilonova inside a supernova).
2. **It is where PBHs could still be all of the dark matter.** Over most masses they are excluded as the bulk of
   it: microlensing of stars excludes ~10⁻¹⁰ to ~10 M☉ as the dominant component, GW merger rates cap them at
   ~10⁻³ of the dark matter for 10–30 M☉, and the gamma rays of evaporating black holes exclude masses below
   ~10⁻¹⁶ M☉. The window left open is **10⁻¹⁶ to 10⁻¹⁰ M☉** (10¹⁴–10²⁰ kg, asteroid masses, black holes
   smaller than an atom to 0.3 µm): sub-solar by ten orders of magnitude or more
   (Carr & Kühnel 2020 [\[104\]](../references.md#ref-104)). Such binaries would emit far above the band of
   ground-based detectors (the last orbit of two masses M radiates at ~4.4 kHz × M☉/(m₁ + m₂)), so
   LIGO–Virgo cannot test that window.

So the two meet only partly: LIGO–Virgo probe sub-solar masses down to ~0.1–0.2 M☉, where a detection would be
striking but PBHs can only be part of the dark matter; the masses where they could be all of it are out of
reach of any GW detector today.

**What the GW data say.**

| Mass range | Data | Bound on \(f_\text{PBH}\) |
|---|---|---|
| 0.2–1 M☉ (sub-solar) | LVK O4a search, no candidate; rates below 110–10 000 Gpc⁻³ yr⁻¹ [\[101\]](../references.md#ref-101) | ≤ 7% at 1 M☉, ≤ 40% at 0.35 M☉ (binaries formed in the early Universe) |
| 1–100 M☉ | the observed BBH merger rate, which no PBH population may exceed (O4a, [\[102\]](../references.md#ref-102)) | ~10⁻² at 1 M☉, ~10⁻³ at 10–30 M☉, ~10⁻⁴ at 100 M☉ |

The second bound rests on the measured BBH rate, the quantity of the [rates](../modes/rates.md) mode: whatever
fraction of the observed mergers is primordial, PBHs cannot produce more mergers than are seen. The stochastic
background gives a bound about 1000 times weaker.

## Spins: small going in, fast coming out

The black holes of the catalog spin little before they merge, but the black hole a merger leaves behind spins
fast. The spin is measured as the dimensionless \(\chi = cJ/(GM^2)\), from 0 (not rotating) to 1 (the maximum
for a [Kerr black hole](glossary.md)).

- **Going in.** The effective spin \(\chi_\text{eff}\), the mass-weighted spin along the orbit, clusters near
  zero: the [spin_population](../modes/spin-population.md) mode finds a mean of about 0.04 for equal masses and
  0.1–0.3 for unequal ones, with a narrow spread (σ ≈ 0.08). Individual spin magnitudes are mostly below about 0.5.
- **Coming out.** Merger remnants have \(\chi_f \approx 0.6\)–0.8; GW150914 left a 63 M☉ black hole with
  \(\chi_f = 0.69\) (+0.05 / −0.04) (GWTC-1 [\[1\]](../references.md#ref-1)).

**The remnant's spin comes from the orbit.** Just before the merger the two black holes orbit each other at
about half the speed of light. Part of that orbital angular momentum is radiated in the gravitational waves; the
rest cannot disappear and becomes the spin of the remnant. For two equal, non-spinning black holes, numerical
relativity gives \(\chi_f \approx 0.69\), set by the geometry of the last orbit (final-spin fits: Rezzolla et
al. 2008 [\[76\]](../references.md#ref-76), used by `catalog_statistics`). Unequal masses give less (the light
companion brings little angular momentum, and \(\chi_f \to 0\) as the mass ratio goes to 0); spins aligned with
the orbit give more, up to 0.9 and above, and spins against it less.

**In everyday terms**, a 60 M☉ remnant with \(\chi = 0.7\) has a horizon of about 150 km radius, which turns
about 110 times per second; in a naive picture its equator moves at about a third of the speed of light.

**Where this enters `gwtc_analysis`:**

- **Hierarchical mergers.** A black hole that is itself the remnant of an earlier merger carries
  \(\chi \approx 0.7\) into its next merger, far above the small spins of the population. Spins near 0.7, with
  masses above the pair-instability gap, are the signature looked for in dense clusters and AGN disks, where
  remnants can merge again. In the [spin_population](../modes/spin-population.md) mode they would appear as a
  tail of large \(|\chi_\text{eff}|\), as many negative as positive when the spins are isotropic; the
  Gaussian χ_eff model does not describe such a subpopulation, and GW231123, the most massive and one of the
  fastest-spinning binaries, is left out of its default event selection.
- **The area law.** The horizon area \(A = 8\pi (GM/c^2)^2 (1 + \sqrt{1-\chi^2})\) shrinks with the spin at
  given mass, so the fast remnant spin works against the area increase, and so does the mass radiated (about
  5% of the total). Area still grows because it scales as the mass squared: one black hole of mass close to
  \(m_1 + m_2\) has more area than two of masses \(m_1\) and \(m_2\). The [area_law](../modes/area-law.md)
  mode measures both sides for GW250114.
- **The remnant of each event.** [parameters_estimation](../modes/parameters-estimation.md) summarizes, for
  every PE label, the remnant mass and spin stored in the samples.

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

## The stochastic background of compact binaries

The catalogs list the mergers loud enough to be detected one by one. Most mergers in the universe are too
distant and too faint: summed over the whole universe, a binary black hole merges every few minutes and a binary
neutron star every few tens of seconds, far more often than the detectors resolve. Their signals add up, with
random phases and from every direction, into a **stochastic background**: not a signal with a shape, but an extra,
persistent noise-like strain common to all detectors. This background of mergers must exist, given the merger
rates the catalogs measure, and it is the one this page and the [stochastic](../modes/stochastic.md) mode deal with.

!!! note "Stochastic: the kind of signal, not its origin"
    "Stochastic" means random, persistent and unresolved, and also applies to gravitational waves from the early
    universe (inflation, phase transitions, cosmic strings). Those are a different subject: the background from
    standard inflation is about a million times too weak for ground-based detectors, and primordial gravitational
    waves are sought mainly in the polarization of the cosmic microwave background.

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
