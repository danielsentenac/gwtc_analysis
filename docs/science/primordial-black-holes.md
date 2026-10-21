# Primordial black holes and sub-solar masses

Could some of the black holes seen by LIGO–Virgo–KAGRA be primordial, and part of the dark matter? This page
explains what primordial black holes are, why sub-solar masses are the cleanest test, what a counterpart would tell,
how PBH binaries form and merge, and which bounds the GW data already set. The masses themselves are explained in
[What a GW signal measures](gw-signals.md#masses-primary-secondary-and-what-the-signal-measures); the merger rate
behind the last bound in [Merger rates](merger-rates.md).

## What they are

A primordial black hole (PBH) would form in the first fraction of a second after the Big Bang,
when an exceptionally dense region of the hot, radiation-dominated Universe collapses under its own gravity,
long before the first stars. Its mass is roughly the mass inside the horizon at that moment,

\[
M \sim \frac{c^3 t}{G} \approx 2 \times 10^5\, M_\odot \left(\frac{t}{1\ \text{s}}\right),
\]

so the formation time sets the mass and **any mass is possible**: about 1 M☉ at \(t \approx 10^{-5}\) s (the
epoch of the quark–hadron transition), asteroid masses at \(10^{-18}\) s, 10⁵ M☉ at 1 s.

## Why they can be dark matter

The evidence for dark matter is only gravitational: something cold, dark and
collisionless, about five times as abundant as ordinary matter. Ordinary (baryonic) matter is fixed at ~5% of
the cosmic budget by the light elements made in the first minutes and by the microwave background, so no object
made from baryons afterwards (faint stars, planets, the black holes of the GW catalogs) can be the dark matter.
A PBH escapes that argument: it forms *before* the light elements, from radiation, and is not counted as
baryons. Once formed it is massive, dark, collisionless and stable (if heavier than ~10¹² kg, below which it
has evaporated by Hawking radiation). Their share is written \(f_\text{PBH}\), the fraction of the dark matter
they make up.

## Why "sub-solar" comes up: two different reasons

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

## Counterparts of sub-solar mergers

Whether a sub-solar binary shines depends on its origin, and the
search for a counterpart is what separates the two:

| Origin | Counterpart | Why |
|---|---|---|
| Primordial black holes | none: the binary is dark | two black holes with no gas or disk around them emit no light |
| Neutron-star fragments of a collapsar disk [\[103\]](../references.md#ref-103) | yes, inside an explosion | they form and merge inside a collapsing massive star: a long gamma-ray burst and a supernova with a kilonova inside it (a "superkilonova") |

A transient in the localization points to the collapsar channel, a new kind of neutron-star merger rather than
dark matter; a dark field in a well-covered map points to primordial black holes. The dark-matter question thus
becomes partly an optical one. Two practical points for the follow-up:

- **They are seen only nearby.** The detection range scales roughly as \(\mathcal{M}^{5/6}\): against 160 Mpc for
  1.4 + 1.4 M☉ in O4, about 70 Mpc for 0.5 + 0.5 M☉ and 30 Mpc for 0.2 + 0.2 M☉. A counterpart would be close,
  hence bright, but such events are rare and none is confirmed.
- **The alerts flag them.** Since O4c every LVK alert carries HasSSM, the probability that a component is below
  1 M☉, and a chirp-mass bin whose lowest bin (up to 0.87 M☉) marks sub-solar systems. S251112cm, from the
  dedicated sub-solar search, was followed up by Rubin (GCN 43257); a dark result is informative only if the
  map was well covered.

## How PBH binaries form and merge

Bounds on f_PBH from GW data compare the merger rate a PBH population would have with the rate allowed by the data.
That rate depends on how the binaries form. The bounds on this page use two channels: **early binaries**, formed in
the radiation era, and **late binaries**, formed by capture in PBH clusters today (Phukon et al. [\[125\]](../references.md#ref-125)). The early
ones dominate the bounds below about 50 M☉, and their history is unusual.

**Close, yet wide.** Two PBHs born at a distance x first move apart with the expansion of the Universe. They stop,
turn around and fall towards each other only when their own gravity wins, i.e. when the density of the pair, ~m/x³,
exceeds that of the radiation around it. The pair's density falls as (1 + z)³ and the radiation's as (1 + z)⁴, so
a closer pair decouples earlier, when the Universe was smaller and had stretched it less. The separation at
decoupling therefore grows as x⁴, and the orbit that results has the semi-major axis
(Sasaki et al. 2016 [\[127\]](../references.md#ref-127); Raidal et al. 2019 [\[123\]](../references.md#ref-123))

\[
a \approx \frac{\alpha}{f_\text{PBH}}\,\frac{x^4}{\bar x^{3}},\qquad \alpha \approx 0.1,
\]

where x̄ is the mean PBH spacing at matter–radiation equality. "Close" and "wide" refer to two different scales:

- **"close"** means closer than the *average* distance between PBHs;
- **"wide"** means much wider than the ~0.2 AU (1 AU = 1.5 × 10⁸ km, the Earth–Sun distance) from which two 30 M☉
  black holes on a circular orbit merge within the age of the Universe. That limit comes from the Peters merger time,
  t = (5/256) c⁵a⁴ / (G³m₁m₂M) [\[99\]](../references.md#ref-99).

The average spacing is enormous. For 30 M☉ PBHs with f_PBH = 10⁻³, the number density today is ~10⁶ per Mpc³,
giving a spacing of ~10 kpc comoving, or ~3 pc ≈ 600 000 AU at equality (z ≈ 3400). With the fourth power:

| Initial distance x / x̄ | Orbit a (30 M☉, f_PBH = 10⁻³) |
|---|---|
| 0.1 | ~6000 AU |
| 0.03 | ~50 AU, a bit wider than Neptune's orbit and still 250 times the circular limit |

**Eccentric.** The pair falls almost straight in. Only the small tidal pull of the neighbouring PBHs and of the
surrounding matter gives it some sideways motion, so the orbit is nearly radial, with eccentricity e very close to 1.

**How a 50 AU orbit merges.** The size of the orbit is not its closest distance: the pericentre is a(1 − e) and the
apocentre a(1 + e). Gravitational-wave power grows roughly as 1/r⁵, so nearly all the energy leaves in a burst at
each pericentre passage. Each burst lowers the apocentre much more than the pericentre, so the orbit shrinks and
becomes round over hundreds of millions of orbits (45 years each at the start for 30 + 30 M☉ at 50 AU). For a very
eccentric orbit the merger time is

\[
t \approx \frac{3}{85}\,\frac{c^5 a^4}{G^3 m_1 m_2 M}\,(1-e^2)^{7/2}.
\]

Going from 0.2 to 50 AU multiplies a⁴ by ~3 × 10⁹, which is compensated if (1 − e²)^(7/2) ≈ 2 × 10⁻¹⁰, i.e.
**e ≈ 0.9991**. That orbit swings from ~100 AU in to a pericentre of ≈ 0.04 AU (6 million km), inside the 0.2 AU
circular limit. At a = 100 AU, a merger within the age of the Universe needs e ≳ 0.9996.

**Circular in the band.** By 10–20 Hz the orbit has shrunk to a few hundred kilometres and lost its eccentricity
(at small e, e ∝ f^(−19/18)). The catalog analyses assume circular orbits anyway: their posterior samples carry
neither eccentricity nor pericentre. An eccentric analysis would measure only the residual e at a reference
frequency. The e ≈ 0.999 birth orbit leaves no trace, so eccentricity cannot tell an early PBH binary from a
stellar one; a measured in-band eccentricity would point instead to a recent capture in a cluster or an AGN disk.

### The suppression factor S

Without correction, the early-binary rate (Hütsi et al. 2021 [\[124\]](../references.md#ref-124), eq. 2.1),

\[
R \approx 1.6\times10^{6}\ \text{Gpc}^{-3}\,\text{yr}^{-1}\, f_\text{PBH}^{53/37}
\left(\frac{t}{t_0}\right)^{-34/37} \left(\frac{M}{M_\odot}\right)^{-32/37} \eta^{-34/37}\, S,
\]

would count every binary that forms as surviving untouched until it merges. S < 1 is the fraction that does. Two
effects reduce it:

- **S₁, disruption at formation.** The other PBHs around the pair, and the ordinary matter fluctuations (rescaled
  variance σ_M ≈ 0.004), pull on it. They change its orbit, which delays the merger, or break it before it merges.
  The more PBHs there are (larger f), the more neighbours each pair has.
- **S₂, disruption later.** Once structures form, PBHs gather in clusters, where a binary can meet a third black
  hole and be broken apart or changed. This too is stronger at large f.

At f_PBH = 1 only about 2 binaries in 1000 survive (S ≈ 0.002); at small f the binaries are isolated and S is closer
to 1. S moves the bounds directly: the limit on f is where the predicted rate equals the rate the data allow. A
smaller S means fewer mergers for a given f, hence a weaker limit; a factor 3 in S moves the limit by a factor of 3
to 6.

The LVK sub-solar papers [\[101\]](../references.md#ref-101), [\[122\]](../references.md#ref-122) use the fit S = 2.3 × 10⁻³ f^−0.65. It equals S₁S₂ at f = 1 but keeps
only the growth of S₂ as f falls, without the drop of S₁ at small f:

| f_PBH | S₁S₂ / LVK fit | Consequence |
|---|---|---|
| 1 | 1.0 | identical |
| 0.1 | 1.3 | close |
| 0.01 | 2.9 | the LVK limit is ~6 times weaker at 1 M☉ (0.018–0.032 against 0.003–0.006) |
| ≲ 10⁻⁴ | — | the fit exceeds 1, more binaries surviving than formed: unphysical |

The [primordial_black_holes](../modes/primordial-black-holes.md) mode offers both (`--suppression`):

- `lvk` reproduces the published sub-solar limits;
- `hutsi` (the default) is the only valid choice above 1 M☉, where the limits lie at f ~ 10⁻⁵–10⁻³.

The choice changes no data, only the model of how many early binaries survive: it is the largest theoretical
uncertainty of these bounds.

## What the GW data say

| Mass range | Data | Bound on \(f_\text{PBH}\) |
|---|---|---|
| 0.2–1 M☉ (sub-solar) | LVK O4b search combined with O4a, no significant candidate; rates below 50–5000 Gpc⁻³ yr⁻¹ [\[122\]](../references.md#ref-122) (O4a alone: 110–10 000 [\[101\]](../references.md#ref-101)) | ≤ 1.8–3.2% at 1 M☉, ≤ 10–18% at 0.35 M☉ (binaries formed in the early Universe; the range spans the three search pipelines; O4a alone: 7% and 40%); binaries formed late: f_PBH < 1 only above 0.62–0.74 M☉ |
| 1–100 M☉ | the observed BBH merger rate, which no PBH population may exceed (O4a, [\[102\]](../references.md#ref-102)) | ~10⁻² at 1 M☉, ~10⁻³ at 10–30 M☉, ~10⁻⁴ at 100 M☉ |

The second bound rests on the measured BBH rate, the quantity of the [rates](../modes/rates.md) mode: whatever
fraction of the observed mergers is primordial, PBHs cannot produce more mergers than are seen. The stochastic
background gives a bound about 1000 times weaker.

The [primordial_black_holes](../modes/primordial-black-holes.md) mode computes both kinds of bound from the public
data: the sub-solar searches (it reproduces the LVK limits above) and the GWTC injections and events from 1 to
150 M☉.
