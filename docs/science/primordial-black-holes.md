# Primordial black holes and sub-solar masses

Could some of the black holes seen by LIGO–Virgo–KAGRA be primordial, and part of the dark matter? This page
explains what primordial black holes are, why sub-solar masses are the cleanest test, what a counterpart would tell,
and which bounds the GW data already set. The masses themselves are explained in
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

## What the GW data say

| Mass range | Data | Bound on \(f_\text{PBH}\) |
|---|---|---|
| 0.2–1 M☉ (sub-solar) | LVK O4b search combined with O4a, no significant candidate; rates below 50–5000 Gpc⁻³ yr⁻¹ [\[122\]](../references.md#ref-122) (O4a alone: 110–10 000 [\[101\]](../references.md#ref-101)) | ≤ 1.8–3.2% at 1 M☉, ≤ 10–18% at 0.35 M☉ (binaries formed in the early Universe; the range spans the three search pipelines; O4a alone: 7% and 40%); binaries formed late: f_PBH < 1 only above 0.62–0.74 M☉ |
| 1–100 M☉ | the observed BBH merger rate, which no PBH population may exceed (O4a, [\[102\]](../references.md#ref-102)) | ~10⁻² at 1 M☉, ~10⁻³ at 10–30 M☉, ~10⁻⁴ at 100 M☉ |

The second bound rests on the measured BBH rate, the quantity of the [rates](../modes/rates.md) mode: whatever
fraction of the observed mergers is primordial, PBHs cannot produce more mergers than are seen. The stochastic
background gives a bound about 1000 times weaker.
