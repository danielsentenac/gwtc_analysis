# gwtc_analysis explained

Notes on the physics behind several modes of `gwtc_analysis`: the bright-siren (counterpart) analysis
(`counterpart.py`, `h0_bright.py`, sections 1–6), the area-law test (`area_law.py`, sections 7–10) and
the scalar-cloud test (`scalar_cloud.py`, `scalar_cloud_mode.py`, section 11), with the sources of the
formulas in section 12.

## 1. Viewing angle vs inclination

They describe the same geometry. The viewing angle is a **folded** version of the inclination
(`counterpart.py`, `viewing_angle()`):

```python
def viewing_angle(samples: dict) -> Optional[np.ndarray]:
    for k in ("theta_jn", "iota"):
        if k in samples:
            th = np.degrees(np.asarray(samples[k], float))
            return np.minimum(th, 180 - th)
    return None
```

**Inclination** is the raw parameter-estimation sample, from 0 to 180°. The posterior files provide
one of two parameters:

- `theta_jn`: the angle between the **total** angular momentum J and the line of sight. Used first.
- `iota`: the angle between the **orbital** angular momentum L and the line of sight, at the
  reference frequency. Used only if `theta_jn` is missing.

For low-spin binary neutron stars like GW170817, J and L point almost the same way, so the two are
nearly equal. When spins precess, they differ.

**Viewing angle** is `min(θ, 180° − θ)`, from 0 to 90° (0 = face-on, 90 = edge-on). Folding drops the
difference between seeing the binary from above (θ near 0°) and from below (θ near 180°). That
difference matters for the gravitational-wave polarization, but not for:

- the distance–inclination degeneracy, since the signal amplitude depends on |cos θ|;
- electromagnetic constraints such as a jet opening angle, quoted relative to the jet axis on either
  side.

So `--viewing-angle MEAN SIGMA` takes a value from 0 to 90° (checked in `h0_bright.py`). It is
compared with the folded samples, which become weights in `viewing_angle_weights()`.

The "distance–inclination degeneracy" wording (plot titles, `distance_inclination_{event}.png`, notes)
is the standard name from the literature; the plotted axis is the folded viewing angle.

## 2. Recession velocity

The **recession velocity** is how fast the host galaxy appears to move away from us, measured from the
redshift of its spectrum: v_r ≈ c·z. In the code it is the CMB-frame velocity of the host
(`Counterpart.v_recession`).

It is the sum of two parts:

- **Hubble flow** (v_H): the part due to the expansion of the universe. Only this part follows
  Hubble's law, v_H = H0 · d.
- **Peculiar velocity** (v_p): the galaxy's own motion through space, caused by the gravity of
  nearby structure (groups, clusters).

For a distant host the peculiar velocity is negligible. For a nearby host it is not: NGC 4993, the
host of GW170817, is about 40 Mpc away, so its peculiar velocity is roughly 10% of its total velocity.
`hubble_flow_redshift()` removes it, adding the uncertainties in quadrature:

```python
vh, sv = vr - vp, float(np.hypot(sr, sp))
return vh / C_KMS, sv / C_KMS, ...
```

With the values registered for NGC 4993, which match the GW170817 standard-siren H0 paper
(Abbott et al. 2017, Nature):

| Quantity | Value |
|---|---|
| v_r (recession, CMB frame) | 3327 ± 72 km/s |
| ⟨v_p⟩ (peculiar) | 310 ± 150 km/s |
| v_H = v_r − v_p | 3017 ± 166 km/s, so z_H ≈ 0.01006 |

This Hubble-flow redshift is combined with the gravitational-wave distance to give H0. The
peculiar-velocity uncertainty (±150 km/s) is a large part of the redshift error, more than the
measured velocity's (±72 km/s).

**CLI options**

- `--v-recession V SIGMA` and `--v-peculiar V SIGMA` override the registered values, in km/s.
- `--redshift Z SIGMA` skips the velocity correction and uses the given redshift directly. Use it for
  a distant host with a spectroscopic redshift, where the peculiar velocity doesn't matter.
- `--v-peculiar` on its own fails if the host is registered by its redshift, because there is no
  recession velocity to subtract it from.

## 3. The CMB frame

### Definition

The **CMB frame** is the reference frame where the cosmic microwave background looks the same in
every direction, apart from its tiny intrinsic fluctuations. It is the closest thing to a cosmic rest
frame: an observer at rest there is moving only with the expansion of the universe.

We are not at rest in it. The Sun moves at about **370 km/s** relative to the CMB, toward galactic
(l, b) ≈ (264°, 48°), in Crater/Leo. We see this as the **CMB dipole**: the CMB is about 3.4 mK hotter
in that direction and colder in the opposite one.

### It is a state of motion, not a place

The CMB frame is defined at every point in the universe. At any location it is the frame of a
**comoving observer**: one with no velocity of their own, only carried along by the expansion.

Why the CMB looks the same in all directions to such an observer:

- The CMB fills the whole universe almost uniformly. It was released everywhere at once, about
  380,000 years after the Big Bang (z ≈ 1100), when the universe became transparent.
- Wherever you are, you receive photons from a sphere around you, the *last-scattering surface*.
  Each observer has their own sphere. The early universe was very nearly uniform (the cosmological
  principle), so every such sphere looks nearly identical: 2.725 K in every direction, with
  fluctuations of about 10⁻⁵.

When you move, the Doppler effect changes this. Photons from ahead are blueshifted and those from
behind redshifted, so T(θ) ≈ T₀ (1 + (v/c) cos θ). For v ≈ 370 km/s, v/c ≈ 1.23 × 10⁻³ and
ΔT ≈ 3.4 mK. That is the measured dipole.

Analogy: rain on a windless day falls vertically and evenly when you stand still. Run, and it seems
to come from the front and hits you harder. The rain hasn't changed, only your motion through it.

An observer at rest in the CMB frame at our position would move at about 370 km/s *relative to the
Sun*, in the direction opposite to the apex. For that observer the dipole would vanish.

### For galaxy velocities

- Each galaxy has its own local CMB frame. Its **peculiar velocity** is its motion relative to it.
- A galaxy with no peculiar velocity moves only with the Hubble flow.
- A **CMB-frame recession velocity** is the cz we would measure if we were comoving observers
  ourselves, so our own 370 km/s is removed.
- Subtracting the galaxy's peculiar velocity (v_H = v_r − v_p) then leaves only the Hubble flow,
  which is what H0 relates to distance.

The code uses the CMB frame because Hubble's law describes motion relative to the cosmic rest frame,
and because the peculiar velocity it subtracts (from velocity-field reconstructions) is also given in
that frame. v_H = v_r − v_p only makes sense if both are in the same frame. **If you pass
`--v-recession` yourself, give a CMB-frame value**: a heliocentric velocity biases H0 by a few
hundred km/s.

### Where to find CMB-frame velocities

Catalogs and tools:

- **NED** (see section 6): each galaxy's "Velocity corrections" table has a "CMB" row, and its
  velocity-correction calculator converts a heliocentric velocity to every frame.
- **HyperLEDA** (leda.univ-lyon1.fr): the `v3k` parameter is the velocity relative to the CMB.
- **Cosmicflows-4 / Extragalactic Distance Database** (edd.ifa.hawaii.edu): CMB-frame velocities plus
  measured and reconstructed peculiar velocities, useful for `--v-peculiar`.
- **Group catalogs** (2MRS groups of Crook et al. 2007, Tully 2015): CMB-frame velocities averaged
  over a galaxy's group.

To compute it yourself, use the Planck 2018 dipole: 369.82 ± 0.11 km/s toward galactic
(l, b) = (264.02°, 48.25°). To first order:

```python
from astropy.coordinates import SkyCoord
import astropy.units as u, numpy as np

gal = SkyCoord(197.4487*u.deg, -23.3839*u.deg)              # NGC 4993
apex = SkyCoord(l=264.021*u.deg, b=48.253*u.deg, frame="galactic")
v_helio = 2930.0                                            # km/s, from the spectrum
v_cmb = v_helio + 369.82 * np.cos(gal.separation(apex))
```

For NGC 4993 this gives cos ≈ 0.84 and **v_CMB ≈ 3241 km/s** for the galaxy itself. The registered
3327 ± 72 km/s is, as far as I recall, the CMB-frame velocity of **NGC 4993's galaxy group** from the
GW170817 H0 paper; averaging over the group reduces the effect of the galaxy's own motion within it
(check the paper's references for the exact source). The value passed to `--v-recession` and the
`--v-peculiar` subtracted from it should be in the same frame and describe the same object, either
the galaxy or its group.

## 4. How Planck corrects for its own motion

**Planck's motion around the Sun (orbital dipole).** Planck orbited L2, so it moved with the Earth
around the Sun at about 30 km/s. This adds a dipole of about 270 μK whose direction turns over the
year. The orbit is known to cm/s, so the signal is predicted exactly. Planck used it as its **main
absolute calibrator**, converting detector output into kelvin (2015 and 2018 analyses, for both LFI
and HFI).

**The Sun's motion relative to the CMB (solar dipole).** The 370 km/s dipole (about 3.36 mK) is
fitted from the maps, and the standard maps are released with it subtracted. Its amplitude and
direction are the values used for the CMB-frame correction.

**Subtler effects.** Moving at v/c ≈ 10⁻³ also:

- **aberrates** the sky slightly, pushing apparent positions toward the apex;
- **Doppler-boosts** the small-scale fluctuations, making them look slightly brighter and smaller in
  that direction;
- adds a small kinematic quadrupole, of order (v/c)², frequency-dependent and corrected in the
  processing.

In 2013 Planck measured the aberration and boosting in the fluctuations themselves and found
v ≈ 384 ± 78 km/s in the same direction: an independent check that the dipole comes from our motion.

**Open question.** The kinematic reading of the dipole is an assumption; an intrinsic CMB dipole would
look the same. Counts of distant quasars and radio sources should show the same dipole, but some
studies (e.g. Secrest et al. 2021) find an amplitude about twice as large. This tension is not
resolved. At about 40 Mpc, everyone uses the standard kinematic interpretation.

## 5. The chain of motions (Milky Way, Local Group)

Our 370 km/s relative to the CMB is the vector sum of several nested motions:

| Motion | Speed | Direction |
|---|---|---|
| Earth around the Sun | ~30 km/s | changes over the year (Planck's calibrator) |
| Sun relative to its local stars (Local Standard of Rest) | ~20 km/s | roughly toward Hercules |
| Sun around the Galactic centre | ~230–250 km/s | l = 90°, the direction of Galactic rotation |
| Milky Way relative to the Local Group barycentre | ~80–120 km/s | roughly toward M31 |
| **Local Group relative to the CMB** | **~620 km/s** | (l, b) ≈ (272°, 30°) |

These don't simply add up. The Sun's rotation around the Galaxy (toward l = 90°) points roughly
*against* the Local Group's motion (toward l ≈ 272°), so it cancels part of the 620 km/s, leaving the
~370 km/s of the Sun relative to the CMB.

**No separate Milky Way correction is needed.** The CMB dipole is measured directly from the Sun's
position (once Planck's orbital motion is removed), so it already includes all of these motions. A
single heliocentric → CMB correction covers them. The intermediate frames are still used, and NED
lists them all: heliocentric, LSR, Galactocentric (GSR), Local Group and CMB.

**Why the Local Group moves at 620 km/s.** It is the gravity of large-scale structure:

- the pull of the Great Attractor / Norma cluster region, near the centre of the Laniakea
  supercluster;
- the pull of the more distant Shapley supercluster;
- a push away from the underdense Local Void.

NGC 4993's 310 km/s peculiar velocity has the same origin; both numbers come from the same
velocity-field reconstructions of the local universe.

## 6. NED

**NED**, the **NASA/IPAC Extragalactic Database** (https://ned.ipac.caltech.edu), is run by IPAC at
Caltech, funded by NASA, and public since 1990. It brings together the published data on objects
**outside the Milky Way**: galaxies, quasars, clusters, radio and X-ray sources.

For each object it collects and cross-matches data from thousands of catalogs and papers:

- **Names and positions**: all aliases (NGC 4993, ESO 508-G018, PGC 45657…) resolved to one object.
- **Redshifts and velocities**: the measured heliocentric redshift, plus a **"Velocity corrections"**
  table converting it to the LSR, Galactocentric, Local Group and **CMB** frames. This is where to
  look up a `--v-recession` value.
- **Redshift-independent distances**: Cepheids, tip of the red giant branch, Tully–Fisher,
  surface-brightness fluctuations, supernovae… Useful to compare with a gravitational-wave distance.
- **Photometry** from radio to gamma rays, images, morphology, Galactic extinction along the line of
  sight, and the bibliography for the object.

Useful tools:

- **Object search** by name or by position (cone search).
- A **velocity-correction calculator**: coordinates and a heliocentric velocity in, each frame out.
- A **cosmology calculator** (distances, lookback time).
- An API, also reachable from Python with `astroquery.ipac.ned`:

  ```python
  from astroquery.ipac.ned import Ned
  t = Ned.query_object("NGC 4993")   # position, redshift, velocity, type...
  ```

Related databases:

- **SIMBAD** (CDS, Strasbourg): objects inside and outside the Galaxy, strongest on stars. NED covers
  extragalactic objects only.
- **HyperLEDA** (Lyon): homogenized galaxy parameters, such as the `v3k` CMB velocity.

Positions and redshifts are commonly checked in both NED and SIMBAD.

## 7. The area-law mode (`area_law.py`)

The area-law mode does **not** compute the remnant mass itself. The LVK release for GW250114
(Zenodo 16877102, arXiv:2509.08054) already provides **area samples**, computed by the collaboration
for each time window. `area_law.py` loads them, pools the two ringdown codes, compares the areas and
scans over the windows. `kerr_area()` is only used for the optional full-signal (`--with-imr`)
consistency check.

**Which steps come from the paper.** Almost everything: the pooling, the Gaussian significance and
the scans. `gwtc_analysis` adds the direct probability P(A_f < A_i), and approximates the
area-change distribution with random pairs. Each step below is marked.

### 7.0 The physics

The horizon area of a Kerr black hole of mass M and spin χ is A = 8πM²[1 + √(1 − χ²)]. Hawking's area
law requires A_f ≥ A_i = A₁ + A₂.

For this to be a real test, the two areas must come from **separate parts of the signal**:

- **A_i, from the inspiral**: parameter estimation (NRSur7dq4) on data **truncated at t<** before the
  peak gives m₁, m₂, χ₁, χ₂. There are 16 truncations, from −250 M to 0 M (M = GM_total/c³ with the
  detector-frame total mass; see section 8).
- **A_f, from the ringdown**: data **starting at t>** after the peak are fit with damped sinusoids, the
  quasinormal modes. In Kerr, a mode's frequency and damping time depend only on (M_f, χ_f), so the
  fit gives the remnant mass and spin without any merger model. Ringdown times are in units of M_f,
  the remnant's mass. Two models are fit:
  - **220**: the fundamental mode;
  - **220 + 221**: adding the first overtone, which decays fast and is only valid from about 6 M_f.

The full-signal fit is not a test: its remnant comes from numerical-relativity fits, which obey the
law by construction.

### 7.1 Pooling the two ringdown codes (`pooled_remnant`)

**Origin: the paper.** Pooling in equal numbers, with the paper's seed 250114.

Two independent ringdown codes analysed the data:

- **`ringdown`**: provided at many start times and for both models;
- **pyRing**: provided only at 10.5 M_f, with the 220 mode.

The code takes n = min(N_ringdown, N_pyRing) = 2000 samples from each code, drawn **without
replacement** with the paper's seed 250114, and concatenates them into 4000 samples. Equal numbers
give the two codes equal weight, so neither dominates just because it has more samples. The fixed seed
reproduces the paper's exact subsample.

### 7.2 Headline result

**Origin: partly the paper.** The Gaussian significance S and its 4.4σ are the paper's; the direct probability P(A_f < A_i) is added by `gwtc_analysis`.

A_i at t< = −40 M is compared with the pooled A_f at 10.5 M_f, in two ways.

**Gaussian significance**: S = (⟨A_f⟩ − ⟨A_i⟩) / √(σ_f² + σ_i²). If both areas are Gaussian and
independent, A_f − A_i is Gaussian with mean ⟨A_f⟩ − ⟨A_i⟩ and width √(σ_f² + σ_i²); S says how many
widths the difference lies above 0. Only the means and standard deviations are used, because with a
few thousand samples the far tails are poorly sampled. At 4.4σ the one-sided tail is about 5 × 10⁻⁶,
so a few thousand samples hold essentially none of it.

**Direct probability**: P(A_f < A_i). The code draws 400,000 random pairs (one A_f, one A_i, **with**
replacement) and counts the fraction where the law is violated. This measures the tails directly, with
no Gaussian assumption. Its resolution is about 1/400,000 ≈ 2.5 × 10⁻⁶, close to the Gaussian 4.4σ
level, so it is reported alongside S rather than replacing it.

### 7.3 Relative area change

**Origin: partly the paper.** The distribution of (A_f − A_i)/A_i is the paper's; drawing 200,000 random pairs instead of all the pairs is an approximation of `gwtc_analysis`.

The code computes (A_f − A_i)/A_i over 200,000 random pairs and reports the median and the 5–95%
interval. The paper uses every possible pair (N_f × N_i); random pairs give nearly the same
percentiles for much less computation. The histogram is plotted against the prior on this ratio
implied by the PE priors (`area_change_prior.dat`), to show that the data, not the prior, drive the
result.

Both areas are in detector-frame masses, so the (1+z)² factor cancels in the ratio and in S.

### 7.4 The two scans

**Origin: the paper.** Both scans and the three numbers drawn from them (section 7.5).

**Scan over the inspiral truncation**: A_f is fixed (pooled, 10.5 M_f) and A_i varies over the 16
values of t<. A later truncation keeps more of the inspiral, so A_i is measured more precisely and S
rises.

**Scan over the ringdown start**: A_i is fixed at −40 M, and A_f varies with t> for each model.

- These scans use the `ringdown` code alone, since pyRing only has 10.5 M_f. So the 220 value at
  10.5 M_f in this scan differs slightly from the pooled headline.
- A later start loses signal (the ringdown decays exponentially), so A_f widens and S falls.
- An earlier start with 220 alone risks bias from overtones that haven't decayed yet. That's why
  220 + 221 is used from 6 M_f.

### 7.5 Summary numbers, compared with the paper

| Quantity | How it's computed | Published |
|---|---|---|
| Headline | S at −40 M vs the pooled A_f | 4.4σ |
| Minimum over the truncations | min of S in the truncation scan | 3.4σ |
| Earliest truncation above 5σ | the most negative t< with S ≥ 5 | −10 M |
| Overtone | S of 220 + 221 at 6 M_f vs A_i at −40 M | 3.6σ |

The "earliest truncation above 5σ" is the first t< that reaches 5σ; the code doesn't check that every
later truncation stays above it.

**Outputs**:

- `area_law.tsv`: the comparison table above;
- `area_law.truncation.tsv` and `area_law.ringdown.tsv`: the two scans;
- two plots (the area distributions with the ratio, and the significance scans) and an HTML report.

## 8. The M unit (t = −40 M, 10.5 M_f…)

### 8.1 What M is

In GR it's common to set G = c = 1, so a mass can be expressed as a time or a length:

- **as a time**: GM/c³ ≈ **4.93 μs per solar mass**;
- **as a length**: GM/c² ≈ **1.477 km per solar mass**. This is the `MSUN_KM` constant in
  `area_law.py`, which converts areas to km².

So "1 M" means GM/c³ for the binary's total mass M, and **t = −40 M means 40 × GM/c³ before the
peak**. The peak, t = 0, is the time of maximum strain amplitude, around the merger.

### 8.2 "G = c = 1" does not mean G equals c

G and c don't even have the same dimensions (G = 6.674 × 10⁻¹¹ m³ kg⁻¹ s⁻², c = 2.998 × 10⁸ m s⁻¹).
"G = c = 1" is a **choice of units**: the units of length, time and mass are picked so that each
constant has the numerical value 1. Each constant then becomes a conversion factor:

1. **c = 1 links time and length.** Time is measured in meters, the distance light travels in that
   time: 1 s ↔ 299,792,458 m. Speeds become dimensionless, and light has speed 1.
2. **G = 1 links mass and length.** The factor G/c² converts a mass into a length:
   1 kg ↔ G/c² = 7.42 × 10⁻²⁸ m. For the Sun: 1.989 × 10³⁰ kg × 7.42 × 10⁻²⁸ m/kg ≈ **1477 m**.
   Dividing that length by c gives a time: 1477 m / c ≈ **4.93 μs**.

With G, c and a mass M, the only time you can build is GM/c³, and the only length GM/c². Writing
"40 M" means 40 × GM/c³: in SI, nothing is set to 1, it's just dimensional analysis.

Analogy: "the town is 2 hours away" uses a speed to express a distance as a time; nobody concludes
that hours equal kilometers. (Imperfect, since c is a true universal constant, but the logic of a
constant as a conversion factor is the same.)

Setting G = c = 1 fixes two of the three base units. The third is fixed by the problem: lengths and
times are expressed in units of the binary's mass M.

### 8.3 What 4.93 μs corresponds to

It isn't a physical event on its own: it is the **natural time scale of gravity for one solar mass**,
the time light takes to cross the gravitational radius GM/c² ≈ 1.477 km. Strong-gravity processes
around a mass M happen on simple multiples of GM/c³:

- **Black hole size**: the Schwarzschild radius is 2GM/c² ≈ 2.95 km per solar mass; light crosses it
  in ≈ 9.9 μs.
- **Orbits**: Kepler's law T = 2π√(r³/GM), with r = n × GM/c², gives T = 2π n^{3/2} × GM/c³.

For a non-spinning black hole:

| Process | In units of M | For 1 M_☉ | For GW250114 (~70 M_☉) |
|---|---|---|---|
| Light crossing GM/c² | 1 M | 4.93 μs | ≈ 0.35 ms |
| Orbit at the ISCO (r = 6M) | 2π·6^{3/2} ≈ 92 M | 0.45 ms | ≈ 32 ms, so f_GW ≈ 63 Hz |
| 220 ringdown period | 2π/0.374 ≈ 16.8 M_f | 83 μs, so f ≈ 12 kHz | ≈ 4 ms, so f ≈ 250 Hz* |
| 220 damping time | ≈ 11.2 M_f | 55 μs | ≈ 4 ms |

\* For GW250114's spinning remnant (χ_f ≈ 0.7), the 220 frequency is higher in units of M_f, about
0.53/M_f instead of 0.374/M_f. Approximate values from memory; the remnant rings around 250 Hz.

- **ISCO** is the innermost stable circular orbit. The gravitational-wave frequency is twice the
  orbital frequency, so f_GW,ISCO ≈ 4400 Hz × (M_☉/M).
- **Everything scales linearly with mass**: a 10 M_☉ black hole rings at about 1.2 kHz, a 100 M_☉
  one at about 120 Hz, a 10⁶ M_☉ one at about 0.01 Hz (the LISA band). This is why LIGO and Virgo,
  sensitive to tens to hundreds of Hz, see stellar-mass black-hole mergers.

### 8.4 In seconds for GW250114

M is the **detector-frame** total mass: about 66 M_☉ in the source frame, times (1 + z) with z ≈ 0.09,
so about 70 M_☉ (approximate values from memory of the paper). M_f is a few percent smaller, since
about 5% of the mass is radiated as gravitational waves.

| Time | In ms |
|---|---|
| 1 M | ≈ 0.35 ms |
| −40 M (headline truncation) | ≈ 14 ms before the peak |
| −250 M (earliest truncation) | ≈ 88 ms before the peak |
| 1 M_f | ≈ 0.33 ms |
| 10.5 M_f (220 start) | ≈ 3.5 ms after the peak |
| 6 M_f (220 + 221 start) | ≈ 2 ms after the peak |

40 M ≈ 14 ms is less than half an ISCO orbit period (~32 ms): t< = −40 M cuts the data in the very
last cycles before the merger.

### 8.5 Why use M instead of seconds

**Vacuum GR has no built-in scale.** For a given mass ratio and spins, a binary of 10 M_☉ and one of
100 M_☉ produce the **same waveform shape**; the heavier one is just stretched by a factor of 10 in
time. In units of M, the waveform depends only on the mass ratio and spins:

- Numerical-relativity simulations, and surrogates like NRSur7dq4, are computed with M = 1 and
  rescaled to the physical mass.
- Times in M mean the same thing for every event: −40 M is always at the same stage of the late
  inspiral, whatever the mass.
- The ringdown scales with M_f. The 220 damping time is about 12 M_f for this remnant's spin (value
  from memory; it depends on the spin), so starting at 10.5 M_f means about one damping time after the
  peak. By then the overtones, which decay several times faster, have mostly died out. That's why the
  220 + 221 model is needed to start as early as 6 M_f.

## 9. What σ (sigma) means in the results

In `gwtc_analysis`, σ is always a **standard deviation**, but it's used in two ways.

### 9.1 σ as a width: the uncertainty of one measurement

σ_f and σ_i in the area-law formula are the standard deviations of the posterior samples
(`final.std()` and `initial.std()` in `gaussian_significance`):

  σ = √( (1/N) Σ (A_k − ⟨A⟩)² )

The samples are draws from the **Bayesian posterior**, the probability distribution of the parameter
given the data. σ measures the width of that distribution, so it is the measurement uncertainty: a
small σ means the data pin the parameter down tightly.

For a Gaussian distribution:

| Range | Probability inside |
|---|---|
| ±1σ | 68.3% |
| ±1.645σ | 90% |
| ±2σ | 95.4% |
| ±3σ | 99.7% |

GW papers, and most outputs of `gwtc_analysis`, quote the **median with a 90% credible interval**
instead: the 5th to 95th percentiles of the samples, like `ratio_q` in the area-law mode or
`viewing_angle_low_90`/`high_90` in the counterpart mode. These intervals use no Gaussian assumption,
so they stay correct for skewed posteriors such as distances and spins. For a Gaussian, the 90%
interval is ±1.645σ.

### 9.2 "4.4σ" as a significance: how far a result is from a hypothesis

  S = (⟨A_f⟩ − ⟨A_i⟩) / √(σ_f² + σ_i²)

The difference A_f − A_i has its own width √(σ_f² + σ_i²), because the two uncertainties are
independent and add in quadrature. S counts how many of these widths separate the measured difference
from 0, the boundary of a violation of the area law.

The value is turned into a probability using the one-sided tail of a Gaussian, P(A_f < A_i):

| Significance | One-sided probability of being wrong | |
|---|---|---|
| 1σ | 0.16 | |
| 3σ | 1.3 × 10⁻³ | |
| 3.4σ | 3.4 × 10⁻⁴ | minimum over the truncations |
| 3.6σ | 1.6 × 10⁻⁴ | 220 + 221 from 6 M_f |
| **4.4σ** | **5.4 × 10⁻⁶** | headline result |
| 5σ | 2.9 × 10⁻⁷ | the usual "discovery" threshold |

So "4.4σ" means that, if the distributions were Gaussian, the probability that A_f < A_i would be
about 5 in a million. It is one-sided because only one direction violates the law.

### 9.3 Caveats

- The conversion assumes Gaussian distributions. Real posteriors have non-Gaussian tails, which is why
  the code also gives the direct P(A_f < A_i) from random pairs of samples.
- A Bayesian σ is not the same as the **detection significance** of a GW search, which is a
  false-alarm rate measured on background noise. The 4.4σ here is about a physical law, not about
  whether the event is real.

## 10. The time windows of the area-law test

The windows are **set by the paper's analyses**, not by `gwtc_analysis`: the release holds samples only
for the windows the LVK analysed, and the code has no option to change them. Another value, e.g.
t< = −60 M, would need a new PE run on the strain data.

| Window | Values in the release | Paper's reference (constant in `area_law.py`) |
|---|---|---|
| **t<** (inspiral) | −250, −200, −150, −100, −75, −50, −45, −40, −35, −30, −25, −20, −15, −10, −5, 0 M | −40 M (`REFERENCE_CUT`) |
| **t>** (ringdown code, 220 and 220 + 221) | 6, 7, 8, 9, 9.5, 10, 10.5, 11, 12, …, 20 M_f | 10.5 M_f for 220 (`REFERENCE_START`), 6 M_f for 220 + 221 (`OVERTONE_START`) |
| **pyRing** | only 220 at 10.5 M_f | pooled with the ringdown code at 10.5 M_f |

The scans already compute every released combination: every t< against the reference A_f, and every t>
and model against the reference A_i.

### Inspiral truncation: t< from −250 M to 0 M

- **Upper end, 0 M (the peak):** the latest truncation the test allows. After the peak, the data would
  include the ringdown, used to measure A_f, and the two measurements would no longer be independent.
- **Lower end, −250 M (≈ 88 ms before the peak):** truncating earlier throws away the loud last cycles.
  Most of the signal-to-noise ratio is in those final cycles, since the amplitude grows toward the
  merger. With less data, A_i gets wide and the test loses power: only ≈ 3.4σ at −250 M.
- **Spacing:** coarse at the far end (every 50 M), dense near the peak (every 5 M from −50 M). The
  question is how close to the merger the inspiral-only description can be pushed, and the significance
  changes fastest there.

### Ringdown start: t> from 6 M_f to 20 M_f

- **Lower end, 6 M_f:** just after the peak, the signal is not yet a sum of a few quasinormal modes. The
  merger is strongly nonlinear, and many overtones contribute. Even 220 + 221 is only considered reliable
  from about 6 M_f, and 220 alone needs a later start: hence 10.5 M_f for the 220 headline.
- **Upper end, 20 M_f (≈ 7 ms):** the ringdown decays exponentially, with a damping time of about
  12 M_f. Between 10.5 and 20 M_f the amplitude falls by roughly e^(−9.5/12) ≈ 0.45, so A_f gets wide and
  the significance drops to about 1.5σ. Later starts would add little.
- **Spacing:** steps of 1 M_f, plus 9.5 and 10.5 around the reference start.

### The gap between the two windows

With the headline choice, the data between −40 M and +10.5 M_f (about 14 ms before the peak and 3.5 ms
after it) are used by **neither** measurement. That stretch is the merger itself, where neither
description (inspiral or quasinormal modes) is valid, and skipping it makes A_i and A_f truly
independent.

The trade-off shows in the two panels of the scan plot: pushing either window toward the peak adds
signal and raises the significance, but moves into the region where the model (inspiral-only, or a few
quasinormal modes) is less trustworthy.

### What the plot shows

In `area_law_scans_GW250114.png`, the **lines and dots are computed** by `gwtc_analysis` from the
release samples; the **orange markers are the published values** (the `PUBLISHED` dict): the 4.4σ
diamonds, the 3.4σ dashed line, the 3.6σ square. The blue dot at 10.5 M_f in the right panel is slightly
below the diamond because the start-time scan uses the ringdown code alone, while the published 4.4σ
pools it with pyRing.

The exact reasons for the paper's reference choices (−40 M, 10.5 M_f) are given in the paper; the
account above is the general reasoning behind such choices.

## 11. The scalar_cloud mode, explained step by step

This section explains in plain terms what `scalar_cloud` does and how it reproduces the paper of
Roy, Vicente, Aurrekoetxea, Clough & Ferreira (arXiv:2510.17967). The detailed reference is
`docs/modes/scalar-cloud.md`; the code is `scalar_cloud.py` (physics, fits, sampling) and
`scalar_cloud_mode.py` (stages, plots, report).

### 11.1 The question

**Does the signal of GW190728_064510 look as if the two black holes merged inside a cloud of an
ultralight scalar field?** The paper reports tentative evidence for it. The mode refits the public
strain data with two models and compares them:

- **vacuum**: the standard catalog waveform, IMRPhenomXPHM (two black holes in empty space);
- **cloud**: the same waveform, modified by the effect of a scalar cloud on the inspiral
  (IMRPhenomXPHM_Scalar, the paper's model).

Unlike `area_law`, this mode **runs its own inference** on the strain: it does not just read results
from a release. That's why it is expensive and needs a cluster for the full result.

**Which tests come from the paper.** The paper's claim rests on the Bayes factor. `gwtc_analysis`
reproduces it and adds a frequentist check: how often noise alone gives such a gain, and what a real
cloud would give. This attribution is based on the code and `docs/modes/scalar-cloud.md` (which
quote the paper's values), not on a fresh reading of the paper.

| Test | Origin |
|---|---|
| 1. Best-fit gain Λ (11.4) | partly the paper: it reports a median likelihood gain over its posterior and the cloud mode; Λ, the best fits and the seeded fits are `gwtc_analysis`'s |
| 2. Background (11.5) | added by `gwtc_analysis` |
| 3. Sensitivity injections (11.6) | added by `gwtc_analysis` |
| 4. Bayes factor (11.7) | the paper: its actual claim |
| Band analysis (11.9) | added by `gwtc_analysis` |

### 11.2 The physics in four steps

**1. An ultralight boson.** Some extensions of the Standard Model predict very light particles
(axion-like), with masses around 10⁻¹³–10⁻¹¹ eV. Such a light particle has a quantum wavelength of
kilometres, comparable to the size of a stellar-mass black hole. The ratio of the two is the key
parameter:

  α = G M m_φ / (ħ c)  (black-hole size / boson wavelength)

For M = 20 M_☉, α = 0.1 corresponds to m_φ ≈ 6.5 × 10⁻¹³ eV.

**2. A cloud around the binary.** The field forms a bound structure, like an atom's electron cloud
around the nucleus ("gravitational atom"), of radius

  r_φ = GM / (c² α²)

The paper and the code write it r_φ = M/α², in units where G = c = 1: there M stands for the length
GM/c² (section 8), and α is dimensionless, so r_φ is a length. In SI the formula is the **Bohr radius of
a gravitational atom**: replacing α by its definition,

  r_φ = ħ² / (G M m_φ²) = λ̄_φ / α, with λ̄_φ = ħ/(m_φ c) the boson's reduced Compton wavelength,

the same as the hydrogen Bohr radius a₀ = ħ²/(m_e k e²) with the electric coupling k e² replaced by the
gravitational one G M m_φ. For α = 0.19 and M ≈ 20 M_☉: GM/c² ≈ 29.5 km, so
r_φ ≈ 29.5 km / 0.036 ≈ 28 GM/c² ≈ 820 km. The cloud is described by its mean density ρ̄_φ (g/cm³). The best
fit for GW190728 has ρ̄_φ ≈ 6 × 10⁶ g/cm³: very dense, but the cloud is so small that its total mass is
only about 10⁻⁴ of the binary's mass.

**The notation ρ̄_φ.** ρ is the density, the subscript φ says it's the density of the scalar field φ
(not of ordinary matter), and the **bar means "mean"**: ρ̄_φ is the mean density of the field inside the
radius r_φ. A bound cloud's density varies with distance, like the electron density of an atom, peaking
near the centre and falling off outside r_φ. The model doesn't use the profile in detail: it describes
the cloud by this single number, and the torque formula is written in terms of it. The docs use the same
notation (`\bar\rho_\phi`) as the paper, which defines it after its Eq. 2 as "the average mass density of
scalars in the region r ≲ r_φ" (checked against the PDF). In the code the parameter has no
bar: `rho_phi` in the sampling parameters, `rho_cgs` / `rho` in the functions, in g/cm³.

"Mean" is meant loosely. The code computes the cloud mass as ρ̄_φ r_φ³ (`cloud_mass_fraction()` in
`scalar_cloud.py`), without the 4π/3 of a sphere's volume, so ρ̄_φ works as a density scale for the cloud
rather than an exact volume average. That doesn't change the fit, since ρ̄_φ is just the amplitude of
the torque. It matters when comparing the cloud mass, or the 0.1 M cap of the agnostic prior, with other
papers. A factor 4π/3 ≈ 4.2 would also be the size of the "2 to 4 times lower density" difference with the
paper (section 11.8); how the paper normalizes it hasn't been checked, so this is only a possibility.

How would such a cloud form? A **spinning** black hole can amplify waves of a bosonic field and transfer
its spin to them (**superradiance**), growing a cloud spontaneously when α is around 0.1. This motivates
the "superradiance" prior below.

**3. The binary loses energy to the cloud.** As the two black holes orbit inside the cloud, they stir
the field and kick part of it out (**ionization**, like ionizing an atom). This removes orbital angular
momentum, on top of what gravitational waves remove. The ratio of the two losses is δ:

  δ = (cloud loss) / (GW loss), which grows as about v^(−13.5)

v is the orbital velocity (v = (πGMf/c³)^(1/3) for GW frequency f). Since δ grows hugely when v is small,
**the effect is strongest at low frequency**, early in the observed inspiral, when the orbit is wide and
slow. Near the merger the GW loss dominates and the cloud doesn't matter.

**4. The signal gets a phase shift.** With the extra loss, the binary spirals in faster: it spends less
time at each frequency. The waveform's phase then differs from vacuum by ΔΨ(f), the **dephasing**. The
model fixes ΔΨ = 0 at the end of the inspiral (M f_c = 0.018), so the merger and ringdown are unchanged,
and the difference grows toward low frequencies. Each mode (ℓ, m) of the waveform is multiplied by
e^{i ΔΨ_ℓm}. With ρ̄_φ = 0, the cloud waveform is exactly the vacuum one (tested in the code).

Example: ρ̄_φ = 10⁶ g/cm³ and α = 0.1 dephase GW190728 by about 1.6 rad at 20 Hz and 0.13 rad at 30 Hz.
The best fit has 38 rad at 20 Hz.

So the cloud model has **two extra parameters**, ρ̄_φ and α, on top of the ~15 binary parameters.

### 11.3 Background notions used below

- **Likelihood L**: how well a waveform matches the data, given the noise. The code uses the
  **log-likelihood ratio** ln L (signal against noise only); ln L ≈ SNR²/2 for the best template. For
  GW190728, ln L ≈ 92, i.e. SNR ≈ 13.5.
- **PSD** (power spectral density): the noise level of each detector against frequency, taken from the
  PE file. It weights each frequency in the likelihood.
- **Best fit** (maximum likelihood): the parameters that maximize L. Found with a local optimizer,
  Nelder–Mead (scipy).
- **Nested sampling** (bilby + dynesty): explores the whole prior to give the posterior samples **and
  the evidence** Z = ∫ L × prior. "Live points" set its resolution; more is slower and more accurate.
- **Chirp mass** 𝓜: the mass combination that controls the inspiral speed; the best-measured mass.

### 11.4 Test 1: does the cloud fit better? (the statistic Λ)

**Origin: partly the paper.** The paper reports the **median** likelihood gain over its posterior samples (Δln L ≈ 3.7, Supplemental Material, Fig. 8 inset) and a cloud posterior at a lower chirp mass (Fig. 8); it gives no best fit. The statistic Λ = 2Δln L between best fits and the seeded best-fit procedure are `gwtc_analysis`'s.

  Λ = 2 (ln L_cloud − ln L_vacuum) at the two best fits

**Λ is always ≥ 0**: the cloud model contains the vacuum one (ρ̄_φ = 0), so it can always do at least as
well. Two extra parameters will always absorb some noise. So a positive Λ means nothing by itself: you
need to know **what Λ noise alone gives** (Test 2). The usual χ² rule doesn't apply either, because
ρ̄_φ = 0 sits at the edge of its range and α does nothing when ρ̄_φ = 0.

**The difficulty: finding the cloud best fit.** The likelihood has two separate peaks:

- the vacuum peak;
- a cloud peak. The cloud posterior sits at a **lower chirp mass**: 9.69 against 10.14 M_☉ in the paper
  (Fig. 8), 9.89 against 10.14 M_☉ in the mode's samplings (detector frame). A lower chirp mass makes the
  inspiral slower, the cloud makes it faster, so the two compensate and the waveform still matches the
  data near the merger while changing at low frequency. The mode's cloud *best fit* (10.11 M_☉) is not
  in the bulk of that posterior: a fit finds the highest point, the posterior weighs the volume.

An optimizer started at the vacuum best fit stays on the vacuum peak. The first version of the mode did
that (the `profile` step: a grid of α and of dephasings around the vacuum fit) and found only Λ = 0.96:
it missed the cloud peak. The current version takes its starting points from the **nested sampling with
the superradiance prior**, which explores the whole space and does find the cloud peak:

1. take the 25 highest-likelihood samples of that sampling;
2. for each, adjust the merger time (the sampling marginalizes over it, the fits need a fixed one) with
   a fine scan;
3. keep the 4 best as **seeds**, and run two Nelder–Mead passes from each with all the parameters free.

This is the `seeded` stage. It gives **Λ = 8.98** for GW190728:

| | ln L ratio | |
|---|---|---|
| vacuum best fit | 91.81 | 𝓜 = 10.07 M_☉ (detector frame), q = 0.55 |
| cloud best fit | 96.30 | 𝓜 = 10.11 M_☉, ρ̄_φ = 5.8 × 10⁶ g/cm³, α = 0.19 (m_φ ≈ 1.1 × 10⁻¹² eV) |

### 11.5 Test 2: how often does noise give such a Λ? (the background)

**Origin: added by gwtc_analysis.** Not in the paper: no published value to compare with.

1. Take 150 stretches of real detector noise, **off-source**: within an hour of the event, at least
   64 s away from it, with the same length.
2. In each, **inject** the vacuum best fit: a pure vacuum signal, by construction without a cloud.
3. Run **exactly the same procedure** (vacuum fit, seeded cloud fits) and compute Λ.
4. The p-value is (k + 1)/(n + 1), with k of the n segments giving Λ ≥ 8.98.

Result: 2 of 151 (the event and one segment with a loud noise transient, Λ = 112), so **p ≈ 0.013**
(0.007 without that segment). 147 of 150 segments give Λ < 4.1. So the event's Λ is unusual for pure
noise around a vacuum signal, at about the 1% level.

**Caveat:** the seeds come from the event's own posterior. The background therefore asks how often
noise produces such an improvement **in that region of parameter space**, not anywhere. The p-value is
in that sense optimistic.

### 11.6 Test 3: what would a real cloud give? (sensitivity)

**Origin: added by gwtc_analysis.** Not in the paper: no published value to compare with.

The reverse check: inject the **cloud** best fit into the same 150 noise segments and rerun.

Result: Λ ≈ 25 (5–95%: 9.9–44), with ρ̄_φ and α recovered without bias. Only 4.7% of these injections
give Λ ≤ 8.98.

So the event sits **between two distributions**: above what noise gives around vacuum (p ≈ 0.01), and
below what this cloud would give (≈ 5%). A weaker cloud, or a noise fluctuation on a vacuum signal,
fits both better than the best-fit cloud.

### 11.7 Test 4: the Bayes factor (the paper's actual claim)

**Origin: the paper.** Its Bayes factors, priors (SM Eqs. 9–10) and sampler settings.

Nested sampling computes the evidence Z of each model, and the Bayes factor B = Z_cloud / Z_vacuum.
Unlike Λ, the evidence **penalizes extra parameters**: it averages the likelihood over the prior, so a
model that fits better only in a tiny corner of a wide prior is penalized (Occam's razor). That's why the
Bayes factor depends on the prior on (ρ̄_φ, α). Two priors are used:

- **agnostic**: α uniform, ρ̄_φ uniform up to 10⁹ g/cm³. The paper doesn't give the upper end, and the
  result depends on it strongly.
- **superradiance**: only clouds that could have grown by superradiance before the merger (taking
  τ_d = 10⁶ yr between black-hole formation and merger) and not have decayed, with ρ̄_φ below the
  saturation density. This physical prior is much narrower, so it penalizes the cloud less.

| | here | paper |
|---|---|---|
| ln B, superradiance | **1.34 ± 0.30** | 3.5 |
| ln B, agnostic (ρ̄_φ ≤ 10⁹) | **−0.82 ± 0.31** | 0.4 |

ln B > 0 favours the cloud; ln B ≈ 1 is weak, ln B ≈ 3 moderate. The agnostic value rises as the upper
end of the prior is lowered: −0.82 at 10⁹, 0.34 at 3 × 10⁸, 1.90 at 5 × 10⁷ g/cm³. **The paper's 0.4
matches an upper end near 3 × 10⁸ g/cm³.**

### 11.8 What is reproduced, and what is not

**Reproduced:** a cloud mode with ρ̄_φ ~ 10⁷ g/cm³ and α of a few 0.1; the agnostic Bayes factor, for a
prior up to about 3 × 10⁸ g/cm³; the vacuum chirp-mass posterior exactly (10.14 M_☉, detector frame).

The likelihood gain, compared like for like (difference of the posterior medians of the time-marginalized
ln L of the samplings, read from the CC-IN2P3 work directory):

| | gwtc_analysis | paper |
|---|---|---|
| median Δln L, superradiance | 4.20 | ≈ 3.7 (SM, Fig. 8 inset) |
| median Δln L, agnostic | 4.05 | |
| Δln L between best fits (Λ/2) | 4.49 | not given |

The paper reports no best fit, no Λ, no background and no p-value.

Chirp-mass posteriors (detector frame, median and 90%):

| | gwtc_analysis | paper (Fig. 8) |
|---|---|---|
| vacuum | 10.14 (10.05–10.24) | 10.14 (+0.10/−0.09) |
| cloud, superradiance | 9.89 (9.73–10.18) | 9.69 (+0.27/−0.28) |
| cloud, agnostic | 9.88 (9.72–10.12) | |

The identical vacuum value places the paper's numbers in the detector frame. The shift is reproduced,
smaller: 2.5% here against 4.4% in the paper.

**Not reproduced:** a smaller superradiance Bayes factor (1.3 against 3.5), a density 2 to 4 times
lower, and only part of the chirp-mass shift (2.5% against 4.4%). The three go together: a weaker cloud
needs less chirp-mass compensation. Possible causes:

- the frame where the correction is applied: here each mode in the **J frame** (fixed to the total
  angular momentum), in the paper the **co-precessing frame** (following the orbital plane). They're
  identical for aligned spins, and differ when the spins precess;
- the likelihood: this one gives ln L = 89.6 at the catalog's best sample against 96.8 in the PE file
  (different PSD, window and calibration handling; it affects both models equally);
- the cloud samplings were run with 500 live points instead of the paper's 1000;
- the normalisation of the torque, which sets which density gives a given dephasing.

### 11.9 Where the improvement comes from (a warning sign)

**Origin: added by gwtc_analysis.** A diagnostic of the mode, not in the paper.

If a cloud were really there, the gain in ln L should **accumulate over the low-frequency inspiral**,
where the dephasing grows. The code splits the gain by detector and frequency band:

| | 20–30 Hz | 30–35 Hz | 35–60 Hz | > 60 Hz | total |
|---|---|---|---|---|---|
| all detectors | +1.45 | **+5.17** | −0.57 | −1.56 | +4.49 |

The gain sits almost entirely **in a few bins between 29.8 and 31.9 Hz, mostly in the L1 detector**. The
PSDs have no line there. A noise-free cloud signal would spread the gain more evenly, and give some in H1
and below 30 Hz too. So whether the data around 30 Hz hold a noise fluctuation or a real effect is the
open question.

### 11.10 The stages of the mode

| Stage | What it does | Cost (GW190728) |
|---|---|---|
| `prepare` | download the PE file and the strain (event + 150 off-source segments), store the settings, PSDs, best samples | once |
| `fit` | vacuum best fit, then the first cloud fits (`profile` grid + `polish`): the old Λ = 0.96 | ~2 CPU-h |
| `background` | Test 2 with those first fits | ~2 CPU-h per segment |
| `sample` | nested sampling: vacuum, agnostic, superradiance → Bayes factors and the seeds | ~1 day on 64 CPUs each (cluster) |
| `seeded` | cloud fits from the 4 seeds, on the event and the 150 segments → Λ = 8.98 and its p-value | ~1 CPU-h per segment |
| `report` | HTML report, TSV, plots | seconds |

The default run (`prepare fit report`) works on a laptop but gives only the first, incomplete result.
The headline needs `sample`, then `seeded`, on a Slurm cluster. The sensitivity injections (Test 3) are
run by hand: `python -m gwtc_analysis.scalar_cloud seeded --inject cloud`.

### 11.11 The bottom line

- The cloud model fits GW190728 better than vacuum (Λ ≈ 9), more than noise usually allows
  (p ≈ 0.01, optimistic).
  This frequentist check is the mode's own: the paper has no background or p-value.
- But the Bayes factor is weak (ln B ≈ 1.3 with the physical prior, ≈ −0.8 to +1.9 with the agnostic one,
  depending on its upper end).
- The gain comes from a narrow band around 30 Hz in one detector, not from the whole inspiral.
- And Λ is a screening statistic: eccentricity, waveform errors, calibration or a noise feature at low
  frequency could also make vacuum fit the inspiral worse. The claim rests on the Bayes factor under a
  physical prior.

Same boson elsewhere: see 11.12 for the continuous-wave link and the other constraints on the boson,
including the tension with the spin of GW190517 that the paper itself discusses.

### 11.12 Superradiance

**Superradiance is not a counterpart.** A *counterpart* (`counterpart` and `h0_bright` modes) is an
electromagnetic signal from the same event, light seen by telescopes, as for GW170817 (kilonova,
gamma-ray burst, afterglow), which identified the host galaxy and gave the redshift. *Superradiance* is a
physical process around a spinning black hole. A scalar cloud is made of a field that doesn't interact
with light, so it is dark, and GW190728 is a binary black hole merger, not expected to emit light anyway.
The only other signal a cloud gives is gravitational: the continuous wave at 2 m_φ c²/h (section 11.12.3).

**The process.** A wave of a bosonic field passing near a rotating black hole can come out amplified,
taking energy and spin from it. This happens when the wave rotates more slowly than the horizon,
ω < m Ω_H (Ω_H the horizon's angular velocity, m the wave's azimuthal number). For an ultralight boson
with α ≈ 0.1 the waves stay bound to the black hole, so the amplification repeats and a cloud grows
spontaneously, like a laser with the black hole's rotation as the energy source.

**It slows the black hole's spin.** The cloud grows by taking the black hole's rotational energy and
angular momentum. As the black hole slows, Ω_H falls, and the growth stops when Ω_H ≈ ω/m: the black hole
keeps a reduced **critical spin**, set by α, and the cloud reaches its **saturation mass**, up to about 10%
of the black hole's mass.

| Phase | Time scale | What happens |
|---|---|---|
| growth | τ_SR, short (from minutes to ~Myr depending on α and the mass) | the cloud grows exponentially and the black hole spins down |
| saturation | at τ_SR | the spin reaches its critical value, the cloud its maximum mass |
| decay | τ_ann, long | the bosons annihilate into gravitons (the continuous wave), and the cloud slowly disappears |

In `scalar_cloud`, superradiance enters only through the **superradiance prior** (SM Eqs. 9–10 of the
paper). It keeps only clouds with **τ_SR ≤ τ_d ≤ τ_ann**, τ_d being the delay between the black hole's
formation and the merger: the cloud had time to grow, and not yet to decay. That gives the window on α,
0.02 x^(1/9) ≤ α/(1+q) ≤ 0.11 x^(1/15); the saturation mass gives the upper bound on the density.

#### 11.12.1 Constraints on the boson outside scalar_cloud

`scalar_cloud` looks for the effect of a cloud on one binary's inspiral. Other analyses look for the
**consequences superradiance would have** wherever black holes exist. The references and numbers below
are from memory: check them before quoting them.

**1. Black hole spins.** If a boson of mass m_φ exists, every black hole with α around 0.1 spins down
quickly. So a black hole measured to spin fast, in the right mass range and old enough for the cloud to
have grown, **excludes that boson mass**.

- **X-ray binaries**: stellar black holes accreting from a companion; their spins come from X-ray
  spectra (inner edge of the disk, shape of the iron line). Cygnus X-1, for example, spins above 0.95.
  Arvanitaki, Baryakhtar & Huang 2015 excluded scalars of about **6 × 10⁻¹³ to 2 × 10⁻¹¹ eV**.
- **GW events**: the spins measured in binary black hole mergers. Ng et al. 2021 used the GWTC-2 spins
  to exclude about 1.3–2.7 × 10⁻¹³ eV (from memory). The paper points to **GW190517**, whose spin
  measurement disfavours boson masses below about **1.1 × 10⁻¹² eV** (see 11.12.2).
- **Supermassive black holes** (e.g. M87*, imaged by the EHT) constrain much lighter bosons, around
  10⁻²¹–10⁻²⁰ eV.

**2. Continuous gravitational waves from a single cloud.** A saturated cloud emits a nearly
monochromatic wave at f = 2 m_φ c²/h for thousands to millions of years. All-sky searches look for a
nearby unknown black hole with a cloud (the LVK O3 boson-cloud search, 20–610 Hz, shown by the
`continuous_waves` mode); directed searches look toward known black holes such as Cygnus X-1. Nothing
has been detected: the upper limits become excluded boson masses only **for an assumed black hole**
(distance, mass, initial spin, age).

**3. Stochastic background from all the clouds.** The clouds around all the black holes of the Universe
would add up to a background of gravitational waves. LVK upper limits on the isotropic background
exclude ranges around **10⁻¹³ eV** (Tsukada et al. 2019 for scalars, 2021 for vectors), again under
assumptions about the black-hole population. That is the topic of the `stochastic` mode, whose docs
don't treat bosons.

#### 11.12.2 Where the GW190728 boson falls

`scalar_cloud` gives m_φ ≈ 0.7–1.2 × 10⁻¹² eV for GW190728:

| Constraint | Excluded range | GW190728 boson |
|---|---|---|
| continuous waves (O3 all-sky) | depends on the assumed black hole | not excluded: no young cloud around a light, fast-spinning black hole within a few kpc; compatible if such systems are rare or old (`docs/modes/continuous-waves.md`) |
| stochastic background | around 10⁻¹³ eV | below the band: not excluded (not discussed by the paper) |
| X-ray binary spins | ~6 × 10⁻¹³ to 2 × 10⁻¹¹ eV | inside, but the paper argues these spins are model-dependent |
| **GW190517's spin** | below ~1.1 × 10⁻¹² eV (as the paper states it) | **most of the band: the paper's main tension** |

**What the paper says** (Discussion section of arXiv:2510.17967; quotes checked against the PDF):

- "The high spins of BHs in X-ray binaries, inferred from continuum-fitting and reflection methods, have
  been interpreted as excluding scalars in this mass range" (citing Arvanitaki et al. 2015, Cardoso et
  al. 2018, Hoof et al. 2024, Witte & Mummery 2025);
- "These spin estimates, however, are highly model-dependent" (Zdziarski et al. 2024, *ApJ* 962, 101,
  "Black hole spin measurements in LMC X-1 and Cyg X-1 are highly model dependent", and Zdziarski et al.
  2024, *ApJL* 967, L9, "What is the black hole spin in Cyg X-1?"): "accounting for a warm Comptonization
  layer in Cyg X-1 can reduce the inferred spin parameter to a/M ≲ 0.1" (the second paper). A footnote
  adds that such systematics may also reconcile the high X-ray spins with the low spins of most GW
  sources;
- "Constraints from few well-measured BH spins **through GWs** [Aswathi et al. 2025, Ng et al. 2021,
  Caputo et al. 2025, on GW231123] do not exclude the full range favored by our analysis of GW190728";
- "the main tension arises from GW190517, whose spin measurement disfavors m_φ ≲ 1.1 × 10⁻¹² eV";
- "Notably, joint spin measurements from the full GWTC-2 [Ng et al. 2021] showed weak evidence for a light
  scalar in the mass range suggested by our analysis";
- the non-detection of continuous waves is cited as a constraint, without a comparison with GW190728;
  the stochastic background is not discussed (Yuan et al. 2022 is only in the references).

So the paper itself identifies the tension: the GW190728 band (0.70–1.17 × 10⁻¹² eV here, "about
10⁻¹² eV" in the paper) lies mostly below the ~1.1 × 10⁻¹² eV that GW190517 disfavours, leaving only its
upper end. The X-ray exclusions it considers weaker, because the spins are model-dependent. Other ways
out, not argued by the paper: the black hole must be old enough for the cloud to have grown, and a boson
with **self-interactions** (as an axion with a small decay constant) can stop the growth early, leaving
the spin high.

#### 11.12.3 The continuous-wave link

The posterior boson mass, m_φ ≈ 0.89 (0.70–1.17) × 10⁻¹² eV with the superradiance prior, would make an
isolated spinning black hole's cloud emit at 2 m_φ c²/h ≈ **340–570 Hz**, inside the band of the LVK
all-sky search for boson clouds. The `continuous_waves` mode draws its exclusions with this band.

## 12. Sources of the formulas and numbers

For each formula or number used in this file: where it comes from. "Code" means the formula is
implemented in `gwtc_analysis`, with the file. References marked **(from memory)** were cited from my
knowledge rather than checked against the paper itself, so check them before quoting them.

### Sections 1 (viewing angle) and 2 (recession velocity)

| Formula / number | Source |
|---|---|
| viewing angle = min(θ, 180° − θ), `theta_jn` before `iota` | code: `counterpart.py`, `viewing_angle()` |
| `theta_jn` (J to line of sight), `iota` (L to line of sight) | parameter conventions of LALSuite / bilby, used in the GWTC PE files |
| v_H = v_r − ⟨v_p⟩, σ added in quadrature | code: `counterpart.py`, `hubble_flow_redshift()`; method of Abbott et al. 2017, *Nature* 551, 85 (GW170817 standard siren) |
| NGC 4993: v_r = 3327 ± 72 km/s, v_p = 310 ± 150 km/s | registered in `counterpart.py`; values of Abbott et al. 2017 (above) |

### Section 3 (CMB frame)

| Formula / number | Source |
|---|---|
| solar dipole 369.82 ± 0.11 km/s toward (l, b) = (264.021°, 48.253°) | Planck Collaboration 2018 results I, *A&A* 641, A1 (2020) (from memory) |
| T(θ) ≈ T₀ (1 + (v/c) cos θ), dipole ≈ 3.36 mK | first-order Doppler effect; classic derivation in Peebles & Wilkinson 1968, *Phys. Rev.* 174, 2168 (from memory) |
| v_CMB ≈ v_helio + v_☉ cos(angle to apex) | first-order projection of the solar velocity; NED's velocity calculator applies the same correction |
| NGC 4993 v_CMB ≈ 3241 km/s | computed in this session with astropy, assuming v_helio = 2930 km/s (heliocentric value from memory; check NED) |
| T₀ = 2.725 K | Fixsen 2009, *ApJ* 707, 916 (COBE/FIRAS) (from memory) |

### Sections 4 (Planck) and 5 (chain of motions)

| Formula / number | Source |
|---|---|
| orbital dipole as Planck's absolute calibrator | Planck 2015 and 2018 results, HFI and LFI calibration papers (from memory) |
| aberration and Doppler boosting, v ≈ 384 ± 78 km/s | Planck Collaboration 2013 results XXVII, *A&A* 571, A27 (2014) (from memory) |
| quasar dipole about twice the CMB one | Secrest et al. 2021, *ApJL* 908, L51 (from memory) |
| Local Group ≈ 620 km/s toward (l, b) ≈ (272°, 30°) | Kogut et al. 1993, *ApJ* 419, 1 (from memory) |
| Sun relative to the LSR, ~20 km/s | Schönrich, Binney & Dehnen 2010, *MNRAS* 403, 1829 (from memory) |
| Sun around the Galactic centre, ~230–250 km/s | range of recent measurements (e.g. Gaia, GRAVITY); not a single source |
| Laniakea, Great Attractor | Tully et al. 2014, *Nature* 513, 71 (from memory) |

### Section 7 (area law) and section 10 (windows)

| Formula / number | Source |
|---|---|
| A = 8π (GM/c²)² [1 + √(1 − χ²)] | horizon area of the Kerr metric (Kerr 1963); code: `area_law.py`, `kerr_area()` |
| area theorem A_f ≥ A_1 + A_2 | Hawking 1971, *Phys. Rev. Lett.* 26, 1344 |
| S = (⟨A_f⟩ − ⟨A_i⟩)/√(σ_f² + σ_i²), pooling, seed 250114 | LVK, GW250114 discovery paper, arXiv:2509.08054, and its release (Zenodo 16877102); code: `area_law.py` |
| P(A_f < A_i) over random pairs, ratio over 200,000 pairs | code: `area_law.py` (choices of `gwtc_analysis`, not the paper's) |
| published 4.4σ, 3.4σ, −10 M, 3.6σ; the windows −250…0 M, 6…20 M_f | arXiv:2509.08054; the windows as listed in the release files |
| NRSur7dq4 | Varma et al. 2019, *Phys. Rev. Research* 1, 033015 (from memory) |
| `ringdown` code | Isi & Farr 2021, arXiv:2107.05609 (from memory) |
| pyRing | Carullo, Del Pozzo & Veitch 2019, *Phys. Rev. D* 99, 123029 (from memory) |
| GW150914 area-law test | Isi et al. 2021, *Phys. Rev. Lett.* 127, 011103 |
| GW250114 masses (~66 M_☉ source, z ≈ 0.09) | from memory of arXiv:2509.08054; can be computed from the release samples |

### Section 8 (M unit)

| Formula / number | Source |
|---|---|
| G = 6.674 × 10⁻¹¹ m³ kg⁻¹ s⁻², c = 299,792,458 m/s | CODATA 2018 (c is exact by definition of the metre) |
| GM_☉/c² = 1.4766 km, GM_☉/c³ = 4.925 μs | IAU 2015 Resolution B3, nominal GM_☉ = 1.3271244 × 10²⁰ m³ s⁻²; code: `MSUN_KM` in `area_law.py` |
| Schwarzschild radius 2GM/c² | Schwarzschild 1916 |
| ISCO at r = 6 GM/c², f_GW ≈ 4400 Hz (M_☉/M) | Bardeen, Press & Teukolsky 1972, *ApJ* 178, 347 (from memory) |
| 220 mode, Schwarzschild: M ω = 0.3737 − 0.0890 i (period 16.8 M, damping 11.2 M) | Leaver 1985, *Proc. R. Soc. A* 402, 285; tables in Berti, Cardoso & Starinets 2009, *CQG* 26, 163001 (from memory) |
| 220 for χ ≈ 0.7: M_f ω ≈ 0.53, τ ≈ 12 M_f | same tables / fits of Berti, Cardoso & Will 2006, *Phys. Rev. D* 73, 064030 (approximate, from memory) |

### Section 9 (σ)

| Formula / number | Source |
|---|---|
| σ = sample standard deviation | code: `area_law.py`, `gaussian_significance()` (`numpy.std`) |
| tail probabilities (5.4 × 10⁻⁶ at 4.4σ, etc.) | computed in this session with `scipy.stats.norm.sf` |
| 5σ discovery convention | convention of particle physics, also used in GW astronomy |

### Section 11 (scalar cloud)

| Formula / number | Source |
|---|---|
| model, α = GMm_φ/(ħc), r_φ = GM/(c²α²) (M/α² with G = c = 1) | Roy, Vicente, Aurrekoetxea, Clough & Ferreira, arXiv:2510.17967; code: `scalar_cloud.py` |
| torque L̇_φ = 4π ρ̄_φ M³ I/v⁴ and the sum I(v, α, q) | arXiv:2510.17967, Eq. 2 and Eq. 6 of its Supplemental Material (as cited in `docs/modes/scalar-cloud.md`); code: `torque_integral()` |
| GW loss L̇_GW = (32/5) η² M v⁷ | quadrupole formula for a circular orbit, Peters 1964, *Phys. Rev.* 136, B1224 |
| δ = (5π/8) ρ̄_φ M² I/(η² v¹¹), ~ v^(−13.5) | ratio of the two losses above; code: `dephasing_22()` |
| ΔΨ₂₂(f) = 2π ∫ (f′ − f) Δ(dt/df) df′ | stationary-phase approximation of an inspiral (e.g. Cutler & Flanagan 1994, *Phys. Rev. D* 49, 2658, from memory), as used in the paper; code: `dephasing_22()` |
| M f_c = 0.018, end of the inspiral | end of the inspiral region of IMRPhenomX (code comment `MF_INSPIRAL_END`) |
| (ℓ, m) modes: (|m|/2) ΔΨ₂₂(2f/|m|) | the paper; code: `scalar_cloud_binary_black_hole()` |
| ionization of gravitational atoms | Baumann, Bertone, Stout & Tomaselli 2022, *Phys. Rev. Lett.* 128, 221102 (from memory) |
| superradiance | Brito, Cardoso & Pani, *Superradiance*, Lect. Notes Phys. 906 (2015); Arvanitaki & Dubovsky 2011, *Phys. Rev. D* 83, 044026 (from memory) |
| superradiance window 0.02 x^(1/9) ≤ α/(1+q) ≤ 0.11 x^(1/15) | arXiv:2510.17967, Supplemental Material Eq. 9 (code: `superradiance_window()`) |
| saturation density 1.3 × 10⁸ g/cm³ … | arXiv:2510.17967, Supplemental Material Eq. 10 (code: `superradiance_density_max()`) |
| cloud mass ρ̄_φ r_φ³ < 0.1 M | the paper's agnostic prior; code: `cloud_mass_fraction()` |
| Λ = 2(ln L_cloud − ln L_vacuum), background p = (k+1)/(n+1) | choices of `gwtc_analysis` (`scalar_cloud.py`, `scalar_cloud_mode.py`), not the paper's |
| ln L ≈ SNR²/2 at the best fit | property of the matched-filter likelihood for a well-fitting template |
| evidence Z = ∫ L π, nested sampling | Skilling 2006, *Bayesian Analysis* 1, 833; dynesty: Speagle 2020, *MNRAS* 493, 3132; bilby: Ashton et al. 2019, *ApJS* 241, 27 (from memory) |
| sampler settings (1000 live points, naccept 60, dlogz 0.1) | the paper's, as recorded in `SAMPLER_SETTINGS` |
| IMRPhenomXPHM | Pratten et al. 2021, *Phys. Rev. D* 103, 104056 (from memory) |
| GW190728 results (Λ = 8.98, ln B, ρ̄_φ, band table) | runs of this mode, as reported in `docs/modes/scalar-cloud.md` |
| paper's ln B = 3.5 and 0.4, ρ̄_φ ≈ 0.7 × 10⁸ | arXiv:2510.17967, as quoted in `docs/modes/scalar-cloud.md` |
| m_φ ≈ 6.5 × 10⁻¹³ eV for α = 0.1, M = 20 M_☉ | m_φ = α ħc³/(GM); code: `scalar_mass_ev()` |
| r_φ ≈ 800 km, cloud mass ~10⁻⁴ M for the best fit | computed in this session from the formulas above |
| continuous wave at 2 m_φ c²/h | boson annihilation into gravitons (Arvanitaki & Dubovsky 2011, above); LVK O3 all-sky boson search band 20–610 Hz as quoted in the docs |
| superradiance condition ω < m Ω_H, spin-down, saturation | Zel'dovich 1971; Brito, Cardoso & Pani, *Superradiance*, Lect. Notes Phys. 906 (2015) (from memory) |
| τ_SR ≤ τ_d ≤ τ_ann (superradiance prior) | arXiv:2510.17967, SM Eqs. 9–10 |
| X-ray binary spins exclude ~6 × 10⁻¹³–2 × 10⁻¹¹ eV | Arvanitaki, Baryakhtar & Huang 2015, *Phys. Rev. D* 91, 084011 (from memory) |
| GW-event spins exclude ~1.3–2.7 × 10⁻¹³ eV | Ng, Vitale, Hannuksela & Li 2021, *Phys. Rev. Lett.* 126, 151102 (from memory) |
| GW190517 disfavours m_φ ≲ 1.1 × 10⁻¹² eV; X-ray spins model-dependent; GWTC-2 spins weak evidence in the range | arXiv:2510.17967, Discussion (checked against the PDF) |
| Cyg X-1 / LMC X-1 spins model-dependent; Cyg X-1 down to a/M ≲ 0.1 | Zdziarski et al. 2024, *ApJ* 962, 101; Zdziarski et al. 2024, *ApJL* 967, L9 (references [153], [154] of the paper) |
| stochastic-background exclusions around 10⁻¹³ eV | Tsukada et al. 2019, *Phys. Rev. D* 99, 103015; Tsukada et al. 2021, *Phys. Rev. D* 103, 083005 (from memory) |
| O3 all-sky boson-cloud search, 20–610 Hz | LVK, *Phys. Rev. D* 105, 102001 (2022) (from memory); `docs/modes/continuous-waves.md` [107] |

`docs/references.md` holds the project's own reference list: [79] GW250114, [80] Hawking 1971,
[81] Isi et al. 2021, [100] Roy et al.
