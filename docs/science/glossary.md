# Glossary

The terms and the physical parameters used in these pages, in two parts: [Terms](#terms) (data, methods, software,
astrophysical objects) and [Physical parameters](#physical-parameters) (the quantities and their symbols). Each part
is grouped by theme.

## Terms

### Detectors, catalogs and data

| Term | Definition |
|---|---|
| **GW** | Gravitational wave |
| **LVK** | LIGO–Virgo–KAGRA, the collaboration of the ground-based GW detectors (H1 Hanford, L1 Livingston, V1 Virgo, K1 KAGRA) |
| **Strain** | The relative change of length a passing gravitational wave produces between the test masses of a detector, \(h = \Delta L / L\), of order 10⁻²¹: the time series the detectors record and that parameter estimation analyses |
| **GWTC** | Gravitational-Wave Transient Catalog: GWTC-1, -2.1, -3, -4.0, -5.0, and the update GWTC-4.1 |
| **GWOSC** | Gravitational Wave Open Science Center, the public GW data portal: strain, event lists and their parameters |
| **Confident / marginal lists** | GWOSC event lists. A confident list holds the candidates with p_astro ≥ 0.5; a marginal list (e.g. `GWTC-3-marginal`) less significant candidates, some of which the LVK population and cosmology analyses still use |
| **O1, O2, O3a, O3b, O4a, O4b** | Observing runs, 2015–2025; O3 and O4 are split in halves |
| **ER15** | Engineering run 15, just before O4a |
| **GPS time** | Seconds since 6 January 1980, the time stamp of GW data and events (GW170817: 1187008882.43) |
| **Superevent, GraceDB** | GraceDB is the LVK database of low-latency candidates and their alerts; a superevent is a candidate there, named like S231113bw, before it gets its catalog name GWyymmdd_hhmmss |
| **GCN** | General Coordinates Network, the NASA system that distributes transient alerts and circulars |
| **Zenodo record, version, concept ID** | Zenodo is the CERN open-data repository where the LVK publishes PE files, sky maps and injections. A release is a record with numbered versions (v1, v2…); the concept ID names all the versions together, and the package resolves the latest one at run time ([Data sources](../data-sources.md#zenodo-release-versions)) |
| **LIGO DCC** | The LIGO Document Control Center (dcc.ligo.org), whose public entries (e.g. LIGO-P1800370) hold papers and data releases: the GWTC-1 GW170817 products, the continuous-wave upper limits |
| **S3 mirror, Galaxy history** | Two mirrors of the public files the package downloads: the `gwtc` bucket of a MinIO server (S3 is the object-storage protocol, MinIO a server implementing it), and the public "GWTC" history on usegalaxy.org. Both fall back on the original source ([Data sources](../data-sources.md)) |
| **HEASARC TAP, VizieR** | HEASARC is NASA's high-energy archive, queried through TAP (Table Access Protocol) for the burst and neutrino catalogs; VizieR is the catalog service of the CDS (Strasbourg), queried through its ASU-TSV interface for GLADE+ |

### Sources and astrophysics

| Term | Definition |
|---|---|
| **BBH / BNS / NSBH** | Binary black hole / binary neutron star / neutron star–black hole binary (written NS-NS and BH-NS by catalog_statistics) |
| **Kerr black hole** | A rotating black hole of general relativity (Kerr 1963), the form an isolated, uncharged black hole settles into. Two numbers describe it completely, its mass M and its dimensionless spin χ = cJ/(GM²) from 0 (not rotating, the Schwarzschild black hole) to 1 (maximal): the no-hair theorem. Its horizon area is A = 8π (GM/c²)² (1 + √(1 − χ²)), used by [area_law](../modes/area-law.md); merger remnants typically have χ ≈ 0.7 |
| **Primordial black hole (PBH)** | A black hole formed in the first fraction of a second after the Big Bang from the collapse of a dense region, not from a star; any mass is possible, set by the formation time. A dark-matter candidate, because it forms before the light elements and is not counted as baryons; \(f_\text{PBH}\) is the fraction of the dark matter it would make up ([Primordial black holes](primordial-black-holes.md)) |
| **Sub-solar mass (SSM)** | A compact object below ~1 M☉, which stars do not leave behind as neutron stars or black holes: the cleanest signature of a PBH merger. LVK alerts carry HasSSM, the probability that a component is below 1 M☉ |
| **Mass gaps** | Lower mass gap, about 3–5 M☉, between the heaviest neutron stars and the lightest black holes; pair-instability gap, above about 50 M☉, where stars should leave no black hole |
| **Hierarchical merger** | A merger whose component is itself the remnant of an earlier merger, expected in dense star clusters and AGN disks; hinted by a mass in the pair-instability gap, a spin near 0.7 or a negative χ_eff ([event_selection](../modes/event-selection.md)) |
| **Formation channels** | The ways a binary comes into existence, from two separate objects to a bound pair that eventually merges; they lead to the same kind of event but predict different rates, masses, spins and eccentricities. Stellar black holes: isolated binary evolution (spins tend to be aligned with the orbit, χ_eff > 0), the dynamical channel in dense clusters (isotropic spins, as many negative χ_eff as positive, [spin_population](../modes/spin-population.md)), the AGN-disk channel, hierarchical mergers. Primordial black holes: early binaries, formed in the radiation era, and late binaries, captured in PBH clusters today ([Primordial black holes](primordial-black-holes.md#how-pbh-binaries-form-and-merge)); each gives its own limit on f_PBH |
| **Equation of state (EOS)** | The relation between pressure and density of neutron-star matter beyond nuclear density. All neutron stars share it; it sets their radius and tidal deformability at each mass, a "stiff" EOS giving large, easily deformed stars ([Tidal deformability](tidal-deformability.md)) |
| **Electromagnetic counterpart** | Light (or neutrinos) from the source of a GW event, such as a kilonova or a flare: its position identifies the host galaxy, whose redshift the GW signal alone does not give |
| **Kilonova, AT** | A kilonova is the optical and infrared glow, lasting days, of the radioactive decay of heavy elements in the matter ejected by a neutron-star merger. AT (astronomical transient) is the prefix of the IAU Transient Name Server: AT2017gfo is the kilonova of GW170817, in the galaxy NGC 4993 |
| **GRB; Fermi-GBM, Swift** | A gamma-ray burst is a short flash of gamma rays; short GRBs come from neutron-star mergers. GBM is the burst monitor of the Fermi satellite (degree-scale positions); Swift localizes bursts with its BAT, XRT and UVOT instruments, to arcseconds for an afterglow |
| **IceCube Gold / Bronze alerts** | Public real-time alerts of IceCube, the neutrino telescope at the South Pole, for track-like high-energy neutrino events; Gold and Bronze are their two purity classes |
| **AGN, ZTF** | An active galactic nucleus is the bright accretion disk around a galaxy's central black hole; a binary merging inside it could cause a flare. ZTF (Zwicky Transient Facility) is an optical survey; ZTF19abanrhr is the AGN flare proposed as the counterpart of GW190521, a candidate only |
| **Continuous wave (CW)** | A long-lived, nearly monochromatic GW signal, from a spinning, slightly deformed neutron star or from a boson cloud around a black hole; none has been detected ([continuous_waves](../modes/continuous-waves.md)) |
| **Millisecond pulsar, LMC** | A millisecond pulsar is a neutron star spinning hundreds of times per second (J0437−4715, at 160 pc); the LMC (Large Magellanic Cloud) is a satellite galaxy about 50 kpc away, host of J0537−6910 |
| **Stochastic background** | The sum of all the GW signals too weak or too many to be resolved one by one, described by its energy-density spectrum Ω_GW(f) ([stochastic](../modes/stochastic.md)) |
| **Ultralight boson, superradiance** | A hypothetical particle far lighter than any known one (m_φ ~ 10⁻¹² eV here), a dark-matter candidate. By superradiance its field draws rotational energy from a spinning black hole and grows into a bound cloud when its Compton wavelength matches the hole's size (α ~ 0.1); the superradiance prior requires the cloud to have grown before the merger and not to have decayed ([scalar_cloud](../modes/scalar-cloud.md)) |
| **Strong lensing, images** | A galaxy near the line of sight splits a GW into two or four images arriving days to months apart, with different magnifications; a saddle-point (type II) image also gets a phase shift of π/2, the Morse phase ([lens_pairs](../modes/lens-pairs.md)) |
| **Optical depth, magnification bias** | The optical depth is the probability that a source at a given redshift is strongly lensed; the magnification bias is the extra weight of the lensed configurations whose images are bright enough to be detected |
| **Kroupa distribution** | A stellar initial mass function; here the power law m^−2.3 of the simulated Galactic black holes of the O3 boson-cloud search |

### Signal and waveforms

| Term | Definition |
|---|---|
| **IMR** | Inspiral–merger–ringdown, the three phases of a binary signal: the two bodies spiral in (the chirp), their horizons fuse (the peak of the signal), the remnant black hole settles down. An IMR waveform model covers all three in one template, and an IMR analysis, the standard catalog one, fits it to the whole signal; in that analysis the remnant comes from numerical-relativity fits, not from the ringdown alone (see [area_law](../modes/area-law.md)) |
| **Polarizations, antenna pattern** | A gravitational wave has two independent polarizations, plus and cross, h₊(t) and h×(t). Each detector records a combination of the two, weighted by its antenna-pattern factors F₊, F×, which depend on the sky position and the polarization angle ψ |
| **Precession** | When the spins are tilted against the orbital angular momentum L, L turns around the nearly fixed total angular momentum J and the orbital plane wobbles, modulating the amplitude and phase of the signal (the P of IMRPhenomXPHM) |
| **Higher multipoles (HM)** | The harmonics of the signal beyond the dominant (ℓ, m) = (2, 2) quadrupole, e.g. (3, 3), (4, 4), (2, 1): strong for unequal masses and inclined orbits, they help to measure the inclination and the mass ratio (the HM of IMRPhenomXPHM) |
| **Co-precessing frame, J frame** | Frames for the modes of a precessing binary: the J frame is fixed to the total angular momentum, the co-precessing frame follows the orbital plane; they coincide for aligned spins |
| **Post-Newtonian (PN) expansion** | The expansion of the inspiral in powers of the orbital speed v/c; an "n PN" term is suppressed by (v/c)^2n relative to the Newtonian one (spins enter at 1.5PN, tides at 5PN). TaylorF2 and SpinTaylor are PN approximants ([Waveform models](waveforms.md)) |
| **Stationary-phase approximation** | The approximation that gives the frequency-domain phase of an inspiral from its time–frequency evolution dt/df; it turns the time-domain effect of a cloud or of tides into a phase shift |
| **Numerical relativity (NR)** | The numerical solution of Einstein's equations for the last orbits and the merger of two compact objects: too expensive for parameter estimation, it builds, calibrates and tests the waveform models and the remnant fits |
| **Phenomenological model (IMRPhenom)** | Analytic amplitude and phase, mostly in the frequency domain, fitted to numerical-relativity and effective-one-body waveforms: the fastest models, the workhorse of the catalogs |
| **Effective-one-body model (SEOBNR, TEOBResumS)** | Maps the two-body problem onto one body in a deformed black-hole spacetime; its equations of motion are integrated in time and calibrated to numerical relativity: slower than phenomenological models, very accurate |
| **Surrogate model (NRSur7dq4), ROM** | A surrogate interpolates a library of numerical-relativity simulations, the most faithful model within its range (NRSur7dq4: 7 parameters, mass ratios up to 4); a reduced-order model (_ROM) is a fast frequency-domain surrogate of an EOB model |
| **IMRPhenomXPHM** | Frequency-domain waveform model with precession and higher modes (Pratten et al. 2021 [\[51\]](../references.md#ref-51)); the other models and their names: [Waveform models](waveforms.md) |
| **XPHM-SpinTaylor** | IMRPhenomXPHM with numerically evolved spin precession (Colleoni et al. 2024 [\[52\]](../references.md#ref-52)) |
| **NRTidal** | Tidal terms calibrated on numerical relativity, added to a binary-black-hole model for neutron stars (IMRPhenomPv2_NRTidal) |
| **Quasinormal modes** | The damped oscillations of the ringdown, like those of a struck bell; for a Kerr black hole their frequencies and damping times depend only on its mass and spin, so a ringdown-only fit measures the remnant. They are labelled (ℓ, m, n): 220 is the fundamental quadrupole mode, 221 its first overtone, which decays faster and matters only near the peak |
| **Dephasing** | The change of the GW phase (radians) caused by an extra loss of energy or angular momentum compared with vacuum, e.g. by a scalar cloud |
| **PSD, ASD** | Power spectral density, the detector noise spectrum, used to whiten the data; the amplitude spectral density is its square root |
| **Calibration envelope** | Uncertainty of the detector calibration in amplitude and phase, marginalized in PE |
| **Whitening** | Dividing the data by the noise amplitude spectrum, so that all frequencies have equal noise |
| **q-transform** | Time–frequency representation of the strain, showing the chirp |
| **Glitch, gating** | A glitch is a short non-Gaussian noise transient in a detector; gating zeroes the data around a loud one, with a smooth window, before filtering |
| **Matched filter, template bank** | The matched filter correlates the data with a waveform template, weighted by the noise spectrum; a search runs it with a template bank, a grid of waveforms covering the masses and spins |
| **Search pipeline** | A program that scans the strain for signals (matched filtering, e.g. GstLAL, PyCBC) and assigns each candidate a FAR; GWOSC lists each event's lowest FAR over the pipelines |
| **SNR** | Signal-to-noise ratio; the matched-filter SNR compares the data with a template |
| **FAR** | False-alarm rate: how often noise alone produces a candidate at least this significant (per year) |
| **p_astro** | Probability that a candidate is astrophysical |

### Parameter estimation and statistics

| Term | Definition |
|---|---|
| **PE** | Parameter estimation: the Bayesian inference of one event's masses, spins, distance, sky position… |
| **PE label** | One analysis in a PE file (waveform, settings), e.g. `C00:IMRPhenomXPHM-SpinTaylor`; `Mixed` combines several. The prefix C00, C01, C02 is the version of the detector calibration the analysis used |
| **LowSpin / HighSpin** | The two spin priors of neutron-star analyses: dimensionless spins up to 0.05, as observed in Galactic neutron stars, or up to 0.89 |
| **PESummary** | The LVK file format and library for PE results (Hoy & Raymond 2021 [\[62\]](../references.md#ref-62)) |
| **`nocosmo` / `cosmo` files** | The two PE files per event of GWTC-2.1 and GWTC-3: `nocosmo` holds the samples under their sampling prior ∝ D_L²; `cosmo` holds them reweighted to uniform in comoving volume (the posteriors of the papers), while still recording the D_L² prior ([Which PE file is read](../data-sources.md#which-pe-file-is-read)) |
| **PE prior, distance prior** | The prior of the released samples: masses uniform in the detector frame, isotropic spins and orientation, and a distance prior, ∝ D_L² (uniform in Euclidean volume, O1–O3) or uniform in comoving volume and source-frame time (O4). A hierarchical analysis divides it out |
| **Maximum-likelihood (maxL) sample** | The posterior sample of highest likelihood; its waveform, projected on each detector, gives the strain overlay and the matched filter of [parameters_estimation](../modes/parameters-estimation.md) |
| **Intrinsic / extrinsic parameters** | Intrinsic: those of the binary itself (masses, spins, tidal deformabilities). Extrinsic: how it is seen (sky position, distance, inclination, polarization, phase, coalescence time) |
| **Distance–inclination degeneracy** | An inclined binary is fainter than a face-on one at the same distance, so the GW amplitude measures distance and inclination together; an independent viewing angle (e.g. from a radio jet) breaks it |
| **Hierarchical inference** | Inference of population parameters from many events, each with its own uncertain parameters |
| **Hyperprior** | Prior on population (and cosmological) parameters |
| **Selection effects** | The detected events are a biased sample: heavy and nearby binaries are over-represented. Population and H₀ inference corrects for it with the detectable fraction ξ (written β(H₀) for the bright siren), measured with injections ([Selection effects](selection-effects.md)) |
| **Injection** | A simulated signal added to the data (or evaluated semi-analytically) to measure the detection probability |
| **Real / semi-analytic injections** | Real injections are added to the O3 and O4 data and searched by the pipelines; semi-analytic ones (O1, O2) are found when their SNR computed on the noise spectra exceeds a threshold. A mixture file combines both, with mixture weights |
| **Found injection** | An injection that passes the detection criteria (FAR or SNR threshold) |
| **Importance reweighting** | Estimating an expectation under one density from samples drawn from another, weighting each by the ratio of the two (e.g. p_pop/p_draw): one set of injections or of posterior samples then serves any population or setting |
| **n_eff, effective sample size** | Effective number of samples of a weighted Monte Carlo sum, \((\sum x)^2/\sum x^2\); a small one makes the sum unreliable |
| **Jeffreys prior, scale-free likelihood** | Uninformative priors on a rate: with p(R) ∝ R^−1/2 and N detections, R⟨VT⟩ follows a Gamma(N + ½) distribution; with p(R) ∝ 1/R the rate integrates out, leaving a likelihood of the population shape only ([Merger rates](merger-rates.md)) |
| **Nested sampling** | Bayesian sampling that computes the evidence and the posterior together (dynesty) |
| **Live points** | The set of samples nested sampling evolves; more live points, finer exploration |
| **ln Z, Bayes factor** | ln Z is the log Bayesian evidence; the Bayes factor of two models is the ratio of their evidences, ln B = ln Z₁ − ln Z₂, and depends on their priors |
| **MCMC, emcee** | Markov-chain Monte Carlo; emcee is an ensemble MCMC sampler whose walkers explore the parameter space together; the autocorrelation time is the number of steps after which its samples are independent |
| **Seed** | Initialization of the random generator of one sampler run |
| **Kernel density estimate (KDE)** | A smooth density built by placing a kernel (e.g. a Gaussian) on each sample; Silverman's rule sets its width, 1.06 σ n^−1/5 for n samples |
| **MAP, HDI, equal-tailed interval** | The maximum a posteriori is the peak of the posterior; the highest-density interval is the narrowest interval holding a given probability; the equal-tailed interval leaves the same probability (e.g. 5%) on each side |
| **P–P test** | A check of the credible intervals on simulated events: the fraction of true values inside each credible level should equal that level, tested with a Kolmogorov–Smirnov p-value |
| **Off-source background, p-value** | The same analysis run on stretches of noise near the event (with a vacuum signal injected) gives what noise alone produces; the p-value is the fraction of that background at least as extreme as the event, (k + 1)/(n + 1) for k of n segments ([scalar_cloud](../modes/scalar-cloud.md)) |
| **False-positive probability (FPP)** | In [lens_pairs](../modes/lens-pairs.md), the fraction of the pairs of the selection with a statistic at least that of a given pair (itself included): a flagged pair is a candidate, not a detection |

### Sky localization

| Term | Definition |
|---|---|
| **Sky map, HEALPix, nside** | A sky map is the probability of the source position over the sky, stored as a FITS file on HEALPix (Hierarchical Equal Area isoLatitude Pixelization) pixels: 12·nside² pixels of equal area, 0.84 deg² at nside 64, 0.23° wide at nside 256 |
| **Multi-order map (MOC)** | A HEALPix map with fine pixels only where the probability is and coarse ones elsewhere, the format of the LVK sky maps and of the threshold map of galaxy_catalog |
| **3D sky map** | A sky map whose pixels also hold the distribution of the luminosity distance along that line of sight (`DISTMU`, `DISTSIGMA`, `DISTNORM`), i.e. a probability per unit volume |
| **Credible level, credible region** | Sorting the pixels by decreasing probability and summing, the greedy credible level of a pixel is the cumulative sum at that pixel; the 90% credible region is the set of pixels at or below 90%, the smallest region holding 90% of the probability |
| **Credible area (A90), credible volume** | The area (deg²) of a credible region, 90% by default (A50, A95 likewise), and the volume (Mpc³) of the 3D one |
| **Searched probability** | The credible level of the smallest region that contains a given position, e.g. that of a counterpart: small when the position sits near the peak of the localization |
| **RAVEN, overlap I_Ω** | RAVEN is the LVK low-latency search for GRB and neutrino coincidences; `search_skymaps --counterparts` reproduces it offline. I_Ω = N_pix Σ p_GW p_ext is the GW probability covered by a counterpart's localization, relative to a uniform sky ([Multi-messenger coincidences](../modes/search-skymaps.md#multi-messenger-coincidences)) |

### Cosmology, populations and galaxies

| Term | Definition |
|---|---|
| **Standard siren** | A GW source used as a distance indicator |
| **Bright / dark / spectral siren** | Redshift from an electromagnetic counterpart / a galaxy catalog / features of the mass distribution |
| **ΛCDM** | The standard cosmological model: a cosmological constant Λ, cold dark matter, flat space |
| **Planck 2015 cosmology** | Flat ΛCDM with the Planck 2015 parameters, the reference cosmology of the GW catalogs: `Planck15` (astropy: H₀ = 67.74, Ω_m = 0.3075) or `Planck15_LAL` (LAL: H₀ = 67.90, Ω_m = 0.3065) (Planck Collaboration 2016 [\[32\]](../references.md#ref-32)); see [How the redshift is obtained](redshift.md#the-planck-2015-cosmology) |
| **Planck (H₀)** | H₀ = 67.4 ± 0.5 km/s/Mpc, inferred from the cosmic microwave background measured by the Planck satellite, assuming flat ΛCDM (Planck Collaboration 2020 [\[33\]](../references.md#ref-33)). It reflects the early Universe: the size of the sound waves in the primordial plasma, carried forward to today by the model |
| **Distance ladder** | Distances built rung by rung, each calibrating the next: geometry (parallaxes, a water maser, eclipsing binaries) → Cepheids → type Ia supernovae → galaxies in the Hubble flow, whose distance against redshift gives H₀ |
| **Cepheid** | A pulsating supergiant star whose brightness varies with a period of days to about 100 days. The longer the period, the more luminous the star (Leavitt law, 1912), so the period gives its true luminosity and, compared with its apparent brightness, its distance. Bright enough to be seen in galaxies up to about 40 Mpc with the Hubble and James Webb telescopes: the middle rung of the distance ladder |
| **Type Ia supernova** | The thermonuclear explosion of a white dwarf. Its peak luminosity is nearly the same each time and is corrected with the width of its light curve (brighter ones fade more slowly), so it is a standard candle visible to gigaparsecs. Cepheids in the same host galaxies calibrate it |
| **SH0ES** | "Supernovae and H₀ for the Equation of State of dark energy" (written with a zero), the team led by A. G. Riess that measures H₀ with the distance ladder: H₀ = 73.04 ± 1.04 km/s/Mpc (Riess et al. 2022 [\[34\]](../references.md#ref-34)). It reflects the late, local Universe |
| **Hubble tension** | The 5σ disagreement between the Planck (67.4) and SH0ES (73.0) values of H₀: either an unknown systematic error in one of them, or physics beyond ΛCDM. GW sirens measure distance without the ladder and without an early-Universe model, so enough of them can arbitrate; GW170817 alone (about 70, with a wide interval) cannot yet |
| **PLP** | Power Law + Peak: BBH primary-mass model, a power law with a Gaussian peak and a smooth low-mass turn-on [\[15\]](../references.md#ref-15) [\[11\]](../references.md#ref-11) |
| **MLTP** | Multi Peak model: power law with two Gaussian peaks, found near 9 and 27 M☉ in the GWTC-4.0 cosmology paper [\[29\]](../references.md#ref-29) and near 10 and 30 M☉ in the reproduction run |
| **FullPop-4.0** | GWTC-4.0 model of the full compact-binary population (neutron stars, mass gap, black holes) [\[13\]](../references.md#ref-13) [\[29\]](../references.md#ref-29) |
| **Madau–Dickinson** | Shape of the merger rate with redshift: rises as (1+z)^γ, peaks near z_p, falls as (1+z)^−κ [\[16\]](../references.md#ref-16) |
| **Galaxy catalog: GLADE+, DES-Y6** | GLADE+ is the all-sky galaxy compilation for GW host searches (VizieR `VII/291/gladep`), whose K band the GWTC-4.0 dark siren uses; DES-Y6 is the deeper Dark Energy Survey catalog, over 12% of the sky ([galaxy_catalog](../modes/galaxy-catalog.md)) |
| **Line-of-sight redshift prior** | For each sky pixel, the prior on the redshift of a GW source in that direction, from the galaxies of the catalog plus the completeness correction: the input of the dark-siren method |
| **In-catalog / out-of-catalog term, completeness** | The in-catalog term is the redshift prior from the galaxies of a pixel brighter than its magnitude threshold; the out-of-catalog term adds the fainter, missed galaxies, uniform in comoving volume. The completeness is the fraction of the galaxies (weighted by luminosity) the catalog holds at a given redshift |
| **Schechter function** | The luminosity function of galaxies: a power law of slope α cut off exponentially above the characteristic magnitude M* (K band of GLADE+: M* = −23.39, α = −1.09) |
| **Spectroscopic / photometric redshift, K-correction** | A spectroscopic redshift comes from spectral lines (σ_z down to 10⁻⁴), a photometric one from colours; the K-correction converts a magnitude observed in a redshifted band to the galaxy's rest frame |

### Searches for continuous waves and lensed pairs

| Term | Definition |
|---|---|
| **Targeted / narrowband / all-sky search** | A targeted search follows the phase of a known pulsar from its radio or X-ray timing; a narrowband search lets the frequency and spin-down wander slightly around the timing values; an all-sky search scans every sky position, frequency and spin-down for unknown sources |
| **Semi-coherent search** | Analyses short stretches of data coherently and adds their power incoherently: far cheaper than one coherent analysis of months of data, at some loss of sensitivity |
| **CW pipelines** | CWInPy (heterodyne to remove the source's phase evolution, then Bayesian sampling of the amplitude), 5n-vector, the F/G/D-statistic, PowerFlux, FrequencyHough, SOAP (a hidden Markov model, HMM, tracking a wandering frequency without templates), and the Band Sampled Data (BSD) framework of the boson-cloud search, whose detection statistic is the critical ratio ([Method](../modes/continuous-waves.md#method)) |
| **Posterior overlap** | ∫ p₁ p₂ / π dθ over masses, spins and sky position: large when the posteriors of two events coincide compared with the prior volume, the first statistic of [lens_pairs](../modes/lens-pairs.md) |
| **Lensing search tiers; GOLUM, Hanabi, Phazap** | Stages of the LVK lensing search: tier 1 is a fast screen (posterior overlap, Phazap for phase consistency); tiers 2–3 are joint PE of the two signals (GOLUM, Hanabi) |
| **IAS, OGC** | Independent catalogs of GW events from public data (the Institute for Advanced Study pipeline, the Open Gravitational-wave Catalog); the package does not read them |

### Software

| Term | Definition |
|---|---|
| **icarogw** | Python package for population and cosmology inference with GW events (Mastrogiovanni et al. 2024 [\[57\]](../references.md#ref-57)) |
| **gwcosmo** | The other LVK code for GW cosmology, run alongside icarogw; the published results combine the samples of both ([icarogw and gwcosmo](icarogw-gwcosmo.md)) |
| **bilby / dynesty** | Bayesian inference library (Ashton et al. 2019 [\[58\]](../references.md#ref-58)) / its nested sampler (Speagle 2020 [\[60\]](../references.md#ref-60)) |
| **ligo.skymap, healpy** | ligo.skymap is the LVK package that reads, plots and cross-matches sky maps; healpy (and mhealpy for multi-order maps) handles HEALPix maps |
| **gwpy, PyCBC, LAL** | gwpy fetches, filters and whitens strain and computes the q-transform; PyCBC is the search library used for the matched filter; LAL (LIGO Algorithm Library) supplies waveforms, detector geometry and time delays |
| **TDinf, ringdown, pyRing** | Time-domain inference codes of the LVK: TDinf fits the full waveform (NRSur7dq4) to data cut before the merger; ringdown and pyRing fit quasinormal modes to the data after it ([area_law](../modes/area-law.md)) |

## Physical parameters

Some symbols mean different things on different pages; the context says which:

- **α**: the boson–black-hole coupling (scalar_cloud), the slope of the χ_eff–q relation (spin_population), the slope of the primary-mass power law (PLP), the slope of the Schechter function (galaxy_catalog);
- **β**: the pairing slope of the secondary mass (PLP, FullPop-4.0), the detectable fraction as a function of H₀ (bright siren);
- **κ**: the slope of the rate at low redshift, R(z) = R₀ (1+z)^κ with κ = 2.9 (rates), the high-redshift fall of the Madau–Dickinson shape (hubble_constant), the spin-induced quadrupole coefficient (tidal deformability);
- **λ**: the fraction in the Gaussian peak (PLP), the dimensional tidal deformability (Q_ij = −λ E_ij);
- **Λ**: the population hyperparameters (spectral siren, selection effects), the tidal deformability (neutron_star_eos), the likelihood-ratio statistic of scalar_cloud;
- **μ, σ**: the mean and width of a mass peak or of the χ_eff distribution; also a lens magnification (μ±) and the velocity dispersion of a lens galaxy (σ);
- **ψ**: the polarization angle, the Madau–Dickinson shape ψ(z);
- **ε**: the ellipticity of a neutron star (continuous_waves), the exponent of the luminosity weight of galaxies (galaxy_catalog).

### Masses

| Parameter | Definition |
|---|---|
| **Primary / secondary mass** | The masses of the heavier (\(m_1\)) and the lighter (\(m_2 \le m_1\)) object of a binary; the labels follow the mass, not the history. Population mass models describe the distribution of \(m_1\) and the secondary through the mass ratio ([Masses](gw-signals.md#masses-primary-secondary-and-what-the-signal-measures)) |
| **Total mass** | \(M = m_1 + m_2\); for heavy binaries, whose merger and ringdown dominate the signal, better measured than the components |
| **Mass ratio** | \(q = m_2 / m_1\), from 0 to 1 (1 for equal masses); some PE files also give \(1/q\). Measured from higher-order effects, partly degenerate with the aligned spins |
| **Chirp mass** | \(\mathcal{M} = (m_1 m_2)^{3/5}/(m_1+m_2)^{1/5}\), the mass combination that sets the inspiral, and the best-measured mass (0.1% for GW170817) |
| **Source / detector frame** | Physical masses / redshifted masses \(m_\text{det} = (1+z)\, m_\text{src}\) seen by the detector |
| **Remnant (final) mass** | Mass of the black hole left by the merger, \(M_f = m_1 + m_2 - E_\text{rad}/c^2\): about 5% less than the total mass (GW150914: 66 → 63.1 M☉) |
| **E_rad, peak luminosity** | The energy radiated in gravitational waves (in M☉c², 1 M☉c² = 1.79 × 10⁵⁴ erg) and the maximum radiated power (10⁵⁶ erg/s; 3.6 × 10⁵⁶ erg/s for GW150914), both from numerical-relativity fits to the masses and spins |

### Spins

| Parameter | Definition |
|---|---|
| **Spin magnitudes a₁, a₂ (χ₁, χ₂), tilts θ₁, θ₂** | The dimensionless spin of each component, from 0 to 1 (below 0.4 for the neutron stars of the population models), and the angle between it and the orbital angular momentum. Isotropic spins have directions uniform on the sphere |
| **Effective spin χ_eff** | The mass-weighted spin along the orbital angular momentum, \(\chi_\text{eff} = (m_1 a_1 \cos\theta_1 + m_2 a_2 \cos\theta_2)/(m_1 + m_2)\), from −1 to 1: spins aligned with the orbit (χ_eff > 0) lengthen the inspiral. The best-measured spin combination |
| **Precession spin χ_p** | The spin component in the orbital plane that drives precession, from 0 to 1 |
| **phi_jl, phi_12** | Azimuthal spin angles: that of the orbital angular momentum around the total one, and that between the in-plane components of the two spins |
| **Final spin χ_f** | The dimensionless spin of the remnant black hole, about 0.69 for equal, non-spinning masses, from numerical-relativity fits |
| **Spin-population parameters μ₀, σ, α** | The model χ_eff \| q ~ N(μ₀ + α(q − 0.5), σ), truncated to [−1, 1]: μ₀ is the mean at q = 0.5, α the slope with the mass ratio (α < 0: unequal-mass binaries spin more), σ the width ([spin_population](../modes/spin-population.md)) |

### Orientation, position and time

| Parameter | Definition |
|---|---|
| **Inclination θ_JN, ι; viewing angle** | The angle between the line of sight and the total (θ_JN) or orbital (ι) angular momentum, from 0° (face-on) through 90° (edge-on) to 180°; the viewing angle folds it into 0–90° |
| **Polarization angle ψ, phase φ** | ψ is the orientation of the binary's plus and cross axes on the sky, in [0, π); φ is the orbital phase at the reference frequency, in [0, 2π) |
| **ra, dec** | Right ascension and declination of the source, in radians in the PE files |
| **geocent_time, <det>_time** | The GPS time of the merger at the Earth's centre, and its arrival time at each detector, from it, ra and dec with the detector delays |
| **Truncation time t_<, ringdown start t_>** | The time before the peak at which the inspiral-only data are cut, and the time after it at which the ringdown fit starts, in units of GM/c³ of the total detector-frame mass M or of the remnant mass M_f ([area_law](../modes/area-law.md)) |

### Distance, redshift and cosmology

| Parameter | Definition |
|---|---|
| **z** | Redshift: the stretch of wavelengths and time scales by cosmic expansion |
| **D_L** | Luminosity distance, measured by the GW amplitude |
| **H₀** | Hubble constant, the present expansion rate of the Universe (km/s/Mpc) |
| **Ω_m** | Present matter density as a fraction of the critical density |
| **E(z)** | The expansion rate in units of its present value, H(z)/H₀; in flat ΛCDM \(\sqrt{\Omega_m (1+z)^3 + 1 - \Omega_m}\) |
| **Comoving volume V_c, dV_c/dz** | Volume that expands with the Universe; merger rates are given per unit comoving volume. dV_c/dz is the comoving volume of the shell between z and z + dz over the whole sky (Gpc³) |
| **Detector-frame Jacobian** | \((1+z)^2\, dD_L/dz\): the factor that carries a density from source-frame variables (m₁, m₂, z) to the measured ones (m₁,det, m₂,det, D_L) |
| **Recession, peculiar and Hubble-flow velocities v_r, v_p, v_H** | The measured velocity of a host galaxy away from us, the part from its own motion, and the part from the expansion, v_H = v_r − v_p (km/s); for a nearby host z = v_H/c (NGC 4993: 3327 − 310 = 3017 km/s) |
| **σ_z; zcmb, zhelio** | σ_z is the redshift uncertainty of a galaxy or host; zhelio is the heliocentric redshift, zcmb the redshift in the rest frame of the cosmic microwave background |
| **Ks magnitude** | The near-infrared (2.2 μm) apparent magnitude of the 2MASS survey (Vega system): it traces stellar mass, the band of the GWTC-4.0 dark siren |
| **Magnitude threshold m_thr** | The apparent magnitude down to which a pixel's catalog is taken as complete, by default the median of its galaxies' magnitudes |
| **Luminosity weight ε** | The exponent of the weight L^ε of each galaxy as a possible host: 1, proportional to luminosity; 0, all galaxies equally likely |

### Populations, rates and selection

| Parameter | Definition |
|---|---|
| **Population model p_pop(θ \| Λ)** | The distribution of source parameters θ (masses, redshift, spins) in the Universe, set by the hyperparameters Λ |
| **PLP parameters α, β, m_min, m_max, δ_m, λ, μ, σ** | α is the slope of the power law in m₁, p(m₁) ∝ m₁^−α, and β that of the pairing, p(m₂ \| m₁) ∝ m₂^β; m_min and m_max bound the black-hole masses, δ_m is the width of the smooth turn-on above m_min; λ is the fraction in the Gaussian peak, μ and σ its position and width (M☉). GWTC-3 values: 3.4, 1.1, 5.1, 87, 4.8, 0.04, 34, 3.6 |
| **MLTP parameters λ_g, λ_low, μ_low, σ_low, μ_high, σ_high** | λ_g is the fraction in the two Gaussian peaks of the Multi Peak model, λ_low the share of the lower one; μ and σ are the positions and widths of the two peaks (M☉) |
| **Madau–Dickinson γ, κ, z_p** | The shape \(\psi(z) \propto (1+z)^\gamma / [1 + ((1+z)/(1+z_p))^{\gamma+\kappa}]\), normalized to 1 at z = 0: it rises as (1+z)^γ, peaks at \((1+z_p)(\gamma/\kappa)^{1/(\gamma+\kappa)} - 1\), close to z_p, and falls as (1+z)^−κ; the star-formation history of Madau & Dickinson 2014 has γ = 2.7, κ = 2.9, z_p = 1.9 |
| **R, R₀, R(z)** | The merger rate density: mergers per unit comoving volume and source-frame time (Gpc⁻³ yr⁻¹); R₀ its local value, R(z) its evolution with redshift, quoted at z = 0.2 for BBHs |
| **VT, ⟨VT⟩** | Sensitive volume × time of a population: the expected number of detections is R ⟨VT⟩ |
| **ξ** | Detectable fraction of a population |
| **p_draw** | The known density from which the injections were drawn |
| **N_gen, T, mixture weight** | The number of injections generated, the analysis time they cover, and the weight that combines the injection sets of several runs into one file |

### Neutron stars and tides

| Parameter | Definition |
|---|---|
| **Tidal deformability Λ, Λ₁, Λ₂; Love number k₂** | How much a star is deformed by its companion's tidal field, \(\Lambda = \tfrac{2}{3} k_2 (Rc^2/Gm)^5\), with k₂ the Love number and R the radius: 0 for a black hole, about 100–1000 for a 1.4 M☉ neutron star; Λ₁, Λ₂ for the two stars ([Tidal deformability](tidal-deformability.md)) |
| **Λ̃, δΛ̃** | Λ̃ is the mass-weighted combination of Λ₁ and Λ₂ that enters the phase at leading (5PN) order, the one the data constrain; δΛ̃ is the second combination, essentially unconstrained |
| **Λ₁.₄, R₁.₄** | The tidal deformability and the radius (km) of a 1.4 M☉ neutron star, which summarize the EOS, related by Λ₁.₄ ≈ 2.88 × 10⁻⁶ (R₁.₄/km)^7.5 ([neutron_star_eos](../modes/neutron-star-eos.md)) |

### Continuous waves and boson clouds

| Parameter | Definition |
|---|---|
| **h₀** | The dimensionless strain amplitude of a continuous wave; searches give 95% upper limits on it, at twice the rotation frequency of a pulsar |
| **Spin-down limit h₀^sd** | The strain a pulsar would emit if all the rotational energy it loses went into gravitational waves, \(h_0^\text{sd} = \frac{1}{d}\left(\frac{5GI_{zz}}{2c^3}\frac{\lvert\dot f_\text{rot}\rvert}{f_\text{rot}}\right)^{1/2}\); an upper limit below it is the first physically meaningful one |
| **Ellipticity ε** | The fractional deformation of a neutron star about its rotation axis, (I_xx − I_yy)/I_zz (a "mountain"); a non-zero ε makes it radiate at twice its rotation frequency |
| **Boson mass m_φ** | The rest mass of an ultralight boson (eV); its cloud annihilates into a continuous wave at f = 2 m_φ c²/h = 483.6 Hz × (m_φ / 10⁻¹² eV) |
| **Coupling α** | α = G M m_φ/(ħc), the size of the black hole over the boson's wavelength; the cloud radius is r_φ = M/α² (G = c = 1) |
| **Cloud density ρ̄_φ** | The mean density of the scalar field within r_φ (g/cm³); 0 is vacuum |
| **τ_d** | The delay between the formation of the black hole and the merger, used by the superradiance prior (default 10⁶ yr) |
| **Λ (scalar_cloud)** | 2 (ln L_cloud − ln L_vacuum), the likelihood-ratio statistic of the cloud model over vacuum (≥ 0), calibrated on off-source noise: not the tidal Λ |

### Stochastic background

| Parameter | Definition |
|---|---|
| **Ω_GW(f)** | The energy density of the background per logarithmic frequency, in units of the critical density \(\rho_c c^2 = 3H_0^2c^2/(8\pi G)\); dimensionless, it grows as f^2/3 in the inspiral band |
| **⟨dE/df_s⟩, ISCO** | The GW energy one merger radiates per unit source-frame frequency, averaged over the population; the inspiral spectra of BNS and NSBH are cut at the ISCO, the innermost stable circular orbit |
| **z_h** | The redshift of the farthest detected events, beyond which the merger rate is extrapolated with the star-formation shape ([stochastic](../modes/stochastic.md)) |

### Strong lensing

| Parameter | Definition |
|---|---|
| **Magnification μ, μ±** | The factor by which a lens amplifies the energy flux of an image, so that its apparent distance is off by √μ; for a singular isothermal sphere μ± = 1/y ± 1, μ− = (1 − y)/y for the fainter image |
| **Singular isothermal sphere: σ, θ_E, y** | The lens model, a galaxy of density ∝ r⁻² set by its velocity dispersion σ (km/s); its Einstein radius is θ_E = 4π(σ/c)² D_ls/D_s, and y is the source position in units of θ_E (two images for y < 1) |
| **Time delay Δt** | The difference of the arrival times of two lensed images, in days |
