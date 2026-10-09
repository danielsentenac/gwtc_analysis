# Glossary

| Term | Definition |
|---|---|
| **GW** | Gravitational wave |
| **LVK** | LIGO–Virgo–KAGRA, the collaboration of the ground-based GW detectors (H1 Hanford, L1 Livingston, V1 Virgo, K1 KAGRA) |
| **GWTC** | Gravitational-Wave Transient Catalog: GWTC-1, -2.1, -3, -4.0, -5.0 |
| **GWOSC** | Gravitational Wave Open Science Center, the public GW data portal |
| **O1, O2, O3a, O3b, O4a, O4b** | Observing runs, 2015–2025; O3 and O4 are split in halves |
| **ER15** | Engineering run 15, just before O4a |
| **BBH / BNS / NSBH** | Binary black hole / binary neutron star / neutron star–black hole binary |
| **PE** | Parameter estimation: the Bayesian inference of one event's masses, spins, distance, sky position… |
| **PE label** | One analysis in a PE file (waveform, settings), e.g. `C00:IMRPhenomXPHM-SpinTaylor`; `Mixed` combines several |
| **PESummary** | The LVK file format and library for PE results (Hoy & Raymond 2021 [\[62\]](../references.md#ref-62)) |
| **PSD** | Power spectral density, the detector noise spectrum, used to whiten the data |
| **Calibration envelope** | Uncertainty of the detector calibration in amplitude and phase, marginalized in PE |
| **FAR** | False-alarm rate: how often noise alone produces a candidate at least this significant (per year) |
| **p_astro** | Probability that a candidate is astrophysical |
| **SNR** | Signal-to-noise ratio; the matched-filter SNR compares the data with a template |
| **Primary / secondary mass** | The masses of the heavier (\(m_1\)) and the lighter (\(m_2 \le m_1\)) object of a binary; the labels follow the mass, not the history. Population mass models describe the distribution of \(m_1\) and the secondary through the mass ratio ([Masses](gw-signals.md#masses-primary-secondary-and-what-the-signal-measures)) |
| **Total mass** | \(M = m_1 + m_2\); for heavy binaries, whose merger and ringdown dominate the signal, better measured than the components |
| **Mass ratio** | \(q = m_2 / m_1\), from 0 to 1 (1 for equal masses); some PE files also give \(1/q\). Measured from higher-order effects, partly degenerate with the aligned spins |
| **Chirp mass** | \(\mathcal{M} = (m_1 m_2)^{3/5}/(m_1+m_2)^{1/5}\), the mass combination that sets the inspiral, and the best-measured mass (0.1% for GW170817) |
| **Remnant (final) mass** | Mass of the black hole left by the merger, \(M_f = m_1 + m_2 - E_\text{rad}/c^2\): about 5% less than the total mass (GW150914: 66 → 63.1 M☉) |
| **Primordial black hole (PBH)** | A black hole formed in the first fraction of a second after the Big Bang from the collapse of a dense region, not from a star; any mass is possible, set by the formation time. A dark-matter candidate, because it forms before the light elements and is not counted as baryons; \(f_\text{PBH}\) is the fraction of the dark matter it would make up ([Primordial black holes](gw-signals.md#primordial-black-holes-and-sub-solar-masses)) |
| **Sub-solar mass (SSM)** | A compact object below ~1 M☉, which stars do not leave behind as neutron stars or black holes: the cleanest signature of a PBH merger. LVK alerts carry HasSSM, the probability that a component is below 1 M☉ |
| **Mass gaps** | Lower mass gap, about 3–5 M☉, between the heaviest neutron stars and the lightest black holes; pair-instability gap, above about 50 M☉, where stars should leave no black hole |
| **z** | Redshift: the stretch of wavelengths and time scales by cosmic expansion |
| **D_L** | Luminosity distance, measured by the GW amplitude |
| **Source / detector frame** | Physical masses / redshifted masses \(m_\text{det} = (1+z)\, m_\text{src}\) seen by the detector |
| **H₀** | Hubble constant, the present expansion rate of the Universe (km/s/Mpc) |
| **Planck 2015 cosmology** | Flat ΛCDM with the Planck 2015 parameters, the reference cosmology of the GW catalogs: `Planck15` (astropy: H₀ = 67.74, Ω_m = 0.3075) or `Planck15_LAL` (LAL: H₀ = 67.90, Ω_m = 0.3065) (Planck Collaboration 2016 [\[32\]](../references.md#ref-32)); see [How the redshift is obtained](redshift.md#the-planck-2015-cosmology) |
| **Planck (H₀)** | H₀ = 67.4 ± 0.5 km/s/Mpc, inferred from the cosmic microwave background measured by the Planck satellite, assuming flat ΛCDM (Planck Collaboration 2020 [\[33\]](../references.md#ref-33)). It reflects the early Universe: the size of the sound waves in the primordial plasma, carried forward to today by the model |
| **Distance ladder** | Distances built rung by rung, each calibrating the next: geometry (parallaxes, a water maser, eclipsing binaries) → Cepheids → type Ia supernovae → galaxies in the Hubble flow, whose distance against redshift gives H₀ |
| **Cepheid** | A pulsating supergiant star whose brightness varies with a period of days to about 100 days. The longer the period, the more luminous the star (Leavitt law, 1912), so the period gives its true luminosity and, compared with its apparent brightness, its distance. Bright enough to be seen in galaxies up to about 40 Mpc with the Hubble and James Webb telescopes: the middle rung of the distance ladder |
| **Type Ia supernova** | The thermonuclear explosion of a white dwarf. Its peak luminosity is nearly the same each time and is corrected with the width of its light curve (brighter ones fade more slowly), so it is a standard candle visible to gigaparsecs. Cepheids in the same host galaxies calibrate it |
| **SH0ES** | "Supernovae and H₀ for the Equation of State of dark energy" (written with a zero), the team led by A. G. Riess that measures H₀ with the distance ladder: H₀ = 73.04 ± 1.04 km/s/Mpc (Riess et al. 2022 [\[34\]](../references.md#ref-34)). It reflects the late, local Universe |
| **Hubble tension** | The 5σ disagreement between the Planck (67.4) and SH0ES (73.0) values of H₀: either an unknown systematic error in one of them, or physics beyond ΛCDM. GW sirens measure distance without the ladder and without an early-Universe model, so enough of them can arbitrate; GW170817 alone (about 70, with a wide interval) cannot yet |
| **ΛCDM, Ω_m** | Standard cosmological model; present matter density as a fraction of the critical density |
| **Comoving volume V_c** | Volume that expands with the Universe; merger rates are given per unit comoving volume |
| **Standard siren** | A GW source used as a distance indicator |
| **Bright / dark / spectral siren** | Redshift from an electromagnetic counterpart / a galaxy catalog / features of the mass distribution |
| **Injection** | A simulated signal added to the data (or evaluated semi-analytically) to measure the detection probability |
| **p_draw** | The known density from which the injections were drawn |
| **Found injection** | An injection that passes the detection criteria (FAR or SNR threshold) |
| **VT, ⟨VT⟩** | Sensitive volume × time of a population: the expected number of detections is R ⟨VT⟩ |
| **ξ** | Detectable fraction of a population |
| **n_eff** | Effective number of samples of a weighted Monte Carlo sum, \((\sum x)^2/\sum x^2\) |
| **Hierarchical inference** | Inference of population parameters from many events, each with its own uncertain parameters |
| **Hyperprior** | Prior on population (and cosmological) parameters |
| **PLP** | Power Law + Peak: BBH primary-mass model, a power law with a Gaussian peak and a smooth low-mass turn-on [\[15\]](../references.md#ref-15) [\[11\]](../references.md#ref-11) |
| **MLTP** | Multi Peak model: power law with two Gaussian peaks, found near 9 and 27 M☉ in the GWTC-4.0 cosmology paper [\[29\]](../references.md#ref-29) and near 10 and 30 M☉ in the reproduction run |
| **FullPop-4.0** | GWTC-4.0 model of the full compact-binary population (neutron stars, mass gap, black holes) [\[13\]](../references.md#ref-13) [\[29\]](../references.md#ref-29) |
| **Madau–Dickinson** | Shape of the merger rate with redshift: rises as (1+z)^γ, peaks near z_p, falls as (1+z)^−κ [\[16\]](../references.md#ref-16) |
| **Nested sampling** | Bayesian sampling that computes the evidence and the posterior together (dynesty) |
| **Live points** | The set of samples nested sampling evolves; more live points, finer exploration |
| **ln Z** | Log Bayesian evidence |
| **Seed** | Initialization of the random generator of one sampler run |
| **icarogw** | Python package for population and cosmology inference with GW events (Mastrogiovanni et al. 2024 [\[57\]](../references.md#ref-57)) |
| **bilby / dynesty** | Bayesian inference library (Ashton et al. 2019 [\[58\]](../references.md#ref-58)) / its nested sampler (Speagle 2020 [\[60\]](../references.md#ref-60)) |
| **IMR** | Inspiral–merger–ringdown, the three phases of a binary signal: the two bodies spiral in (the chirp), their horizons fuse (the peak of the signal), the remnant black hole settles down. An IMR waveform model covers all three in one template, and an IMR analysis, the standard catalog one, fits it to the whole signal; in that analysis the remnant comes from numerical-relativity fits, not from the ringdown alone (see [area_law](../modes/area-law.md)) |
| **Kerr black hole** | A rotating black hole of general relativity (Kerr 1963), the form an isolated, uncharged black hole settles into. Two numbers describe it completely, its mass M and its dimensionless spin χ = cJ/(GM²) from 0 (not rotating, the Schwarzschild black hole) to 1 (maximal): the no-hair theorem. Its horizon area is A = 8π (GM/c²)² (1 + √(1 − χ²)), used by [area_law](../modes/area-law.md); merger remnants typically have χ ≈ 0.7 |
| **Quasinormal modes** | The damped oscillations of the ringdown, like those of a struck bell; for a Kerr black hole their frequencies and damping times depend only on its mass and spin, so a ringdown-only fit measures the remnant |
| **IMRPhenomXPHM** | Frequency-domain waveform model with precession and higher modes (Pratten et al. 2021 [\[51\]](../references.md#ref-51)); the other models and their names: [Waveform models](waveforms.md) |
| **XPHM-SpinTaylor** | IMRPhenomXPHM with numerically evolved spin precession (Colleoni et al. 2024 [\[52\]](../references.md#ref-52)) |
| **q-transform** | Time–frequency representation of the strain, showing the chirp |
| **Whitening** | Dividing the data by the noise amplitude spectrum, so that all frequencies have equal noise |
