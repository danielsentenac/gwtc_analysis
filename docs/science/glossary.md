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
| **Chirp mass** | \(\mathcal{M} = (m_1 m_2)^{3/5}/(m_1+m_2)^{1/5}\), the mass combination that sets the inspiral |
| **z** | Redshift: the stretch of wavelengths and time scales by cosmic expansion |
| **D_L** | Luminosity distance, measured by the GW amplitude |
| **Source / detector frame** | Physical masses / redshifted masses \(m_\text{det} = (1+z)\, m_\text{src}\) seen by the detector |
| **H₀** | Hubble constant, the present expansion rate of the Universe (km/s/Mpc) |
| **Planck 2015 cosmology** | Flat ΛCDM with the Planck 2015 parameters, the reference cosmology of the GW catalogs: `Planck15` (astropy: H₀ = 67.74, Ω_m = 0.3075) or `Planck15_LAL` (LAL: H₀ = 67.90, Ω_m = 0.3065) (Planck Collaboration 2016 [\[32\]](../references.md#ref-32)); see [How the redshift is obtained](redshift.md#the-planck-2015-cosmology) |
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
| **IMRPhenomXPHM** | Frequency-domain waveform model with precession and higher modes (Pratten et al. 2021 [\[51\]](../references.md#ref-51)) |
| **XPHM-SpinTaylor** | IMRPhenomXPHM with numerically evolved spin precession (Colleoni et al. 2024 [\[52\]](../references.md#ref-52)) |
| **q-transform** | Time–frequency representation of the strain, showing the chirp |
| **Whitening** | Dividing the data by the noise amplitude spectrum, so that all frequencies have equal noise |
