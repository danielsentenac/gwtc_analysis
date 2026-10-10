# scalar_cloud

Is a binary black hole better fitted with a cloud of an ultralight scalar field around it? The mode refits the
public strain of an event with the **IMRPhenomXPHM_Scalar** waveform of Roy, Vicente, Aurrekoetxea, Clough and
Ferreira [\[100\]](../references.md#ref-100), against the vacuum IMRPhenomXPHM of the catalog, and calibrates the
improvement on off-source noise. It reproduces their analysis of **GW190728_064510**, the event for which they report
tentative evidence for a scalar environment.

```bash
gwtc_analysis scalar_cloud                                     # GW190728: prepare, best fit, report (laptop, ~2 h)
gwtc_analysis scalar_cloud --stages fit background sample --executor slurm --submit \
    --python /path/to/env/bin/python --slurm-option=--partition=htc --slurm-option=--licenses=sps
gwtc_analysis scalar_cloud --stages report                     # once the jobs are done
```

## The model

An ultralight boson of mass \(m_\phi\) around a binary of total mass \(M\) forms a bound structure of radius
\(r_\phi = M/\alpha^2\), with \(\alpha = G M m_\phi / (\hbar c)\) (the black-hole size over the boson's wavelength;
\(\alpha = 0.1\) is \(m_\phi \approx 6.5 \times 10^{-13}\) eV for \(M = 20\,M_\odot\)). In the late inspiral the
orbit is well inside it (\(\beta = R / r_\phi \ll 1\)) and the binary ionizes it: it loses angular momentum to the
field at the rate (the paper's Eq. 2 and Eq. 6 of its Supplemental Material, \(G = c = 1\))

\[
\dot L_\phi = 4\pi \bar\rho_\phi M^3 \frac{I(v, \alpha, q)}{v^4},
\qquad
I = \sum_{\ell, m > 0} \pi^2 \beta^{\frac14 + \frac{\ell}{2}}
\left[\frac{2(1+q_2)}{m}\right]^{\frac12 - \ell}
\left[1 + (-1)^\ell q_2^{\ell - 1}\right]^2
\left[\frac{Y_{\ell m}(\pi/2, 0)}{\Gamma(\ell + 3/2)}\right]^2 ,
\]

with \(\bar\rho_\phi\) the mean density of the field within \(r_\phi\), \(v = (M\Omega)^{1/3}\),
\(\beta = \alpha^2/v^2\), \(q_2 = m_2/M\), and the sum up to \(\ell = 4\) (the \(\ell = 2\) term dominates). Against
the gravitational-wave loss \(\dot L_\text{GW} = \frac{32}{5}\eta^2 M v^7\) the orbit shrinks faster by the factor
\(1 + \delta\), \(\delta = \frac{5\pi}{8}\,\bar\rho_\phi M^2 I/(\eta^2 v^{11})\), which grows as about
\(v^{-13.5}\): the effect sits at the lowest frequencies. In the stationary-phase approximation the phase of the
(2, 2) mode changes by

\[
\Delta\Psi_{22}(f) = 2\pi \int_f^{f_c} (f' - f)\, \Delta\!\left(\frac{dt}{df}\right)\!(f')\, df',
\qquad
\Delta\!\left(\frac{dt}{df}\right) = \left(\frac{dt}{df}\right)_\text{GW}\left[\frac{1}{1+\delta} - 1\right],
\]

zero, with its derivative, at the end of the inspiral \(f_c\) (\(M f_c = 0.018\)): only the inspiral changes. Each
\((\ell, m)\) mode gets \(\frac{|m|}{2}\Delta\Psi_{22}(2f/|m|)\).

| Choice | Here | Note |
|---|---|---|
| vacuum model | IMRPhenomXPHM of bilby, as the catalog | \(\bar\rho_\phi = 0\) reproduces it exactly (tested) |
| where the correction is applied | each mode of `SimInspiralChooseFDModes` (J frame, seen at θ_JN), rotated into the polarization frame of `SimInspiralChooseFDWaveform` by a fitted angle | the paper uses the co-precessing frame: identical for aligned spins, an approximation for precessing ones |
| \(1/(1+\delta)\) | kept, not linearised | identical while \(\delta \ll 1\); finite when the torque exceeds the GW loss |
| masses in \(\delta\) | source frame (\(\bar\rho_\phi\) is a physical density), redshift from the distance (Planck15) | the kinematics use the detector-frame masses |

The dephasing scales roughly as \(\bar\rho_\phi\, \alpha^{2.5}\), i.e. as \(\bar\rho_\phi\, m_\phi^{2.5}\) whatever the
mass of the binary. For GW190728 (\(M \approx 20\,M_\odot\)), \(\bar\rho_\phi = 10^6\) g/cm³ and \(\alpha = 0.1\)
dephase the (2, 2) mode by about 1.6 rad at 20 Hz and 0.13 rad at 30 Hz.

## The test

**Best fit.** The vacuum model is fitted from the catalog's best samples (Nelder–Mead in coordinates whitened by
the covariance of the 2000 highest-likelihood samples: the likelihood is too narrow for a plain optimizer). The cloud
model is then profiled on a grid of α and of the dephasing at the lowest frequency (0.5, 2, 8 and 30 rad, which sets
\(\bar\rho_\phi\)), re-fitting the binary at each point, and the best points are refined with the cloud free. The
statistic is

\[
\Lambda = 2\,(\ln L_\text{cloud} - \ln L_\text{vacuum}) \;\ge\; 0 .
\]

The cloud model contains vacuum, so Λ is positive even without a cloud, and its distribution is not the textbook
χ² with two degrees of freedom: the density sits at the edge of its range, and α does nothing when it vanishes.

**Background.** The vacuum best fit is injected into off-source GWOSC noise, within an hour of the event, and the
same procedure is run on each segment. The fraction of the background with Λ at least the event's is its p-value.

**Bayes factor.** Nested sampling (bilby, dynesty, acceptance-walk, naccept 60, dlogz 0.1) of the vacuum model
and of the cloud model with two priors on the field:

| Prior | α | \(\bar\rho_\phi\) |
|---|---|---|
| agnostic | uniform, \(0.008 \le \alpha \le (1+q)/2\) | uniform up to \(10^9\) g/cm³ (the paper does not give its upper end), cloud mass \(\bar\rho_\phi r_\phi^3 < 0.1\,M\) |
| superradiance, delay \(\tau_d\) | \(0.02\,x^{1/9} \le \alpha/(1+q) \le 0.11\,x^{1/15}\), \(x = \frac{m_1}{10 M_\odot}\frac{1\,\text{Myr}}{\tau_d}\): the cloud grew before the merger and did not decay | below the saturation mass, \(1.3 \times 10^8\,\text{g/cm}^3\, \frac{(1+q)^3}{8} \left[\frac{\alpha/(1+q)}{0.07}\right]^7 \left[\frac{10 M_\odot}{m_1}\right]^2\) |

The binary priors are those of the GWTC-2.1 analysis. The merger time is marginalized; the distance is sampled
(the correction depends on the redshift, so the waveform is not simply proportional to \(1/d_L\)).

## Limits

- The torque is Newtonian and valid for \(\beta \ll 1\); close to the merger the model is extrapolated, and the
  correction is switched off there.
- A large Λ means that the vacuum model fits the inspiral worse than usual. Eccentricity, waveform systematics,
  calibration or a noise feature at low frequency can do the same: Λ is a screening statistic, the Bayes factor
  under a physical prior the actual claim.
- This likelihood gives ln L = 89.6 at the catalog's best sample of GW190728, against 96.8 in the PE file: the phases
  agree in every detector, the difference comes from the noise-spectrum, window and calibration handling. Both
  models use the same likelihood, so the comparison is not affected.

## Options

| Option | Default | Meaning |
|---|---|---|
| `--event` | `GW190728_064510` | the event (name of the GWTC PE releases) |
| `--pe-file`, `--pe-label` | the event's Zenodo file, `C01:IMRPhenomXPHM` | settings, noise spectra and starting points |
| `--workdir` | `scalar_cloud` | work directory |
| `--stages` | `prepare fit report` | among `prepare fit background sample report` |
| `--sampling-frequency` | 2048 | Hz |
| `--n-background` | 150 | off-source injections |
| `--jobs` | 1 | local parallel processes (fit, background), sampler pool (sample) |
| `--sample-models` | `vacuum agnostic superradiance` | models of the sample stage |
| `--tau-d` | 10⁶ | superradiance prior: years between black-hole formation and merger |
| `--executor` | `local` | `slurm`: the heavy stages as sbatch scripts in `<workdir>/slurm` and `<workdir>/slurm_sample` |
| `--python`, `--slurm-option`, `--slurm-env-setup`, `--submit`, `--npool` | | Slurm: the jobs' interpreter, extra #SBATCH options, setup lines, submission, CPUs of each sampling job |
| `--out-report`, `--out-summary` | `scalar_cloud.html`, `scalar_cloud.tsv` | report and table of the fits |

Cost for GW190728 (16 s at 2048 Hz): a vacuum likelihood 0.035 s, a cloud likelihood 0.2 s (the modes of LAL). The
best fit takes about 2 CPU-hours, each background injection about 2, the nested sampling of the cloud models days of
CPU (a cluster job).
