# Selection effects and injections

## The catalogs are biased samples

Detections are not a fair sample of the binaries in the Universe. A heavy binary is louder and is
detected out to much larger distances than a light one, so heavy binaries are strongly
over-represented, and at a given mass nearby sources are favored. A population inferred from the raw
catalog, a rate obtained by counting events, or a redshift inferred from the population, would all be
wrong.

The correction needs the probability \(P_\text{det}(\theta)\) that the whole detection chain
(detectors, noise, search pipelines, thresholds) detects a source of parameters \(\theta\). For real
searches on real, non-stationary data with glitches, it cannot be computed analytically: the LVK
measures it by simulation.

## Search-sensitivity injections

The LVK draws a large number of simulated signals, the **injections**, with masses, spins, redshift,
sky position and orientation drawn from a known density \(p_\text{draw}(\theta)\). For O3 and O4 the
signals are added to the real detector data, and the same search pipelines that found the real events
analyse them; each injection is recorded with the false-alarm rate (FAR) assigned by every pipeline.
For O1 and O2 the LVK uses a **semi-analytic** estimate: the SNR each signal would have had in the
measured noise spectrum of the time. The campaigns and their file format are described in
Essick et al. 2025 [\[35\]](../references.md#ref-35); the files are on Zenodo
([Data sources](../data-sources.md#search-sensitivity-injections)).

Each injection file holds, per injection:

- the source parameters;
- the draw density \(\ln p_\text{draw}\) (in source-frame masses, redshift and Cartesian spin
  components);
- a mixture weight: each observing period has its own draw, and the weights combine them into one set;
- the FAR of each search, or the semi-analytic SNR;

and, globally, the number of signals drawn \(N_\text{gen}\) (including the many too weak to be
recorded) and the total analysis time \(T\).

An injection is **found** if it passes the same criteria as the real events: lowest FAR below a
threshold for the real injections, network SNR above a threshold for the semi-analytic ones.

## Reweighting: one set of injections for every population

For a population model \(p_\text{pop}(\theta \mid \Lambda)\), the detectable fraction is

\[
\xi(\Lambda) = \int P_\text{det}(\theta)\, p_\text{pop}(\theta \mid \Lambda)\, d\theta
\;\approx\; \frac{1}{N_\text{gen}} \sum_{j \,\in\, \text{found}} \frac{p_\text{pop}(\theta_j \mid \Lambda)}{p_\text{draw}(\theta_j)} .
\]

The injections were drawn from \(p_\text{draw}\), not from the population being tested; weighting each
found injection by \(p_\text{pop}/p_\text{draw}\) corrects for this (importance sampling). One set of
injections thus serves for every population and every cosmology that an analysis tries, without
running the searches again. With the mixture weights \(w_j\), the draw density of injection \(j\) is
\(p_\text{draw}(\theta_j)/w_j\).

The **sensitive volume-time** of a population with a rate density \(R\) is

\[
\langle VT \rangle = T \cdot \frac{1}{N_\text{gen}} \sum_{j \,\in\, \text{found}} \frac{w_j\, p_\text{pop}(\theta_j)}{p_\text{draw}(\theta_j)},
\]

where \(p_\text{pop}\) includes the redshift distribution \(\frac{dV_c}{dz}\frac{1}{1+z}\) (see
[What a GW signal measures](gw-signals.md#comoving-volume-and-the-rate-of-mergers)), so that the expected
number of detections is \(N_\text{exp} = R \langle VT \rangle\).

## The effective number of injections

The sum is a Monte Carlo estimate. If a few injections dominate it (the population puts its weight
where few injections were found), the estimate is noisy. Its reliability is measured by the effective
sample size

\[
n_\text{eff} = \frac{\left(\sum_j x_j\right)^2}{\sum_j x_j^2}, \qquad x_j = \frac{w_j\, p_\text{pop}(\theta_j)}{p_\text{draw}(\theta_j)} .
\]

In a population analysis of \(N\) events, the selection term enters as \(\xi^N\), so its relative
error is multiplied by \(N\). Farr (2019) [\[34\]](../references.md#ref-34) showed that
\(n_\text{eff} > 4N\) keeps the resulting bias small; icarogw rejects the points of parameter space
where this fails. The same statistic, computed over the PE samples of each event, measures the
reliability of the per-event sums of the [hierarchical likelihood](spectral-siren.md#the-hierarchical-likelihood).

## A consequence: rare detections, common sources

Selection effects can invert the picture. Over O3–O4b there are 248 BBH detections for one BNS, yet
the two merger rates per unit volume are comparable: a BNS is detectable only within a few hundred
Mpc, a heavy BBH out to several Gpc, and the sensitive volume grows as the cube of the reach. In the
`rates` mode, with a constant rate per comoving volume for both, ⟨VT⟩ is 0.045 Gpc³ yr for BNS and
6.0 Gpc³ yr for BBH, a factor ~130 that offsets most of the factor 248 in the counts
([Merger rates](merger-rates.md)). The LVK population analyses even find the BNS rate likely higher
than the BBH one (GWTC-3: 10–1700 against 17.9–44 Gpc⁻³ yr⁻¹,
Abbott et al. 2023 [\[12\]](../references.md#ref-12)).
