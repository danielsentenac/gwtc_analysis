# area_law

Hawking's area law tested with **GW250114**, the loudest binary black hole so far (network SNR 80),
reproduced from the LVK data release of its discovery paper [\[79\]](../references.md#ref-79).

```bash
gwtc_analysis area_law               # downloads 114 MB once
gwtc_analysis area_law --with-imr    # also the full-signal PE, for comparison (27 MB)
```

## The test

In classical general relativity the total area of black-hole horizons never decreases (Hawking 1971
[\[80\]](../references.md#ref-80)). A merger should therefore leave a remnant whose area exceeds the sum
of the initial ones, \(A_f \ge A_1 + A_2\). The horizon area of a Kerr black hole of mass \(M\) and spin
\(\chi\) is

\[
A = 8\pi \left(\frac{GM}{c^2}\right)^2 \left(1 + \sqrt{1 - \chi^2}\right).
\]

A real test measures the two sides from **different parts of the signal**:

- **before:** \(A_i = A_1 + A_2\) from parameter estimation (NRSur7dq4) on data truncated at a time
  \(t_<\) before the peak, from −250 M to 0 M, with M the total detector-frame mass;
- **after:** \(A_f\) from fits of the ringdown quasinormal modes to the data from a time \(t_>\) after
  the peak (the fundamental 220 mode, or 220 and its first overtone 221). The mode frequencies and
  damping times give the remnant mass and spin through the Kerr spectrum alone, with no merger model.

The full-signal PE cannot test the law: its remnant comes from fits to numerical-relativity simulations,
which obey it by construction. `--with-imr` shows it only as a consistency check.

The masses are detector-frame on both sides, so the redshift cancels in \((A_f - A_i)/A_i\). The
significance is the paper's, \((\langle A_f\rangle - \langle A_i\rangle)/\sqrt{\sigma_f^2 + \sigma_i^2}\)
(a *z-score*: the difference of the means in units of its standard deviation, read as a one-sided Gaussian tail
probability, 4.4σ ≈ 5 × 10⁻⁶),
with the remnant areas of two ringdown codes pooled in equal numbers; the report also gives the direct
probability \(P(A_f < A_i)\) over random pairs of samples.

## Method

**The mode runs no parameter estimation and no ringdown fit.** Those analyses are the LVK's, run on the
strain of GW250114. The data release of the paper gives their result as **posterior samples of
the horizon areas**; the mode reads them and computes the test from them, in seconds once downloaded.

### What the LVK analyses did

| Side | Code | Data | What the release holds |
|---|---|---|---|
| before, \(A_i\) | **TDinf**, a time-domain inference code (emcee sampler), with the NRSur7dq4 waveform of a quasicircular binary (15 parameters) | the strain cut sharply at a time \(t_<\) before the peak, 16 times from −250 M to 0 M | for each \(t_<\), samples of \(A_1 + A_2\), from the masses and spins of the two black holes |
| after, \(A_f\) | **ringdown** and **pyRing**, two time-domain codes that fit damped sinusoids: the Kerr quasinormal modes 220, or 220 and 221 | the strain from a time \(t_>\) after the peak, in units of the remnant mass M\(_f\) | samples of \(A_f\), from the remnant mass and spin given by the mode frequencies and damping times: from ringdown for each model at 17 start times from 6 to 20 M\(_f\); from pyRing for the 220 mode from 10.5 M\(_f\) (a reweighted file, `remnant_area_pyring_reweighted.npy`) |
| prior | | | the prior of \((A_f - A_i)/A_i\) that these analyses imply |

### The time windows are the paper's

**Every time window is set by the LVK analyses, not by the mode.** The release holds samples only for the
windows the paper analysed, and the mode has no option to change them: another window, such as
\(t_< = -60\) M, would need a new parameter estimation on the strain. Times are in units of \(GM/c^3\)
(about 0.35 ms for the detector-frame total mass of GW250114, about 70 M\(_\odot\)) before the peak, and
of \(GM_f/c^3\) after it.

| Window | Values in the release | Reference of the paper (used by the mode) |
|---|---|---|
| \(t_<\), inspiral truncation | −250, −200, −150, −100, −75, −50, −45, −40, −35, −30, −25, −20, −15, −10, −5, 0 M | −40 M, the headline result |
| \(t_>\), ringdown start, ringdown code (220 and 220 + 221) | 6, 7, 8, 9, 9.5, 10, 10.5, 11, 12, …, 20 M\(_f\) | 10.5 M\(_f\) for the 220 mode; 6 M\(_f\) for 220 + 221 |
| \(t_>\), pyRing | 10.5 M\(_f\) (220 only) | 10.5 M\(_f\), pooled with ringdown |

The ranges follow the limits of each description of the signal:

- \(t_<\) stops at the peak (0 M): later, the inspiral data would contain the ringdown, and the two areas
  would no longer be independent. Earlier than −250 M (about 88 ms), too little of the signal is left: its
  loudest cycles are the last ones. The truncations are dense near the peak (every 5 M from −50 M), where
  the significance changes fastest.
- \(t_>\) starts at 6 M\(_f\): just after the peak the signal is not yet a sum of a few quasinormal modes,
  and even 220 + 221 is considered valid only from about 6 M\(_f\), the 220 mode alone later. After about
  20 M\(_f\) the ringdown, damped in about 12 M\(_f\), has decayed too much to measure \(A_f\).
- With the reference windows, the data from −40 M to +10.5 M\(_f\), the merger, are used by neither
  analysis.

The constants `REFERENCE_CUT`, `REFERENCE_START` and `OVERTONE_START` of `area_law.py` hold the
paper's reference windows; the scans run over all the released ones. See the paper
[\[79\]](../references.md#ref-79) for the reasons of its choices.

### Independence of the two sides

Neither analysis sees the merger, and neither uses the other: the test needs no model of the merger. Their one
shared assumption is that both parts of the signal are redshifted by the same factor, since the areas are compared
in the detector frame; the recoil kick of the remnant could break it, by at most 10% and only for large spins,
inconsistent with GW250114 (paper, Supplement IV).

### What the mode does

1. **Download** the 114 MB archive of the release once (Zenodo 16877102) and keep the 37 files of the area law.
2. **Pool the two ringdown codes** for the remnant area: as many samples of ringdown (220 from 10.5 M\(_f\)) as of
   pyRing, 2000 of each, drawn without replacement with the paper's random seed.
3. **Headline result**: \(A_i\) at \(t_< = -40\) M (the paper's reference window) against the pooled \(A_f\), with the significance
   \((\langle A_f\rangle - \langle A_i\rangle)/\sqrt{\sigma_f^2 + \sigma_i^2}\) (the paper's: it uses the means and
   widths only, since the far tails of the samples are unreliable), and the direct probability \(P(A_f < A_i)\)
   over 400 000 random pairs of samples.
4. **Area change**: \((A_f - A_i)/A_i\) over 200 000 random pairs (the paper takes all the pairs), with its
   median and 90% interval.
5. **Scans**: the same significance for each of the 16 truncation times \(t_<\) (against the pooled \(A_f\)), and
   for each ringdown start time \(t_>\) and model (against \(A_i\) at −40 M). pyRing is released only at 10.5 M\(_f\),
   so the start-time scans use ringdown alone. From them: the lowest significance over the truncations, the
   earliest truncation above 5σ, and the 220 + 221 result from 6 M\(_f\).
6. **Units**: the areas are in \((G M_\odot/c^2)^2\) in the release, shown in km² (detector frame).

Only `--with-imr` computes areas itself: from the masses and spins of the full-signal NRSur7dq4 samples, with
the Kerr formula above, as the consistency check described in [The test](#the-test).

## Results

| Result | gwtc_analysis | Published [\[79\]](../references.md#ref-79) |
|---|---|---|
| \(A_f > A_i\), \(t_< = -40\) M, 220 mode from \(t_> = 10.5\) M\(_f\) | **4.45σ** | 4.4σ |
| minimum over the truncations −250 M to 0 M | 3.36σ | 3.4σ |
| earliest truncation above 5σ | −10 M | −10 M |
| 220 + 221 from 6 M\(_f\), its earliest time of validity | 3.63σ | 3.6σ |

The initial area is \(2.82 \times 10^5\) km² and the remnant area \(4.85 \times 10^5\) km² (medians,
detector frame): \((A_f - A_i)/A_i = 0.73\) (90%: 0.41–1.16). The direct probability of a decrease is
\(10^{-4}\).

![Initial and remnant areas, and the fractional area change](../img/modes/area_law_GW250114.png)

*Left: the initial and remnant areas. Right: the area change against the prior; the full-signal PE
(black outline) is far narrower, but tells nothing about the law.*

![Significance against the truncation and the ringdown start](../img/modes/area_law_scans_GW250114.png)

*Keeping more of the inspiral (later \(t_<\)) measures \(A_i\) better; starting the ringdown later loses
signal, and the two-mode model costs precision.*

## Why only GW250114

The test needs inspiral-only and ringdown-only analyses, which the catalog PE files do not contain, and
enough SNR in both parts. Of the earlier results:

- **GW150914** (Isi et al. 2021 [\[81\]](../references.md#ref-81), 97% with overtones): no data release
  of the inspiral side.
- The **GWTC-3 tests of GR** release has inspiral (frequency cut, IMR consistency test) and pyRing
  ringdown posteriors for 11 O1–O3 events. Combining them is possible but is a new analysis: with a
  frequency cut instead of a time cut GW150914 gives 74% instead of 97%, and the weak events need quality
  cuts (an uninformative inspiral or ringdown can mimic a violation). Not implemented.

## Outputs

- `--out-report` (HTML): the result, the comparison with the paper, the scans, the plots;
- `--out-summary` (TSV): the comparison with the paper; `<name>.truncation.tsv` and `<name>.ringdown.tsv`:
  significance and median areas against \(t_<\) and \(t_>\);
- `--plots-dir`: `area_law_GW250114.png`, `area_law_scans_GW250114.png`.

The release is extracted (37 files) into the Zenodo cache, or `--cache-dir`.

## Options

| Option | Default | Meaning |
|---|---|---|
| `--src-name` | `GW250114` | the event; only GW250114 has a public release of inspiral-only and ringdown-only analyses ([why](#why-only-gw250114)) |
| `--with-imr` | off | also show the area change of the full-signal analysis, made with an **IMR** waveform (see below): a consistency check, not a test; downloads 27 MB more |
| `--cache-dir` | the Zenodo cache | where the release is extracted |
| `--out-report`, `--out-summary`, `--plots-dir` | `area_law.html`, `area_law.tsv`, `area_law_plots` | outputs ([above](#outputs)) |
| `--data-repo {zenodo,s3,galaxy}` | `zenodo` | where the GW250114 data release are downloaded from: their public sources, or a mirror holding the same files, `s3` (the MinIO bucket, faster) or `galaxy` (the public usegalaxy.org history), else the public source ([S3 mirror of the public files](../data-sources.md#s3-mirror-of-the-public-files)) |

**IMR** stands for *inspiral–merger–ringdown*, the three phases of a compact-binary coalescence:

1. **inspiral:** the two bodies orbit and draw closer as they radiate; frequency and amplitude rise (the chirp);
2. **merger:** the horizons touch and fuse, the peak of the signal, described accurately only by numerical
   relativity;
3. **ringdown:** the distorted remnant settles into a Kerr black hole, emitting damped quasinormal modes whose
   frequencies and damping times depend only on its mass and spin.

The term is used in two ways:

- an **IMR waveform** models all three phases in one template (the IMR of IMRPhenomXPHM; NRSur7dq4 and SEOBNRv5PHM
  are IMR models too), and the standard catalog analysis fits one to the whole signal;
- the **IMR consistency test** of general relativity estimates the masses and spins separately from the
  low-frequency part (mostly inspiral) and the high-frequency part (merger and ringdown), and checks that they
  predict the same remnant.

For the area law, the full-signal (IMR) analysis cannot be the test: its remnant mass and spin come from fits to
numerical-relativity simulations, which obey the law by construction. The test compares the inspiral-only and
ringdown-only analyses; `--with-imr` only shows that the full-signal remnant agrees with them. See also the
[glossary](../science/glossary.md) and [waveforms](../science/waveforms.md).
