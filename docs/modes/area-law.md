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
significance is the paper's, \((\langle A_f\rangle - \langle A_i\rangle)/\sqrt{\sigma_f^2 + \sigma_i^2}\),
with the remnant areas of two ringdown codes pooled in equal numbers; the report also gives the direct
probability \(P(A_f < A_i)\) over random pairs of samples.

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
