"""Hawking's area law with GW250114, reproduced from the LVK data release of its discovery paper.

The area of a Kerr black hole of mass M and spin chi is A = 8 pi M^2 [1 + sqrt(1 - chi^2)] (G = c = 1). The
area law requires the remnant area A_f to exceed the sum A_i = A_1 + A_2 of the initial ones. A real test
measures the two from different parts of the signal:

- A_i from parameter estimation on data truncated before the peak (NRSur7dq4, 16 truncation times from
  -250 M to 0 M, M the total detector-frame mass);
- A_f from fits of the ringdown quasinormal modes to the data after the peak (the 220 mode, or 220 + 221),
  which give the remnant mass and spin through the Kerr spectrum alone, with no merger model.

The release (Zenodo 16877102, LVK 2025, arXiv:2509.08054) provides both, in detector-frame masses: the
redshift cancels in (A_f - A_i) / A_i. The significance is the paper's Gaussian one,
(<A_f> - <A_i>) / sqrt(sigma_f^2 + sigma_i^2), with the remnant areas of the two ringdown codes pooled; the
direct probability P(A_f < A_i) over random pairs of samples is given alongside.

Published: 4.4 sigma for t< = -40 M with the 220 mode from 10.5 M; at least 3.4 sigma for every
truncation from -250 M (the minimum of the paper's script, over all of them); above 5 sigma for t< >= -10 M;
3.6 sigma with 220 + 221 from 6 M.

The full-signal (IMR) PE gives A_f from fits to numerical relativity, which obey the law by construction:
it is shown only as a consistency check (`with_imr`).
"""
from __future__ import annotations

import re
import tarfile
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .report import write_simple_html_report

RECORD = "16877102"
ARCHIVE = "GW250114_data_release.tar.gz"
IMR_FILE = "posterior_samples_NRSur7dq4.h5"
MSUN_KM = 1.4766250614046494          # G M_sun / c^2 in km
REFERENCE_CUT = -40                   # t< of the paper's headline result (M)
REFERENCE_START = 10.5                # t> of its 220 ringdown (M_f)
OVERTONE_START = 6                    # earliest t> of the 220 + 221 model (M_f)
SEED = 250114                         # the paper's subsampling seed
PUBLISHED = dict(main=4.4, min_all_cuts=3.4, five_sigma_from=-10, overtone=3.6)
_MEMBERS = re.compile(r"(area_law_inspiral_data\.hdf5|ringdown_areas/.+\.hdf5|remnant_area_pyring_reweighted\.npy|"
                      r"area_change_prior\.dat)$")


def _log(msg: str) -> None:
    print(f"[area_law] {msg}", flush=True)


def kerr_area(mass, spin):
    """Horizon area of a Kerr black hole, in (G M_sun / c^2)^2 for a mass in M_sun."""
    spin = np.clip(np.asarray(spin, float), 0.0, 1.0)
    return 8 * np.pi * np.asarray(mass, float) ** 2 * (1 + np.sqrt(1 - spin ** 2))


def gaussian_significance(final: np.ndarray, initial: np.ndarray) -> float:
    """The paper's significance: (<A_f> - <A_i>) / sqrt(sigma_f^2 + sigma_i^2)."""
    return float((final.mean() - initial.mean()) / np.hypot(final.std(), initial.std()))


def pair_probability(final: np.ndarray, initial: np.ndarray, n: int = 400_000, seed: int = 1) -> float:
    """P(A_f < A_i) over random pairs of samples."""
    rng = np.random.default_rng(seed)
    return float(np.mean(rng.choice(final, n) < rng.choice(initial, n)))


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def fetch_release(cache_dir: Optional[str | Path] = None, with_imr: bool = False) -> Path:
    """Download the release once and extract the area-law files; returns the directory holding them."""
    from .data_repo import zenodo_cache_dir
    from .parameters_estimation import _download_http_with_progress

    root = Path(cache_dir).expanduser() if cache_dir else zenodo_cache_dir() / f"zenodo_{RECORD}_area_law"
    data = root / "data"
    if not (data / "area_law_inspiral_data.hdf5").exists():
        root.mkdir(parents=True, exist_ok=True)
        tgz = root / ARCHIVE
        if not tgz.exists():
            _log(f"downloading {ARCHIVE} (114 MB) from Zenodo {RECORD}")
            tmp = tgz.with_suffix(".part")
            _download_http_with_progress(f"https://zenodo.org/records/{RECORD}/files/{ARCHIVE}?download=1", tmp,
                                         desc=ARCHIVE)
            tmp.replace(tgz)
        with tarfile.open(tgz) as tar:
            members = [m for m in tar.getmembers() if m.isfile() and _MEMBERS.search(m.name)]
            for m in members:
                m.name = "data/" + m.name.split("data/", 1)[-1]
            try:
                tar.extractall(root, members=members, filter="data")
            except TypeError:                     # Python < 3.10.12 / 3.11.4: no extraction filters
                tar.extractall(root, members=members)
        tgz.unlink()
        _log(f"extracted {len(members)} files to {data}")
    if with_imr and not (data / IMR_FILE).exists():
        _log(f"downloading {IMR_FILE} (27 MB)")
        tmp = (data / IMR_FILE).with_suffix(".part")
        _download_http_with_progress(f"https://zenodo.org/records/{RECORD}/files/{IMR_FILE}?download=1", tmp,
                                     desc=IMR_FILE)
        tmp.replace(data / IMR_FILE)
    return data


def load(data: Path) -> dict:
    """Inspiral areas by truncation time, ringdown areas by model and start time, pyRing areas, prior."""
    import h5py

    out: dict = {"inspiral": {}, "ringdown": {"220": {}, "220+221": {}}}
    with h5py.File(data / "area_law_inspiral_data.hdf5", "r") as h:
        for i, t in enumerate(h["times"][:]):
            out["inspiral"][int(t)] = h[f"area_insp_{i}"][:]
    for p in sorted((data / "ringdown_areas").glob("*_final_mass_spin_area.hdf5")):
        m = re.match(r"(220(?:\+221)?)_([\d.]+)M_final_mass_spin_area\.hdf5$", p.name)
        if not m:
            continue
        with h5py.File(p, "r") as h:
            out["ringdown"][m.group(1)][float(m.group(2))] = h["Area_f"][:]
    out["pyring"] = np.load(data / "remnant_area_pyring_reweighted.npy")
    prior = np.loadtxt(data / "area_change_prior.dat")
    out["prior"] = (prior[0], prior[1])
    return out


def pooled_remnant(d: dict, start: float = REFERENCE_START, seed: int = SEED) -> np.ndarray:
    """Remnant areas of the 220 ringdown at `start`, pooled with pyRing in equal numbers (as the paper)."""
    rd = d["ringdown"]["220"][start]
    rng = np.random.default_rng(seed)
    n = min(len(rd), len(d["pyring"]))
    return np.concatenate([rng.choice(rd, n, replace=False), rng.choice(d["pyring"], n, replace=False)])


def imr_areas(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """(A_i, A_f) of the full-signal NRSur7dq4 PE, the remnant from its NR fits (circular)."""
    import h5py

    with h5py.File(path, "r") as h:
        g = next(h[k] for k in h if isinstance(h[k], h5py.Group) and "posterior_samples" in h[k])
        ps = g["posterior_samples"][()]
    ai = kerr_area(ps["mass_1"], ps["a_1"]) + kerr_area(ps["mass_2"], ps["a_2"])
    return ai, kerr_area(ps["final_mass"], ps["final_spin"])


# ---------------------------------------------------------------------------
# analysis
# ---------------------------------------------------------------------------
def analyse(d: dict) -> dict:
    """Headline result, truncation-time scan, ringdown start-time scans, overtone result."""
    ai = d["inspiral"][REFERENCE_CUT]
    af = pooled_remnant(d)
    rng = np.random.default_rng(1)
    ratio = rng.choice(af, 200_000) / rng.choice(ai, 200_000) - 1
    cut_scan = pd.DataFrame([dict(t_cut=t, significance=gaussian_significance(af, a), p_violation=pair_probability(af, a),
                                  initial_area_km2=float(np.median(a)) * MSUN_KM ** 2)
                             for t, a in sorted(d["inspiral"].items())])
    start_scan = pd.DataFrame([dict(model=m, t_start=t, significance=gaussian_significance(a, ai),
                                    p_violation=pair_probability(a, ai), final_area_km2=float(np.median(a)) * MSUN_KM ** 2)
                               for m in ("220", "220+221") for t, a in sorted(d["ringdown"][m].items())])
    five = cut_scan[cut_scan["significance"] >= 5]
    ot = d["ringdown"]["220+221"].get(OVERTONE_START)
    return dict(
        main=gaussian_significance(af, ai), main_p=pair_probability(af, ai),
        initial_km2=float(np.median(ai)) * MSUN_KM ** 2, final_km2=float(np.median(af)) * MSUN_KM ** 2,
        ratio=ratio, ratio_q=np.percentile(ratio, [5, 50, 95]),
        min_all_cuts=float(cut_scan["significance"].min()), five_sigma_from=int(five["t_cut"].min()) if len(five) else None,
        overtone=gaussian_significance(ot, ai) if ot is not None else None,
        cut_scan=cut_scan, start_scan=start_scan, initial=ai, final=af)


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------
_COL = dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", grid="#e4e3df", blue="#2a78d6", orange="#eb6834",
            green="#1a9e77", grey="#9b9a96")


def _style(ax):
    ax.set_facecolor(_COL["surface"])
    ax.grid(color=_COL["grid"], lw=0.7, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(_COL["grid"])
    ax.tick_params(colors=_COL["ink2"], labelsize=9)


def plot_areas(r: dict, prior: tuple, imr: Optional[tuple], out_png: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    c = _COL
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), dpi=150)
    fig.patch.set_facecolor(c["surface"])
    ax = axes[0]
    k2 = MSUN_KM ** 2 / 1e5
    bins = np.linspace(min(r["initial"].min(), r["final"].min()) * k2, max(np.percentile(r["final"], 99.8),
                       np.percentile(r["initial"], 99.8)) * k2, 80)
    ax.hist(r["initial"] * k2, bins=bins, density=True, color=c["blue"], alpha=0.8, lw=0, zorder=2,
            label=f"initial A₁ + A₂ (inspiral, t< = {REFERENCE_CUT} M)")
    ax.hist(r["final"] * k2, bins=bins, density=True, color=c["orange"], alpha=0.8, lw=0, zorder=2,
            label=f"remnant A_f (ringdown 220, t> = {REFERENCE_START:g} M_f)")
    ax.set_xlabel("Horizon area, detector frame (10⁵ km²)", color=c["ink2"])
    ax.set_yticks([])
    ax.legend(frameon=False, fontsize=8, labelcolor=c["ink2"], loc="upper right")
    ax = axes[1]
    x = np.linspace(-0.6, 2.0, 300)
    h, _, _ = ax.hist(r["ratio"], bins=x, density=True, color=c["orange"], alpha=0.85, lw=0, zorder=2,
                      label=f"inspiral vs ringdown: {r['main']:.1f}σ")
    if imr is not None:
        ax.hist(imr[1] / imr[0] - 1, bins=x, density=True, histtype="step", color=c["ink2"], lw=1.3, zorder=3,
                label="full-signal PE (NR fits: circular)")
    ax.plot(prior[0], prior[1], color=c["grey"], lw=1.2, ls=(0, (3, 2)), zorder=3, label="prior")
    ax.axvline(0, color=c["ink"], lw=1, zorder=4)
    ax.annotate("area law: A_f ≥ A_i", xy=(0, 1), xycoords=("data", "axes fraction"), xytext=(4, -12),
                textcoords="offset points", fontsize=8, color=c["ink2"])
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(0, 1.6 * h.max())          # the full-signal posterior is much narrower: clipped
    ax.set_xlabel("(A_f − A_i) / A_i", color=c["ink2"])
    ax.set_yticks([])
    ax.legend(frameon=False, fontsize=8, labelcolor=c["ink2"], loc="upper right")
    for a in axes:
        _style(a)
    fig.suptitle("GW250114: Hawking's area law from the inspiral and the ringdown", color=c["ink"], fontsize=11,
                 x=0.01, ha="left")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, facecolor=c["surface"])
    plt.close(fig)
    return out_png


def plot_scans(r: dict, out_png: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    c = _COL
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0), dpi=150, sharey=True)
    fig.patch.set_facecolor(c["surface"])
    cs = r["cut_scan"]
    ax = axes[0]
    ax.plot(cs["t_cut"], cs["significance"], "o-", color=c["blue"], ms=4, zorder=3, label="gwtc_analysis")
    ax.plot([REFERENCE_CUT], [PUBLISHED["main"]], "D", color=c["orange"], ms=7, zorder=4, label="published (4.4σ)")
    ax.axhline(PUBLISHED["min_all_cuts"], color=c["orange"], lw=0.9, ls=(0, (3, 2)), zorder=2,
               label="published minimum (3.4σ)")
    ax.set_xlabel("Inspiral truncation t< (M before the peak)", color=c["ink2"])
    ax.set_ylabel("Significance of A_f > A_i (σ)", color=c["ink2"])
    ax.legend(frameon=False, fontsize=8, labelcolor=c["ink2"], loc="upper left")
    ax = axes[1]
    ss = r["start_scan"]
    for m, col in (("220", c["blue"]), ("220+221", c["green"])):
        s = ss[ss["model"] == m]
        ax.plot(s["t_start"], s["significance"], "o-", color=col, ms=4, zorder=3, label=f"ringdown {m}")
    ax.plot([REFERENCE_START], [PUBLISHED["main"]], "D", color=c["orange"], ms=7, zorder=4, label="published 220 (4.4σ)")
    ax.plot([OVERTONE_START], [PUBLISHED["overtone"]], "s", color=c["orange"], ms=7, zorder=4,
            label="published 220+221 (3.6σ)")
    ax.set_xlabel(f"Ringdown start t> (M_f after the peak; inspiral at {REFERENCE_CUT} M)", color=c["ink2"])
    ax.legend(frameon=False, fontsize=8, labelcolor=c["ink2"], loc="upper right")
    for a in axes:
        _style(a)
    fig.suptitle("Area-law significance against the truncation of the inspiral and the start of the ringdown",
                 color=c["ink"], fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, facecolor=c["surface"])
    plt.close(fig)
    return out_png


# ---------------------------------------------------------------------------
# mode
# ---------------------------------------------------------------------------
def run_area_law(
    src_name: str = "GW250114",
    cache_dir: Optional[str | Path] = None,
    with_imr: bool = False,
    out_report_html: Optional[str | Path] = "area_law.html",
    out_summary_tsv: Optional[str | Path] = "area_law.tsv",
    plots_dir: str | Path = "area_law_plots",
) -> dict:
    """Reproduce the area-law test of GW250114 from the LVK release; see the module docstring."""
    if src_name.split("_")[0] != "GW250114":
        raise ValueError("the area-law test needs inspiral-only and ringdown-only analyses: only GW250114 has a "
                         "public release of them")
    data = fetch_release(cache_dir, with_imr=with_imr)
    d = load(data)
    r = analyse(d)
    d_prior = d["prior"]
    imr = imr_areas(data / IMR_FILE) if with_imr else None
    q = r["ratio_q"]
    _log(f"A_i = {r['initial_km2']:.3g} km², A_f = {r['final_km2']:.3g} km², (A_f - A_i)/A_i = {q[1]:.2f} "
         f"[{q[0]:.2f}, {q[2]:.2f}]; significance {r['main']:.2f}σ (published {PUBLISHED['main']}σ)")

    rows = [
        dict(result=f"A_f > A_i, t< = {REFERENCE_CUT} M, 220 from {REFERENCE_START:g} M_f (pooled)",
             gwtc_analysis=round(r["main"], 2), published=PUBLISHED["main"]),
        dict(result="minimum over the truncations -250 M to 0 M", gwtc_analysis=round(r["min_all_cuts"], 2),
             published=PUBLISHED["min_all_cuts"]),
        dict(result="earliest truncation above 5σ (M)", gwtc_analysis=r["five_sigma_from"],
             published=PUBLISHED["five_sigma_from"]),
        dict(result=f"220+221 from {OVERTONE_START} M_f", gwtc_analysis=None if r["overtone"] is None else
             round(r["overtone"], 2), published=PUBLISHED["overtone"]),
    ]
    table = pd.DataFrame(rows)
    if out_summary_tsv:
        Path(out_summary_tsv).parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(out_summary_tsv, sep="\t", index=False)
        r["cut_scan"].to_csv(Path(out_summary_tsv).with_suffix(".truncation.tsv"), sep="\t", index=False,
                             float_format="%.4g")
        r["start_scan"].to_csv(Path(out_summary_tsv).with_suffix(".ringdown.tsv"), sep="\t", index=False,
                               float_format="%.4g")
    if out_report_html:
        images = [plot_areas(r, d_prior, imr, Path(plots_dir) / "area_law_GW250114.png"),
                  plot_scans(r, Path(plots_dir) / "area_law_scans_GW250114.png")]
        paras = [
            f"<b>A<sub>f</sub> &gt; A<sub>i</sub> at {r['main']:.1f}σ</b> (published {PUBLISHED['main']}σ): the sum of "
            f"the initial horizon areas, {r['initial_km2']:.3g} km² (median), measured from the data up to "
            f"{-REFERENCE_CUT} M before the peak, against the remnant area, {r['final_km2']:.3g} km², from the 220 "
            f"ringdown mode starting {REFERENCE_START:g} M<sub>f</sub> after it. (A<sub>f</sub> − A<sub>i</sub>)/A<sub>i</sub> "
            f"= {q[1]:.2f} (90%: {q[0]:.2f} to {q[2]:.2f}). Areas are in the detector frame; the redshift cancels in "
            "the ratio.",
            "The two areas come from different parts of the signal: the inspiral analyses stop before the peak, and "
            "the ringdown fits use the Kerr quasinormal-mode spectrum alone, with no merger model. This is what "
            "makes it a test: the remnant of the full-signal PE comes from fits to numerical relativity, which obey "
            "the law by construction" + (" (shown in the plot for comparison only)." if imr is not None else "."),
            f"Significance as in the paper: (⟨A<sub>f</sub>⟩ − ⟨A<sub>i</sub>⟩)/√(σ<sub>f</sub>² + σ<sub>i</sub>²), "
            "with the remnant areas of the two ringdown codes pooled in equal numbers. The direct probability "
            f"P(A<sub>f</sub> &lt; A<sub>i</sub>) over random pairs of samples is {r['main_p']:.1e}, the tail of the "
            "distributions being less Gaussian.",
            f"Over the inspiral truncations from −250 M to 0 M the significance stays above {r['min_all_cuts']:.1f}σ and "
            f"exceeds 5σ from {r['five_sigma_from']} M on; with the 220 + 221 model from {OVERTONE_START} M<sub>f</sub> "
            f"(its earliest time of validity) it is {r['overtone']:.1f}σ. Later truncations keep more of the "
            "inspiral and measure A<sub>i</sub> better; later ringdown starts lose signal.",
            f"Data: LVK, GW250114 discovery paper release (Zenodo {RECORD}, arXiv:2509.08054).",
        ]
        tables = [("Comparison with the paper", table.to_html(index=False, na_rep="")),
                  ("Significance against the inspiral truncation", r["cut_scan"].to_html(
                      index=False, float_format=lambda x: f"{x:.3g}")),
                  ("Significance against the ringdown start", r["start_scan"].to_html(
                      index=False, float_format=lambda x: f"{x:.3g}"))]
        Path(out_report_html).parent.mkdir(parents=True, exist_ok=True)
        write_simple_html_report(out_report_html, title="Hawking's area law with GW250114", paragraphs=paras,
                                 images=images, tables=tables)
        _log(f"report written to {out_report_html}")
    return dict(table=table, **{k: r[k] for k in ("main", "min_all_cuts", "five_sigma_from", "overtone")})
