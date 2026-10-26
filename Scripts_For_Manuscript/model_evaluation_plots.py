"""Patient-level diagnostic plots from saved predictions; never fits a model.

The same small module is copied into each independent analysis repository.
Inputs retain their original scores, split roles, decisions and HIV metadata.
"""

from html import escape
from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


COLORS = {0: "#4C78A8", 1: "#C94B50"}
HIV_MARKERS = {"Negative": "o", "Positive": "D", "Unknown": "^", "Conflict": "X"}
HIV_SUFFIX = {"Negative": "", "Positive": " [HIV+]", "Unknown": " [HIV?]", "Conflict": " [HIV conflict]"}
REFERENCE = "https://pmc.ncbi.nlm.nih.gov/articles/PMC12377778/"


def normalize_hiv(value):
    """Normalize recorded HIV infection only; exposure is not infection."""
    if pd.isna(value):
        return "Unknown"
    text = str(value).strip().lower()
    if text in {"positive", "pos", "hiv+", "1", "1.0", "true"}:
        return "Positive"
    if text in {"negative", "neg", "hiv-", "0", "0.0", "false"}:
        return "Negative"
    if text == "conflict":
        return "Conflict"
    if text in {"", "na", "n/a", "nan", "none", "unknown", "missing", "not known", "not tested"}:
        return "Unknown"
    raise ValueError(f"Unrecognized HIV infection status: {value!r}; resolve its meaning explicitly")


def prepare_predictions(frame):
    """Fail on silent duplication, invalid labels or unknown status encodings."""
    needed = {"model", "split", "sample_id", "y_true", "score", "y_pred", "hiv_status"}
    missing = needed - set(frame)
    if missing:
        raise ValueError(f"Missing plot columns: {sorted(missing)}")
    out = frame.copy()
    if out.empty:
        raise ValueError("No predictions available for this panel")
    for key in ["model", "split", "sample_id"]:
        if out[key].isna().any() or out[key].astype(str).str.strip().eq("").any():
            raise ValueError(f"Missing {key}")
        out[key] = out[key].astype(str)
    for key in ["y_true", "score", "y_pred"]:
        out[key] = pd.to_numeric(out[key], errors="raise")
    if not out.y_true.isin([0, 1]).all() or not np.isfinite(out.score).all():
        raise ValueError("Plots require known binary TB labels and finite scores")
    if not out.y_pred.dropna().isin([0, 1]).all():
        raise ValueError("Decisions must be 0, 1 or missing")
    out["hiv_status"] = out.hiv_status.map(normalize_hiv)
    out["score_kind"] = out.get("score_kind", pd.Series("probability", index=out.index))
    out["threshold"] = pd.to_numeric(out.get("threshold", np.nan), errors="raise")
    probability = out.score_kind.eq("probability")
    if not out.loc[probability, "score"].between(0, 1).all():
        raise ValueError("A probability must be between zero and one; retain margins as decision scores")
    keys = ["model", "split", "sample_id"]
    if "repeat_id" in out:
        if out.repeat_id.isna().any():
            raise ValueError("Every repeated prediction needs an explicit repeat_id")
        keys.append("repeat_id")
    if out.duplicated(keys).any():
        raise ValueError("Duplicate predictions; do not count duplicate rows as repeated evaluations")
    if "subject_id" in out:
        subject_keys = ["model", "split", "subject_id"] + (["repeat_id"] if "repeat_id" in out else [])
        if out.subject_id.isna().any() or out.duplicated(subject_keys).any():
            raise ValueError("Use one library per participant in each evaluation")
    for _, group in out.groupby(["model", "split", "sample_id"], sort=False):
        if group.y_true.nunique() != 1 or group.hiv_status.nunique() != 1:
            raise ValueError("Clinical label or HIV status changes across repeated predictions")
    out["correct"] = np.where(out.y_pred.notna(), out.y_pred.eq(out.y_true).astype(float), np.nan)
    out["sample_label"] = out.sample_id + out.hiv_status.map(HIV_SUFFIX)
    return out


def _metrics(group):
    known = group.loc[group.y_pred.notna()]
    negative = known.loc[known.y_true.eq(0)]
    positive = known.loc[known.y_true.eq(1)]
    return dict(n_participants=group.sample_id.nunique(), n_evaluations=len(group),
                n_decisions=len(known), tb_positive=int(group.y_true.sum()),
                accuracy=known.correct.mean(), sensitivity=positive.y_pred.mean(),
                specificity=1 - negative.y_pred.mean(),
                roc_auc=roc_auc_score(group.y_true, group.score) if group.y_true.nunique() == 2 else np.nan)


def metric_table(frame):
    """Descriptive metrics with explicit denominators, HIV and age strata."""
    rows = []
    for (model, split), subset in frame.groupby(["model", "split"], sort=False):
        strata = [("All", "All", subset)]
        for status in HIV_MARKERS:
            strata.append(("All", status, subset.loc[subset.hiv_status.eq(status)]))
        if "age_group" in subset:
            for age, age_subset in subset.groupby("age_group", dropna=False, sort=False):
                strata.append((str(age), "All", age_subset))
                strata.extend((str(age), status, age_subset.loc[age_subset.hiv_status.eq(status)])
                              for status in HIV_MARKERS)
        for age, status, group in strata:
            rows.append(dict(model=model, split=split, age_group=age, hiv_status=status, **_metrics(group)))
    return pd.DataFrame(rows)


def _percent(value):
    return "NA" if pd.isna(value) else f"{value:.0%}"


def _score_axis(axis, group):
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", color="0.92", linewidth=0.6, zorder=0)
    if group.score_kind.eq("probability").all():
        axis.set_ylim(-0.035, 1.04)
        axis.yaxis.set_major_formatter(PercentFormatter(1))
        axis.set_ylabel("TB prediction score (%)")
    else:
        axis.set_ylabel("TB decision score")


def _threshold(axis, group):
    values = group.threshold.dropna().unique()
    if len(values) == 1:
        axis.axhline(values[0], color="0.2", linestyle="--", linewidth=1, zorder=2)


def _violin(axis, group, split):
    if group.empty:
        axis.axis("off")
        return
    rng = np.random.default_rng(87)
    counts = []
    for label in [0, 1]:
        subset = group.loc[group.y_true.eq(label)]
        counts.append(subset.sample_id.nunique())
        values = subset.score.to_numpy()
        if len(values) >= 3 and np.ptp(values) > 1e-10:
            violin = axis.violinplot(values, positions=[label], widths=0.78,
                                     showmeans=False, showextrema=False, points=80)
            for body in violin["bodies"]:
                body.set(facecolor=COLORS[label], edgecolor=COLORS[label], alpha=0.18)
        for status, marker in HIV_MARKERS.items():
            status_rows = subset.loc[subset.hiv_status.eq(status)]
            x = label + rng.uniform(-0.17, 0.17, len(status_rows))
            axis.scatter(x, status_rows.score, s=20 if status == "Negative" else 30,
                         c=COLORS[label], marker=marker, alpha=0.7, linewidths=0.5,
                         edgecolors="black" if status == "Positive" else COLORS[label], zorder=3)
    _score_axis(axis, group)
    _threshold(axis, group)
    axis.set_xlim(-0.55, 1.55)
    axis.set_xticks([0, 1], [f"TB negative\nn={counts[0]}", f"TB positive\nn={counts[1]}"])
    m = _metrics(group)
    stats = f"Accuracy {_percent(m['accuracy'])} | Sensitivity {_percent(m['sensitivity'])}\nSpecificity {_percent(m['specificity'])}"
    axis.set_title(textwrap.fill(str(split), 30) + "\n" + stats, fontsize=9.5, pad=9)


def _ranked_rows(group):
    repeated = group.duplicated("sample_id").any()
    if repeated:
        if "repeat_id" not in group or group.y_pred.isna().any():
            raise ValueError("Success fractions require real repeats with recorded decisions")
        rows = group.groupby("sample_id", sort=False).agg(
            value=("correct", "mean"), appearances=("repeat_id", "nunique"),
            y_true=("y_true", "first"), hiv_status=("hiv_status", "first"),
            sample_label=("sample_label", "first"), correct=("correct", "mean")).reset_index()
    else:
        rows = group.assign(value=group.score, appearances=1).copy()
    return rows.sort_values(["value", "sample_id"], ascending=[False, True]), repeated


def _ranked(axis, group, show_ids=False, aligned_annotations=False):
    rows, repeated = _ranked_rows(group)
    _draw_ranked(axis, rows, group, repeated, show_ids, aligned_annotations)


def _annotation_rows(axis, rows):
    """Align one annotation cell to each unchanged patient bar; no aggregation."""
    error_colors = [to_rgba("#242424" if value == 0 else "#dddddd" if pd.isna(value) else "white")
                    for value in rows.correct]
    hiv_colors = {"Positive": "#8056a3", "Negative": "white", "Unknown": "#c58c23", "Conflict": "#b32f3b"}
    tracks = [("Error", error_colors),
              ("HIV", [to_rgba(hiv_colors[status]) for status in rows.hiv_status])]
    for index, (label, colors) in enumerate(tracks):
        bottom = -.13 - index * .075
        track = axis.inset_axes([0, bottom, 1, .045])
        track.imshow(np.asarray(colors)[None, :, :], aspect="auto", interpolation="nearest",
                     extent=(-.5, len(rows) - .5, 0, 1))
        track.set_xlim(axis.get_xlim())
        track.set_axis_off()
        axis.text(-.015, bottom + .0225, label, transform=axis.transAxes,
                  ha="right", va="center", fontsize=9)


def _draw_ranked(axis, rows, group, repeated, show_ids, aligned_annotations=False):
    positions = np.arange(len(rows))
    axis.bar(positions, rows.value, color=[COLORS[int(v)] for v in rows.y_true], width=0.92, alpha=0.7)
    if not aligned_annotations:
        for status in ["Positive", "Unknown", "Conflict"]:
            selected = rows.hiv_status.eq(status).to_numpy()
            axis.scatter(positions[selected], rows.value.to_numpy()[selected], marker=HIV_MARKERS[status],
                         s=23, facecolors="white", edgecolors="black", linewidths=0.7, zorder=4)
    if repeated:
        axis.set_ylim(-0.035, 1.04)
        axis.set_ylabel("Correct test predictions (%)")
        axis.yaxis.set_major_formatter(PercentFormatter(1))
        counts_axis = axis.twinx()
        counts_axis.plot(positions, rows.appearances, color="0.2", linewidth=1)
        counts_axis.set_ylabel("Evaluation appearances")
        counts_axis.set_ylim(0, max(1, rows.appearances.max()) * 1.2)
    else:
        incorrect = rows.correct.eq(0).to_numpy()
        if not aligned_annotations:
            axis.scatter(positions[incorrect], rows.value.to_numpy()[incorrect], marker="x", color="black", s=18, zorder=5)
        _score_axis(axis, group)
        _threshold(axis, group)
    axis.spines[["top", "right"]].set_visible(False)
    axis.set_xlim(-1, max(1, len(rows)))
    if aligned_annotations:
        if repeated:
            raise ValueError("Aligned annotation rows require one prediction per person")
        _annotation_rows(axis, rows)
    if show_ids:
        axis.set_xticks(positions, rows.sample_label, rotation=75, ha="right", fontsize=7.5)
        for tick, status in zip(axis.get_xticklabels(), rows.hiv_status):
            if status == "Positive":
                tick.set_fontweight("bold")
    else:
        axis.set_xticks([])
    axis.set_xlabel("Sample ID; [HIV+] positive, [HIV?] unknown" if show_ids else
                    "Samples ranked by score; full IDs in the sample pages / table")
    if aligned_annotations:
        axis.set_xlabel("One bar per patient, ordered from highest to lowest score", labelpad=48)


def _focus(splits, requested):
    if requested is not None:
        if requested not in splits:
            raise ValueError(f"Requested split {requested!r} is not available")
        return requested
    for word in ["test", "new participant", "validation", "val", "oof", "tuning"]:
        match = [split for split in splits if word in split.lower() and "train" not in split.lower()]
        if match:
            return match[-1]
    return splits[-1]


def _save_details(frame, prefix, focus_split, title):
    """Export all exact values; PDF pages keep every focused sample ID readable."""
    frame.to_csv(str(prefix) + "_predictions.tsv", sep="\t", index=False)
    metrics = metric_table(frame)
    metrics.to_csv(str(prefix) + "_metrics.tsv", sep="\t", index=False)
    with PdfPages(str(prefix) + "_samples.pdf") as pdf:
        for model, group in frame.loc[frame.split.eq(focus_split)].groupby("model", sort=False):
            rows, repeated = _ranked_rows(group)
            for start in range(0, len(rows), 48):
                fig, axis = plt.subplots(figsize=(15, 7.5))
                _draw_ranked(axis, rows.iloc[start:start + 48], group, repeated, True)
                axis.set_title(f"{model} | {focus_split}\nRanked samples {start + 1}–{min(start + 48, len(rows))} / {len(rows)}", fontsize=11)
                fig.subplots_adjust(left=0.07, right=0.94, top=0.87, bottom=0.49)
                fig.text(0.07, 0.015, "Color: recorded TB class (blue negative / red positive). Diamond: HIV+. Triangle: HIV unknown. ×: wrong decision.\n"
                         "Single-run scores are not repeated-test success rates. Recorded scores/labels are unchanged; decisions use the documented boundary.", fontsize=9)
                pdf.savefig(fig)
                plt.close(fig)
    columns = [name for name in ["model", "split", "sample_label", "subject_id", "age_group", "hiv_status", "y_true", "score", "score_kind", "y_pred", "correct", "threshold", "threshold_operator", "classification_rule", "repeat_id"] if name in frame]
    # Offline, searchable exact-value companion; no external assets or patient-data upload.
    content = ("<!doctype html><meta charset='utf-8'><title>Patient evaluation details</title>"
               "<style>body{font:15px system-ui;margin:24px;color:#243746}table{border-collapse:collapse}"
               "th,td{padding:6px 10px;border-bottom:1px solid #ddd;text-align:left}th{position:sticky;top:0;background:#edf2f5}"
               "input{padding:8px;width:440px;max-width:90%}tr.hiv-positive{background:#fff5d8}</style>"
               f"<h1>{escape(title)}</h1><p>1 = TB positive; 0 = TB negative. [HIV+] = positive; [HIV?] = unknown. "
               "No HIV suffix = recorded negative. Correct: 1 yes, 0 no, blank decision unavailable. "
               "Scores are probabilities only when marked as such in the TSV. Missing HIV is never treated as negative.</p>"
               "<label>Find a sample, model or HIV status: <input id='filter' type='search'></label>"
               + frame[columns].to_html(index=False, border=0, table_id="predictions", na_rep="", escape=True)
               + "<script>const rows=[...document.querySelectorAll('#predictions tbody tr')];"
               "for(const r of rows){if(r.textContent.includes('[HIV+]'))r.classList.add('hiv-positive')}"
               "document.querySelector('#filter').addEventListener('input',e=>{const q=e.target.value.toLowerCase();"
               "for(const r of rows)r.hidden=!r.textContent.toLowerCase().includes(q)});</script>")
    Path(str(prefix) + "_samples.html").write_text(content)


def evaluation_panel(frame, *, title, output_prefix=None, split_order=None, focus_split=None,
                     aligned_annotations=False):
    """One combined figure per evaluated model set, with exact-value companions.

    Each model has split-specific score violins and one ranked evaluation panel.
    Actual repeated held-out predictions, if supplied, yield correct/repeat counts;
    one prediction per person yields scores plus wrong-decision markers instead.
    """
    data = prepare_predictions(frame)
    models = list(data.model.unique())
    splits = list(split_order) if split_order is not None else list(data.split.unique())
    if set(splits) != set(data.split.unique()):
        raise ValueError("split_order must contain every available split exactly once")
    focused = _focus(splits, focus_split)
    with plt.rc_context({"font.size": 9, "svg.fonttype": "none", "pdf.fonttype": 42}):
        figure, axes = plt.subplots(len(models), len(splits) + 1, squeeze=False,
                                   figsize=(4.0 * len(splits) + (14 if aligned_annotations else 10),
                                            (4.5 if aligned_annotations else 4.1) * len(models) + 1.4),
                                   gridspec_kw={"width_ratios": [1] * len(splits) + [3.7 if aligned_annotations else 2.4]})
        for row, model in enumerate(models):
            subset = data.loc[data.model.eq(model)]
            for col, split in enumerate(splits):
                group = subset.loc[subset.split.eq(split)]
                _violin(axes[row, col], group, split)
            group = subset.loc[subset.split.eq(focused)]
            if len(group):
                _ranked(axes[row, -1], group, show_ids=False, aligned_annotations=aligned_annotations)
                repeat = group.duplicated("sample_id").any()
                axes[row, -1].set_title(f"{focused}: " + ("correct / repeated evaluations" if repeat else "one prediction per person"), fontsize=10)
            else:
                axes[row, -1].axis("off")
            if not subset.score_kind.eq("probability").all():
                limits = np.r_[subset.score.to_numpy(), subset.threshold.dropna().to_numpy()]
                padding = max(0.1, np.ptp(limits) * 0.08)
                for axis in axes[row]:
                    axis.set_ylim(limits.min() - padding, limits.max() + padding)
        legend = [Line2D([], [], marker="o", linestyle="", color=COLORS[0], label="TB negative"),
                  Line2D([], [], marker="o", linestyle="", color=COLORS[1], label="TB positive"),
                  Line2D([], [], marker="D", linestyle="", color="black", markerfacecolor="white",
                         label="HIV positive (points)" if aligned_annotations else "HIV positive"),
                  Line2D([], [], marker="^", linestyle="", color="black", markerfacecolor="white",
                         label="HIV unknown (points)" if aligned_annotations else "HIV unknown"),
                  Line2D([], [], marker="|" if aligned_annotations else "x", linestyle="", color="black",
                         label="Wrong decision (error row)" if aligned_annotations else "Wrong decision (single run)"),
                  Line2D([], [], linestyle="--", color="black", label="Decision boundary (if fixed)")]
        if aligned_annotations:
            legend[4:4] = [Line2D([], [], marker="s", linestyle="", color="#8056a3", label="HIV row: positive"),
                           Line2D([], [], marker="s", linestyle="", color="#c58c23", label="HIV row: unknown")]
        figure.suptitle(title, fontsize=13, y=1 - 0.10 / figure.get_figheight())
        figure.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, 1 - 0.38 / figure.get_figheight()),
                      ncol=4 if aligned_annotations else 3, frameon=False)
        figure.subplots_adjust(top=1 - 1.8 / figure.get_figheight(),
                               bottom=(1.35 if aligned_annotations else 1.0) / figure.get_figheight(),
                               left=3.2 / figure.get_figwidth(), right=0.98, wspace=0.45, hspace=0.85)
        for row, model in enumerate(models):
            position = axes[row, 0].get_position()
            figure.text(position.x0 - 1.0 / figure.get_figwidth(), (position.y0 + position.y1) / 2,
                        textwrap.fill(model, 22), fontsize=10, weight="bold", va="center", ha="right")
        note = ("Below bars: Error black = wrong; HIV purple = positive, gold = unknown; blank = correct / HIV negative.\n"
                "Every patient and score is retained. Full sample IDs and HIV labels are in the linked details."
                if aligned_annotations else
                "Recorded scores/labels are unchanged. Decisions use the documented boundary and split roles.\n"
                "HIV comparisons are descriptive, not adjusted effects. Full IDs and denominators are exported.")
        figure.text(0.065, 0.12 / figure.get_figheight(), note, fontsize=8.5)
        if output_prefix is not None:
            prefix = Path(output_prefix)
            prefix.parent.mkdir(parents=True, exist_ok=True)
            for suffix in ["png", "pdf"]:
                figure.savefig(str(prefix) + "." + suffix, dpi=120, bbox_inches="tight")
            _save_details(data, prefix, focused, title)
    return figure
