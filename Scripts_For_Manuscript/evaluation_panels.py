"""Attach clinical metadata to saved/live model scores without fitting any model."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from model_evaluation_plots import evaluation_panel


SPLIT_LABELS = {"train": "Train (apparent)", "test": "Test", "val": "External validation"}
MARGIN_MODELS = {"LogisticRegression", "RidgeClassifier", "LinearSVC", "SVC", "GaussianNB"}
BROAD_INTERPRETATION_MODELS = MARGIN_MODELS | {
    "KNeighborsClassifier", "RandomForestClassifier", "ExtraTreesClassifier",
}


def metadata_by_split(root):
    """Read recorded HIV status; missing status is never recoded as negative."""
    result = {}
    for split in SPLIT_LABELS:
        frame = pd.read_csv(Path(root) / f"Files_For_Manuscript/metadata.{split}.tsv", sep="\t")
        if frame.sample_id.duplicated().any():
            raise ValueError(f"Duplicate library IDs in {split} metadata")
        result[split] = frame.set_index("sample_id")
    return result


def score_rows(predictions, metadata):
    """Convert feature-set/model/split arrays from the preceding notebook cell."""
    tables = []
    for feature_set, models in predictions.items():
        for model, splits in models.items():
            for split, values in splits.items():
                rows = metadata[split]
                tables.append(pd.DataFrame({
                    "feature_set": feature_set, "model": model, "split": split,
                    "sample_id": rows.index, "tb_positive": rows.tb_binary.to_numpy(),
                    "score": np.asarray(values),
                }))
    return pd.concat(tables, ignore_index=True)


def saved_scores(root, section):
    """Read the existing prediction exports for figure-only execution."""
    root = Path(root)
    source = root / "Figure_Output/Gene_Interactions/focused"
    if section == "paper":
        frame = pd.read_csv(root / "Files_For_Manuscript/gene_cutoff_reference/paper_model/r_predictions.csv",
                            float_precision="round_trip")
        return frame.rename(columns={"tb_binary": "tb_positive", "probability": "score"}).assign(
            feature_set="Six genes", model="Paper replica GFS")
    if section.startswith("gfs_") or section == "xgb_selected":
        frame = pd.read_csv(source / "gfs_prediction_scores.csv", float_precision="round_trip")
        features = {
            "gfs_paper": ["GFS: paper rule"],
            "gfs_cv": ["GFS: CV AUC"],
            "gfs_revised": ["GFS: paper rule + revisions", "GFS: CV AUC + revisions"],
            "xgb_selected": ["GFS: CV AUC"],
        }[section]
        frame = frame.loc[frame.feature_set.isin(features)]
        if section == "xgb_selected":
            frame = frame.loc[frame.model.eq("XGBClassifier")]
        return frame.copy()
    if section == "six_gene_panels":
        return pd.read_csv(source / "six_gene_panel_prediction_scores.csv", float_precision="round_trip")
    frame = pd.read_csv(source / "prediction_scores.csv", float_precision="round_trip")
    if section in {"six_genes", "logistic", "cutoff_logistic"}:
        frame = frame.loc[frame.feature_set.eq("Six genes")]
        if section != "six_genes":
            frame = frame.loc[frame.model.str.startswith("Logistic GLM:") | frame.model.eq("Paper replica GFS")]
    elif section in {"broad", "cutoff_broad"}:
        frame = frame.loc[frame.feature_set.isin(["DE genes", "All variable genes"])]
        if section == "cutoff_broad":
            frame = frame.loc[frame.model.isin(BROAD_INTERPRETATION_MODELS)]
    else:
        raise ValueError(f"Unknown evaluation section: {section}")
    if section.startswith("cutoff_"):
        frame = frame.loc[frame.split.eq("test")]
    return frame.copy()


def clinical_scores(root, frame):
    """Keep score units and apply original/default boundaries without tuning."""
    metadata = metadata_by_split(root)
    reference = json.loads((Path(root) / "Files_For_Manuscript/gene_cutoff_reference/paper_model/model.json").read_text())
    merged = []
    for split, rows in frame.groupby("split", sort=False):
        fields = metadata[split].loc[:, ["hiv", "tb_binary"]]
        joined = rows.merge(fields, left_on="sample_id", right_index=True, validate="many_to_one", how="left")
        if joined.tb_binary.isna().any() or not joined.tb_positive.eq(joined.tb_binary).all():
            raise ValueError(f"Prediction IDs or labels disagree with {split} metadata")
        merged.append(joined)
    frame = pd.concat(merged, ignore_index=True)
    recorded_hiv = frame.hiv.fillna("").astype(str).str.strip().str.lower()
    if not set(recorded_hiv) <= {"positive", "negative", "", "unknown"}:
        raise ValueError(f"Unrecognized HIV categories: {sorted(set(recorded_hiv))}")
    frame["hiv_status"] = recorded_hiv.map({"positive": "Positive", "negative": "Negative"}).fillna("Unknown")
    margin = frame.model.isin(MARGIN_MODELS) | frame.model.str.startswith("Logistic GLM:")
    frame["score_kind"] = np.where(margin, "decision", "probability")
    frame["threshold"] = np.where(margin, 0.0, 0.5)
    paper = frame.model.eq("Paper replica GFS") | (
        frame.model.eq("GFS") & frame.feature_set.eq("Paper GFS genes"))
    frame.loc[paper, "threshold"] = reference["threshold"]
    frame["classification_rule"] = "classifier default; not optimized for this plot"
    frame.loc[paper, "classification_rule"] = reference["threshold_rule"]
    frame["threshold_operator"] = ">"
    frame["y_pred"] = frame.score.gt(frame.threshold).astype(int)
    if "predicted" in frame and frame.predicted.notna().any():
        observed = frame.predicted.notna()
        np.testing.assert_array_equal(frame.loc[observed, "y_pred"], frame.loc[observed, "predicted"])
    frame["y_true"] = frame.tb_positive.astype(int)
    frame["split"] = frame.split.map(SPLIT_LABELS)
    if frame.feature_set.nunique() > 1:
        frame["model"] = frame.feature_set + " | " + frame.model
    return frame


def plot_comparison(root, section, predictions=None):
    """Create one combined panel for all models in the preceding evaluation cell."""
    titles = {
        "paper": "Paper six-gene GLM",
        "logistic": "Six-gene logistic models",
        "six_genes": "All classifiers on the paper six-gene panel",
        "broad": "Classifiers using DE genes and all variable genes",
        "gfs_paper": "Model-specific GFS: paper selection rule",
        "gfs_cv": "Model-specific GFS: fivefold training CV selection",
        "gfs_revised": "Model-specific GFS with feature revisions",
        "six_gene_panels": "Same classifiers on paper and XGBoost gene panels",
        "cutoff_logistic": "Test predictions of the five logistic models",
        "cutoff_broad": "Test predictions with DE genes and all variable genes",
        "xgb_selected": "XGBoost model used for the SHAP interpretation",
    }
    raw = saved_scores(root, section) if predictions is None else score_rows(predictions, metadata_by_split(root))
    frame = clinical_scores(root, raw)
    destination = Path(root) / "Figure_Output/Model_Evaluation"
    destination.mkdir(exist_ok=True)
    return evaluation_panel(
        frame, title=titles[section], output_prefix=destination / section,
        split_order=[label for label in SPLIT_LABELS.values() if label in set(frame.split)],
        focus_split="Test",
        aligned_annotations=True,
    )
