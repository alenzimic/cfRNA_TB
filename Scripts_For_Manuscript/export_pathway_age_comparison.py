"""Publish age-pathway model results and rendered patient distributions.

Read the already executed overview analysis and manuscript pathway results.
No models are fitted and no participant-level data tables are copied to the export.
"""

import argparse
import hashlib
import json
import shutil
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

AGGREGATE_FILES = [
    "primary_tests.tsv",
    "primary_curves.tsv",
    "coverage.tsv",
    "score_scaling.tsv",
    "test_summary.tsv",
    "sensitivity_tests.tsv",
    "sensitivity_summary.tsv",
    "model_diagnostics.json",
    "sensitivity_diagnostics.json",
    "age_summary.tsv",
    "age_support.tsv",
    "participant_flow.tsv",
    "clinical_missingness.tsv",
    "run_summary.txt",
    "pathway_timeline_index.tsv",
    "summary_tests.tsv",
    "summary_curves.tsv",
    "summary_model_diagnostics.json",
    "summary_leave_one_study_out_curves.tsv",
    "summary_study_sensitivity.tsv",
]
FIGURES = [
    "age_support",
    "main_difference_curves",
    "main_group_curves",
    "other_study_difference_curves",
    "pathway_age_atlas",
    "pathway_coverage",
]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def timeline_figure_names(source):
    """Select only the focus figures declared by the executed timeline index."""
    index = pd.read_csv(source / "pathway_timeline_index.tsv", sep="\t")
    if index.pathway.duplicated().any() or not index.pathway.str.fullmatch(
        r"HALLMARK_[A-Z0-9_]+"
    ).all():
        raise ValueError("Unexpected or duplicate pathway in timeline index")
    shown = index.shown_in_notebook.astype(str).str.lower().map(
        {"true": True, "false": False}
    )
    if shown.isna().any():
        raise ValueError("Unexpected shown_in_notebook flag in timeline index")
    return [
        "timeline_" + pathway.removeprefix("HALLMARK_").lower()
        for pathway in index.loc[shown, "pathway"]
    ]


def compare_saved_results(ROOT, OVERVIEW, OUT):
    """Compare numerical results and verify participant overlap locally."""
    OLD = ROOT / "Figure_Output/Gene_Interactions/focused"
    NEW = OVERVIEW / "results/pathway_age"
    old_score = pd.read_csv(OLD / "pathway_score_tb_associations.csv")
    old_gsea = pd.read_csv(OLD / "pathway_enrichment.csv")
    new_test = pd.read_csv(NEW / "primary_tests.tsv", sep="\t")
    new_curve = pd.read_csv(NEW / "primary_curves.tsv", sep="\t")
    people = pd.read_csv(NEW / "participants.tsv", sep="\t")
    old_meta = pd.concat(
        [
            pd.read_csv(
                ROOT / f"Files_For_Manuscript/metadata.{split}.tsv", sep="\t"
            ).assign(split=split)
            for split in ["train", "test", "val"]
        ],
        ignore_index=True,
    )
    old_meta["sample_id_normalized"] = old_meta.sample_id.str.strip()
    people["old_id_normalized"] = people.original_sample_id.str.strip()
    assert old_meta.sample_id_normalized.is_unique
    assert people.loc[people.old_id_normalized.notna(), "old_id_normalized"].is_unique
    matched = old_meta.merge(
        people.loc[
            people.old_id_normalized.notna(),
            ["old_id_normalized", "subject_id", "cohort", "tb", "original_split"],
        ],
        left_on="sample_id_normalized",
        right_on="old_id_normalized",
        how="left",
        validate="one_to_one",
        suffixes=("_old", "_new"),
        indicator=True,
    )
    assert matched._merge.eq("both").all()
    assert matched.subject_id.is_unique
    assert matched.split.eq(matched.original_split).all()
    assert matched.tb_old.map({"positive": 1, "negative": 0}).eq(matched.tb_new).all()
    overlap_rows = (
        matched.groupby(["split", "cohort_new"])
        .size()
        .rename("overlapping_participants")
        .reset_index()
        .rename(columns={"cohort_new": "cohort"})
    )
    overlap_rows.to_csv(OUT / "participant_overlap_summary.csv", index=False)

    focus = old_score.pathway.drop_duplicates().tolist()
    base = old_score.loc[old_score.adjustment.eq("Cohort + composition")]
    rows = []
    for p in focus:
        row = {
            "pathway": p,
            "manuscript_score_definition": "Mean standardized gene expression; excludes XGBoost six; training-reference SD",
            "age_score_definition": "Single-sample directional rank enrichment; all set genes; pooled 812-person SD"
            if p in set(new_test.pathway)
            else "Not tested (Reactome set absent from Hallmark-only atlas)",
        }
        for split in ["train", "test", "val"]:
            r = base.loc[base.pathway.eq(p) & base.split.eq(split)].iloc[0]
            for field in ["tb_difference", "ci_lower", "ci_upper", "bh_fdr_30"]:
                row[f"old_{split}_{field}"] = r[field]
        relevant = new_test.loc[new_test.pathway.eq(p)]
        row["new_studies_with_any_TB_q_lt_005"] = ";".join(
            relevant.loc[relevant.any_TB_q.lt(0.05), "cohort"]
        ) or ("None" if len(relevant) else "Not tested")
        row["new_min_any_TB_q"] = relevant.any_TB_q.min() if len(relevant) else np.nan
        row["new_min_age_interaction_q"] = (
            relevant.age_interaction_q.min() if len(relevant) else np.nan
        )
        row["new_age_interaction_discoveries"] = (
            int(relevant.age_interaction_q.lt(0.05).sum()) if len(relevant) else np.nan
        )
        for study in ["PROTECT", "ADTB", "GHL", "R2D2", "TBSQUARE"]:
            test = relevant.loc[relevant.cohort.eq(study)]
            c = new_curve.loc[
                new_curve.pathway.eq(p)
                & new_curve.cohort.eq(study)
                & new_curve.supported
            ]
            row[f"new_{study}_any_TB_q"] = (
                test.any_TB_q.iloc[0] if len(test) else np.nan
            )
            row[f"new_{study}_supported_difference_min"] = (
                c.difference.min() if len(c) else np.nan
            )
            row[f"new_{study}_supported_difference_max"] = (
                c.difference.max() if len(c) else np.nan
            )
        rows.append(row)
    focus_table = pd.DataFrame(rows)
    focus_table.to_csv(OUT / "focused_pathway_comparison.csv", index=False)

    # Exhaustive shared-set comparison under a declared old sensitivity specification.
    robust_gsea = old_gsea.loc[
        old_gsea.adjustment.eq("Cohort + composition")
        & old_gsea.panel.eq("Exclude XGBoost six")
    ]
    shared = sorted(set(robust_gsea.pathway) & set(new_test.pathway))
    assert len(shared) == 49
    rows = []
    for p in shared:
        row = {
            "pathway": p,
            "old_adjustment": "Cohort + composition",
            "old_gene_panel": "Exclude XGBoost six",
        }
        for split in ["train", "test", "val"]:
            r = robust_gsea.loc[
                robust_gsea.pathway.eq(p) & robust_gsea.split.eq(split)
            ].iloc[0]
            row[f"old_{split}_GSEA_NES"] = r.nes
            row[f"old_{split}_GSEA_FDR"] = r.gsea_fdr
            row[f"old_{split}_genes_tested"] = int(r.pathway_genes_tested)
        relevant = new_test.loc[new_test.pathway.eq(p)]
        row["old_GSEA_significant_splits"] = (
            ";".join(
                robust_gsea.loc[
                    robust_gsea.pathway.eq(p) & robust_gsea.gsea_fdr.lt(0.05), "split"
                ]
            )
            or "None"
        )
        row["old_GSEA_all_three_positive"] = bool(
            robust_gsea.loc[robust_gsea.pathway.eq(p), "nes"].gt(0).all()
        )
        row["new_any_TB_significant_studies"] = (
            ";".join(relevant.loc[relevant.any_TB_q.lt(0.05), "cohort"]) or "None"
        )
        row["new_min_age_interaction_q"] = relevant.age_interaction_q.min()
        for study in ["PROTECT", "ADTB", "GHL", "R2D2", "TBSQUARE"]:
            test = relevant.loc[relevant.cohort.eq(study)].iloc[0]
            c = new_curve.loc[
                new_curve.pathway.eq(p)
                & new_curve.cohort.eq(study)
                & new_curve.supported
            ]
            row[f"new_{study}_any_TB_q"] = test.any_TB_q
            row[f"new_{study}_age_interaction_q"] = test.age_interaction_q
            row[f"new_{study}_supported_difference_min"] = c.difference.min()
            row[f"new_{study}_supported_difference_max"] = c.difference.max()
        rows.append(row)
    all_shared = pd.DataFrame(rows)
    all_shared.to_csv(OUT / "all_shared_hallmark_comparison.csv", index=False)
    new_hits = new_test.loc[new_test.any_TB_q.lt(0.05)].copy()
    new_hits["in_old_eligible_GSEA"] = new_hits.pathway.isin(shared)
    new_hits.to_csv(OUT / "age_analysis_any_TB_hits.csv", index=False)
    paths = [
        OLD / "pathway_score_tb_associations.csv",
        OLD / "pathway_enrichment.csv",
        OLD / "pathway_gene_program_tb_associations.csv",
        NEW / "primary_tests.tsv",
        NEW / "primary_curves.tsv",
        NEW / "participants.tsv",
    ] + [
        ROOT / f"Files_For_Manuscript/metadata.{s}.tsv"
        for s in ["train", "test", "val"]
    ]
    audit = {
        "inputs": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        "overlap": {
            "old_participants": len(old_meta),
            "new_participants": len(people),
            "old_in_new": len(matched),
            "new_not_old": len(people) - len(matched),
            "one_to_one_matching": "Trim whitespace from old sample_id and new original_sample_id; exact join",
            "all_split_labels_match": True,
            "all_binary_labels_match": True,
            "new_GHL_R2D2_TBSQUARE_exactly_old": bool(
                set(
                    people.loc[
                        people.cohort.isin(["GHL", "R2D2", "TBSQUARE"]), "subject_id"
                    ]
                )
                == set(matched.subject_id)
            ),
        },
        "comparison": {
            "old_focused_pathways": len(focus),
            "shared_eligible_Hallmarks": len(shared),
            "old_Reactome_neutrophil_absent_new": True,
            "new_pancreas_beta_cells_absent_old_eligible_GSEA": True,
            "age_interaction_tests": len(new_test),
            "age_interaction_discoveries_FDR05": int(
                new_test.age_interaction_q.lt(0.05).sum()
            ),
            "age_interaction_min_q": float(new_test.age_interaction_q.min()),
            "any_TB_discoveries_FDR05": int(new_test.any_TB_q.lt(0.05).sum()),
            "new_nonoverlap_studies_any_TB_discoveries": int(
                new_test.loc[new_test.cohort.isin(["PROTECT", "ADTB"]), "any_TB_q"]
                .lt(0.05)
                .sum()
            ),
            "age_omnibus_tests_do_not_have_a_single_signed_coefficient": True,
            "supported_curve_extrema_are_descriptive_not_tested": True,
        },
        "scope": "No new model fitting, p values, gene-set permutations or participant-level public outputs were generated. Tables compare previously saved estimates and multiplicity results. Old and new scores and tests have different definitions; no effect-size equivalence or replication claim.",
    }
    (OUT / "comparison_audit.json").write_text(json.dumps(audit, indent=2))

    return audit, focus_table


def write_readable_comparison(out, focus):
    rows = []
    for _, record in focus.iterrows():
        tested = pd.notna(record.new_min_age_interaction_q)
        rows.append(
            {
                "Pathway": record.pathway.removeprefix("HALLMARK_")
                .removeprefix("REACTOME_")
                .replace("_", " ")
                .title(),
                "Earlier TB difference, train / test / validation (training SD)": " / ".join(
                    f"{record[f'old_{split}_tb_difference']:+.2f}"
                    for split in ["train", "test", "val"]
                ),
                "Earlier adjusted P, same order": " / ".join(
                    f"{record[f'old_{split}_bh_fdr_30']:.3g}"
                    for split in ["train", "test", "val"]
                ),
                "Age analysis: studies with general TB evidence, FDR<0.05": record.new_studies_with_any_TB_q_lt_005,
                "New age-dependent TB evidence": (
                    f"{int(record.new_age_interaction_discoveries)} discoveries; minimum q={record.new_min_age_interaction_q:.3f}"
                    if tested
                    else "Not tested"
                ),
            }
        )
    pd.DataFrame(rows).to_csv(out / "readable_focused_comparison.csv", index=False)


def write_readme(out, audit):
    overlap, comparison = audit["overlap"], audit["comparison"]
    (out / "README.md").write_text(f"""# Pathway analysis and continuous-age comparison

This is an archived **aggregate result snapshot** from the executed local
`cfRNA_TB_overview/notebooks/02_pathway_enrichment_by_age.ipynb`, together with a
comparison against `Scripts_For_Manuscript/gene_interactions_pathways.ipynb`.
The snapshot contains model estimates, tests, coverage, aggregate counts,
figures and provenance. Faint patient dots show observed distributions in the
figures; the interactive explorer embeds them as raster images. It does not
include participant identifiers, individual score tables, raw counts or
sample-level QC records.

## Findings and interpretation

- The earlier manuscript analysis supports general TB-associated expression,
  most consistently for interferon-gamma after excluding the selected genes.
  The focused comparison table retains all five previously specified programs,
  including uncertain and untested comparisons. The exhaustive table retains
  all {comparison["shared_eligible_Hallmarks"]} shared eligible Hallmark sets.
- The age analysis found **{comparison["age_interaction_discoveries_FDR05"]} of
  {comparison["age_interaction_tests"]} age-interaction tests at FDR <0.05**
  (smallest q={comparison["age_interaction_min_q"]:.6f}). Its separate omnibus
  test for any TB difference identified {comparison["any_TB_discoveries_FDR05"]}
  pathway–study associations. General TB association does not demonstrate an
  age-dependent TB difference.
- All **{overlap["old_in_new"]} earlier participants overlap** the
  {overlap["new_participants"]}-person age analysis. GHL, R2D2 and TBSQUARE are
  exactly this overlapping group; the {overlap["new_not_old"]} additional people
  are ADTB and PROTECT. Of the {comparison["any_TB_discoveries_FDR05"]} general-TB
  discoveries, {comparison["any_TB_discoveries_FDR05"] - comparison["new_nonoverlap_studies_any_TB_discoveries"]}
  occur in overlapping adult studies and
  {comparison["new_nonoverlap_studies_any_TB_discoveries"]} occur in ADTB or PROTECT.
  Agreement is **same-participant analytical
  concordance**, not independent replication or confirmed transfer.
- Children and adults were sampled in different studies. The primary age curves are
  fitted separately by study and do not infer a continuous biological transition
  across unsampled adolescence. Sparse ages are masked, including internal gaps.

## Different scores and questions

| Feature | Earlier manuscript pathways | Continuous-age analysis |
|---|---|---|
| Participants | 251, divided into training/test/validation | 812, divided into five studies |
| Expression | Training-reference DESeq2 VST | Within-person ranks of length-adjusted counts |
| Participant pathway score | Average standardized expression; six selected genes excluded in focused comparison | Directional single-sample rank enrichment; all measured set members included |
| SD scale | Training-reference SD | Pooled SD across 812 people |
| Primary comparison | TB–nonTB level difference | Whether TB–nonTB difference changes with age |
| Adjustment in the comparison | Recruiting study + estimated RNA sources | Smooth age, sex, HIV, country and count-read length as estimable within study |
| Multiple testing | 30 focused pathway-score comparisons; GSEA FDR separately per run | 250 pathway–study tests for each primary question |

The old GSEA normalized enrichment score, the old participant mean-expression
score and the new single-sample rank score are different quantities. Their effect
sizes and adjusted p-values are not interchangeable. New curve extrema in the
comparison CSV are descriptive summaries of supported ages, not additional tests.
A GSEA output of zero reflects finite-permutation resolution, not a literal zero
probability. Neither analysis establishes biochemical activity or causation.

## Contents

- `focused_pathway_comparison.csv`: five prespecified older programs, their old
  effects/intervals/FDR and the corresponding new study tests/curve ranges.
- `all_shared_hallmark_comparison.csv`: all shared eligible Hallmarks under the
  old specification adjusting for study and estimated RNA sources and excluding
  the six XGBoost genes. This prevents selecting only concordant pathways.
- `participant_overlap_summary.csv`: aggregate verified overlap by old split and
  study. Participant identities were checked locally but are not exported here.
- `age_analysis_any_TB_hits.csv`: all new general-TB discoveries, with a flag
  indicating whether each set was eligible in the old GSEA.
- `primary_tests.tsv` and `primary_curves.tsv`: complete new test family and
  supported age-grid estimates. Missing curves are unsupported, not zero.
- `sensitivity_*`, `model_diagnostics.json`: robustness variants and explicit
  reasons for any skipped model fits.
- `pathway_timelines_all.pdf`: all eligible pathways, one per page, with two
  group curves on one age axis. Blue solid shows TB; orange dashed shows nonTB.
  Children retain the original PROTECT model. Adults use an exploratory summary
  with fixed study/clinical composition, study-specific baseline and TB levels,
  and shared smooth age shapes. These are new summary estimates, not averages of
  whichever studies happen to cover an age. The adolescent gap stays unconnected.
  Whiskers show pointwise model intervals conditional on the model assumptions;
  study-omission sensitivity is reported. At ages covered by fewer studies, the
  summary relies on the shared-shape assumption to represent the other studies.
  Faint dots show actual patient ages and observed standardized scores, before
  covariate adjustment. They remain visible where sparse support withholds curves.
- `summary_*`: separate summary-model tests, complete predictions, diagnostics,
  and leave-one-adult-study-out checks. Its FDR families differ from the original
  study-specific primary analysis, which remains the basis of the manuscript
  comparison tables. The shared-age-shape assumption is not validated: both
  full-range and supported-age formal heterogeneity checks fail the leverage
  safeguards and are explicitly unavailable. Omission checks include changes in
  the reference population and do not constitute independent validation.
- `pathway_age_explorer.html`: self-contained interactive view with one age axis,
  pathway/source selectors, optional confidence bands and a patient-dot toggle.
  Dots follow the study/omission selection and are embedded as transparent raster
  images, without patient identifiers or machine-readable patient rows. Numerical
  model estimates are aggregate results. The file requires no external scripts.
  The primary figure
  design follows two-group lifespan displays such as
  [Bethlehem et al., 2022](https://www.nature.com/articles/s41586-022-04554-y),
  while preserving the limitations of the much smaller, disjoint-age TB cohorts.
- `pathway_timeline_index.tsv`: pathway-to-page lookup, with a flag identifying
  the focus pathways shown in the notebook. Their `timeline_*.png` and `.pdf`
  figures are included individually.
- `age_pathway_analysis.html`: archived HTML report with embedded figures and
  aggregate outputs. It loads MathJax/RequireJS from a CDN for rendering and is
  not an independently rerunnable analysis package.
- `age_analysis_provenance.json`: original overview source/software/output hashes.
- `comparison_audit.json`: local input hashes, matching checks and comparison scope.
- `snapshot_manifest.json`: exact SHA-256 values of every published snapshot file.

## Rebuild

From the `cfRNA_TB` repository root, after executing the original overview and
manuscript notebooks:

```bash
.venv/bin/python Scripts_For_Manuscript/export_pathway_age_comparison.py \
  --overview ../cfRNA_TB_overview
```

The exporter reads `../cfRNA_TB_overview/results/pathway_age/`, the manuscript
`Figure_Output/Gene_Interactions/focused/pathway*.csv` results, and old/new
metadata locally to verify overlap. These source paths and hashes are recorded
without copying participant records. Use the original overview notebook and its
local data/resources to refit the age models; this exporter only copies and
compares completed results.
""")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    inferred_root = Path(__file__).resolve().parents[1]
    parser.add_argument("--root", type=Path, default=inferred_root)
    parser.add_argument(
        "--overview", type=Path, default=inferred_root.parent / "cfRNA_TB_overview"
    )
    args = parser.parse_args()
    root, overview = args.root.resolve(), args.overview.resolve()
    out = root / "Figure_Output/Gene_Interactions/pathway_age_comparison"
    source = overview / "results/pathway_age"
    out.mkdir(parents=True, exist_ok=True)
    copies = {name: name for name in AGGREGATE_FILES}
    copies.update(
        {
            f"{name}.{suffix}": f"{name}.{suffix}"
            for name in FIGURES + timeline_figure_names(source)
            for suffix in ["pdf", "png"]
        }
    )
    copies.update(
        {
            "02_pathway_enrichment_by_age.html": "age_pathway_analysis.html",
            "provenance.json": "age_analysis_provenance.json",
            "pathway_timelines_all.pdf": "pathway_timelines_all.pdf",
            "pathway_age_explorer.html": "pathway_age_explorer.html",
        }
    )
    copied_sources = {}
    for source_name, output_name in copies.items():
        path = source / source_name
        if not path.is_file():
            raise FileNotFoundError(path)
        shutil.copyfile(path, out / output_name)
        copied_sources[output_name] = {"source_path": str(path), "sha256": sha256(path)}
    audit, focus = compare_saved_results(root, overview, out)
    write_readable_comparison(out, focus)
    write_readme(out, audit)
    forbidden = {
        "participants.tsv",
        "sample_qc.tsv",
        "scores_raw.tsv",
        "scores_standardized.tsv",
    }
    if forbidden.intersection(p.name for p in out.iterdir()):
        raise ValueError("Participant-level export found in aggregate snapshot")
    files = {
        p.name: sha256(p)
        for p in sorted(out.iterdir())
        if p.is_file() and p.name != "snapshot_manifest.json"
    }
    manifest = {
        "created_utc": datetime.now(UTC).isoformat(),
        "description": "Aggregate model outputs and rendered patient distributions; no participant-level data tables or identifiers",
        "source_roots": {"manuscript": str(root), "overview": str(overview)},
        "exporter": {"path": str(Path(__file__).resolve()), "sha256": sha256(__file__)},
        "copied_sources": copied_sources,
        "comparison_inputs": audit["inputs"],
        "files": files,
    }
    (out / "snapshot_manifest.json").write_text(json.dumps(manifest, indent=2))
    for name, expected in files.items():
        assert sha256(out / name) == expected, name
    print(f"Exported {len(files)} checksummed aggregate artifacts to {out}")
    print(
        f"Verified overlap: {audit['overlap']['old_in_new']} / {audit['overlap']['old_participants']} earlier participants"
    )


if __name__ == "__main__":
    main()
