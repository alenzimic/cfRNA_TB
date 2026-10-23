"""Regenerate the fixed training folds with the notebook's pinned Python environment."""

from pathlib import Path

import pandas as pd
from sklearn.model_selection import StratifiedKFold


reference_dir = Path(__file__).resolve().parent
metadata = pd.read_csv(reference_dir.parent / "metadata.train.tsv", sep="\t")
folds = metadata[["sample_id", "tb_binary"]].copy()
folds["validation_fold"] = -1
splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=87)
for fold, (_, held_out) in enumerate(splitter.split(metadata, metadata["tb_binary"])):
    folds.loc[held_out, "validation_fold"] = fold
folds.to_csv(reference_dir / "folds.csv", index=False)
