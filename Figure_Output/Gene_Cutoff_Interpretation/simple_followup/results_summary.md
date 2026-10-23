**1. Cutoff sensitivity.** Mean fold AUC across the four DE cutoffs spans **0.812–0.876**. At 0.05 it is **0.812**, compared with **0.881** for abundance only. These are descriptive comparisons from five training folds; no winning cutoff is selected.

**2. Filter stability.** **42/142** full-training DE genes pass the 0.05 filter in all five fitting subsets. The published genes are retained as follows: GBP5 5/5, BNIP3L 5/5, KLF6 2/5, DYSF 5/5, LASP1 5/5, PCBP1 3/5. Inclusion is not predictive importance.

**3. Gene removal.** The largest observed AUC loss is **+0.107 after removing GBP5** on test and **+0.129 after removing GBP5** on validation. The same gene has the largest point-estimate loss in both cohorts. Read the paired intervals before interpreting these point estimates; the models and six-gene list are already fixed.

**Scope:** these analyses explain sensitivity of the available models and data. They do not recover the missing 113-gene input pool, establish a superior classifier, or demonstrate causal gene interactions.
