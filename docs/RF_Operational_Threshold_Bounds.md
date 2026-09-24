# RF operational threshold point estimates and stability bounds

All four rainfall products use a prospective upper-tail activation rule: **RF >= threshold**. The RF feature uses t-4..t-1 and excludes the disease-trigger week. A1/A2 remains the independent true/false standard: an RF-supported Constant TA or Outbreak Threshold anchor is true iff its trigger week lies in A1 OR A2; otherwise it is false. Products 1-2 use rolling four-week RF_HDX accumulation; Products 3-4 use the empirically selected antecedent rainfall STA/LTA R(t) feature.

## Point estimate

For each product and spatial scope, candidate upper-tail thresholds are evaluated using sensitivity, specificity, balanced accuracy, PPV, trigger rate, false-alarm rate, and false-alarm weeks per evaluable unit. Thresholds within 0.01 of the maximum Youden J are retained; the adopted point estimate is selected by the lowest false-alarm rate, then lower false-alarm burden per evaluable unit, higher PPV, sensitivity, specificity, and finally the higher threshold. The complete candidate search is exported as `RF_Operational_Threshold_Candidate_Grid.csv`.

## Scope-specific fixed thresholds

NCR thresholds are derived only from NCR historical observations and remain fixed across NCR evaluation years. Their stability bounds are obtained from leave-one-year-out re-derivation. Regional thresholds are derived from the pooled 17-region dataset and remain fixed across all 17 regions; their stability bounds are obtained from leave-one-region-out re-derivation. For Products 3-4, each held-out fold repeats the rainfall STA/LTA candidate-window selection before deriving that fold's threshold, so the bounds include both operating-threshold and RF-window selection instability.

## Bounds and activation status

For every product, the lower and upper values are empirical stability bounds around the fixed point estimate, with `lower <= point <= upper`. They are not annual thresholds. Operational status is monotonic in rainfall:

- RF < lower: No RF activation
- lower <= RF < point: Precautionary RF trigger
- point <= RF < upper: Formal RF trigger
- RF >= upper: Strong RF trigger

The point estimate is the threshold used for headline RF-gated metrics. The lower/upper analyses are sensitivity and operational-uncertainty analyses.
