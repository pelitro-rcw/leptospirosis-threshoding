# Independent per-region and RF-to-anchor revision

## Independent per-region analysis

Stage I runs the NCR-style analytical pipeline independently in all 17 canonical Philippine regions. This is not the pooled Regional analysis. For each region, Stage I produces a disease-only run plus RF Products 1-4.

Stage C additionally derives four independent RF thresholds for each region:

- P1: RF mm preceding Constant TA
- P2: RF mm preceding Outbreak Threshold
- P3: RF STA/LTA R(t) preceding Constant TA
- P4: RF STA/LTA R(t) preceding Outbreak Threshold

The per-region thresholds use the same Stage C first-A1/A2-qualifying trigger rule and A1/A2 false-alarm comparator definition used in the NCR analysis. Stability bounds are obtained by region-specific leave-one-year-out re-derivation.

## Pooled Regional RF products

The pooled Regional RF reruns now expose a dedicated RF-to-anchor layer. For each RF product the pipeline exports all target anchor trigger weeks, their antecedent RF feature, threshold status, and A1/A2 true/false alarm classification. A 17-region summary and geographic map report the performance of the RF qualification of the intended disease anchor.

P1/P3 are always tied to Constant TA. P2/P4 are always tied to Outbreak Threshold. Figure B5 six-detector products remain available as diagnostic context and are not the primary RF-product visualization.
