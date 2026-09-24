# RF bounds QA

The current pipeline enforces the following contracts for all Products 1-4:

1. Operational direction is `>=` only.
2. Low RF below the lower stability bound is never labelled a strong trigger.
3. Each scope/product has one finite ordered threshold triplet: lower <= point <= upper.
4. NCR thresholds are NCR-derived and use LOYO stability bounds.
5. Regional thresholds are pooled-17-region-derived and use LORO stability bounds.
6. Candidate point estimates are selected using upper-tail discrimination with explicit false-alarm-rate tie-breaking.
7. The full candidate threshold grid is exported for reproducibility.
8. Point estimates remain the primary formal trigger; bound-based classifications are operational sensitivity information.
