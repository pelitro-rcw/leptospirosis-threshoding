# Stage C RF threshold derivation: primary A1/A2 qualifying trigger

For each evaluable region-year and each disease detector, Stage C first enumerates every detector-active trigger week.

The positive threshold-derivation anchor is the **first detector-active trigger week that falls within A1 or A2**. This yields at most one primary true anchor per region-year. Detector-active trigger weeks outside both A1 and A2 are retained as false-alarm comparators. Later detector-active weeks that also fall in A1/A2 are retained in the audit inventory but are excluded from threshold fitting so a prolonged detector episode cannot contribute repeated positive observations from the same region-year.

The rule is applied identically to both rainfall representations:

- P1 and P3 use the Constant TA trigger inventory.
- P2 and P4 use the Outbreak Threshold trigger inventory.
- P1/P2 use cumulative RF_HDX during t-4 through t-1.
- P3/P4 use the selected rainfall STA/LTA R(t) feature during t-4 through t-1.
- The disease-trigger week is excluded from the RF feature.

Stage C exports the complete active-trigger audit inventory, the first A1/A2-true anchor inventory, the threshold-fitting inventory, and an explicit support table distinguishing disease-anchor qualifying years from RF-eligible qualifying years.

Stage C figures are exported separately for NCR and Regional scopes. No mixed-scope threshold-product, STA/LTA-window, operational-bound, or Stage C head-to-head figure is produced.
