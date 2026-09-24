# Legend Layout QA

This revision standardizes figure legends to prevent clipping, overlap, excessive lateral dead space, and visually disoriented key placement.

## Global rule

All figures based on the shared publication theme use bottom-positioned legends with horizontal reading order, centered alignment, compact key widths, compact x/y spacing, and additional bottom plot margin.

## Long legends

- NCR outbreak-threshold drift: the three requested entries wrap horizontally rather than extending beyond the canvas.
- NCR and Regional RF products: environmental series, adopted RF threshold, RF trigger, disease anchor, and the outbreak-threshold key when applicable are constrained to compact horizontal rows.
- Regional detector map: up to six detector keys wrap into two horizontal rows.
- NCR 2019/2024 comparison: the previous 26-point inter-item spacing was removed and replaced with compact publication spacing.
- Multipanel NCR figures retain a shared collected legend where applicable.

## Analytical scope

This is a display-layout revision only. It does not alter detector calculations, RF thresholds, RF operators, region inclusion, trigger weeks, metrics, or exported numerical tables. The no-RF leptospirosis figures retain their disease-only analytical plotting structure.


## Deterministic shared-legend correction

- NCR multipanel and 2019/2024 comparison figures no longer use data-dependent `patchwork` guide collection.
- One complete legend is extracted once with `cowplot::get_legend()` and placed below a legend-free panel grid.
- Trigger keys are trained independently of whether a trigger occurs in a particular year, preventing missing or duplicated entries.
- The disease-only timing legend is organized as two compact horizontal rows: line-series keys and detector-trigger keys.
- The regional detector map always presents all six detector categories in a fixed two-row order.
