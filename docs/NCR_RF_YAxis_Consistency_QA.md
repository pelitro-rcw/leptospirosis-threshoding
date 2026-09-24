# NCR RF Weekly-Case Y-Axis Consistency QA

The NCR RF-product timing figures and the corresponding disease-only NCR timing figures use the same left-axis weekly leptospirosis scale.

- `Y_MAX` is derived once from the full NCR weekly leptospirosis series.
- `NCR_WEEKLY_CASE_AXIS_MAX <- Y_MAX * 1.10` is the common visible upper limit.
- Disease-only NCR panels use `ylim = c(0, NCR_WEEKLY_CASE_AXIS_MAX)`.
- RF Product 1-4 NCR panels set `case_upper <- NCR_WEEKLY_CASE_AXIS_MAX` and use the same 0-to-upper-limit range.
- RF products retain their own secondary environmental axis: rainfall mm for Products 1/2 and rainfall STA/LTA R(t) for Products 3/4.
- No analytical values are truncated or rescaled in the underlying tables or trigger calculations; this is a display-scale consistency rule only.
