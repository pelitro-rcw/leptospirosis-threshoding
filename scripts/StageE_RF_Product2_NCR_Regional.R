# Full NCR + Regional rerun using RF Product 2 as an antecedent gate.
.run_rf_product2 <- function() {
  old_p <- Sys.getenv("TA_RF_GATE_PRODUCT", unset="")
  old_s <- Sys.getenv("TA_OUTPUT_SUFFIX", unset="")
  on.exit(Sys.setenv(TA_RF_GATE_PRODUCT=old_p, TA_OUTPUT_SUFFIX=old_s), add=TRUE)
  Sys.setenv(TA_RF_GATE_PRODUCT="2", TA_OUTPUT_SUFFIX="_RF_Product2")
  sys.source(file.path("scripts","StageA_NCR_Lepto_Analysis.R"), envir=new.env(parent=globalenv()))
  sys.source(file.path("scripts","StageB_Regional_Lepto_Analysis.R"), envir=new.env(parent=globalenv()))
}
.run_rf_product2()
