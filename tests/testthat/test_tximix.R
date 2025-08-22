context("tximix")
test_that("tximix works as expected", {

  dir <- system.file("extdata/oarfish", package="tximportData")
  names <- paste0("rep", 2:4)
  files <- file.path(dir, paste0("sgnex_h9_", names, ".quant.gz"))
  coldata <- data.frame(files, names)
  se <- tximeta(coldata, type="oarfish", skipMeta=TRUE)

  se <- tximix(coldata, type="oarfish")
  
})
