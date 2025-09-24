context("alevin")

test_that("tximeta can import alevin", {

  dir <- system.file("extdata", package="tximportData")
  samps <- list.files(file.path(dir, "alevin"))
  files <- file.path(dir,"alevin",samps[1],"alevin/quants_mat.gz")
  file.exists(files)
  coldata <- data.frame(files, names="neurons")

  #se <- tximeta(coldata, type="alevin")
  se <- tximeta(coldata, type="alevin", skipMeta=TRUE)
  expect_true(metadata(se)$tximetaInfo$type == "alevin")
  
})
