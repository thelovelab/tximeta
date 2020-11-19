.onLoad <- function(libname, pkgname, ...) {

  envBFCLoc <- Sys.getenv("TXIMETA_HUB_CACHE")
  if (envBFCLoc != "") {
    bfcLoc <- getTximetaBFC()
    if (is.null(bfcLoc) || envBFCLoc != bfcLoc) {
      setTximetaBFC(envBFCLoc, quiet=TRUE)
    }
  }
  
}
