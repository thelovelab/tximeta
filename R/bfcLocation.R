#' @rdname getTximetaBFC
#' 
#' @export
setTximetaBFC <- function() {
  message("Which BiocFileCache directory should tximeta use? (press Enter to cancel)")
  bfcloc <- file.choose()
  tximetaDir <- user_cache_dir("tximeta")
  bfclocFile <- file.path(tximetaDir, "bfcloc.txt")
  if (!file.exists(tximetaDir)) dir.create(tximetaDir)
  write(bfcloc, file=bfclocFile)
  bfcloc
}

#' Get or set the directory of the BiocFileCache used by tximeta
#'
#' Running \code{getTximetaBFC} will report the saved directory,
#' if it has been determined, or will return NULL.
#' Running \code{setTximetaBFC} will ask the user to specify a
#' BiocFileCache directory for accessing and saving TxDb sqlite files.
#'
#' @return the directory of the BiocFileCache used by tximeta
#'
#' @rdname getTximetaBFC
#' 
#' @export
getTximetaBFC <- function() {
  tximetaDir <- user_cache_dir("tximeta")
  bfclocFile <- file.path(tximetaDir, "bfcloc.txt")
  if (!file.exists(bfclocFile)) {
    message("tximeta's BiocFileCache location has not yet been set")
    NULL
  } else {
    scan(bfclocFile, what="char", quiet=TRUE)
  }
}

# an internal function for getting the BFC location
# (above functions are user-facing for getting/setting)
getBFCLoc <- function() {
  defaultDir <- user_cache_dir(appname="BiocFileCache")

  prompt <- paste("",
  "tximeta needs a BiocFileCache directory to access and save TxDb objects.",
  paste0("Do you wish to use the default directory: '",defaultDir,"'?"),
  "If not, a temporary directory that is specific to this R session will be used.","",
  "You can always change this directory later by running: setTximetaBFC()",
  "Or enter [0] to exit and set this directory manually now.",
  sep="\n")
  
  # this file tells us which BFC dir has been previously chosen use with tximeta
  tximetaDir <- user_cache_dir("tximeta")
  bfclocFile <- file.path(tximetaDir, "bfcloc.txt")
  if (!file.exists(bfclocFile)) {
    if (interactive()) {
      
      # exception: temp dir already specified as BFC location in this session
      if (file.exists(file.path(tempdir(),"BiocFileCache.sqlite"))) {
        bfcloc <- tempdir()
        return(bfcloc)
      }
      
      ans <- menu(c("Yes (use default)", "No (use temp)"), title=prompt)
      if (ans == 0) stop("no BiocFileCache directory choice made at this time")
      if (ans == 1) {
        bfcloc <- defaultDir
        if (!file.exists(tximetaDir)) dir.create(tximetaDir)
        write(bfcloc, file=bfclocFile)
      } else {
        bfcloc <- tempdir()
      }
    } else {
      bfcloc <- tempdir()
    }
  } else {
    bfcloc <- scan(bfclocFile, what="char", quiet=TRUE)
  }
  bfcloc
}
