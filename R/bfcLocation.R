#' Get or set the directory of the BiocFileCache used by tximeta
#'
#' Running \code{getTximetaBFC} will report the saved directory,
#' if it has been determined, or will return NULL.
#' Running \code{setTximetaBFC} will ask the user to specify a
#' BiocFileCache directory for accessing and saving TxDb sqlite files.
#'
#' @return the directory of the BiocFileCache used by tximeta
#' (or nothing, in the case of \code{setTximetaBFC})
#'
#' @rdname getTximetaBFC
#' 
#' @export
getTximetaBFC <- function() {
  bfclocFile <- bfclocFile()
  if (!file.exists(bfclocFile)) {
    message("tximeta's BiocFileCache location has not yet been set")
    NULL
  } else {
    readBFCLocFile(bfclocFile)
  }
}

#' @rdname getTximetaBFC
#' 
#' @export
setTximetaBFC <- function() {
  message("Which BiocFileCache directory should tximeta use? (press Enter to cancel)")
  bfcloc <- file.choose()  
  bfclocFile <- bfclocFile()
  writeBFCLocFile(bfcloc)
  invisible()
}

# not exported:

bfclocFile <- function() {
  tximetaDir <- user_cache_dir("tximeta")
  file.path(tximetaDir, "bfcloc.json")
}

writeBFCLocFile <- function(bfcloc) {
  tximetaDir <- user_cache_dir("tximeta")
  if (!file.exists(tximetaDir)) dir.create(tximetaDir)
  bfclocFile <- bfclocFile()
  write(toJSON(bfcloc, pretty=TRUE), file=bfclocFile)
}

readBFCLocFile <- function(bfclocFile) {
  fromJSON(bfclocFile)
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

  # this is the JSON file where we store the location of the tximeta BiocFileCache
  bfclocFile <- bfclocFile()
  
  # this file tells us which BFC dir has been previously chosen use with tximeta
  if (!file.exists(bfclocFile)) {
    if (interactive()) {

      # exception: temp dir already specified as BFC location in this session
      if (file.exists(file.path(tempdir(),"BiocFileCache.sqlite"))) {
        bfcloc <- tempdir()
        return(bfcloc)
      }
      
      # otherwise ask user:
      ans <- menu(c("Yes (use default)", "No (use temp)"), title=prompt)
      if (ans == 0) stop("no BiocFileCache directory choice made at this time")

      # user wants to use default dir:
      if (ans == 1) {
        bfcloc <- defaultDir
        writeBFCLocFile(bfcloc)
        
        # user wants to use temp dir:
      } else if (ans == 2) {
        bfcloc <- tempdir()
      }
      
      # not interactive, use temp dir:
    } else {
      bfcloc <- tempdir()
    }
    
    # file already exists, read BFC loc:
  } else {
    bfcloc <- readBFCLocFile(bfclocFile)
  }
  
  bfcloc
}
