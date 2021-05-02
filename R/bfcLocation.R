#' Get or set the directory of the BiocFileCache used by tximeta
#'
#' Running \code{getTximetaBFC} will report the saved directory,
#' if it has been determined, or will return NULL.
#' Running \code{setTximetaBFC} will ask the user to specify a
#' BiocFileCache directory for accessing and saving TxDb sqlite files.
#' Note that tximeta's BiocFileCache can be set by the environmental
#' variable \code{TXIMETA_HUB_CACHE}, which will reset the cache location.
#' 
#' @param dir the location for tximeta's BiocFileCache. can be missing
#' in which case the function will call \code{file.choose} for choosing
#' location interactively
#' @param quiet whether to suppress feedback message
#' 
#' @return the directory of the BiocFileCache used by tximeta
#' (or nothing, in the case of \code{setTximetaBFC})
#'
#' @rdname getTximetaBFC
#'
#' @examples
#'
#' # getting the BiocFileCache used by tximeta
#' # (may not be set, which uses BiocFileCache default or temp directory)
#' getTximetaBFC()
#'
#' # don't want to actually change user settings so this is not run:
#' # setTximetaBFC()
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
setTximetaBFC <- function(dir, quiet=FALSE) {
  if (missing(dir)) {
    message("which BiocFileCache directory should tximeta use? (press Enter to cancel)")
    bfcloc <- file.choose()
  } else {
    stopifnot(is(dir, "character"))
    bfcloc <- dir
  }
  bfclocFile <- bfclocFile()
  writeBFCLocFile(bfcloc)
  if (!quiet) message("for group use, set the permissions of this directory to allow group write (g+w)")
  invisible()
}

# not exported:

# functions to read or edit the BiocFileCache location for tximeta (bfcloc)
# this information is stored in the 'tximeta' default R_user_dir
# under a file 'bfcloc.json'
bfclocFile <- function() {
  tximetaDir <- R_user_dir("tximeta", which="config")
  file.path(tximetaDir, "bfcloc.json")
}

writeBFCLocFile <- function(bfcloc) {
  tximetaDir <- R_user_dir("tximeta", which="config")
  if (!file.exists(tximetaDir)) dir.create(tximetaDir, recursive=TRUE)
  bfclocFile <- bfclocFile()
  write(toJSON(bfcloc, pretty=TRUE), file=bfclocFile)
}

readBFCLocFile <- function(bfclocFile) {
  fromJSON(bfclocFile)
}

# an internal function for getting the BFC location
# the logic here depends on whether a location has been set before,
# and whether tximport is being run interactively or not
getBFCLoc <- function() {
  defaultDir <- R_user_dir("BiocFileCache", which="cache")

  prompt <- paste("",
  "tximeta needs a BiocFileCache directory to access and save TxDb objects.",
  paste0("Do you wish to use the default directory: '",defaultDir,"'?"),
  "If not, a temporary directory that is specific to this R session will be used.","",
  "You can always change this directory later by running: setTximetaBFC()",
  "Or enter [0] to exit and set this directory manually now.",
  "This location can also be set by environmental variable TXIMETA_HUB_CACHE.",
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
