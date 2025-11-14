# Get or set the directory of the BiocFileCache used by tximeta

Running `getTximetaBFC` will report the saved directory, if it has been
determined, or will return NULL. Running `setTximetaBFC` will ask the
user to specify a BiocFileCache directory for accessing and saving TxDb
sqlite files. Note that tximeta's BiocFileCache can be set by the
environmental variable `TXIMETA_HUB_CACHE`, which will reset the cache
location.

## Usage

``` r
getTximetaBFC()

setTximetaBFC(dir, quiet = FALSE)
```

## Arguments

- dir:

  the location for tximeta's BiocFileCache. can be missing in which case
  the function will call `file.choose` for choosing location
  interactively

- quiet:

  whether to suppress feedback message

## Value

the directory of the BiocFileCache used by tximeta (or nothing, in the
case of `setTximetaBFC`)

## Examples

``` r

# getting the BiocFileCache used by tximeta
# (may not be set, which uses BiocFileCache default or temp directory)
getTximetaBFC()
#> tximeta's BiocFileCache location has not yet been set
#> NULL

# don't want to actually change user settings so this is not run:
# setTximetaBFC()
```
