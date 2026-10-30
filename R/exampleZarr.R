#' Path to an example Zarr store
#'
#' Returns the path to a small example Zarr store containing simulated counts
#' for 200 cells and 100 genes, stored like the \code{X} group of an AnnData
#' \code{.zarr} store (CSR layout, cells as rows).
#'
#' The store is shipped as a zip archive, because Zarr v2 metadata files
#' (\code{.zarray}, \code{.zattrs}, \code{.zgroup}) are hidden files, which
#' are not allowed in R packages. The first call unzips it to
#' \code{tempdir()}, and later calls in the same session reuse that copy.
#'
#' The script that creates the store is
#' \code{system.file("scripts", "make-example-data.R",
#' package="ZarrSparseMatrix")}.
#'
#' @return The path to the \code{example.zarr} directory.
#'
#' @seealso \link{ZarrSparseMatrix}
#'
#' @examples
#' store <- exampleZarr()
#' list.files(store, recursive=TRUE, all.files=TRUE)
#'
#' X <- ZarrSparseMatrix(store, "X")
#' X
#'
#' @importFrom utils unzip
#' @export
exampleZarr <- function()
{
    exdir <- file.path(tempdir(), "ZarrSparseMatrix-example")
    store <- file.path(exdir, "example.zarr")
    if (!dir.exists(store)) {
        zipfile <- system.file("extdata", "example.zarr.zip",
                               package="ZarrSparseMatrix", mustWork=TRUE)
        unzip(zipfile, exdir=exdir)
    }
    store
}
