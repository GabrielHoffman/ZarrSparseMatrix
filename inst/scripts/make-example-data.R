### Create inst/extdata/example.zarr: a small simulated counts matrix
### (200 cells x 100 genes) stored like the X group of an AnnData .zarr
### store, i.e. in CSR layout with cells as rows.
###
### Run from the package root:
###   Rscript inst/scripts/make-example-data.R

library(Matrix)
library(Rarr)
library(jsonlite)

### Write a 1D Zarr array the way zarr-python does, i.e. with the last chunk
### padded to full size (Rarr::write_zarr_array() writes a short last chunk).
write_1d_zarr <- function(v, path, chunk)
{
    pad <- (-length(v)) %% chunk
    write_zarr_array(array(c(v, rep(v[1L], pad))), path, chunk_dim=chunk)
    f <- file.path(path, ".zarray")
    meta <- read_json(f)
    meta$shape <- list(length(v))
    write_json(meta, f, auto_unbox=TRUE)
}

set.seed(42)
n_cells <- 200L
n_genes <- 100L
mu <- rlnorm(n_genes, meanlog=-1, sdlog=1.5)
counts <- matrix(rnbinom(n_cells * n_genes, mu=rep(mu, each=n_cells),
                         size=0.5),
                 n_cells, n_genes)
counts <- as(counts, "CsparseMatrix")   # cells x genes

store <- file.path("inst", "extdata", "example.zarr")
unlink(store, recursive=TRUE)
g <- file.path(store, "X")
dir.create(g, recursive=TRUE)
writeLines('{"zarr_format": 2}', file.path(store, ".zgroup"))
writeLines('{"zarr_format": 2}', file.path(g, ".zgroup"))
write_json(list(`encoding-type`="csr_matrix", `encoding-version`="0.1.0",
                shape=c(n_cells, n_genes)),
           file.path(g, ".zattrs"), auto_unbox=TRUE)

### CSR of 'counts' == CSC of its transpose.
S <- as(t(counts), "CsparseMatrix")
write_1d_zarr(as.integer(S@x), file.path(g, "data"), 1000L)
write_1d_zarr(S@i, file.path(g, "indices"), 1000L)
write_1d_zarr(S@p, file.path(g, "indptr"), 1000L)
