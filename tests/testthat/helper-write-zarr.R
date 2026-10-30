### Write a 1D Zarr array the way zarr-python does, i.e. with the last chunk
### padded to full size. Rarr::write_zarr_array() writes a short last chunk,
### which Rarr::read_zarr_array() then fails to read with an 'index'.
write_1d_zarr <- function(v, path, chunk)
{
    pad <- (-length(v)) %% chunk
    Rarr::write_zarr_array(array(c(v, rep(v[1L], pad))), path,
                           chunk_dim=chunk)
    f <- file.path(path, ".zarray")
    meta <- jsonlite::read_json(f)
    meta$shape <- list(length(v))
    jsonlite::write_json(meta, f, auto_unbox=TRUE)
}

### Write sparse matrix 'M' to group "X" of a new Zarr v2 store, in AnnData
### layout ("csr" or "csc").
write_sparse_zarr <- function(store, M, layout=c("csr", "csc"))
{
    layout <- match.arg(layout)
    g <- file.path(store, "X")
    dir.create(g, recursive=TRUE)
    writeLines('{"zarr_format": 2}', file.path(store, ".zgroup"))
    writeLines('{"zarr_format": 2}', file.path(g, ".zgroup"))
    writeLines(sprintf(paste0('{"encoding-type": "%s_matrix", ',
                              '"encoding-version": "0.1.0", ',
                              '"shape": [%d, %d]}'),
                       layout, nrow(M), ncol(M)),
               file.path(g, ".zattrs"))
    if (layout == "csr") {
        S <- as(Matrix::t(M), "CsparseMatrix")
    } else {
        S <- as(M, "CsparseMatrix")
    }
    write_1d_zarr(S@x, file.path(g, "data"), 7L)
    write_1d_zarr(S@i, file.path(g, "indices"), 7L)
    write_1d_zarr(S@p, file.path(g, "indptr"), 5L)
    invisible(store)
}
