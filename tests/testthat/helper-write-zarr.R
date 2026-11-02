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

### Write integer-valued doubles 'v' as an uncompressed 1D Zarr v2 int64
### ("<i8") or uint64 ("<u8") array, as zarr-python would. The bytes of each
### value are written directly as its low and high 32-bit words.
write_int64_zarr <- function(v, path, chunk, dtype="<i8")
{
    dir.create(path, recursive=TRUE)
    writeLines(sprintf(paste0('{"zarr_format": 2, "shape": [%d], ',
                              '"chunks": [%d], "dtype": "%s", ',
                              '"compressor": null, "fill_value": 0, ',
                              '"filters": null, "order": "C"}'),
                       length(v), chunk, dtype),
               file.path(path, ".zarray"))
    pad <- (-length(v)) %% chunk
    v <- c(v, rep(0, pad))
    hi <- floor(v / 2^32)
    lo <- v - hi * 2^32
    ## A word of 0x80000000 becomes NA_integer_, which has that bit pattern,
    ## so the coercion warning can be ignored.
    to_int32 <- function(w)
        suppressWarnings(as.integer(ifelse(w >= 2^31, w - 2^32, w)))
    words <- as.vector(rbind(to_int32(lo), to_int32(hi %% 2^32)))
    for (k in seq_len(length(v) %/% chunk)) {
        w <- words[seq(2L * (k - 1L) * chunk + 1L, length.out=2L * chunk)]
        writeBin(w, file.path(path, k - 1L), size=4L, endian="little")
    }
}
