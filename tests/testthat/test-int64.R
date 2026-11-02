dir <- tempfile()
dir.create(dir)

test_that("int64 values beyond the R integer range are read exactly", {
    v <- c(0, 1, 5, 2^31 - 1, 2^31, 2^32 - 1, 2^32, 2^32 + 7, 3e12,
           2^53 - 1, -1, -2^31, -2^31 - 1, -2^40)
    path <- file.path(dir, "i8")
    write_int64_zarr(v, path, 4L)
    read_1d <- ZarrSparseMatrix:::.read_zarr_1d
    expect_identical(read_1d(path), v)
    expect_identical(read_1d(path, c(10L, 2L, 5L, 5L)), v[c(10L, 2L, 5L, 5L)])
    expect_identical(read_1d(path, integer(0)), double(0))
})

test_that("uint64 values are read exactly", {
    v <- c(0, 2^31, 2^32 + 1, 2^53 - 1)
    path <- file.path(dir, "u8")
    write_int64_zarr(v, path, 3L, dtype="<u8")
    expect_identical(ZarrSparseMatrix:::.read_zarr_1d(path), v)
})

test_that("int64 values beyond 2^53 are an error", {
    path <- file.path(dir, "big")
    write_int64_zarr(2^60, path, 1L)
    expect_error(ZarrSparseMatrix:::.read_zarr_1d(path), "2\\^53")
})

test_that("a matrix with an int64 indptr can be read", {
    set.seed(5)
    M <- Matrix::rsparsematrix(30L, 20L, density=0.2)
    for (layout in c("csr", "csc")) {
        store <- file.path(tempfile(), "s.zarr")
        write_sparse_zarr(store, M, layout)
        indptr <- file.path(store, "X", "indptr")
        p <- as.numeric(ZarrSparseMatrix:::.read_zarr_1d(indptr))
        unlink(indptr, recursive=TRUE)
        write_int64_zarr(p, indptr, 5L)
        X <- ZarrSparseMatrix(store, "X")
        expect_equal(as.matrix(X), as.matrix(M), ignore_attr=TRUE)
        expect_equal(as.matrix(X[c(3, 1), c(2, 20)]),
                     as.matrix(M[c(3, 1), c(2, 20)]), ignore_attr=TRUE)
    }
})

test_that("positions beyond 2^31 - 1 are computed exactly", {
    pos <- ZarrSparseMatrix:::.ranges_to_positions(c(1, 2^31 + 5, 3e9),
                                                    c(2L, 3L, 0L))
    expect_identical(pos, c(1, 2, 2^31 + 5, 2^31 + 6, 2^31 + 7))
})

test_that("columns stored beyond position 2^31 - 1 can be read", {
    ## A (2^31 - 1) x 2 CSC matrix: column 1 is full, so column 2 starts at
    ## position 2^31 of 'data' and 'indices'. Only the chunk holding column 2
    ## exists on disk (Zarr treats missing chunks as fill values), and only
    ## column 2 is read.
    nrow <- 2^31 - 1
    N <- nrow + 2
    chunk <- 50
    store <- file.path(tempfile(), "big.zarr")
    g <- file.path(store, "X")
    dir.create(g, recursive=TRUE)
    write_last_chunk <- function(name, dtype, values, size) {
        path <- file.path(g, name)
        dir.create(path)
        writeLines(sprintf(paste0('{"zarr_format": 2, "shape": [%s], ',
                                  '"chunks": [%d], "dtype": "%s", ',
                                  '"compressor": null, "fill_value": 0, ',
                                  '"filters": null, "order": "C"}'),
                           format(N, scientific=FALSE), chunk, dtype),
                   file.path(path, ".zarray"))
        ## 'values' go at positions N - 1 and N, in the last chunk, which
        ## is padded to full size.
        cid <- (N - 1) %/% chunk
        before <- N - cid * chunk - length(values)
        zero <- vector(typeof(values), 1L)
        v <- c(rep(zero, before), values,
               rep(zero, chunk - before - length(values)))
        writeBin(v, file.path(path, format(cid, scientific=FALSE)),
                 size=size, endian="little")
    }
    writeLines('{"zarr_format": 2}', file.path(store, ".zgroup"))
    writeLines('{"zarr_format": 2}', file.path(g, ".zgroup"))
    writeLines(sprintf(paste0('{"encoding-type": "csc_matrix", ',
                              '"encoding-version": "0.1.0", ',
                              '"shape": [%s, 2]}'),
                       format(nrow, scientific=FALSE)),
               file.path(g, ".zattrs"))
    write_last_chunk("data", "<f8", c(1.5, -2), 8L)
    write_last_chunk("indices", "<i4", c(0L, 2L), 4L)
    write_int64_zarr(c(0, nrow, N), file.path(g, "indptr"), 3L)

    seed <- ZarrSparseMatrixSeed(store, "X")
    expect_identical(nzcount(seed), N)
    trip <- ZarrSparseMatrix:::.load_CSC_triplets(seed, j=2L)
    expect_identical(trip$dim, c(as.integer(nrow), 1L))
    expect_identical(trip$row, c(1L, 3L))
    expect_identical(trip$val, c(1.5, -2))
    trip <- ZarrSparseMatrix:::.load_CSC_triplets(seed, i=c(3L, 2L), j=2L)
    expect_identical(trip$row, 1L)
    expect_identical(trip$val, -2)
})
