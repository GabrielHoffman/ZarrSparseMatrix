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
