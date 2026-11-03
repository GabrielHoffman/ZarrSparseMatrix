set.seed(1)
M <- Matrix::rsparsematrix(30, 20, 0.2)
M@x <- round(M@x * 10)
D <- as.matrix(M)

for (layout in c("csr", "csc")) {

    store <- file.path(tempfile(), paste0(layout, ".zarr"))
    write_sparse_zarr(store, M, layout)
    X <- ZarrSparseMatrix(store, "X")
    seed_class <- paste0(toupper(layout), "_ZarrSparseMatrixSeed")

    test_that(paste("constructor works on", layout), {
        expect_s4_class(X, "ZarrSparseMatrix")
        expect_s4_class(DelayedArray::seed(X), seed_class)
        expect_identical(dim(X), dim(M))
        expect_true(DelayedArray::is_sparse(X))
        expect_equal(nzcount(DelayedArray::seed(X)), length(M@x))
        expect_identical(ZarrSparseMatrix(DelayedArray::seed(X)), X)
    })

    test_that(paste("full and subset reads match on", layout), {
        expect_equal(as.matrix(X), D, ignore_attr=TRUE)
        i <- c(5, 2, 2, 30)
        j <- c(20, 1, 3)
        expect_equal(as.matrix(X[i, j]), D[i, j], ignore_attr=TRUE)
        expect_equal(as.matrix(t(X)[1:3, ]), t(D)[1:3, ], ignore_attr=TRUE)
        expect_equal(as.matrix(X[integer(0), 1:2]),
                     D[integer(0), 1:2, drop=FALSE], ignore_attr=TRUE)
    })

    test_that(paste("summaries and coercion match on", layout), {
        expect_equal(as.vector(DelayedArray::colSums(X)), colSums(D))
        expect_equal(as.vector(DelayedArray::rowSums(X)), rowSums(D))
        expect_equal(as.matrix(as(X, "dgCMatrix")), D, ignore_attr=TRUE)
    })
}

test_that("sparse.layout and dim arguments override the attributes", {
    store <- file.path(tempfile(), "csc.zarr")
    write_sparse_zarr(store, M, "csc")
    ## Reading CSC data as CSR gives the transpose.
    Xt <- ZarrSparseMatrix(store, "X", dim=rev(dim(M)), sparse.layout="CSR")
    expect_equal(as.matrix(Xt), t(D), ignore_attr=TRUE)
})

test_that("dimnames are propagated", {
    store <- file.path(tempfile(), "csr.zarr")
    write_sparse_zarr(store, M, "csr")
    dn <- list(paste0("r", seq_len(nrow(M))), paste0("c", seq_len(ncol(M))))
    X <- ZarrSparseMatrix(store, "X", dimnames=dn)
    expect_identical(dimnames(X), dn)
    expect_identical(rownames(X[3:4, ]), dn[[1L]][3:4])
})

test_that("bad input is rejected", {
    store <- file.path(tempfile(), "csr.zarr")
    write_sparse_zarr(store, M, "csr")
    expect_error(ZarrSparseMatrix(store, "Y"))
    expect_error(ZarrSparseMatrix(store, "X", dim=c(5L, 5L)), "indptr")
})

test_that("row and column selections match the in-memory matrix", {
    set.seed(3)
    for (layout in c("csr", "csc")) {
        store <- file.path(tempfile(), paste0(layout, ".zarr"))
        write_sparse_zarr(store, M, layout)
        X <- ZarrSparseMatrix(store, "X")
        seed <- DelayedArray::seed(X)
        idx_list <- list(
            list(c(30L, 1L, 15L), NULL),           # unsorted rows
            list(c(4L, 4L, 9L, 4L), c(2L, 2L)),    # duplicated rows and cols
            list(integer(0), 1:5),                 # no rows
            list(1:5, integer(0)),                 # no columns
            list(NULL, c(20L, 1L)),                # all rows
            list(sample(30L, 12L), sample(20L, 7L)))
        for (index in idx_list) {
            i <- if (is.null(index[[1L]])) seq_len(nrow(D)) else index[[1L]]
            j <- if (is.null(index[[2L]])) seq_len(ncol(D)) else index[[2L]]
            a <- extract_array(seed, index)
            expect_equal(a, D[i, j, drop=FALSE], ignore_attr=TRUE)
            sa <- extract_sparse_array(seed, index)
            expect_identical(typeof(a), type(sa))
            expect_identical(dim(sa), c(length(i), length(j)))
            expect_equal(as.matrix(sa), D[i, j, drop=FALSE],
                         ignore_attr=TRUE)
        }
    }
})

test_that("t() transposes the seed instead of delaying the transpose", {
    dn <- list(paste0("r", seq_len(nrow(D))), paste0("c", seq_len(ncol(D))))
    for (layout in c("csr", "csc")) {
        store <- file.path(tempfile(), paste0(layout, ".zarr"))
        write_sparse_zarr(store, M, layout)
        X <- ZarrSparseMatrix(store, "X", dimnames=dn)
        tX <- t(X)
        expect_s4_class(tX, "ZarrSparseMatrix")
        expect_s4_class(seed(tX), if (layout == "csr") "CSC_ZarrSparseMatrixSeed"
                                  else "CSR_ZarrSparseMatrixSeed")
        expect_identical(dimnames(tX), rev(dn))
        expect_equal(as.matrix(tX), t(D), ignore_attr=TRUE)
        expect_equal(as.matrix(tX[c(5, 2, 2), c(30, 1)]),
                     t(D)[c(5, 2, 2), c(30, 1)], ignore_attr=TRUE)
        expect_identical(class(seed(t(tX))), class(seed(X)))
        ## After subsetting, t() is delayed as for any DelayedMatrix.
        Y <- t(X[1:4, ])
        expect_false(is(Y, "ZarrSparseMatrix"))
        expect_equal(as.matrix(Y), t(D[1:4, ]), ignore_attr=TRUE)
    }
})

test_that("consecutive selections with empty rows and columns are read", {
    E <- M
    E[, c(1L, 2L, 7L, 20L)] <- 0
    E[c(1L, 5L, 30L), ] <- 0
    E <- Matrix::drop0(E)
    DE <- as.matrix(E)
    for (layout in c("csr", "csc")) {
        store <- file.path(tempfile(), paste0(layout, ".zarr"))
        write_sparse_zarr(store, E, layout)
        X <- ZarrSparseMatrix(store, "X")
        for (index in list(list(3:20, 1:8), list(c(30L, 2L, 5L), 6:20),
                           list(1:30, 2:2), list(NULL, 1:20),
                           list(10:12, NULL), list(c(1L, 5L), 1:20))) {
            i <- index[[1L]] %||% seq_len(nrow(DE))
            j <- index[[2L]] %||% seq_len(ncol(DE))
            a <- extract_array(DelayedArray::seed(X), index)
            expect_equal(a, DE[i, j, drop=FALSE], ignore_attr=TRUE)
            sa <- extract_sparse_array(DelayedArray::seed(X), index)
            expect_equal(as.matrix(sa), DE[i, j, drop=FALSE],
                         ignore_attr=TRUE)
        }
        expect_equal(as.matrix(t(X)[3:15, 2:25]), t(DE)[3:15, 2:25],
                     ignore_attr=TRUE)
    }
})

test_that("rows are filtered chunk by chunk", {
    by_chunk <- ZarrSparseMatrix:::.load_CSC_rows_by_chunk
    E <- M
    E[, c(1L, 2L, 7L, 20L)] <- 0
    E <- Matrix::drop0(E)
    DE <- as.matrix(E)
    store <- file.path(tempfile(), "csc.zarr")
    write_sparse_zarr(store, E, "csc")
    x <- DelayedArray::seed(ZarrSparseMatrix(store, "X"))
    for (index in list(list(c(30L, 2L, 5L), 1:20), list(3:20, c(1L, 4:9)),
                       list(1:30, 2L), list(4L, integer(0)))) {
        i <- index[[1L]]
        j <- index[[2L]]
        lookup <- integer(nrow(DE))
        lookup[i] <- seq_along(i)
        trip <- by_chunk(x, x@indptr_ranges[j, , drop=FALSE], lookup)
        expect_false(is.null(trip))
        A <- matrix(0, length(i), length(j))
        A[cbind(trip$row, trip$col)] <- trip$val
        expect_equal(A, DE[i, j, drop=FALSE], ignore_attr=TRUE)
        expect_identical(length(trip$row), sum(DE[i, j] != 0))
    }
})
