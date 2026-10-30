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
