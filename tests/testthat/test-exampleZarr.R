test_that("exampleZarr() returns a readable store", {
    store <- exampleZarr()
    expect_true(dir.exists(store))
    expect_true(file.exists(file.path(store, "X", ".zattrs")))
    ## A second call reuses the unzipped copy.
    expect_identical(exampleZarr(), store)

    X <- ZarrSparseMatrix(store, "X")
    expect_s4_class(DelayedArray::seed(X), "CSR_ZarrSparseMatrixSeed")
    expect_identical(dim(X), c(200L, 100L))
    expect_identical(DelayedArray::type(X), "integer")
})
