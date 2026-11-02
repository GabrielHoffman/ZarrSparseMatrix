read_1d <- ZarrSparseMatrix:::.read_zarr_1d
read_by_chunk <- ZarrSparseMatrix:::.read_zarr_1d_by_chunk

set.seed(2)
dbl <- round(rnorm(103), 3)
int <- sample.int(1000L, 103L)

dir <- tempfile()
dir.create(dir)
dbl_path <- file.path(dir, "dbl")
int_path <- file.path(dir, "int")
write_1d_zarr(dbl, dbl_path, 10L)
write_1d_zarr(int, int_path, 10L)

test_that("reading by chunk matches the stored values", {
    for (i in list(1:103, 5L, 103L, c(1L, 10L, 11L, 20L, 21L, 103L),
                   95:103, 37:64)) {
        expect_identical(read_by_chunk(dbl_path, i), dbl[i])
        expect_identical(read_by_chunk(int_path, i), int[i])
    }
    expect_identical(read_by_chunk(dbl_path), dbl)
    expect_identical(read_by_chunk(int_path), int)
})

test_that("unsorted and duplicated indices are returned in order", {
    i <- c(50L, 3L, 103L, 3L, 12L, 50L)
    expect_identical(read_by_chunk(dbl_path, i), dbl[i])
    expect_identical(read_by_chunk(int_path, i), int[i])
    expect_identical(read_1d(dbl_path, i), dbl[i])
    ## As many elements as a chunk has, but with a duplicate.
    i <- c(1:9, 9L)
    expect_identical(read_by_chunk(int_path, i), int[i])
})

test_that("empty selections have the right type without reading data", {
    expect_identical(read_by_chunk(dbl_path, integer(0)), double(0))
    expect_identical(read_by_chunk(int_path, integer(0)), integer(0))
    expect_identical(read_1d(int_path, integer(0)), integer(0))
})

test_that("Rarr's short last chunk is handled", {
    ## write_zarr_array() writes the last chunk shorter than the others.
    short_path <- file.path(dir, "short")
    Rarr::write_zarr_array(array(int), short_path, chunk_dim=10L)
    expect_identical(read_1d(short_path, c(2L, 101L, 103L)),
                     int[c(2L, 101L, 103L)])
    expect_identical(read_1d(short_path), int)
})

test_that("out-of-range indices are rejected", {
    expect_error(read_by_chunk(dbl_path, 104L))
    expect_error(read_by_chunk(dbl_path, 0L))
})
