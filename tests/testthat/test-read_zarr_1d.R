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
    rarr_write(array(int), short_path, chunk_dim=10L)
    expect_identical(read_1d(short_path, c(2L, 101L, 103L)),
                     int[c(2L, 101L, 103L)])
    expect_identical(read_1d(short_path), int)
})

test_that("every chunk reader gives the same values", {
    i <- c(50L, 3L, 103L, 3L, 12L, 50L, 11L)
    expect_identical(read_by_chunk(dbl_path, i,
        reader=ZarrSparseMatrix:::.public_chunk_reader(dbl_path)), dbl[i])
    expect_identical(read_by_chunk(int_path, i,
        reader=ZarrSparseMatrix:::.public_chunk_reader(int_path)), int[i])
    expect_identical(read_by_chunk(int_path, integer(0),
        reader=ZarrSparseMatrix:::.public_chunk_reader(int_path)), integer(0))
})

test_that("Zarr v3 arrays are read (Rarr 2.x only)", {
    skip_if_not("zarr_version" %in% names(formals(Rarr::write_zarr_array)))
    v3_path <- file.path(dir, "v3")
    rarr_write(array(dbl), v3_path, chunk_dim=10L, zarr_version=3L)
    i <- c(103L, 1L, 55L, 55L)
    expect_identical(read_by_chunk(v3_path, i), dbl[i])
    expect_identical(read_1d(v3_path), dbl)
})

test_that("out-of-range indices are rejected", {
    expect_error(read_by_chunk(dbl_path, 104L))
    expect_error(read_by_chunk(dbl_path, 0L))
})

test_that("reading ranges matches reading their positions", {
    read_ranges <- ZarrSparseMatrix:::.read_zarr_1d_ranges
    positions <- ZarrSparseMatrix:::.ranges_to_positions
    cases <- list(
        list(start=1, width=103),                  # whole array
        list(start=c(3, 8, 15), width=c(5, 7, 2)), # adjacent: one run
        list(start=c(2, 25, 61), width=c(3, 30, 1)),  # across chunks
        list(start=c(11, 31), width=c(10, 10)),    # whole chunks
        list(start=c(5, 9, 9, 40), width=c(4, 0, 3, 64)),  # empty range
        list(start=c(1, 50), width=c(0, 0)),       # nothing selected
        list(start=c(20, 5), width=c(3, 3)),       # unsorted
        list(start=c(5, 6), width=c(4, 4))         # overlapping
    )
    for (r in cases) {
        idx <- positions(r$start, r$width)
        expect_identical(read_ranges(dbl_path, r$start, r$width), dbl[idx])
        expect_identical(read_ranges(int_path, r$start, r$width), int[idx])
    }
    short_path <- file.path(dir, "short_ranges")
    rarr_write(array(int), short_path, chunk_dim=10L)
    expect_identical(read_ranges(short_path, c(95, 101), c(4, 3)),
                     int[c(95:98, 101:103)])
    runs <- ZarrSparseMatrix:::.read_zarr_1d_runs_by_chunk
    expect_identical(runs(int_path, c(3, 45), c(12, 101),
        reader=ZarrSparseMatrix:::.public_chunk_reader(int_path)),
        int[c(3:12, 45:101)])
    expect_error(runs(int_path, 100, 104))
})
