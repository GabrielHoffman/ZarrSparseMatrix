### =========================================================================
### ZarrSparseMatrix objects
### -------------------------------------------------------------------------
###
### A DelayedMatrix backend for a CSR or CSC sparse matrix stored in a Zarr
### group as 3 one-dimensional arrays: <group>/data, <group>/indices, and
### <group>/indptr (e.g. the X group of an AnnData .zarr store).
###
### Modeled after HDF5Array::H5SparseMatrix. Reading is done with Rarr.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers for reading Zarr metadata and 1D arrays
###

### Join path components with "/". All components are single strings except
### the last one, which can be a vector (the result is then vectorized).
.zarr_path <- function(filepath, ...)
{
    parts <- list(sub("/+$", "", filepath), ...)
    n <- length(parts)
    if (n == 1L)
        return(parts[[1L]])
    prefix <- paste(unlist(parts[-n]), collapse="/")
    paste(prefix, parts[[n]], sep="/")
}

### Read the JSON metadata of a Zarr node. Supports Zarr v2 (.zarray/.zgroup
### + .zattrs) and Zarr v3 (zarr.json).
.read_zarr_json <- function(path, file)
{
    f <- .zarr_path(path, file)
    if (!file.exists(f))
        return(NULL)
    read_json(f, simplifyVector=TRUE)
}

.read_zarr_node_attrs <- function(path)
{
    v3 <- .read_zarr_json(path, "zarr.json")
    if (!is.null(v3))
        return(v3$attributes)
    .read_zarr_json(path, ".zattrs")
}

.zarr_node_exists <- function(path)
{
    any(file.exists(.zarr_path(path, c(".zarray", ".zgroup", "zarr.json"))))
}

### Length of a 1D Zarr array, without reading its data.
.zarr_length <- function(path)
{
    meta <- .read_zarr_json(path, ".zarray")
    if (is.null(meta))
        meta <- .read_zarr_json(path, "zarr.json")
    if (is.null(meta))
        stop(wmsg("no Zarr array found at \"", path, "\""))
    shape <- meta$shape
    if (length(shape) != 1L)
        stop(wmsg("Zarr array at \"", path, "\" is not 1-dimensional"))
    as.numeric(shape)
}

### Read elements 'idx' (1-based) of a 1D Zarr array. Reads the whole array
### if 'idx' is NULL.
###
### For a local array, chunks are read one at a time (see
### .read_zarr_1d_by_chunk() below), with Rarr's internal chunk decoder or,
### if that fails (e.g. Rarr's internals changed), with one public
### Rarr::read_zarr_array() call per chunk. For a remote array, or if both
### fail, this falls back to a single Rarr::read_zarr_array() call.
.read_zarr_1d <- function(path, idx=NULL)
{
    if (!is.null(.zarr_int64_info(path)))
        return(.read_zarr_int64_1d(path, idx))
    if (!grepl("^[a-z]+://", path)) {
        ans <- tryCatch(.read_zarr_1d_by_chunk(path, idx),
                        error=function(e) NULL)
        if (is.null(ans))
            ans <- tryCatch(
                .read_zarr_1d_by_chunk(path, idx,
                                       reader=.public_chunk_reader(path)),
                error=function(e) NULL)
        if (!is.null(ans))
            return(ans)
    }
    if (is.null(idx))
        return(as.vector(read_zarr_array(path)))
    if (length(idx) == 0L) {
        ## Rarr can't read zero elements: read one to get the type, then drop.
        return(as.vector(read_zarr_array(path, index=list(1L)))[0L])
    }
    as.vector(read_zarr_array(path, index=list(idx)))
}

### Read the elements in ranges start[k], ..., start[k] + width[k] - 1 (1-based)
### of a 1D Zarr array, in order. Same as
### .read_zarr_1d(path, .ranges_to_positions(start, width)), but without
### building that vector of positions when the ranges are sorted and
### disjoint, as for sorted columns of a CSC matrix. Adjacent ranges (e.g.
### consecutive columns) are merged into one run, and each run is read as
### whole chunks, trimmed at its ends.
.read_zarr_1d_ranges <- function(path, start, width)
{
    nonempty <- width > 0
    start <- start[nonempty]
    end <- start + width[nonempty] - 1
    disjoint <- length(start) <= 1L ||
                all(start[-1L] > end[-length(end)])
    if (disjoint && is.null(.zarr_int64_info(path)) &&
        !grepl("^[a-z]+://", path)) {
        runs <- .merge_runs(start, end)
        start <- runs$start
        end <- runs$end
        ans <- tryCatch(.read_zarr_1d_runs_by_chunk(path, start, end),
                        error=function(e) NULL)
        if (is.null(ans))
            ans <- tryCatch(
                .read_zarr_1d_runs_by_chunk(path, start, end,
                                            reader=.public_chunk_reader(path)),
                error=function(e) NULL)
        if (!is.null(ans))
            return(ans)
    }
    .read_zarr_1d(path, .ranges_to_positions(start, end - start + 1))
}

### 64-bit integer arrays
###
### AnnData writes "indptr" as int64 when the matrix has more than 2^31 - 1
### nonzero values. Rarr converts 64-bit integers to R integers, which turns
### values above .Machine$integer.max into NA. To read them exactly (up to
### 2^53, as doubles), we let Rarr read the array as float64, which gives
### each value's 8 bytes unchanged, and then convert those bytes back.

### NULL if the Zarr array at 'path' isn't a local 64-bit integer array.
### Otherwise a list with the metadata file name and whether the type is
### signed.
.zarr_int64_info <- function(path)
{
    if (grepl("^[a-z]+://", path))
        return(NULL)
    meta <- .read_zarr_json(path, ".zarray")
    if (!is.null(meta)) {
        m <- regmatches(meta$dtype, regexec("^[<>|=]?([iu])8$", meta$dtype))
        if (length(m[[1L]]) == 0L)
            return(NULL)
        return(list(file=".zarray", signed=m[[1L]][[2L]] == "i"))
    }
    meta <- .read_zarr_json(path, "zarr.json")
    if (!is.null(meta) && identical(meta$node_type, "array") &&
        meta$data_type %in% c("int64", "uint64"))
        return(list(file="zarr.json", signed=meta$data_type == "int64"))
    NULL
}

.read_zarr_int64_1d <- function(path, idx=NULL)
{
    info <- .zarr_int64_info(path)
    path <- sub("/+$", "", normalizePath(path, mustWork=TRUE))

    ## A copy of the array that declares float64 instead of (u)int64: edited
    ## metadata, and symlinks to the original chunks (copies if symlinks
    ## aren't supported). The JSON is edited as text so that nothing else
    ## in it changes.
    tmp <- tempfile("zarr_int64_")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive=TRUE))
    json <- readLines(file.path(path, info$file), warn=FALSE)
    if (info$file == ".zarray") {
        json <- sub('("dtype"\\s*:\\s*"[<>|=]?)[iu]8"', '\\1f8"', json)
    } else {
        json <- sub('("data_type"\\s*:\\s*")u?int64"', '\\1float64"', json)
    }
    writeLines(json, file.path(tmp, info$file))
    entries <- setdiff(list.files(path, all.files=TRUE, no..=TRUE), info$file)
    from <- file.path(path, entries)
    to <- file.path(tmp, entries)
    linked <- suppressWarnings(file.symlink(from, to))
    if (!all(linked))
        file.copy(from[!linked], tmp, recursive=TRUE)

    x <- .read_zarr_1d(tmp, idx)

    ## Split each value's 8 bytes into its low and high 32-bit words, read
    ## as signed integers. Fix up the words that are negative or NA (R's NA
    ## integer is the bit pattern 0x80000000, i.e. -2^31).
    words <- readBin(writeBin(x, raw(), size=8L, endian="little"), "integer",
                     n=2L * length(x), size=4L, endian="little")
    unsigned <- function(w) {
        w <- as.numeric(w)
        w[is.na(w)] <- 2^31
        w[w < 0] <- w[w < 0] + 2^32
        w
    }
    lo <- unsigned(words[seq.int(1L, by=2L, length.out=length(x))])
    hi <- words[seq.int(2L, by=2L, length.out=length(x))]
    if (info$signed) {
        hi <- as.numeric(hi)
        hi[is.na(hi)] <- -2^31
    } else {
        hi <- unsigned(hi)
    }
    if (any(hi >= 2^21 | hi < -2^21))
        stop(wmsg("64-bit integer array \"", path, "\" has values that ",
                  "can't be represented exactly as doubles (beyond 2^53)"))
    hi * 2^32 + lo
}

### Reading 1D arrays chunk by chunk
###
### Rarr::read_zarr_array() is slow for a long 'index': for each chunk, it
### scans the whole index to find that chunk's elements, so its cost grows
### with (number of chunks) x length(idx). .read_zarr_1d_by_chunk() instead
### sorts 'idx' once, finds each chunk's elements with findInterval(), and
### decodes each needed chunk once. On a 1.2 billion element AnnData store,
### reading 43M elements went from 5 s to 0.7 s, and 240M from 86 s to 4 s.
###
### Chunks are decoded by a "chunk reader": a list with the array length 'n',
### the chunk length 'chunk', the R type of the values 'type', and a function
### 'read(cid)' that returns the values of chunk 'cid' (0-based). The fastest
### readers use Rarr's internal chunk decoder, whose interface differs
### between Rarr versions, so there is one per version. The public reader,
### which works with any version, is the fallback.

### Get an internal Rarr function. Using get() rather than ':::' keeps R CMD
### check from flagging the functions missing in the installed Rarr version.
.rarr_internal <- function(name) get(name, envir=asNamespace("Rarr"))

.r_type_of_zarr_base_type <- function(base_type)
{
    switch(base_type, bool="logical", int="integer", uint="integer",
           float="double", stop("unsupported data type"))
}

### The chunk reader for the installed version of Rarr, chosen from the
### arguments of its internal read_chunk(). 'path' ends with "/".
.rarr_chunk_reader <- function(path)
{
    read_chunk <- .rarr_internal("read_chunk")
    args <- names(formals(read_chunk))
    if ("zarr_array_path" %in% args) {
        .rarr1_chunk_reader(path, read_chunk)
    } else if ("chunk_path" %in% args) {
        .rarr2_0_chunk_reader(path, read_chunk)
    } else if ("chunk_bytes" %in% args) {
        .rarr2_1_chunk_reader(path, read_chunk)
    } else {
        stop("unknown Rarr version")
    }
}

### Rarr 1.x: Zarr v2 only. read_chunk() reads and decodes a chunk file.
.rarr1_chunk_reader <- function(path, read_chunk)
{
    meta <- .rarr_internal(".read_array_metadata")(path, ".zarray",
                                                   s3_client=NULL)
    if (length(meta$shape) != 1L)
        stop("not a 1D array")
    n <- as.numeric(meta$shape[[1L]])
    chunk <- as.numeric(meta$chunks[[1L]])
    ## For a C-order array, read_chunk() calls aperm() on each chunk, which
    ## copies it. C and Fortran order are the same in 1D, so skip that.
    meta$order <- "F"
    read <- function(cid) {
        ## zarr-python always writes full-size chunks, but Rarr's writer
        ## makes the last one shorter: 'alt_chunk_dim' accepts both.
        read_chunk(path, chunk_id=cid, metadata=meta,
                   alt_chunk_dim=min(chunk, n - cid * chunk))$chunk_data
    }
    list(n=n, chunk=chunk,
         type=.r_type_of_zarr_base_type(meta$datatype$base_type), read=read)
}

### Metadata of a Zarr v2 or v3 array in Rarr 2.x, with its decoders.
.rarr2_metadata <- function(path, metadata_file=NULL)
{
    read_meta <- .rarr_internal(".read_array_metadata")
    if ("metadata_file" %in% names(formals(read_meta))) {
        file <- if (file.exists(paste0(path, "zarr.json"))) "zarr.json"
                else ".zarray"
        meta <- read_meta(path, file, s3_client=NULL)
    } else {
        meta <- read_meta(path, s3_client=NULL)
    }
    if (length(meta$shape) != 1L)
        stop("not a 1D array")
    meta$configured_decoders <-
        .rarr_internal(".configure_codecs")(meta$codecs, operation="decode")
    meta
}

### File of chunk 'cid', e.g. "<path>0" (Zarr v2) or "<path>c/0" (Zarr v3).
.rarr2_chunk_file <- function(path, meta, cid)
{
    key <- matrix(as.integer(cid), ncol=1L)
    paste0(path, .rarr_internal(".create_chunk_names")(as.data.frame(key),
                                                        meta))
}

### Rarr 2.0: read_chunk() reads and decodes a chunk file.
.rarr2_0_chunk_reader <- function(path, read_chunk)
{
    meta <- .rarr2_metadata(path)
    n <- as.numeric(meta$shape[[1L]])
    chunk <- as.numeric(meta$chunk_grid$configuration$chunk_shape[[1L]])
    read <- function(cid) {
        len <- min(chunk, n - cid * chunk)
        f <- .rarr2_chunk_file(path, meta, cid)
        if (!file.exists(f))
            return(rep(meta$fill_value, len))
        read_chunk(f, metadata=meta, alt_chunk_dim=len)
    }
    list(n=n, chunk=chunk,
         type=.r_type_of_zarr_base_type(meta$datatype$base_type), read=read)
}

### Rarr 2.1: read_chunk() decodes the bytes of a chunk.
.rarr2_1_chunk_reader <- function(path, read_chunk)
{
    meta <- .rarr2_metadata(path)
    n <- as.numeric(meta$shape[[1L]])
    chunk <- as.numeric(meta$chunk_grid$configuration$chunk_shape[[1L]])
    decode <- function(bytes, dim)
        read_chunk(chunk_bytes=bytes, chunk_dim=dim,
                   decoders=meta$configured_decoders,
                   datatype=meta$datatype, fill_value=meta$fill_value)
    read <- function(cid) {
        len <- min(chunk, n - cid * chunk)
        f <- .rarr2_chunk_file(path, meta, cid)
        size <- file.size(f)  # NA if missing
        if (is.na(size))
            return(rep(meta$fill_value, len))
        bytes <- readBin(f, "raw", n=size)
        ## Full-size chunk as written by zarr-python, or a short last chunk
        ## as written by Rarr 1.x.
        tryCatch(decode(bytes, chunk), error=function(e) decode(bytes, len))
    }
    list(n=n, chunk=chunk,
         type=.r_type_of_zarr_base_type(meta$datatype$base_type), read=read)
}

### Any Rarr version: one Rarr::read_zarr_array() call per chunk, which is
### fast because each index is short. Zarr metadata is read with jsonlite.
.public_chunk_reader <- function(path)
{
    meta <- .read_zarr_json(path, ".zarray")
    if (!is.null(meta)) {
        chunk <- meta$chunks
        base <- sub("^[<>|=]?([biuf]).*$", "\\1", meta$dtype)
    } else {
        meta <- .read_zarr_json(path, "zarr.json")
        chunk <- meta$chunk_grid$configuration$chunk_shape
        base <- substr(meta$data_type, 1L, 1L)
    }
    if (length(meta$shape) != 1L)
        stop("not a 1D array")
    n <- as.numeric(meta$shape)
    chunk <- as.numeric(chunk)
    type <- switch(base, b="logical", i="integer", u="integer", f="double",
                   stop("unsupported data type"))
    read <- function(cid) {
        from <- cid * chunk + 1
        to <- min((cid + 1) * chunk, n)
        as.vector(read_zarr_array(path, index=list(from:to)))
    }
    list(n=n, chunk=chunk, type=type, read=read)
}

### Read elements 'idx' (1-based, any order, duplicates allowed) of a local
### 1D Zarr array with 'reader', which defaults to the fastest one for the
### installed Rarr.
.read_zarr_1d_by_chunk <- function(path, idx=NULL, reader=NULL)
{
    path <- paste0(sub("/+$", "", normalizePath(path, mustWork=TRUE)), "/")
    if (is.null(reader))
        reader <- .rarr_chunk_reader(path)
    n <- reader$n
    chunk <- reader$chunk

    ## Empty selection: no need to read anything, just get the type right.
    if (!is.null(idx) && length(idx) == 0L)
        return(vector(reader$type, 0L))

    perm <- NULL
    if (is.null(idx)) {
        idx <- seq_len(n)
        strict <- TRUE
    } else {
        idx <- as.numeric(idx)
        ## One pass in the common case of strictly increasing positions.
        ## is.unsorted() is NA if 'idx' has NAs. 'strict' (no duplicates)
        ## tells whether a chunk is fully selected when the number of
        ## positions in it equals its length.
        unsorted <- is.unsorted(idx, strictly=TRUE)
        if (is.na(unsorted))
            stop("subscript out of bounds")
        strict <- !unsorted
        if (unsorted && is.unsorted(idx)) {
            perm <- order(idx)
            idx <- idx[perm]
        }
        if (idx[[1L]] < 1 || idx[[length(idx)]] > n)
            stop("subscript out of bounds")
    }

    ## For each chunk overlapping [idx[1], idx[length(idx)]], the positions
    ## in 'idx' of its first and last selected elements: lo[k]..hi[k].
    ## 'idx' is known to be sorted and free of NAs, so findInterval() can
    ## skip checking it.
    first <- (idx[[1L]] - 1) %/% chunk
    last <- (idx[[length(idx)]] - 1) %/% chunk
    bounds <- findInterval(seq(first, last + 1) * chunk, idx,
                           checkSorted=FALSE, checkNA=FALSE) + 1
    lo <- bounds[-length(bounds)]
    hi <- bounds[-1L] - 1
    cids <- seq(first, last)
    needed <- which(hi >= lo)

    ## Decode each needed chunk once (Rarr returns 1D arrays; unlist() drops
    ## their "dim" attribute) and keep the selected elements.
    pieces <- lapply(needed, function(k) {
        cid <- cids[[k]]
        chunk_data <- reader$read(cid)
        if (strict && hi[[k]] - lo[[k]] + 1 == length(chunk_data))
            return(chunk_data)  # whole chunk selected
        chunk_data[idx[lo[[k]]:hi[[k]]] - cid * chunk]
    })
    ans <- unlist(pieces, use.names=FALSE)

    if (!is.null(perm))
        ans[perm] <- ans
    ans
}

### Merge ranges start[k]..end[k] (sorted and disjoint) that start right
### after the previous one ends, e.g. consecutive columns of a CSC matrix.
.merge_runs <- function(start, end)
{
    if (length(start) <= 1L)
        return(list(start=start, end=end))
    brk <- c(TRUE, start[-1L] != end[-length(end)] + 1)
    list(start=start[brk], end=end[c(brk[-1L], TRUE)])
}

### Call FUN(values, from, to, offset) on each chunk that overlaps the runs
### start[k]..end[k] (1-based, nonempty, sorted and disjoint) of local 1D
### Zarr arrays of the same length and chunking, read with 'readers' (see
### .read_zarr_1d_by_chunk()). 'values' are the decoded values of the chunk
### in each array, 'from'..'to' the selected positions in it (1-based, one
### pair per run), and 'offset' the position of the chunk in the arrays.
### Returns the results in order.
.map_runs_by_chunk <- function(readers, start, end, FUN)
{
    n <- readers[[1L]]$n
    chunk <- readers[[1L]]$chunk
    for (reader in readers)
        if (reader$n != n || reader$chunk != chunk)
            stop("arrays differ in length or chunking")
    if (length(start) == 0L)
        return(list())
    if (start[[1L]] < 1 || end[[length(end)]] > n)
        stop("subscript out of bounds")

    ## Split the runs at chunk boundaries into segments, each within one
    ## chunk: chunk 'cid' (0-based, nondecreasing) and positions 'from'..'to'
    ## in it (1-based).
    nseg <- (end - 1) %/% chunk - (start - 1) %/% chunk + 1
    run <- rep.int(seq_along(start), nseg)
    cid <- (start[run] - 1) %/% chunk + sequence(nseg) - 1
    offset <- cid * chunk
    from <- pmax(start[run] - offset, 1)
    to <- pmin(end[run] - offset, chunk)

    ## Decode each chunk once.
    groups <- split(seq_along(cid), cumsum(c(TRUE, diff(cid) != 0)))
    lapply(groups, function(s) {
        c1 <- cid[[s[[1L]]]]
        values <- lapply(readers, function(reader) reader$read(c1))
        FUN(values, from[s], to[s], c1 * chunk)
    })
}

### The selected positions 'from'..'to' (see .map_runs_by_chunk()) in a
### chunk of length 'len'. NULL if the whole chunk is selected.
.chunk_positions <- function(from, to, len)
{
    if (length(from) == 1L) {
        if (from == 1 && to >= len)
            return(NULL)
        return(from:to)
    }
    sequence(to - from + 1, from=from)
}

### Read the elements in runs start[k]..end[k] (1-based, nonempty, sorted
### and disjoint) of a local 1D Zarr array, in order, with 'reader' as for
### .read_zarr_1d_by_chunk().
.read_zarr_1d_runs_by_chunk <- function(path, start, end, reader=NULL)
{
    path <- paste0(sub("/+$", "", normalizePath(path, mustWork=TRUE)), "/")
    if (is.null(reader))
        reader <- .rarr_chunk_reader(path)
    pieces <- .map_runs_by_chunk(list(reader), start, end,
        function(values, from, to, offset) {
            v <- values[[1L]]
            p <- .chunk_positions(from, to, length(v))
            if (is.null(p)) v else v[p]
        })
    if (length(pieces) == 0L)
        return(vector(reader$type, 0L))
    unlist(pieces, use.names=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ZarrSparseMatrixSeed class
###

#' ZarrSparseMatrixSeed objects
#'
#' ZarrSparseMatrixSeed is a low-level helper class for representing a
#' pointer to a sparse matrix stored in a Zarr group in CSR or CSC layout,
#' as three 1D arrays \code{data}, \code{indices} and \code{indptr}
#' (e.g. the \code{X} group of an AnnData \code{.zarr} store).
#' It is the Zarr counterpart of \code{HDF5Array::H5SparseMatrixSeed}.
#'
#' Note that a ZarrSparseMatrixSeed object is not intended to be used
#' directly. Most end users will typically create and manipulate a
#' higher-level \link{ZarrSparseMatrix} object instead.
#'
#' @param filepath Path to the \code{.zarr} store.
#' @param group Name of the group in the store that contains the sparse
#'   matrix, e.g. \code{"X"}.
#' @param subdata Name of the array in \code{group} to use instead of
#'   \code{"data"} (e.g. a layer of alternative values that share the same
#'   \code{indices} and \code{indptr}). \code{NULL} means \code{"data"}.
#' @param dim Dimensions of the matrix (\code{c(nrow, ncol)}). By default
#'   they are read from the \code{shape} array or attribute of \code{group}.
#' @param sparse.layout \code{"CSR"} or \code{"CSC"}. By default it is read
#'   from the \code{encoding-type} (or \code{h5sparse_format}) attribute of
#'   \code{group}, falling back to \code{"CSR"}.
#' @param dimnames Optional list of 2 elements: row names and column names.
#'
#' @return A \code{CSC_ZarrSparseMatrixSeed} or \code{CSR_ZarrSparseMatrixSeed}
#'   object.
#'
#' @seealso \link{ZarrSparseMatrix}
#'
#' @examples
#' store <- exampleZarr()
#' seed <- ZarrSparseMatrixSeed(store, "X")
#' seed
#' dim(seed)
#' nzcount(seed)
#' extract_array(seed, list(1:3, 1:5))
#'
#' ## A CSR seed transposes to a CSC seed without reading any data
#' t(seed)
#'
#' @importClassesFrom S4Arrays Array
#' @importClassesFrom BiocGenerics OutOfMemoryObject
#' @importClassesFrom S4Vectors character_OR_NULL
#' @exportClass ZarrSparseMatrixSeed CSC_ZarrSparseMatrixSeed CSR_ZarrSparseMatrixSeed
#' @rdname ZarrSparseMatrixSeed
setClass("ZarrSparseMatrixSeed",
  contains=c("Array", "OutOfMemoryObject"),
  representation(
    "VIRTUAL",
    filepath="character",       # Path to the .zarr store
    group="character",          # Name of the group in the store
    subdata="character_OR_NULL",
    dim="integer",
    indptr_ranges="data.frame", # 1-based start + width of each
          # column (CSC) or row (CSR)
    dimnames="list"
  ),
  prototype(
    subdata=NULL,
    dimnames=list(NULL, NULL)
  )
)

#' @rdname ZarrSparseMatrixSeed
setClass("CSC_ZarrSparseMatrixSeed", contains="ZarrSparseMatrixSeed")

#' @rdname ZarrSparseMatrixSeed
setClass("CSR_ZarrSparseMatrixSeed", contains="ZarrSparseMatrixSeed")

.get_data_name <- function(subdata)
{
  if (is.null(subdata)) "data" else subdata
}

.data_path <- function(x)
  .zarr_path(x@filepath, x@group, .get_data_name(x@subdata))

.indices_path <- function(x) .zarr_path(x@filepath, x@group, "indices")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and basic methods
###

#' @importFrom BiocGenerics path
#' @exportMethod path
#' @rdname ZarrSparseMatrixSeed
#' @param object,x A ZarrSparseMatrixSeed object.
setMethod("path", "ZarrSparseMatrixSeed", function(object) object@filepath)

#' @exportMethod dim
#' @rdname ZarrSparseMatrixSeed
setMethod("dim", "ZarrSparseMatrixSeed", function(x) x@dim)

#' @exportMethod dimnames
#' @rdname ZarrSparseMatrixSeed
setMethod("dimnames", "ZarrSparseMatrixSeed",
  function(x)
  {
    if (all(vapply(x@dimnames, is.null, logical(1L))))
      return(NULL)
    x@dimnames
  }
)

#' @importFrom S4Arrays is_sparse
#' @exportMethod is_sparse
#' @rdname ZarrSparseMatrixSeed
setMethod("is_sparse", "ZarrSparseMatrixSeed", function(x) TRUE)

#' @importFrom SparseArray nzcount
#' @exportMethod nzcount
#' @rdname ZarrSparseMatrixSeed
setMethod("nzcount", "ZarrSparseMatrixSeed",
    function(x) .zarr_length(.data_path(x))
)

### Transposing a CSR seed gives a CSC seed (and vice versa) at no cost.

#' @importFrom BiocGenerics t
#' @exportMethod t
#' @rdname ZarrSparseMatrixSeed
setMethod("t", "CSR_ZarrSparseMatrixSeed",
  function(x)
  {
    x@dim <- rev(x@dim)
    x@dimnames <- rev(x@dimnames)
    class(x) <- class(new("CSC_ZarrSparseMatrixSeed"))
    x
  }
)

#' @rdname ZarrSparseMatrixSeed
setMethod("t", "CSC_ZarrSparseMatrixSeed",
  function(x)
  {
    x@dim <- rev(x@dim)
    x@dimnames <- rev(x@dimnames)
    class(x) <- class(new("CSR_ZarrSparseMatrixSeed"))
    x
  }
)

#' @importFrom methods show
#' @importFrom BiocGenerics type
#' @exportMethod show
#' @rdname ZarrSparseMatrixSeed
setMethod("show", "ZarrSparseMatrixSeed",
  function(object)
  {
    cat(sprintf("<%s> sparse %s object of type \"%s\":\n",
                paste(dim(object), collapse=" x "),
                class(object), type(object)))
    cat("# dirname: ", dirname(object@filepath), "\n", sep="")
    cat("# basename: ", basename(object@filepath), "\n", sep="")
    cat("# group: ", object@group, "\n", sep="")
  }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

.read_zarr_sparse_dim <- function(group_path)
{
  ## A 'shape' array inside the group (10x style) takes precedence.
  shape_path <- .zarr_path(group_path, "shape")
  if (.zarr_node_exists(shape_path))
    return(as.integer(.read_zarr_1d(shape_path)))
  attrs <- .read_zarr_node_attrs(group_path)
  shape <- attrs$shape
  if (is.null(shape))
    shape <- attrs$h5sparse_shape
  if (is.null(shape))
    stop(wmsg("Zarr group \"", group_path, "\" contains no 'shape' ",
      "array and has no 'shape' or 'h5sparse_shape' ",
      "attribute. As a consequence, the dimensions of the ",
      "sparse matrix can't be determined. Please supply ",
      "them via the 'dim' argument."))
  ## Unlike HDF5 attributes read by rhdf5, JSON attributes are already
  ## in row-major (Python) order, i.e. c(nrow, ncol).
  as.integer(shape)
}

.read_zarr_sparse_layout <- function(group_path)
{
  if (.zarr_node_exists(.zarr_path(group_path, "shape")))
    return("CSR")
  attrs <- .read_zarr_node_attrs(group_path)
  layout <- attrs[["encoding-type"]]
  if (is.null(layout))
    layout <- attrs[["h5sparse_format"]]
  if (is.null(layout))
    return("CSR")
  ans <- toupper(substr(layout, 1L, 3L))
  if (!(ans %in% c("CSR", "CSC")))
    stop(wmsg("sparse matrix in Zarr group \"", group_path, "\" is ",
      "stored in unsupported layout \"", layout, "\""))
  ans
}

#' @importFrom methods new
#' @importFrom S4Vectors wmsg
#' @importFrom jsonlite read_json
#' @importFrom Rarr read_zarr_array
#' @export
#' @rdname ZarrSparseMatrixSeed
ZarrSparseMatrixSeed <- function(
  filepath, group, subdata=NULL,
  dim=NULL, sparse.layout=NULL,
  dimnames=NULL)
{
if (!(is.character(filepath) && length(filepath) == 1L && !is.na(filepath)))
    stop(wmsg("'filepath' must be a single string"))
  if (!grepl("^[a-z]+://", filepath))
    filepath <- normalizePath(filepath, mustWork=TRUE)
  if (!(is.character(group) && length(group) == 1L && ! is.na(group)))
    stop(wmsg("'group' must be a single string"))
  group <- gsub("^/+|/+$", "", group)
  group_path <- .zarr_path(filepath, group)

  for (comp in c(.get_data_name(subdata), "indices", "indptr"))
    if (!.zarr_node_exists(.zarr_path(group_path, comp)))
      stop(wmsg("Zarr group \"", group_path, "\" has no \"", comp, "\" array"))

  if (is.null(dim)) {
    dim <- .read_zarr_sparse_dim(group_path)
  } else {
    dim <- as.integer(dim)
    stopifnot(length(dim) == 2L, !anyNA(dim), all(dim >= 0L))
  }
  if (is.null(sparse.layout)) {
    layout <- .read_zarr_sparse_layout(group_path)
  } else {
    layout <- match.arg(toupper(sparse.layout), c("CSC", "CSR"))
  }
  if (layout == "CSC") {
    expected_indptr_len <- dim[[2L]] + 1L
    ans_class <- "CSC_ZarrSparseMatrixSeed"
  } else {
    expected_indptr_len <- dim[[1L]] + 1L
    ans_class <- "CSR_ZarrSparseMatrixSeed"
  }

  data_len <- .zarr_length(.zarr_path(group_path, .get_data_name(subdata)))
  indices_len <- .zarr_length(.zarr_path(group_path, "indices"))
  if (indices_len != data_len)
    stop(wmsg("\"data\" and \"indices\" arrays have different lengths"))

  indptr <- as.numeric(.read_zarr_1d(.zarr_path(group_path, "indptr")))
  if (length(indptr) != expected_indptr_len ||
    indptr[[1L]] != 0 || indptr[[length(indptr)]] != data_len)
    stop(wmsg("\"indptr\" array is inconsistent with the ",
                "dimensions or layout of the sparse matrix"))
  indptr_ranges <- data.frame(start=indptr[-length(indptr)] + 1,
                              width=as.integer(diff(indptr)))

  if (is.null(dimnames)) {
    dimnames <- list(NULL, NULL)
  } else {
    stopifnot(is.list(dimnames), length(dimnames) == 2L)
  }

  new(ans_class, filepath=filepath, group=group, subdata=subdata, dim=dim, indptr_ranges=indptr_ranges, dimnames=dimnames)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Data extraction
###

### Load the columns 'j' (sorted, unique, or NULL for all) and the rows 'i'
### (unique, any order, or NULL for all) of a CSC seed as triplets: a list
### with the dimensions length(i) x length(j), and the row, column and value
### of each nonzero.
### The positions start[k], ..., start[k] + width[k] - 1 for each k, as
### doubles. sequence(width, from=start) would do the same, but it coerces
### 'start' to integer, which gives NA for positions beyond 2^31 - 1 (in
### matrices with more nonzero values than that).
.ranges_to_positions <- function(start, width)
{
  rep.int(start - 1, width) + sequence(width)
}

.load_CSC_triplets <- function(x, i=NULL, j=NULL)
{
  ranges <- x@indptr_ranges
  if (!is.null(j))
    ranges <- ranges[j, , drop=FALSE]
  widths <- ranges$width
  ncol <- nrow(ranges)

  nrow <- x@dim[[1L]]
  if (!is.null(i)) {
    ## Keep only the selected rows, before building the COO_SparseArray.
    ## 'lookup' maps each row to its position in 'i' (0 if not selected),
    ## which is much faster than match() on long 'row_idx'.
    lookup <- integer(nrow)
    lookup[i] <- seq_along(i)
    dim <- c(length(i), ncol)
    trip <- .load_CSC_rows_by_chunk(x, ranges, lookup)
    if (!is.null(trip))
      return(c(list(dim=dim), trip))
  }

  if (is.null(j)) {
    ## All columns: read the whole data and indices arrays.
    read <- function(path) .read_zarr_1d(path)
  } else {
    ## The ranges in data/indices of the selected columns.
    read <- function(path) .read_zarr_1d_ranges(path, ranges$start, widths)
  }
  nzdata <- read(.data_path(x))
  row_idx <- as.integer(read(.indices_path(x))) + 1L

  if (is.null(i)) {
    col_idx <- rep.int(seq_len(ncol), widths)
  } else {
    row_idx <- lookup[row_idx]
    keep <- which(row_idx != 0L)
    row_idx <- row_idx[keep]
    nzdata <- nzdata[keep]
    ## The column of each kept value, from its position: column k holds
    ## positions cumsum(widths)[k-1] + 1, ..., cumsum(widths)[k]. Cheaper
    ## than a column for every loaded value when most rows are dropped.
    col_idx <- findInterval(keep, cumsum(as.numeric(widths)) + 1) + 1L
    nrow <- length(i)
  }
  list(dim=c(nrow, ncol), row=row_idx, col=col_idx, val=nzdata)
}

### The row, column and value of the nonzeros in the rows that 'lookup'
### selects (see .load_CSC_triplets()) and the columns 'ranges' (sorted) of
### a CSC seed. 'data' and 'indices' are read together chunk by chunk, and
### the other rows are dropped from each chunk as it is decoded, so the
### values of the dropped rows are never all held in memory. NULL if the
### arrays can't be read this way (remote, 64-bit or chunked differently).
.load_CSC_rows_by_chunk <- function(x, ranges, lookup)
{
  paths <- c(.data_path(x), .indices_path(x))
  if (any(grepl("^[a-z]+://", paths)) ||
      !all(vapply(paths, function(p) is.null(.zarr_int64_info(p)), TRUE)))
    return(NULL)
  readers <- tryCatch(lapply(paths, function(p) {
      .rarr_chunk_reader(paste0(normalizePath(p, mustWork=TRUE), "/"))
    }), error=function(e) NULL)
  if (is.null(readers))
    return(NULL)

  cols <- which(ranges$width > 0)
  start <- ranges$start[cols]
  runs <- .merge_runs(start, start + ranges$width[cols] - 1)
  select <- function(v, from, to) {
    p <- .chunk_positions(from, to, length(v))
    if (is.null(p)) v else v[p]
  }
  pieces <- tryCatch(.map_runs_by_chunk(readers, runs$start, runs$end,
      function(values, from, to, offset) {
        row <- lookup[select(values[[2L]], from, to) + 1L]
        keep <- which(row != 0L)
        p <- .chunk_positions(from, to, length(values[[2L]]))
        list(row=row[keep], val=select(values[[1L]], from, to)[keep],
             pos=offset + (if (is.null(p)) keep else p[keep]))
      }), error=function(e) NULL)
  if (is.null(pieces))
    return(NULL)
  if (length(pieces) == 0L)
    return(list(row=integer(0), col=integer(0),
                val=vector(readers[[1L]]$type, 0L)))
  get <- function(name) unlist(lapply(pieces, `[[`, name), use.names=FALSE)
  ## The column of each kept value, from its position in data/indices:
  ## the last nonempty selected column that starts at or before it.
  pos <- get("pos")
  list(row=get("row"), col=cols[findInterval(pos, start)], val=get("val"))
}

### The triplets of a CSC seed as a COO_SparseArray.
.load_CSC_ZarrSparseMatrixSeed <- function(x, i=NULL, j=NULL)
{
  trip <- .load_CSC_triplets(x, i=i, j=j)
  nzcoo <- cbind(trip$row, trip$col, deparse.level=0L)
  COO_SparseArray(trip$dim, nzcoo, trip$val, check=FALSE)
}

### The triplets of a CSC seed as an ordinary matrix, filled directly
### instead of through a sparse object. With 'transpose=TRUE', the
### transposed matrix is filled instead, which is how a CSR seed is read.
.load_CSC_as_dense <- function(x, i=NULL, j=NULL, transpose=FALSE)
{
  trip <- .load_CSC_triplets(x, i=i, j=j)
  zero <- vector(typeof(trip$val), 1L)
  if (transpose) {
    ans <- matrix(zero, trip$dim[[2L]], trip$dim[[1L]])
    ans[cbind(trip$col, trip$row, deparse.level=0L)] <- trip$val
  } else {
    ans <- matrix(zero, trip$dim[[1L]], trip$dim[[2L]])
    ans[cbind(trip$row, trip$col, deparse.level=0L)] <- trip$val
  }
  ans
}

### Normalize 'index' of a CSC seed for loading: the columns are sorted and
### deduplicated and the rows deduplicated (rows are filtered while loading,
### in any order). 'i0'/'j0' are the indices to reorder or expand the loaded
### result with afterwards, or NULL if not needed.
.normalize_CSC_index <- function(index)
{
  j <- index[[2L]]
  j0 <- NULL
  if (!is.null(j)) {
    if (!is.integer(j))
      j <- as.integer(j)
    if (!isStrictlySorted(j)) {
      ju <- unique(sort(j))
      j0 <- match(j, ju)
      j <- ju
    }
  }
  i <- index[[1L]]
  i0 <- NULL
  if (!is.null(i)) {
    if (!is.integer(i))
      i <- as.integer(i)
    if (anyDuplicated(i) != 0L) {
      iu <- unique(i)
      i0 <- match(i, iu)
      i <- iu
    }
  }
  list(i=i, j=j, i0=i0, j0=j0)
}

### Dense extraction from a CSC seed, or from a CSR seed via its CSC
### transpose 't(x)' and reversed index.
.extract_array_from_CSC_ZarrSparseMatrixSeed <- function(
    x, index, transpose=FALSE)
{
  idx <- .normalize_CSC_index(index)
  ans <- .load_CSC_as_dense(x, i=idx$i, j=idx$j, transpose=transpose)
  if (is.null(idx$i0) && is.null(idx$j0))
    return(ans)
  ## Dimensions of the loaded result in CSC orientation.
  d <- if (transpose) rev(dim(ans)) else dim(ans)
  i0 <- idx$i0 %||% seq_len(d[[1L]])
  j0 <- idx$j0 %||% seq_len(d[[2L]])
  if (transpose) ans[j0, i0, drop=FALSE] else ans[i0, j0, drop=FALSE]
}

#' @importFrom S4Vectors isStrictlySorted
#' @importFrom methods as
#' @importClassesFrom SparseArray SVT_SparseArray
.extract_sparse_array_from_CSC_ZarrSparseMatrixSeed <- function(x, index)
{
  ## Duplicated rows and columns are loaded once and expanded afterwards.
  idx <- .normalize_CSC_index(index)
  coo <- .load_CSC_ZarrSparseMatrixSeed(x, i=idx$i, j=idx$j)
  index2 <- list(idx$i0, idx$j0)
  if (is.null(index2[[1L]]) && is.null(index2[[2L]]))
    return(coo)
  ## Subset an SVT_SparseArray, not the COO_SparseArray: in SparseArray
  ## 1.10.10, extract_sparse_array() on a COO_SparseArray returns zeros for
  ## repeated columns, e.g. index2 = list(NULL, c(1L, 1L)).
  extract_sparse_array(as(coo, "SVT_SparseArray"), index2)
}

#' @importFrom SparseArray extract_sparse_array COO_SparseArray
#' @exportMethod extract_sparse_array
#' @rdname ZarrSparseMatrixSeed
#' @param index A list of 2 elements (\code{NULL} or integer vector),
#'   as for \code{S4Arrays::extract_array()}.
setMethod("extract_sparse_array", "CSC_ZarrSparseMatrixSeed",
  function(x, index)
    .extract_sparse_array_from_CSC_ZarrSparseMatrixSeed(x, index)
)

#' @rdname ZarrSparseMatrixSeed
setMethod("extract_sparse_array", "CSR_ZarrSparseMatrixSeed",
    function(x, index)
        t(.extract_sparse_array_from_CSC_ZarrSparseMatrixSeed(t(x),rev(index)))
)

#' @importFrom S4Arrays extract_array
#' @exportMethod extract_array
#' @rdname ZarrSparseMatrixSeed
setMethod("extract_array", "CSC_ZarrSparseMatrixSeed",
  function(x, index) .extract_array_from_CSC_ZarrSparseMatrixSeed(x, index)
)

#' @rdname ZarrSparseMatrixSeed
setMethod("extract_array", "CSR_ZarrSparseMatrixSeed",
  function(x, index)
    .extract_array_from_CSC_ZarrSparseMatrixSeed(t(x), rev(index),
                                                 transpose=TRUE)
)

### Chunk geometry: a CSC matrix is best read by column, a CSR by row.

#' @importFrom DelayedArray chunkdim
#' @exportMethod chunkdim
#' @rdname ZarrSparseMatrixSeed
setMethod("chunkdim", "CSC_ZarrSparseMatrixSeed",
    function(x) c(x@dim[[1L]], 1L)
)

#' @rdname ZarrSparseMatrixSeed
setMethod("chunkdim", "CSR_ZarrSparseMatrixSeed",
    function(x) c(1L, x@dim[[2L]])
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ZarrSparseMatrix class
###

#' Sparse matrices stored in a Zarr store
#'
#' The ZarrSparseMatrix class is a \link[DelayedArray]{DelayedMatrix}
#' subclass for representing and operating on a sparse matrix stored in
#' CSR or CSC layout in a Zarr group, as three 1D arrays \code{data},
#' \code{indices} and \code{indptr} (e.g. the \code{X} group of an AnnData
#' \code{.zarr} store). It is the Zarr counterpart of
#' \code{HDF5Array::H5SparseMatrix}.
#'
#' All the operations available for \link[DelayedArray]{DelayedMatrix}
#' objects work on ZarrSparseMatrix objects. Data is only read from disk
#' when needed, and subsetting by column (CSC) or by row (CSR) only reads
#' the corresponding elements of \code{data} and \code{indices}.
#'
#' \code{t()} on a ZarrSparseMatrix transposes its seed, so the result is
#' another ZarrSparseMatrix: a CSR seed becomes a CSC seed and vice versa,
#' and no data is read. Reading from the result is as fast as reading from
#' the original, whereas a delayed transpose (what \code{t()} does on other
#' DelayedMatrix objects) transposes every block after reading it. Once a
#' ZarrSparseMatrix is subset or modified, it becomes a DelayedMatrix and
#' \code{t()} is delayed as usual.
#'
#' @param filepath Path to the \code{.zarr} store, or a
#'   \link{ZarrSparseMatrixSeed} object.
#' @param group Name of the group in the store that contains the sparse
#'   matrix, e.g. \code{"X"}.
#' @param ... Further arguments passed to \code{\link{ZarrSparseMatrixSeed}()}
#'   (\code{subdata}, \code{dim}, \code{sparse.layout}, \code{dimnames}).
#' @param seed A ZarrSparseMatrixSeed object.
#' @param x A ZarrSparseMatrix object.
#'
#' @return A ZarrSparseMatrix object.
#'
#' @seealso \link{ZarrSparseMatrixSeed}, \code{HDF5Array::H5SparseMatrix}
#'
#' @examples
#' ## Simulated counts for 200 cells x 100 genes, stored like the X group
#' ## of an AnnData .zarr store (CSR layout)
#' store <- exampleZarr()
#' X <- ZarrSparseMatrix(store, "X")
#' X
#'
#' X[1:5, 1:10]
#' colSums(X)
#'
#' ## Genes as rows, cells as columns: t() gives a ZarrSparseMatrix
#' ## backed by a CSC seed, without reading any data
#' tX <- t(X)
#' tX
#' seed(tX)
#'
#' ## Load into memory
#' as(X, "dgCMatrix")
#'
#' @importClassesFrom DelayedArray DelayedMatrix
#' @exportClass ZarrSparseMatrix
#' @rdname ZarrSparseMatrix
setClass("ZarrSparseMatrix",
    contains="DelayedMatrix",
    representation(seed="ZarrSparseMatrixSeed")
)

#' @importFrom DelayedArray DelayedArray new_DelayedArray
#' @exportMethod DelayedArray
#' @rdname ZarrSparseMatrix
setMethod("DelayedArray", "ZarrSparseMatrixSeed",
    function(seed) new_DelayedArray(seed, Class="ZarrSparseMatrix")
)

#' @importFrom DelayedArray seed
#' @exportMethod seed
#' @rdname ZarrSparseMatrix
setMethod("seed", "ZarrSparseMatrix",
  function(x) x@seed
)

### Transpose the seed rather than delaying the transpose. This reads no
### data (a CSC seed becomes a CSR seed and vice versa), and the result is
### read directly in the transposed layout. A delayed transpose would make
### DelayedArray call aperm() on every block it reads, which copies it.
#' @exportMethod t
#' @rdname ZarrSparseMatrix
setMethod("t", "ZarrSparseMatrix",
  function(x) DelayedArray(t(seed(x)))
)

#' @importFrom methods is
#' @export
#' @rdname ZarrSparseMatrix
ZarrSparseMatrix <- function(filepath, group, ...)
{
  if (is(filepath, "ZarrSparseMatrixSeed")) {
    if (!missing(group) || length(list(...)) != 0L)
      stop(wmsg("ZarrSparseMatrix() must be called with a single ", "argument when passed a ZarrSparseMatrixSeed object"))
    seed <- filepath
  } else {
    seed <- ZarrSparseMatrixSeed(filepath, group, ...)
  }
  DelayedArray(seed)
}
