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
### For a local Zarr v2 array, chunks are read one at a time with Rarr's
### internal read_chunk() (see .read_zarr_1d_by_chunk() below). For anything
### else, or if that fails (e.g. Rarr's internals changed), this falls back
### to Rarr::read_zarr_array().
.read_zarr_1d <- function(path, idx=NULL)
{
    if (!grepl("^[a-z]+://", path)) {
        ans <- tryCatch(.read_zarr_1d_by_chunk(path, idx),
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

### Rarr::read_zarr_array() is slow for a long 'index': for each chunk, it
### scans the whole index to find that chunk's elements, so its cost grows
### with (number of chunks) x length(idx). Here we sort 'idx' once, find each
### chunk's elements with findInterval(), and decompress each needed chunk
### once with Rarr's read_chunk(). On a 1.2 billion element AnnData store,
### reading 43M elements went from 5 s to 0.7 s, and 240M from 86 s to 4 s.
.read_zarr_1d_by_chunk <- function(path, idx=NULL)
{
    path <- paste0(sub("/+$", "", normalizePath(path, mustWork=TRUE)), "/")
    meta <- Rarr:::.read_array_metadata(path, ".zarray", s3_client=NULL)
    if (length(meta$shape) != 1L)
        stop("not a 1D array")
    n <- as.numeric(meta$shape[[1L]])
    chunk <- as.numeric(meta$chunks[[1L]])
    ## For a C-order array, read_chunk() calls aperm() on each chunk, which
    ## copies it. C and Fortran order are the same in 1D, so skip that.
    meta$order <- "F"

    ## Empty selection: no need to read anything, just get the type right.
    if (!is.null(idx) && length(idx) == 0L) {
        return(switch(meta$datatype$base_type,
                      bool=logical(0), int=integer(0), uint=integer(0),
                      float=double(0),
                      stop("unsupported data type")))
    }

    if (is.null(idx)) {
        idx <- seq_len(n)
        perm <- NULL
    } else {
        idx <- as.numeric(idx)
        if (anyNA(idx) || min(idx) < 1 || max(idx) > n)
            stop("subscript out of bounds")
        perm <- NULL
        if (is.unsorted(idx)) {
            perm <- order(idx)
            idx <- idx[perm]
        }
    }

    ## With no duplicates, a chunk is fully selected when the number of
    ## selected elements equals its length.
    strict <- isStrictlySorted(idx)
    first <- (idx[[1L]] - 1) %/% chunk
    last <- (idx[[length(idx)]] - 1) %/% chunk
    ## bounds[k] is the position in 'idx' of the first element in chunk
    ## (first + k - 1); bounds[k + 1] - 1 is the position of its last.
    bounds <- findInterval(seq(first, last + 1) * chunk, idx) + 1L
    pieces <- vector("list", last - first + 1)
    for (k in seq_along(pieces)) {
        lo <- bounds[[k]]
        hi <- bounds[[k + 1L]] - 1L
        if (hi < lo)
            next
        cid <- first + k - 1
        ## zarr-python always writes full-size chunks, but Rarr's writer
        ## makes the last one shorter: 'alt_chunk_dim' accepts both.
        chunk_data <- Rarr:::read_chunk(path, chunk_id=cid, metadata=meta,
                                        alt_chunk_dim=min(chunk, n - cid * chunk)
                      )$chunk_data
        chunk_data <- as.vector(chunk_data)
        ## Whole chunk selected: no need to subset it.
        if (strict && hi - lo + 1 == length(chunk_data)) {
            pieces[[k]] <- chunk_data
        } else {
            pieces[[k]] <- chunk_data[idx[lo:hi] - cid * chunk]
        }
    }
    ans <- unlist(pieces, use.names=FALSE)

    if (!is.null(perm))
        ans[perm] <- ans
    ans
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
.load_CSC_triplets <- function(x, i=NULL, j=NULL)
{
  ranges <- x@indptr_ranges
  if (!is.null(j))
    ranges <- ranges[j, , drop=FALSE]
  widths <- ranges$width
  ncol <- nrow(ranges)

  if (is.null(j)) {
    ## All columns: read the whole data and indices arrays.
    idx <- NULL
  } else {
    ## Positions (1-based) in data/indices of the selected columns.
    idx <- sequence(widths, from=ranges$start)
  }
  nzdata <- .read_zarr_1d(.data_path(x), idx)
  row_idx <- as.integer(.read_zarr_1d(.indices_path(x), idx)) + 1L
  col_idx <- rep.int(seq_len(ncol), widths)

  nrow <- x@dim[[1L]]
  if (!is.null(i)) {
    ## Keep only the selected rows, before building the COO_SparseArray.
    ## 'lookup' maps each row to its position in 'i' (0 if not selected),
    ## which is much faster than match() on long 'row_idx'.
    lookup <- integer(nrow)
    lookup[i] <- seq_along(i)
    row_idx <- lookup[row_idx]
    keep <- row_idx != 0L
    row_idx <- row_idx[keep]
    col_idx <- col_idx[keep]
    nzdata <- nzdata[keep]
    nrow <- length(i)
  }
  list(dim=c(nrow, ncol), row=row_idx, col=col_idx, val=nzdata)
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
.extract_array_from_CSC_ZarrSparseMatrixSeed <- function(x, index,
                                                          transpose=FALSE)
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
#' ## Genes as rows, cells as columns: transposing the seed gives a
#' ## ZarrSparseMatrix backed by a CSC seed, without reading any data
#' ZarrSparseMatrix(t(seed(X)))
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
