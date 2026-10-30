# ZarrSparseMatrix

A `DelayedArray` backend for sparse matrices stored in [Zarr](https://zarr.dev).

`ZarrSparseMatrix` reads a sparse matrix saved in a Zarr group as the three 1D
arrays `data`, `indices` and `indptr`, in CSR or CSC layout. This is how
AnnData's `write_zarr()` stores its `X` matrix. The matrix stays on disk and
is read with [Rarr](https://bioconductor.org/packages/Rarr) only when needed.
It is the Zarr counterpart of `HDF5Array::H5SparseMatrix`.

Because a `ZarrSparseMatrix` is a `DelayedMatrix`, it works with any code
written for `DelayedArray` objects: subsetting, `rowSums()`/`colSums()`,
block processing, delayed arithmetic, and assays in a `SingleCellExperiment`.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Rarr", "DelayedArray", "SparseArray"))

# Install from Bioconductor
BiocManager::install("ZarrSparseMatrix")

# Install from GitHub
# remotes::install_github("GabrielHoffman/ZarrSparseMatrix")
```

## Usage

```r
library(ZarrSparseMatrix)

# Example store: 200 cells x 100 genes, in AnnData's CSR layout
store <- system.file("extdata", "example.zarr", package = "ZarrSparseMatrix")
X <- ZarrSparseMatrix(store, "X")
X
#> <200 x 100> sparse ZarrSparseMatrix object of type "integer":

X[1:5, 1:10]          # delayed subsetting
colSums(X)            # computed block by block
as(X, "dgCMatrix")    # load into memory

# Genes as rows: transposing the seed turns CSR into CSC without reading data
counts <- ZarrSparseMatrix(t(seed(X)))

library(SingleCellExperiment)
sce <- SingleCellExperiment(assays = list(counts = counts))
```

The dimensions and layout are read from the group's `shape` and
`encoding-type` attributes. If they are missing, pass `dim` and
`sparse.layout`. Row and column names can be given with `dimnames`.

See the vignette for more examples:

```r
vignette("ZarrSparseMatrix")
```

## Writing a store from Python

Use Zarr format 2. Rarr reads Zarr v2 arrays, and recent anndata versions are
switching their default to v3:

```python
import anndata as ad

ad.settings.zarr_write_format = 2
adata.write_zarr("data.zarr")
```

```r
X <- ZarrSparseMatrix("data.zarr", "X")
```

## Limitations

- Only Zarr v2 stores have been tested.
- Only local stores have been tested.
- `Rarr::write_zarr_array()` writes the last chunk of an array shorter than
  the other chunks. Rarr's own reader then fails on indexed reads that touch
  that chunk. Stores written by zarr-python and anndata are not affected.
  `inst/scripts/make-example-data.R` shows how to write a compatible store
  from R.

## License

Artistic-2.0
