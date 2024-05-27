# create_seurat_object.R

# Install and load necessary packages
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
if (!requireNamespace("hdf5r", quietly = TRUE)) {
  install.packages("hdf5r")
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  install.packages("Matrix")
}

library(Seurat)
library(hdf5r)
library(Matrix)

# Define the path to your HDF5 file
h5_file = "/home/student2/seurat/file.h5"

# Open the HDF5 file
h5 = H5File$new(h5_file, mode = "r")

# Read the data from the HDF5 file
data = h5[["X/data"]]$read()
indices = h5[["X/indices"]]$read()
indptr = h5[["X/indptr"]]$read()
var = h5[["var"]]$read()

# Ensure data is in the correct format
data = as.numeric(data)
indices = as.integer(indices)
indptr = as.integer(indptr)

# Get the dimensions of the matrix
num_genes = nrow(var)
num_cells = length(indptr) - 1

# Adjust indices for 0-based indexing used in HDF5
i = indices + 1
p = indptr
x = data

# Create a sparse matrix
sparse_matrix = sparseMatrix(i = i, p = p, x = x, dims = c(num_genes, num_cells))

# Create the Seurat object
pbmc = CreateSeuratObject(counts = sparse_matrix, project = "pbmc_h5", min.cells = 3, min.features = 200)

# Print the Seurat object
print(pbmc)

# Close the HDF5 file
h5$close()
