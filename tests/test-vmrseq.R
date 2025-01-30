# Load required libraries
library(testthat)
library(SummarizedExperiment)
library(GenomicRanges)
library(vmrseq)  # Ensure the package is loaded

# Mock data for testing
create_mock_SE <- function() {
  # Create a simple GRanges object
  gr <- GRanges(seqnames = Rle(c("chr1", "chr2"), c(3, 2)),
                ranges = IRanges(start = c(1, 200, 400, 100, 300),
                                 width = 100),
                strand = Rle(strand("*"), 5))
  
  # Create binary methylation matrix (mock assay data)
  methylation_data <- matrix(c(1, 0, 1, 0, 1, 0, 1, NA, 1, 0), nrow = 5)
  
  # Create SummarizedExperiment object
  SummarizedExperiment(assays = list(M_mat = methylation_data), rowRanges = gr)
}

# Test 1: Check that the function runs without errors on valid input
test_that("vmrseqSmooth runs without errors on valid input", {
  mock_SE <- create_mock_SE()
  result <- vmrseqSmooth(mock_SE)
  
  expect_s4_class(result, "GRanges")  # Check that the output is a GRanges object
})

# Test 2: Check that the result contains the expected metadata columns
test_that("vmrseqSmooth output contains expected metadata columns", {
  mock_SE <- create_mock_SE()
  result <- vmrseqSmooth(mock_SE)
  
  # Check for columns 'meth', 'total', and 'var' in the result metadata
  expect_true("meth" %in% names(mcols(result)))
  expect_true("total" %in% names(mcols(result)))
  expect_true("var" %in% names(mcols(result)))
})

# Test 3: Check that the function handles missing values correctly
test_that("vmrseqSmooth handles missing values correctly", {
  mock_SE <- create_mock_SE()
  assays(mock_SE)$M_mat[1, 2] <- NA  # Introduce an NA value
  
  result <- vmrseqSmooth(mock_SE, sparseNAdrop = FALSE)
  
  # Expect that missing values do not cause errors
  expect_s4_class(result, "GRanges")
  expect_true(!is.na(mcols(result)$meth[1]))  # meth column should not be NA
})

# Test 4: Check that the function throws an error when positions are too close
test_that("vmrseqSmooth throws error for CpG sites with position differences < 2", {
  gr <- GRanges(seqnames = "chr1",
                ranges = IRanges(start = c(1, 2), width = 1))  # Positions differ by 1 bp
  
  methylation_data <- matrix(c(1, 0), nrow = 2)
  mock_SE <- SummarizedExperiment(assays = list(M_mat = methylation_data), rowRanges = gr)
  
  expect_error(vmrseqSmooth(mock_SE), "position difference less than 2 bp")
})

# Test 5: Check that parallel execution is handled correctly
test_that("vmrseqSmooth handles parallel execution", {
  mock_SE <- create_mock_SE()
  
  result <- vmrseqSmooth(mock_SE, BPPARAM = BiocParallel::MulticoreParam(2), verbose = TRUE)
  
  # Check that parallelization completed without errors
  expect_s4_class(result, "GRanges")
})
