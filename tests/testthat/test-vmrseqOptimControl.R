library(testthat)

# Test 1: Check that the function returns a list with the correct structure
test_that("vmrseqOptimControl returns a list with the correct structure", {
  result <- vmrseqOptimControl()
  
  # Check that the result is a list
  expect_type(result, "list")
  
  # Check that the list contains all expected elements
  expect_named(result, c("inits", "epsilon", "backtrack", "eta", "maxIter"))
  
  # Check types of each element
  expect_type(result$inits, "double")    # inits should be numeric
  expect_type(result$epsilon, "double")  # epsilon should be numeric
  expect_type(result$backtrack, "logical")  # backtrack should be logical
  expect_type(result$eta, "double")      # eta should be numeric
  expect_type(result$maxIter, "double")  # maxIter should be numeric
})

# Test 2: Check default values
test_that("vmrseqOptimControl returns correct default values", {
  result <- vmrseqOptimControl()
  
  expect_equal(result$inits, c(0.2, 0.5, 0.8))
  expect_equal(result$epsilon, 1e-3)
  expect_true(result$backtrack)
  expect_equal(result$eta, 0.05)  # Default eta when backtrack = TRUE
  expect_equal(result$maxIter, 100)
})

# Test 3: Check that custom values are used correctly
test_that("vmrseqOptimControl uses custom parameter values", {
  custom_inits <- c(0.1, 0.4, 0.9)
  custom_epsilon <- 1e-5
  custom_backtrack <- FALSE
  custom_eta <- 0.01
  custom_maxIter <- 200
  
  result <- vmrseqOptimControl(
    inits = custom_inits,
    epsilon = custom_epsilon,
    backtrack = custom_backtrack,
    eta = custom_eta,
    maxIter = custom_maxIter
  )
  
  expect_equal(result$inits, custom_inits)
  expect_equal(result$epsilon, custom_epsilon)
  expect_false(result$backtrack)
  expect_equal(result$eta, custom_eta)
  expect_equal(result$maxIter, custom_maxIter)
})

# Test 4: Check that eta defaults based on backtrack
test_that("vmrseqOptimControl adjusts eta based on backtrack", {
  # When backtrack is TRUE, eta should be 0.05
  result_with_backtrack <- vmrseqOptimControl(backtrack = TRUE)
  expect_equal(result_with_backtrack$eta, 0.05)
  
  # When backtrack is FALSE, eta should default to 0.005
  result_without_backtrack <- vmrseqOptimControl(backtrack = FALSE)
  expect_equal(result_without_backtrack$eta, 0.005)
})

# Test 5: Check that inits values outside the range (0, 1) trigger an error
test_that("vmrseqOptimControl throws an error for invalid inits values", {
  # Create invalid inits vectors
  invalid_inits_low <- c(-0.1, 0.5, 0.8)     # Contains a value < 0
  invalid_inits_high <- c(0.2, 0.5, 1.1)      # Contains a value > 1
  invalid_inits_edge <- c(0, 0.5, 0.8)        # Contains a value exactly at 0
  
  # Expect an error when inits values are outside the valid range (0, 1)
  expect_error(vmrseqOptimControl(inits = invalid_inits_low), 
               "All values in inits has to between 0 and 1!")
  
  expect_error(vmrseqOptimControl(inits = invalid_inits_high), 
               "All values in inits has to between 0 and 1!")
  
  expect_error(vmrseqOptimControl(inits = invalid_inits_edge), 
               "All values in inits has to between 0 and 1!")
  
  # Valid inits should not throw an error
  valid_inits <- c(0.1, 0.5, 0.9)
  expect_silent(vmrseqOptimControl(inits = valid_inits))
})
