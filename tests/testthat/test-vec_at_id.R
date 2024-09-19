test_that("returns the right shape and class", {
    # load("data/test_data.RData")
    data("test_pd_mats")
    test_pd_mats[[1]] |> Matrix::unpack() |> Matrix::symmpart() |> Matrix::pack() |> 
    vec_at_id() |> (\(x) c(length(x), inherits(x, what = c("vector", "numeric"))))() |> 
    all() |> expect_true()
    
})

test_that("returns as expected when the input is all 1s", {
    Matrix::Matrix(rep(1, 4), ncol = 2, nrow = 2) |> Matrix::symmpart() |> 
    Matrix::pack() |> vec_at_id() |> expect_equal(object = _, expected = c(1, sqrt(2), 1))
})

test_that("Inputs are correctly validated", {
    Matrix::Matrix(1:4, nrow = 2, ncol = 2) |> vec_at_id() |> expect_error()
})
