test_that("Returns the right shape and class", {
    load("data/test_data.RData")
    result <- do.call(airm_exp, args = test_pd_mats)
    expect_equal(result@Dim, c(10,10)) 
    inherits(result, "dppMatrix") |> expect_true()
})

test_that("Returns expected results when the reference point is the identity", {
    load("data/test_data.RData")
    Sigma <- diag(10) |> methods::as(object =_, Class ="dpoMatrix") |> Matrix::pack()
    Lambda <- test_pd_mats[[1]]; comp <- expm::expm(Lambda)
    result <- airm_exp(Sigma, Lambda); 
    expect_lt(Matrix::norm(comp - result, "F")/Matrix::norm(comp, "F"), 1e-3)
})

test_that("Inputs are correctly validated", {
    Sigma <- Matrix::Matrix(c(1,2,3,4), nrow = 2)
    Lambda <- Sigma |> Matrix::symmpart() |> Matrix::pack()
    Lambda2 <- Matrix::Matrix(c(1,2,3,4,5,6,7,8,9), nrow = 3) |>
        Matrix::symmpart() |> Matrix::pack()
    expect_error(airm_exp(Sigma, Lambda))
    expect_error(airm_exp(Lambda, Sigma))
    expect_error(airm_exp(Sigma, Lambda2))
})

