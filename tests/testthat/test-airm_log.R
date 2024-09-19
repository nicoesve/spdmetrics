test_that("Returns the right shape and class", {
    # load("data/test_data.RData")
    data("test_pd_mats")
    result <- do.call(airm_log, args = test_pd_mats)
    expect_equal(result@Dim, c(10,10)) 
    inherits(result, "dspMatrix") |> expect_true()
})

test_that("Returns expected results when the reference point is the identity", {
    # load("data/test_data.RData")
    data("test_pd_mats")
    Sigma <- diag(10) |> methods::as(object =_, Class ="dpoMatrix") |> Matrix::pack()
    Lambda <- test_pd_mats[[1]]; comp <- expm::logm(Lambda)
    result <- airm_log(Sigma, Lambda); 
    expect_lt(Matrix::norm(comp - result, "F")/Matrix::norm(comp, "F"), 1e-3)
})

test_that("Riemannian log and exp are mutual inverses", {
    # load("data/test_data.RData")
    data("test_pd_mats")
    comp1 <- airm_log(test_pd_mats[[1]], airm_exp(test_pd_mats[[1]], test_pd_mats[[2]]))
    comp2 <- airm_exp(test_pd_mats[[1]], airm_log(test_pd_mats[[1]], test_pd_mats[[2]]))
    list(comp1, comp2) |> purrr::map(.x= _, .f = \(x) Matrix::norm(test_pd_mats[[2]] - x, "F")) |>
        purrr::walk(\(x) expect_lt(object = x/Matrix::norm(test_pd_mats[[2]], "F"), 1e-3))
})

test_that("Inputs are correctly validated", {
    Sigma <- Matrix::Matrix(c(1,2,3,4), nrow = 2)
    Lambda <- Sigma |> Matrix::nearPD() |> _$mat |> Matrix::pack()
    Lambda2 <- Matrix::Matrix(c(1,2,3,4,5,6,7,8,9), nrow = 3) |> 
        Matrix::nearPD() |> _$mat |> Matrix::pack()
    expect_error(airm_log(Sigma, Lambda))
    expect_error(airm_log(Lambda, Sigma))
    expect_error(airm_log(Sigma, Lambda2))
})
    




