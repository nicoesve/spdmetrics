test_that("Method returns right length and class", {
    load("data/test_data.RData")
    Lambda <- test_pd_mats |> 
	(\(l) list(l[[1]], l[[2]] |> Matrix::unpack() |> 
	    methods::as(object=_, "dsyMatrix") |> Matrix::pack()))() |>
	(\(l) list(l, do.call(airm_vec, args = l)))() |> 
	(\(l) airm_unvec(l[[1]][[1]], l[[2]]))()
    Lambda |> inherits(x = _, what="dspMatrix") |> expect_true()
    Lambda |> _@Dim |> expect_equal(object=_, c(10, 10))
})

test_that("arguments are correctly validated", {
    load("data/test_data.RData")
    Sigma2 <- matrix(1:100, ncol = 10)
    w2 <- 1:55
    Lambda2 <- Matrix::Matrix(rep(1, 4), nrow = 2) |> Matrix::symmpart() |> 
	Matrix::pack()
    airm_unvec(test_pd_mats[[1]], Sigma2) |> expect_error()
    airm_unvec(Sigma2, w2) |> expect_error()
    airm_unvec(Lambda2, w2) |> expect_error()
}) 
