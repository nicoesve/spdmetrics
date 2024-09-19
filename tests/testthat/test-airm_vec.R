test_that("Method returns right length and class", {
    load("data/test_data.RData")
    w <- test_pd_mats |> do.call(airm_vec, args = _)
    w |> inherits(x = _, what=c("vector", "numeric")) |> expect_true()
    w |> length() |> expect_equal(object=_, 55)
})

test_that("arguments are correctly validated", {
    load("data/test_data.RData")
    Sigma2 <- matrix(1:100, ncol = 10)
    Lambda2 <- Matrix::Matrix(rep(1, 4), nrow = 2) |> Matrix::symmpart() |> 
	Matrix::pack()
    airm_vec(test_pd_mats[[1]], Sigma2) |> expect_error()
    airm_vec(Sigma2, test_pd_mats[[1]]) |> expect_error()
    airm_vec(Sigma, Lambda2) |> expect_error()
})

test_that("returns reasonable values in simple cases", {
    load("data/test_data.RData")
    list(list(diag(10) |> Matrix::nearPD() |> _$mat |> Matrix::pack(), 
	    test_pd_mats[[1]] |> Matrix::unpack() |> 
		methods::as(object=_, Class="dsyMatrix") |> Matrix::pack()) |> 
	do.call(airm_vec, args=_),
	vec_at_id(test_pd_mats[[1]])) |> 
	(\(l) norm((l[[1]] - l[[2]]) |> as.matrix(x = _, ncol=1), "F")/
		Matrix::norm(l[[2]] |> as.matrix(x = _, ncol=1), "F"))() |>
	expect_lt(object=_, expected=1e-3)
})

test_that("vec and unvec are mutual inverses", {
    load("data/test_data.RData")
    list(test_pd_mats[[1]], 
	test_pd_mats[[2]] |> Matrix::unpack() |> 
	    methods::as(object =_, Class="dsyMatrix") |> Matrix::pack()) |> 
    (\(l) list(l, do.call(airm_vec, args = l)))() |>  
    (\(l) list(l[[1]][[2]], airm_unvec(l[[1]][[1]], l[[2]])))() |>
    (\(l) Matrix::norm(l[[1]] - l[[2]], "F")/Matrix::norm(l[[1]], "F"))() |>
    expect_lt(object=_, expected=1e-3)

    list(test_pd_mats[[1]], rep(1, 55)) |> 
    (\(l) list(l, do.call(airm_unvec, args = l)))() |> 
    (\(l) list(l[[1]][[2]], airm_vec(l[[1]][[1]], l[[2]])))() |>
    (\(l) norm((l[[1]] - l[[2]]) |> as.matrix(x = _, ncol=1))/
	norm(l[[1]] |> as.matrix(x =_, ncol=1)))() |>
    expect_lt(object=_, expected=1e-3)
})
