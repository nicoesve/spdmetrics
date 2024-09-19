test_that("Initialization using connectomes works", {
    # load("data/test_data.RData")
    # load("data/airm.RData")
	data("test_pd_mats")
	data("airm")
    CSample$new(conns = test_pd_mats, geom = airm) |>  
	(\(s) list(expect_identical(object = s$connectomes, expected = test_pd_mats),
	    expect_equal(s$sample_size, 2),
	    expect_equal(s$matrix_size, 10),
	    expect_equal(s$mfd_dim, 55), 
	    list(s$frechet_mean, s$tangent_images, s$vector_images, s$is_centered,
		s$variation, s$sample_cov) |> 
		purrr::map_lgl(.x =_, .f = is.null) |> all() |> expect_true()))()
})

 test_that("Initialization using tangent images works", {
    # load("data/test_data.RData")
    # load("data/airm.RData")
	data("test_pd_mats")
	data("airm")
    ref <- diag(10) |> Matrix::nearPD() |> _$mat |> Matrix::pack()
    ts <- test_pd_mats |> purrr::map(.x =_, .f = \(x) airm$log(ref, Lambda=x))
    list(ref, ts) |> CSample$new(tan_imgs = _, geom = airm, centered = FALSE) |>  
	(\(s) list(expect_identical(object = s$tangent_images, expected = list(ref, ts)),
	    expect_equal(s$sample_size, 2),
	    expect_equal(s$matrix_size, 10),
	    expect_equal(s$mfd_dim, 55), 
	    list(s$frechet_mean, s$connectomes, s$vector_images,
		s$variation, s$sample_cov) |> 
		purrr::map_lgl(.x =_, .f = is.null) |> all() |> expect_true()))()
})

test_that("Initialization using vecs works", {
    # load("data/test_data.RData")
    # load("data/airm.RData")
	data("test_pd_mats")
	data("airm")
    ref <- diag(10) |> Matrix::nearPD() |> _$mat |> Matrix::pack()
    vs <- test_pd_mats |> purrr::map(.x =_, .f = \(x) airm$log(ref, Lambda=x)) |>
	purrr::map(.x =_, .f = \(x) airm$vec(ref, v=x)) |>
	Reduce(rbind, x =_)
    list(ref, vs) |> CSample$new(vec_imgs = _, geom = airm, centered = FALSE) |>  
	(\(s) list(expect_identical(object = s$vector_images, expected = list(ref, vs)),
	    expect_equal(s$sample_size, 2),
	    expect_equal(s$matrix_size, 10),
	    expect_equal(s$mfd_dim, 55), 
	    list(s$frechet_mean, s$tangent_images, s$connectomes,
		s$variation, s$sample_cov) |> 
		purrr::map_lgl(.x =_, .f = is.null) |> all() |> expect_true()))()
})
 
