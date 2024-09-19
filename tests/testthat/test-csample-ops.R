test_that("computation of tangent images works", {
    # load("data/test_data.RData")
    # load("data/airm.RData")
	data("test_pd_mats")
	data("airm")
    s <- CSample$new(conns = test_pd_mats, geom = airm)
	s$compute_tangents() 
	s$tangent_images |> (\(x) list(x |> is.null() |> expect_false(),
	    x |> inherits(x = _, "list") |> expect_true(),
	    x |> length() |> expect_equal(object=_, 2)))()
	s$tangent_images[[1]] |> (\(m) list(
			m |> inherits(x=_, "dppMatrix") |> expect_true(),
			m@Dim |> expect_equal(object=_, c(10,10))))()
	s$tangent_images[[2]] |> (\(l) list(
			l |> length() |> expect_equal(object=_, 2),
			l |> purrr::map_lgl(.x = _, .f= \(m) (m@Dim == c(10,10)) |> all()) |>
			all() |> expect_true()))()
})

test_that("computation of vectorized images works", {
	# load("data/test_data.RData")
    # load("data/airm.RData")
	data("test_pd_mats")
	data("airm")
    s <- CSample$new(conns = test_pd_mats, geom = airm)
	s$compute_tangents() 
	s$compute_vecs()
	s$vector_images[[1]] |> (\(m) list(
		m |> inherits(x=_, "dppMatrix") |> expect_true(),
		m |> _@Dim |> (\(v) v == c(10,10))() |> all() |> expect_true()
	))()
	s$vector_images[[2]] |> (\(m) list(
		m |> inherits(x=_, "matrix") |> expect_true(),
		m |> nrow() |> expect_equal(object=_, 2),
		m |> ncol() |> expect_equal(object=_, 55),
		m |> as.vector() |> is.null() |> any() |> expect_false()
	))()
})

test_that("inverse of vectorization works", {
	# load("data/test_data.RData")
    # load("data/airm.RData")
	data("test_pd_mats")
	data("airm")
    s <- CSample$new(conns = test_pd_mats, geom = airm)
	s$compute_tangents() 
	old_tan <- s$tangent_images
	s$compute_vecs()
	s$compute_unvecs()
	expect_equal(s$tangent_images, old_tan)
})	    	     

test_that("changing reference points works", {
    # load("data/test_data.RData")
    # load("data/airm.RData")
	data("test_pd_mats")
	data("airm")
    s <- CSample$new(conns = test_pd_mats, geom = airm)
    s$compute_tangents()
    new_ref <- 2*diag(10) |> Matrix::nearPD() |> _$mat |> Matrix::pack()
    new_ref |> s$change_ref_pt()
    s$tangent_images[[1]] |> expect_equal(object=_, new_ref)
    s$tangent_images[[2]] |> (\(l) list(l |> length() |> expect_equal(object=_, 2),
		l |> purrr::map_lgl(.x =_, .f = \(x) inherits(x, "dspMatrix")) |> 
			all() |> expect_true(),
		l |> purrr::map_lgl(.x = _, .f = \(x) (x@Dim == c(10,10)) |> all()) |> 
			all() |> expect_true()))()
})

test_that("computation of frechet mean works", {
    # load("data/test_data.RData")
    # load("data/airm.RData")
	data("test_pd_mats")
	data("airm")
    s <- CSample$new(conns = test_pd_mats, geom = airm) 
	s$compute_fmean() |> suppressWarnings()
	s$frechet_mean |> (\(m) list(
		m |> is.null() |> expect_false(),
	    m |> inherits(x = _, "dppMatrix") |> expect_true(),
	    m |> _@Dim |> expect_equal(object=_, c(10,10))
	))()
})

test_that("centering works", {
	# load("data/test_data.RData")
    # load("data/airm.RData")
	data("test_pd_mats")
	data("airm")
    s <- CSample$new(conns = test_pd_mats, geom = airm)
	s$center() |> expect_error()
	s$compute_tangents(); s$center()
	s$is_centered |> is.null() |> expect_false()
	s$compute_vecs()
	s$vector_images[[2]] |> apply(X=_, 2, mean) |> abs() |> 
		(\(r) r < 0.2)() |> all() |> expect_true()
})

test_that("computing variation works", {
	# load("data/test_data.RData")
    # load("data/airm.RData")
	data("test_pd_mats")
	data("airm")
    s <- CSample$new(conns = test_pd_mats, geom = airm)
	s$compute_variation() 
	s$variation |> (\(x) list(
		x |> is.null() |> expect_false(),
		x |> inherits(x=_, "numeric") |> expect_true(),
		x |> expect_gt(object=_, 0) 
	))()
})

test_that("computing sample covariance works", {
	# load("data/test_data.RData")
    # load("data/airm.RData")
	data("test_pd_mats")
	data("airm")
    s <- CSample$new(conns = test_pd_mats, geom = airm)
	s$compute_sample_cov()
	s$sample_cov |> (\(x) list(
		x |> is.null() |> expect_false(),
		x |> inherits(x=_, "matrix") |> expect_true(),
		x |> as.vector() |> is.null() |> any() |> expect_false(),
		x |> isSymmetric() |> expect_true()
	))()
})
    