CSample <- R6::R6Class(
    classname = "CSample",
    private = list(
        conns = NULL, tan_imgs = NULL, vec_imgs = NULL, 
        n = NULL, p = NULL, d = NULL, centered = NULL, 
        f_mean = NULL, metric = NULL, var = NULL, s_cov = NULL
    ),
    public = list(
        initialize = function(conns = NULL, tan_imgs = NULL, 
                              vec_imgs = NULL, centered = NULL, metric) {
            if (is.null(metric)) stop("metric must be specified.")

            if (!is.null(conns)) {
                if (!is.null(tan_imgs) || !is.null(vec_imgs)) {
                    stop("When initializing, if conns is not NULL, 
                          tan_imgs and vec_imgs must be NULL.")
                }

                if (!is.null(centered)) { 
                    warning("If conns is not NULL, centered is ignored")
                }

                class_flag <- conns |> 
                    purrr::map_lgl(\(x) inherits(x, "dppMatrix")) |> all()
                if (!class_flag) stop("conns must be a list of dppMatrix objects.")

                n <- length(conns); p <- nrow(conns[[1]]); d <- p * (p + 1) / 2
                private$conns <- conns; private$tan_imgs <- NULL
                private$vec_imgs <- NULL; private$n <- n; private$p <- p
                private$d <- d; private$centered <- NULL; private$f_mean <- NULL
                private$metric <- metric; private$var <- NULL; private$s_cov <- NULL

            } else if (!is.null(tan_imgs)) {
                if (!is.null(vec_imgs)) {
                    stop("If tan_imgs is not NULL, conns and vec_imgs must be NULL.")
                }
                if (is.null(centered)) {
                    stop("If tan_imgs is not NULL, centered must be specified.")
                }
                if (!is.logical(centered)) stop("centered must be a logical.")
                if (!inherits(tan_imgs[[1]], "dppMatrix")) {
                    stop("The first element of tan_imgs must be a dppMatrix object.")
                }
                if (!is.list(tan_imgs[[2]])) {
                    stop("The second element of tan_imgs must be a list.")
                }
                class_flag <- tan_imgs[[2]] |> 
                    purrr::map_lgl(\(x) inherits(x, "dspMatrix")) |> all()
                if (!class_flag) {
                    stop("The second element of tan_imgs must be a list of dspMatrix objects.")
                }

                n <- length(tan_imgs[[2]]); p <- nrow(tan_imgs[[1]]); d <- p * (p + 1) / 2
                frechet_mean <- if (centered) tan_imgs[[1]] else NULL
                private$conns <- NULL; private$tan_imgs <- tan_imgs; private$vec_imgs <- NULL
                private$n <- n; private$p <- p; private$d <- d; private$centered <- centered
                private$f_mean <- frechet_mean; private$metric <- metric; private$var <- NULL
                private$s_cov <- NULL
            } else {
                if (is.null(vec_imgs)) {
                    stop("At least one of conns, tan_imgs, or vec_imgs must be specified.")
                }
                if (is.null(centered)) {
                    stop("If vec_imgs is not NULL, centered must be specified.")
                }
                if (!is.logical(centered)) stop("centered must be a logical.")
                if (!inherits(vec_imgs[[1]], "dppMatrix")) {
                    stop("The first element of vec_imgs must be a dppMatrix object.")
                }
                if (!is.matrix(vec_imgs[[2]])) {
                    stop("The second element of vec_imgs must be a matrix.")
                }

                n <- nrow(vec_imgs[[2]])
                p <- nrow(vec_imgs[[1]])
                d <- p * (p + 1) / 2

                dims_flag <- d == (vec_imgs[[2]] |> ncol())
                if (!dims_flag) stop("Dimensions don't match")

                frechet_mean <- if (centered) vec_imgs[[1]] else NULL

                private$conns <- NULL; private$tan_imgs <- NULL
                private$vec_imgs <- vec_imgs; private$n <- n; private$p <- p
                private$d <- d; private$centered <- centered; 
                private$f_mean <- frechet_mean; private$metric <- metric; 
                private$var <- NULL; private$s_cov <- NULL
            }
        },
        compute_tangents = function(
            ref_pt = diag(private$p) |> Matrix::nearPD() |> 
                _$mat |> Matrix::pack()) {
            if (!inherits(ref_pt, "dppMatrix")) { 
                stop("ref_pt must be a dppMatrix object.")
            }
            if (is.null(private$conns)) stop("conns must be specified.")
            private$tan_imgs <- private$conns |> 
                purrr::map(\(conn) private$metric$log(ref_pt, conn)) |>
                (\(x) list(ref_pt, x))()
        },
        compute_conns = function() {
            if (is.null(private$tan_imgs)) stop("tan_imgs must be specified.")
            private$conns <- private$tan_imgs[[2]] |> 
                purrr::map(\(tan) private$metric$exp(private$tan_imgs[[1]], tan))
        },
        compute_vecs = function() {
            if (is.null(private$tan_imgs)) stop("tan_imgs must be specified.")
            private$vec_imgs <- private$tan_imgs[[2]] |> 
                purrr::map(\(tan) self$metric$vec(private$tan_imgs[[1]], tan)) |> 
                do.call(rbind, args = _) |> (\(x) list(private$tan_imgs[[1]], x))()
        },
        compute_unvecs = function() {
            if (is.null(private$vec_imgs)) stop("vec_imgs must be specified.")
            private$tan_imgs <- 1:nrow(private$vec_imgs[[2]]) |> 
                purrr::map(\(i) self$metric$unvec(private$vec_imgs[[1]], 
                        private$vec_imgs[[2]][i, ])) |> 
                (\(x) list(private$vec_imgs[[1]], x))()
        },
        compute_fmean = function(tol = 0.05, max_iter = 20, lr = 0.2) {
            if (!is.null(private$f_mean)) { 
                warning("The Frechet mean has already been computed.")
            }
            if (is.null(private$tan_imgs)) { 
                message("tangent images were null, so they will be computed")
                self$compute_tangents()
            }
            if (!is.numeric(tol)) stop("tol must be a numeric.")
            if (max_iter < 1) stop("max_iter must be at least 1.")

            aux_sample <- self; delta <- Inf; iter <- 0

            while ((delta > tol) && (iter < max_iter)) {
                old_tan <- aux_sample$tangent_images; iter <- iter + 1 
                old_ref_pt <- old_tan[[1]]

                if (iter > max_iter) {
                    warning("Computation of Frechet mean exceeded maximum 
                        number of iterations.")
                }

                tan_step <- lr * Reduce(`+`, old_tan[[2]]) / 
                    aux_sample$sample_size
                new_ref_pt <- self$metric$exp(old_ref_pt, tan_step)

                delta <- Matrix::norm(new_ref_pt - old_ref_pt, "F") / 
                        Matrix::norm(old_ref_pt, "F")

                new_tan_imgs <- relocate(old_ref_pt, new_ref_pt, old_tan[[2]], 
                    self$metric) 

                aux_sample <- CSample$new(
                    tan_imgs = list(new_ref_pt, new_tan_imgs),
                    centered = FALSE, metric = self$metric)
            }
            private$f_mean <- aux_sample$tangent_images[[1]]
        },
        change_ref_pt = function(new_ref_pt) {
            if (is.null(private$tan_imgs)) {
                stop("tangent images have not been computed")
            }
            if (!inherits(new_ref_pt, "dppMatrix")) { 
                stop("new_ref_pt must be a dppMatrix object.")
            }
            new_tan_imgs <- relocate(private$tan_imgs[[1]], new_ref_pt, 
                private$tan_imgs[[2]], self$metric)
            private$tan_imgs <- list(new_ref_pt, new_tan_imgs)
        },
        center = function() {
            if (is.null(private$tan_imgs)) stop("tan_imgs must be specified.")
            if (!(private$centered |> is.null())) {
                if (private$centered) stop("The sample is already centered.")
            }
            if (is.null(private$f_mean)) self$compute_fmean()
            private$centered <- TRUE; self$change_ref_pt(private$f_mean)
        },
        compute_variation = function() {
            if (self$vector_images |> is.null()) {
                if (self$tangent_images |> is.null()) { 
                    self$compute_tangents()
                } 
                self$compute_vecs()
            }

            if (is.null(private$vec_imgs)) stop("vec_imgs must be specified.")
            if (is.null(private$centered) || !private$centered) {
                self$compute_unvecs(); self$center(); self$compute_vecs()
            }
            private$var <- private$vec_imgs[[2]] |>
                apply(X = _, MARGIN = 1, FUN = function(x) sum(x^2)) |> mean()
        },
        compute_sample_cov = function() {
            if (self$vector_images |> is.null()) {
                if (self$tangent_images |> is.null()) { 
                    self$compute_tangents()
                } 
                self$compute_vecs()
            } 

            if (is.null(private$vec_imgs)) stop("vec_imgs must be specified.")
            private$s_cov <- cov(private$vec_imgs[[2]])
        }
    ),
    active = list(
        connectomes = function() private$conns,
        tangent_images = function() private$tan_imgs,
        vector_images = function() private$vec_imgs,
        sample_size = function() private$n,
        matrix_size = function() private$p,
        mfd_dim = function() private$d,
        is_centered = function() private$centered,
        frechet_mean = function() private$f_mean,
        metric = function() private$metric,
        variation = function() private$var,
        sample_cov = function() private$s_cov
    )
)

