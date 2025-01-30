CSample <- R6::R6Class(
    classname = "CSample",
    private = list(
        conns = NULL, vec_imgs = NULL,
        n = NULL, p = NULL, d = NULL, centered = NULL,
        f_mean = NULL, metric = NULL, var = NULL, s_cov = NULL,
        tangent_handler = NULL
    ),
    public = list(
        initialize = function(conns = NULL, tan_imgs = NULL,
                              vec_imgs = NULL, centered = NULL, metric) {
            validate_metric(metric)
            private$metric <- metric
            private$tangent_handler <- TangentImageHandler$new(metric)

            if (!is.null(conns)) {
                validate_conns(conns, tan_imgs, vec_imgs, centered)

                n <- length(conns)
                p <- nrow(conns[[1]])
                d <- p * (p + 1) / 2
                private$conns <- conns
                private$vec_imgs <- NULL
                private$n <- n
                private$p <- p
                private$d <- d
                private$centered <- NULL
                private$f_mean <- NULL
                private$var <- NULL
                private$s_cov <- NULL
            } else if (!is.null(tan_imgs)) {
                validate_tan_imgs(tan_imgs, vec_imgs, centered)

                n <- length(tan_imgs[[2]])
                p <- nrow(tan_imgs[[1]])
                d <- p * (p + 1) / 2
                frechet_mean <- if (centered) tan_imgs[[1]] else NULL
                private$conns <- NULL
                private$vec_imgs <- NULL
                private$n <- n
                private$p <- p
                private$d <- d
                private$centered <- centered
                private$f_mean <- frechet_mean
                private$var <- NULL
                private$s_cov <- NULL
                private$tangent_handler$set_tangent_images(
                    tan_imgs[[1]], tan_imgs[[2]]
                )
            } else {
                validate_vec_imgs(vec_imgs, centered)

                n <- nrow(vec_imgs[[2]])
                p <- nrow(vec_imgs[[1]])
                d <- p * (p + 1) / 2

                dims_flag <- d == (vec_imgs[[2]] |> ncol())
                if (!dims_flag) stop("Dimensions don't match")

                frechet_mean <- if (centered) vec_imgs[[1]] else NULL

                private$conns <- NULL
                private$vec_imgs <- vec_imgs
                private$n <- n
                private$p <- p
                private$d <- d
                private$centered <- centered
                private$f_mean <- frechet_mean
                private$var <- NULL
                private$s_cov <- NULL
            }
            default_ref_pt <- diag(p) |>
                methods::as("dpoMatrix") |>
                Matrix::pack()
        },
        compute_tangents = function(ref_pt = default_ref_pt) {
            if (!inherits(ref_pt, "dppMatrix")) {
                stop("ref_pt must be a dppMatrix object.")
            }
            if (is.null(private$conns)) stop("conns must be specified.")
            private$tangent_handler$set_reference_point(ref_pt)
            private$tangent_handler$compute_tangents(private$conns)
        },
        compute_conns = function() {
            if (is.null(private$tangent_handler$tangent_images)) {
                stop("tangent images must be specified.")
            }
            private$conns <- private$tangent_handler$compute_conns()
        },
        compute_vecs = function() {
            if (is.null(private$tangent_handler$tangent_images)) {
                stop("tangent images must be specified.")
            }
            private$vec_imgs <- private$tangent_handler$compute_vecs()
        },
        compute_unvecs = function() {
            if (is.null(private$vec_imgs)) stop("vec_imgs must be specified.")
            private$tangent_handler$set_tangent_images(
                private$vec_imgs[[1]],
                1:nrow(private$vec_imgs[[2]]) |>
                    purrr::map(\(i) {
                        private$metric$unvec(
                            private$vec_imgs[[1]],
                            private$vec_imgs[[2]][i, ]
                        )
                    })
            )
        },
        compute_fmean = function(tol = 0.05, max_iter = 20, lr = 0.2) {
            compute_fmean(self, tol, max_iter, lr)
        },
        change_ref_pt = function(new_ref_pt) {
            if (is.null(private$tangent_handler$tangent_images)) {
                stop("tangent images have not been computed")
            }
            if (!inherits(new_ref_pt, "dppMatrix")) {
                stop("new_ref_pt must be a dppMatrix object.")
            }
            private$tangent_handler$relocate_tangents(new_ref_pt)
        },
        center = function() {
            if (is.null(private$tangent_handler$tangent_images)) {
                stop("tangent images must be specified.")
            }
            if (!(private$centered |> is.null())) {
                if (private$centered) stop("The sample is already centered.")
            }
            if (is.null(private$f_mean)) self$compute_fmean()
            private$centered <- TRUE
            self$change_ref_pt(private$f_mean)
        },
        compute_variation = function() {
            if (self$vector_images |> is.null()) {
                if (private$tangent_handler$tangent_images |> is.null()) {
                    self$compute_tangents()
                }
                self$compute_vecs()
            }

            if (is.null(private$vec_imgs)) stop("vec_imgs must be specified.")
            if (is.null(private$centered) || !private$centered) {
                self$compute_unvecs()
                self$center()
                self$compute_vecs()
            }
            private$var <- private$vec_imgs[[2]] |>
                apply(X = _, MARGIN = 1, FUN = function(x) sum(x^2)) |>
                mean()
        },
        compute_sample_cov = function() {
            if (self$vector_images |> is.null()) {
                if (private$tangent_handler$tangent_images |> is.null()) {
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
        tangent_images = function() private$tangent_handler$tangent_images,
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
