#' CSample Class
#'
#' @description
#' This class represents a sample of connectomes, with various properties and methods to handle their tangent and vectorized images. # nolint: line_length_linter.
#'
#' @field conns A list of connectomes.
#' @field vec_imgs A list of vector images.
#' @field n Number of connectomes.
#' @field p Number of rows (equivalently, of columns) in each connectome.
#' @field d Dimensionality of the tangent space of connectomes.
#' @field centered Whether the data is centered.
#' @field f_mean Frechet mean of the sample.
#' @field metric Riemannian metric used.
#' @field var Variation of the data.
#' @field s_cov Sample covariance of the data.
#' @field tangent_handler Handler for tangent images.
#'
#' @export
CSample <- R6::R6Class( # nolint: cyclocomp_linter
    classname = "CSample",
    private = list(
        conns = NULL, vec_imgs = NULL,
        n = NULL, p = NULL, d = NULL, centered = NULL,
        f_mean = NULL, metric = NULL, var = NULL, s_cov = NULL,
        tangent_handler = NULL
    ),
    public = list(
        #' Initialize a CSample object
        #'
        #' @param conns A list of connectomes (default is NULL).
        #' @param tan_imgs A list of tangent images (default is NULL).
        #' @param vec_imgs A list of vectorized images (default is NULL).
        #' @param centered Boolean indicating whether tangent or vectorized images are centered (default is NULL). # nolint: line_length_linter
        #' @param metric Object of class `rmetric` representing the Riemannian metric used. # nolint: line_length_linter
        #'
        #' @return A new `CSample` object.
        #' @examples
        #' data(airm)
        #' conn_sample <- CSample$new(
        #' conns = list_of_conns,
        #' metric = airm
        #' )
        initialize = function(conns = NULL, tan_imgs = NULL,
                              vec_imgs = NULL, centered = NULL, metric) {
            # Validate and set the metric
            validate_metric(metric)
            private$metric <- metric
            private$tangent_handler <- TangentImageHandler$new(metric)

            # If connectomes are provided
            if (!is.null(conns)) {
                validate_conns(conns, tan_imgs, vec_imgs, centered)

                # Set dimensions and initialize properties
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

                # If tangent images are provided
            } else if (!is.null(tan_imgs)) {
                validate_tan_imgs(tan_imgs, vec_imgs, centered)

                # Set dimensions and initialize properties
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

                # Set tangent images in the handler
                private$tangent_handler$set_tangent_images(
                    tan_imgs[[1]], tan_imgs[[2]]
                )

                # If vector images are provided
            } else {
                validate_vec_imgs(vec_imgs, centered)

                # Set dimensions and initialize properties
                n <- nrow(vec_imgs[[2]])
                p <- nrow(vec_imgs[[1]])
                d <- p * (p + 1) / 2

                # Check if dimensions match
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

            # Set default reference point
            default_ref_pt <- diag(p) |>
                methods::as("dpoMatrix") |>
                Matrix::pack()
        },

        #' Compute Tangent Images
        #'
        #' @description
        #' This function computes the tangent images from the connectomes.
        #'
        #' @param ref_pt A reference point, which must be a `dppMatrix` object (default is `default_ref_pt`). # nolint: line_length_linter
        #'
        #' @return None
        #' @throws Error if `ref_pt` is not a `dppMatrix` object or if `conns` is not specified. # nolint: line_length_linter
        #' @examples
        #' \dontrun{
        #'   CSample$compute_tangents(ref_pt = some_ref_pt)
        #' }
        compute_tangents = function(ref_pt = default_ref_pt) {
            if (!inherits(ref_pt, "dppMatrix")) {
                stop("ref_pt must be a dppMatrix object.")
            }
            if (is.null(private$conns)) stop("conns must be specified.")
            private$tangent_handler$set_reference_point(ref_pt)
            private$tangent_handler$compute_tangents(private$conns)
        },

        #' Compute connectomes
        #'
        #' @description
        #' This function computes the connectomes from the tangent images.
        #'
        #' @return None
        #' @throws Error if tangent images are not specified.
        #' @examples
        #' \dontrun{
        #'   CSample$compute_conns()
        #' }
        compute_conns = function() {
            if (is.null(private$tangent_handler$tangent_images)) {
                stop("tangent images must be specified.")
            }
            private$conns <- private$tangent_handler$compute_conns()
        },

        #' Compute Vectorized Tangent Images
        #'
        #' @description
        #' This function computes the vectorized tangent images from the tangent images. # nolint: line_length_linter
        #'
        #' @return None
        #' @throws Error if tangent images are not specified.
        #' @examples
        #' \dontrun{
        #'   CSample$compute_vecs()
        #' }
        compute_vecs = function() {
            if (is.null(private$tangent_handler$tangent_images)) {
                stop("tangent images must be specified.")
            }
            private$vec_imgs <- private$tangent_handler$compute_vecs()
        },

        #' Compute Unvectorized Tangent Images
        #'
        #' @description
        #' This function computes the tangent images from the vector images.
        #'
        #' @return None
        #' @throws Error if `vec_imgs` is not specified.
        #' @examples
        #' \dontrun{
        #'   CSample$compute_unvecs()
        #' }
        compute_unvecs = function() {
            if (is.null(private$vec_imgs)) stop("vec_imgs must be specified.")
            private$tangent_handler$set_tangent_images(
                private$vec_imgs[[1]],
                1:nrow(private$vec_imgs[[2]]) |> # nolint: seq_linter
                    purrr::map(\(i) {
                        private$metric$unvec(
                            private$vec_imgs[[1]],
                            private$vec_imgs[[2]][i, ]
                        )
                    })
            )
        },

        #' Compute Frechet Mean
        #'
        #' @description
        #' This function computes the Frechet mean of the sample.
        #'
        #' @param tol Tolerance for the convergence of the mean (default is 0.05). # nolint: line_length_linter
        #' @param max_iter Maximum number of iterations for the computation (default is 20). # nolint: line_length_linter
        #' @param lr Learning rate for the optimization algorithm (default is 0.2). # nolint: line_length_linter
        #'
        #' @return None
        #' @examples
        #' \dontrun{
        #'   CSample$compute_fmean(tol = 0.01, max_iter = 30, lr = 0.1)
        #' }
        compute_fmean = function(tol = 0.05, max_iter = 20, lr = 0.2) {
            compute_frechet_mean(self, tol, max_iter, lr)
        },

        #' Change Reference Point
        #'
        #' @description
        #' This function changes the reference point for the tangent images.
        #'
        #' @param new_ref_pt A new reference point, which must be a `dppMatrix` object. # nolint: line_length_linter
        #'
        #' @return None
        #' @throws Error if tangent images have not been computed or if `new_ref_pt` is not a `dppMatrix` object. # nolint: line_length_linter
        #' @examples
        #' \dontrun{
        #'   CSample$change_ref_pt(new_ref_pt)
        #' }
        change_ref_pt = function(new_ref_pt) {
            if (is.null(private$tangent_handler$tangent_images)) {
                stop("tangent images have not been computed")
            }
            if (!inherits(new_ref_pt, "dppMatrix")) {
                stop("new_ref_pt must be a dppMatrix object.")
            }
            private$tangent_handler$relocate_tangents(new_ref_pt)
        },

        #' Center the sample
        #'
        #' This function centers the sample by computing the Frechet mean if it is not already computed, and then changing the reference point to the computed Frechet mean. # nolint: line_length_linter
        #'
        #' @return None. This function is called for its side effects.
        #' @throws error if tangent images are not specified.
        #' @throws error if the sample is already centered.
        #' @examples
        #' \dontrun{
        #' obj$center()
        #' }
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

        #' Compute Sample Covariance
        #'
        #' This function computes the sample covariance matrix for the vector images. It first checks if the vector images are null, and if so, it computes the vectors, computing first the tangent images if necessary.  # nolint: line_length_linter
        #'
        #' @return None. This function is called for its side effects.
        #' @examples
        #' \dontrun{
        #'   obj <- CSample$new()
        #'   obj$compute_sample_cov()
        #' }
        compute_sample_cov = function() {
            if (self$vector_images |> is.null()) {
                if (private$tangent_handler$tangent_images |> is.null()) {
                    self$compute_tangents()
                }
                self$compute_vecs()
            }

            private$s_cov <- cov(private$vec_imgs[[2]])
        }
    ),
    active = list(
        #' @field connectomes Connectomes data
        connectomes = function() private$conns,

        #' @field tangent_images Tangent images data
        tangent_images = function() private$tangent_handler$tangent_images,

        #' @field vector_images Vector images data
        vector_images = function() private$vec_imgs,

        #' @field sample_size Sample size
        sample_size = function() private$n,

        #' @field matrix_size Matrix size
        matrix_size = function() private$p,

        #' @field mfd_dim Manifold dimension
        mfd_dim = function() private$d,

        #' @field is_centered Centering status
        is_centered = function() private$centered,

        #' @field frechet_mean Frechet mean
        frechet_mean = function() private$f_mean,

        #' @field metric Riemannian Metric used
        riem_metric = function() private$metric,

        #' @field variation Variation of the sample
        variation = function() private$var,

        #' @field sample_cov Sample covariance
        sample_cov = function() private$s_cov
    )
)
