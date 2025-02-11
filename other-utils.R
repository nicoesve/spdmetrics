#' Metric Object Constructor
#'
#' Constructs a metric object that contains the necessary functions for Riemannian operations.# nolint: line_length_linter.
#'
#' @param log A function representing the Riemannian logarithmic map. This function should accept a `dppMatrix` (the reference point) and another `dppMatrix` (the matrix whose logarithm is to be computed), and it outputs a `dspMatrix` (the tangent image).# nolint: line_length_linter.
#' @param exp A function representing the Riemannian exponential map. This function should accept a `dppMatrix` (the reference point) and a `dspMatrix` (the matrix whose exponential is to be computed) and return a `dppMatrix` (the image on the manifold).# nolint: line_length_linter.
#' @param vec A function representing the vectorization operation for tangent spaces. This function should accept a `dppMatrix` (the reference point) and a `dspMatrix` (the tangent image) and return a vector (the vectorized image).# nolint: line_length_linter.
#' @param unvec A function representing the inverse of the vectorization operation. This function should accept a `dppMatrix` (the reference point) and a vector (the vectorized image), and it returns a `dspMatrix` (the tangent image).# nolint: line_length_linter.
#'
#' @return An object of class `rmetric` containing the specified functions.
#' @export
metric <- function(log, exp, vec, unvec) {
    met <- list(log = log, exp = exp, vec = vec, unvec = unvec)
    class(met) <- "rmetric"
    return(met)
}

#' Generate Random Samples from a Riemannian Normal Distribution
#'
#' Simulates random samples from a Riemannian normal distribution on symmetric positive definite matrices.# nolint: line_length_linter.
#'
#' @param n Number of samples to generate.
#' @param refpt Reference point on the manifold, represented as a symmetric positive definite matrix. Must be an object of class `dppMatrix` from the Matrix package.# nolint: line_length_linter.
#' @param disp Dispersion matrix defining the spread of the distribution. Must be an object of class `dppMatrix` from the Matrix package.# nolint: line_length_linter.
#' @param met A metric object of class `rmetric`.
#'
#' @return An object of class `CSample` containing the generated samples.
#' @export
rspdnorm <- function(n, refpt, disp, met) {
    p <- refpt@Dim[1]
    d <- p * (p + 1) / 2
    mu <- rep(0, d)
    Sigma <- as.matrix(disp)
    smat <- MASS::mvrnorm(n, mu, Sigma)
    CSample$new(vec_imgs = list(refpt, smat), met = met, centered = FALSE)
}

#' Relocate Tangent Representations to a New Reference Point
#'
#' Changes the reference point for tangent space representations on a Riemannian manifold.# nolint: line_length_linter.
#'
#' @param old_ref A reference point on the manifold to be replaced. Must be an object of class `dppMatrix` from the Matrix package.# nolint: line_length_linter.
#' @param new_ref The new reference point on the manifold. Must be an object of class `dppMatrix` from the Matrix package.# nolint: line_length_linter.
#' @param images A list of tangent representations relative to the old reference point. Each element in the list must be an object of class `dspMatrix`.# nolint: line_length_linter.
#' @param met A metric object of class `rmetric`, containing functions for Riemannian operations (logarithmic map, exponential map, vectorization, and inverse vectorization).# nolint: line_length_linter.
#'
#' @return A list of tangent representations relative to the new reference point. Each element in the returned list will be an object of class `dspMatrix`.# nolint: line_length_linter.
#' @export
relocate <- function(old_ref, new_ref, images, met) {
    images |> furrr::future_map(
        \(tan) met$exp(old_ref, tan) |> met$log(Sigma = new_ref, Lambda = _)
    )
}

#' Compute the Frechet Mean
#'
#' This function computes the Frechet mean of a sample using an iterative algorithm.# nolint: line_length_linter.
#'
#' @param sample An object of class `CSample` containing the sample data.
#' @param tol A numeric value specifying the tolerance for convergence. Default is 0.05.# nolint: line_length_linter.
#' @param max_iter An integer specifying the maximum number of iterations. Default is 20.# nolint: line_length_linter.
#' @param lr A numeric value specifying the learning rate. Default is 0.2.
#'
#' @return The function updates the `frechet_mean` field of the `sample` object with the computed Frechet mean.# nolint: line_length_linter.
#'
#' @details
#' The function iteratively updates the reference point of the sample until the change in the reference point is less than the specified tolerance or the maximum number of iterations is reached. If the tangent images are not already computed, they will be computed before starting the iterations.# nolint: line_length_linter.
#'
#' @examples
#' \dontrun{
#' sample <- conn_sample$new(
#'     conns = list_of_dppMatrix_objects,
#'     metric = some_metric
#' )
#' compute_frechet_mean(sample, tol = 0.01, max_iter = 50, lr = 0.1)
#' }
compute_frechet_mean <- function(sample, tol = 0.05, max_iter = 20, lr = 0.2) {
    # Validating parameters
    if (!is.null(sample$frechet_mean)) {
        warning("The Frechet mean has already been computed.")
    }
    if (is.null(sample$tangent_images)) {
        message("tangent images were null, so they will be computed")
        sample$compute_tangents()
    }
    if (!is.numeric(tol)) stop("tol must be a numeric.")
    if (max_iter < 1) stop("max_iter must be at least 1.")

    aux_sample <- sample
    delta <- Inf
    iter <- 0

    while ((delta > tol) && (iter < max_iter)) {
        old_tan <- sample$tangent_images
        iter <- iter + 1
        old_ref_pt <- sample$tangent_handler$reference_point

        if (iter > max_iter) {
            warning("Computation of Frechet mean exceeded maximum
                number of iterations.")
        }

        # Computing the step
        tan_step <- lr * Reduce(`+`, old_tan) /
            aux_sample$sample_size
        new_ref_pt <- sample$metric$exp(old_ref_pt, tan_step)

        delta <- Matrix::norm(new_ref_pt - old_ref_pt, "F") /
            Matrix::norm(old_ref_pt, "F")

        # Mapping tangent images to the new step
        new_tan_imgs <- relocate(
            old_ref_pt, new_ref_pt, old_tan,
            sample$metric
        )

        aux_sample <- CSample$new(
            tan_imgs = list(new_ref_pt, new_tan_imgs),
            centered = FALSE, metric = sample$metric
        )
    }
    sample$frechet_mean <- sample$tangent_handler$reference_point
}

#' Validate Metric
#'
#' Validates that the metric is not NULL.
#'
#' @param metric The metric to validate.
#' @return None. Throws an error if the metric is NULL.
validate_metric <- function(metric) {
    if (is.null(metric)) stop("metric must be specified.")
}

#' Validate Connections
#'
#' Validates the connections input.
#'
#' @param conns List of connection matrices.
#' @param tan_imgs List of tangent images.
#' @param vec_imgs List of vector images.
#' @param centered Logical indicating if the data is centered.
#' @return None. Throws an error if the validation fails.
validate_conns <- function(conns, tan_imgs, vec_imgs, centered) {
    if (!is.null(conns)) {
        if (!is.null(tan_imgs) || !is.null(vec_imgs)) {
            stop("When initializing, if conns is not NULL, tan_imgs and vec_imgs must be NULL.") # nolint: line_length_linter.
        }
        if (!is.null(centered)) {
            warning("If conns is not NULL, centered is ignored")
        }
        class_flag <- conns |>
            purrr::map_lgl(\(x) inherits(x, "dppMatrix")) |>
            all()
        if (!class_flag) stop("conns must be a list of dppMatrix objects.")
    }
}

#' Validate Tangent Images
#'
#' Validates the tangent images input.
#'
#' @param tan_imgs List of tangent images.
#' @param vec_imgs List of vector images.
#' @param centered Logical indicating if the data is centered.
#' @return None. Throws an error if the validation fails.
validate_tan_imgs <- function(tan_imgs, vec_imgs, centered) {
    if (!is.null(tan_imgs)) {
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
            purrr::map_lgl(\(x) inherits(x, "dspMatrix")) |>
            all()
        if (!class_flag) {
            stop("The second element of tan_imgs must be a list of dspMatrix objects.") # nolint: line_length_linter.
        }
    }
}

#' Validate Vector Images
#'
#' Validates the vector images input.
#'
#' @param vec_imgs List of vector images.
#' @param centered Logical indicating if the data is centered.
#' @return None. Throws an error if the validation fails.
validate_vec_imgs <- function(vec_imgs, centered) {
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
}

#' Validate arguments for Riemannian logarithms
#'
#' @param sigma A dppMatrix object
#' @param lambda A dppMatrix object
#' @throws Error if sigma and lambda are not of the same dimensions
#' @return None
validate_log_args <- function(sigma, lambda) {
    inheritance_flag <- list(sigma, lambda) |>
        purrr::map(\(x) inherits(x, "dppMatrix")) |>
        all()

    if (!inheritance_flag) {
        stop("Both arguments should be of class dppMatrx")
    }

    dim_flag <- list(sigma, lambda) |>
        purrr::map(\(x) x@Dim) |>
        (\(l) identical(l[[1]], l[[2]]))()

    if (!dim_flag) {
        stop("Arguments should be matrices of the same dimension")
    }
}

#' Validate arguments for Riemannian logarithms
#'
#' @param sigma A dppMatrix object
#' @param lambda A dspMatrix object
#' @throws Error if sigma and lambda are not of the same dimensions
#' @return None
validate_exp_args <- function(sigma, v) {
    inheritance_flag <- c(
        sigma |> inherits("dppMatrix"),
        v |> inherits("dspMatrix")
    ) |>
        all()

    if (!inheritance_flag) {
        stop("sigma should be of class dppMatrx and v should be of class dspMatrix") # nolint: line_length_linter
    }

    dim_flag <- list(sigma, v) |>
        purrr::map(\(x) x@Dim) |>
        (\(l) identical(l[[1]], l[[2]]))()

    if (!dim_flag) {
        stop("Arguments should be matrices of the same dimension")
    }
}

#' Validate arguments for vectorization
#'
#' @param sigma A dppMatrix object
#' @param v A dspMatrix object
#' @throws Error if sigma and v are not of the same dimensions
#' @return None
validate_vec_args <- function(sigma, v) {
    validate_exp_args(sigma, v)
}

#' Validate arguments for inverse vectorization
#'
#' @param sigma A dppMatrix object
#' @param w A numeric vector
#' @throws Error if the dimensionalities don't match
#' @return None
validate_unvec_args <- function(sigma, w) {
    inheritance_flag <- c(
        inherits(sigma, "dppMatrix"),
        inherits(w, what = c("numeric", "vector"))
    ) |>
        all()
    if (!inheritance_flag) {
        stop("sigma should be of class dppMatrix and v should be a numeric vector") # nolint: line_length_linter
    }

    dim_flag <- list(
        sigma@Dim[1] * (sigma@Dim[1] + 1) / 2,
        length(w) |> as.numeric()
    ) |>
        do.call(identical, args = _)
    if (!dim_flag) {
        stop("Dimensions of sigma and v don't match")
    }
}

#' Create an Identity Matrix

#' @param sigma A matrix.
#' @return An identity matrix of the same dimensions as `sigma`.
id_matr <- function(sigma) {
    sigma |>
        nrow() |>
        diag() |>
        methods::as("dpoMatrix") |>
        Matrix::pack()
}

#' Differential of Matrix Logarithm Map
#'
#' Computes the differential of the matrix logarithm map at a point Sigma, evaluated at H # nolint: line_length_linter
#'
#' @param sigma A symmetric positive definite matrix of class dspMatrix
#' @param H A symmetric matrix representing tangent vector of class dsyMatrix
#' @return A symmetric matrix representing the differential evaluated at H of class dsyMatrix # nolint: line_length_linter
#' @export
dlog <- function(sigma, h) {
    if (!inherits(sigma, "dppMatrix")) {
        stop("sigma must be a symmetric positive definite matrix of class dppMatrix") # nolint: line_length_linter
    }
    if (!inherits(h, "dspMatrix")) {
        stop("H must be a symmetric matrix of class dspMatrix")
    }

    aux_matr <- sigma |> id_matr()
    n <- sigma |> nrow()
    t_vals <- seq(0, 1, length.out = 100) # Integration points
    dt <- t_vals[2] - t_vals[1]

    result <- Matrix::Matrix(0, n, n, sparse = FALSE)
    for (t in t_vals) {
        gamma_t <- t * sigma + (1 - t) * aux_matr
        gamma_t_inv <- Matrix::solve(gamma_t)
        result <- result + gamma_t_inv %*% h %*% gamma_t_inv * dt
    }

    result |>
        Matrix::symmpart() |>
        Matrix::pack()
}

#' Differential of Matrix Exponential Map
#'
#' Computes the differential of the matrix exponential map located at a point a, evaluated at x # nolint: line_length_linter
#'
#' @param a A symmetric matrix of class dspMatrix
#' @param x A symmetric matrix representing tangent vector of class dspMatrix
#' @return A positive definite symmetric matrix representing the differential located at a and evaluated at x, of class dppMatrix # nolint: line_length_linter
#' @export
dexp <- function(a, x) {
    if (!inherits(a, "dspMatrix")) {
        stop("a must be a symmetric matrix of class dspMatrix")
    }
    if (!inherits(x, "dspMatrix")) {
        stop("x must be a symmetric matrix of class dspMatrix")
    }

    n <- a |> nrow()
    t_vals <- seq(0, 1, length.out = 100) # Integration points
    dt <- t_vals[2] - t_vals[1]

    result <- Matrix::Matrix(0, n, n, sparse = FALSE)
    for (t in t_vals) {
        gamma_left <- (1 - t) * a |>
            as.matrix() |>
            expm::expm(method = "hybrid_Eigen_Ward")
        gamma_right <- t * a |>
            as.matrix() |>
            expm::expm(method = "hybrid_Eigen_Ward")
        result <- result + gamma_left %*% x %*% gamma_right * dt
    }

    result |>
        Matrix::nearPD() |>
        _$mat |>
        Matrix::pack()
}

#' Wrapper for the matrix logarithm
#'
#' @param x A matrix
#' @return Its matrix logarithm
safe_logm <- function(x) {
    tryCatch(
        expm::logm(x, method = "Eigen"),
        error = function(e) {
            print(e)
            expm::logm(x, method = "Higham08")
        }
    )
}

#' Default reference point
#'
#' @param p the dimension
#' @return A diagonal matrix of the desired dimension
default_ref_pt <- function(p) {
    diag(p) |>
        methods::as("dpoMatrix") |>
        Matrix::pack()
}
