#' Metric Object Constructor
#'
#' Constructs a metric object that contains the necessary functions for Riemannian operations.
#'
#' @param log A function representing the Riemannian logarithmic map. This function should accept a `dppMatrix` (the reference point) and another `dppMatrix` (the matrix whose logarithm is to be computed), and it outputs a `dspMatrix` (the tangent image).
#' @param exp A function representing the Riemannian exponential map. This function should accept a `dppMatrix` (the reference point) and a `dspMatrix` (the matrix whose exponential is to be computed) and return a `dppMatrix` (the image on the manifold).
#' @param vec A function representing the vectorization operation for tangent spaces. This function should accept a `dppMatrix` (the reference point) and a `dspMatrix` (the tangent image) and return a vector (the vectorized image).
#' @param unvec A function representing the inverse of the vectorization operation. This function should accept a `dppMatrix` (the reference point) and a vector (the vectorized image), and it returns a `dspMatrix` (the tangent image).
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
#' Simulates random samples from a Riemannian normal distribution on symmetric positive definite matrices.
#'
#' @param n Number of samples to generate.
#' @param refpt Reference point on the manifold, represented as a symmetric positive definite matrix. Must be an object of class `dppMatrix` from the Matrix package.
#' @param disp Dispersion matrix defining the spread of the distribution. Must be an object of class `dppMatrix` from the Matrix package.
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
#' Changes the reference point for tangent space representations on a Riemannian manifold.
#'
#' @param old_ref A reference point on the manifold to be replaced. Must be an object of class `dppMatrix` from the Matrix package.
#' @param new_ref The new reference point on the manifold. Must be an object of class `dppMatrix` from the Matrix package.
#' @param images A list of tangent representations relative to the old reference point. Each element in the list must be an object of class `dspMatrix`.
#' @param met A metric object of class `rmetric`, containing functions for Riemannian operations (logarithmic map, exponential map, vectorization, and inverse vectorization).
#'
#' @return A list of tangent representations relative to the new reference point. Each element in the returned list will be an object of class `dspMatrix`.
#' @export
relocate <- function(old_ref, new_ref, images, met) {
    images |> furrr::future_map(
        \(tan) met$exp(old_ref, tan) |> met$log(Sigma = new_ref, Lambda = _)
    )
}

#' Compute the Frechet Mean
#'
#' This function computes the Frechet mean of a sample using an iterative algorithm.
#'
#' @param sample An object of class `CSample` containing the sample data.
#' @param tol A numeric value specifying the tolerance for convergence. Default is 0.05.
#' @param max_iter An integer specifying the maximum number of iterations. Default is 20.
#' @param lr A numeric value specifying the learning rate. Default is 0.2.
#'
#' @return The function updates the `frechet_mean` field of the `sample` object with the computed Frechet mean.
#'
#' @details
#' The function iteratively updates the reference point of the sample until the change in the reference point is less than the specified tolerance or the maximum number of iterations is reached. If the tangent images are not already computed, they will be computed before starting the iterations.
#'
#' @examples
#' \dontrun{
#' sample <- CSample$new(conns = list_of_dppMatrix_objects, metric = some_metric)
#' compute_fmean(sample, tol = 0.01, max_iter = 50, lr = 0.1)
#' }
compute_fmean <- function(sample, tol = 0.05, max_iter = 20, lr = 0.2) {
    if (!is.null(sample$frechet_mean)) { 
        warning("The Frechet mean has already been computed.")
    }
    if (is.null(sample$tangent_images)) { 
        message("tangent images were null, so they will be computed")
        sample$compute_tangents()
    }
    if (!is.numeric(tol)) stop("tol must be a numeric.")
    if (max_iter < 1) stop("max_iter must be at least 1.")

    aux_sample <- sample; delta <- Inf; iter <- 0

    while ((delta > tol) && (iter < max_iter)) {
        old_tan <- sample$tangent_images; iter <- iter + 1 
        old_ref_pt <- sample$tangent_handler$reference_point

        if (iter > max_iter) {
            warning("Computation of Frechet mean exceeded maximum 
                number of iterations.")
        }

        tan_step <- lr * Reduce(`+`, old_tan) / 
            aux_sample$sample_size
        new_ref_pt <- sample$metric$exp(old_ref_pt, tan_step)

        delta <- Matrix::norm(new_ref_pt - old_ref_pt, "F") / 
                Matrix::norm(old_ref_pt, "F")

        new_tan_imgs <- relocate(old_ref_pt, new_ref_pt, old_tan, 
            sample$metric) 

        aux_sample <- CSample$new(
            tan_imgs = list(new_ref_pt, new_tan_imgs),
            centered = FALSE, metric = sample$metric)
    }
    sample$frechet_mean <- sample$tangent_handler$reference_point
}



