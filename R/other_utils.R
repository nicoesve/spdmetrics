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



