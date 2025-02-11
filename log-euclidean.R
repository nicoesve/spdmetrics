#' Compute the Log-Euclidean Logarithm
#'
#' @param ref_pt A reference point.
#' @param mfd_pt A point on the manifold.
#'
#' @return The tangent space image of `mfd_pt` at `ref_pt`.
#' @export
log_euclidean_log <- function(ref_pt, mfd_pt) {
    validate_log_args(ref_pt, mfd_pt)

    aux_matr_1 <- ref_pt |> safe_logm()
    aux_matr_2 <- mfd_pt |> safe_logm()

    aux_matr_3 <- aux_matr_1 - aux_matr_2

    dexp(ref_pt, aux_matr_3)
}

#' Compute the Log-Euclidean Exponential
#'
#' This function computes the Euclidean exponential map.
#'
#' @param ref_pt A reference point.
#' @param tangent A tangent vector to be mapped back to the manifold at `ref_pt`. # nolint: line_length_linter
#'
#' @return The point on the manifold corresponding to the tangent vector at `ref_pt`. # nolint: line_length_linter
#' @export
log_euclidean_exp <- function(ref_pt, tangent) {
    validate_exp_args(ref_pt, tangent)

    # compute the functions that compose the exponential
    aux_matr_1 <- ref_pt |>
        as.matrix() |>
        expm::logm()
    aux_matr_2 <- tangent |>
        (\(x) dlog(ref_pt, x))()

    aux_matr_3 <- aux_matr_1 + aux_matr_2

    aux_matr_3 |>
        expm::expm(method = "hybrid_Eigen_Ward") |>
        Matrix::nearPD() |>
        _$mat |>
        Matrix::pack()
}

#' Vectorize at Identity Matrix (Euclidean)
#'
#' Converts a symmetric matrix into a vector representation. # nolint: line_length_linter
#'
#' @param sigma A symmetric matrix.
#' @param v A vector.
#'
#' @return A numeric vector, representing the vectorized tangent image.
#' @export
log_euclidean_vec <- function(sigma, v) {
    airm_vec(sigma |> id_matr(), v)
}

#' Compute the Inverse Vectorization (Euclidean)
#'
#' Converts a vector back into a tangent matrix relative to a reference point using Euclidean metric. # nolint: line_length_linter
#'
#' @param sigma A symmetric matrix.
#' @param w A numeric vector, representing the vectorized tangent image.
#'
#' @return A symmetric matrix, representing the tangent vector.
#' @export
log_euclidean_unvec <- function(sigma, w) {
    airm_unvec(sigma |> id_matr(), w)
}
