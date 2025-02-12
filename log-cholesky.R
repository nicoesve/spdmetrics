#' Compute the Log-Cholesky Logarithm
#'
#' This function computes the Riemannian logarithmic map using the Log-Cholesky metric for symmetric positive-definite matrices. The Log-Cholesky metric operates by transforming matrices via their Cholesky decomposition.
#'
#' @param ref_pt A symmetric positive-definite matrix of class `dppMatrix`, representing the reference point.
#' @param mfd_pt A symmetric positive-definite matrix of class `dppMatrix`, representing the target point.
#'
#' @return A symmetric matrix of class `dspMatrix`, representing the tangent space image of `mfd_pt` at `ref_pt`.
#' @export
log_cholesky_log <- function(ref_pt, mfd_pt) {
    validate_log_args(ref_pt, mfd_pt)

    # Compute Cholesky decompositions - get lower triangular factors
    l_ref <- ref_pt |>
        chol() |>
        t()
    l_mfd <- mfd_pt |>
        chol() |>
        t()

    # Compute off-diagonal difference and diagonal terms
    lower_diff <- l_mfd - l_ref
    diag_ratio <- diag(l_mfd) / diag(l_ref)
    diag_terms <- diag(l_ref) * log(diag_ratio)

    # Set diagonal terms
    diag(lower_diff) <- diag_terms

    # Project to SPD tangent space and return
    result <- l_ref %*% t(lower_diff) + lower_diff %*% t(l_ref)
    result |>
        Matrix::symmpart() |>
        Matrix::pack()
}

#' Compute the Log-Cholesky Exponential
#'
#' This function computes the Riemannian exponential map using the Log-Cholesky metric for symmetric positive-definite matrices. The map operates by transforming the tangent vector via Cholesky decomposition of the reference point.
#'
#' @param ref_pt A symmetric positive-definite matrix of class `dppMatrix`, representing the reference point.
#' @param tangent A symmetric matrix of class `dspMatrix`, representing the tangent vector to be mapped.
#'
#' @return A symmetric positive-definite matrix of class `dppMatrix`, representing the point on the manifold.
#' @export
log_cholesky_exp <- function(ref_pt, tangent) {
    validate_exp_args(ref_pt, tangent)

    # Compute Cholesky decomposition - get lower triangular factor
    l_ref <- ref_pt |>
        chol() |>
        t()

    # Transform tangent vector to Cholesky space
    l_inv <- solve(l_ref)
    temp <- l_inv %*% tangent %*% t(l_inv)
    temp_under_half <- temp |> half_underscore()
    x_l <- l_ref %*% temp_under_half
    x_l <- Matrix::tril(x_l)

    # Compute off-diagonal difference and diagonal terms
    lower_sum <- x_l + l_ref
    diag_ratio <- diag(x_l) / diag(l_ref)
    diag_terms <- diag(l_ref) * exp(diag_ratio)

    # Set diagonal terms
    diag(lower_sum) <- diag_terms

    # Return SPD matrix
    aux <- (lower_sum %*% t(lower_sum))
    aux_2 <- aux |>
        Matrix::nearPD() |>
        _$mat
    aux_2 |> Matrix::pack()
}

#' Isometry from tangent space at P to tangent space at identity
#'
#' @param sigma A symmetric positive-definite matrix of class `dppMatrix`
#' @param v A symmetric matrix of class `dspMatrix`
#' @return A symmetric matrix of class `dspMatrix`
#' @export
spd_isometry_to_identity <- function(sigma, v) {
    validate_vec_args(sigma, v)

    # Cholesky decomposition
    l_sigma <- sigma |>
        chol() |>
        t()
    l_inv <- solve(l_sigma)

    # Transform to Cholesky space
    x <- l_inv %*% v %*% t(l_inv)
    x_sym <- Matrix::symmpart(x)

    # Split into parts and apply isometry
    x_lower <- x_sym * Matrix::tril(matrix(1, nrow(x), ncol(x)), -1)
    x_diag <- diag(x_sym)

    y_lower <- x_lower
    y_diag <- x_diag / diag(l_sigma)

    # Reconstruct symmetric matrix
    y <- Matrix::Matrix(0, nrow(x), ncol(x))
    y[lower.tri(y, diag = FALSE)] <- y_lower[lower.tri(y, diag = FALSE)]
    diag(y) <- y_diag

    y |>
        Matrix::symmpart() |>
        Matrix::pack()
}

#' Compute the Log-Cholesky Vectorization
#'
#' @param sigma A symmetric positive-definite matrix of class `dppMatrix`
#' @param v A symmetric matrix of class `dspMatrix`
#' @return A numeric vector representing the vectorized tangent image
#' @export
log_cholesky_vec <- function(sigma, v) {
    validate_vec_args(sigma, v)

    # Apply isometry then vectorize at identity
    v |>
        spd_isometry_to_identity(sigma = sigma, v = _) |>
        vec_at_id()
}

#' Reverse isometry from tangent space at identity to tangent space at P
#'
#' @param sigma A symmetric positive-definite matrix of class `dppMatrix`
#' @param v A symmetric matrix of class `dspMatrix`
#' @return A symmetric matrix of class `dspMatrix`
spd_isometry_from_identity <- function(sigma, v) {
    validate_vec_args(sigma, v)

    # Get Cholesky decomposition
    l_sigma <- sigma |>
        chol() |>
        t()

    # Split input matrix
    v_sym <- Matrix::symmpart(v)
    v_lower <- v_sym * Matrix::tril(matrix(1, nrow(v), ncol(v)), -1)
    v_diag <- diag(v_sym)

    # Reverse the isometry
    x_lower <- v_lower
    x_diag <- v_diag * diag(l_sigma)

    # Reconstruct matrix
    x <- Matrix::Matrix(0, nrow(v), ncol(v))
    x[lower.tri(x, diag = FALSE)] <- x_lower[lower.tri(x, diag = FALSE)]
    diag(x) <- x_diag

    # Transform back using Cholesky
    x |>
        Matrix::symmpart() |>
        (\(mat) l_sigma %*% mat %*% t(l_sigma))() |>
        Matrix::symmpart() |>
        Matrix::pack()
}

#' Compute the Log-Cholesky Inverse Vectorization
#'
#' @param sigma A symmetric positive-definite matrix of class `dppMatrix`
#' @param w A numeric vector representing the vectorized tangent image
#' @return A symmetric matrix of class `dspMatrix`
#' @export
log_cholesky_unvec <- function(sigma, w) {
    validate_unvec_args(sigma, w)

    # First undo vec_at_id like in airm_unvec
    w_scaled <- w
    for (i in 1:sigma@Dim[1]) {
        w_scaled[i * (i + 1) / 2] <- w_scaled[i * (i + 1) / 2] * sqrt(2)
    }
    w_scaled <- w_scaled / sqrt(2)

    # Create matrix and reverse isometry
    methods::new(
        "dsyMatrix",
        x = w_scaled,
        Dim = as.integer(c(sigma@Dim[1], sigma@Dim[1]))
    ) |>
        Matrix::pack() |>
        spd_isometry_from_identity(sigma = sigma, v = _)
}
