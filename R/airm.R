#' Compute the AIRM Logarithm
#'
#' This function computes the Riemannian logarithmic map for the Affine-Invariant Riemannian Metric (AIRM). # nolint: line_length_linter
#'
#' @param sigma A symmetric positive-definite matrix of class `dppMatrix`, representing the reference point. # nolint: line_length_linter
#' @param lambda A symmetric positive-definite matrix of class `dppMatrix`, representing the target point. # nolint: line_length_linter
#'
#' @return A symmetric matrix of class `dspMatrix`, representing the tangent space image of `lambda` at `sigma`. # nolint: line_length_linter
#' @examples
#' sigma <- diag(2) |> Matrix::nearPD()$mat |> Matrix::pack()
#' lambda <- diag(c(2, 3)) |> Matrix::nearPD()$mat |> Matrix::pack()
#' airm_log(sigma, lambda)
#' @export
airm_log <- function(sigma, lambda) {
    inheritance_flag <- list(sigma, lambda) |>
        purrr::map(\(x) inherits(x, "dppMatrix")) |>
        all()

    if (!inheritance_flag) {
        stop("Both arguments should be of class dppMatrx")
    }

    dim_flag <- list(sigma, lambda) |>
        purrr::map(\(x) x@Dim) |>
        (\(l) identical(l[[1]], l[[2]]))()

    if (!dim_flag){
        stop("Arguments should be matrices of the same dimension")
    }

    sigma_sqrt <- expm::sqrtm(sigma) |>
        Matrix::nearPD() |>
        _$mat

    sigma_sqrt_inv <- Matrix::solve(sigma_sqrt)

    lambda |> (\(x) sigma_sqrt_inv %*% x %*% sigma_sqrt_inv)() |>
        Matrix::symmpart() |>
        as.matrix() |>
        (
            \(x){
                tryCatch(
                    expm::logm(x, method = "Eigen"),
                    error = function(e) {
                        print(e)
                        expm::logm(x, method = "Higham08")
                    }
                )
            }
        )() |>
        (\(x) sigma_sqrt %*% x %*% sigma_sqrt)() |>
        Matrix::Matrix(sparse = FALSE, doDiag=FALSE) |>
        Matrix::symmpart() |>
        Matrix::pack()
}

#' Compute the AIRM Exponential
#'
#' This function computes the Riemannian exponential map for the Affine-Invariant Riemannian Metric (AIRM). # nolint: line_length_linter
#'
#' @param sigma A symmetric positive-definite matrix of class `dppMatrix`, representing the reference point. # nolint: line_length_linter
#' @param v A tangent vector of class `dspMatrix`, to be mapped back to the manifold at `sigma`. # nolint: line_length_linter
#'
#' @return A symmetric positive-definite matrix of class `dppMatrix`.
#' @examples
#' sigma <- diag(2) |> Matrix::nearPD()$mat |> Matrix::pack()
#' v <- diag(c(1, 0.5)) |> Matrix::nearPD()$mat |> Matrix::pack()
#' airm_exp(sigma, v)
#' @export
airm_exp <- function(sigma, v) {
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

    sigma_sqrt <- expm::sqrtm(sigma) |> Matrix::nearPD() |> _$mat
    sigma_sqrt_inv <- Matrix::solve(sigma_sqrt)
    v |> (\(x) sigma_sqrt_inv %*% x %*% sigma_sqrt_inv)() |>
        Matrix::symmpart() |> 
        as.matrix() |>
        expm::expm(method = "hybrid_Eigen_Ward") |>
        (\(x) sigma_sqrt %*% x %*% sigma_sqrt)() |> 
        Matrix::nearPD() |> 
        _$mat |> 
        Matrix::pack()
}

#' Vectorize at Identity Matrix
#'
#' Converts a symmetric matrix into a vector representation specific to operations at the identity matrix. # nolint: line_length_linter
#'
#' @param v A symmetric matrix of class `dspMatrix`.
#'
#' @return A numeric vector, representing the vectorized tangent image.
#' @examples
#' v <- diag(c(1, sqrt(2))) |> Matrix::nearPD()$mat |> Matrix::pack()
#' vec_at_id(v)
#' @export
vec_at_id <- function(v) {
    if (!inherits(v, "dspMatrix")) {
        stop("v should be an object of class dspMatrix")
    }

    w <- v@x; w <- sqrt(2) * w
    for (i in 1:v@Dim[1]) {
        w[i * (i + 1) / 2] <- w[i * (i + 1) / 2] / sqrt(2)
    }
    upper_part <- vector("numeric", length = v@Dim[1] * (v@Dim[1] + 1) / 2)
    for (i in 1:v@Dim[1]) {
        for (j in 1:i) {
            upper_part[j + i * (i - 1) / 2] <- w[j + i * (i - 1) / 2]
        }
    }
    return(upper_part)
}

#' Compute the AIRM Vectorization of Tangent Space
#'
#' Vectorizes a tangent matrix into a vector in Euclidean space using AIRM.
#'
#' @param sigma A symmetric positive-definite matrix of class `dppMatrix`, representing the reference point. # nolint: line_length_linter
#' @param v A symmetric matrix of class `dspMatrix`, representing a tangent vector. # nolint: line_length_linter
#'
#' @return A numeric vector, representing the vectorized tangent image.
#' @examples
#' sigma <- diag(2) |> Matrix::nearPD()$mat |> Matrix::pack()
#' v <- diag(c(1, 0.5)) |> Matrix::nearPD()$mat |> Matrix::pack()
#' airm_vec(sigma, v)
#' @export
airm_vec <- function(sigma, v) {
    inheritance_flag <- c(
        inherits(sigma, "dppMatrix"),
        inherits(v, "dspMatrix")
    ) |> all()
    if (!inheritance_flag) {
        stop("sigma should be of class dppMatrix and v should be of class dspMatrix") # nolint: line_length_linter
    }

    dim_flag <- list(sigma, v) |>
        purrr::map(\(x) x@Dim) |>
        (\(x) identical(x[[1]], x[[2]]))()
    if (!dim_flag) {
        stop("Dimensions of sigma and v don't match")
    }

    sigma_sqrt <- expm::sqrtm(sigma) |>
        Matrix::nearPD() |>
        _$mat
    sigma_sqrt_inv <- Matrix::solve(sigma_sqrt)
    v |> (\(x) sigma_sqrt_inv %*% x %*% sigma_sqrt_inv)() |>
        Matrix::Matrix(sparse = FALSE, doDiag = FALSE) |>
        Matrix::symmpart() |>
        Matrix::pack() |>
        vec_at_id()
}

#' Compute the Inverse Vectorization (AIRM)
#'
#' Converts a vector back into a tangent matrix relative to a reference point using AIRM. # nolint: line_length_linter
#'
#' @param sigma A symmetric positive-definite matrix of class `dppMatrix`, representing the reference point. # nolint: line_length_linter
#' @param w A numeric vector, representing the vectorized tangent image.
#'
#' @return A symmetric matrix of class `dspMatrix`, representing the tangent vector. # nolint: line_length_linter
#' @examples
#' sigma <- diag(2) |> Matrix::nearPD()$mat |> Matrix::pack()
#' w <- c(1, sqrt(2), 2)
#' airm_unvec(sigma, w)
#' @export
airm_unvec <- function(sigma, w) {
    inheritance_flag <- c(
        inherits(sigma, "dppMatrix"),
        inherits(w, what = c("numeric", "vector"))
    ) |>
        all()
    if (!inheritance_flag) {
        stop("sigma should be of class dppMatrix and v should be a numeric vector") # nolint: line_length_linter
    }

    dim_flag <- list(
        sigma@Dim[1] * ( sigma@Dim[1] + 1 ) / 2,
        length(w) |> as.numeric()
    ) |>
        do.call(identical, args = _)
    if (!dim_flag) {
        stop("Dimensions of sigma and v don't match")
    }

    sigma_sqrt <- expm::sqrtm(sigma) |>
        Matrix::nearPD() |>
        _$mat
    for (i in 1:sigma@Dim[1]) {
        w[i * (i + 1) / 2] <- w[i * (i + 1) / 2] * sqrt(2)
    }
    w <- w / sqrt(2)
    methods::new(
        "dspMatrix",
        x = w,
        Dim = as.integer(c(sigma@Dim[1], sigma@Dim[1]))
    ) |>
        (\(x) sigma_sqrt %*% x %*% sigma_sqrt)() |>
        Matrix::Matrix(sparse = FALSE, doDiag = FALSE) |>
        Matrix::symmpart() |>
        Matrix::pack()
}
