#' Compute the AIRM Logarithm 
#'
#' This function computes the Riemannian logarithmic map for the Affine-Invariant Riemannian Metric (AIRM).
#'
#' @param Sigma A symmetric positive-definite matrix of class `dppMatrix`, representing the reference point.
#' @param Lambda A symmetric positive-definite matrix of class `dppMatrix`, representing the target point.
#' 
#' @return A symmetric matrix of class `dspMatrix`, representing the tangent space image of `Lambda` at `Sigma`.
#' @examples
#' Sigma <- diag(2) |> Matrix::nearPD()$mat |> Matrix::pack()
#' Lambda <- diag(c(2, 3)) |> Matrix::nearPD()$mat |> Matrix::pack()
#' airm_log(Sigma, Lambda)
#' @export
airm_log <- function(Sigma, Lambda) {
    inheritance_flag <- list(Sigma, Lambda) |> 
	    purrr::map_lgl(\(x) inherits(x, "dppMatrix")) |> all()

    if(!inheritance_flag) stop("Both arguments should be of class dppMatrx")
    
    dim_flag <- list(Sigma, Lambda) |> purrr::map(\(x) x@Dim) |> 
        (\(l) identical(l[[1]],l[[2]]))()

    if(!dim_flag){stop("Arguments should be matrices of the same dimension")} 

    Sigma.sqrt <- expm::sqrtm(Sigma) |> Matrix::nearPD() |> _$mat
    Sigma.sqrt.inv <- Matrix::solve(Sigma.sqrt)
    Lambda |> (\(x) Sigma.sqrt.inv %*% x %*% Sigma.sqrt.inv)() |> 
    	Matrix::symmpart() |> as.matrix() |> 
	(\(x){
        tryCatch(expm::logm(x, method = "Eigen"), 
        error = function(e){print(e); expm::logm(x, method = "Higham08")})
    }
	)() |> (\(x) Sigma.sqrt %*% x %*% Sigma.sqrt)() |> 
    	Matrix::Matrix(sparse = FALSE, doDiag=FALSE) |> Matrix::symmpart() |> 
    	Matrix::pack() 
}

#' Compute the AIRM Exponential
#'
#' This function computes the Riemannian exponential map for the Affine-Invariant Riemannian Metric (AIRM).
#'
#' @param Sigma A symmetric positive-definite matrix of class `dppMatrix`, representing the reference point.
#' @param v A tangent vector of class `dspMatrix`, to be mapped back to the manifold at `Sigma`.
#'
#' @return A symmetric positive-definite matrix of class `dppMatrix`.
#' @examples
#' Sigma <- diag(2) |> Matrix::nearPD()$mat |> Matrix::pack()
#' v <- diag(c(1, 0.5)) |> Matrix::nearPD()$mat |> Matrix::pack()
#' airm_exp(Sigma, v)
#' @export
airm_exp <- function(Sigma, v) {
    inheritance_flag <- c(Sigma |> inherits("dppMatrix"), v |> inherits("dspMatrix")) |> all()
        
    if(!inheritance_flag) stop("Sigma should be of class dppMatrx and v should be of class dspMatrix")
    
    dim_flag <- list(Sigma, v) |> purrr::map(\(x) x@Dim) |> 
        (\(l) identical(l[[1]],l[[2]]))()

    if(!dim_flag){stop("Arguments should be matrices of the same dimension")} 

    Sigma.sqrt <- expm::sqrtm(Sigma) |> Matrix::nearPD() |> _$mat
    Sigma.sqrt.inv <- Matrix::solve(Sigma.sqrt)
    v |> (\(x) Sigma.sqrt.inv %*% x %*% Sigma.sqrt.inv)() |>
        Matrix::symmpart() |> as.matrix() |>
        expm::expm(method = "hybrid_Eigen_Ward") |>
        (\(x) Sigma.sqrt %*% x %*% Sigma.sqrt)() |> 
        Matrix::nearPD() |> _$mat |> Matrix::pack()
}

#' Vectorize at Identity Matrix
#'
#' Converts a symmetric matrix into a vector representation specific to operations at the identity matrix.
#'
#' @param v A symmetric matrix of class `dspMatrix`.
#'
#' @return A numeric vector, representing the vectorized tangent image.
#' @examples
#' v <- diag(c(1, sqrt(2))) |> Matrix::nearPD()$mat |> Matrix::pack()
#' vec_at_id(v)
#' @export
vec_at_id <- function(v) {
    if(!inherits(v, "dspMatrix")) stop("v should be an object of class dspMatrix")

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
#' @param Sigma A symmetric positive-definite matrix of class `dppMatrix`, representing the reference point.
#' @param v A symmetric matrix of class `dspMatrix`, representing a tangent vector.
#'
#' @return A numeric vector, representing the vectorized tangent image.
#' @examples
#' Sigma <- diag(2) |> Matrix::nearPD()$mat |> Matrix::pack()
#' v <- diag(c(1, 0.5)) |> Matrix::nearPD()$mat |> Matrix::pack()
#' airm_vec(Sigma, v)
#' @export
airm_vec <- function(Sigma, v) {
    inheritance_flag <- c(inherits(Sigma, "dppMatrix"), inherits(v, "dspMatrix")) |> all()
    if(!inheritance_flag){
	stop("Sigma should be of class dppMatrix and v should be of class dspMatrix")
    }

    dim_flag <- list(Sigma, v) |> purrr::map(\(x) x@Dim) |> 
        (\(x) identical(x[[1]], x[[2]]))() 
    if(!dim_flag) stop("Dimensions of Sigma and v don't match")
    
    Sigma.sqrt <- expm::sqrtm(Sigma) |> Matrix::nearPD() |> _$mat
    Sigma.sqrt.inv <- Matrix::solve(Sigma.sqrt)
    v |> (\(x) Sigma.sqrt.inv %*% x %*% Sigma.sqrt.inv)() |>
	    Matrix::Matrix(sparse=FALSE, doDiag=FALSE) |>
        Matrix::symmpart() |> Matrix::pack() |> vec_at_id()
}

#' Compute the Inverse Vectorization (AIRM)
#'
#' Converts a vector back into a tangent matrix relative to a reference point using AIRM.
#'
#' @param Sigma A symmetric positive-definite matrix of class `dppMatrix`, representing the reference point.
#' @param w A numeric vector, representing the vectorized tangent image.
#'
#' @return A symmetric matrix of class `dspMatrix`, representing the tangent vector.
#' @examples
#' Sigma <- diag(2) |> Matrix::nearPD()$mat |> Matrix::pack()
#' w <- c(1, sqrt(2), 2)
#' airm_unvec(Sigma, w)
#' @export
airm_unvec <- function(Sigma, w) {
    inheritance_flag <- c(inherits(Sigma, "dppMatrix"), inherits(w, what=c("numeric","vector"))) |> all()
    if(!inheritance_flag){
	stop("Sigma should be of class dppMatrix and v should be a numeric vector")
    }

    dim_flag <- list(Sigma@Dim[1]*(Sigma@Dim[1] +1)/2, length(w) |> as.numeric()) |>
	    do.call(identical, args = _)
    if(!dim_flag){stop("Dimensions of Sigma and v don't match")}
    
    Sigma.sqrt <- expm::sqrtm(Sigma) |> Matrix::nearPD() |> _$mat
    for (i in 1:Sigma@Dim[1]) {
        w[i * (i + 1) / 2] <- w[i * (i + 1) / 2] * sqrt(2)
    }
    w <- w / sqrt(2)
    methods::new("dspMatrix", x = w, Dim = as.integer(c(Sigma@Dim[1], Sigma@Dim[1]))) |>
        (\(x) Sigma.sqrt %*% x %*% Sigma.sqrt)() |>
        Matrix::Matrix(sparse = FALSE, doDiag=FALSE) |>
        Matrix::symmpart() |> Matrix::pack()
}
