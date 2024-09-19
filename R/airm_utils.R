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
