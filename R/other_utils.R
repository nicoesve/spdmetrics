geometry <- function(log, exp, vec, unvec) {
    geom <- list(log = log, exp = exp, vec = vec, unvec = unvec)
    class(geom) <- "rgeom"
    return(geom)
}

rspdnorm <- function(n, refpt, disp, geom) {
    p <- refpt@Dim[1]
    d <- p * (p + 1) / 2
    mu <- rep(0, d)
    Sigma <- as.matrix(disp)
    smat <- MASS::mvrnorm(n, mu, Sigma)
    CSample$new(vec_imgs = list(refpt, smat), geom = geom, centered=FALSE)
}

relocate <- function(old_ref, new_ref, images, geom){
    images |> furrr::future_map(\(tan) geom$exp(old_ref, tan) |> 
        geom$log(Sigma=new_ref, Lambda=_))
}

