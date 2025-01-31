euclidean_log <- function(ref_pt, mfd_pt) {
    if (ref_pt |> inherits("dppMatrix") |> (\(x) !x)()) {
        stop("reference point should be a positive definite matrix (dppMatrix object)") # nolint: line_length_linter
    }
    if (mfd_pt |> inherits("dspMatrix") |> (\(x) !x)()) {
        stop("tangent point should be a symmetric matrix (dspMatrix object)") # nolint: line_length_linter
    }

    ref_pt - mfd_pt
}

euclidean_exp <- function(ref_pt, tangent) {
    if (ref_pt |> inherits("dppMatrix") |> (\(x) !x)()) {
        stop("reference point should be a positive definite matrix (dppMatrix object)") # nolint: line_length_linter
    }
    if (tangent |> inherits("dppMatrix") |> (\(x) !x)()) {
        stop("tangent should be a symmetric matrix (dspMatrix object)") # nolint: line_length_linter
    }
    tryCatch(
        {
            chol(ref_pt + tangent)
            ref_pt + tangent
        },
        error = function(e) {
            stop("Exponential map is not defined for those arguments")
        }
    )
}
