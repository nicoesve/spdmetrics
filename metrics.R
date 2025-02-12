source("airm.R", "log-euclidean.R", "euclidean.R", "log-cholesky.R", "bures-wasserstein.R")

# Create the metrics objects
airm <- metric(
    log = airm_log,
    exp = airm_exp,
    vec = airm_vec,
    unvec = airm_unvec
)

log_euclidean <- metric(
    log = log_euclidean_log,
    exp = log_euclidean_exp,
    vec = log_euclidean_vec,
    unvec = log_euclidean_unvec
)

euclidean <- metric(
    log = euclidean_log,
    exp = euclidean_exp,
    vec = euclidean_vec,
    unvec = euclidean_unvec
)

log_cholesky <- metric(
    log = log_cholesky_log,
    exp = log_cholesky_exp,
    vec = log_cholesky_vec,
    unvec = log_cholesky_unvec
)

bures_wasserstein <- metric(
    log = bures_wasserstein_log,
    exp = bures_wasserstein_exp,
    vec = bures_wasserstein_vec,
    unvec = bures_wasserstein_unvec
)

# Save the metrics as package data
usethis::use_data(airm, log_euclidean, euclidean, log_cholesky, 
                  bures_wasserstein, internal = FALSE, overwrite = TRUE)