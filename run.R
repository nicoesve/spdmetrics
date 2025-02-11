# Script to run all R files in a directory
library(here)

#' Run all R files in a directory
#' 
#' @param dir_path Path to the directory containing R files. Defaults to current directory.
#' @param pattern Pattern to match files. Defaults to "\\.R$" for R files.
#' @param recursive Whether to search subdirectories. Defaults to FALSE.
#' @param exclude Vector of filenames to exclude. Defaults to empty vector.
#' @return Vector of successfully executed file names
run_all_r_files <- function(dir_path = ".", 
                           pattern = "\\.R$", 
                           recursive = FALSE,
                           exclude = c()) {
    
    # Get list of R files
    r_files <- list.files(
        path = dir_path,
        pattern = pattern,
        full.names = TRUE,
        recursive = recursive
    )
    
    # Remove excluded files
    r_files <- r_files[!basename(r_files) %in% exclude]
    
    # Sort files alphabetically
    r_files <- sort(r_files)
    
    # Keep track of successfully run files
    successful_files <- character()
    
    # Run each file
    for (file in r_files) {
        tryCatch({
            message(sprintf("\nExecuting %s...", basename(file)))
            source(file, echo = TRUE)
            successful_files <- c(successful_files, basename(file))
            message(sprintf("Successfully executed %s", basename(file)))
        }, error = function(e) {
            message(sprintf("Error in %s: %s", basename(file), e$message))
        }, warning = function(w) {
            message(sprintf("Warning in %s: %s", basename(file), w$message))
        })
    }
    
    # Print summary
    message("\nExecution Summary:")
    message(sprintf("Total files: %d", length(r_files)))
    message(sprintf("Successfully executed: %d", length(successful_files)))
    message(sprintf("Failed: %d", length(r_files) - length(successful_files)))
    
    if (length(successful_files) > 0) {
        message("\nSuccessfully executed files:")
        message(paste("-", successful_files, collapse = "\n"))
    }
    
    if (length(successful_files) < length(r_files)) {
        failed_files <- setdiff(basename(r_files), successful_files)
        message("\nFailed files:")
        message(paste("-", failed_files, collapse = "\n"))
    }
    
    return(invisible(successful_files))
}

# Example usage:
# To run all R files in current directory:
# run_all_r_files()

# To run all R files in a specific directory:
# run_all_r_files("path/to/directory")

# To exclude certain files:
# run_all_r_files(exclude = c("setup.R", "cleanup.R"))

# To run files recursively including subdirectories:
# run_all_r_files(recursive = TRUE)
