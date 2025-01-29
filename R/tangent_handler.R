TangentImageHandler <- R6::R6Class(
    classname = "TangentImageHandler",
    private = list(
        reference_point = NULL,  # The current reference point on the manifold
        tangent_images = NULL,   # List of tangent images (dspMatrix objects)
        metric = NULL            # Metric object for operations
    ),
    public = list(
        initialize = function(metric, reference_point = NULL) {
            if (is.null(metric)) stop("Metric object must be provided.")
            private$metric <- metric
            private$reference_point <- reference_point
            private$tangent_images <- list()
        },
        
        set_reference_point = function(new_ref_pt) {
            if (!inherits(new_ref_pt, "dppMatrix")) {
                stop("Reference point must be of class dppMatrix.")
            }
            if (!is.null(private$reference_point) && !is.null(private$tangent_images)) {
                private$tangent_images <- private$tangent_images |> 
                    purrr::map(\(tan) private$metric$exp(private$reference_point, tan)) |>
                    purrr::map(\(point) private$metric$log(new_ref_pt, point))
            }
            private$reference_point <- new_ref_pt
        },
        
        compute_tangents = function(manifold_points) {
            if (is.null(private$reference_point)) {
                stop("Reference point must be set before computing tangents.")
            }
            private$tangent_images <- manifold_points |> 
                purrr::map(\(point) private$metric$log(private$reference_point, point))
        },
        
        set_tangent_images = function(reference_point, tangent_images) {
            if (!inherits(reference_point, "dppMatrix")) {
                stop("Reference point must be of class dppMatrix.")
            }
            class_flag <- tangent_images |> purrr::map_lgl(\(x) inherits(x, "dspMatrix")) |> all()
            if (!class_flag) {
                stop("All tangent images must be of class dspMatrix.")
            }
            private$reference_point <- reference_point
            private$tangent_images <- tangent_images
        },
        
        add_tangent_image = function(image) {
            if (!inherits(image, "dspMatrix")) {
                stop("Tangent image must be of class dspMatrix.")
            }
            private$tangent_images <- c(private$tangent_images, list(image))
        },
        
        get_tangent_images = function() {
            private$tangent_images
        },
        
        relocate_tangents = function(new_ref_pt) {
            self$set_reference_point(new_ref_pt)
        },
        
        clear_tangents = function() {
            private$tangent_images <- list()
        }
    ),
    active = list(
        reference_point = function(value) {
            if (missing(value)) return(private$reference_point)
            self$set_reference_point(value)
        },
        tangent_images = function() {
            self$get_tangent_images()
        }
    )
)
