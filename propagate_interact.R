# Load necessary library
library(purrr)

# Step 1: Define your character vectors and their associated logical statuses
vector1 <- c("A", "F", "C")  # Vector 1
vector2 <- c("B", "D", "E")  # Vector 2
vector3 <- c("A", "E", "F")  # Vector 3

# Initial logical values for each vector
statuses <- c(TRUE, F, TRUE)  # TRUE/FALSE status for vector1, vector2, vector3

# Step 2: Create a named list of vectors
vectors <- list(vector1, vector2, vector3)

# Step 3: Function to propagate FALSE status across vectors
propagate_false <- function(vectors, statuses) {
    # Initialize a flag to check if further propagation is needed
    updated <- TRUE

    # Continue until no more changes
    while (updated) {
        updated <- FALSE

        # Identify all elements from vectors marked as FALSE
        false_elements <- unique(unlist(vectors[statuses == FALSE]))

        # Check each vector and update its status if it contains any false element
        for (i in seq_along(vectors)) {
            if (statuses[i] == TRUE && any(vectors[[i]] %in% false_elements)) {
                statuses[i] <- FALSE
                updated <- TRUE  # Set flag to TRUE to indicate changes made
            }
        }
    }

    return(statuses)
}

# Step 4: Apply the function to propagate FALSE statuses
final_statuses <- propagate_false(vectors, statuses)

# Combine vectors with their final statuses
result <- map2(vectors, final_statuses, ~ list(elements = .x, status = .y))

# Display the results
print(result)
