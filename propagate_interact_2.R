# Load necessary libraries
library(purrr)

# Step 1: Define your character vectors and associated logical statuses
vector1 <- c("X", "Z", "C")  # Vector 1
vector2 <- c("B", "D", "E")  # Vector 2
vector3 <- c("A", "E", "F")  # Vector 3

# Initial logical values for each vector
statuses <- c(TRUE, F, TRUE)  # TRUE/FALSE status for vector1, vector2, vector3

# Step 2: Create a list of vectors
vectors <- list(vector1, vector2, vector3)

# Step 3: Function to propagate FALSE status across vectors
propagate_false_status <- function(vectors, statuses) {
    # Repeat until no more changes
    repeat {
        # Collect all elements from vectors marked as FALSE
        false_elements <- unique(unlist(vectors[statuses == FALSE]))

        # Check each vector and update status if it contains any false element
        new_statuses <- map_lgl(seq_along(vectors), function(i) {
            if (statuses[i] == TRUE && any(vectors[[i]] %in% false_elements)) {
                return(FALSE)
            } else {
                return(statuses[i])
            }
        })

        # Stop if no changes in statuses
        if (all(new_statuses == statuses)) break

        # Update statuses
        statuses <- new_statuses
    }

    # Create a named vector of unique elements with their final statuses
    all_elements <- unique(unlist(vectors))
    element_status <- setNames(rep(TRUE, length(all_elements)), all_elements)

    # Set the final status of each element based on the propagated statuses
    for (i in seq_along(vectors)) {
        element_status[vectors[[i]]] <- statuses[i]
    }

    return(element_status)
}

# Step 4: Apply the function to get the propagated status of each element
final_element_status <- propagate_false_status(vectors, statuses)

# Display the final status of each unique element
print(final_element_status)
