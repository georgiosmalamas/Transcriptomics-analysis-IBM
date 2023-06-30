
#For microarrays platform
platform <- "microarrays"

if (platform == "microarrays") {
# Run the R script for microarrays
source("test_microarrays_analysis.R")
} else if (platform == "RNA-seq") {
# Run the R script for RNA-seq
source("RNA-Seq_analysis.R")
} else {
# Print an error message if the condition is not valid
stop("Invalid condition. Please choose either microarrays or RNA-seq.")
}
