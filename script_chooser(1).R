#!/usr/bin/env Rscript

# Retrieve the condition input
condition <- commandArgs(trailingOnly = TRUE)[1]

# Condition

if (condition == "rna_seq") {
  t <-
} else if (condition == "microarrays") {
  result <- FALSE
} else {
  stop("Invalid condition provided")
}

# Print
write.table(result, "chooser.txt")
