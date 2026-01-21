#!/usr/bin/env Rscript

# R script to convert Table S1 (spike-in sequences) to FASTA format
# Wilson et al. 2022 - Cell Reports Methods

# Install required package if needed
if (!require("readxl", quietly = TRUE)) {
    install.packages("readxl", repos = "https://cloud.r-project.org/")
}

library(readxl)

# Function to write FASTA format
write_fasta <- function(names, sequences, output_file) {
    # Open file connection
    con <- file(output_file, "w")
    
    # Write each sequence in FASTA format
    for (i in seq_along(names)) {
        # Write header line (starts with >)
        writeLines(paste0(">", names[i]), con)
        # Write sequence line
        writeLines(sequences[i], con)
    }
    
    # Close connection
    close(con)
    
    cat("FASTA file written to:", output_file, "\n")
    cat("Number of sequences:", length(names), "\n")
}

# Main script
main <- function() {
    # Input and output file paths
    input_file <- "TableS1_spike_in_sequences.xlsx"
    output_file <- "spike_in_controls.fa"
    
    # Read the Excel file - Table S1 sheet contains the actual data
    cat("Reading Excel file:", input_file, "\n")
    spike_in_data <- read_excel(input_file, sheet = "Table S1")
    
    # Display information about the data
    cat("\nData summary:\n")
    cat("  Total sequences:", nrow(spike_in_data), "\n")
    cat("  Columns:", paste(colnames(spike_in_data), collapse = ", "), "\n")
    
    # Show length distribution
    cat("\nSequence length distribution:\n")
    print(table(spike_in_data$len))
    
    # Show methylation status
    cat("\nMethylation status:\n")
    print(table(spike_in_data$methylated))
    
    # Extract name and sequence columns (first two columns)
    seq_names <- spike_in_data$name
    seq_sequences <- spike_in_data$seq
    
    # Verify no missing data
    if (any(is.na(seq_names)) || any(is.na(seq_sequences))) {
        stop("ERROR: Missing data in name or sequence columns!")
    }
    
    # Write to FASTA file
    cat("\nWriting FASTA file...\n")
    write_fasta(seq_names, seq_sequences, output_file)
    
    # Verify output
    cat("\nVerification:\n")
    cat("  Output file:", output_file, "\n")
    cat("  File size:", file.size(output_file), "bytes\n")
    
    # Show first few sequences as preview
    cat("\nFirst 3 sequences in FASTA format:\n")
    cat(paste(readLines(output_file, n = 6), collapse = "\n"), "\n")
    
    cat("\nDone!\n")
}

# Run the main function
main()
