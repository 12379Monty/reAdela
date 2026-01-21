# Simple R code to convert Table S1 to FASTA format
# Just the essential lines

library(readxl)

# Read the Excel file
spike_in <- read_excel("TableS1_spike_in_sequences.xlsx", sheet = "Table S1")

# Extract name and sequence columns
names <- spike_in$name
sequences <- spike_in$seq

# Write FASTA file
output <- file("spike_in_controls.fa", "w")
for (i in 1:length(names)) {
    writeLines(paste0(">", names[i]), output)
    writeLines(sequences[i], output)
}
close(output)

cat("Created spike_in_controls.fa with", length(names), "sequences\n")
