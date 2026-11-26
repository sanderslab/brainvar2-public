local_libs = "/home/ec2-user/R_libraries"
.libPaths(local_libs)

library(data.table)

# Set your working directory to the folder containing the files
setwd("/home/ec2-user/pseudobulk/means")

# Get a list of all files in the directory with a specific file extension
file_list <- list.files(pattern = "BrainVar2_period_\\d+_level_1_pseudobulk_mean.txt")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each file and rbind them to the combined_data
for (i in 1:length(file_list)) {
  # Read the data from the file
  current_data <- fread(file_list[i], data.table=F)

  period_number <- gsub(".*_period_(\\d+)_.*", "\\1", file_list[i])

  current_data$V1 <- paste0(current_data$V1, "_p", period_number)

  # Combine the data
  combined_data <- rbind(combined_data, current_data)
}

rownames(combined_data) <- combined_data$V1
combined_data <- combined_data[2:ncol(combined_data)]

tmp <- t(combined_data)

gene <- rownames(tmp)
final_df <- cbind(gene, tmp)

# Save the combined data to a new file
write.table(final_df, "Brainvar2_combined_pseudobulk_data_mean.txt", col.names=T, row.names = F, quote=F, sep="\t")
