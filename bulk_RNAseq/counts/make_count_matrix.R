local_libs = "/home/ec2-user/R_libraries"
.libPaths(local_libs)

library(data.table)

# Get list of input files
input_files <- list.files(pattern = "_counts.csv")

# Initialize an empty list to store data frames
dfs <- list()

# Read each file into a data frame and store in the list
for (file in input_files) {
  df <- fread(file, data.table=F, header=FALSE, col.names=c("gene", gsub("_counts.csv", "", file)))
  dfs[[length(dfs) + 1]] <- df
}

# Merge data frames based on 'gene' column
merged_df <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), dfs)

# Write the merged data frame to a tab-separated text file
fwrite(merged_df, file = "merged_count_matrix.txt", sep = "\t", row.names = FALSE, col.names=T, quote=F)