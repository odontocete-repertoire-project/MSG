library(fs)

src_dir <- "Z:/MaiaProjects/Manzanilloanalysis/contours_update"  # source directory
dest_dir <- file.path(src_dir, "bigger_random_subset") # destination directory for new subset

# Create the destination directory if it does not exist
if (!dir.exists(dest_dir)) {
  dir.create(dest_dir)
}

# Define prefixes
prefixes <- c("msg")#, "tt", "sg")
num_files_to_select <- 200  # Number of files per prefix
set.seed(123)

for (prefix in prefixes) {
  files <- list.files(src_dir, pattern = paste0("^", prefix), full.names = TRUE)
  
  # Check if there are enough files
  num_files <- min(num_files_to_select, length(files))
  
  if (num_files > 0) {
    # Randomly select files
    selected_files <- sample(files, num_files)
    
    # Copy selected files to the new subdirectory
    file.copy(selected_files, dest_dir)
    
    cat("Copied", num_files, "files with prefix", prefix, "to", dest_dir, "\n")
  } else {
    cat("No files found for prefix", prefix, "\n")
  }
}


#Bocas -add prefix



source_dir <- "Z:/MaiaProjects/Bocas2023analysis/smm/smm_contour_subset"
dest_dir <- "Z:/MaiaProjects/Manzanilloanalysis/contours_w_outgroups/bigger_random_subset"
prefix <- "bocas."

if (!dir_exists(dest_dir)) {
  dir_create(dest_dir)
}

# List all files in source directory
all_files <- dir_ls(source_dir, type = "file")

# Select a random sample of 200 files (or fewer if not enough files)
set.seed(123)  # Set seed for reproducibility (optional)
sampled_files <- sample(all_files, min(200, length(all_files)))

# Copy files with new prefixed names
for (file in sampled_files) {
  file_name <- path_file(file)
  new_name <- paste0(prefix, file_name)
  new_path <- path(dest_dir, new_name)
  file_copy(file, new_path)
}
