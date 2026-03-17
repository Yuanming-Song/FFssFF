# Function to find trajectory files
find_CG_trajectory_files <- function(conc_dir) {
    # Find all gro files that are part of the production run
    gro_files <- list.files(conc_dir, pattern = ".*md.*\\.gro$", full.names = TRUE)
    # Find corresponding xtc files
    xtc_files <- list.files(conc_dir, pattern = ".*md.*\\.xtc$", full.names = TRUE)
    
    # Sort files to ensure they're in the correct order
    gro_files <- sort(gro_files)
    xtc_files <- sort(xtc_files)
    
    return(list(gro = gro_files, xtc = xtc_files))
} 