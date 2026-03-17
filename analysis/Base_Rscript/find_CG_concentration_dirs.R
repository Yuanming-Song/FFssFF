# Function to find concentration directories
find_CG_concentration_dirs <- function(base_dir) {
    # Look for directories matching concentration patterns (e.g., 100mM, 40mM)
    dirs <- list.dirs(base_dir, recursive = FALSE)
    conc_dirs <- dirs[grep("[0-9]+mM", basename(dirs))]
    return(conc_dirs)
} 