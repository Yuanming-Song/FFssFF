#!/usr/bin/env Rscript

# Contact analysis for concentration-varied Base simulations
# New sims: Base_low (~14 mM), Base_high1 (~53 mM), Base_high2 (~87 mM)

step <- 1  # Process every frame

# Set the R library path (HPC location)
.libPaths("/dfs9/tw/yuanmis1/R_libs/")

# Define main directory
main_dir <- "/dfs9/tw/yuanmis1/mrsec/FFssFF"

# List all required packages
requiredpackages <- c(
    "sna",
    "ggplot2",
    "doParallel",
    "dplyr",
    "plotly",
    "bio3d",
    "geometry",
    "bigmemory",
    "SOMMD",
    "Rcpp"
)

# Load (and if needed, install) all required packages
for (packagename in requiredpackages) {
    if (!requireNamespace(packagename, quietly = TRUE)) {
        install.packages(packagename)
    }
    library(packagename, character.only = TRUE)
}

# Register parallel backend
registerDoParallel(cores = detectCores())

# Source the individual function scripts
source(file.path(main_dir, "Base_Rscript/base_fun.R"))
source(file.path(main_dir, "Base_Rscript/getContact_AA_by_type_base.R"))
source(file.path(main_dir, "Base_Rscript/build_AA_molecule_dict.R"))

# List of concentration-varied systems
systems <- list(
    Base_High1 = list(
        dir = "/dfs9/tw/yuanmis1/mrsec/FFssFF/AA/Mix/Base_high1/",
        pdb = "FFssFF_Base_mix_water_ions.pdb",
        prefix = "FFssFF_Base_AA_high1_contact",
        uniseltextlist = c("FSSF")
    ),
    Base_High2 = list(
        dir = "/dfs9/tw/yuanmis1/mrsec/FFssFF/AA/Mix/Base_high2/",
        pdb = "FFssFF_Base_mix_water_ions.pdb",
        prefix = "FFssFF_Base_AA_high2_contact",
        uniseltextlist = c("FSSF")
    ),
    Base_High3 = list(
        dir = "/dfs9/tw/yuanmis1/mrsec/FFssFF/AA/Mix/Base_high3/",
        pdb = "FFssFF_Base_mix_water_ions.pdb",
        prefix = "FFssFF_Base_AA_high3_contact",
        uniseltextlist = c("FSSF")
    )
)

# Create output directory
output_dir <- file.path(main_dir, "AA/Analysis/contact_analysis/output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Process each system
for (system_name in names(systems)) {
    cat(sprintf("\nAnalyzing system: %s\n", system_name))

    # Get configuration for this system
    config <- systems[[system_name]]

    # Check if PDB exists (simulation may not be done yet)
    pdb_path <- file.path(config$dir, config$pdb)
    if (!file.exists(pdb_path)) {
        cat(sprintf("PDB not found: %s — skipping\n", pdb_path))
        next
    }

    # Read PDB file
    mypdb <- read.pdb(pdb_path, hex = TRUE)

    # Create molecule dictionary
    molecule_dict <- build_AA_molecule_dict(mypdb$atom, "FSSF", is_acidic = FALSE)

    # Find all DCD files
    dcd_files <- list.files(config$dir, pattern = ".*\\.dcd$", full.names = TRUE)

    # Filter out equilibrium files
    dcd_files <- dcd_files[!grepl("eq\\.dcd$", dcd_files)]

    if (length(dcd_files) == 0) {
        cat(sprintf("No production DCD files found in %s — skipping\n", config$dir))
        next
    }

    # Sort DCD files by number (e.g., npt01.dcd, npt02.dcd, etc.)
    dcd_files <- dcd_files[order(as.numeric(gsub(".*?([0-9]+)\\.dcd$", "\\1", basename(dcd_files))))]

    # Initialize contact data dataframe
    contact_data <- data.frame(
        frame = integer(),
        pi_pi_contacts = integer(),
        bb_bb_contacts = integer(),
        total_contacts = integer(),
        stringsAsFactors = FALSE
    )

    processed_frame <- 0

    # Process each DCD file
    for (dcd_file in dcd_files) {
        cat(sprintf("\nProcessing trajectory file: %s\n", basename(dcd_file)))

        # Read trajectory
        mydcd <- read.dcd(dcd_file, verbose = FALSE, big = TRUE)
        mydcdcell <- read.dcd(dcd_file, verbose = FALSE, cell = TRUE)

        # Get number of frames
        tot_frame <- nrow(mydcdcell)
        cat(sprintf("Found %d frames\n", tot_frame))

        frame_indices <- seq(1, tot_frame, step)

        # Process frames in parallel
        contact_results <- foreach(frame = frame_indices) %dopar% {
            # Get coordinates for this frame
            frame_coords <- mydcd[frame,]

            # Analyze contacts using the molecule dictionary
            analyze_AA_contacts_by_type(frame_coords, molecule_dict, mydcdcell[frame,])
        }

        # Combine results into dataframe
        for (i in seq_along(contact_results)) {
            frame_num <- frame_indices[i] + processed_frame
            result <- contact_results[[i]]
            contact_data <- rbind(contact_data, data.frame(
                frame = frame_num,
                pi_pi_contacts = result$pi_pi_contacts,
                bb_bb_contacts = result$bb_bb_contacts,
                total_contacts = result$total_contacts,
                stringsAsFactors = FALSE
            ))
        }

        # Update processed frame count
        processed_frame <- processed_frame + tot_frame

        # Clean up memory
        rm(mydcd, mydcdcell, contact_results)
        gc()
    }

    # Save final results for this system
    output_file <- file.path(output_dir,
                           paste0(gsub("\\.pdb$", "", config$pdb),
                                 "_", tolower(system_name),
                                 "_step", step,
                                 ".rda"))
    save(contact_data, file = output_file)
    cat(sprintf("Saved contact data to: %s\n", output_file))
}
