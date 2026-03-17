# Function to build a dictionary of molecule information from CG topology
# This function processes the topology to identify:
# - Number of molecules
# - For each molecule: residue numbers, backbone and sidechain bead indices
# - Special handling for PHE (BB + 3 SC beads) and CYS (no BB/SC beads)
build_CG_molecule_dict <- function(topology) {
    # Filter out water rows first
    topology <- topology[topology$resid != "W",]
    
    # Initialize dictionary
    molecule_dict <- list()
    current_mol <- 1
    current_resno <- 0
    current_mol_res <- list()
    
    # Process topology line by line
    for (i in 1:nrow(topology)) {
        row <- topology[i,]
        
        # Check if we're starting a new molecule
        if (row$resno == 1 && current_resno == 6) {
            # Save previous molecule
            molecule_dict[[current_mol]] <- current_mol_res
            # Start new molecule
            current_mol <- current_mol + 1
            current_mol_res <- list()
        }
        
        # Update current residue number
        current_resno <- row$resno
        
        # Process residue based on type
        if (row$resid == "PHE") {
            if (row$elety == "BB") {
                current_mol_res[[paste0("PHE", row$resno, "_BB")]] <- row$eleno
            } else if (grepl("SC", row$elety)) {
                current_mol_res[[paste0("PHE", row$resno, "_", row$elety)]] <- row$eleno
            }
        } else if (row$resid == "CYS") {
            # For CYS, just store the bead number
            current_mol_res[[paste0("CYS", row$resno, "_", row$elety)]] <- row$eleno
        }
    }
    
    # Add the last molecule
    if (length(current_mol_res) > 0) {
        molecule_dict[[current_mol]] <- current_mol_res
    }
    
    return(molecule_dict)
} 