# Generic C++ function for efficient pairwise contact analysis between molecules
# This function organizes coordinates by molecules and analyzes contacts between different molecules
# Input:
#   - coords1: First set of coordinates (Nx3 matrix)
#   - coords2: Second set of coordinates (Mx3 matrix)
#   - cutoff: Distance cutoff for contact definition
#   - tempcell: Box dimensions for PBC [x, y, z]
#   - mol_dict: Dictionary mapping molecule indices to coordinate indices
# Output:
#   - Integer count of contacts between molecules
cppFunction('
int contact_analysis_counter_AA(
    NumericMatrix coords1,
    NumericMatrix coords2,
    double cutoff,
    NumericVector tempcell,
    List mol_dict
) {
    int counter = 0;
    int n_mols = mol_dict.length();
    
    // Loop through each molecule pair (avoiding self-comparison)
    for (int i = 0; i < n_mols; ++i) {
        // Get indices of coordinates for first molecule
        IntegerVector mol1_indices = mol_dict[i];
        
        // Extract coordinates for first molecule into a new matrix
        NumericMatrix mol1_coords(mol1_indices.length(), 3);
        for (int k = 0; k < mol1_indices.length(); ++k) {
            mol1_coords(k, _) = coords1(mol1_indices[k] - 1, _);
        }
        
        // Compare with all other molecules (j > i to avoid duplicates)
        for (int j = i + 1; j < n_mols; ++j) {
            // Get indices of coordinates for second molecule
            IntegerVector mol2_indices = mol_dict[j];
            
            // Extract coordinates for second molecule into a new matrix
            NumericMatrix mol2_coords(mol2_indices.length(), 3);
            for (int k = 0; k < mol2_indices.length(); ++k) {
                mol2_coords(k, _) = coords2(mol2_indices[k] - 1, _);
            }
            
            // Apply periodic boundary conditions (PBC) correction
            // For each atom in first molecule
            for (int k = 0; k < mol1_coords.nrow(); ++k) {
                // Extract coordinates of reference atom from first molecule for PBC correction
                NumericVector ref_atom = mol1_coords(k, _);
                
                // Correct coordinates in x, y, z dimensions
                for (int dim = 0; dim < 3; ++dim) {
                    // Get box size for current dimension (x, y, or z) from tempcell vector
                    double box_size = tempcell[dim];
                    
                    // For each atom in second molecule
                    for (int l = 0; l < mol2_coords.nrow(); ++l) {
                        // Calculate distance in current dimension
                        double diff = mol2_coords(l, dim) - ref_atom[dim];
                        double abs_diff = std::abs(diff);
                        
                        // If distance is greater than half box size, apply PBC correction
                        if (abs_diff > box_size / 2) {
                            mol2_coords(l, dim) -= (diff / abs_diff) * box_size * 
                                std::ceil((abs_diff - box_size / 2) / box_size);
                        }
                    }
                }
            }
            
            // Calculate distances and count contacts between molecules
            for (int k = 0; k < mol1_coords.nrow(); ++k) {
                for (int l = 0; l < mol2_coords.nrow(); ++l) {
                    // Calculate Euclidean distance
                    double dist = 0;
                    for (int dim = 0; dim < 3; ++dim) {
                        dist += pow(mol1_coords(k, dim) - mol2_coords(l, dim), 2);
                    }
                    dist = sqrt(dist);
                    
                    // Count contact if within cutoff distance
                    if (dist <= cutoff) {
                        counter += 1;
                    }
                }
            }
        }
    }
    
    return counter;
}
')

# C++ function specifically for N-O contact analysis
cppFunction('
int contact_analysis_NO_AA(
    NumericMatrix n_coords,
    NumericMatrix o_coords,
    double cutoff,
    NumericVector tempcell,
    List n_mol_dict,
    List o_mol_dict
) {
    int counter = 0;
    int n_mols = n_mol_dict.length();
    
    // Loop through each molecule
    for (int i = 0; i < n_mols; ++i) {
        // Get indices of N coordinates for this molecule
        IntegerVector n_indices = n_mol_dict[i];
        
        // Extract N coordinates for this molecule
        NumericMatrix n_mol_coords(n_indices.length(), 3);
        for (int k = 0; k < n_indices.length(); ++k) {
            n_mol_coords(k, _) = n_coords(n_indices[k] - 1, _);
        }
        
        // Compare with O atoms from other molecules
        for (int j = 0; j < n_mols; ++j) {
            if (i == j) continue;  // Skip same molecule
            
            // Get indices of O coordinates for other molecule
            IntegerVector o_indices = o_mol_dict[j];
            
            // Extract O coordinates for other molecule
            NumericMatrix o_mol_coords(o_indices.length(), 3);
            for (int k = 0; k < o_indices.length(); ++k) {
                o_mol_coords(k, _) = o_coords(o_indices[k] - 1, _);
            }
            
            // Apply periodic boundary conditions (PBC) correction
            for (int k = 0; k < n_mol_coords.nrow(); ++k) {
                NumericVector ref_atom = n_mol_coords(k, _);
                
                for (int dim = 0; dim < 3; ++dim) {
                    double box_size = tempcell[dim];
                    
                    for (int l = 0; l < o_mol_coords.nrow(); ++l) {
                        double diff = o_mol_coords(l, dim) - ref_atom[dim];
                        double abs_diff = std::abs(diff);
                        
                        if (abs_diff > box_size / 2) {
                            o_mol_coords(l, dim) -= (diff / abs_diff) * box_size * 
                                std::ceil((abs_diff - box_size / 2) / box_size);
                        }
                    }
                }
            }
            
            // Calculate distances and count contacts
            for (int k = 0; k < n_mol_coords.nrow(); ++k) {
                for (int l = 0; l < o_mol_coords.nrow(); ++l) {
                    double dist = 0;
                    for (int dim = 0; dim < 3; ++dim) {
                        dist += pow(n_mol_coords(k, dim) - o_mol_coords(l, dim), 2);
                    }
                    dist = sqrt(dist);
                    
                    if (dist <= cutoff) {
                        counter += 1;
                    }
                }
            }
        }
    }
    
    return counter;
}
')

# Function to analyze backbone-backbone interactions
# This function handles contacts between backbone atoms (N-O pairs) in the AA model
# Input:
#   - frame_coords: Full frame coordinates
#   - tempcell: Box dimensions for PBC
#   - molecule_dict: Dictionary of molecule information
#   - cutoff: Distance cutoff for backbone interactions (default 3.4 Å)
# Output:
#   - Integer count of backbone-backbone contacts
analyze_AA_bb_bb_contacts_by_type <- function(frame_coords, tempcell, molecule_dict, cutoff = 3.4) {
    # Initialize lists to store N and O coordinates separately
    n_coords <- list()
    o_coords <- list()
    n_counter <- 1
    o_counter <- 1
    
    # Process each molecule
    for (mol_idx in seq_along(molecule_dict)) {
        mol <- molecule_dict[[mol_idx]]
        
        # Get N and O coordinates
        n_indices <- as.numeric(mol[["N"]])
        o_indices <- as.numeric(mol[["O"]])
        
        # Store N coordinates
        if (length(n_indices) > 0) {
            n_coords[[n_counter]] <- frame_coords[n_indices,]
            n_counter <- n_counter + 1
        }
        
        # Store O coordinates
        if (length(o_indices) > 0) {
            o_coords[[o_counter]] <- frame_coords[o_indices,]
            o_counter <- o_counter + 1
        }
    }
    
    # Convert lists to matrices
    n_coords_matrix <- do.call(rbind, n_coords)
    o_coords_matrix <- do.call(rbind, o_coords)
    
    # Create molecule dictionaries for N and O atoms
    n_mols <- length(molecule_dict)
    n_mol_dict <- lapply(1:n_mols, function(i) {
        start_idx <- (i-1)*6 + 1  # 6 N atoms per molecule
        end_idx <- i*6
        start_idx:end_idx
    })
    
    o_mol_dict <- lapply(1:n_mols, function(i) {
        start_idx <- (i-1)*4 + 1  # 4 O atoms per molecule
        end_idx <- i*4
        start_idx:end_idx
    })
    
    # Use C++ function for efficient pairwise contact analysis between N and O
    contacts <- contact_analysis_NO_AA(n_coords_matrix, o_coords_matrix, cutoff, tempcell, n_mol_dict, o_mol_dict)
    
    return(contacts)
}

# Function to analyze pi-pi interactions between PHE rings
# This function handles contacts between aromatic rings in the AA model
# Input:
#   - frame_coords: Full frame coordinates (flat format: x1,y1,z1,x2,y2,z2,...)
#   - tempcell: Box dimensions for PBC
#   - molecule_dict: Dictionary of molecule information
#   - cutoff: Distance cutoff for pi-pi interactions (default 5.4 Å)
# Output:
#   - Integer count of pi-pi contacts
analyze_AA_pi_pi_contacts_by_type <- function(frame_coords, tempcell, molecule_dict, cutoff = 5.4) {

    # Initialize list to store ring coordinates
    ring_coords <- list()
    ring_counter <- 1
    
    # Process each molecule
    for (mol_idx in seq_along(molecule_dict)) {
        mol <- molecule_dict[[mol_idx]]
        
        # Get coordinates for each ring and calculate COM
        for (ring_name in c("ring1", "ring2", "ring3", "ring4")) {
            ring_indices <- as.numeric(mol[[ring_name]])
            ring_atom_coords <- frame_coords[ring_indices,]
            ring_coords[[ring_counter]] <- colMeans(ring_atom_coords)  # Calculate COM
            ring_counter <- ring_counter + 1
        }
    }
    
    # Convert list of ring coordinates to matrix
    ring_coords_matrix <- do.call(rbind, ring_coords)
    
    # Create molecule dictionary for rings
    # Each molecule has 4 rings (2 per PHE)
    n_mols <- length(molecule_dict)
    ring_mol_dict <- lapply(1:n_mols, function(i) {
        start_idx <- (i-1)*4 + 1
        end_idx <- i*4
        start_idx:end_idx
    })
    
    # Use C++ function for efficient pairwise contact analysis
    contacts <- contact_analysis_counter_AA(ring_coords_matrix, ring_coords_matrix, cutoff, tempcell, ring_mol_dict)
    
    return(contacts)
}

# Main function to analyze AA contacts by type
# Input:
#   - frame_coords: Full frame coordinates
#   - molecule_dict: Dictionary of molecule information
#   - mycell: Box dimensions for PBC
# Output:
#   - List containing:
#     * pi_pi_contacts: Count of pi-pi interactions
#     * bb_bb_contacts: Count of backbone-backbone interactions
#     * total_contacts: Sum of all contacts
analyze_AA_contacts_by_type <- function(frame_coords, molecule_dict, mycell) {
        # Reshape frame coordinates into nx3 matrix
    frame_coords <- matrix(frame_coords, ncol = 3, byrow = TRUE)
    # Analyze different types of contacts
    pi_pi_contacts <- analyze_AA_pi_pi_contacts_by_type(frame_coords, mycell, molecule_dict)
    bb_bb_contacts <- analyze_AA_bb_bb_contacts_by_type(frame_coords, mycell, molecule_dict)
    
    # Calculate total contacts
    total_contacts <- pi_pi_contacts + bb_bb_contacts
    
    # Return results as a list
    return(list(
        pi_pi_contacts = pi_pi_contacts,
        bb_bb_contacts = bb_bb_contacts,
        total_contacts = total_contacts
    ))
} 