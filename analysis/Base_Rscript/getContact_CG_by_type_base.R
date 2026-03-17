#assuming x number copies of molecules
#each molecuels is similar to peptide, with 6 residues
#it's CG version, so only heavy atoms exists and we care about them all
#Break interactions into 2 types, pi-pi int, backbone-backbone int, each with a distinct distance cutoff
#each residue has 1 backbone beads, phe has 3 sidechain beads
#during contact analysis, backbone is easy, loop through them all
#for side chain, first get COM of each ring, then loop through them all
#for pair wise interactions looping, between 2 backbone beads or 2 phe COM,  should be handled with a c++ function
#now work on it based of getContact_base.R

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
int contact_analysis_counter_CG(
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
                            // Move the atom to the nearest periodic image
                            // Formula explanation:
                            // 1. (diff / abs_diff) gives the sign (-1 or +1) of the displacement
                            // 2. (abs_diff - box_size / 2) / box_size calculates how many box lengths we need to move
                            // 3. std::ceil() rounds up to nearest integer number of box lengths
                            // 4. Multiply by box_size to get actual distance to move
                            // 5. Apply the sign to determine direction
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

# Function to calculate center of mass for a group of beads
# This utility function computes the center of mass for a set of coordinates
# Input:
#   - coords: Matrix of coordinates (Nx3)
# Output:
#   - Vector of COM coordinates (length 3)
calculate_com <- function(coords) {
    colMeans(coords)
}

# Function to analyze pi-pi interactions between PHE sidechains
# This function specifically handles aromatic ring interactions in the CG model
# It first calculates the COM of each PHE ring, then analyzes contacts between rings
# Input:
#   - frame_coords: Full frame coordinates
#   - tempcell: Box dimensions for PBC
#   - molecule_dict: Dictionary of molecule information
#   - cutoff: Distance cutoff for pi-pi interactions (default 5.4 Å)
# Output:
#   - Integer count of pi-pi contacts
analyze_CG_pi_pi_contacts_by_type <- function(frame_coords, tempcell, molecule_dict, cutoff = 0.6) {
    # Initialize list to store PHE ring COMs
    phe_coms <- list()
    com_counter <- 1
    
    # Process each molecule
    for (mol_idx in seq_along(molecule_dict)) {
        mol <- molecule_dict[[mol_idx]]
        
        # Find all PHE residues in this molecule
        phe_res <- grep("^PHE[0-9]+_SC", names(mol), value = TRUE)
        
        # Group SC beads by PHE residue number
        phe_groups <- split(phe_res, sub("_SC[0-9]+$", "", phe_res))
        
        # For each PHE residue, calculate one COM using all its SC beads
        for (phe_group in phe_groups) {
            # Get all SC bead indices for this PHE
            sc_indices <- as.numeric(mol[phe_group])
            
            # Extract coordinates for these beads
            sc_coords <- frame_coords[sc_indices,]
            
            # Calculate COM for this PHE ring
            phe_coms[[com_counter]] <- calculate_com(sc_coords)
            com_counter <- com_counter + 1
        }
    }
    
    # Convert list of COMs to matrix
    phe_coms_matrix <- do.call(rbind, phe_coms)
    
    # Create molecule dictionary for PHE COMs
    # Each molecule has 4 PHE residues, so we group COMs accordingly
    n_mols <- length(molecule_dict)
    phe_mol_dict <- lapply(1:n_mols, function(i) {
        start_idx <- (i-1)*4 + 1
        end_idx <- i*4
        start_idx:end_idx
    })
    
    # Use C++ function for efficient pairwise contact analysis
    contacts <- contact_analysis_counter_CG(phe_coms_matrix, phe_coms_matrix, cutoff, tempcell, phe_mol_dict)
    
    return(contacts)
}

# Function to analyze backbone-backbone interactions
# This function handles contacts between backbone beads in the CG model
# It identifies backbone beads and analyzes their pairwise interactions
# Input:
#   - frame_coords: Full frame coordinates
#   - tempcell: Box dimensions for PBC
#   - molecule_dict: Dictionary of molecule information
#   - cutoff: Distance cutoff for backbone interactions (default 4.6 Å)
# Output:
#   - Integer count of backbone-backbone contacts
analyze_CG_bb_bb_contacts_by_type <- function(frame_coords, tempcell, molecule_dict, cutoff = 0.5) {
    # Initialize list to store backbone coordinates
    bb_coords <- list()
    bb_counter <- 1
    
    # Process each molecule
    for (mol_idx in seq_along(molecule_dict)) {
        mol <- molecule_dict[[mol_idx]]
        
        # Find all BB beads in this molecule
        bb_res <- grep("_BB$", names(mol), value = TRUE)
        
        # Get coordinates for each BB bead
        for (res in bb_res) {
            bb_index <- as.numeric(mol[res])
            bb_coords[[bb_counter]] <- frame_coords[bb_index,]
            bb_counter <- bb_counter + 1
        }
    }
    
    # Convert list of BB coordinates to matrix
    bb_coords_matrix <- do.call(rbind, bb_coords)
    
    # Create molecule dictionary for BB beads
    # Each molecule has 4 BB beads (one per PHE)
    n_mols <- length(molecule_dict)
    bb_mol_dict <- lapply(1:n_mols, function(i) {
        start_idx <- (i-1)*4 + 1
        end_idx <- i*4
        start_idx:end_idx
    })
    
    # Use C++ function for efficient pairwise contact analysis
    contacts <- contact_analysis_counter_CG(bb_coords_matrix, bb_coords_matrix, cutoff, tempcell, bb_mol_dict)
    
    return(contacts)
}

# Function to calculate box dimensions from frame coordinates
# This function computes the box dimensions by finding the range of coordinates in each dimension
# Input:
#   - frame_coords: Matrix of coordinates (Nx3)
# Output:
#   - Vector of box dimensions [x, y, z]
calculate_box_dimensions <- function(frame_coords) {
    # Calculate min and max for each dimension
    x_range <- range(frame_coords[,1])
    y_range <- range(frame_coords[,2])
    z_range <- range(frame_coords[,3])
    
    # Return box dimensions as vector
    c(x_range[2] - x_range[1],
      y_range[2] - y_range[1],
      z_range[2] - z_range[1])
}

# Main function to analyze CG contacts by type
# This is the primary function that coordinates the analysis of different contact types
# It combines pi-pi and backbone-backbone contact analysis and calculates overall metrics
# Input:
#   - frame_coords: Full frame coordinates
#   - molecule_dict: Dictionary of molecule information
# Output:
#   - List containing:
#     * pi_pi_contacts: Count of pi-pi interactions
#     * bb_bb_contacts: Count of backbone-backbone interactions
#     * total_contacts: Sum of all contacts
analyze_CG_contacts_by_type <- function(frame_coords, molecule_dict) {
    # Calculate box dimensions from frame coordinates
    tempcell <- calculate_box_dimensions(frame_coords)
    
    # Analyze different types of contacts
    pi_pi_contacts <- analyze_CG_pi_pi_contacts_by_type(frame_coords, tempcell, molecule_dict)
    bb_bb_contacts <- analyze_CG_bb_bb_contacts_by_type(frame_coords, tempcell, molecule_dict)
    
    # Calculate total contacts
    total_contacts <- pi_pi_contacts + bb_bb_contacts
    
    # Return results as a list
    return(list(
        pi_pi_contacts = pi_pi_contacts,
        bb_bb_contacts = bb_bb_contacts,
        total_contacts = total_contacts
    ))
} 