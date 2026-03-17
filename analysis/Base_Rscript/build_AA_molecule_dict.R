# Function to build a dictionary of molecule information from AA topology
# This function processes the topology to identify:
# - Number of molecules
# - For each molecule: residue numbers, atom indices for rings and backbone
build_AA_molecule_dict <- function(topology, unresid, is_acidic = FALSE) {
    # Get unique molecule numbers (each FSSF has a unique resno)
    mol_numbers <- unique(topology$resno[topology$resid == unresid])
    n_mols <- length(mol_numbers)
    
    # Define ring atoms based on pH condition
    if (is_acidic) {
        # Acidic condition ring definitions
        ring1_atoms <- c("C20", "C17", "C19", "C21", "C18", "C15")  # PHE1 ring 1
        # VMD selection: name C20 C17 C19 C21 C18 C15
        
        ring2_atoms <- c("C27", "C22", "C26", "C23", "C25", "C24")  # PHE2 ring 1
        # VMD selection: name C27 C22 C26 C23 C25 C24
        
        ring3_atoms <- c("C36", "C35", "C34", "C33", "C32", "C29")  # PHE2 ring 2
        # VMD selection: name C36 C35 C34 C33 C32 C29
        
        ring4_atoms <- c("C38", "C37", "C39", "C40", "C31", "C13")  # PHE1 ring 2
        # VMD selection: name C38 C37 C39 C40 C31 C13
    } else {
        # Basic condition ring definitions
        ring1_atoms <- c("C14", "C16", "C20", "C17", "C19", "C18")  # PHE1 ring 1
        # VMD selection: name C14 C16 C20 C17 C19 C18
        
        ring2_atoms <- c("C28", "C31", "C35", "C32", "C34", "C33")  # PHE2 ring 1
        # VMD selection: name C28 C31 C35 C32 C34 C33
        
        ring3_atoms <- c("C30", "C36", "C40", "C37", "C39", "C38")  # PHE1 ring 2
        # VMD selection: name C30 C36 C40 C37 C39 C38
        
        ring4_atoms <- c("C21", "C22", "C26", "C23", "C25", "C24")  # PHE2 ring 2
        # VMD selection: name C21 C22 C26 C23 C25 C24
    }
    
    # Initialize dictionary
    molecule_dict <- list()
    
    # Process each molecule
    for (mol in 1:n_mols) {
        # Get indices for this molecule
        mol_indices <- which(topology$resid == unresid & topology$resno == mol_numbers[mol])
        
        # Get atoms for this molecule
        mol_atoms <- topology[mol_indices,]
        
        # Initialize molecule dictionary
        mol_dict <- list()
        
        # Store ring atoms
        mol_dict[["ring1"]] <- mol_atoms$eleno[mol_atoms$elety %in% ring1_atoms]
        mol_dict[["ring2"]] <- mol_atoms$eleno[mol_atoms$elety %in% ring2_atoms]
        mol_dict[["ring3"]] <- mol_atoms$eleno[mol_atoms$elety %in% ring3_atoms]
        mol_dict[["ring4"]] <- mol_atoms$eleno[mol_atoms$elety %in% ring4_atoms]
        
        # Store backbone atoms (N1-N6 and O1-O4)
        mol_dict[["N"]] <- mol_atoms$eleno[startsWith(mol_atoms$elety, "N")]
        mol_dict[["O"]] <- mol_atoms$eleno[startsWith(mol_atoms$elety, "O")]
        
        # Add molecule to dictionary
        molecule_dict[[mol]] <- mol_dict
    }
    
    return(molecule_dict)
} 