# Define atom types
atom_types <- c("N", "C", "O", "S")

# Initialize a named cutoff matrix with default value 4.6 Å for all pairs
cutoff_matrix <- matrix(4.6, nrow = 4, ncol = 4, dimnames = list(atom_types, atom_types))

# Set specific cutoffs for Carbon-Carbon and Sulfur-Sulfur
cutoff_matrix["C", "C"] <- 5.4  # Carbon-Carbon
cutoff_matrix["S", "S"] <- 6.3  # Sulfur-Sulfur


# AtomSelTextList for selection of atoms
AtomSelTextList <- {
  list(
    "CSSC" = list(
      "AM1OC" = c("C1", "C20"),
      "AM1aC" = c("C2", "C15"),
      "THIC" = c("C3", "C11"),
      "AM2OC" = c("C4", "C17"),
      "PHEaC" = c("C5", "C16"),
      "PHEbC" = c("C6", "C10", "C13", "C19"),
      "PHEgC" = c("C7", "C9", "C12", "C18"),
      "PHEdC" = c("C8", "C14"),
      "AM1O" = c("O1", "O4"),
      "AM2O" = c("O2", "O3"),
      "AM1N" = c("N1", "N4"),
      "AM2N" = c("N2", "N3"),
      "S" = c("S1", "S2")
    ),
    "CSH" = list(
      "AM1OC" = c("C1"),
      "AM1aC" = c("C2"),
      "THIC" = c("C3"),
      "AM2OC" = c("C4"),
      "PHEaC" = c("C5"),
      "PHEbC" = c("C6", "C10"),
      "PHEgC" = c("C7", "C9"),
      "PHEdC" = c("C8"),
      "AM1O" = c("O1"),
      "AM2O" = c("O2"),
      "AM1N" = c("N1"),
      "AM2N" = c("N2"),
      "S" = c("S")
    )
  )
}

# AtomTypeLib for atom types
AtomTypeLib <- {
  list(
    "AM1OC" = "C",
    "AM1aC" = "C",
    "THIC" = "C",
    "AM2OC" = "C",
    "PHEaC" = "C",
    "PHEbC" = "C",
    "PHEgC" = "C",
    "PHEdC" = "C",
    "AM1O" = "O",
    "AM2O" = "O",
    "AM1N" = "N",
    "AM2N" = "N",
    "S" = "S"
  )
}
heavy_atom_names<-names(AtomTypeLib)
# Create an empty AtomIndexListPerRes
AtomIndexListPerRes <- list()

# Loop through each Resname (CSH, CSSC) and AtomName (AM2N, AM1N, etc.)
for (Resname in names(AtomSelTextList)) {
  AtomIndexListPerRes[[Resname]] <- list()
  
  for (AtomName in names(AtomSelTextList[[Resname]])) {
    # Use atomselect on mypdb to get the index
    AtomIndexListPerRes[[Resname]][[AtomName]] <- atom.select(
      mypdb, 
      resid = Resname, 
      elety = AtomSelTextList[[Resname]][[AtomName]]
    )$atom
  }
}

# Turn AtomIndexListPerRes into AtomIndexList by merging entries with the same AtomName, regardless of Resname
AtomIndexList <- list()

for (Resname in names(AtomIndexListPerRes)) {
  for (AtomName in names(AtomIndexListPerRes[[Resname]])) {
    if (!is.null(AtomIndexList[[AtomName]])) {
      # Merge with existing entries
      AtomIndexList[[AtomName]] <- c(AtomIndexList[[AtomName]], AtomIndexListPerRes[[Resname]][[AtomName]])
    } else {
      # Create new entry if it doesn't exist
      AtomIndexList[[AtomName]] <- AtomIndexListPerRes[[Resname]][[AtomName]]
    }
  }
}
AtomIndexLib<-list()
for (resno in 1:2000) {
  # Use atomselect to select and store the index
  selection <- atom.select(
    mypdb,
    resid = c("CSH","CSSC"),
    resno = resno,
    "noh"
  )
  # Only store the selection if there are atoms selected
  if (length(selection$atom)!=0) {
    # Ensure the AtomIndexLib list has an entry for the current resno
    
    AtomIndexLib[[resno]] <- selection$atom
  }
}

# Assuming particle, cutoff_matrix, AtomIndexLib, and AtomTypeLib are already defined


cppFunction('
    NumericMatrix find_contacts_for_residue( 
                                            List AtomIndexLib, 
                                            NumericMatrix particle,
                                            NumericMatrix cutoff_matrix, 
                                            List AtomIndexList, 
                                            List AtomTypeLib) {
                                            // set max num of res
                                            int max_resno = AtomIndexLib.size();
      // Initialize an empty edge matrix with two columns
      NumericMatrix edge_matrix(0, 2);
      
      // Map atom types to matrix indices
      std::map<std::string, int> type_to_index = {{"N", 0}, {"C", 1}, {"O", 2}, {"S", 3}};
      
      // Loop over each subsequent residue for resno_i
      for (int resno_i = 1; resno_i <= max_resno; ++resno_i) {
        
        IntegerVector atoms_i_list = AtomIndexLib[resno_i - 1];
        IntegerVector atoms_i = Rcpp::as<IntegerVector>(Rcpp::wrap(atoms_i_list)); // Flattened atoms_i list
        
        
        // Loop over each subsequent residue for resno_j
        for (int resno_j = resno_i + 1; resno_j <= max_resno; ++resno_j) {
          
          IntegerVector atoms_j_list = AtomIndexLib[resno_j - 1];
          IntegerVector atoms_j = Rcpp::as<IntegerVector>(Rcpp::wrap(atoms_j_list)); // Flattened atoms_j list
          
          // Loop over each atom in atoms_i
          bool contact_found = false;
          for (int a = 0; a < atoms_i.size(); ++a) {
            int atom_i = atoms_i[a] - 1; // Convert to 0-based index
            
            // Check bounds for atom_i in particle
            if (atom_i < 0 || atom_i >= particle.nrow()) {
              continue;
            }
            
            // Identify type of atom_i
            std::string type_i;
            CharacterVector atom_names = AtomIndexList.names();
            for (int k = 0; k < atom_names.size(); ++k) {
              std::string name = Rcpp::as<std::string>(atom_names[k]);
              IntegerVector atom_list = AtomIndexList[name];
              if (std::find(atom_list.begin(), atom_list.end(), atoms_i[a]) != atom_list.end()) {
                type_i = Rcpp::as<std::string>(AtomTypeLib[name]);
                break;
              }
            }
            
            // Convert type_i to index
            if (type_to_index.find(type_i) == type_to_index.end()) {
              continue;
            }
            int type_i_index = type_to_index[type_i];
            
            
            // Loop over each atom in atoms_j
            for (int b = 0; b < atoms_j.size(); ++b) {
              int atom_j = atoms_j[b] - 1;  // Convert to 0-based index
              
              // Check bounds for atom_j in particle
              if (atom_j < 0 || atom_j >= particle.nrow()) {
                continue;
              }
              
              // Identify type of atom_j
              std::string type_j;
              CharacterVector atom_names = AtomIndexList.names();
              for (int k = 0; k < atom_names.size(); ++k) {
                std::string name = Rcpp::as<std::string>(atom_names[k]);
                IntegerVector atom_list = AtomIndexList[name];
                if (std::find(atom_list.begin(), atom_list.end(), atoms_j[b]) != atom_list.end()) {
                  type_j = Rcpp::as<std::string>(AtomTypeLib[name]);
                  break;
                }
              }
              
              // Convert type_j to index
              if (type_to_index.find(type_j) == type_to_index.end()) {
                continue;
              }
              int type_j_index = type_to_index[type_j];
              
              // Determine the cutoff for the atom pair
              double cutoff = cutoff_matrix(type_i_index, type_j_index);
              
              // Calculate the Euclidean distance between atom_i and atom_j
              double dist = 0.0;
              for (int d = 0; d < 3; ++d) {
                dist += pow(particle(atom_i, d) - particle(atom_j, d), 2.0);
              }
              dist = sqrt(dist);
              
              // Check if the distance is within the cutoff
              if (dist <= cutoff) {
                // Print debug message for contact detection
                //Rcpp::Rcout << "Contact found between residues " << resno_i << " and " << resno_j << std::endl;
                
                // Contact found, add new row to edge_matrix
                NumericMatrix new_edge_matrix(edge_matrix.nrow() + 1, 2);
                for (int row = 0; row < edge_matrix.nrow(); ++row) {
                  new_edge_matrix(row, 0) = edge_matrix(row, 0);
                  new_edge_matrix(row, 1) = edge_matrix(row, 1);
                }
                new_edge_matrix(edge_matrix.nrow(), 0) = resno_i;
                new_edge_matrix(edge_matrix.nrow(), 1) = resno_j;
                edge_matrix = new_edge_matrix; // Update edge_matrix with the new row
                
                contact_found = true;
                break;
              }
            }
            if (contact_found) break;
          }
        }
      }
        return edge_matrix;
      
    }
')

# Run the function with specific parameters to test contact detection and edge list population




# Main R function to call the C++ function and prepare data
record_contact_edges <- function(particle, cutoff_matrix, AtomIndexLib, AtomTypeLib) {
  # Ensure the particle data is in an appropriate format (e.g., matrix with atom types as a column)
  # Call the C++ function
  edge_list <-  find_contacts_for_residue(AtomIndexLib, particle, cutoff_matrix, AtomIndexList, AtomTypeLib)
  return(edge_list)
}

getEdge_per_frame <- function(frame) {
  tempdcd<-mydcd[frame,] #get coor for current frame 
  tempcell<-mydcdcell[frame,] #get box size
  tempdcd<-pbcwrap_fibre(tempdcd,tempcell) #wrap coordinates
  #edges<-getedge_cutoff(frame)
  particle<-matrix(tempdcd,ncol = 3,byrow = TRUE) #convert coord to a 3 column matrix (x, y, z) 
  edges<-record_contact_edges(particle, cutoff_matrix, AtomIndexLib, AtomTypeLib)
  edges<-as.matrix(cbind(edges,1))
  edges<-rbind(edges,edges[,c(2,1,3)])
  attr(edges,"n")<-length(AtomIndexLib)
  edges
}
