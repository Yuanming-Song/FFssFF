# Define atom types
atom_types <- c("N", "C", "O", "S")

# Initialize a named cutoff matrix with default value 4.6 Å for all pairs
cutoff_matrix <- matrix(4.6, nrow = 4, ncol = 4, dimnames = list(atom_types, atom_types))

# Set specific cutoffs for Carbon-Carbon and Sulfur-Sulfur
cutoff_matrix["C", "C"] <- 5.4  # Carbon-Carbon
cutoff_matrix["S", "S"] <- 6.3  # Sulfur-Sulfur

cppFunction('
    NumericMatrix find_contacts_for_residue(List AtomIndexLib, 
                                          NumericMatrix particle,
                                          NumericMatrix cutoff_matrix, 
                                          List AtomIndexList, 
                                          List AtomTypeLib,
                                          NumericVector tempcell) {
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
              
              // Calculate the Euclidean distance between atom_i and atom_j with PBC
              double dist = 0.0;
              for (int d = 0; d < 3; ++d) {
                double diff = particle(atom_i, d) - particle(atom_j, d);
                double box_size = tempcell[d];
                double half_box = box_size / 2.0;
                if (std::abs(diff) > half_box) {
                  diff -= box_size * std::round(diff / box_size);
                }
                dist += diff * diff;
              }
              dist = sqrt(dist);
              
              // Check if the distance is within the cutoff
              if (dist <= cutoff) {
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

# Main R function to call the C++ function and prepare data
record_contact_edges <- function(particle, cutoff_matrix, AtomIndexLib, AtomTypeLib, tempcell) {
  # Call the C++ function
  edge_list <- find_contacts_for_residue(AtomIndexLib, particle, cutoff_matrix, AtomIndexList, AtomTypeLib, tempcell)
  return(edge_list)
}

getEdge_per_frame <- function(frame) {
  tempdcd<-mydcd[frame,] #get coor for current frame 
  tempcell<-mydcdcell[frame,] #get box size
  #tempdcd<-pbcwrap(tempdcd,tempcell) #wrap coordinates
  particle<-matrix(tempdcd,ncol = 3,byrow = TRUE) #convert coord to a 3 column matrix (x, y, z) 
  edges<-record_contact_edges(particle, cutoff_matrix, AtomIndexLib, AtomTypeLib, tempcell)
  edges<-as.matrix(cbind(edges,1))
  edges<-rbind(edges,edges[,c(2,1,3)])
  attr(edges,"n")<-length(AtomIndexLib)
  edges
}
