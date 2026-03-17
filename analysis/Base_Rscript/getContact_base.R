if (0) {# carbon
  # 
  # AM1OC
  # CSH C1
  # CSSC C1 C20
  # AM1aC
  # CSH C2
  # CSSC C2 C15
  # THIC
  # CSH C3
  # CSSC C3 C11
  # AM2OC
  # CSH C4
  # CSSC C4 C17
  # PHEaC
  # CSH C5
  # CSSC C5 C16
  # PHEbC
  # CSH C6 C10
  # CSSC C6 C10 C13 C19
  # PHEgC
  # CSH C7 C9 
  # CSSC C7 C9 C12 C18
  # PHEgC
  # CSH C8 
  # CSSC C8 C14
  # 
  # Oxygen
  # AM1O 
  # CSH O1
  # CSSC O1 O4
  # AM2O 
  # CSH O2
  # CSSC O2 O3
  # 
  # Nitrogen
  # AM1N
  # CSH N1
  # CSSC N1 N4
  # AM2N
  # CSH N2
  # CSSC N2 N3
  # 
  # Sulfur
  # CSH S
  # CSSC S1 S2
}
# Define atom types
atom_types <- c("N", "C", "O", "S")

# Initialize a named cutoff matrix with default value 4.6 Å for all pairs
cutoff_matrix <- matrix(4.6, nrow = 4, ncol = 4, dimnames = list(atom_types, atom_types))

# Set specific cutoffs for Carbon-Carbon and Sulfur-Sulfur
cutoff_matrix["C", "C"] <- 5.4  # Carbon-Carbon
cutoff_matrix["S", "S"] <- 6.3  # Sulfur-Sulfur



# Function to return all indices from the same resno & resid
get_atom_indices <- function(index) {
  # Get the resno and resid corresponding to the index
  atom_info <- mypdb$atom[index,]
  resno <- atom_info$resno
  resid <- atom_info$resid
  out<-list()
  # Return all indices for the same resno & resid
  
  out[[1]]<-AtomIndexLibPerRes[[resno]][[resid]]
  out[[2]]<-c(resid,resno)
  out
}
# C++ code to speed up the loop
{cppFunction('int contact_analysis_counter_nonfib(
  IntegerVector RefIndex, 
  List AtomIndexList, 
  NumericMatrix particle, 
  double cutoff, 
  std::string atom_type_j, 
  Function get_atom_indices, 
  int j, 
  int i,
  NumericVector tempdcd, 
  NumericVector tempcell
) {
  int counter = 0;  // Initialize counter before ref_idx loop
  
  // Loop through each index in RefIndex
  for (int ref = 0; ref < RefIndex.size(); ++ref) {
    int ref_idx = RefIndex[ref];  // Set current ref_idx
    
    IntegerVector CompIndex = AtomIndexList[atom_type_j];  // Equivalent to R: CompIndex <- AtomIndexList[[atom_type_j]]
    List ref_inx_residinfo = get_atom_indices(ref_idx);    // Equivalent to R: ref_inx_residinfo <- get_atom_indices(ref_idx)
    
    
    int tempdcdsize = tempdcd.size();  // Total number of points
    int nrows = tempdcdsize / 3;  // Calculate number of rows for the 3-column matrix
    NumericVector origin = particle(ref_idx - 1, _);
    
    // PBC correction
    for (int offset = 0; offset <= 2; ++offset) {
      double tempbox = tempcell[offset];  // Get box size for current axis
      for (int i = offset; i < tempdcdsize; i += 3) {  // Iterate over each coordinate axis
        double diff = tempdcd[i] - origin[offset];  // Difference for PBC
        double absdiff = std::abs(diff);  // Absolute difference
        
        // Apply PBC correction if necessary
        if (absdiff > tempbox / 2) {
          tempdcd[i] -= (diff / absdiff) * tempbox * 
            std::ceil((absdiff - tempbox / 2) / tempbox);
        }
      }
    }
    
    // Convert corrected tempdcd vector into nx3 matrix (particle)
    NumericMatrix particle(nrows, 3);
    for (int i = 0; i < nrows; ++i) {
      for (int j = 0; j < 3; ++j) {
        particle(i, j) = tempdcd[i * 3 + j];
      }
    }
    
    // Compute CloseIndex: indices within cutoff distance
    std::vector<int> CloseIndex;
    for (int k = 0; k < particle.nrow(); ++k) {
      bool within_cutoff = true;
      for (int d = 0; d < particle.ncol(); ++d) {
        if (std::abs(particle(k, d) - particle(ref_idx - 1, d)) > cutoff) {
          within_cutoff = false;
          break;
        }
      }
      if (within_cutoff) CloseIndex.push_back(k + 1);  // Store 1-based index to match R
    }
    
    // Avoid self-counting by checking residue indices
    IntegerVector AvoidIndex;
    if (j == i) {
      IntegerVector residinfo = ref_inx_residinfo[0];
      int max_resid = *std::max_element(residinfo.begin(), residinfo.end());
      for (int idx = 0; idx < CompIndex.size(); ++idx) {
        if (CompIndex[idx] <= max_resid) AvoidIndex.push_back(CompIndex[idx]);
      }
    } else {
      IntegerVector residinfo = ref_inx_residinfo[0];
      for (int idx = 0; idx < CompIndex.size(); ++idx) {
        if (std::find(residinfo.begin(), residinfo.end(), CompIndex[idx]) != residinfo.end()) {
          AvoidIndex.push_back(CompIndex[idx]);
        }
      }
    }
    
    // Filter CompIndex to exclude AvoidIndex and retain only indices in CloseIndex
    std::vector<int> final_CompIndex;
    for (int idx = 0; idx < CompIndex.size(); ++idx) {
      if (std::find(AvoidIndex.begin(), AvoidIndex.end(), CompIndex[idx]) == AvoidIndex.end() &&
          std::find(CloseIndex.begin(), CloseIndex.end(), CompIndex[idx]) != CloseIndex.end()) {
        final_CompIndex.push_back(CompIndex[idx]);
      }
    }
    
    // Loop through each index in final_CompIndex and calculate distance
    for (int ci = 0; ci < final_CompIndex.size(); ci++) {
      int comp_idx = final_CompIndex[ci] - 1;  // Adjust to 0-based index for C++
        
        // Calculate Euclidean distance
      double dist = 0;
      for (int d = 0; d < particle.ncol(); d++) {
        dist += pow(particle(ref_idx - 1, d) - particle(comp_idx, d), 2);  // Equivalent to dist <- euclidean(particle[ref_idx, ], particle[comp_idx, ])
      }
      dist = sqrt(dist);
      
      // Check if within cutoff
      if (dist <= cutoff) {
        counter += 1;
      }
    }
  }
  
  return counter;  // Return the total counter after looping through RefIndex
}
  '
)}
contact_analysis_nonfib <- function(framenum,fiber=TRUE) {
  # Extract DCD for current frame and wrap coordinates based on fiber center
  tempdcd<-mydcd[framenum,] #get coor for current frame 
  tempcell<-mydcdcell[framenum,] #get box size
  #tempdcd<-pbcwrap_fibre(tempdcd,tempcell)
  particle <- matrix(tempdcd, ncol = 3, byrow = TRUE) # Convert coordinates to nx3 matrix
  
  # Initiate an empty matrix, where each row/column corresponds to a unique heavy atom
  num_atoms <- length(heavy_atom_names)
  contact_matrix <- matrix(0, nrow = num_atoms, ncol = num_atoms)
  #debugcounter<-0
  #num_atoms<-1
  
  # Loop through each unique heavy atom as reference atoms
  for (i in 1:num_atoms) {
    atom_type_i <- heavy_atom_names[i] # Assuming fibresel contains atom types
    RefIndex<-AtomIndexList[[atom_type_i]]
    
    for (j in i:num_atoms) { # Avoid double counting by starting loop at i for j
      
      # Determine atom types involved and get cutoff from cutoff_matrix
      atom_type_j <- heavy_atom_names[j]
      cutoff <- cutoff_matrix[AtomTypeLib[[atom_type_i]],AtomTypeLib[[atom_type_j]]] # Extract appropriate cutoff

      # Store the final count in the contact matrix
      contact_matrix[i, j] <- contact_analysis_counter_nonfib(RefIndex, AtomIndexList, particle, cutoff, atom_type_j, get_atom_indices, j, i,tempdcd,tempcell) 
      #contact_matrix[j, i] <- counter # Symmetric matrix
    }
  }
  colnames(contact_matrix)<- names(AtomTypeLib)
  rownames(contact_matrix)<- names(AtomTypeLib)
  # Return the final contact matrix
  return(contact_matrix)
}

#old part of looping in R, converted to cpp for efficiency
contact_analysis_nonfib_part_R<-function(num_atoms=num_atoms,contact_matrix=contact_matrix,particle=particle,i=i,j=j,atom_type_j=atom_type_j,atom_type_i=atom_type_i,RefIndex=RefIndex,cutoff=cutoff)  {
  counter <- 0
      # Loop through each index for reference and comparison atoms
      for (ref_idx in RefIndex) {
        # Call up on index to be compared against
        CompIndex<-AtomIndexList[[atom_type_j]]
        ref_inx_residinfo<-get_atom_indices(ref_idx)
        tempdcd<-pbcwrap(tempdcd,tempcell,particle[ref_idx,])
        particle <- matrix(tempdcd, ncol = 3, byrow = TRUE) # Convert coordinates to nx3 matrix
        CloseIndex <- which(apply(particle, 1, function(row) all(abs(row - particle[ref_idx, ]) <= cutoff)))
        # Avoid self-counting by skipping same residue indices
        if (j == i) {
          #avoid index smaller the max of index in the same resid as ref index 
          AvoidIndex<-which(CompIndex<=max(ref_inx_residinfo[[1]]))
        } else {
          #avoid index  in the same resid as ref index 
          AvoidIndex<-which(CompIndex %in% ref_inx_residinfo[[1]])
        }
        CompIndex<-CompIndex[-AvoidIndex]
        CompIndex<-CompIndex[which(CompIndex %in% CloseIndex)]
        
        for (comp_idx in CompIndex) {
          #debugcounter<-debugcounter+1
          
          # Calculate distance and check if within cutoff
          dist <- euclidean(particle[ref_idx, ], particle[comp_idx, ])
          if (dist <= cutoff) {
            counter <- counter + 1
          }
        }
      }
      counter
      }