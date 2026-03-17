# Function to analyze contacts in a single frame
# This function processes one frame of the trajectory and returns contact data
# Input:
#   - frame_coords: Matrix of coordinates for the current frame
#   - molecule_dict: Dictionary mapping molecule indices to coordinate indices
# Output:
#   - Data frame with contact information for this frame
analyze_CG_frame <- function(frame_coords, molecule_dict) {
    # Analyze contacts by type
    contact_results <- analyze_CG_contacts_by_type(frame_coords, molecule_dict)
    
    # Return results as a data frame
    return(data.frame(
        pi_pi_contacts = contact_results$pi_pi_contacts,
        bb_bb_contacts = contact_results$bb_bb_contacts,
        total_contacts = contact_results$total_contacts
    ))
}

# Main function to analyze contacts in a trajectory
# This function processes the entire trajectory and returns contact data for all frames
# Input:
#   - simtraj: Trajectory object from read.trj
# Output:
#   - Data frame with contact information for all frames
analyze_CG_contacts <- function(simtraj) {
    # Get number of frames
    n_frames <- dim(simtraj$coord)[3]
    
    # Build molecule dictionary from topology
    molecule_dict <- build_CG_molecule_dict(simtraj$top)
    
    # Process frames in parallel
    contact_data <- foreach(frame = 1:n_frames, .combine = rbind) %dopar% {
        # Extract coordinates for current frame
        frame_coords <- simtraj$coord[,,frame]
        
        # Analyze contacts in this frame
        frame_results <- analyze_CG_frame(frame_coords, molecule_dict)
        
        # Add frame number
        frame_results$frame <- frame
        
        return(frame_results)
    }
    
    return(contact_data)
} 