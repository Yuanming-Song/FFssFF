#resname to change elety
tempreslist<-c("ILE", "PHE")
tempreslist_elety <- c("I", "F")
for (temprescount in 1:length(tempreslist)){
  tempresind <- which(mypdb$atom$resid == tempreslist[temprescount])
  mypdb$atom$elety[tempresind] <- paste0(mypdb$atom$elety[tempresind],tempreslist_elety[temprescount])
}
# Get indices for IF residues ("ILE" or "PHE") and sort them
IF_ind <- which(mypdb$atom$resid %in% tempreslist)
IF_indices_sorted <- sort(IF_ind)

# Calculate the expected number of atoms per IF group:
# atoms_per_IF = 2 * (total number of IF atoms) / (number of unique IF residue numbers)
num_unique <- length(unique(mypdb$atom$resno[IF_ind]))
atoms_per_IF <- as.integer(2 * length(IF_ind) / num_unique)

# Set new residue numbers to start from one more than the maximum FSSF resno
start <- max(mypdb$atom$resno[mypdb$atom$resid == "FSSF"]) + 1

# Divide the sorted IF indices into blocks of size 'atoms_per_IF'
n_groups <- ceiling(length(IF_indices_sorted) / atoms_per_IF)
for (i in seq_len(n_groups)) {
  block_start <- (i - 1) * atoms_per_IF + 1
  block_end <- min(i * atoms_per_IF, length(IF_indices_sorted))
  indices_block <- IF_indices_sorted[block_start:block_end]
  # Assign a new residue number for this block
  mypdb$atom$resno[indices_block] <- start + i - 1
}
#change resname ILE PHE to IF
mypdb$atom$resid[IF_ind] <- "IF"
# AtomSelTextList for selection of atoms
AtomSelTextList <- {
  list(
    "FSSF"=list(
      "PHE1aC" = c("C1", "C11"),
      "PHE1bC" = c("C13", "C12"),
      "PHE1gC" = c("C14", "C30"),
      "PHE1dC" = c("C16", "C20", "C36", "C40"),
      "PHE1eC" = c("C37", "C39", "C17", "C19"),
      "PHE1zC" = c("C18", "C38"),
      "PHE2aC" = c("C10", "C3"),
      "PHE2bC" = c("C15", "C27"),
      "PHE2gC" = c("C28", "C21"),
      "PHE2dC" = c("C31", "C35", "C22", "C26"),
      "PHE2eC" = c("C32", "C34", "C23", "C25"),
      "PHE2zC" = c("C33", "C24"),
      "CYSC1" = c("C5", "C8"),
      "CYSC2" = c("C6", "C7"),
      "CYSS" = c("S1", "S2"),
      "PHE1BOC" = c("C2", "C11"),
      "PHE2BOC" = c("C4", "C9"),
      "PHE1BO" = c("O1", "O4"),
      "PHE2BO" = c("O2", "O3"),
      "PHE1N" = c("N1", "N6"),
      "PHE2N" = c("N2", "N5"),
      "CYSN" = c("N3", "N4")
    ),
    # Define the mapping list for the "IF" group
   #need to make sure these names are corrected 
    "IF" = list(
      "PHEaC"  = c("CAF"),
      "PHEbC"  = c("CBF"),
      "PHEgC"  = c("CGF"),
      "PHEdC"  = c("CD1F", "CD2F"),
      "PHEeC"  = c("CE1F", "CE2F"),
      "PHEzC"  = c("CZF"),
      "PHEBOC" = c("CF"),
      "PHEBO"  = c("OT1F", "OT2F"),
      "PHEN"   = c("NF"),
      "ILEN"   = c("NI"),
      "ILEBO"  = c("OI"),
      "ILEBOC" = c("CI"),
      "ILEaC"  = c("CAI"),
      "ILEbC" = c("CBI"),
      "ILEgC1" = c("CG1I"),
      "ILEgC2" = c("CG2I"),
      "ILEdC"  = c("CDI")
    )
  )
  
}
FFssFF_type <- {
  list(
    "SideChain" = c("PHE1bC", "PHE1gC", "PHE1dC", "PHE1eC", "PHE1zC", "PHE2bC", "PHE2gC", "PHE2dC", "PHE2eC", "PHE2zC", "CYSC1", "CYSC2", "CYSS"),
    "Backbone" = c("PHE1BOC", "PHE2BOC", "PHE1BO", "PHE2BO", "PHE1N", "PHE2N", "CYSN","PHE1aC","PHE2aC")
  )
}
FI_type <- {
  list(
    "SideChain" = c( "PHEbC", "PHEgC", "PHEdC", "PHEeC", "PHEzC", "ILEbC","ILEgC1", "ILEgC2", "ILEdC"),
    "Backbone"  = c("PHEBOC", "PHEBO", "PHEN", "ILEN", "ILEBO", "ILEBOC","ILEaC", "PHEaC")
  )
}
# AtomTypeLib for atom types
AtomTypeLib <- list()
for (resname in names(AtomSelTextList)) {
  for (name in names(AtomSelTextList[[resname]])) {
    if (grepl("C$", name) || grepl("C[0-9]+$", name)) {
      AtomTypeLib[[name]] <- "C"
    } else if (grepl("O$", name)) {
      AtomTypeLib[[name]] <- "O"
    } else if (grepl("N$", name)) {
      AtomTypeLib[[name]] <- "N"
    } else if (grepl("S$", name)) {
      AtomTypeLib[[name]] <- "S"
    }
  }
}
# Get all heavy atom names from AtomTypeLib for all residue types
heavy_atom_names <- names(AtomTypeLib)
# Create an empty AtomIndexListPerRes
AtomIndexListPerRes <- list()

# Loop through each Resname  and AtomName (AM2N, AM1N, etc.)
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

# Create an empty list AtomIndexLib
AtomIndexLib<-list()
for (resno in 1:2000) {
  # Use atomselect to select and store the index
  selection <- atom.select(
    mypdb,
    resid = c("IF","FSSF"),
    resno = resno,
    "noh"
  )
  # Only store the selection if there are atoms selected
  if (length(selection$atom)!=0) {
    # Ensure the AtomIndexLib list has an entry for the current resno
    
    AtomIndexLib[[resno]] <- selection$atom
  }
}

# Create AtomIndexLibPerRes (old version) for compatibility
AtomIndexLibPerRes <- list()

# Loop through 1:2000 as resno
for (resno in 1:200) {
  # Loop through each resid (FSSF, IF)
  for (resid in names(AtomSelTextList)) {
    # Use atomselect to select and store the index
    selection <- atom.select(
      mypdb,
      resid = resid,
      resno = resno,
      elety = unlist(AtomSelTextList[[resid]])
    )
    # Only store the selection if there are atoms selected
    if (!is.null(selection$atom)) {
      # Ensure the AtomIndexLibPerRes list has an entry for the current resno
      if (!exists(as.character(resno), where = AtomIndexLibPerRes)) {
        AtomIndexLibPerRes[[as.character(resno)]] <- list()
      }
      # Store the atom indices in AtomIndexLibPerRes for the corresponding resno and resid
      AtomIndexLibPerRes[[as.character(resno)]][[resid]] <- selection$atom
    }
  }
}



