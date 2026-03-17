#base function for conducting pbc correction, chunk and cores optional for parallel computation
pbcwrap_base<-function(tempdcd_base,origin,tempcell,chunk=1,cores=1) { 
  tempdcd_basesize<-length(tempdcd_base) #how many points in total
  chunksize<-round(tempdcd_basesize/3/cores) #for parallel comp, makesure each chunck include complete xyz coord
  ini<-((chunk-1)*chunksize)*3+1 #first atom
  if(chunk<cores) {
    fini<-chunk*chunksize*3 #last atom in this chunk
  } else {
    fini<-tempdcd_basesize #special treat for last chunk
  }
  for (offset in 0:2) { 
    #offset for x y z, x0 y1 z2
    tempbox <- tempcell[offset+1] #get box size
    # Ensure the sequence goes from smaller to larger number
    start_idx <- min(ini+offset, fini-2+offset)
    end_idx <- max(ini+offset, fini-2+offset)
    indices <- seq(start_idx, end_idx, 3) #indx for x, y or z
    diff <- tempdcd_base[indices] - origin[offset+1] #difference in z 
    absdiff <- abs(diff) #absolute value
    update_indices <- absdiff > tempbox/2 #which ones need correction 
    #pbc correction
    tempdcd_base[indices[update_indices]] <- tempdcd_base[indices[update_indices]] -
      (diff[update_indices]/absdiff[update_indices]) * 
      tempbox * 
      ceiling(
        (absdiff[update_indices] - tempbox/2) / 
          tempbox
      )
  }
  tempdcd_base[seq(ini,fini,1)] #output for this chunk
}

pbcwrap<-function(tempdcd_in,tempcell,origin=c(0,0,0)) { 
  tempdcd_in<-pbcwrap_base(tempdcd_in,origin,tempcell) #wrap everything around COM of fibre

  tempdcd_in #return wrapped coordinates
}

cross_manual <- function(v1, v2,matrix=TRUE) {
  if (TRUE) {
    cbind(v1[,2] * v2[,3] - v1[,3] * v2[,2],
          v1[,3] * v2[,1] - v1[,1] * v2[,3],
          v1[,1] * v2[,2] - v1[,2] * v2[,1])
  } else {
    cbind(v1[2] * v2[3] - v1[3] * v2[2],
          v1[3] * v2[1] - v1[1] * v2[3],
          v1[1] * v2[2] - v1[2] * v2[1])
  }
}
#get distance between two vector
euclidean <- function(a, b) {
  sqrt(sum((a - b)^2))  # Calculate Euclidean distance
}

pbc <- function(comi, cur, currbox) {
  # Loop through x, y, and z coordinates
  for (l in 1:3) {
    diff <- cur[l] - comi[l]  # Calculate the difference between the current coordinate and the center of mass coordinate
    absdiff <- abs(diff)  # Calculate the absolute value of the difference
    # Update coordinate only when necessary to account for periodic boundary conditions
    if (absdiff > currbox[l]/2) {
      cur[l] <- cur[l] - (diff/abs(diff)) * currbox[l] * ceiling((abs(diff) - currbox[l]/2) / currbox[l])
    } else {
      cur[l] <- cur[l]
    }
  }
  # Return updated current coordinates
  cur
}

#base function for conducting pbc correction, chunk and cores optional for parallel computation
pbcwrap_fibre_base<-function(tempdcd_base,z,tempcell,chunk=1,cores=1) { 
  origin<-c(0,0,z)
  tempdcd_basesize<-length(tempdcd_base) #how many points in total
  chunksize<-round(tempdcd_basesize/3/cores) #for parallel comp, makesure each chunck include complete xyz coord
  ini<-((chunk-1)*chunksize)*3+1 #first atom
  if(chunk<cores) {
    fini<-chunk*chunksize*3 #last atom in this chunk
  } else {
    fini<-tempdcd_basesize #special treat for last chunk
  }
  for (offset in 0:2) { #offset for x y z, x0 y1 z2
    tempbox <- tempcell[offset+1] #get box size
    indices <- seq(ini+offset, fini-2+offset, 3) #indx for x, y or z
    diff <- tempdcd_base[indices] - origin[offset+1] #difference in z 
    absdiff <- abs(diff) #absolute value
    update_indices <- absdiff > tempbox/2 #which ones need correction 
    #pbc correction
    tempdcd_base[indices[update_indices]] <- tempdcd_base[indices[update_indices]] -
      (diff[update_indices]/absdiff[update_indices]) * 
      tempbox * 
      ceiling(
        (absdiff[update_indices] - tempbox/2) / 
          tempbox
      )
  }
  tempdcd_base[seq(ini,fini,1)] #output for this chunk
}

#include z dir correction for fibre
pbcwrap_fibre<-function(tempdcd_in,tempcell) { 
  fibresel<-atom.select(mypdb,"noh",resid=c("CSSC","CSH")) #heavy atoms only
  fibrecoor<-pbcwrap_fibre_base(tempdcd_in[fibresel$xyz],0,tempcell) #wrap fibre heavy atoms for now
  com<-com.xyz(fibrecoor) #get center of mass
  tempdcd_in<-pbcwrap_fibre_base(tempdcd_in,com[3],tempcell) #wrap everything around COM of fibre
  tempdcd_in[seq(3,length(tempdcd_in),3)]<-  tempdcd_in[seq(3,length(tempdcd_in),3)]-com[3] #shift to origin
  tempdcd_in #return wrapped coordinates
}


#get distribution by moiety
getcoordist_fibre_by_moiety<-function(frame) { 
  tempdcd<-mydcd[frame,] #get coor for current frame 
  tempcell<-mydcdcell[frame,] #get box size
  tempdcd<-pbcwrap_fibre(tempdcd,tempcell) #wrap coordinates
  particle<-matrix(tempdcd,ncol = 3,byrow = TRUE) #convert coord to a 3 column matrix (x, y, z) 
  zhisout<-c() #empty output name 
  rhisout<-c()
  for (name in names(MoietyCoorList)) {
    MoietyResid<-MoietyCoorList[[name]]

    FibreCoorX<-particle[MoietyResid,1] #x
    FibreCoorY<-particle[MoietyResid,2] #y
    fibrecoorR<-sqrt(FibreCoorX^2+FibreCoorY^2) #radial distance
    FibreCoorZ<-particle[MoietyResid,3] #z 
    
    FibreCoorZfib<-FibreCoorZ#[which(fibrecoorR<27)]  #in case only need to fibre 
    fibrecoorRfib<-fibrecoorR#[which(abs(FibreCoorZ)<100)]
    
    coorzhis<-table(cut(FibreCoorZfib,breaks = seq(-180,180,1))) #get raw dist
    coorrhis<-table(cut(fibrecoorRfib,breaks = seq(0,40,1))) #get raw dist
    coorrhis<-table(cut(fibrecoorRfib,breaks = seq(0,71,1))) #get raw dist
    zhisout<-cbind(zhisout,coorzhis) #combine output 
    rhisout<-cbind(rhisout,coorrhis) #combine output 
  }
  out<-list() 
  out[[1]]<-zhisout
  out[[2]]<-rhisout
  return(out)
  }


#get density for fibre
getrho_fibre <- function(frame) {
  tempdcd<-mydcd[frame,] #get coor for current frame 
  tempcell<-mydcdcell[frame,] #get box size
  tempdcd<-pbcwrap_fibre(tempdcd,tempcell) #wrap coordinates
  fibresel<-atom.select(mypdb,resid=c("CSH","CSSC")) #select heavy atoms only
  particle<-matrix(tempdcd[fibresel$xyz],ncol = 3,byrow = TRUE) #convert coord to a 3 column matrix (x, y, z) 
  c(frame,length(fibresel$atom)/convhulln(particle,"FA")$vol) #calculate heavy atom density 
}

#get dcd name list
getdcdlist <- function(firstrun,midfix="") { # Function to get list of DCD files
  dcdname_list <- character(0)
  for (runindex in firstrun:lastrun) {   
    dcdname <- sprintf("%s%02d.%sdcd", dcdprefix, runindex, midfix)#construct file name
    dcdfile <- file.path(dcddir, dcdname)  # Construct file name with dir 
    if (file.exists(dcdfile)) {     # Check if file exists
      dcdname_list <- c(dcdname_list, dcdname)       # Add to the list
    } else {
      break       # Break if file doesn't exist
    }
  }
  return(dcdname_list)   # Return DCD name list without dcddir in front
}

#get cos of angle between two vector
getcostheta <- function(ref, comp) {
  sum(ref * comp) / (sqrt(sum(ref^2)) * sqrt(sum(comp^2)))  # Calculate dot product and normalize the vectors to obtain the cosine of the angle
}

#get resid depending on z and r cut off provided
getFiberSectionResid<-function(tempdcd,frame,FibreCoorPercentCutZ,PURE,MoietyCut=0.75,TIPS=FALSE) {
  out<-list()
  fibresel<-atom.select(mypdb,"noh",resid=UniqResnameList) #heavy atoms only
  FibreCoor<-tempdcd[fibresel$xyz] #get coor
  FibreCoorX<-FibreCoor[seq(1,length(FibreCoor)-2,3)] #x
  FibreCoorY<-FibreCoor[seq(2,length(FibreCoor)-1,3)] #y
  FibreCoorZ<-FibreCoor[seq(3,length(FibreCoor)-0,3)] #z
  FibreCoorR<-sqrt(FibreCoorX^2+FibreCoorY^2) #radial distance
  FibreCoorHighZ<-quantile(FibreCoorZ, FibreCoorPercentCutZ) #find upper limit in z
  FibreCoorLowZ<-quantile(FibreCoorZ, 1-FibreCoorPercentCutZ) #find lower limit in z
  if (TIPS) {
    FibreIndexZ<-which(FibreCoorZ>FibreCoorHighZ & FibreCoorZ<FibreCoorLowZ) #get index outside z range
    FibreIndexR<-which(FibreCoorR>FibreRangeR) #get index outside r range
  } else {
    FibreIndexZ<-which(FibreCoorZ<FibreCoorHighZ & FibreCoorZ>FibreCoorLowZ) #get index within z range
    FibreIndexR<-which(FibreCoorR<FibreRangeR) #get index within r range
  }  
  FibreIndexZ<-fibresel$atom[FibreIndexZ] #get atomic index corresponding to pdb 
  FibreIndexR<-fibresel$atom[FibreIndexR] #get atomic index corresponding to pdb 
  FibreIndex <- intersect(FibreIndexZ, FibreIndexR) #get index for the inner fiber
  FibreResidCSSC<-mypdb$atom[which(mypdb$atom$eleno %in% FibreIndex & mypdb$atom$resid=="CSSC"),]$resno #get resid for CSSC
  FibreResidCSSC <- table(FibreResidCSSC) #count occurance
  FibreResidCSSC <- as.numeric(names(FibreResidCSSC[FibreResidCSSC > 15 * MoietyCut])) #only certain percent counts
  if (PURE==TRUE) {
    FibreResidCSH<-mypdb$atom[which(mypdb$atom$eleno %in% FibreIndex & mypdb$atom$resid=="CSH"),]$resno #get resid for CSH
    FibreResidCSH <- table(FibreResidCSH) #count occurance
    FibreResidCSH <- as.numeric(names(FibreResidCSH[FibreResidCSH > 15 * MoietyCut])) #only certain percent counts
    out[["CSH"]][["Inside"]]<-FibreResidCSH
    out[["CSH"]][["Outside"]]<-setdiff(seq(1, maxCSH, 1), FibreResidCSH)
  }
  out[["CSSC"]][["Inside"]]<-FibreResidCSSC
  out[["CSSC"]][["Outside"]]<-setdiff(seq(1, maxCSSC, 1), FibreResidCSSC)
  out[["tempdcd"]]<-tempdcd
  out
}

#get orientation of CSH equvalents relative to the fibre axis
getOrient<-function(frame) {
  out<-list()
  tempdcd<-mydcd[frame,] #get coor for current frame 
  tempcell<-mydcdcell[frame,] #get box size
  tempdcd<-pbcwrap_fibre(tempdcd,tempcell) #pbd correction
  particle<-matrix(tempdcd,ncol = 3,byrow = TRUE) #convert coord to a 3 column matrix (x, y, z) 
  tempout<-do.call(rbind,sapply(1:2,function(i) {
    do.call(rbind,sapply(seq(1,length(NameIndexList[[i*2]]),1),function(Index) {
      Coor1<-particle[NameIndexList[[2*i]][Index],]
      Coor2<-particle[NameIndexList[[2*i-1]][Index],]
      CSHvec<-Coor2-Coor1
      CSHcom<-(Coor2+Coor1)/2
      CSHcomR<-sqrt(CSHcom[1]^2+CSHcom[2]^2)
      if (CSHcomR<28) {
        c(getcostheta(c(0,0,1),CSHvec),CSHcomR,CSHcom[3])
      }
    }, simplify = FALSE))
  }, simplify = FALSE))
  hist_r <- cut(tempout[,2], breaks = seq(0, 30, by = 1), include.lowest = TRUE)
  hist_z <- cut(tempout[,3], breaks = seq(-160, 160, by = 1), include.lowest = TRUE)
  hist_ang <- cut(tempout[,1], breaks = seq(-1, 1.1, by = 0.1), include.lowest = TRUE)
  abs_hist_ang <- cut(abs(tempout[,1]), breaks = seq(0, 1.1, by = 0.1), include.lowest = TRUE)
  out<-list()
  out[["hist_ang_r"]] <- as.matrix(table(hist_r, hist_ang)) #cut 2 column out data into 2d histogram, range from 0 to 30 for both column.
  out[["hist_ang_z"]] <- as.matrix(table(hist_z, hist_ang)) #cut 2 column out data into 2d histogram, range from 0 to 30 for both column.
  out[["abs_hist_ang_r"]] <- as.matrix(table(hist_r, abs_hist_ang)) #cut 2 column out data into 2d histogram, range from 0 to 30 for both column.
  out[["abs_hist_ang_z"]] <- as.matrix(table(hist_z, abs_hist_ang)) #cut 2 column out data into 2d histogram, range from 0 to 30 for both column.
  return(out)
}

getOrientOld<-function(frame,PURE=FALSE) {
  #Resids<-getFiberSectionResid(tempdcd,frame,FibreCoorPercentCutZ,PURE)
  
  for (Position in c("Inside","Outside")) {
    for (Resname in UniqResnameList) {
      # for (resid in Resids[[Resname]][[Position]][which(Resids[[Resname]][[Position]]%%2!=0)]) {
      #   OrientSelCom<-list()
      #   for (i in 1:2) {
      #     OrientSel<-atom.select(mypdb,resid=gsub("[^[:alpha:]]", "", Resname),resno=resid,elety=ifelse(resid%%2==0,OrientSelText[[Resname]][[2]][i],OrientSelText[[Resname]][[1]][i]))
      #     print(OrientSel$atom)
      #     OrientSelCom[[i]]<-Resids$tempdcd[OrientSel$xyz]
      #   }
      #   CosOrient<-getcostheta(c(0,0,1),OrientSelCom[[2]]-OrientSelCom[[1]])
      #   out[[gsub("[^[:alpha:]]", "", Resname)]][[Position]]<-c(out[[gsub("[^[:alpha:]]", "", Resname)]][[Position]],CosOrient)
      # }
      if (Resname == "CSSC") {
        Moiety1<-sapply(Resids[[Resname]][[Position]][which(Resids[[Resname]][[Position]]%%2!=0)], function(resid) {
          zindex1<-(resid*26-12)*3
          zindex2<-(resid*26-22)*3
          getcostheta(c(0,0,1),Resids$tempdcd[seq(zindex1-2,zindex1,1)]-Resids$tempdcd[seq(zindex2-2,zindex2,1)])
        })
        Moiety2<-sapply(Resids[[Resname]][[Position]][which(Resids[[Resname]][[Position]]%%2==0)], function(resid) {
          zindex1<-(resid*26-16)*3
          zindex2<-(resid*26-14)*3
          getcostheta(c(0,0,1),Resids$tempdcd[seq(zindex1-2,zindex1,1)]-Resids$tempdcd[seq(zindex2-2,zindex2,1)])
        }
        )
        out[[gsub("[^[:alpha:]]", "", Resname)]][[Position]]<-c(Moiety1,Moiety2)
        
      } else {
        out[[gsub("[^[:alpha:]]", "", Resname)]][[Position]]<-sapply(Resids[[Resname]][[Position]], function(resid) {
          BaseInd<-resid*27+26*maxCSSC
          zindex1<-(BaseInd-12)*3
          zindex2<-(BaseInd-23)*3
          getcostheta(c(0,0,1),Resids$tempdcd[seq(zindex1-2,zindex1,1)]-Resids$tempdcd[seq(zindex2-2,zindex2,1)])
        })
      }
    }
    for (Resname in unique(gsub("[^[:alpha:]]", "",UniqResnameList))) {
      if (abs==FALSE) {
        out[[Resname]][[Position]]<-table(cut(out[[Resname]][[Position]],breaks = c(-Inf,seq(0.1,1,0.9),Inf)))
      } else {
        out[[Resname]][[Position]]<-table(cut(abs(out[[Resname]][[Position]]),breaks = c(-Inf,seq(-0.9,1,0.9),Inf),include.lowest = TRUE))
      }
    }
  }
  out
}

#get orientation of CSH equvalents relative to the fibre axis for CSSC only
getFiberSectionResidCSSC<-function(tempdcd,frame,FibreCoorPercentCutZ,MoietyCut=0.75,TIPS=FALSE) {
  out<-list()
  fibresel<-atom.select(mypdb,"noh",resid=UniqResnameList) #heavy atoms only
  particle<-matrix(tempdcd[fibresel$xyz],ncol = 3,byrow = TRUE) #convert coord to a 3 column matrix (x, y, z) 
  FibreCoorR<-sqrt(particle[,1]^2+particle[,2]^2) #radial distance
  FibreCoorHighZ<-quantile(particle[,3], FibreCoorPercentCutZ) #find upper limit in z
  FibreCoorLowZ<-quantile(particle[,3], 1-FibreCoorPercentCutZ) #find lower limit in z
  if (TIPS) {
    FibreIndexZ<-which(particle[,3]>FibreCoorHighZ & particle[,3]<FibreCoorLowZ) #get index outside z range
    FibreIndexR<-which(FibreCoorR>FibreRangeR) #get index outside r range
  } else {
    FibreIndexZ<-which(particle[,3]<FibreCoorHighZ & particle[,3]>FibreCoorLowZ) #get index within z range
    FibreIndexR<-which(FibreCoorR<FibreRangeR) #get index within r range
  }  
  FibreIndexZ<-fibresel$atom[FibreIndexZ] #get atomic index corresponding to pdb 
  FibreIndexR<-fibresel$atom[FibreIndexR] #get atomic index corresponding to pdb 
  FibreIndex <- intersect(FibreIndexZ, FibreIndexR) #get index for the inner fiber
  FibreResidCSSC<-mypdb$atom[which(mypdb$atom$eleno %in% FibreIndex & mypdb$atom$resid=="CSSC"),]$resno #get resid for CSSC
  FibreResidCSSC <- table(FibreResidCSSC) #count occurance
  FibreResidCSSC <- as.numeric(names(FibreResidCSSC[FibreResidCSSC > 30 * MoietyCut])) #only certain percent counts
  out[["CSSC"]][["Inside"]]<-FibreResidCSSC
  out[["CSSC"]][["Outside"]]<-setdiff(seq(1, maxCSSC, 1), FibreResidCSSC)
  out[["tempdcd"]]<-tempdcd
  out
}

#get PHE PHE distance depending on output of getFiberSectionResidCSSC
getConfor<-function(frame) {
  out<-list()
  tempdcd<-mydcd[frame,] #get coor for current frame 
  tempcell<-mydcdcell[frame,] #get box size
  tempdcd<-pbcwrap_fibre(tempdcd,tempcell) #pbd correction
  particle<-matrix(tempdcd,ncol = 3,byrow = TRUE) #convert coord to a 3 column matrix (x, y, z) 
  Resids<-getFiberSectionResidCSSC(tempdcd,frame,FibreCoorPercentCutZ)
  for (Position in c("Inside","Outside")) {
    out[["CSSC"]][[Position]]<-sapply(Resids[["CSSC"]][[Position]], function(resid) {
      euclidean(particle[PHEClist[resid,1],],particle[PHEClist[resid,2],])
    })
    out[["CSSC"]][[Position]]<-table(cut(out[["CSSC"]][[Position]],breaks = c(-Inf,PHEdisBreak,Inf)))
  }
  out
}




getnorm <- function(v1, v2, v3) {
  va <- v2 - v1   # Calculate difference vector between v1 and v2
  vb <- v3 - v1   # Calculate difference vector between v1 and v3
  vc <- cross_manual(va, vb)  # Calculate cross product of va and vb
  t(apply(vc, 1, function(row) {
    row / sqrt(sum(row^2))
  }))
  #vc / sqrt(sum(vc^2))  # Normalize vc and return it
}



#get the COM of particular moiety
getCOM <- function(startresi,finalresi,tempdcd,moiety,resname="CSSC",maxcssc=maxcssc) {
  return(##### create an empty com matrix
    do.call(rbind,
            sapply(seq(startresi,finalresi,1),
                   function(i) {
                     TempMoietyXYZ<-moietyxyz[[moiety]][[1]][[i]]
                     #####   get coordinate
                     refmoietyxyz<-tempdcd[TempMoietyXYZ]
                     #####   calculate COM of res i
                     com1<-as.vector(com.xyz(refmoietyxyz))
                     if (resname=="CSSC") {
                       TempMoietyXYZ<-moietyxyz[[moiety]][[2]][[i]]
                       #####   get coordinate
                       refmoietyxyz<-tempdcd[TempMoietyXYZ]
                       #####   calculate COM of res i
                       com2<-as.vector(com.xyz(refmoietyxyz))
                       out<-rbind(com1,com2)
                       out
                     } else {
                     }
                   }
                   ,simplify = FALSE)
    )
  )
}




getBBvec <- function(startresi,finalresi,tempdcd,maxcssc=maxcssc) {
  do.call(rbind,
          sapply(seq(startresi,finalresi,1),
                 function(i) {
                   C1COM<-tempdcd[BBxyz$C1[[i]]]
                   C20COM<-tempdcd[BBxyz$C20[[i]]]
                   Vref<-C20COM-C1COM
                   rbind(Vref,Vref*-1)
                 },simplify = FALSE)
  )
}

getBBCOM <- function(tempdcd) {
  out<-do.call(rbind,
          sapply(1:maxcssc,
                 function(j) {
                   C1COM<-tempdcd[BBxyz$C1[[j]]]
                   C20COM<-tempdcd[BBxyz$C20[[j]]]
                   BBCOM<-(C20COM+C1COM)/2
                   BBCOM
                   },simplify = FALSE)
  )
  return(out)
}

getBBangHist<-function(framenum) {
  #wrap cell 
  tempdcd<-mydcd[framenum,] #get coor for current frame 
  tempcell<-mydcdcell[framenum,] #get box size
  tempdcd<-pbcwrap_fibre(tempdcd,tempcell)
  #get COM of each  moiety
  COMlistPHE<-{
    getCOM(mincssc,maxcssc,tempdcd,"phe") 
  }  
  COMlistAM<-{
    getCOM(mincssc,maxcssc,tempdcd,"am") 
  }  
  COMlistTHI<-{
    getCOM(mincssc,maxcssc,tempdcd,"thi") 
  }  
  #get backbone vector
  BBvec<-getBBvec(mincssc,maxcssc,tempdcd)
  #loop through CSSC, and caluclate angle based on COM based contact
  out<-do.call(rbind,sapply(1:maxCSSC,
                            function (i) {
                              #ignore self interaction
                              AvoidList<-c(i*2-1,i*2)
                              #create an empty vector
                              TempOut<-c()
                              #loop thorugh all the moiety
                              for (moietycounter1 in 1:length(moietylist)) {
                                #avoid double counting
                                for (moietycounter2 in moietycounter1:length(moietylist)) {
                                  #get moiety name 
                                  COM1<-moietylist[moietycounter1]
                                  COM2<-moietylist[moietycounter2]
                                  #get moiety COM
                                  COMlist1<-get(paste0("COMlist",COM1))
                                  COMlist2<-get(paste0("COMlist",COM2))
                                  #get C1 & C20 COM if res i
                                  C1COMi<-tempdcd[BBxyz$C1[[i]]]
                                  C20COMi<-tempdcd[BBxyz$C20[[i]]]
                                  #####get cutoff
                                  contactcutoff<-getcutoff("CSSC","CSSC",COM1,COM2)
                                  #get corresponding CSX index
                                  for (IndexI in c(2*i-1,2*i)) {
                                    #####   get COM 
                                    comi<-COMlist1[IndexI,]
                                    #make it faster by screening out moieties too far away
                                    CloseIndex <- which(apply(COMlist2, 1, function(row) all(abs(row - comi) <= contactcutoff))) #screen out atoms that are too far
                                    #avoid double count by only looking for bigger resid
                                    CloseIndex<-CloseIndex[which(CloseIndex>i*2)]
                                    #also compare to avoidlist to avoid double counting
                                    CloseIndex<-CloseIndex[which(CloseIndex%in%AvoidList==FALSE)]
                                    #loop thorugh moieties that are close enough
                                    for (ComIndex in CloseIndex) {
                                      #double check if it should be avoided
                                      if (ComIndex%in%AvoidList) {
                                        next
                                      }
                                      comComp<-COMlist1[ComIndex,]
                                      
                                      #there's contact if they are within a cutoff
                                      if (euclidean(comComp,comi)<=contactcutoff) {
                                        #get coordinates of C1 and C20
                                        C1COMcomp<-tempdcd[BBxyz$C1[[ceiling(ComIndex/2)]]]
                                        C20COMcomp<-tempdcd[BBxyz$C20[[ceiling(ComIndex/2)]]]
                                        #setup a comparison matrix
                                        BBmatrix<-outer(c("C1COMi","C20COMi"),c("C1COMcomp","C20COMcomp"), function(x,y) {
                                          paste(x,y)
                                        }
                                        )
                                        #setup a corresponding matrix for BBvec selection
                                        VecSelIndexMatrix<-outer(c(i*2-1,i*2),c(ceiling(ComIndex/2)*2-1,ceiling(ComIndex/2)*2), function(x,y) {
                                          paste(x,y)
                                        }
                                        )
                                        #calculate backbone end differences in a matrix
                                        DisMatrix<-mapply(function(i) {
                                          split_entry <- strsplit(i, " ")[[1]]
                                          x = split_entry[1]
                                          y = split_entry[2]
                                          euclidean(get(x),get(y))
                                        }
                                        ,BBmatrix)             
                                        #select closest ends
                                        VecIndex<-as.numeric(strsplit(VecSelIndexMatrix[which(DisMatrix==min(DisMatrix))], " ")[[1]])
                                        #get average position of closest contact
                                        CloseBBcom<-strsplit(BBmatrix[which(DisMatrix==min(DisMatrix))], " ")[[1]]
                                        #average xyz of closest contact
                                        BBcontCoor<-(get(CloseBBcom[1])+get(CloseBBcom[2]))/2
                                        BBcontCoor<-c(sqrt(BBcontCoor[1]^2+sqrt(BBcontCoor[2]^2)),BBcontCoor[3])
                                        #compute angles
                                        BBang<-getcostheta(BBvec[VecIndex[1],],BBvec[VecIndex[2],])
                                        #record angle 
                                        #TempOut<-c(TempOut,getcostheta(Vcomp,Vref))
                                        TempOut<-rbind(TempOut,c(BBang,BBcontCoor))
                                        #for debugging
                                        #,i,ceiling(ComIndex/2),COM1,COM2,IndexI,ComIndex))
                                        #stop()
                                        #avoid comparison between the same CSSC pair
                                        AvoidList<-c(AvoidList,c(ceiling(ComIndex/2)*2,ceiling(ComIndex/2)*2-1))
                                      }
                                    }
                                  }
                                }
                              }
                              #return output
                              TempOut
                            },simplify = FALSE
  )
  )
  hist_r <- cut(out[,2], breaks = seq(0, 30, by = 1), include.lowest = TRUE)
  hist_z <- cut(out[,3], breaks = seq(-160, 160, by = 1), include.lowest = TRUE)
  hist_ang_BB <- cut(out[,1], breaks = seq(-1, 1+Bin, by = Bin), include.lowest = TRUE)
  abs_hist_ang_BB <- cut(abs(out[,1]), breaks = seq(1, 1+Bin, by = Bin), include.lowest = TRUE)
  out<-list()
  out[["hist_ang_r"]] <- as.matrix(table(hist_r, hist_ang_BB)) #cut 2 column out data into 2d histogram, range from 0 to 30 for both column.
  out[["hist_ang_z"]] <- as.matrix(table(hist_z, hist_ang_BB)) #cut 2 column out data into 2d histogram, range from 0 to 30 for both column.
  out[["abs_hist_ang_r"]] <- as.matrix(table(hist_r, abs_hist_ang_BB)) #cut 2 column out data into 2d histogram, range from 0 to 30 for both column.
  out[["abs_hist_ang_z"]] <- as.matrix(table(hist_z, abs_hist_ang_BB)) #cut 2 column out data into 2d histogram, range from 0 to 30 for both column.
  out[["abs_hist_ang"]]<-table( abs_hist_ang_BB)
  out[["hist_ang"]]<-table(hist_ang_BB)
  
  out  
}

getBBaxialAngHist<-function(framenum) {
  # #wrap cell 
  tempdcd<-mydcd[framenum,] #get coor for current frame 
  tempcell<-mydcdcell[framenum,] #get box size
  tempdcd<-pbcwrap_fibre(tempdcd,tempcell)
  
  #get backbone vector
  BBvec<-getBBvec(mincssc,maxcssc,tempdcd)
  BBcom<-getBBCOM(tempdcd)
  
  #loop through CSSC, and caluclate angle based on COM based contact
  tempout<-do.call(rbind,
               sapply(1:maxCSSC,
                            function (i) {
                              tempBBvec<-BBvec[2*i,]
                              
                              tempBBcom<-BBcom[i,]
                              tempBBcomz<-tempBBcom[3]
                              tempBBcomr<-sqrt(tempBBcom[1]^2+tempBBcom[2]^2)
                              tempBBcom[3]<-0
                              c(getcostheta(tempBBvec,c(0,0,1)),tempBBcomr,tempBBcomz)
                            },simplify = FALSE
               )
  )
  hist_r <- cut(tempout[,2], breaks = seq(0, 30, by = 1), include.lowest = TRUE)
  hist_z <- cut(tempout[,3], breaks = seq(-160, 160, by = 1), include.lowest = TRUE)
  hist_ang_BB <- cut(tempout[,1], breaks = seq(-1, 1+Bin, by = Bin), include.lowest = TRUE)
  abs_hist_ang_BB <- cut(abs(tempout[,1]), breaks = seq(0, 1+Bin, by = Bin), include.lowest = TRUE)
  OrdPar<-(3*tempout[,1]^2-1)/2
  # Create a new column for the bins based on the third column
  tempbin <- cut(tempout[, 3], breaks = seq(-160, 160, by = 1), include.lowest = TRUE, labels = FALSE)
  
  # Calculate the average of the first column within each bin
  OrdPar <- aggregate(OrdPar, by = list(tempbin), FUN = mean)
  

  out<-list()
  out[["hist_ang_r"]] <- as.matrix(table(hist_r, hist_ang_BB)) #cut 2 column out data into 2d histogram, range from 0 to 30 for both column.
  out[["hist_ang_z"]] <- as.matrix(table(hist_z, hist_ang_BB)) #cut 2 column out data into 2d histogram, range from 0 to 30 for both column.
  out[["abs_hist_ang_r"]] <- as.matrix(table(hist_r, abs_hist_ang_BB)) #cut 2 column out data into 2d histogram, range from 0 to 30 for both column.
  out[["abs_hist_ang_z"]] <- as.matrix(table(hist_z, abs_hist_ang_BB)) #cut 2 column out data into 2d histogram, range from 0 to 30 for both column.
  out[["abs_hist_ang"]]<-table( abs_hist_ang_BB)
  out[["hist_ang"]]<-table(hist_ang_BB)
  out[["OrdPar"]]<-OrdPar[,2]
  out  
}
if (0) {function(framenum) {
  # #wrap cell 
  tempdcd<-mydcd[framenum,] #get coor for current frame 
   tempcell<-mydcdcell[framenum,] #get box size
   tempdcd<-pbcwrap_fibre(tempdcd,tempcell)
  
  #get backbone vector
  BBvec<-getBBvec(mincssc,maxcssc,tempdcd)
  #loop through CSSC, and caluclate angle based on COM based contact
  do.call(rbind,sapply(1:maxCSSC,
                       function (i) {
                         
                         Vref<-BBvec[2*i,]
                         
                         getcostheta(c(0,0,1),Vref)
                         
                         
                       },simplify = FALSE
  )
  )
}
}
getBBradialAngHist<-function(framenum) {
  # #wrap cell 
  tempdcd<-mydcd[framenum,] #get coor for current frame 
   tempcell<-mydcdcell[framenum,] #get box size
   tempdcd<-pbcwrap_fibre(tempdcd,tempcell)
  
  #get backbone vector
  BBvec<-getBBvec(mincssc,maxcssc,tempdcd)
  BBcom<-getBBCOM(tempdcd)

  #loop through CSSC, and caluclate angle based on COM based contact
  out<-do.call(rbind,sapply(1:maxCSSC,
                            function (i) {
                              tempBBvec<-BBvec[2*i,]
                              tempBBvec[3]<-0
                              tempBBcom<-BBcom[i,]
                              tempBBcomz<-tempBBcom[3]
                              tempBBcomr<-sqrt(tempBBcom[1]^2+tempBBcom[2]^2)
                              tempBBcom[3]<-0
                              c(getcostheta(tempBBvec,tempBBcom),tempBBcomr,tempBBcomz)
                            },simplify = FALSE
  )
  )
  hist_r <- cut(out[,2], breaks = seq(0, 30, by = 1), include.lowest = TRUE)
  hist_z <- cut(out[,3], breaks = seq(-160, 160, by = 1), include.lowest = TRUE)
  hist_ang_BB <- cut(out[,1], breaks = seq(-1, 1+Bin, by = Bin), include.lowest = TRUE)
  abs_hist_ang_BB <- cut(abs(out[,1]), breaks = seq(0, 1+Bin, by = Bin), include.lowest = TRUE)
  out<-list()
  out[["hist_ang_r"]] <- as.matrix(table(hist_r, hist_ang_BB)) #cut 2 column out data into 2d histogram, range from 0 to 30 for both column.
  out[["hist_ang_z"]] <- as.matrix(table(hist_z, hist_ang_BB)) #cut 2 column out data into 2d histogram, range from 0 to 30 for both column.
  out[["abs_hist_ang_r"]] <- as.matrix(table(hist_r, abs_hist_ang_BB)) #cut 2 column out data into 2d histogram, range from 0 to 30 for both column.
  out[["abs_hist_ang_z"]] <- as.matrix(table(hist_z, abs_hist_ang_BB)) #cut 2 column out data into 2d histogram, range from 0 to 30 for both column.
  out[["abs_hist_ang"]]<-table( abs_hist_ang_BB)
  out[["hist_ang"]]<-table(hist_ang_BB)
  
  out  
}


getcutoff<- function(pair1,pair2,COM1,COM2) {
  if (pair1==pair2) {
    cutofflist[[paste0(pair1,"-",pair2)]][COM1,COM2]
  } else {
    cutofflist[["CSH-CSSC"]][COM1,COM2]
  }
}

getCOM_perResID <- function(framenum,wrapping=TRUE) {
  out<-list()
  #wrap cell 
  tempdcd<-mydcd[framenum,] #get coor for current frame 
  tempcell<-mydcdcell[framenum,] #get box size
  if (wrapping) {
    tempdcd<-pbcwrap_fibre(tempdcd,tempcell)
  }
  for (resname in UniqResnameList)   {
    out[[resname]]<-do.call(rbind,
                            sapply(seq(1,length(moietyxyz[[resname]]),1),
                                   function(i) {
                                     TempMoietyXYZ<-moietyxyz[[resname]][[i]]
                                     #####   get coordinate
                                     refmoietyxyz<-tempdcd[TempMoietyXYZ]
                                     #####   calculate COM of res i
                                     as.vector(com.xyz(refmoietyxyz))
                                   }
                                   ,simplify = FALSE)
    )
  }
  return(out)
}

getMSD<-function(ntau,resname) {
  XcoorTemp<-list()
  YcoorTemp<-list()
  ZcoorTemp<-list()
  for (resid in 1:get(paste0("max",resname))) {
    YcoorTemp[[resid]]<-0
    XcoorTemp[[resid]]<-0
    ZcoorTemp[[resid]]<-0
    for (frame in (ntau-1)*tau+2:(ntau*tau)) {
      if (frame > length(COM)) {
        break
      }
      XcoorTemp[[resid]]<-XcoorTemp[[resid]]+(COM[[frame]][[resname]][resid,1]-COM[[(ntau-1)*tau+1]][[resname]][resid,1])^2
      YcoorTemp[[resid]]<-YcoorTemp[[resid]]+(COM[[frame]][[resname]][resid,2]-COM[[(ntau-1)*tau+1]][[resname]][resid,2])^2
      ZcoorTemp[[resid]]<-ZcoorTemp[[resid]]+(COM[[frame]][[resname]][resid,3]-COM[[(ntau-1)*tau+1]][[resname]][resid,3])^2
    }
    XcoorTemp[[resid]]<-XcoorTemp[[resid]]/(tau)
    YcoorTemp[[resid]]<-YcoorTemp[[resid]]/(tau)
    ZcoorTemp[[resid]]<-ZcoorTemp[[resid]]/(tau)
  }
  XcoorTemp<-mean(unlist(XcoorTemp))
  YcoorTemp<-mean(unlist(YcoorTemp))
  ZcoorTemp<-mean(unlist(ZcoorTemp))
  out<-list()
  out[[1]]<-XcoorTemp
  out[[2]]<-YcoorTemp
  out[[3]]<-ZcoorTemp
  out
}


#basic moiety break down
moietyname <- c("N1", "C1", "O1", "C2", "C3", "S", "N2", "C4", "O2", "C5", "C6", "C7", "C8", "C9", "C10","S1")
moietyname2 <- c("S2", "C11", "C12", "C13", "C14", "N3", "C15", "C16", "C17", "C18", "C19", "C20", "O3", "O4", "N4")
phemoietyname <- c("C5", "C6", "C7", "C8", "C9", "C10")
ammoietyname <- c("N1", "C1", "O1", "C3", "S", "N2", "C4", "O2", "S1")
thimoietyname <- c("C3", "S","S1")
phemoietyname2 <- c("C12", "C13", "C14", "C16","C18","C19")
ammoietyname2 <- c("S2", "C11", "N3", "C17", "C20", "O3", "O4", "N4")
thimoietyname2 <- c("C11", "S2")

# Define the moiety name vectors
phemoietyname <- c("C5", "C6", "C7", "C8", "C9", "C10")
ammoietyname <- c("N1", "C1", "O1", "N2", "C4", "O2", "C2")
thimoietyname <- c("C3", "S", "S1")
phemoietyname2 <- c("C12", "C13", "C14", "C16", "C18", "C19")
ammoietyname2 <- c("N3", "C17", "C20", "O3", "O4", "N4","C15")
thimoietyname2 <- c("C11", "S2")

# Combine all moiety names into a list
MoietyNameList <- list(
  phemoietyname = phemoietyname,
  ammoietyname = ammoietyname,
  thimoietyname = thimoietyname,
  phemoietyname2 = phemoietyname2,
  ammoietyname2 = ammoietyname2,
  thimoietyname2 = thimoietyname2
)

convert_COMlist_to_df <- function(com_list) {
  temp_data <- data.frame(x = numeric(), y = numeric(), z = numeric(), resid = integer(), moietyind = integer())
  
  for (resid in 1:length(com_list)) {
    for (moietyind in 1:length(com_list[[resid]])) {
      # Extract the vector
      vector <- com_list[[resid]][[moietyind]]
      # Append the vector and identifiers to the temp_data
      temp_data <- rbind(temp_data, data.frame(
        x = vector[1],
        y = vector[2],
        z = vector[3],
        resid = as.numeric(resid),
        moietyind = as.numeric(moietyind)
      ))
      
    }
  }
  return(temp_data)
}



