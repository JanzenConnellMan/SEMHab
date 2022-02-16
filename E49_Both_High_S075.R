filepath="~/Validation/SEMHab"
setwd(filepath)


############## create Spatial Grid of Environment Types ####################
dm <- 499

set.seed(100)
Habs<-read.csv("Hab_Med.csv")
Habs2 <- unlist(Habs)
Habs3 <- unname(Habs2)
Habs3 <- sample(Habs3)
Mat.S.H <- matrix(Habs3,nrow=dm,ncol=dm)


################################## Species Fitness + Habitat Type ###################################

#### function for fitness ####

Sigma <- .075   # niche width 

Habitat_Survival <- function(x,y){ # Gaussian Function of Fitness + habtitat Type 
  Hab.Pr <- exp(  -(x-y)^2/(2*Sigma^2)    )
  return(Hab.Pr)
}


###### assigning Habitat Values of species #### 
S <- 300

ratio <- 1/S
H_Opt <- seq(ratio-(ratio/2),ratio*S,ratio)



mat.torus <- function(Matrix,Total.Neighbors,xcord,ycord){      # Torus of Space
  
  if(Total.Neighbors==0){return(Matrix[xcord,ycord])}
  
  dm    <- nrow(Matrix)                               # dimension of matrix (n x n marix, so rows = col)
  
  Crown.Pot       <-  c(seq(1,10,by=2)^2)-1    # arbitrarily set to 10, because that many neighbors gets silly biologically, just used for computation 
  Crowns          <-  seq(0,length(Crown.Pot))
  Crown           <-  Crowns[min(which(Crown.Pot >= Total.Neighbors))]  #returns crown exension total (even if not whole crown)
  Sample          <-  Total.Neighbors - (Crown.Pot[Crown])   # Returns how many to sample from the extra plots in addtion to crown
  
  rows  <- c(  (xcord-Crown):(xcord+Crown)       )    # figure out which rows of matrix crown falls in
  cols  <- c(  (ycord-Crown):(ycord+Crown)       )    # figure out which columns of matrix crown falls in
  
  rows[which(rows<1)]    <- rows[which(rows<1)] + dm  # if crown extends to a row value less than 1, go to opposite side of torus 
  rows[which(rows>dm)]   <- rows[which(rows>dm)] - dm # if crown extends to a row value greater than dm, go to opposite side of torus
  
  cols[which(cols<1)]    <- cols[which(cols<1)] + dm  # if crown extends to a column value less than 1, go to opposite side of torus 
  cols[which(cols>dm)]   <- cols[which(cols>dm)] - dm # if crown extends to a column value greater than dm, go to opposite side of torus
  
  JC_Matrix              <- Matrix[rows,cols ]        # returns subset of matrix / trees that are in JCE zone + extras 
  JC_Matrixb             <- JC_Matrix                # Copy matrix 
  
  JC_Matrix[c(1,nrow(JC_Matrix)),]    <- 0             
  JC_Matrix[,c(1,ncol(JC_Matrix))]    <- 0 
  
  vec.JC_Matrix <- c(JC_Matrix)
  vec.JC_Matrix <- vec.JC_Matrix[vec.JC_Matrix>0]
  
  
  JC_Matrixb[2:(nrow(JC_Matrix)-1),2:(ncol(JC_Matrix)-1)]    <- 0 
  
  vec.JC_Matrixb <- c(JC_Matrixb)
  vec.JC_Matrixb <- JC_Matrixb[JC_Matrixb>0]
  vec.JC_Matrixb <- sample(vec.JC_Matrixb,Sample, replace = FALSE)
  
  return(unique(c(vec.JC_Matrix,vec.JC_Matrixb )  ))
  
}
'%!in%' <- function(x,y)!('%in%'(x,y))
mat.torus <- function(Matrix,Rad,xcord,ycord){      # Torus of Space
  
  dm    <- nrow(Matrix)                               # dimension of matrix (n x n marix, so rows = col)
  
  Crown.Pot       <-  c(seq(1,Rad,by=2)^2)-1    # arbitrarily set to 10, because that many neighbors gets silly biologically, just used for computation 
  Crowns          <-  seq(0,length(Crown.Pot))
  Crown           <-  Crowns[min(which(Crown.Pot >= (Rad-1)^2))]  #returns crown exension total (even if not whole crown)
  
  rows  <- c(  (xcord-Crown):(xcord+Crown)       )    # figure out which rows of matrix crown falls in
  cols  <- c(  (ycord-Crown):(ycord+Crown)       )    # figure out which columns of matrix crown falls in
  
  rows[which(rows<1)]    <- rows[which(rows<1)] + dm  # if crown extends to a row value less than 1, go to opposite side of torus 
  rows[which(rows>dm)]   <- rows[which(rows>dm)] - dm # if crown extends to a row value greater than dm, go to opposite side of torus
  
  cols[which(cols<1)]    <- cols[which(cols<1)] + dm  # if crown extends to a column value less than 1, go to opposite side of torus 
  cols[which(cols>dm)]   <- cols[which(cols>dm)] - dm # if crown extends to a column value greater than dm, go to opposite side of torus
  
  JC_Matrix              <- Matrix[rows,cols ]        # returns subset of matrix / trees that are in JCE zone + extras 
  
  return(JC_Matrix)
  
}




#dm <- 100
S.list <- seq(1:S)
TimeSteps <- 50000

dd <- sample(1:S, dm*dm, replace = TRUE)
Mat.S <- matrix(dd,nrow=dm,ncol=dm)

df.Props                <- data.frame( matrix(NA,ncol=S+1,nrow=TimeSteps) )
df.Props[,1]            <- seq(1:TimeSteps)
df.Props[1,2:(S+1)]     <- c(table(Mat.S))/(dm*dm)

A <- rep(.5,S)

set.seed(100)
Y <- rlnorm(S,meanlog=0,sd=1)
#Y <- seq(1,100,length=S)

d <- 1
Dist.Rate <- .0025
Rad <- 7

for(mm in 1:TimeSteps){
  
  Mat.S2 <- Mat.S  
  P      <- c(table(factor(Mat.S2, levels = 1:S)))/(dm*dm)   # P goes to proportion of each species in environment
  
  Prob.Dist                             <- matrix(runif(dm*dm),ncol=dm) # Matrix that defines hte probability that each species is disurbed
  Prob.Dist[Prob.Dist >= Dist.Rate]     <- NA # Spcies with draw less than Dist.Rate are distriubed 
  
  df.Rep                      <- which(!is.na(Prob.Dist), arr.ind=TRUE) # saves the indexes of each location that is distrubed
  
  x.val  <- df.Rep[,1]   # x coordinates of disturbance
  y.val  <- df.Rep[,2]   # y coordinates of distriubance
  
  Replaceb <- length(x.val)  # total number of distrubances
  
  
  
  Replacements <- sapply(1:Replaceb, function(x){  # function that determines which speices replaces disturbed patches (apply function loops over all distrubed)
    
    
    Habitat_Type  <- Mat.S.H[x.val[x],y.val[x]]  # Returns Habitat Type 
    
    JC.Victims   <- mat.torus(Mat.S2,Rad,x.val[x],y.val[x])
    JC.Pred       <- c(table(factor(JC.Victims, levels = 1:S))) 
    Predation     <- exp(-A*JC.Pred)
    H_Effect      <- sapply(1:S, function(z) Habitat_Survival(H_Opt[z],  Habitat_Type)  )
    
    Seeds         <-  Y*H_Effect*Predation*P
    Total.Seeds   <- as.numeric(sum(Seeds)) # total scaled number of seeds in local patch
    Vec.Probs     <- c(Seeds)/Total.Seeds  # probability that each species wins lottery         
    
    
    Vec.Probs.Ordred    <- Vec.Probs[order(as.numeric(names(Vec.Probs)))]  # order previous probability values (1 to S)
    Vec.Sum             <- cumsum(Vec.Probs.Ordred) # creates probability intervals to determine which species wins 
    prob.rep <- runif(1) # draw from uniform distribution to determine the winner of the lottery
    
    Replacement <- as.numeric(names(Vec.Sum[min(which(Vec.Sum > prob.rep))])) # store winner of lottery
    return(Replacement) # return winner
  }
  )
  
  Mat.S[df.Rep] <- Replacements # put winner of lottery into correct location
  df.Props[mm+1,2:(S+1)]     <- c(table(factor(Mat.S, levels = 1:S)))/(dm*dm) # Store proportion of each species at each time step
  
} # Code that runs the simulation (notes inside)


df.PropsM <- as.matrix(df.Props)


filepath="~/Validation/Outputs2"
setwd(filepath)


write.csv(df.PropsM,"TS_E49_RAND_BOTH_S075.csv",quote=F,row.names=F)
write.csv(Mat.S,"DIST_E49_RAND_BOTH_S075.csv",quote=F,row.names=F)






