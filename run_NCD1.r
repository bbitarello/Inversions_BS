##########################################################################################
#       Use Phase 1 1000G data to verify allele frequencies within inversion in chr4
#       Bárbara Bitarello
#       Created: 21.04.2017
#       Last modified: 25.04.2017
##########################################################################################

#do: source('run_NCD1.R') #1529.50  seconds! for 4 pops
#preamble
library(pegas);library(dplyr)
library(plyr);library(data.table)
library(parallel);library(lattice)
library(SOAR);Sys.setenv(R_LOCAL_CACHE="inversions")
library(ggplot2);library(splitstackshape)
library(pryr)
library(doMC)
#source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/bedtools_inR.R')
#source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/SFS_script.r')
#source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/mclapply2.R')
#source('NCD_func.R')
#my_cores<-detectCores()/2
registerDoMC(11);
Objects()

##################################################################################

#test:

#victor's code with modifications by me:


NCD1 <- function(X, W = 100000, S = 50000) {
  
  windows_dt <- 
    data.table(POS = seq(X[1, POS], X[nrow(X), POS], S))[
      , POS2 := POS + W][
        -length(POS)] #

    print (paste0('Finished setting up coordinates for chr ', unique(X$CHR)))
  
    setkey(windows_dt, POS, POS2) #
  X[, POS2 := POS] #
    
  X_windows <-
    foverlaps(X, windows_dt, type = "within", nomatch = 0L)[ #this is not ideal but for now it's fine
      , window := .GRP, by = .(POS, POS2)][
        order(window, i.POS)][
          , .(Win.ID = paste(CHR, POS, POS2, sep = "|"), MAF)][,N_Raw:=.N,by = Win.ID]
  
  print (paste0('Finished selecting SNPs per window for chr ', unique(X$CHR)))
  
  X_NCD <-
    X_windows[MAF != 0 & MAF != 1][
      , .(N_Raw= N_Raw,
          N_SNPs_cor = .N,
          NCD1_tf0.5 = sqrt(sum((MAF-0.5)^2)/.N),  #.N is the equivalent of n() in dplyr.
          NCD1_tf0.4 = sqrt(sum((MAF-0.4)^2)/.N),
          NCD1_tf0.3=sqrt(sum((MAF-0.3)^2)/.N)), 
      by = Win.ID]
  
  print (paste0('NCD1 calculations done for chr ', unique(X$CHR)))
  
  unique_windows <-
    
    X_windows[, .(N_SNPs_cor = sum(MAF != 0 & MAF != 1)), by = Win.ID]
  
  setkey(unique_windows, Win.ID, N_SNPs_cor)
  setkey(X_NCD, Win.ID, N_SNPs_cor)
  
  unique(X_NCD[unique_windows])
}

x.1<-vector('list', 26)
for(j in 1:26){
system.time(x.1[[j]]<-foreach(x=1:22, .combine="rbind", .packages=c("data.table")) %dopar%
         NCD1(X=POPS_AF[[j]][[x]], W=3000, S=1500)); # very fast
print(pops[j])}
#test:
mclapply2(x.1, function(X) na.omit(X))-> pops_NCD1_gen
mclapply2(pops_NCD1_gen, function(X) setDT(filter(X, N_SNPs_cor>=10)))-> pops_NCD1_gen_IS;
mclapply2(pops_NCD1_gen_IS, function(X) X[, c('Chr', 'POS1','POS2') := tstrsplit(Win.ID, "|", fixed=TRUE)])-> pops_NCD1_gen_IS;
mclapply2(pops_NCD1_gen_IS, function(X) X[order(as.numeric(Chr), as.numeric(POS1))])-> pops_NCD1_gen_IS
Store(x.1, pops_NCD1_gen, pops_NCD1_gen_IS)
#####
#The End








