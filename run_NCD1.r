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
registerDoMC(22);
Objects()

##################################################################################


x.1<-vector('list', 26)
for(j in 1:26){
system.time(x.1[[j]]<-foreach(x=1:22, .combine="rbind", .packages=c("data.table")) %dopar%
         NCD1(X=POPS_AF[[j]][[x]], W=2000, S=1000)); # very fast
print(pops[j])}
#test:
mclapply2(x.1, function(X) na.omit(X))-> pops_NCD1_gen
mclapply2(pops_NCD1_gen, function(X) setDT(filter(X, N_SNPs_cor>=10)))-> pops_NCD1_gen_IS;
mclapply2(pops_NCD1_gen_IS, function(X) X[, c('Chr', 'POS1','POS2') := tstrsplit(Win.ID, "|", fixed=TRUE)])-> pops_NCD1_gen_IS;
mclapply2(pops_NCD1_gen_IS, function(X) X[order(as.numeric(Chr), as.numeric(POS1))])-> pops_NCD1_gen_IS
Store(x.1, pops_NCD1_gen, pops_NCD1_gen_IS)
#####
#The End








