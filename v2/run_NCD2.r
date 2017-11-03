##########################################################################################
#       Use Phase 1 1000G data to verify allele frequencies within inversion in chr4
#       Bárbara Bitarello
#       Created: 21.04.2017
#       Last modified: 27.07.2017
##########################################################################################

#do: system.time(source('run_NCD2.r')) 808.963 seconds for four populations!
#preamble
library(pegas);library(dplyr)
library(plyr);library(data.table)
library(parallel);library(lattice)
library(SOAR);Sys.setenv(R_LOCAL_CACHE="inversions")
library(ggplot2);library(splitstackshape)
library(pryr)
library(doMC)
source('/mnt/sequencedb/PopGen/barbara/NCD-Statistics/scripts/NCD_func.R')
registerDoMC(11)
Objects()
##################################################################################

x.2<-vector('list', 7)
for(j in 1:7){
system.time(x.2[[j]]<-foreach(x=c(1:22,"X"), .combine="rbind", .packages=c("data.table")) %dopar%
         NCD2(X=POPS_AF[[j]][CHR==x],Y=FD_list[[x]], W=2000, S=1000)); # very fast
print(pops[j])}
 
mclapply2(x.2, function(X) na.omit(X))-> pops_NCD2_gen
mclapply2(pops_NCD2_gen, function(X) X[,IS:=N_FDs_cor+N_SNPs_cor][IS>=8])-> pops_NCD2_gen_IS; #less stringent  filter
mclapply2(pops_NCD2_gen_IS, function(X) X[, c('Chr', 'POS1','POS2') := tstrsplit(Win.ID, "|", fixed=TRUE)])-> pops_NCD2_gen_IS;
mclapply2(pops_NCD2_gen_IS, function(X) X[order(as.numeric(Chr), as.numeric(POS1))])-> pops_NCD2_gen_IS
Store(x.2, pops_NCD2_gen, pops_NCD2_gen_IS)



