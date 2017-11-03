##########################################################################################
#       Bárbara Bitarello
#       Created: 21.04.2017
#       Last modified: 01.11.2017
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
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/mclapply2.R')
source('/mnt/sequencedb/PopGen/barbara/NCD-Statistics/scripts/NCD_func.R')
registerDoMC(11);
Objects()

##################################################################################


x.1<-vector('list', 7)
for(j in 1:7){
	
system.time(x.1[[j]]<-foreach(x=c(1:22,"X"), .combine="rbind", .packages=c("data.table")) %dopar%
         #NCD1(X=POPS_AF[[j]][,CHR==x], W=2000, S=1000)); # very fast
		#POPS_AF[[j]][,NCD1(W=2000, S=1000), by='CHR']); #need to test first
		NCD1(POPS_AF[[j]][CHR==x], W=2000, S=1000))
print(pops[j])}
#test:
mclapply2(x.1, function(X) na.omit(X))-> pops_NCD1_gen
mclapply2(pops_NCD1_gen, function(X) setDT(filter(X, N_SNPs_cor>=8)))-> pops_NCD1_gen_IS; #using less strict filter
mclapply2(pops_NCD1_gen_IS, function(X) X[, c('Chr', 'POS1','POS2') := tstrsplit(Win.ID, "|", fixed=TRUE)])-> pops_NCD1_gen_IS;
mclapply2(pops_NCD1_gen_IS, function(X) X[order(as.numeric(Chr), as.numeric(POS1))])-> pops_NCD1_gen_IS
Store(x.1, pops_NCD1_gen, pops_NCD1_gen_IS)
#####
#The End








