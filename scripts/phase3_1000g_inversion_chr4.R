##########################################################################################
#	Use Phase 1 1000G data to verify allele frequencies within inversion in chr4
#	Bárbara Bitarello
#	Created: 14.02.2017
#	Last modified: 09.08.2017
##########################################################################################
#preamble #######################################################
library(pegas);library(dplyr)
library(plyr);library(data.table)
library(parallel);library(lattice)
library(SOAR);Sys.setenv(R_LOCAL_CACHE="inversions")
library(ggplot2);library(splitstackshape)
library(pryr)
library(doMC); library(purrr)
registerDoMC(4);library(bigmemory);
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/mclapply2.R')
ecdf_fun <- function(x,perc) ecdf(x)(perc) #http://stats.stackexchange.com/questions/50080/estimate-quantile-of-value-in-a-vector
source('/mnt/sequencedb/PopGen/barbara/NCD-Statistics/scripts/NCD_func.R')
fread('input_files/RegionsToAnalyse_45inversions_hg19.txt')-> all_inv
all_inv[,CHR:= Chromosome]
all_inv[,Chromosome:=NULL]
fread('input_files/VariantsClassified_HsInv0102_using434Ph3Samples.txt')-> var_dt  #in light of coordinates Carla sent on 21.07, I should filter this to contain only variables within HsInv102.
setDT(filter(var_dt, POS>= as.numeric(((all_inv %>% dplyr::filter(Inversion=="HsInv0102"))$StartNonRecombining)), POS<=as.numeric(((all_inv %>% dplyr::filter(Inversion=="HsInv0102"))$EndNonRecombining))))-> var_dt
#check to see if there are any SNPs between the outside and inside coordiantes, which should be excluded (see email):
c(filter(var_dt, POS>=40235026, POS<=40235029)$POS, filter(var_dt, POS>=40237058, POS<=40237061)$POS)
#there aren't any.
pops<-c("ACB","ASW","BEB","CDX","CEU","CHB","CHS","CLM","ESN","FIN","GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL","PEL","PJL","PUR","STU","TSI","YRI")
Objects()
fread('zcat /mnt/sequencedb/PopGen/cesare/hg19/bedfiles/hg19.pantro2.Map50_100.TRF.SDs.bed.gz')-> H_C_cov #addes in 8.9.17
#######################################
##########can be skipped ##############
#######################################
FD_list<-factor('list',22)
mclapply2(1:22, function(i)
paste0('zcat < /mnt/sequencedb/PopGen/barbara/NCV_dir_package/input_data/outgroup_files/fds.chr', i, '.hg19_pantro2.Map50_100.TRF.SDs.bed.gz') %>%
data.table::fread(sep='\t')) -> FD_list

for (i in 1:22){
colnames(FD_list[[i]])<-c('CHR', 'POS', 'REF', 'Chimp_REF')
FD_list[[i]] %>% mutate(ID=paste0(CHR, "|",POS)) %>% as.data.table -> FD_list[[i]]
}
#use uppercase
mclapply2(1:22, function(i)
#HC_div1<-HC_div %>% dplyr::mutate(REF=toupper(REF), Chimp_REF=toupper(Chimp_REF)) %>% dplyr::filter(POS>=40224426 & POS<= 40247234) %>% as.data.table  #select range of inversions +- 10 kb
FD_list[[i]] %>% dplyr::mutate(REF=toupper(REF), Chimp_REF=toupper(Chimp_REF)) %>% as.data.table) -> FD_list
Store(FD_list)
#SFS
source('DAF_SFS.R')# this script plots DAF SFS for the inversion in chr4

#the following was done wihtout removing indels. Might need to redo properly.
source('run_vcftools.R') #run vctools and clean data
#
source('read_1000gph3_data.R')
################################
################################
pops<-c("ACB","ASW","BEB","CDX","CEU","CHB","CHS","CLM","ESN","FIN","GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL","PEL","PJL","PUR","STU","TSI","YRI")
names(POPS_AF)<- pops

system.time(source('run_NCD1.r')) #working on this one
system.time(source('run_NCD2.r')) #need to re-run 09.08.17. Having problems with NCD2. Need to fix and re-run every NCD2 instance.

Objects()
gc()
# need to re-run from this point forward using W=2000 and S=1000 08.08.17
HsInv0102_list<-vector('list',2)
names(HsInv0102_list)<-c('NCD1','NCD2')

mclapply2(1:26, function(Z) NCD1(X=POPS_AF[[Z]][[4]][POS %in% var_dt$POS], W=2000, S=1000))-> HsInv0102_list[[1]]
mclapply2(1:26, function(Z) NCD2(X=POPS_AF[[Z]][[4]][POS %in% var_dt$POS],Y=FD_list[[4]], W=2000, S=1000))-> HsInv0102_list[[2]] #NCD2 yielding very strange results!! NCD2>0.5! Ok, the problem is I didn't remove inversions etc.
#
mclapply2(HsInv0102_list[[1]], function(Z) Z[,c("Chr", "POS1","POS2") := tstrsplit(Win.ID, "|", fixed=TRUE)])->HsInv0102_list[[1]]
mclapply2(HsInv0102_list[[2]], function(Z) Z[,c("Chr", "POS1","POS2") := tstrsplit(Win.ID, "|", fixed=TRUE)])->HsInv0102_list[[2]]
#
mclapply2(HsInv0102_list[[1]], function(Z) Z[order(as.numeric(Chr), as.numeric(POS1))])-> HsInv0102_list[[1]]
mclapply2(HsInv0102_list[[2]], function(Z) Z[order(as.numeric(Chr), as.numeric(POS1))])-> HsInv0102_list[[2]]
#
HsInv0102_privstd_list<-vector('list',2)
names(HsInv0102_privstd_list)<-c('NCD1','NCD2')
mclapply2(1:26, function(Z) NCD1(X=POPS_AF[[Z]][[4]][POS %in% var_dt[Type=='PrivateStd']$POS], W=2000, S=1000))-> HsInv0102_privstd_list[[1]]
mclapply2(1:26, function(Z) NCD2(X=POPS_AF[[Z]][[4]][POS %in% var_dt[Type=='PrivateStd']$POS], Y=FD_list[[4]], W=2000, S=1000))-> HsInv0102_privstd_list[[2]]
#
mclapply2(HsInv0102_privstd_list[[1]], function(Z) Z[,c("Chr", "POS1","POS2") := tstrsplit(Win.ID, "|", fixed=TRUE)])-> HsInv0102_privstd_list[[1]]
mclapply2(HsInv0102_privstd_list[[2]], function(Z) Z[,c("Chr", "POS1","POS2") := tstrsplit(Win.ID, "|", fixed=TRUE)])-> HsInv0102_privstd_list[[2]]
#
mclapply2(HsInv0102_privstd_list[[1]], function(Z) Z[order(as.numeric(Chr), as.numeric(POS1))])-> HsInv0102_privstd_list[[1]]
mclapply2(HsInv0102_privstd_list[[2]], function(Z) Z[order(as.numeric(Chr), as.numeric(POS1))])-> HsInv0102_privstd_list[[2]]
#
Store(HsInv0102_list, HsInv0102_privstd_list)

#tf0,5
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD1_gen_IS[[x]]$NCD1_tf0.5, min(HsInv0102_list[[1]][[x]]$NCD1_tf0.5, na.rm=T))[1])) <=0.005)] #this is a very naive approach but shwos some pops have low ncd for this inversion region.
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD2_gen_IS[[x]]$NCD2_tf0.5, min(HsInv0102_list[[2]][[x]]$NCD2_tf0.5, na.rm=T))[1])) <=0.005)] 
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD1_gen_IS[[x]]$NCD1_tf0.5, min(HsInv0102_list[[1]][[x]]$NCD1_tf0.5, na.rm=T))[1])) <=0.0005)]
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD2_gen_IS[[x]]$NCD2_tf0.5, min(HsInv0102_list[[2]][[x]]$NCD2_tf0.5, na.rm=T))[1])) <=0.0005)]

pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD1_gen_IS[[x]]$NCD1_tf0.5, min(HsInv0102_privstd_list[[1]][[x]]$NCD1_tf0.5, na.rm=T))[1])) <=0.005)] #this is a very naive approach but shwos some pops have low ncd for this inversion
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD2_gen_IS[[x]]$NCD2_tf0.5, min(HsInv0102_privstd_list[[2]][[x]]$NCD2_tf0.5, na.rm=T))[1])) <=0.005)]
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD1_gen_IS[[x]]$NCD1_tf0.5, min(HsInv0102_privstd_list[[1]][[x]]$NCD1_tf0.5, na.rm=T))[1])) <=0.0005)] #this is a very naive approach but shwos some pops have low ncd for this inversion
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD2_gen_IS[[x]]$NCD2_tf0.5, min(HsInv0102_privstd_list[[2]][[x]]$NCD2_tf0.5, na.rm=T))[1])) <=0.0005)]

#tf0.4
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD1_gen_IS[[x]]$NCD1_tf0.4, min(HsInv0102_list[[1]][[x]]$NCD1_tf0.4, na.rm=T))[1])) <=0.005)] #this is a very naive approach but shwos some pops have low ncd for this inversion region.
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD2_gen_IS[[x]]$NCD2_tf0.4, min(HsInv0102_list[[2]][[x]]$NCD2_tf0.4, na.rm=T))[1])) <=0.005)]
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD1_gen_IS[[x]]$NCD1_tf0.4, min(HsInv0102_list[[1]][[x]]$NCD1_tf0.4, na.rm=T))[1])) <=0.0005)] #this is a very naive approach but shwos some pops have low ncd for this inversion region.
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD2_gen_IS[[x]]$NCD2_tf0.4, min(HsInv0102_list[[2]][[x]]$NCD2_tf0.4, na.rm=T))[1])) <=0.0005)]

pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD1_gen_IS[[x]]$NCD1_tf0.4, min(HsInv0102_privstd_list[[1]][[x]]$NCD1_tf0.4, na.rm=T))[1])) <=0.005)] #this is a very naive approach but shwos some pops have low ncd for this inversion
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD2_gen_IS[[x]]$NCD2_tf0.4, min(HsInv0102_privstd_list[[2]][[x]]$NCD2_tf0.4, na.rm=T))[1])) <=0.005)]
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD1_gen_IS[[x]]$NCD1_tf0.4, min(HsInv0102_privstd_list[[1]][[x]]$NCD1_tf0.4, na.rm=T))[1])) <=0.0005)] #
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD2_gen_IS[[x]]$NCD2_tf0.4, min(HsInv0102_privstd_list[[2]][[x]]$NCD2_tf0.4, na.rm=T))[1])) <=0.0005)]

#tf0,3
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD1_gen_IS[[x]]$NCD1_tf0.3, min(HsInv0102_list[[1]][[x]]$NCD1_tf0.3, na.rm=T))[1])) <=0.005)] #
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD2_gen_IS[[x]]$NCD2_tf0.3, min(HsInv0102_list[[2]][[x]]$NCD2_tf0.3, na.rm=T))[1])) <=0.005)]
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD1_gen_IS[[x]]$NCD1_tf0.3, min(HsInv0102_list[[1]][[x]]$NCD1_tf0.3, na.rm=T))[1])) <=0.0005)] #this is a very naive approach but shwos some pops have low ncd for this inversion region.
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD2_gen_IS[[x]]$NCD2_tf0.3, min(HsInv0102_list[[2]][[x]]$NCD2_tf0.3, na.rm=T))[1])) <=0.0005)]

pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD1_gen_IS[[x]]$NCD1_tf0.3, min(HsInv0102_privstd_list[[1]][[x]]$NCD1_tf0.3, na.rm=T))[1])) <=0.005)] #this is a very naive approach but shwos some pops have low ncd for this inversion
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD2_gen_IS[[x]]$NCD2_tf0.3, min(HsInv0102_privstd_list[[2]][[x]]$NCD2_tf0.3, na.rm=T))[1])) <=0.005)]
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD1_gen_IS[[x]]$NCD1_tf0.3, min(HsInv0102_privstd_list[[1]][[x]]$NCD1_tf0.3, na.rm=T))[1])) <=0.0005)] #this is a very naive approach but shwos some pops have low ncd for this inversion
pops[which(unlist(mclapply2(1:26, function(x) ecdf_fun(pops_NCD2_gen_IS[[x]]$NCD2_tf0.3, min(HsInv0102_privstd_list[[2]][[x]]$NCD2_tf0.3, na.rm=T))[1])) <=0.0005)]

#summary

nrow(FD_list[[4]][POS>=var_dt[1,POS] & POS<= var_dt[nrow(var_dt),POS]]) # 34
sapply(1:26, function(Z) nrow(POPS_AF[[Z]][[4]][POS>=var_dt[1,POS] & POS<= var_dt[nrow(var_dt),POS]])) #57 SNPs in inversion
sapply(1:26, function(Z) nrow(POPS_AF[[Z]][[4]][POS>=var_dt[1,POS] & POS<= var_dt[nrow(var_dt),POS]][MAF!=0 & MAF!=1])) # between  12 and 19 SNPs actually segregating in each pop in inversion.

## of segregating SNPs per pop.
#
sapply(1:26, function(Z) nrow(POPS_AF[[Z]][[4]][POS %in% var_dt[Type=='PrivateStd']$POS])) #27 SNPs are privatestd
sapply(1:26, function(Z) nrow(POPS_AF[[Z]][[4]][POS %in% var_dt[Type=='PrivateStd']$POS][MAF!=0 & MAF!=1])) #13-17 actually segregating in each pop SNPs

#actually segregating
sapply(1:26, function(Z) setkey(POPS_AF[[Z]][[4]], ID))
setkey(FD_list[[4]], ID)

#PtoDper pop
sapply(1:26, function(Z) HsInv0102_list[[2]][[Z]][,PtoD]) #all PtoDs are <1 which explains why the NCD2 is not low, although NCD1 is.

#So, corrrected PtoD for privateStd
sapply(1:26, function(Z) HsInv0102_privstd_list[[2]][[Z]][,PtoD]) #the same applies. All <1.

#Project: consider extending all of the above to the other inversions as well...

#fread('SummaryInfo_45Invs_Frequencies_v3.1.csv')-> all_inv #on 21.07 Carla sent a new file with inversion coordinates so I am updating this. See preamble, above.
##########################################################
##########################################################
# # # # # START HERE # # # # # # # # START HERE # # # #  # 
Objects()

POPS_AF[[1]][[1]]

#names(POPS_AF)<- pops

############ ############## ############### ################## #######################
pops<-c("ACB","ASW","BEB","CDX","CEU","CHB","CHS","CLM","ESN","FIN","GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL","PEL","PJL","PUR","STU","TSI","YRI")

all_inv[CHR!='X' & CHR!="Y"]-> all_inv_notXY
all_inv_notXY[,.(.N),.(CHR)]

all_inv_notXY[order(as.numeric(CHR))]-> all_inv_notXY

a<-as.numeric(unique(all_inv_notXY$CHR))

all_inv_notXY[,Len:=as.numeric(EndNonRecombining)-as.numeric(StartNonRecombining)]
na.omit(all_inv_notXY)[, Type:='long']-> all_inv_notXY_b
all_inv_notXY[is.na(Len)==T][,Len:= SizeInside][,Type:='short']-> all_inv_notXY_a
rbind(all_inv_notXY_a, all_inv_notXY_b)-> all_inv_notXY

all_inv_notXY[Len>=2000][order(as.numeric(CHR))]->inv_2000bp #to relate to the scan

inv_2000bp[,Inversion]-> tmp_all_inv
####
#make a function
func_inv<-function(inv, dt,pop){
type<-dt[Inversion==inv][,Type]
as.numeric(dt[Inversion==inv][,CHR])-> chr
if(type=='long'){
as.numeric(dt[Inversion==inv][,StartNonRecombining])-> st
as.numeric(dt[Inversion==inv][,EndNonRecombining])-> end
POPS_AF[[J]][[chr]][POS>=st & POS<=end]-> res
tmp<-unique(c(res[POS>= dt[Inversion==inv,EndOutside] & POS<=dt[Inversion==inv,StartInside]], res[POS>=dt[Inversion==inv,EndInside] & POS<=dt[Inversion==inv,StartOutside]]))
if(length(tmp)!=0){
res[!(POS %in% tmp)]-> res
}
else{
res<-res
}
}
if(type=='short'){
as.numeric(dt[Inversion==inv][,StartInside])-> st
as.numeric(dt[Inversion==inv][,EndInside])-> end
POPS_AF[[J]][[chr]][POS>=st & POS<=end]-> res
}
NCD1(X=res, W=2000, S=1000)-> ncd1_res
NCD2(X=res, Y=FD_list[[chr]],  W=2000, S=1000)-> ncd2_res
return(list(res, ncd1_res, ncd2_res))
}
#######################
res_all_inv<-vector('list', 26)
system.time(for(J in 1:26){
mclapply2(1:33, function(Z) func_inv(tmp_all_inv[Z], inv_2000bp, pop=pops[J]))-> res_all_inv[[J]]})

names(res_all_inv)<-pops
for(J in 1:26){
names(res_all_inv[[J]])<- inv_2000bp[,Inversion]
}

for(J in 1:26){
res_all_inv[[J]][inv_2000bp[,Inversion]]
}
Store(res_all_inv)
#this is where the problem is!!!!
for(J in 1:26){
for (K in 1:length(res_all_inv[[J]])){
res_all_inv[[J]][[K]][[3]][,IS:=N_SNPs_cor+N_FDs_cor]-> res_all_inv[[J]][[K]][[3]]
}}
#HsInv0102_privstd_list[[2]][[J]][,IS:=N_SNPs_cor+N_FDs_cor]-> HsInv0102_privstd_list[[2]][[J]]
#for(J in 1:26){
#res_all_inv[[J]][[33]][[3]] %>% dplyr::select(Win.ID, PtoD, N_FDs_Raw, N_FDs_cor, N_Raw, N_SNPs_cor, NCD2_tf0.5:NCD2_tf0.3, IS) %>% as.data.table ->  res_all_inv[[J]][[33]][[3]]
#}

Store(HsInv0102_privstd_list)
Store(res_all_inv)


#stopped here 17.08.2017
dt_ncd<-vector('list',26)
names(dt_ncd)<-pops
for(J in 1:26){ #attention 
	setDT(do.call(rbind,mclapply2(names(res_all_inv[[J]]), function(K) cbind((res_all_inv[[J]][[K]][[2]] %>% dplyr::summarise(ncd1_rawp_0.5=ecdf_fun(pops_NCD1_gen_IS[[J]]$NCD1_tf0.5, min(res_all_inv[[J]][[K]][[2]][N_SNPs_cor>=10]$NCD1_tf0.5, na.rm=T)), 
	ncd1_rawp_0.4= ecdf_fun(pops_NCD1_gen_IS[[J]]$NCD1_tf0.4, min(res_all_inv[[J]][[K]][[2]][N_SNPs_cor>=10]$NCD1_tf0.4, na.rm=T)),
	 ncd1_rawp_0.3= ecdf_fun(pops_NCD1_gen_IS[[J]]$NCD1_tf0.3, min(res_all_inv[[J]][[K]][[2]][N_SNPs_cor>=10]$NCD1_tf0.3, na.rm=T)))),res_all_inv[[J]][[K]][[3]] %>% dplyr::summarise(ncd2_rawp_0.5= ecdf_fun(pops_NCD2_gen_IS[[J]]$NCD2_tf0.5, min(res_all_inv[[J]][[K]][[3]][IS>=10]$NCD2_tf0.5, na.rm=T)), ncd2_rawp_0.4 = ecdf_fun(pops_NCD2_gen_IS[[J]]$NCD2_tf0.4, min(res_all_inv[[J]][[K]][[3]][IS>=10]$NCD2_tf0.4, na.rm=T)), ncd2_rawp_0.3 = ecdf_fun(pops_NCD2_gen_IS[[J]]$NCD2_tf0.3, min(res_all_inv[[J]][[K]][[3]][IS>=10]$NCD2_tf0.3, na.rm=T)), Inversion=K) %>% as.data.table))))-> dt_ncd[[J]];
dt_ncd[[J]][, POP:=pops[J]]
print(pops[J])}

names(dt_ncd)<-pops
Store(dt_ncd)

#first, add length of inversions to dtu_ncd
var1a<-lapply(names(res_all_inv[[1]]), function(X) which(dt_ncd[[1]][,Inversion]==X))


dt_ncd_3<-dt_ncd
for(J in 1:26){
setkey(dt_ncd[[J]], Inversion)
setkey(inv_2000bp, Inversion)
dt_ncd[[J]][inv_2000bp]-> dt_ncd_3[[J]]
}

#have inversion in the same order everywhere

var1a<-sapply(names(res_all_inv[[1]]), function(X) which(dt_ncd[[1]][,Inversion]==X))

lapply(1:26, function(X)  dt_ncd[[X]][var1a])-> dt_ncd

lapply(1:26, function(X)  dt_ncd_3[[X]][var1a])-> dt_ncd_3

Store(dt_ncd, dt_ncd_3)
Objects()
#stopped here bionc06
for (J in 1:26){
rbind(dt_ncd_3[[J]], list(ncd1_rawp_0.5=ecdf_fun(pops_NCD1_gen_IS[[J]]$NCD1_tf0.5, min(HsInv0102_privstd_list[[1]][[J]][N_SNPs_cor>=10]$NCD1_tf0.5, na.rm=T)), ncd1_rawp_0.4=ecdf_fun(pops_NCD1_gen_IS[[J]]$NCD1_tf0.4, min(HsInv0102_privstd_list[[1]][[J]][N_SNPs_cor>=10]$NCD1_tf0.4, na.rm=T)), ncd1_rawp_0.3=ecdf_fun(pops_NCD1_gen_IS[[J]]$NCD1_tf0.3, min(HsInv0102_privstd_list[[1]][[J]][N_SNPs_cor>=10]$NCD1_tf0.3, na.rm=T)),ncd2_rawp_0.5=ecdf_fun(pops_NCD2_gen_IS[[J]]$NCD2_tf0.5, min(HsInv0102_privstd_list[[2]][[J]][IS>=10]$NCD2_tf0.5, na.rm=T)), ncd2_rawp_0.4=ecdf_fun(pops_NCD2_gen_IS[[J]]$NCD2_tf0.4, min(HsInv0102_privstd_list[[2]][[J]][IS>=10]$NCD2_tf0.4, na.rm=T)), ncd2_rawp_0.3=ecdf_fun(pops_NCD2_gen_IS[[J]]$NCD2_tf0.3, min(HsInv0102_privstd_list[[2]][[J]][IS>=10]$NCD2_tf0.3, na.rm=T)), Inversion='HsInv0102_PrivStd',POP=pops[J], StartNonRecombining=40235009, EndOutside=40235026, StartInside=40235029, EndInside=40237058, StartOutside=40237061, EndNonRecombining=40237631 ,SizeInside=2030, SizeUsable=2619, CHR=4,  Len=2622, Type='long'))-> dt_ncd_3[[J]]
print(pops[J])}

for (J in 1:26){
rbind(dt_ncd_3[[J]], list(ncd1_rawp_0.5=NA, ncd1_rawp_0.4=NA, ncd1_rawp_0.3=NA,ncd2_rawp_0.5=NA, ncd2_rawp_0.4=NA, ncd2_rawp_0.3=NA, Inversion='HsInv0102_homozPh3',POP=pops[J], StartNonRecombining=40235009, EndOutside=40235026, StartInside=40235029, EndInside=40237058, StartOutside=40237061, EndNonRecombining=40237631 ,SizeInside=2030, SizeUsable=2619, CHR=4,  Len=2622, Type='long'))-> dt_ncd_3[[J]]
print(pops[J])}

names(dt_ncd_3)<-pops[1:26]

Store(res_all_inv)
Store(dt_ncd)
Store(dt_ncd_3)
# corrected p-values
#stopped here on 15.08.2017

#sample a number between 1st and final positions in chromosome. add +- Len_Inv/2 to each side. Thats the fake inversion.
tmp<- mclapply2(1:22, function(x) seq(POPS_AF[[1]][[x]][order(as.numeric(CHR))][2000, POS], POPS_AF[[1]][[x]][order(as.numeric(CHR))][nrow(POPS_AF[[1]][[x]])-2000, POS]), mc.cores=1)
Store(tmp)


#for(i in 1:26){
#setDT(dt_ncd_3[[i]])

#dt_ncd_3[[i]][order(as.numeric(CHR))][Inversion!="HsInv0102 PrivStd"][order(as.numeric(CHR))]-> dt_ncd[[i]]
#}
dt_ncd_3[[1]][,Inversion]-> inversions
Store(inversions)
#####################################
func_inv3<-function(inv, dt, pop, rep){
type<-dt[[pop]][Inversion==inv][,Type]
len<- dt[[pop]][Inversion==inv][,Len]
as.numeric(dt[[pop]][Inversion==inv][,CHR])-> chr
res<-data.table(CHR=chr, Central_Pos=sample(tmp[[chr]], 1))[,POS:=Central_Pos-len/2][,POS2:=Central_Pos+ len/2]
res2<- POPS_AF_sub[[pop]][[chr]][POS>= res$POS & POS<= res$POS2]
NCD1(X=res2, W=2000, S=1000)-> ncd1_res
NCD2(X=res2, Y=FD_list[[chr]],  W=2000, S=1000)-> ncd2_res
setkey(ncd1_res, Win.ID)
setkey(ncd2_res, Win.ID)
merge(ncd1_res, ncd2_res) -> ncd_res
ncd_res[,POP:=pop][,Inversion:=inv][,Rep:=rep]-> ncd_res
ncd_res-> ncd_res
ncd_res
remove(res)
remove(res2)
return(ncd_res)
}
#
func_inv4<-function(inv, dt, pop, rep){
type<-dt[[pop]][Inversion==inv][,Type]
len<- dt[[pop]][Inversion==inv][,Len]
as.numeric(dt[[pop]][Inversion==inv][,CHR])-> chr
res<-data.table(CHR=chr, Central_Pos=sample(tmp[[chr]], 1))[,POS:=Central_Pos-len/2][,POS2:=Central_Pos+ len/2]
res2<- POPS_AF_sub[[pop]][POS>= res$POS & POS<= res$POS2]
NCD1(X=res2, W=2000, S=1000)-> ncd1_res
NCD2(X=res2, Y=FD_list[[chr]],  W=2000, S=1000)-> ncd2_res
setkey(ncd1_res, Win.ID)
setkey(ncd2_res, Win.ID)
merge(ncd1_res, ncd2_res) -> ncd_res
ncd_res[,POP:=pop][,Inversion:=inv][,Rep:=rep]-> ncd_res
ncd_res-> ncd_res
ncd_res
remove(res)
remove(res2)
return(ncd_res)
}
#


#bionc01 after samp3
POPS_AF[c(1,2,3,4,5)]-> POPS_AF_1_2_3_4_5
mem_change(remove(POPS_AF))
Store(POPS_AF_1_2_3_4_5)
POPS_AF_1_2_3_4_5-> POPS_AF_sub
head(tmp[[1]])
samp<-vector('list', 5)
samp[[1]]<- vector('list', 33)
samp[[2]]<- vector('list', 33)
samp[[3]]<- vector('list', 33)
samp[[4]]<- vector('list', 33)
samp[[5]]<- vector('list', 33)
system.time(for(J in 1:5){  #RUNNING BIONC06
for (K in 1:33){
samp[[J]][[K]]<-foreach(X=1:1000) %dopar%
try(func_inv3(dt=dt_ncd, inv=inversions[K], pop=pops[J], rep=X))
}
gc();
print(pops[J]);
print('done\n');}
)
Store(samp)
##not all re-sampling worked

for(J in 1:5){#bionc01 (running second round of this). need to check because the sapply below is not correct I think.
for(K in 1:33){
1000-length(samp[[J]][[K]][-which(sapply(samp[[J]][[K]], is.null))])->a
for(i in a){
samp[[J]][[K]]<-append(samp[[J]][[K]], 
foreach(X=1:i) %dopar%
try(func_inv3(dt=dt_ncd, inv=inversions[K], pop=pops[J], rep=1000+X)))
}
}}
#check first and then
lapply(1:5, function(Y) sapply(1:33, function(Z) sum(sapply(samp[[Y]][[Z]], function(X) length(X)==17))))
#still some missing
for(J in 1:5){#bionc01 (running second round of this). Fixed code in this round:
for(K in 1:33){
1000-sum(sapply(samp[[J]][[K]], function(X) length(X)==17))->a
for(i in a){
if(i>0){
samp[[J]][[K]]<-append(samp[[J]][[K]],
foreach(X=1:i) %dopar%
try(func_inv3(dt=dt_ncd, inv=inversions[K], pop=pops[J], rep=1000+X)))
}
}}

Store(samp)
#still missing a few

for(J in 1:5){#bionc01 (running second round of this). Fixed code in this round:
for(K in 1:33){
1000-sum(sapply(samp[[J]][[K]], function(X) length(X)==17))->a
for(i in a){
if(i>0){
samp[[J]][[K]]<-append(samp[[J]][[K]],
foreach(X=1:i) %dopar%
try(func_inv3(dt=dt_ncd, inv=inversions[K], pop=pops[J], rep=1000+X)))
}
}}

Store(samp)
#
samp<-mclapply2(1:5, function(Z) lapply(1:33, function(K) samp[[Z]][[K]][which(sapply(samp[[Z]][[K]], function(X)  length(X)==17))]))
names(samp)<- pops[1:5]
for(J in 1:5){
names(samp[[J]])<- inversions[1:33]
}

Store(samp)
#
for(J in 1:5){
append(append(samp[[J]],vector('list', 1000)), vector('list',1000))-> samp[[J]]
}
#
for(J in 1:5){#bionc01 (running second round of this). Fixed code in this round:
K<-34
1000-sum(sapply(samp[[J]][[K]], function(X) length(X)==17))->a
for(i in a){
if(i>0){
samp[[J]][[K]]<-append(samp[[J]][[K]],foreach(X=1:i) %dopar%
try(func_inv3(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=1000+X)))
}
}}

samp<-mclapply2(1:5, function(Z) lapply(1:34, function(K) samp[[Z]][[K]][which(sapply(samp[[Z]][[K]], function(X)  length(X)==17))]))
names(samp)<- pops[1:5]
for(J in 1:5){
names(samp[[J]])<- inversions[1:34]
}
#
Store(samp)
remove(POPS_AF_sub)
##bionc09
POPS_AF_sub<-hsinv0102[c(1,2,3,4,5)]
for(J in 1:5){
append(samp[[J]],vector('list', 1000))-> samp[[J]]
}
for(J in 1:5){
K<-35
samp[[J]][[K]]<-append(samp[[J]][[K]],foreach(X=1:1000) %dopar%
try(func_inv4(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=1000+X)))
}
for(J in 1:5){
K<-35
1000-sum(sapply(samp[[J]][[K]], function(X) length(X)==17))->a
for(i in a){
if(i>0){
samp[[J]][[K]]<-append(samp[[J]][[K]],foreach(X=1:i) %dopar%
try(func_inv4(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=1000+X)))
}
samp<-mclapply2(1:5, function(Z) lapply(1:35, function(K) samp[[Z]][[K]][which(sapply(samp[[Z]][[K]], function(X)  length(X)==17))]))
names(samp)<- pops[1:5]
for(J in 1:5){
names(samp[[J]])<- inversions[1:35]


#################################################
#################################################
#bionc06: done
POPS_AF[c(6,7,8,9,10)]-> POPS_AF_6_7_8_9_10
mem_change(remove(POPS_AF))
Store(POPS_AF_6_7_8_9_10)
POPS_AF_6_7_8_9_10-> POPS_AF_sub
head(tmp[[1]])
samp1<-vector('list', 5)
samp1[[1]]<- vector('list', 33)
samp1[[2]]<- vector('list', 33)
samp1[[3]]<- vector('list', 33)
samp1[[4]]<- vector('list', 33)
samp1[[5]]<- vector('list', 33)
system.time(for(J in 6:10){  #RUNNING BIONC06: 29204.84 seconds
for (K in 1:33){
samp1[[J-5]][[K]]<-foreach(X=1:1000) %dopar%
try(func_inv3(dt=dt_ncd, inv=inversions[K], pop=pops[J], rep=X))
}
gc();
print(pops[J]);
print('done\n');} 
)
Store(samp1)
#
##not enough resamplings, need to complete
for(J in 6:10){#bionc03
for(K in 1:33){
1000-sum(sapply(samp1[[J-5]][[K]], function(X) length(X)==17))->a
for(i in a){
samp1[[J-5]][[K]]<-append(samp1[[J-5]][[K]], foreach(X=1:i) %dopar%
try(func_inv3(dt=dt_ncd, inv=inversions[K], pop=pops[J], rep=1000+X)))
}
}}
#check first and then
Store(samp1)
#a few missing
for(J in 6:10){#bionc03
for(K in 1:33){
1000-sum(sapply(samp1[[J-5]][[K]], function(X) length(X)==17))->a
for(i in a){
if(i>0){
samp1[[J-5]][[K]]<-append(samp1[[J-5]][[K]], foreach(X=1:i) %dopar%
try(func_inv3(dt=dt_ncd, inv=inversions[K], pop=pops[J], rep=1000+X)))
}
}}}

Store(samp1)

samp1<-mclapply2(1:5, function(Z) lapply(1:33, function(K) samp1[[Z]][[K]][which(sapply(samp1[[Z]][[K]], function(X)  length(X)==17))]))
names(samp1)<- pops[6:10]
for(J in 6:10){
names(samp1[[J-5]])<- inversions[1:33]
}

#

for(J in 6:10){
append(append(samp1[[J-5]],vector('list', 1000)), vector('list',1000))-> samp1[[J-5]]
}
#
for(J in 6:10){#bionc01 (running second round of this). Fixed code in this round:
K<-34
1000-sum(sapply(samp1[[J-5]][[K]], function(X) length(X)==17))->a
for(i in a){
if(i>0){
samp1[[J-5]][[K]]<-append(samp1[[J-5]][[K]],foreach(X=1:i) %dopar%
try(func_inv3(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=X)))
}
}}

samp1<-mclapply2(1:5, function(Z) lapply(1:34, function(K) samp1[[Z]][[K]][which(sapply(samp1[[Z]][[K]], function(X)  length(X)==17))]))
names(samp1)<- pops[6:10]
for(J in 6:10){
names(samp1[[J-5]])<- inversions[1:34]
}
Store(samp1)
#
POPS_AF_sub<-hsinv0102[c(6,7,8,9,10)]
for(J in 6:10){
append(samp1[[J-5]],vector('list', 1000))-> samp1[[J-5]]
}
for(J in 6:10){
K<-35
samp1[[J-5]][[K]]<-append(samp1[[J-5]][[K]],foreach(X=1:1000) %dopar%
try(func_inv4(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=1000+X)))
}
for(J in 6:10){
K<-35
1000-sum(sapply(samp1[[J-5]][[K]], function(X) length(X)==17))->a
for(i in a){
if(i>0){
samp1[[J-5]][[K]]<-append(samp1[[J-5]][[K]],foreach(X=1:i) %dopar%
try(func_inv4(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=1000+X)))}
samp1<-mclapply2(1:5, function(Z) lapply(1:35, function(K) samp1[[Z]][[K]][which(sapply(samp1[[Z]][[K]], function(X)  length(X)==17))]))
names(samp1)<- pops[6:10]
for(J in 6:10){
names(samp1[[J-5]])<- inversions[1:35]
Store(samp1)
##################################################################################################
#bionc09: done
POPS_AF[c(11,12,13,14,15,16)]-> POPS_AF_11_12_13_14_15_16
mem_change(remove(POPS_AF))
Store(POPS_AF_11_12_13_14_15_16)
gc()
mem_used()≈
POPS_AF_11_12_13_14_15_16-> POPS_AF_sub
head(tmp[[1]])
samp2<-vector('list', 6) #: 34606.56
samp2[[1]]<- vector('list', 33)
samp2[[2]]<- vector('list', 33)
samp2[[3]]<- vector('list', 33)
samp2[[4]]<- vector('list', 33)
samp2[[5]]<- vector('list', 33)
samp2[[6]]<- vector('list', 33)
system.time(for(J in 11:16){
for (K in 1:33){
samp2[[J-10]][[K]]<-foreach(X=1:1000) %dopar%
try(func_inv3(dt=dt_ncd, inv=inversions[K], pop=pops[J], rep=X))
}
gc();
print(pops[J]);
print('done\n');} #21214
)
Store(samp2)
#
###not enough resamplings, need to complete
for(J in 11:16){#bionc09
for(K in 1:33){
1000-length(samp2[[J-10]][[K]][-which(sapply(samp2[[J-10]][[K]], is.null))])->a
for(i in a){
samp2[[J-10]][[K]]<-append(samp2[[J-10]][[K]], foreach(X=1:i) %dopar%
try(func_inv3(dt=dt_ncd, inv=inversions[K], pop=pops[J], rep=1000+X)))
}
}
print(pops[J]);}

#check first and then
samp2<-mclapply2(1:6, function(Z) lapply(1:33, function(K) samp2[[Z]][[K]][which(sapply(samp2[[Z]][[K]], function(X)  length(X)==17))]))
names(samp2)<- pops[11:16]
for(J in 11:16){
names(samp2[[J-10]])<- inversions[1:33]
}
Store(samp2) #done.
#

for(J in 11:16){
append(append(samp2[[J-10]],vector('list', 1000)), vector('list',1000))-> samp2[[J-10]]
}
#
for(J in 11:16){#bionc01 (running second round of this). Fixed code in this round:
K<-34
1000-sum(sapply(samp2[[J-10]][[K]], function(X) length(X)==17))->a
for(i in a){
if(i>0){
samp2[[J-10]][[K]]<-append(samp2[[J-10]][[K]],foreach(X=1:i) %dopar%
try(func_inv3(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=X)))
}
}}

samp2<-mclapply2(1:6, function(Z) lapply(1:34, function(K) samp2[[Z]][[K]][which(sapply(samp2[[Z]][[K]], function(X)  length(X)==17))]))
names(samp2)<- pops[11:16]
for(J in 11:16){
names(samp2[[J-10]])<- inversions[1:34]
}

Store(samp2)
#
POPS_AF_sub<-hsinv0102[c(11,12,13,14,15,16)]
for(J in 11:16){
append(samp2[[J-10]],vector('list', 1000))-> samp2[[J-10]]
}
for(J in 11:16){
K<-35
samp2[[J-10]][[K]]<-append(samp2[[J-10]][[K]],foreach(X=1:1000) %dopar%
try(func_inv4(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=X)))
}
for(J in 11:16){
K<-35
1000-sum(sapply(samp2[[J-10]][[K]], function(X) length(X)==17))->a
for(i in a){
if(i>0){
samp2[[J-10]][[K]]<-append(samp2[[J-10]][[K]],foreach(X=1:i) %dopar%
try(func_inv4(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=1000+X)))}}}

samp2<-mclapply2(1:6, function(Z) lapply(1:35, function(K) samp2[[Z]][[K]][which(sapply(samp2[[Z]][[K]], function(X)  length(X)==17))]))
names(samp2)<- pops[11:16]
for(J in 11:16){
names(samp2[[J-10]])<- inversions[1:35]}
Store(samp2)
##########################################################
POPS_AF[c(17,18,19)]-> POPS_AF_17_18_19
mem_change(remove(POPS_AF))
Store(POPS_AF_17_18_19)
gc()
mem_used()
POPS_AF_17_18_19[[1]]
POPS_AF_17_18_19-> POPS_AF_sub
samp3<-vector('list', 3)  #running bionc01
samp3[[1]]<- vector('list', 33)
samp3[[2]]<- vector('list', 33)
samp3[[3]]<- vector('list', 33)
system.time(for(J in 17:19){
for (K in 1:33){
samp3[[J-16]][[K]]<-foreach(X=1:1000) %dopar%
try(func_inv3(dt=dt_ncd, inv=inversions[K], pop=pops[J], rep=X))
}
gc();
print(pops[J]);
print('done\n');} #
)
Store(samp3)
#
##not enough resamplings, need to complete

for(J in 17:19){# ??
for(K in 1:33){
1000-length(samp3[[J-16]][[K]][-which(sapply(samp3[[J-16]][[K]], is.null))])->a
for(i in a){
if(i>0){
samp3[[J-16]][[K]]<-append(samp3[[J-16]][[K]], foreach(X=1:i) %dopar%
try(func_inv3(dt=dt_ncd, inv=inversions[K], pop=pops[J], rep=1000+X)))
}
}}}
#check first and then
Store(samp3)
#need 

for(J in 17:19){# bion06
for(K in 1:34){
1000-sum(sapply(samp3[[J-16]][[K]], function(X) length(X)==17))->a
for(i in a){
if(i>0){
samp3[[J-16]][[K]]<-append(samp3[[J-16]][[K]], foreach(X=1:i) %dopar%
try(func_inv3(dt=dt_ncd, inv=inversions[K], pop=pops[J], rep=1000+X)))
}
}}}

#
for(J in 17:19){
append(append(samp3[[J-16]],vector('list', 1000)), vector('list',1000))-> samp3[[J-16]]
}

for(J in 17:19){# bion06
for(K in 1:34){
1000-sum(sapply(samp3[[J-16]][[K]], function(X) length(X)==17))->a
for(i in a){
if(i>0){
samp3[[J-16]][[K]]<-append(samp3[[J-16]][[K]], foreach(X=1:i) %dopar%
try(func_inv3(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=1000+X)))
}
}}}
Store(samp3)
#
samp3<-mclapply2(1:3, function(Z) lapply(1:34, function(K) samp3[[Z]][[K]][which(sapply(samp3[[Z]][[K]], function(X)  length(X)==17))]))
names(samp3)<- pops[17:19]
for(J in 17:19){
names(samp3[[J-16]])<- inversions[1:34]
}

Store(samp3)
#

POPS_AF_sub<-hsinv0102[c(17,18,19)]
for(J in 17:19){
append(samp3[[J-16]],vector('list', 1000))-> samp3[[J-16]]
}
for(J in 17:19){
K<-35
samp3[[J-16]][[K]]<-append(samp3[[J-16]][[K]],foreach(X=1:1000) %dopar%
try(func_inv4(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=X)))
}
for(J in 17:19){
K<-35
1000-sum(sapply(samp3[[J-16]][[K]], function(X) length(X)==17))->a
for(i in a){
if(i>0){
samp3[[J-16]][[K]]<-append(samp3[[J-16]][[K]],foreach(X=1:i) %dopar%
try(func_inv4(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=1000+X)))}

samp3<-mclapply2(1:3, function(Z) lapply(1:35, function(K) samp3[[Z]][[K]][which(sapply(samp3[[Z]][[K]], function(X)  length(X)==17))]))
names(samp3)<- pops[17:19]
for(J in 17:19){
names(samp3[[J-16]])<- inversions[1:35]}
##########################################################
#another node: bionc03: done
POPS_AF[c(20,21,22)]-> POPS_AF_20_21_22
mem_change(remove(POPS_AF))
Store(POPS_AF_20_21_22)
gc()
mem_used()
POPS_AF_20_21_22-> POPS_AF_sub
mem_used()
head(tmp[[1]])
mem_used()
samp4<-vector('list', 3)
samp4[[1]]<- vector('list', 33)
samp4[[2]]<- vector('list', 33)
samp4[[3]]<- vector('list', 33)
system.time(for(J in 20:22){ # 24182.99
for (K in 1:33){
samp4[[J-19]][[K]]<-foreach(X=1:1000) %dopar%
try(func_inv3(dt=dt_ncd, inv=inversions[K], pop=pops[J], rep=X))
}
gc();
print(pops[J]);
print('done\n');} 
)
Store(samp4)
#
###not enough resamplings, need to complete
for(J in 20:22){#bionc02
for(K in 1:33){
1000-sum(sapply(samp4[[J-19]][[K]], function(X) length(X)==17))->a
for(i in a){
samp4[[J-19]][[K]]<-append(samp4[[J-19]][[K]], foreach(X=1:i) %dopar%
try(func_inv3(dt=dt_ncd, inv=inversions[K], pop=pops[J], rep=1000+X)))
}
}}
#check
Store(samp4)
#a few missing
for(J in 20:22){#bionc02
for(K in 1:33){
1000-sum(sapply(samp4[[J-19]][[K]], function(X) length(X)==17))->a
for(i in a){
if(i>0){
samp4[[J-19]][[K]]<-append(samp4[[J-19]][[K]], foreach(X=1:i) %dopar%
try(func_inv3(dt=dt_ncd, inv=inversions[K], pop=pops[J], rep=1000+X)))
}
}}}

#

samp4<-mclapply2(1:3, function(Z) lapply(1:33, function(K) samp4[[Z]][[K]][which(sapply(samp4[[Z]][[K]], function(X)  length(X)==17))]))
names(samp4)<- pops[20:22]
for(J in 20:22){
names(samp4[[J-19]])<- inversions[1:33]
}
#

for(J in 20:22){
append(append(samp4[[J-19]],vector('list', 1000)), vector('list',1000))-> samp4[[J-19]]
}
#
for(J in 20:22){#bionc01 (running second round of this). Fixed code in this round:
K<-34
1000-sum(sapply(samp4[[J-19]][[K]], function(X) length(X)==17))->a
for(i in a){
if(i>0){
samp4[[J-19]][[K]]<-append(samp4[[J-19]][[K]],foreach(X=1:i) %dopar%
try(func_inv3(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=X)))
}
}}

samp4<-mclapply2(1:3, function(Z) lapply(1:34, function(K) samp4[[Z]][[K]][which(sapply(samp4[[Z]][[K]], function(X)  length(X)==17))]))
names(samp4)<- pops[20:22]
for(J in 20:22){
names(samp4[[J-19]])<- inversions[1:34]
}

Store(samp4)
#

POPS_AF_sub<-hsinv0102[c(20,21,22)]
for(J in 20:22){
append(samp4[[J-19]],vector('list', 1000))-> samp4[[J-19]]
}
for(J in 20:22){
K<-35
samp4[[J-19]][[K]]<-append(samp4[[J-19]][[K]],foreach(X=1:1000) %dopar%
try(func_inv4(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=X)))
}
for(J in 20:22){
K<-35
1000-sum(sapply(samp4[[J-19]][[K]], function(X) length(X)==17))->a
for(i in a){
if(i>0){
samp4[[J-19]][[K]]<-append(samp4[[J-19]][[K]],foreach(X=1:i) %dopar%
try(func_inv4(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=1000+X)))}
}}
samp4<-mclapply2(1:3, function(Z) lapply(1:35, function(K) samp4[[Z]][[K]][which(sapply(samp4[[Z]][[K]], function(X)  length(X)==17))]))
names(samp4)<- pops[20:22]
for(J in 20:22){
names(samp4[[J-19]])<- inversions[1:35]}

Store(samp4)

##########################################################
POPS_AF[c(23,24,25,26)]->POPS_AF_23_24_25_26 
mem_change(remove(POPS_AF))
gc()
mem_used()
POPS_AF_23_24_25_26 -> POPS_AF_sub
head(tmp[[1]])
samp5<-vector('list', 4)
samp5[[1]]<- vector('list', 33)
samp5[[2]]<- vector('list', 33)
samp5[[3]]<- vector('list', 33)
samp5[[3]]<- vector('list', 33)
system.time(for(J in 23:26){
for (K in 1:33){
samp5[[J-22]][[K]]<-foreach(X=1:1000) %dopar%
try(func_inv3(dt=dt_ncd, inv=inversions[K], pop=pops[J], rep=X))
}
gc();
print(pops[J]);
print('done\n');}
)
Store(samp5)

#need another round
for(J in 23:26){#running in bionc09
for(K in 1:33){
1000-sum(sapply(samp5[[J-22]][[K]], function(X) length(X)==17))->a
for(i in a){
samp5[[J-22]][[K]]<-append(samp5[[J-22]][[K]], foreach(X=1:i) %dopar%
try(func_inv3(dt=dt_ncd, inv=inversions[K], pop=pops[J], rep=1000+X)))
}
}}
Store(samp5)
#check first and then
samp5<-mclapply2(1:4, function(Z) lapply(1:33, function(K) samp5[[Z]][[K]][which(sapply(samp5[[Z]][[K]], function(X)  length(X)==17))]))
names(samp5)<- pops[23:26]
for(J in 23:26){
names(samp5[[J-22]])<- inversions[1:33]
}

Store(samp5)


for(J in 23:26){
append(append(samp5[[J-22]],vector('list', 1000)), vector('list',1000))-> samp5[[J-22]]
}
#
for(J in 23:26){#bionc01 (running second round of this). Fixed code in this round:
K<-34
1000-sum(sapply(samp5[[J-22]][[K]], function(X) length(X)==17))->a
for(i in a){
if(i>0){
samp5[[J-22]][[K]]<-append(samp5[[J-22]][[K]],foreach(X=1:i) %dopar%
try(func_inv3(dt=dt_ncd, inv=inversions[K], pop=pops[J], rep=X)))
}
}}

samp5<-mclapply2(1:4, function(Z) lapply(1:34, function(K) samp5[[Z]][[K]][which(sapply(samp5[[Z]][[K]], function(X)  length(X)==17))]))
names(samp5)<- pops[23:26]
for(J in 23:26){
names(samp5[[J-22]])<- inversions[1:34]
}
Store(samp5)
#
POPS_AF_sub<-hsinv0102[c(23,24,25,26)]
for(J in 23:26){
append(samp5[[J-22]],vector('list', 1000))-> samp5[[J-22]]
}
for(J in 23:26){
K<-35
samp5[[J-22]][[K]]<-append(samp5[[J-22]][[K]],foreach(X=1:1000) %dopar%
try(func_inv4(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=X)))
}
for(J in 23:26){
K<-35
1000-sum(sapply(samp5[[J-22]][[K]], function(X) length(X)==17))->a
for(i in a){
if(i>0){
samp5[[J-22]][[K]]<-append(samp5[[J-22]][[K]],foreach(X=1:i) %dopar%
try(func_inv4(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=1000+X)))}
}}
samp5<-mclapply2(1:4, function(Z) lapply(1:35, function(K) samp5[[Z]][[K]][which(sapply(samp5[[Z]][[K]], function(X)  length(X)==17))]))
names(samp5)<- pops[23:26]
for(J in 23:26){
names(samp5[[J-22]])<- inversions[1:35]}
Store(samp5)

#################################################################

###
#after there is 1,000 resampiings for each inversion:

test<-append(append(append(append(append(samp,samp1),samp2), samp3), samp4),  samp5)

test2<- mclapply2(test, function(X) mclapply2(X, function(Y) Y[1:1000]))

for (J in 1:length(test2)){

lapply(1:35, function(K) lapply(1:1000, function(X) test2[[J]][[K]][[X]][, Rep:=X]))-> test2[[J]]
}

remove(test)

gc()

test3<-lapply(test2, function(X) lapply(X, function(Y) do.call(rbind, Y)))

for(J in length(test3)){

for(K in 1:35){
unique(test3[[J]][[K]]$Rep)-> tp

for(i in 1:1000){
test3[[J]][[K]][Rep==tp[i]][,Rep:= i]
}}}

remove(test2)

test4<- lapply(1:length(test3), function(J) do.call(rbind, test3[[J]])) 


test5<-do.call(rbind,test4)


#now, group by pop and rep, and select lowest p-value for each tf

#dplur

system.time(test5  %>% dplyr::filter(N_SNPs_cor.x>=10) %>% dplyr::group_by(POP,Inversion,Rep) %>% 
dplyr::summarise(min_ncd1_tf0.5=min(NCD1_tf0.5), min_ncd1_tf0.4=min(NCD1_tf0.4), min_ncd1_tf0.3=min(NCD1_tf0.3)) %>%
as.data.table-> test6
)

system.time(test5  %>% dplyr::filter(N_SNPs_cor.x+N_FDs_cor>=10) %>% dplyr::group_by(POP,Inversion,Rep) %>%
dplyr::summarise(min_ncd2_tf0.5=min(NCD2_tf0.5), min_ncd2_tf0.4=min(NCD2_tf0.4), min_ncd2_tf0.3=min(NCD2_tf0.3)) %>%
as.data.table -> test6b)
#for(J in 1:26){names(res_all_inv[[J]])<-inversions[1:33]}

#for(J in 1:26){dt_ncd_3[[J]][,POP:=pops[J]]-> dt_ncd_3[[J]]}
merge(test6b, test6, all.x=T)-> test7
Store(test6b, test6, test7)

dt_ncd_4<-do.call(rbind, dt_ncd_3)
dt_ncd_4 %>% dplyr::select(ncd1_rawp_0.5:POP, StartNonRecombining:Type) %>% as.data.table-> dt_ncd_4
setkey(dt_ncd_4, POP, Inversion)
#test5[[1]][N_SNPs_cor.x>=10][, P_cor_tf0.5:= ecdf_fun(pops_NCD1_gen_IS[[1]][N_SNPs_cor>=10]$NCD1_tf0.5,min(NCD2_tf0.5)), by=c('Inversion','Rep')]


do.call(rbind, 
lapply(pops, function(J) do.call(rbind,
lapply(inversions[1:33],function(i) data.table(Inversion=i,
ncd1_corp_0.5=ecdf_fun(test6[POP==J & Inversion==i][,min_ncd1_tf0.5], min(res_all_inv[[J]][[i]][[2]][N_SNPs_cor>=10][, NCD1_tf0.5], na.rm=T)),ncd1_corp_0.4=ecdf_fun(test6[POP==J & Inversion==i][,min_ncd1_tf0.4], min(res_all_inv[[J]][[i]][[2]][N_SNPs_cor>=10][, NCD1_tf0.4], na.rm=T)),ncd1_corp_0.3=ecdf_fun(test6[POP==J & Inversion==i][,min_ncd1_tf0.3], min(res_all_inv[[J]][[i]][[2]][N_SNPs_cor>=10][, NCD1_tf0.3], na.rm=T)), POP=J)))))->a1


do.call(rbind,
lapply(unique(test7$POP), function(J) do.call(rbind,
lapply(inversions[1:33],function(i) data.table(Inversion=i,
ncd2_corp_0.5=ecdf_fun(test6b[POP==J & Inversion==i][,min_ncd2_tf0.5], min(res_all_inv[[J]][[i]][[3]][IS>=10][, NCD2_tf0.5], na.rm=T)),ncd2_corp_0.4=ecdf_fun(test6b[POP==J & Inversion==i][,min_ncd2_tf0.4], min(res_all_inv[[J]][[i]][[3]][IS>=10][, NCD2_tf0.4], na.rm=T)),ncd2_corp_0.3=ecdf_fun(test6b[POP==J & Inversion==i][,min_ncd2_tf0.3], min(res_all_inv[[J]][[i]][[3]][IS>=10][, NCD2_tf0.3], na.rm=T)),POP=J)))))->a2

merge(a1,a2)-> a

Store(a)
#add 34 and 35 to a
b<-vector('list', 26)
for(J in pops){
data.table(Inversion=inversions[34],
ncd1_corp_0.5=ecdf_fun(test6[POP==J & Inversion==inversions[34]][,min_ncd1_tf0.5], min(HsInv0102_privstd_list[[1]][[J]][N_SNPs_cor>=10][, NCD1_tf0.5], na.rm=T)),
ncd1_corp_0.4=ecdf_fun(test6[POP==J & Inversion==inversions[34]][,min_ncd1_tf0.4], min(HsInv0102_privstd_list[[1]][[J]][N_SNPs_cor>=10][, NCD1_tf0.4], na.rm=T)),
ncd1_corp_0.3=ecdf_fun(test6[POP==J & Inversion==inversions[34]][,min_ncd1_tf0.3], min(HsInv0102_privstd_list[[1]][[J]][N_SNPs_cor>=10][, NCD1_tf0.3], na.rm=T)),
ncd2_corp_0.5=ecdf_fun(test6b[POP==J & Inversion==inversions[34]][,min_ncd2_tf0.5], min(HsInv0102_privstd_list[[2]][[J]][IS>=10][, NCD2_tf0.5], na.rm=T)),
ncd2_corp_0.4=ecdf_fun(test6b[POP==J & Inversion==inversions[34]][,min_ncd2_tf0.4], min(HsInv0102_privstd_list[[2]][[J]][IS>=10][, NCD2_tf0.4], na.rm=T)),
ncd2_corp_0.3=ecdf_fun(test6b[POP==J & Inversion==inversions[34]][,min_ncd2_tf0.3], min(HsInv0102_privstd_list[[2]][[J]][IS>=10][, NCD2_tf0.3], na.rm=T)),POP=J)->b[[J]]
}

do.call(rbind, b)-> b

#

names(hsinv0102_homozph3)<-pops
Store(hsinv0102_homozph3)

test6b[Inversion=="HsInv0102_homozStdPhase3", Inversion:="HsInv0102_homozPh3"]
test6[Inversion=="HsInv0102_homozStdPhase3", Inversion:="HsInv0102_homozPh3"]
d<-vector('list', 26)
for(J in pops){
data.table(Inversion=inversions[35],
ncd1_corp_0.5=ecdf_fun(test6[POP==J & Inversion==inversions[35]][,min_ncd1_tf0.5], min(hsinv0102_homozph3[[J]][[2]][N_SNPs_cor>=10][, NCD1_tf0.5], na.rm=T)),
ncd1_corp_0.4=ecdf_fun(test6[POP==J & Inversion==inversions[35]][,min_ncd1_tf0.4], min(hsinv0102_homozph3[[J]][[2]][N_SNPs_cor>=10][, NCD1_tf0.4], na.rm=T)),
ncd1_corp_0.3=ecdf_fun(test6[POP==J & Inversion==inversions[35]][,min_ncd1_tf0.3], min(hsinv0102_homozph3[[J]][[2]][N_SNPs_cor>=10][, NCD1_tf0.3], na.rm=T)),
ncd2_corp_0.5=ecdf_fun(test6b[POP==J & Inversion==inversions[35]][,min_ncd2_tf0.5], min(hsinv0102_homozph3[[J]][[3]][(N_SNPs_cor+N_FDs_cor)>=10][, NCD2_tf0.5], na.rm=T)),
ncd2_corp_0.4=ecdf_fun(test6b[POP==J & Inversion==inversions[35]][,min_ncd2_tf0.4], min(hsinv0102_homozph3[[J]][[3]][(N_SNPs_cor+N_FDs_cor)>=10][, NCD2_tf0.4], na.rm=T)),
ncd2_corp_0.3=ecdf_fun(test6b[POP==J & Inversion==inversions[35]][,min_ncd2_tf0.3], min(hsinv0102_homozph3[[J]][[3]][(N_SNPs_cor+N_FDs_cor)>=10][, NCD2_tf0.3], na.rm=T)),POP=J)->d[[J]]
}

do.call(rbind, d)-> d
rbind(a,b,d)-> final
#then

setkey(final, POP, Inversion)
etkey(dt_ncd_4, POP, Inversion)

dt_ncd_4[final][order(Inversion)] %>% dplyr::select(POP,Inversion,Type,CHR, StartNonRecombining:SizeUsable, Len,ncd1_rawp_0.5:ncd2_rawp_0.3,ncd1_corp_0.5:ncd2_corp_0.3) %>% as.data.table-> almost_final

remove(test6, test7)

#add calculate and add two additional lines fro HsInv0102:

#TO DO TO DO TO DO

#read table with genotypes

fread("HsInv0102_Genotypes.txt")-> gen_hsinv0102
#select indiviudals with 1000G 0|0 ogenotypes

#save individual's names to file in

for(i in pops){gen_hsinv0102[Ph3Genotype=="0|0"] %>% dplyr::filter(pop==i) %>% dplyr::select(sample) %>% as.data.table->bla
write.table(unique(bla), file=paste0("/mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/", i,"_samples_gen_00_hsinv0102.txt" ), quote=F, row.names=F, col.names=F)
}

foreach(X=1:24) %dopar%
system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/', pops[X],'_samples_gen_00_hsinv0102.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/', pops[X], '_hsinv0102_AF'))


for(X in 1:26){
ta<- paste0('sed \'s/\t/  /\' input_files/', pops[X], '_hsinv0102_AF.frq |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' | sed \'s/  /\t/g\'|sed \'s/ALLELE:/ALLELE_/\'| sed \'s/;/\t/\' > input_files/', pops[X], '_hsinv0102_AF_2.frq')
system(ta)
print(pops[X])
}

#for(X in 1:26){
#ta<-paste0('rm input_files/', pops[X], '_hsinv0102_AF.frq')
#system(ta)
#}
#
#now, read in this data
#start here

hsinv0102<-vector('list', 26)
for(X in 1:26){
fread(paste0('input_files/', pops[X], '_hsinv0102_AF_2.frq'))-> hsinv0102[[X]]
}
names(hsinv0102)<- pops

#clean it up

for (X in 1:26){ #attention
	colnames(hsinv0102[[X]])<-c('CHR','POS','nAL', 'N_Chr','AF');
        hsinv0102[[X]][,ID:=paste0(CHR,"|",POS)][order(POS)] -> hsinv0102[[X]];
	setkey(hsinv0102[[X]], ID); setDT(unique(hsinv0102[[X]]))-> hsinv0102[[X]];
        setkey(hsinv0102[[X]], CHR, POS); setkey(Res_Alt_list[[4]], CHR, POS)
        hsinv0102[[X]][Res_Alt_list[[4]]][order(POS)] -> hsinv0102[[X]];
        setkey(hsinv0102[[X]], ID); unique(hsinv0102[[X]])-> hsinv0102[[X]];
	print(pops[X]);
}
for (X in 1:26){
        hsinv0102[[X]][-(grep("\\b[A-Z]{2,}:\\b",hsinv0102[[X]]$AF)),][order(POS)] -> hsinv0102[[X]];#exclude lines with indels etc
	gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",hsinv0102[[X]]$AF))))-> hsinv0102[[X]]$AF; #clean up and keep only AF
        setDT(hsinv0102[[X]])[, paste0("AF", 1:3) := tstrsplit(AF, ";")]; #split AF into 3 cols
        as.numeric(hsinv0102[[X]]$AF1)-> hsinv0102[[X]]$AF1; as.numeric(hsinv0102[[X]]$AF2)-> hsinv0102[[X]]$AF2; as.numeric(hsinv0102[[X]]$AF3)-> hsinv0102[[X]]$AF3; #make them numeric
        hsinv0102[[X]][,MAF:=pmin(AF1,AF2,AF3, na.rm=T)] -> hsinv0102[[X]];
	hsinv0102[[X]][MAF<=0.5]-> hsinv0102[[X]];
	gc();
	print(pops[X]);
}

for(J in 1:26){ hsinv0102[[J]][, MAF:=as.numeric(MAF)]}

lapply(1:26, function(X) dplyr::select(hsinv0102[[X]], CHR, POS, ID:ALT, MAF) %>% as.data.table)-> hsinv0102
names(hsinv0102)<-pops
Store(hsinv0102)
##select coordinates of inversion:


mclapply2(1:26, function(J) setDT(filter(hsinv0102[[J]], POS>= as.numeric(((all_inv %>% dplyr::filter(Inversion=="HsInv0102"))$StartNonRecombining)), POS<=as.numeric(((all_inv %>% dplyr::filter(Inversion=="HsInv0102"))$EndNonRecombining)))))-> var2_dt


hsinv0102_privstd<-mclapply2(1:26, function(J) func_inv(inv="HsInv0102_PrivStd", dt=dt_ncd_3[[J]], pop=pops[J])) #raw
Store(hsinv0102_privstd)
##

func_inv2<-function(inv, dt,pop){
type<-dt[Inversion==inv][,Type]
as.numeric(dt[Inversion==inv][,CHR])-> chr
#if(type=='long'){
as.numeric(dt[Inversion==inv][,StartNonRecombining])-> st
as.numeric(dt[Inversion==inv][,EndNonRecombining])-> end
hsinv0102[[pop]][POS>=st & POS<=end]-> res
tmp<-unique(c(res[POS>= dt[Inversion==inv,EndOutside] & POS<=dt[Inversion==inv,StartInside]], res[POS>=dt[Inversion==inv,EndInside] & POS<=dt[Inversion==inv,StartOutside]]))
if(length(tmp)!=0){
res[!(POS %in% tmp)]-> res
}
else{
res<-res
}
#if(type=='short'){
#as.numeric(dt[Inversion==inv][,StartInside])-> st#as.numeric(dt[Inversion==inv][,EndInside])-> end#POPS_AF[[J]][[chr]][POS>=st & POS<=end]-> res
NCD1(X=res, W=2000, S=1000)-> ncd1_res
NCD2(X=res, Y=FD_list[[chr]],  W=2000, S=1000)-> ncd2_res
return(list(res, ncd1_res, ncd2_res))
}

##

#nraw p
lapply(1:26, function(X) rbind(rbind(dt_ncd[[X]],dt_ncd[[X]][Inversion=='HsInv0102']), dt_ncd[[X]][Inversion=='HsInv0102']))-> dt_ncd


hsinv0102_homozph3<-mclapply2(1:26, function(J) func_inv2(inv="HsInv0102_homozPh3", dt=dt_ncd_3[[J]], pop=pops[J]))  #raw
Store(hsinv0102_homozph3)

##

for(X in 1:26){
dt_ncd[[X]][34,Inversion:="HsInv0102_PrivStd"];
dt_ncd[[X]][35,Inversion:="HsInv0102_homozStdPh3"];
dt_ncd[[X]][34,ncd1_rawp_0.5:=NA];dt_ncd[[X]][34, ncd2_rawp_0.5:=NA];dt_ncd[[X]][34, ncd1_rawp_0.4:=NA];dt_ncd[[X]][34, ncd2_rawp_0.4:=NA];dt_ncd[[X]][34, ncd1_rawp_0.3:=NA];dt_ncd[[X]][34, ncd2_rawp_0.3:=NA];
dt_ncd[[X]][35, ncd1_corp_0.5:=NA];dt_ncd[[X]][35, ncd2_corp_0.5:=NA];dt_ncd[[X]][35, ncd1_corp_0.4:=NA];dt_ncd[[X]][35, ncd2_corp_0.4:=NA];dt_ncd[[X]][35, ncd1_corp_0.3:=NA];dt_ncd[[X]][35, ncd2_corp_0.3:=NA]
}

for(X in 1:26){
dt_ncd_3[[X]][34, 
ncd1_rawp_0.5:= ecdf_fun(pops_NCD1_gen_IS[[X]]$NCD1_tf0.5, min(HsInv0102_privstd_list[[1]][[X]][N_SNPs_cor>=10]$NCD1_tf0.5, na.rm=T)),
		ncd1_rawp_0.4:= ecdf_fun(pops_NCD1_gen_IS[[X]]$NCD1_tf0.4, min(HsInv0102_privstd_list[[1]][[X]][N_SNPs_cor>=10]$NCD1_tf0.4, na.rm=T)),
		ncd1_rawp_0.3:= ecdf_fun(pops_NCD1_gen_IS[[X]]$NCD1_tf0.3, min(HsInv0102_privstd_list[[1]][[X]][N_SNPs_cor>=10]$NCD1_tf0.3, na.rm=T)),
		ncd2_rawp_0.5:= ecdf_fun(pops_NCD2_gen_IS[[X]]$NCD2_tf0.5, min(HsInv0102_privstd_list[[2]][[X]][(N_SNPs_cor+N_FDs_cor)>=10]$NCD2_tf0.5, na.rm=T)),
		ncd2_rawp_0.4:= ecdf_fun(pops_NCD2_gen_IS[[X]]$NCD2_tf0.4, min(HsInv0102_privstd_list[[2]][[X]][(N_SNPs_cor+N_FDs_cor)>=10]$NCD2_tf0.4, na.rm=T)),
		ncd2_rawp_0.3:= ecdf_fun(pops_NCD2_gen_IS[[X]]$NCD2_tf0.3, min(HsInv0102_privstd_list[[2]][[X]][(N_SNPs_cor+N_FDs_cor)>=10]$NCD2_tf0.3, na.rm=T))]
		print(pops[X]);	
	}	


for(X in 1:26){
dt_ncd_3[[X]][35,
ncd1_rawp_0.5:= ecdf_fun(NCD1_chr4_homozstd[[X]][N_SNPs_cor>=10]$NCD1_tf0.5, min(hsinv0102_homozph3[[X]][[2]][N_SNPs_cor>=10]$NCD1_tf0.5, na.rm=T)),
                ncd1_rawp_0.4:= ecdf_fun(NCD1_chr4_homozstd[[X]][N_SNPs_cor>=10]$NCD1_tf0.4, min(hsinv0102_homozph3[[X]][[2]][N_SNPs_cor>=10]$NCD1_tf0.4, na.rm=T)),
                ncd1_rawp_0.3:= ecdf_fun(NCD1_chr4_homozstd[[X]][N_SNPs_cor>=10]$NCD1_tf0.3, min(hsinv0102_homozph3[[X]][[2]][N_SNPs_cor>=10]$NCD1_tf0.3, na.rm=T)),
                ncd2_rawp_0.5:= ecdf_fun(NCD2_chr4_homozstd[[X]][IS>=10]$NCD2_tf0.5, min(hsinv0102_homozph3[[X]][[3]][(N_SNPs_cor+N_FDs_cor)>=10]$NCD2_tf0.5, na.rm=T)),
                ncd2_rawp_0.4:= ecdf_fun(NCD2_chr4_homozstd[[X]][IS>=10]$NCD2_tf0.4, min(hsinv0102_homozph3[[X]][[3]][(N_SNPs_cor+N_FDs_cor)>=10]$NCD2_tf0.4, na.rm=T)),
                ncd2_rawp_0.3:= ecdf_fun(NCD2_chr4_homozstd[[X]][IS>=10]$NCD2_tf0.3, min(hsinv0102_homozph3[[X]][[3]][(N_SNPs_cor+N_FDs_cor)>=10]$NCD2_tf0.3, na.rm=T))]
                print(pops[X]);
        }





Store(dt_ncd_3)

#ficou igual o dt_ncd_3...
#calculate NCD1 and NCD2 for entire chr4 in hgins0102
NCD1_chr4_homozstd<- mclapply2(1:26, function(X) NCD1(X=hsinv0102[[X]], W=2000, S=1000))
mclapply(1:26, function(X) NCD1_chr4_homozstd[[X]][N_SNPs_cor>=10])-> NCD1_chr4_homozstd
Store(NCD1_chr_homozstd)

NCD2_chr4_homozstd<- mclapply2(1:26, function(X) NCD2(X=hsinv0102[[X]], FD_list[[4]], W=2000, S=1000))

mclapply(1:26, function(X) NCD2_chr4_homozstd[[X]][,IS:=N_SNPs_cor+N_FDs_cor])-> NCD2_chr_homozstd
Store(NCD1_chr4_homozstd)

#calculate NCD1 and NCD2 with func_inv and add results to 34th position in
#res_all_inv

#and dt_ncd

#do the re-samplings and add to 

all_samps



###################
#make a table

## old stuff ## ## old stuff ## ## old stuff ## ## old stuff ##

#do other populations

#now get the lowest NCD

#collect min NCD values from resamplings

system.time(vec_NCD1_tf0.5_sim<-mclapply2(1:length(test_B_2), function(y) unlist(as.numeric(mclapply2(1:1000, function(x) test_B_2[[y]][[x]][test_B_2[[y]][[x]][, which.min(NCD1_tf0.5)], NCD1_tf0.5]))))) #
names(vec_NCD1_tf0.5_sim)<- names(test_A); Store(vec_NCD1_tf0.5_sim)
system.time(vec_NCD1_tf0.4_sim<-mclapply2(1:length(test_B_2), function(y) unlist(as.numeric(mclapply2(1:1000, function(x) test_B_2[[y]][[x]][test_B_2[[y]][[x]][, which.min(NCD1_tf0.4)], NCD1_tf0.4]))))) # 157 seconds
names(vec_NCD1_tf0.4_sim)<- names(test_A);Store(vec_NCD1_tf0.4_sim)
system.time(vec_NCD1_tf0.3_sim<-mclapply2(1:length(test_B_2), function(y) unlist(as.numeric(mclapply2(1:1000, function(x) test_B_2[[y]][[x]][test_B_2[[y]][[x]][, which.min(NCD1_tf0.3)], NCD1_tf0.3]))))) # 474
names(vec_NCD1_tf0.3_sim)<- names(test_A); Store(vec_NCD1_tf0.3_sim)
#
system.time(vec_NCD1_tf0.5_sim_GBR<-mclapply2(1:length(test_B_GBR), function(y) unlist(as.numeric(mclapply2(1:1000, function(x) test_B_GBR[[y]][[x]][test_B_GBR[[y]][[x]][, which.min(NCD1_tf0.5)], NCD1_tf0.5]))))) #
names(vec_NCD1_tf0.5_sim_GBR)<- names(test_A);Store(vec_NCD1_tf0.5_sim_GBR)
system.time(vec_NCD1_tf0.4_sim_GBR<-mclapply2(1:length(test_B_GBR), function(y) unlist(as.numeric(mclapply2(1:1000, function(x) test_B_GBR[[y]][[x]][test_B_GBR[[y]][[x]][, which.min(NCD1_tf0.4)], NCD1_tf0.4]))))) # 
names(vec_NCD1_tf0.4_sim_GBR)<- names(test_A);Store(vec_NCD1_tf0.4_sim_GBR)
system.time(vec_NCD1_tf0.3_sim_GBR<-mclapply2(1:length(test_B_GBR), function(y) unlist(as.numeric(mclapply2(1:1000, function(x) test_B_GBR[[y]][[x]][test_B_GBR[[y]][[x]][, which.min(NCD1_tf0.3)], NCD1_tf0.3]))))) #
names(vec_NCD1_tf0.3_sim_GBR)<- names(test_A);Store(vec_NCD1_tf0.3_sim_GBR)
#
system.time(vec_NCD2_tf0.5_sim<-mclapply2(1:length(test_D), function(y) unlist(as.numeric(mclapply2(1:1000, function(x) test_D[[y]][[x]][test_D[[y]][[x]][, which.min(NCD2_tf0.5)], NCD2_tf0.5])))))
system.time(vec_NCD2_tf0.4_sim<-mclapply2(1:length(test_D), function(y) unlist(as.numeric(mclapply2(1:1000, function(x) test_D[[y]][[x]][test_D[[y]][[x]][, which.min(NCD2_tf0.4)], NCD2_tf0.4])))))
system.time(vec_NCD2_tf0.3_sim<-mclapply2(1:length(test_D), function(y) unlist(as.numeric(mclapply2(1:1000, function(x) test_D[[y]][[x]][test_D[[y]][[x]][, which.min(NCD2_tf0.3)], NCD2_tf0.3])))))

Store(vec_NCD2_tf0.5_sim, vec_NCD2_tf0.4_sim, vec_NCD2_tf0.3_sim)

# add correct p-value to table using ecdf

#test: works, need to do for the others.

res_ncd1[, cor_p_tf0.5:=c(sapply(1:(length(test_B_2)-2), function(i) ecdf_fun(vec_NCD1_tf0.5_sim[[i]], min(na.omit(test_ncd1[[i]][N_SNPs_cor>=10])$NCD1_tf0.5, na.rm=T))), ecdf_fun(vec_NCD1_tf0.5_sim[[23]], min(na.omit(Inv_Chr4_NCD1[[1]][N_SNPs_cor>=10]$NCD1_tf0.5, na.rm=T))), ecdf_fun(vec_NCD1_tf0.5_sim[[24]],  min(na.omit(Inv_Chr4_PrivStd_NCD1[[1]][N_SNPs_cor>=10]$NCD1_tf0.5, na.rm=T))))]

res_ncd1[, cor_p_tf0.4:=c(sapply(1:(length(test_B_2)-2), function(i) ecdf_fun(vec_NCD1_tf0.4_sim[[i]], min(na.omit(test_ncd1[[i]][N_SNPs_cor>=10])$NCD1_tf0.4, na.rm=T))), ecdf_fun(vec_NCD1_tf0.4_sim[[23]], min(na.omit(Inv_Chr4_NCD1[[1]][N_SNPs_cor>=10]$NCD1_tf0.4, na.rm=T))), ecdf_fun(vec_NCD1_tf0.4_sim[[24]],  min(na.omit(Inv_Chr4_PrivStd_NCD1[[1]][N_SNPs_cor>=10]$NCD1_tf0.4, na.rm=T))))]

res_ncd1[, cor_p_tf0.3:=c(sapply(1:(length(test_B_2)-2), function(i) ecdf_fun(vec_NCD1_tf0.3_sim[[i]], min(na.omit(test_ncd1[[i]][N_SNPs_cor>=10])$NCD1_tf0.3, na.rm=T))), ecdf_fun(vec_NCD1_tf0.3_sim[[23]], min(na.omit(Inv_Chr4_NCD1[[1]][N_SNPs_cor>=10]$NCD1_tf0.3, na.rm=T))), ecdf_fun(vec_NCD1_tf0.3_sim[[24]],  min(na.omit(Inv_Chr4_PrivStd_NCD1[[1]][N_SNPs_cor>=10]$NCD1_tf0.3, na.rm=T))))]

Store(res_ncd1)



#not obsolete: add human-chimp coverage filter

colnames(H_C_cov)<-c('chr','POS1', 'POS2')

 H_C_cov[order(chr,POS1)]-> H_C_cov
 split(H_C_cov, by="chr")-> H_C_cov





