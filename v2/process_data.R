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
library(dplyr)

pops<-c("CHB" ,"CEU" ,"LWK" ,"GIH" ,"JPT", "TSI","YRI")

#system.time(source('run_NCD1.r')) #

#system.time(source('run_NCD2.r')) #

#################
#stopped here 03.11
fread('../input_files/RegionsToAnalyse_45inversions_hg19.txt')-> all_inv
fread('../input_files/VariantsClassified_HsInv0102_using434Ph3Samples.txt')-> var_dt

all_inv[,CHR:= Chromosome]
all_inv[,Chromosome:=NULL]

all_inv[CHR!="Y"]-> all_inv_notY
all_inv_notY[,.(.N),.(CHR)]

all_inv_notY[order(as.numeric(CHR))]-> all_inv_notY #check this before continuing.

a<-as.numeric(unique(all_inv_notY$CHR)) #check this before continuing.does not work because of chr X

all_inv_notY[,Len:=as.numeric(EndNonRecombining)-as.numeric(StartNonRecombining)]
na.omit(all_inv_notY)[, Type:='long']-> all_inv_notY_b
all_inv_notY[is.na(Len)==T][,Len:= SizeInside][,Type:='short']-> all_inv_notY_a
rbind(all_inv_notY_a, all_inv_notY_b)-> all_inv_notY

all_inv_notY[Len>=2000][order(as.numeric(CHR))]->inv_2000bp #to relate to the scan #38 inversions

inv_2000bp[,Inversion]-> tmp_all_inv
####
#make a function
func_inv<-function(inv, dt,pop){
type<-dt[Inversion==inv][,Type]
dt[Inversion==inv][,CHR]-> chr
if(type=='long'){
	as.numeric(dt[Inversion==inv][,StartNonRecombining])-> st
	as.numeric(dt[Inversion==inv][,EndNonRecombining])-> end
	POPS_AF[[J]][CHR==chr][POS>=st & POS<=end]-> res
	tmp<-unique(c(res[POS>= dt[Inversion==inv,EndOutside] & POS<=dt[Inversion==inv,StartInside]], res[POS>=dt[Inversion==inv,EndInside] & POS<=dt[Inversion==inv,StartOutside]]))
	if(length(tmp)!=0){
		res[!(POS %in% tmp)]-> res
	}
	else{
	res<-res
	}	
} else if(type=='short'){
	as.numeric(dt[Inversion==inv][,StartInside])-> st
	as.numeric(dt[Inversion==inv][,EndInside])-> end
	POPS_AF[[J]][CHR==chr][POS>=st & POS<=end]-> res
	}
if(nrow(res)>0){
	NCD1(X=res, W=2000, S=1000)-> ncd1_res
	NCD2(X=res, Y=FD_list[[chr]],  W=2000, S=1000)-> ncd2_res
	return(list(res, ncd1_res, ncd2_res))
	} else if(nrow(FD_list[[chr]][POS>=st & POS<=end])>0){
	NCD2(X=res, Y=FD_list[[chr]],  W=2000, S=1000)-> ncd2_res
	return(list(res="No SNPs in this range", ncd1_res="Cannot be calculated", ncd2_res))
	} else{
	return(list(res="No SNPs in this range", ncd1_res="Cannot be calculated", ncd2_res="Cannot be calculated"))
}
}

#####

res_all_inv<-vector('list',7)
system.time(for(J in 1:7){
mclapply2(1:length(tmp_all_inv), function(Z) func_inv(tmp_all_inv[Z], inv_2000bp, pop=pops[J]))-> res_all_inv[[J]]
print(paste0(pops[J],'done\n'))
}) #588.366  #ok, so out of these 38 inversions, two could not have NCD calculated:

inv_2000bp[c(12,33),]

names(res_all_inv)<-pops

for(J in 1:7){
names(res_all_inv[[J]])<- tmp_all_inv
}
for(J in 1:7){
res_all_inv[[J]][inv_2000bp[,Inversion]]
}

for(J in 1:7){
for (K in 1:11){
res_all_inv[[J]][[K]][[3]][,IS:=N_SNPs_cor+N_FDs_cor]-> res_all_inv[[J]][[K]][[3]]
}
for (K in 13:32){
res_all_inv[[J]][[K]][[3]][,IS:=N_SNPs_cor+N_FDs_cor]-> res_all_inv[[J]][[K]][[3]]
}
for (K in 34:38){
res_all_inv[[J]][[K]][[3]][,IS:=N_SNPs_cor+N_FDs_cor]-> res_all_inv[[J]][[K]][[3]]
}}

Store(res_all_inv)

###################################
###################################

dt_ncd<-vector('list',7)
names(dt_ncd)<-pops
for(J in 1:7){ #attention 
	for (K in names(res_all_inv[[J]])[c(1:11,13:32,34:38)]){
		cbind((res_all_inv[[J]][[K]][[2]] %>% dplyr::summarise(
		ncd1_rawp_0.5=ecdf_fun(pops_NCD1_gen_IS[[J]]$NCD1_tf0.5, min(res_all_inv[[J]][[K]][[2]][N_SNPs_cor>=8]$NCD1_tf0.5, na.rm=T)), 
		ncd1_rawp_0.4= ecdf_fun(pops_NCD1_gen_IS[[J]]$NCD1_tf0.4, min(res_all_inv[[J]][[K]][[2]][N_SNPs_cor>=8]$NCD1_tf0.4, na.rm=T)),
		ncd1_rawp_0.3= ecdf_fun(pops_NCD1_gen_IS[[J]]$NCD1_tf0.3, min(res_all_inv[[J]][[K]][[2]][N_SNPs_cor>=8]$NCD1_tf0.3, na.rm=T)))),res_all_inv[[J]][[K]][[3]] %>% dplyr::summarise(
		ncd2_rawp_0.5= ecdf_fun(pops_NCD2_gen_IS[[J]]$NCD2_tf0.5, min(res_all_inv[[J]][[K]][[3]][IS>=8]$NCD2_tf0.5, na.rm=T)), 
		ncd2_rawp_0.4 = ecdf_fun(pops_NCD2_gen_IS[[J]]$NCD2_tf0.4, min(res_all_inv[[J]][[K]][[3]][IS>=8]$NCD2_tf0.4, na.rm=T)), 
		ncd2_rawp_0.3 = ecdf_fun(pops_NCD2_gen_IS[[J]]$NCD2_tf0.3, min(res_all_inv[[J]][[K]][[3]][IS>=8]$NCD2_tf0.3, na.rm=T)), Len=inv_2000bp[Inversion==K,Len], Inversion=K)) %>% 
		as.data.table-> dt_ncd[[J]][[K]];
		print(K)
	}
	do.call(rbind, dt_ncd[[J]])-> dt_ncd[[J]]
	dt_ncd[[J]][, POP:=pops[J]]
	setDT(dt_ncd[[J]])
	print(pops[J])
}

names(dt_ncd)<-pops
Store(dt_ncd)


#priv std
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
Store(HsInv0102_privstd_list)

###

#stopped here 03.11
for (J in 1:7){
	rbind(dt_ncd[[J]], 
	list(ncd1_rawp_0.5=ecdf_fun(pops_NCD1_gen_IS[[J]]$NCD1_tf0.5, min(HsInv0102_privstd_list[[1]][[J]][N_SNPs_cor>=8]$NCD1_tf0.5, na.rm=T)), ncd1_rawp_0.4=ecdf_fun(pops_NCD1_gen_IS[[J]]$NCD1_tf0.4, min(HsInv0102_privstd_list[[1]][[J]][N_SNPs_cor>=8]$NCD1_tf0.4, na.rm=T)), ncd1_rawp_0.3=ecdf_fun(pops_NCD1_gen_IS[[J]]$NCD1_tf0.3, min(HsInv0102_privstd_list[[1]][[J]][N_SNPs_cor>=8]$NCD1_tf0.3, na.rm=T)),ncd2_rawp_0.5=ecdf_fun(pops_NCD2_gen_IS[[J]]$NCD2_tf0.5, min(HsInv0102_privstd_list[[2]][[J][IS>=8]$NCD2_tf0.5, na.rm=T)), ncd2_rawp_0.4=ecdf_fun(pops_NCD2_gen_IS[[J]]$NCD2_tf0.4, min(HsInv0102_privstd_list[[2]][[J]][IS>=8]$NCD2_tf0.4, na.rm=T)), ncd2_rawp_0.3=ecdf_fun(pops_NCD2_gen_IS[[J]]$NCD2_tf0.3, min(HsInv0102_privstd_list[[2]][[J]][IS>=8]$NCD2_tf0.3, na.rm=T)), Inversion='HsInv0102_PrivStd',POP=pops[J], StartNonRecombining=40235009, EndOutside=40235026, StartInside=40235029, EndInside=40237058, StartOutside=40237061, EndNonRecombining=40237631 ,SizeInside=2030, SizeUsable=2619, CHR=4,  Len=2622, Type='long'))-> dt_ncd_3[[J]]
	print(pops[J])
}

for (J in 1:26){
rbind(dt_ncd_3[[J]], list(ncd1_rawp_0.5=NA, ncd1_rawp_0.4=NA, ncd1_rawp_0.3=NA,ncd2_rawp_0.5=NA, ncd2_rawp_0.4=NA, ncd2_rawp_0.3=NA, Inversion='HsInv0102_homozPh3',POP=pops[J], StartNonRecombining=40235009, EndOutside=40235026, StartInside=40235029, EndInside=40237058, StartOutside=40237061, EndNonRecombining=40237631 ,SizeInside=2030, SizeUsable=2619, CHR=4,  Len=2622, Type='long'))-> dt_ncd_3[[J]]
print(pops[J])}


