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

system.time(source('run_NCD1.r')) #

system.time(source('run_NCD2.r')) #

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

mclapply(1:7, function(X) res_all_inv[[X]][-c(12,33)])-> res_all_inv2

Store(res_all_inv)
Store(res_all_inv2)

###################################
###################################

POPS_AF<- readRDS('POPS_AF_v2.RData')

dt_ncd<-vector('list',7)
names(dt_ncd)<-pops

for(J in 1:7){ #attention 
	for (K in names(res_all_inv2[[J]])){
		cbind((res_all_inv2[[J]][[K]][[2]] %>% dplyr::summarise(
		ncd1_rawp_0.5 = ecdf_fun(pops_NCD1_gen_IS[[J]]$NCD1_tf0.5, min(res_all_inv2[[J]][[K]][[2]][N_SNPs_cor>=8]$NCD1_tf0.5, na.rm=T)), 
		ncd1_rawp_0.4 = ecdf_fun(pops_NCD1_gen_IS[[J]]$NCD1_tf0.4, min(res_all_inv2[[J]][[K]][[2]][N_SNPs_cor>=8]$NCD1_tf0.4, na.rm=T)),
		ncd1_rawp_0.3 = ecdf_fun(pops_NCD1_gen_IS[[J]]$NCD1_tf0.3, min(res_all_inv2[[J]][[K]][[2]][N_SNPs_cor>=8]$NCD1_tf0.3, na.rm=T)))),
		res_all_inv2[[J]][[K]][[3]] %>% dplyr::summarise(
		ncd2_rawp_0.5 = ecdf_fun(pops_NCD2_gen_IS[[J]]$NCD2_tf0.5, min(res_all_inv2[[J]][[K]][[3]][IS>=8]$NCD2_tf0.5, na.rm=T)), 
		ncd2_rawp_0.4 = ecdf_fun(pops_NCD2_gen_IS[[J]]$NCD2_tf0.4, min(res_all_inv2[[J]][[K]][[3]][IS>=8]$NCD2_tf0.4, na.rm=T)), 
		ncd2_rawp_0.3 = ecdf_fun(pops_NCD2_gen_IS[[J]]$NCD2_tf0.3, min(res_all_inv2[[J]][[K]][[3]][IS>=8]$NCD2_tf0.3, na.rm=T)), Len=inv_2000bp[Inversion==K,Len], Inversion=K)) %>% 
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
mclapply2(1:7, function(Z) NCD1(X=POPS_AF[[Z]][CHR==4][POS %in% var_dt[Type=='PrivateStd']$POS], W=2000, S=1000))-> HsInv0102_privstd_list[[1]]
mclapply2(1:7, function(Z) NCD2(X=POPS_AF[[Z]][CHR==4][POS %in% var_dt[Type=='PrivateStd']$POS], Y=FD_list[[4]], W=2000, S=1000)[,IS:=N_SNPs_cor+N_FDs_cor])-> HsInv0102_privstd_list[[2]]
#
mclapply2(HsInv0102_privstd_list[[1]], function(Z) Z[,c("Chr", "POS1","POS2") := tstrsplit(Win.ID, "|", fixed=TRUE)])-> HsInv0102_privstd_list[[1]]
mclapply2(HsInv0102_privstd_list[[2]], function(Z) Z[,c("Chr", "POS1","POS2") := tstrsplit(Win.ID, "|", fixed=TRUE)])-> HsInv0102_privstd_list[[2]]
#
mclapply2(HsInv0102_privstd_list[[1]], function(Z) Z[order(as.numeric(Chr), as.numeric(POS1))])-> HsInv0102_privstd_list[[1]]
mclapply2(HsInv0102_privstd_list[[2]], function(Z) Z[order(as.numeric(Chr), as.numeric(POS1))])-> HsInv0102_privstd_list[[2]]
#
Store(HsInv0102_privstd_list)

###
dt_ncd_3<-dt_ncd
for(J in 1:7){
setkey(dt_ncd[[J]], Inversion)
setkey(inv_2000bp, Inversion)
dt_ncd[[J]][inv_2000bp, nomatch=0]-> dt_ncd_3[[J]]
dt_ncd_3[[J]][,i.Len:=NULL]
}

for (J in 1:7){
	rbind(dt_ncd_3[[J]], 
	list(ncd1_rawp_0.5=ecdf_fun(pops_NCD1_gen_IS[[J]]$NCD1_tf0.5, min(HsInv0102_privstd_list[[1]][[J]][N_SNPs_cor>=8]$NCD1_tf0.5, na.rm=T)), 
	ncd1_rawp_0.4=ecdf_fun(pops_NCD1_gen_IS[[J]]$NCD1_tf0.4, min(HsInv0102_privstd_list[[1]][[J]][N_SNPs_cor>=8]$NCD1_tf0.4, na.rm=T)), 
	ncd1_rawp_0.3=ecdf_fun(pops_NCD1_gen_IS[[J]]$NCD1_tf0.3, min(HsInv0102_privstd_list[[1]][[J]][N_SNPs_cor>=8]$NCD1_tf0.3, na.rm=T)),
	ncd2_rawp_0.5=ecdf_fun(pops_NCD2_gen_IS[[J]]$NCD2_tf0.5, min(HsInv0102_privstd_list[[2]][[J]][IS>=8]$NCD2_tf0.5, na.rm=T)), 
	ncd2_rawp_0.4=ecdf_fun(pops_NCD2_gen_IS[[J]]$NCD2_tf0.4, min(HsInv0102_privstd_list[[2]][[J]][IS>=8]$NCD2_tf0.4, na.rm=T)), 
	ncd2_rawp_0.3=ecdf_fun(pops_NCD2_gen_IS[[J]]$NCD2_tf0.3, min(HsInv0102_privstd_list[[2]][[J]][IS>=8]$NCD2_tf0.3, na.rm=T)), 
	Inversion='HsInv0102_PrivStd',POP=pops[J], StartNonRecombining=40235009, EndOutside=40235026, StartInside=40235029, 
	EndInside=40237058, StartOutside=40237061, EndNonRecombining=40237631 ,SizeInside=2030, SizeUsable=2619, CHR=4,  Len=2622, Type='long'))-> dt_ncd_3[[J]]
	print(pops[J])
}

for (J in 1:7){
	rbind(dt_ncd_3[[J]], 
	list(ncd1_rawp_0.5=NA, 
	ncd1_rawp_0.4=NA, 
	ncd1_rawp_0.3=NA,
	ncd2_rawp_0.5=NA, 
	ncd2_rawp_0.4=NA, 
	ncd2_rawp_0.3=NA, 
	Inversion='HsInv0102_homozPh3',
	POP=pops[J], 
	StartNonRecombining=40235009, 
	EndOutside=40235026, 
	StartInside=40235029, 
	EndInside=40237058, 
	StartOutside=40237061,
	EndNonRecombining=40237631 ,SizeInside=2030, SizeUsable=2619, CHR=4,  Len=2622, Type='long'))-> dt_ncd_3[[J]]
	print(pops[J])
}


Store(dt_ncd, dt_ncd_3)
Objects()
#stopped here 05.11
tmp<- mclapply2(c(1:22,"X"), function(x) seq(POPS_AF[[1]][CHR==x][order(as.numeric(CHR))][2000, POS], POPS_AF[[1]][CHR==x][order(as.numeric(CHR))][nrow(POPS_AF[[1]][CHR==x])-2000, POS]), mc.cores=2)
names(tmp)<-c(1:22,"X")
Store(tmp)
Objects()


dt_ncd_3[[1]][,Inversion]-> inversions
Store(inversions)


func_inv2<-function(inv, dt,pop){
type<-dt[Inversion==inv][,Type]
as.numeric(dt[Inversion==inv][,CHR])-> chr
#if(type=='long'){
as.numeric(dt[Inversion==inv][,StartNonRecombining])-> st
as.numeric(dt[Inversion==inv][,EndNonRecombining])-> end
POPS_AF_homoz[[pop]][CHR==chr][POS>=st & POS<=end]-> res
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
####
##
hsinv0102_homozph3<-mclapply2(1:7, function(J) func_inv2(inv="HsInv0102_homozPh3", dt=dt_ncd_3[[J]], pop=pops[J]))  #raw
Store(hsinv0102_homozph3)

#stopped her eon Nov 15. just need to add 37 and 38 to dt ncd 3 and make the final table and edit the report.
for(X in 1:7){
                dt_ncd_3[[X]][37,
                ncd1_rawp_0.5:= ecdf_fun(pops_NCD1_gen_IS[[X]]$NCD1_tf0.5, min(HsInv0102_privstd_list[[1]][[X]][N_SNPs_cor>=8]$NCD1_tf0.5, na.rm=T)),
                ncd1_rawp_0.4:= ecdf_fun(pops_NCD1_gen_IS[[X]]$NCD1_tf0.4, min(HsInv0102_privstd_list[[1]][[X]][N_SNPs_cor>=8]$NCD1_tf0.4, na.rm=T)),
                ncd1_rawp_0.3:= ecdf_fun(pops_NCD1_gen_IS[[X]]$NCD1_tf0.3, min(HsInv0102_privstd_list[[1]][[X]][N_SNPs_cor>=8]$NCD1_tf0.3, na.rm=T)),
                ncd2_rawp_0.5:= ecdf_fun(pops_NCD2_gen_IS[[X]]$NCD2_tf0.5, min(HsInv0102_privstd_list[[2]][[X]][(N_SNPs_cor+N_FDs_cor)>=8]$NCD2_tf0.5, na.rm=T)),
                ncd2_rawp_0.4:= ecdf_fun(pops_NCD2_gen_IS[[X]]$NCD2_tf0.4, min(HsInv0102_privstd_list[[2]][[X]][(N_SNPs_cor+N_FDs_cor)>=8]$NCD2_tf0.4, na.rm=T)),
                ncd2_rawp_0.3:= ecdf_fun(pops_NCD2_gen_IS[[X]]$NCD2_tf0.3, min(HsInv0102_privstd_list[[2]][[X]][(N_SNPs_cor+N_FDs_cor)>=8]$NCD2_tf0.3, na.rm=T))]
                print(pops[X]);
        }
for(X in 1:7){
                dt_ncd_3[[X]][38,
                ncd1_rawp_0.5:= ecdf_fun(pops_NCD1_gen_IS_homoz[[X]][N_SNPs_cor>=8]$NCD1_tf0.5, min(hsinv0102_homozph3[[X]][[2]][N_SNPs_cor>=8]$NCD1_tf0.5, na.rm=T)),
                ncd1_rawp_0.4:= ecdf_fun(pops_NCD1_gen_IS_homoz[[X]][N_SNPs_cor>=8]$NCD1_tf0.4, min(hsinv0102_homozph3[[X]][[2]][N_SNPs_cor>=8]$NCD1_tf0.4, na.rm=T)),
                ncd1_rawp_0.3:= ecdf_fun(pops_NCD1_gen_IS_homoz[[X]][N_SNPs_cor>=8]$NCD1_tf0.3, min(hsinv0102_homozph3[[X]][[2]][N_SNPs_cor>=8]$NCD1_tf0.3, na.rm=T)),
                ncd2_rawp_0.5:= ecdf_fun(pops_NCD2_gen_IS_homoz[[X]][N_SNPs_cor+N_FDs_cor>=8]$NCD2_tf0.5, min(hsinv0102_homozph3[[X]][[3]][(N_SNPs_cor+N_FDs_cor)>=8]$NCD2_tf0.5, na.rm=T)),
                ncd2_rawp_0.4:= ecdf_fun(pops_NCD2_gen_IS_homoz[[X]][N_SNPs_cor+N_FDs_cor>=8]$NCD2_tf0.4, min(hsinv0102_homozph3[[X]][[3]][(N_SNPs_cor+N_FDs_cor)>=8]$NCD2_tf0.4, na.rm=T)),
                ncd2_rawp_0.3:= ecdf_fun(pops_NCD2_gen_IS_homoz[[X]][N_SNPs_cor+N_FDs_cor>=8]$NCD2_tf0.3, min(hsinv0102_homozph3[[X]][[3]][(N_SNPs_cor+N_FDs_cor)>=8]$NCD2_tf0.3, na.rm=T))]
                print(pops[X]);
        }

Store(dt_ncd_3)

#######################################
#######################################
func_inv3<-function(inv, dt, pop, rep){
type<-dt[[pop]][Inversion==inv][,Type]
len<- dt[[pop]][Inversion==inv][,Len]
dt[[pop]][Inversion==inv][,CHR]-> chr
res<-data.table(CHR=chr, Central_Pos=sample(tmp[[chr]], 1))[,POS:=Central_Pos-len/2][,POS2:=Central_Pos+ len/2]
res2<- POPS_AF_sub[[pop]][CHR==chr][POS>= res$POS & POS<= res$POS2]
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
######################################
######################################

POPS_AF[c(1,2,3,4)]-> POPS_AF_1_2_3_4
mem_change(remove(POPS_AF))
Store(POPS_AF_1_2_3_4)
POPS_AF_1_2_3_4-> POPS_AF_sub
head(tmp[[1]])
samp<-vector('list', 4)
samp[[1]]<- vector('list', 38)
samp[[2]]<- vector('list', 38)
samp[[3]]<- vector('list', 38)
samp[[4]]<- vector('list', 38)
system.time(for(J in 1:4){  #RUNNING BIONC06
		for (K in 1:37){ #later need to do hmozph3
			samp[[J]][[K]]<-foreach(X=1:1000) %dopar%
			try(func_inv3(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=X))
		}
		gc();
		print(pops[J]);
		print('done\n');
	}
) #63900.99 


Store(samp) #
Objects()
lapply(1:4, function(Y) sapply(1:37, function(Z) sum(sapply(samp[[Y]][[Z]], function(X) length(X)==17))))
#still some missing
system.time(for(J in 1:4){#repeat until full
	for(K in 1:37){
		1000-sum(sapply(samp[[J]][[K]], function(X) length(X)==17))->a
		for(i in a){
			if(i>0){
				samp[[J]][[K]]<-append(samp[[J]][[K]],
				foreach(X=1:i) %dopar%
				try(func_inv3(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=1000+X)))
			}
		}
	}	
})

#bion09
#############################################
#############################################
POPS_AF[c(5,6,7)]-> POPS_AF_5_6_7
mem_change(remove(POPS_AF))
Store(POPS_AF_5_6_7)
POPS_AF_5_6_7-> POPS_AF_sub
head(tmp[[1]])
samp1<-vector('list', 3)
samp1[[1]]<- vector('list', 38)
samp1[[2]]<- vector('list', 38)
samp1[[3]]<- vector('list', 38)

system.time(for(J in 5:7){  #RUNNING BIONC06
for (K in 1:37){ #later need to do hmozph3
        samp1[[J-4]][[K]]<-foreach(X=1:1000) %dopar%
        try(func_inv3(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=X))
}
gc();
print(pops[J]);
print('done\n');}
) #37107.81

Store(samp1)

Objects()
lapply(1:3, function(Y) sapply(1:37, function(Z) sum(sapply(samp1[[Y]][[Z]], function(X) length(X)==17))))
# filling missing ones
system.time(for(J in 5:7){ #run this until it fills up
        for(K in 1:37){
                1000-sum(sapply(samp1[[J-4]][[K]], function(X) length(X)==17))->a
                for(i in a){
                        if(i>0){
                                samp1[[J-4]][[K]]<-append(samp1[[J-4]][[K]],
                                foreach(X=1:i) %dopar%
                                try(func_inv3(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=1000+X)))
                        }
                }
        }
}) #8665.621 #900

Store(samp1)

#remove empty elements
samp<-mclapply2(1:4, function(Z) lapply(1:37, function(K) samp[[Z]][[K]][which(sapply(samp[[Z]][[K]], function(X)  length(X)==17))]))
names(samp)<- pops[1:4]
for(J in 1:4){
names(samp[[J]])<- inversions[1:37]
}
Store(samp)
#
samp1<-mclapply2(1:3, function(Z) lapply(1:37, function(K) samp1[[Z]][[K]][which(sapply(samp1[[Z]][[K]], function(X)  length(X)==17))]))
names(samp1)<- pops[5:7]
for(J in 5:7){
names(samp1[[J-4]])<- inversions[1:37]
}
Store(samp1)

####
#now for homoz ph3 individuals for hsinv0102

POPS_AF_homoz-> POPS_AF_sub
head(tmp[[1]])
samp3<-vector('list', 7)

system.time(for(J in 1:7){  #RUNNING BIONC06
		K<-38
        	samp3[[J]]<-foreach(X=1:3000) %dopar%
        	try(func_inv3(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=X))
		
	print(pops[J]);
	print('done\n');
	}
) #37107.81

Store(samp3)

Objects()

# filling missing ones
system.time(for(J in 1:7){ #run this until it fills up
        	K<-38
                1000-sum(unlist(sapply(samp3[[J]], function(X) ncol(X)==17)))-> i 
		if(i>0){
                                samp3[[J]]<-append(samp3[[J]],
                                foreach(X=1:i) %dopar%
                                try(func_inv3(dt=dt_ncd_3, inv=inversions[K], pop=pops[J], rep=3000+X)))
                        }
                }

		             )	

Store(samp3)
mclapply2(1:7, function(Z) samp3[[Z]][which(sapply(samp3[[Z]], function(X)  length(X)==17))])-> samp4
names(samp4)<-pops

Store(samp4)
#########################

#combine all pops

test<-append(samp,samp1)

mclapply2(1:length(test), function(X) list.append(test[[X]], samp4[[X]]))-> test2

test2<- mclapply2(test2, function(X) mclapply2(X, function(Y) Y[1:1000]))

for(J in 1:7){
	names(test2[[J]])[[38]]<-inversions[38]
	}

for (J in 1:length(test2)){

lapply(1:38, function(K) lapply(1:1000, function(X) test2[[J]][[K]][[X]][, Rep:=X]))-> test2[[J]]
}

remove(test)

gc()

test3<-lapply(test2, function(X) lapply(X, function(Y) do.call(rbind, Y)))


for(J in length(test3)){

	for(K in 1:38){
		unique(test3[[J]][[K]]$Rep)-> tp

		for(i in 1:1000){
			test3[[J]][[K]][Rep==tp[i]][,Rep:= i]
		}
	}
}

remove(test2)


test4<- lapply(1:length(test3), function(J) do.call(rbind, test3[[J]])) 


test5<-do.call(rbind,test4)

setkey(test5, Win.ID)

for(i in 1:7){
setkey(pops_NCD2_gen_IS_chimpcov[[i]], Win.ID)
}


lapply(1:7, function(X) test5[POP==pops[X]][pops_NCD2_gen_IS_chimpcov[[X]]])-> BLA
BLA<-do.call(rbind, BLA)
setkey(BLA, Win.ID, POP, Rep, Inversion)
unique(BLA)-> BLA
BLA[COV2>=500/3000]-> test5

system.time(test5  %>% dplyr::filter(N_SNPs_cor.x>=8) %>% dplyr::group_by(POP,Inversion,Rep) %>% 
dplyr::summarise(min_ncd1_tf0.5=min(NCD1_tf0.5), min_ncd1_tf0.4=min(NCD1_tf0.4), min_ncd1_tf0.3=min(NCD1_tf0.3)) %>%
as.data.table-> test6
)

system.time(test5  %>% dplyr::filter(N_SNPs_cor.x+N_FDs_cor>=8) %>% dplyr::group_by(POP,Inversion,Rep) %>%
dplyr::summarise(min_ncd2_tf0.5=min(NCD2_tf0.5), min_ncd2_tf0.4=min(NCD2_tf0.4), min_ncd2_tf0.3=min(NCD2_tf0.3)) %>%
as.data.table -> test6b)

#merge(test6b, test6, all.x=T)-> test7

dt_ncd_4<-do.call(rbind, dt_ncd_3)
dt_ncd_4 %>% dplyr::select(ncd1_rawp_0.5:POP, StartNonRecombining:Type) %>% as.data.table-> dt_ncd_4
setkey(dt_ncd_4, POP, Inversion)
#


Store(test6, test6b)		
cmp<-vector('list',7)
names(cmp)<-pops

for(i in 1:7){
cmp[[i]]<-vector('list', 38) 
names(cmp[[i]])<-inversions[1:38]
}


for(i in 1:7){
	 for(j in 1:length(res_all_inv2[[1]])){
	 setkey(res_all_inv2[[i]][[j]][[3]], Win.ID)
}
}

#stopped here on 29.11. need to find a way to check whether inversion windows have chimp cov or not.

for (J in pops){
	for(i in inversions[1:36]){
		if(is.data.table(res_all_inv2[[J]][[i]][[2]]) & length(res_all_inv2[[J]][[i]][[2]][N_SNPs_cor>=8][, NCD1_tf0.5])>=1){
			data.table(Inversion=i,
			ncd1_corp_0.5=ecdf_fun(test6[POP==J & Inversion==i][,min_ncd1_tf0.5], min(res_all_inv2[[J]][[i]][[2]][N_SNPs_cor>=8][, NCD1_tf0.5], na.rm=T)),
			ncd1_corp_0.4=ecdf_fun(test6[POP==J & Inversion==i][,min_ncd1_tf0.4], min(res_all_inv2[[J]][[i]][[2]][N_SNPs_cor>=8][, NCD1_tf0.4], na.rm=T)),
			ncd1_corp_0.3=ecdf_fun(test6[POP==J & Inversion==i][,min_ncd1_tf0.3], min(res_all_inv2[[J]][[i]][[2]][N_SNPs_cor>=8][, NCD1_tf0.3], na.rm=T)), POP=J)-> cmp[[J]][[i]]
			print(i)
		} else if(length(res_all_inv2[[J]][[i]][[2]][N_SNPs_cor>=8][, NCD1_tf0.5])==0){

			data.table(Inversion=i,
                        ncd1_corp_0.5=NA,
                        ncd1_corp_0.4=NA,
                        ncd1_corp_0.3=NA, POP=J)-> cmp[[J]][[i]]
		} else {
		 	data.table(Inversion=i,
                        ncd1_corp_0.5=NA,
                        ncd1_corp_0.4=NA,
                        ncd1_corp_0.3=NA, POP=J)-> cmp[[J]][[i]]
		
		}
	}
}

#names(HsInv0102_privstd_list[[1]])<-pops
#names(HsInv0102_privstd_list[[2]])<-pops
#Store(HsInv0102_privstd_list) 
#names(hsinv0102_homozph3)<-pops
#Store(hsinv0102_homozph3)

for(J in pops){
	i<-inversions[37]
	 data.table(Inversion=i,
                        ncd1_corp_0.5=ecdf_fun(test6[POP==J & Inversion==i][,min_ncd1_tf0.5], min(HsInv0102_privstd_list[[1]][[J]][N_SNPs_cor>=8][, NCD1_tf0.5], na.rm=T)),
                        ncd1_corp_0.4=ecdf_fun(test6[POP==J & Inversion==i][,min_ncd1_tf0.4], min(HsInv0102_privstd_list[[1]][[J]][N_SNPs_cor>=8][, NCD1_tf0.4], na.rm=T)),
                        ncd1_corp_0.3=ecdf_fun(test6[POP==J & Inversion==i][,min_ncd1_tf0.3], min(HsInv0102_privstd_list[[1]][[J]][N_SNPs_cor>=8][, NCD1_tf0.3], na.rm=T)), POP=J)-> cmp[[J]][[i]]
}


for (J in pops){
	i<-inversions[38]
	data.table(Inversion=i,
                       ncd1_corp_0.5=ecdf_fun(test6[POP==J & Inversion==i][,min_ncd1_tf0.5], min(hsinv0102_homozph3[[J]][[2]][N_SNPs_cor>=8][, NCD1_tf0.5], na.rm=T)),
                       ncd1_corp_0.4=ecdf_fun(test6[POP==J & Inversion==i][,min_ncd1_tf0.4], min(hsinv0102_homozph3[[J]][[2]][N_SNPs_cor>=8][, NCD1_tf0.4], na.rm=T)),
                       ncd1_corp_0.3=ecdf_fun(test6[POP==J & Inversion==i][,min_ncd1_tf0.3], min(hsinv0102_homozph3[[J]][[2]][N_SNPs_cor>=8][, NCD1_tf0.3], na.rm=T)), POP=J)-> cmp[[J]][[i]]
}


mclapply(cmp, function(X) do.call(rbind, X))-> YAY
remove(cmp)

cmp2<-vector('list',7)
names(cmp2)<-pops

for(i in 1:7){
	cmp2[[i]]<-vector('list', 38)
	names(cmp2[[i]])<-inversions[1:38]
}

for (J in pops){
        for(i in inversions[1:36]){
                if(is.data.table(res_all_inv2[[J]][[i]][[3]]) & length(res_all_inv2[[J]][[i]][[3]][IS>=8][, NCD2_tf0.5])>=1){
                        data.table(Inversion=i,
                        ncd2_corp_0.5=ecdf_fun(test6b[POP==J & Inversion==i][,min_ncd2_tf0.5], min(res_all_inv2[[J]][[i]][[3]][IS>=8][, NCD2_tf0.5], na.rm=T)),
                        ncd2_corp_0.4=ecdf_fun(test6b[POP==J & Inversion==i][,min_ncd2_tf0.4], min(res_all_inv2[[J]][[i]][[3]][IS>=8][, NCD2_tf0.4], na.rm=T)),
                        ncd2_corp_0.3=ecdf_fun(test6b[POP==J & Inversion==i][,min_ncd2_tf0.3], min(res_all_inv2[[J]][[i]][[3]][IS>=8][, NCD2_tf0.3], na.rm=T)), POP=J)-> cmp2[[J]][[i]]
                        print(i)
                } else if(length(res_all_inv2[[J]][[i]][[3]][IS>=8][, NCD2_tf0.5])==0){
                        data.table(Inversion=i,
                        ncd2_corp_0.5=NA,
                        ncd2_corp_0.4=NA,
                        ncd2_corp_0.3=NA, POP=J)-> cmp2[[J]][[i]]
                } else {
                        data.table(Inversion=i,
                        ncd2_corp_0.5=NA,
                        ncd2_corp_0.4=NA,
                        ncd2_corp_0.3=NA, POP=J)-> cmp2[[J]][[i]]

                }
        }
}



for(J in pops){
        i<-inversions[37]
         data.table(Inversion=i,
                        ncd2_corp_0.5=ecdf_fun(test6b[POP==J & Inversion==i][,min_ncd2_tf0.5], min(HsInv0102_privstd_list[[2]][[J]][N_SNPs_cor+N_FDs_cor>=8][, NCD2_tf0.5], na.rm=T)),
                        ncd2_corp_0.4=ecdf_fun(test6b[POP==J & Inversion==i][,min_ncd2_tf0.4], min(HsInv0102_privstd_list[[2]][[J]][N_SNPs_cor+N_FDs_cor>=8][, NCD2_tf0.4], na.rm=T)),
                        ncd2_corp_0.3=ecdf_fun(test6b[POP==J & Inversion==i][,min_ncd2_tf0.3], min(HsInv0102_privstd_list[[2]][[J]][N_SNPs_cor+N_FDs_cor>=8][, NCD2_tf0.3], na.rm=T)), POP=J)-> cmp2[[J]][[i]]
}


for (J in pops){
        i<-inversions[38]
        data.table(Inversion=i,
                       ncd2_corp_0.5=ecdf_fun(test6b[POP==J & Inversion==i][,min_ncd2_tf0.5], min(hsinv0102_homozph3[[J]][[3]][N_SNPs_cor+N_FDs_cor>=8][, NCD2_tf0.5], na.rm=T)),
                       ncd2_corp_0.4=ecdf_fun(test6b[POP==J & Inversion==i][,min_ncd2_tf0.4], min(hsinv0102_homozph3[[J]][[3]][N_SNPs_cor+N_FDs_cor>=8][, NCD2_tf0.4], na.rm=T)),
                       ncd2_corp_0.3=ecdf_fun(test6b[POP==J & Inversion==i][,min_ncd2_tf0.3], min(hsinv0102_homozph3[[J]][[3]][N_SNPs_cor+N_FDs_cor>=8][, NCD2_tf0.3], na.rm=T)), POP=J)-> cmp2[[J]][[i]]
}

mclapply(cmp2, function(X) do.call(rbind, X))-> YAY2
remove(cmp2)

do.call(rbind,YAY)-> YAY
do.call(rbind,YAY2)-> YAY2

setkey(YAY, Inversion, POP)
setkey(YAY2, Inversion, POP)
dt_ncd_4[YAY][YAY2]-> final

final[, Cont:='BLA']
for(i in 1:nrow(final)){
	if(final[i,POP]==pops[3]|final[i,POP]==pops[7]){
		final[i, Cont:='AFRICA']
	} else if(final[i,POP]==pops[2]|final[i,POP]==pops[6]){
		final[i, Cont:='EUROPE']
	} else{
		final[i,Cont:='ASIA']
	}
}	

final %>% select(Inversion,CHR, StartNonRecombining:SizeUsable, Type, ncd1_rawp_0.5:ncd2_rawp_0.3, ncd1_corp_0.5:ncd2_corp_0.3, POP, Cont) %>% as.data.table-> final
write.table(final, file='final.txt', col.names=T, row.names=F, sep="\t", quote=F)

Store(final)

fread('zcat /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/v2/hg19.pantro2.sort.bed.gz')-> H_C_cov

colnames(H_C_cov)<-c('CHR','POS', 'POS2')
H_C_cov[,CHR:= gsub("chr","",CHR)]
setkey(H_C_cov, CHR, POS, POS2)

H_C_cov[order(chr,POS1)]-> H_C_cov


bed<-pops_NCD2_gen_IS[[1]]
bed[,CHR:=Chr]
bed[,POS:=POS1]
bed %>% select(Win.ID:IS, CHR, POS, POS2) %>% mutate(POS=as.numeric(POS), POS2=as.numeric(POS2)) %>%  as.data.table -> bed



#testing
apes<-c()
#test<-vector('list', 23)
#names(test)<-c(1:22,"X")
	
chimp_cov<-function(chr=22, w=2000){ # cutoff=500/3000){ #chr, window size used in the scan, chimp_cov cutoff desired (this default value is based on the ncd paper)
fread(paste0('zcat ', chr, '.tab.gz'))-> apes	
colnames(apes)<-c('CHR', 'POS', 'hg19', 'pantro4', 'panpan1.1', 'gorgor3', 'ponabe2','rhemac3','Vindija33.19','Altai', 'Denisova')
apes %>% select(hg19:Denisova, CHR, POS) %>% as.data.table -> apes
apes[, POS2:=POS]
apes %>%  mutate(CHR=as.character(CHR)) %>% as.data.table -> apes
bed1<-bed[CHR==chr]
setkey(bed1, POS, POS2)
setkey(apes, POS, POS2)
foverlaps(apes, bed1, type='within')-> test
setDT(test)[is.na(Win.ID)==F][,.(.N), by=.(Win.ID)]-> testA;
setDT(test)[is.na(Win.ID)==F][pantro4 !="N" & pantro4 != "-"][,.(.N), by=.(Win.ID)]-> testB;
setkey(testA, Win.ID);
setkey(testB, Win.ID);		
testA[testB][,COV1:=N/w][, COV2:=i.N/w]-> res
#stats<-nrow(res[COV2>=cutoff])/nrow(res)
remove(apes)
remove(bed1)
gc()
return(list(res=res))
}
#

chimp_cov2<-function(chr=22, w=2000, bed=POPS_AF[[1]]){ # cutoff=500/3000){ #chr, window size used in the scan, chimp_cov cutoff desired (this default value is based on the ncd paper)
fread(paste0('zcat ', chr, '.tab.gz'))-> apes
colnames(apes)<-c('CHR', 'POS', 'hg19', 'pantro4', 'panpan1.1', 'gorgor3', 'ponabe2','rhemac3','Vindija33.19','Altai', 'Denisova')
apes %>% select(hg19:Denisova, CHR, POS) %>% as.data.table -> apes
apes[, POS2:=POS]
apes %>%  mutate(CHR=as.character(CHR)) %>% as.data.table -> apes
bed1<-bed[CHR==chr]
setkey(bed1, POS, POS2)
setkey(apes, POS, POS2)
foverlaps(apes, bed1, type='within')-> test
setDT(test)[is.na(Win.ID)==F][,.(.N), by=.(Win.ID)]-> testA;
setDT(test)[is.na(Win.ID)==F][pantro4 !="N" & pantro4 != "-"][,.(.N), by=.(Win.ID)]-> testB;
setkey(testA, Win.ID);
setkey(testB, Win.ID);
testA[testB][,COV1:=N/w][, COV2:=i.N/w]-> res
#stats<-nrow(res[COV2>=cutoff])/nrow(res)
remove(apes)
remove(bed1)
gc()
return(list(res=res))
}


test<-foreach (i=c(15:22,"X"), .combine = c)  %dopar% 
	chimp_cov(chr=i)
do.call(rbind, test)-> B.test

remove(test)
test<-foreach (i=c(6:14), .combine = c)  %dopar%  
        chimp_cov(chr=i)
	
do.call(rbind, test)-> C.test

remove(test)

test<-foreach (i=c(3:5), .combine = c)  %dopar%      
        chimp_cov(chr=i)

do.call(rbind, test)-> D.test

remove(test)

gc()

test.2<-chimp_cov(chr=2)
rbind(test.2[[1]], D.test, C.test, B.test)-> TMP
Store(TMP)
test.1<-chimp_cov(chr=1)
rbind(test.1[[1]],TMP)-> TMP2
Store(TMP2)

#add this info to scan

pops_NCD2_gen_IS_chimpcov<-list()

for(i in 1:7){
setkey(pops_NCD2_gen_IS[[i]], Win.ID)
setkey(TMP2, Win.ID)

pops_NCD2_gen_IS[[i]][TMP2]-> pops_NCD2_gen_IS_chimpcov[[i]]
}


##############################
#all o the remaining code is about homozph3 and has been ran already.

source('ph3_homoz.R ')
##########################
##########################
####### The End ##########
##########################
##########################

