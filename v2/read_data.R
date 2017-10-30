##################################################################################################
#		Barbara D. Bitarello
#		Re-doing analyses for inversions for only seven populations, and including chrX
#
#################################################################################################

library(plyr);
library(dplyr);
library(data.table);
library(parallel);
library(SOAR);Sys.setenv(R_LOCAL_CACHE="inversions")
library(doMC);
registerDoMC(4)
pops<-c("CHB" ,"CEU" ,"LWK" ,"GIH" ,"JPT", "TSI","YRI")
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/mclapply2.R')
ecdf_fun <- function(x,perc) ecdf(x)(perc) #http://stats.stackexchange.com/questions/50080/estimate-quantile-of-value-in-a-vector
source('/mnt/sequencedb/PopGen/barbara/NCD-Statistics/scripts/NCD_func.R')
##########

fread('../input_files/RegionsToAnalyse_45inversions_hg19.txt')-> all_inv
fread('../input_files/VariantsClassified_HsInv0102_using434Ph3Samples.txt')-> var_dt 


#####


FD_list<-factor('list',23)
mclapply2(c(1:22,"X"), function(i)
paste0('zcat < /mnt/sequencedb/PopGen/cesare/bs_genomescan/outgroup_files/fds.hg19_pantro2.', i, '.tsv.gz') %>%
data.table::fread(sep='\t')) -> FD_list

for (i in 1:23){
colnames(FD_list[[i]])<-c('CHR', 'POS', 'REF', 'Chimp_REF')
FD_list[[i]] %>% mutate(ID=paste0(CHR, "|",POS)) %>% as.data.table -> FD_list[[i]]
}

names(FD_list)<-c(paste0("chr",c(seq(1:22),"X")))

mclapply2(1:23, function(i)
#HC_div1<-HC_div %>% dplyr::mutate(REF=toupper(REF), Chimp_REF=toupper(Chimp_REF)) %>% dplyr::filter(POS>=40224426 & POS<= 40247234) %>% as.data.table  #select range of inversions +- 10 kb
FD_list[[i]] %>% dplyr::mutate(REF=toupper(REF), Chimp_REF=toupper(Chimp_REF)) %>% as.data.table) -> FD_list
Store(FD_list)

#############################################
POPS_AF<-vector('list', 7)
names(POPS_AF)<-pops
for (i in pops){
        POPS_AF[[i]]<-fread(paste0('zcat ../v2/',i,'_final.bed.gz'))
        POPS_AF[[i]]<- POPS_AF[[i]][,.(V2,V2,V4,V5,V6)]
        colnames(POPS_AF[[i]])<-c('CHR','POS','nAL', 'N_chr','AF');
	POPS_AF[[i]][,ID:=paste0(CHR,"|",POS)][order(POS)] -> POPS_AF[[i]];
	setkey(POPS_AF[[i]], ID); setDT(unique(POPS_AF[[i]]))-> POPS_AF[[i]]
print(i)
}

Store(POPS_AF)

for (i in pops){
	gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",POPS_AF[[i]]$AF))))-> POPS_AF[[i]]$AF; #clean up and keep only AF
        setDT(POPS_AF[[i]])[, paste0("AF", 1:3) := tstrsplit(AF, ";")]; #split AF into 3 cols
        as.numeric(POPS_AF[[i]]$AF1)-> POPS_AF[[i]]$AF1; as.numeric(POPS_AF[[i]]$AF2)-> POPS_AF[[i]]$AF2; as.numeric(POPS_AF[[i]]$AF3)-> POPS_AF[[i]]$AF3; #make them numeric
        POPS_AF[[i]][,MAF:=pmin(AF1,AF2,AF3, na.rm=T)] -> POPS_AF[[i]];
	POPS_AF[[j]][[i]][MAF<=0.5]-> POPS_AF[[j]][[i]];
print(i)
}

mclapply2(POPS_AF, function(X) X[,.(CHR,POS,ID,REF,ALT,MAF)])-> POPS_AF_v2

remove(POPS_AF); gc()


system.time(saveRDS(POPS_AF_v2,file="POPS_AF_v2.RData")) #has less columns

#end


