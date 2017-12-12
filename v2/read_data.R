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
#####

#can be skipped

Res_Alt_list<-vector('list',23)
for (i in 1:22) {
fread(paste0("gunzip -c /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr", i, ".phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz|grep \"^[^#;]\"|awk \'{print $1,$2,$4,$5}\'"))-> Res_Alt_list[[i]]
print(i)
}

fread(paste0("gunzip -c /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1a.20130502.genotypes.vcf.gz |grep \"^[^#;]\"|awk \'{print $1,$2,$4,$5}\'"))-> Res_Alt_list[[23]]

do.call('rbind', Res_Alt_list)-> Res_Alt_list;
colnames(Res_Alt_list)<-c("CHR","POS","REF","ALT");
Res_Alt_list[order(CHR, POS)]-> Res_Alt_list
Store(Res_Alt_list);

####

FD_list<-factor('list',23)
mclapply2(c(1:22,"X"), function(i)
paste0('zcat < /mnt/sequencedb/PopGen/cesare/bs_genomescan/outgroup_files/fds.hg19_pantro2.', i, '.tsv.gz') %>%
data.table::fread(sep='\t')) -> FD_list

#pantro4 (10.12.2017)
setwd('/mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/v2/')

FD_list_pantro4<-factor('list',23)
mclapply2(c(1:22,"X"), function(i)
paste0('zcat < ', i, ".tab.gz") %>% data.table::fread(sep="\t"))-> FD_list_pantro4

for(i in 1:23){
	colnames(FD_list_pantro4[[i]])<-c('CHR', 'POS', 'REF', 'Chimp_REF', 'panpan1.1', 'gorgor3', 'ponabe2','rhemac3','Vindija33.19','Altai', 'Denisova')
	FD_list_pantro4[[i]]<-FD_list_pantro4[[i]] %>% dplyr::select('CHR', 'POS', 'REF', 'Chimp_REF') %>% as.data.table
	print(i)
	}

#

for(i in 1:23){
	FD_list_pantro4[[i]][Chimp_REF %in% c("A","C","G","T")]-> FD_list_pantro4[[i]] #remove Ns and other weirdness
	FD_list_pantro4[[i]][, ID:=paste0(CHR, "|",POS)] ->FD_list_pantro4[[i]]
	}

names(FD_list_pantro4)<-c(1:22,"X")
for (i in 1:23){
	FD_list_pantro4[[i]][REF!=Chimp_REF]-> FD_list_pantro4[[i]]
	}

Store(FD_list_pantro4)


for (i in 1:23){
colnames(FD_list[[i]])<-c('CHR', 'POS', 'REF', 'Chimp_REF')
FD_list[[i]] %>% mutate(ID=paste0(CHR, "|",POS)) %>% as.data.table -> FD_list[[i]]
}

names(FD_list)<-c(1:22,"X")
Store(FD_list)

mclapply2(1:23, function(i)
#HC_div1<-HC_div %>% dplyr::mutate(REF=toupper(REF), Chimp_REF=toupper(Chimp_REF)) %>% dplyr::filter(POS>=40224426 & POS<= 40247234) %>% as.data.table  #select range of inversions +- 10 kb
FD_list[[i]] %>% dplyr::mutate(REF=toupper(REF), Chimp_REF=toupper(Chimp_REF)) %>% as.data.table) -> FD_list
names(FD_list)<-c(1:22,"X")
Store(FD_list)


Store(FD_list_pantro4)
#interesting, there are more FDs using pantro2 than pantro4
#############################################
POPS_AF<-vector('list', 7)
names(POPS_AF)<-pops
for (i in pops){
        POPS_AF[[i]]<-fread(paste0('zcat ../v2/',i,'_final.bed.gz'))
        POPS_AF[[i]]<- POPS_AF[[i]][,.(V1,V2,V4,V5,V6)]
        colnames(POPS_AF[[i]])<-c('CHR','POS','nAL', 'N_chr','AF');
	POPS_AF[[i]][,CHR:= gsub("chr","",CHR)]
	POPS_AF[[i]][,ID:=paste0(CHR,"|",POS)][order(POS)] -> POPS_AF[[i]];
	setkey(POPS_AF[[i]], ID); setDT(unique(POPS_AF[[i]]))-> POPS_AF[[i]]
print(i)
}

for (i in pops){
	setkey(POPS_AF[[i]], CHR, POS); 
	setkey(Res_Alt_list, CHR, POS);
	POPS_AF[[i]][order(CHR, POS)]-> POPS_AF[[i]]
	Res_Alt_list[order(CHR, POS)]-> Res_Alt_list
        POPS_AF[[i]][Res_Alt_list, nomatch=0] -> POPS_AF[[i]];
print(i)
}

for(i in pops){
	setkey(POPS_AF[[i]], ID); 
	unique(POPS_AF[[i]])-> POPS_AF[[i]];
}

for(i in pops){
        POPS_AF[[i]][-(grep("\\b[A-Z]{2,}:\\b",POPS_AF[[i]]$AF)),][order(POS)] -> POPS_AF[[i]];#exclude lines with indels etc
	gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",POPS_AF[[i]]$AF))))-> POPS_AF[[i]]$AF; #clean up and keep only AF
        setDT(POPS_AF[[i]])[, paste0("AF", 1:3) := tstrsplit(AF, ";")]; #split AF into 3 cols
        as.numeric(POPS_AF[[i]]$AF1)-> POPS_AF[[i]]$AF1; as.numeric(POPS_AF[[i]]$AF2)-> POPS_AF[[i]]$AF2; as.numeric(POPS_AF[[i]]$AF3)-> POPS_AF[[i]]$AF3; #make them numeric
        POPS_AF[[i]][,MAF:=pmin(AF1,AF2,AF3, na.rm=T)] -> POPS_AF[[i]];
	POPS_AF[[i]][MAF<=0.5]-> POPS_AF[[i]];
print(i)
}

mclapply2(POPS_AF, function(X) X[,.(CHR,POS,ID,REF,ALT,MAF)])-> POPS_AF_v2

remove(POPS_AF); gc()


system.time(saveRDS(POPS_AF_v2,file="POPS_AF_v2.RData")) #has less columns

#end
##add great apes

apes<-list()
for(chr in c(1:22,"X")){
		if(chr %in% c(1:22,"X")[1:22]){
	        fread(paste0('zcat ', chr, '.tab.gz'))-> apes[[chr]]
        	} else {
		fread(paste0('zcat ', chr, '.tab.gz'))-> apes[[23]]
		print(chr)
}
}
saveRDS(apes, file='apes.RData')
#stopped here 24.11.

