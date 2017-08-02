#Last modified: 24.08.2017
library(pegas);library(dplyr)
library(plyr);library(data.table)
library(parallel);library(lattice)
library(SOAR);Sys.setenv(R_LOCAL_CACHE="inversions")
library(ggplot2);library(splitstackshape)
library(pryr)
#biocLite(VariantAnnotation)
#library(vcfR);
library(doMC)
registerDoMC(10);library(bigmemory);
#source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/bedtools_inR.R')
#source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/SFS_script.r')
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/mclapply2.R')
#Note: there is room for improvement of this script if I use more data.table.
####################################
####################################
####################################
Objects()

mclapply2(1:22, function(i)
fread(paste0("gunzip -c /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr", i, ".phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz|grep \"^[^#;]\"|awk \'{print $1,$2,$4,$5}\'")))-> Res_Alt_list

for(i in 1:22){colnames(Res_Alt_list[[i]])<-c("CHR","POS","REF","ALT")};

Store(Res_Alt_list);

###
###

pops<-c("ACB","ASW","BEB","CDX","CEU","CHB","CHS","CLM","ESN","FIN","GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL","PEL","PJL","PUR","STU","TSI","YRI")

POPS_AF<-vector('list', 26);
tmp<-vector('list',22);
#for(a in 1:26){
#POPS_AF[[a]]<- tmp
#}

for(i in 1:26){
print (pops[i]);
mclapply2(1:22, function(x) fread(paste0('input_files/', pops[[i]], '_chr',x, '_AF_2.frq')))-> POPS_AF[[i]]#create list for freq files
}

for(j in 1:26){
        for (i in 1:22){
                colnames(POPS_AF[[j]][[i]])<-c('CHR','POS','nAL', 'N_chr','AF'); #add colnames
        }
}

#remove duplicated positions...
for(j in 1:26){
                for (i in 1:22){
                        POPS_AF[[j]][[i]] %>% dplyr::mutate(ID=paste0(CHR,"|",POS)) %>% arrange(POS) %>% as.data.table-> POPS_AF[[j]][[i]];
                        setkey(POPS_AF[[j]][[i]], ID); unique(POPS_AF[[j]][[i]])-> POPS_AF[[j]][[i]];
                        POPS_AF[[j]][[i]] %>% dplyr::left_join(Res_Alt_list[[i]]) %>% arrange(POS) %>% as.data.table -> POPS_AF[[j]][[i]];
                        setkey(POPS_AF[[j]][[i]], ID); unique(POPS_AF[[j]][[i]])-> POPS_AF[[j]][[i]];
                        POPS_AF[[j]][[i]][-(grep("\\b[A-Z]{2,}:\\b",POPS_AF[[j]][[i]]$AF)),] %>% arrange(POS) %>% as.data.table -> POPS_AF[[j]][[i]];#exclude lines with indels etc
                        print(paste0('chr', i));
                };
        gc(); print (pops[j]);
        }

# allele frequencies:
for (j in 1:26){
                for (i in 1:22){
#               Pops_AF[[j]][[i]][-(grep("\\b[A-Z]{2,}:\\b",Pops_AF[[j]][[i]]$AF)),] %>% arrange(POS) %>% as.data.table -> Pops_AF[[j]][[i]];#exclude lines with indels etc
                        gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",POPS_AF[[j]][[i]]$AF))))-> POPS_AF[[j]][[i]]$AF; #clean up and keep only AF
                        setDT(POPS_AF[[j]][[i]])[, paste0("AF", 1:3) := tstrsplit(AF, ";")]; #split AF into 3 cols
                        as.numeric(POPS_AF[[j]][[i]]$AF1)-> POPS_AF[[j]][[i]]$AF1; as.numeric(POPS_AF[[j]][[i]]$AF2)-> POPS_AF[[j]][[i]]$AF2; as.numeric(POPS_AF[[j]][[i]]$AF3)-> POPS_AF[[j]][[i]]$AF3; #make them numeric
                        POPS_AF[[j]][[i]] %>% dplyr::mutate(MAF=pmin(AF1,AF2,AF3, na.rm=T)) %>%  as.data.table -> POPS_AF[[j]][[i]];
                        print(paste0('chr',i));
                };
        gc(); print (pops[j]);
        }
#there are still some indels left, so

for(j in 1:26){
	for (i in 1:22){
		filter(POPS_AF[[j]][[i]], MAF<=0.5) %>% as.data.table-> POPS_AF[[j]][[i]];
		print(paste0('chr', i));
}
}
system.time(saveRDS(POPS_AF,file="POPS_AF_v2.RData"))

