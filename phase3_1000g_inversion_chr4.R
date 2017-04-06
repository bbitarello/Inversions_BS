##########################################################################################
#	Use Phase 1 1000G data to verify allele frequencies within inversion in chr4
#	Bárbara Bitarello
#	Created: 14.02.2017
#	Last modified: 04.04.2017
##########################################################################################

#preamble
library(pegas)
library(dplyr)
library(plyr)
library(data.table)
library(parallel)
library(lattice)
library(SOAR)
Sys.setenv(R_LOCAL_CACHE="inversions")
library(ggplot2)
library(splitstackshape)
#biocLite(VariantAnnotation)
library(vcfR)
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/bedtools_inR.R')
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/SFS_script.r')
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/mclapply2.R')
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/SFS_script.r')

#load('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/read_scan_data/list.SCAN.3.RData')
####
#read in data 

#Inversions file:
fread('input_files/VariantsClassified_HsInv0102_using434Ph3Samples.txt')-> var_dt

#human-chimp divergence file:
'zcat < /mnt/sequencedb/PopGen/barbara/NCV_dir_package/input_data/outgroup_files/fds.chr4.hg19_pantro2.Map50_100.TRF.SDs.bed.gz' %>%
data.table::fread(sep='\t') -> HC_div
colnames(HC_div)<-c('CHROM', 'POS', 'REF', 'Chimp_REF')

#use uppercase
HC_div1<-HC_div %>% dplyr::mutate(REF=toupper(REF), Chimp_REF=toupper(Chimp_REF)) %>% dplyr::filter(POS>=40224426 & POS<= 40247234) %>% as.data.table  #select range of inversions +- 10 kb
HC_div2<-HC_div %>% dplyr::mutate(REF=toupper(REF), Chimp_REF=toupper(Chimp_REF)) %>% as.data.table  #all FDs for chr4

#######################################################################
##################   SFS SFS SFS SFS SFS SFS SFS ######################
#######################################################################
# SKIP to LINE 115
#TO DO: SO FAR I DID THIS ONLY FOR THE INVERSION IN CHR4. COULD CONSIDER DOING IT FOR ALL (37).

#see in README. how we recoded the chr4 phase 3 vcf file so that ALT==DERIVED
#DAF
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --freq --out chr4_1000gph3_DAF')
#now each pop separately:
#DAF
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/LWK_samples.txt --freq --out input_files/LWK_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/YRI_samples.txt --freq --out input_files/YRI_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/TSI_samples.txt --freq --out input_files/TSI_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/GBR_samples.txt --freq --out input_files/GBR_chr4_DAF')

#read in the DAFs in order to plot SFS

fread('input_files/LWK_chr4_DAF.frq', header=F)->LWK_chr4_DAF
fread('input_files/YRI_chr4_DAF.frq', header=F)->YRI_chr4_DAF
fread('input_files/GBR_chr4_DAF.frq', header=F)->GBR_chr4_DAF
fread('input_files/TSI_chr4_DAF.frq', header=F)->TSI_chr4_DAF

colnames(LWK_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(YRI_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(GBR_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(TSI_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')

as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",LWK_chr4_DAF$DAF)))))-> LWK_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",YRI_chr4_DAF$DAF)))))-> YRI_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",GBR_chr4_DAF$DAF)))))-> GBR_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",TSI_chr4_DAF$DAF)))))-> TSI_chr4_DAF$DAF
#plot SFS

pdf('figures/test.SFS.all_chr4.pdf')
#par(mfrow=c(4,1))
barplot(table(factor(round(LWK_chr4_DAF$DAF,2),levels=seq(0,1,0.01)))/sum(table(factor(round(LWK_chr4_DAF$DAF,2),levels=seq(0,1,0.01)))), main="LWK all SNPs in chr4", xlab="DAF",ylim=c(0,1), cex.names=0.7) #rounding the DAF makes the plot look nicer
barplot(table(factor(round(YRI_chr4_DAF$DAF,2),levels=seq(0,1,0.01)))/sum(table(factor(round(YRI_chr4_DAF$DAF,2),levels=seq(0,1,0.01)))), main="YRI all SNPs in chr4", xlab="DAF",ylim=c(0,1), cex.names=0.7)
barplot(table(factor(round(GBR_chr4_DAF$DAF,2),levels=seq(0,1,0.01)))/sum(table(factor(round(GBR_chr4_DAF$DAF,2),levels=seq(0,1,0.01)))), main="GBR all SNPs in chr4", xlab="DAF",ylim=c(0,1), cex.names=0.7)
barplot(table(factor(round(TSI_chr4_DAF$DAF,2),levels=seq(0,1,0.01)))/sum(table(factor(round(TSI_chr4_DAF$DAF,2),levels=seq(0,1,0.01)))), main="TSI all SNPs in chr4", xlab="DAF",ylim=c(0,1), cex.names=0.7)
dev.off()

as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(LWK_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> LWK_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(YRI_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> YRI_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(GBR_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> GBR_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(TSI_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> TSI_inv_DAF

#now only PrivateStd SNps in the inversion
dplyr::filter(var_dt, Type=="PrivateStd")$POS-> Std
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(LWK_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> LWK_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(YRI_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> YRI_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(GBR_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> GBR_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(TSI_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> TSI_invPrivStd_DAF

pdf('figures/test.SFS.allinv_chr4.pdf')
#par(mfrow=c(4,1))
barplot(table(factor(round(LWK_inv_DAF,2),levels=seq(0,1,0.01)))/sum(table(factor(round(LWK_inv_DAF,2),levels=seq(0,1,0.01)))), main="LWK all SNPs in Inversion", xlab="DAF",ylim=c(0,1),cex.names=0.7)
barplot(table(factor(round(YRI_inv_DAF,2),levels=seq(0,1,0.01)))/sum(table(factor(round(YRI_inv_DAF,2),levels=seq(0,1,0.01)))), main="YRI all SNPs in Inversion", xlab="DAF",ylim=c(0,1),cex.names=0.7)
barplot(table(factor(round(GBR_inv_DAF,2),levels=seq(0,1,0.01)))/sum(table(factor(round(GBR_inv_DAF,2),levels=seq(0,1,0.01)))), main="GBR all SNPs in Inversion", xlab="DAF",ylim=c(0,1),cex.names=0.7)
barplot(table(factor(round(TSI_inv_DAF,2),levels=seq(0,1,0.01)))/sum(table(factor(round(TSI_inv_DAF,2),levels=seq(0,1,0.01)))), main="TSI all SNPs in Inversion", xlab="DAF",ylim=c(0,1),cex.names=0.7)
dev.off()

pdf('figures/test.SFS.PrivStdinv_chr4.pdf')
#par(mfrow=c(4,1))
barplot(table(factor(round(LWK_invPrivStd_DAF,2),levels=seq(0,1,0.01)))/sum(table(factor(round(LWK_invPrivStd_DAF,2),levels=seq(0,1,0.01)))), main="LWK PrivateStd SNPs in Inversion", xlab="DAF",ylim=c(0,1), cex.names=0.7)
barplot(table(factor(round(YRI_invPrivStd_DAF,2),levels=seq(0,1,0.01)))/sum(table(factor(round(YRI_invPrivStd_DAF,2),levels=seq(0,1,0.01)))), main="YRI PrivateStd SNPs in Inversion", xlab="DAF",ylim=c(0,1), cex.names=0.7)
barplot(table(factor(round(GBR_invPrivStd_DAF,2),levels=seq(0,1,0.01)))/sum(table(factor(round(GBR_invPrivStd_DAF,2),levels=seq(0,1,0.01)))), main="GBR PrivateStd SNPs in Inversion", xlab="DAF",ylim=c(0,1), cex.names=0.7)
barplot(table(factor(round(TSI_invPrivStd_DAF,2),levels=seq(0,1,0.01)))/sum(table(factor(round(TSI_invPrivStd_DAF,2),levels=seq(0,1,0.01)))), main="TSI PrivateStd SNPs in Inversion", xlab="DAF",ylim=c(0,1), cex.names=0.7)
dev.off()

#Store

Store(LWK_chr4_DAF, YRI_chr4_DAF, GBR_chr4_DAF, TSI_chr4_DAF, LWK_inv_DAF, YRI_inv_DAF, GBR_inv_DAF, TSI_inv_DAF, LWK_invPrivStd_DAF, YRI_invPrivStd_DAF, GBR_invPrivStd_DAF, TSI_invPrivStd_DAF)

##################################################################################
### END OF SFS #### END OF SFS #### END OF SFS #### END OF SFS ### END OF SFS ####
##################################################################################
##################################################################################
### END OF SFS #### END OF SFS #### END OF SFS #### END OF SFS ### END OF SFS ####
##################################################################################
#################################################################################
### END OF SFS #### END OF SFS #### END OF SFS #### END OF SFS ### END OF SFS ####
##################################################################################
##################################################################################
### END OF SFS #### END OF SFS #### END OF SFS #### END OF SFS ### END OF SFS ####
##################################################################################
#
##################################################################################
############## PtoD ### PtoD #### PToD ### PtoD #### PToD ########################
##################################################################################
#the following was done wihtout removing indels. Might need to redo properly.
#skip to line 155
system.time(mclapply2(1:22, function(i) system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/LWK_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/LWK_chr', i, '_AF'))))

system.time(mclapply2(1:22, function(i) system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/GBR_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/GBR_chr', i, '_AF'))))

system.time(mclapply2(1:22, function(i) system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/YRI_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/YRI_chr', i, '_AF'))))

system.time(mclapply2(1:22, function(i) system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/TSI_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/TSI_chr', i, '_AF'))))  #10369.160 

mclapply2(1:22, function(i)
system(paste0('sed \'s/\t/  /\'  input_files/LWK_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/LWK_chr', i, '_AF_2.frq')))
mclapply2(1:22, function(i)
system(paste0('sed \'s/\t/  /\'  input_files/YRI_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/YRI_chr', i, '_AF_2.frq')))
mclapply2(1:22, function(i)
system(paste0('sed \'s/\t/  /\'  input_files/GBR_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/GBR_chr', i, '_AF_2.frq')))
mclapply2(1:22, function(i)
system(paste0('sed \'s/\t/  /\'  input_files/TSI_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/TSI_chr', i, '_AF_2.frq')))

#remove duplicated files
system('rm input_files/*AF.frq')
system('rm input_files/*AF.log')
#####
#
mclapply2(1:22, function(i)
fread(paste0("gunzip -c /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr", i, ".phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz|grep \"^[^#;]\"|awk \'{print $1,$2,$4,$5}\'")))-> Res_Alt_list

for(i in 1:22){colnames(Res_Alt_list[[i]])<-c("CHR","POS","REF","ALT")};

Store(Res_Alt_list); Objects();
#
tmp<-vector('list',22);Pops_AF<-list(LWK=tmp, YRI=tmp, GBR=tmp, TSI=tmp); #create list for freq files
mclapply2(1:22, function(i) fread(paste0('input_files/LWK_chr',i, '_AF_2.frq')))->Pops_AF[['LWK']];
mclapply2(1:22, function(i) fread(paste0('input_files/YRI_chr',i, '_AF_2.frq')))->Pops_AF[['YRI']];
mclapply2(1:22, function(i) fread(paste0('input_files/GBR_chr',i, '_AF_2.frq')))->Pops_AF[['GBR']];
mclapply2(1:22, function(i) fread(paste0('input_files/TSI_chr',i, '_AF_2.frq')))->Pops_AF[['TSI']];
#
for(j in 1:4){for (i in 1:22){colnames(Pops_AF[[j]][[i]])<-c('CHR','POS','nAL', 'N_chr','AF')}}; #add colnames

# clean up, remove indels, add MAF colum:
for(j in 1:4){
	for (i in 1:22){
		Pops_AF[[j]][[i]] %>% dplyr::inner_join(Res_Alt_list[[i]]) %>% as.data.table -> Pops_AF[[j]][[i]];
		Pops_AF[[j]][[i]] %>% dplyr::mutate(ID=paste0(CHR,"|",POS)) %>% arrange(POS) %>% as.data.table-> Pops_AF[[j]][[i]];
		setkey(Pops_AF[[j]][[i]], ID); unique(Pops_AF[[j]][[i]]);
	}
	print (j);
}

for (j in 1:4){
	for (i in 1:22){
		Pops_AF[[j]][[i]][-(grep("\\b[A-Z]{2,}:\\b",Pops_AF[[j]][[i]]$AF)),] %>% arrange(POS) %>% as.data.table -> Pops_AF[[j]][[i]];#exclude lines with indels etc
		gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",Pops_AF[[j]][[i]]$AF))))-> Pops_AF[[j]][[i]]$AF; #clean up and keep only AF
		setDT(Pops_AF[[j]][[i]])[, paste0("AF", 1:3) := tstrsplit(AF, ";")]; #split AF into 3 cols
		as.numeric(Pops_AF[[j]][[i]]$AF1)-> Pops_AF[[j]][[i]]$AF1; as.numeric(Pops_AF[[j]][[i]]$AF2)-> Pops_AF[[j]][[i]]$AF2; as.numeric(Pops_AF[[j]][[i]]$AF3)-> Pops_AF[[j]][[i]]$AF3; #make them numeric
		Pops_AF[[j]][[i]] %>% dplyr::mutate(MAF=pmin(AF1,AF2,AF3, na.rm=T)) %>%  as.data.table -> Pops_AF[[j]][[i]];
		}
	print (j);
	}

#remove MAF=0 and MAF=1
for(j in 1:4){
	for (i in 1:22){
		Pops_AF[[j]][[i]] %>% dplyr::filter(!(MAF==0|MAF==1)) %>% as.data.table -> Pops_AF[[j]][[i]];
	}
}

#### PLayground with NCD functions #################################
source('NCD_func.R')

system.time(NCD1(X=LWK_chr4_AF_3, W=3000, S=1500)->out) #running NCD1  
gc()
system.time(NCD1(X=YRI_chr4_AF_3, W=3000, S=1500)->out_YRI) #1422.759 
gc()
system.time(NCD1(X=GBR_chr4_AF_3, W=3000, S=1500)->out_GBR)
gc()
system.time(NCD1(X=TSI_chr4_AF_3, W=3000, S=1500)->out_TSI)  # 10130.651 
gc()

#NCD1 within inversion compared to chromosome 4....(i tried smaller, 1 kb windows, but turns out I had several with few SNPs, so I am going abck to 3 kb)
filter(var_dt, Type=="PrivateStd")$POS->a
system.time(NCD1(X=LWK_chr4_AF_3 %>% dplyr::filter(POS %in% var_dt$POS) %>% as.data.table, W=3000, S=1500)-> out_LWK_inv)
system.time(NCD1(X=LWK_chr4_AF_3 %>% dplyr::filter(POS %in% a) %>% as.data.table, W=3000, S=1500)-> out_LWK_inv_PrivStd)

ecdf_fun <- function(x,perc) ecdf(x)(perc) #http://stats.stackexchange.com/questions/50080/estimate-quantile-of-value-in-a-vector

ecdf_fun(out$tf0.5, min(out_LWK_inv$tf0.5))[1]*100 #1.14% lower tail

ecdf_fun(out$tf0.5, min(out_LWK_inv_PrivStd$tf0.5))[1]*100  #0.5920426 %

ecdf_fun(out$tf0.4, min(out_LWK_inv$tf0.4))[1]*100 # 1.744926 %

ecdf_fun(out$tf0.4, min(out_LWK_inv_PrivStd$tf0.4))[1]*100 # 0.6064437 %

ecdf_fun(out$tf0.3, min(out_LWK_inv$tf0.3))[1]*100 # 6.184445 %

ecdf_fun(out$tf0.3, min(out_LWK_inv_PrivStd$tf0.3))[1]*100 # 0.9528686 %
#inversion

w<-40247234-40224426 #22208 bp, pretty big...


### perhaps obsolete (check later)
#system('sed \'s/\t/  /\'  freq_chr4_inv.txt.frq.count  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/\'|sed \'s/\t/;/\'|sed \'s/  /\t/g\' > freq_inv.txt')

#fread('freq_inv.txt')-> freq_inv #global frequencies of the SNPs in the inversion. If I want pop specific it will require some further filtering.
#colnames(freq_inv)[5]<- gsub(":","_",gsub("\\{|\\}","", colnames(freq_inv)[5])) #fix col name
#index vcf

#system('cat header2.txt test2.vcf.recode.vcf > real_vcf')

#fread('real_vcf', header=T)-> chr4_inv_vcf #read in vcf

#adding a header to this data.table
#colnames(chr4_inv_vcf)<-unlist(strsplit(gsub("#","",system('tail -n 1 header.txt', intern=T)), "\t"))
#chr4_simp<- select(chr4_inv_vcf, CHROM:FORMAT)
#join var_dt to frequency info
#left_join(freq_inv, var_dt) %>% as.data.table -> comp_dt
#now use chimp info:

nrow(HC_div) #209 FD in the inversion +- 10 kb
nrow(comp_dt) # 577 SNPs

comp_dt %>% dplyr::group_by(Type) %>% dplyr::summarise(N=n(), PtoD=N/nrow(HC_div))  #PtoD without correcting FDs

(na.omit(left_join(comp_dt, HC_div) %>% as.data.table) %>% dplyr::filter(Chimp_REF==ALT))$POS-> temp #FDs that are not FDs.9 FDs actually have the chimp state segreagating as alternate allele

filter(HC_div, !(POS %in% temp)) #180 FDs

#this works fine
comp_dt %>% dplyr::group_by(Type) %>% dplyr::summarise(N=n(), D=209, Dcor=180, PtoD=N/D, PtoD_cor=N/Dcor, PtoD_all=nrow(var_dt)/209, PtoD_all_cor=nrow(var_dt)/180)  #PtoD without correcting FDs


#first, use vcf tools to filter the region

#then read in that vcf file
#something new

#i got this from http://grch37.ensembl.org/Homo_sapiens/Tools/AlleleFrequency

#fread('calculated_fra.40225087-40246984.ALL.txt', sep="\t")-> ens_inv
#colnames(ens_inv)[1]<-'CHROM'

#add ancestral info

apply(select(chr4_inv_vcf, INFO), 2, function(x) strsplit(x, ";"))[[1]]-> test
ANCESTRAL<-unlist(mclapply(test, function(x) gsub("\\|","",gsub("AA=","",x[[11]]))))

chr4_inv_vcf %>% mutate(ANC=toupper(ANCESTRAL)) %>% as.data.table -> yay
#i love dplyr
ens_inv %>% dplyr::filter(POS %in% chr4_inv_vcf$POS) %>% mutate(ANC=ANCESTRAL) %>% as.data.table -> yay2
#aso far the 57 cases where ALT is NOT derived, we need to redefine frequencies

yay2 %>% dplyr::mutate(RefAnc=ifelse(REF==ANC, 'TRUE', 'FALSE')) %>% as.data.table %>% filter(RefAnc=='FALSE') %>% as.data.table-> yay3

#now I need to fix the frequencies for these 54 lines.
####################################################################
### END OF OBSOLETE #### END OF OBSOLETE #### END OF OBSOLETE #####
###################################################################

#PRoject: consider extending all of the above to the other inversions as well...
