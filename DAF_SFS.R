library(parallel)
library(doMC)
registerDoMC(11)
library(data.table)
library(dplyr)

#######################################################################
##################   SFS SFS SFS SFS SFS SFS SFS ######################
#######################################################################
#TO DO: SO FAR I DID THIS ONLY FOR THE INVERSION IN CHR4. COULD CONSIDER DOING IT FOR ALL (37).

#see in README. how we recoded the chr4 phase 3 vcf file so that ALT==DERIVED
#DAF
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --freq --out chr4_1000gph3_DAF')
#now each pop separately:
#DAF
#read in the DAFs in order to plot SFS
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/ACB_samples.txt --freq --out input_files/ACB_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/ASW_samples.txt --freq --out input_files/ASW_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/BEB_samples.txt --freq --out input_files/BEB_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/CDX_samples.txt --freq --out input_files/CDX_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/CEU_samples.txt --freq --out input_files/CEU_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/CHB_samples.txt --freq --out input_files/CHB_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/CHS_samples.txt --freq --out input_files/CHS_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/CLM_samples.txt --freq --out input_files/CLM_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/ESN_samples.txt --freq --out input_files/ESN_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/FIN_samples.txt --freq --out input_files/FIN_chr4_DAF')

system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/GBR_samples.txt --freq --out input_files/GBR_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/GIH_samples.txt --freq --out input_files/GIH_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/GWD_samples.txt --freq --out input_files/GWD_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/IBS_samples.txt --freq --out input_files/IBS_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/ITU_samples.txt --freq --out input_files/ITU_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/JPT_samples.txt --freq --out input_files/JPT_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/KHV_samples.txt --freq --out input_files/KHV_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/LWK_samples.txt --freq --out input_files/LWK_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/MSL_samples.txt --freq --out input_files/MSL_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/MXL_samples.txt --freq --out input_files/MXL_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/PEL_samples.txt --freq --out input_files/PEL_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/PJL_samples.txt --freq --out input_files/PJL_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/PUR_samples.txt --freq --out input_files/PUR_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/STU_samples.txt --freq --out input_files/STU_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/TSI_samples.txt --freq --out input_files/TSI_chr4_DAF')
system('vcftools --gzvcf chr4_1000gph3_DAF.vcf.gz --keep input_files/YRI_samples.txt --freq --out input_files/YRI_chr4_DAF')


fread('input_files/ACB_chr4_DAF.frq', header=F)->ACB_chr4_DAF
fread('input_files/ASW_chr4_DAF.frq', header=F)->ASW_chr4_DAF
fread('input_files/BEB_chr4_DAF.frq', header=F)->BEB_chr4_DAF
fread('input_files/CDX_chr4_DAF.frq', header=F)->CDX_chr4_DAF
fread('input_files/CEU_chr4_DAF.frq', header=F)->CEU_chr4_DAF
fread('input_files/CHB_chr4_DAF.frq', header=F)->CHB_chr4_DAF
fread('input_files/CHS_chr4_DAF.frq', header=F)->CHS_chr4_DAF
fread('input_files/CLM_chr4_DAF.frq', header=F)->CLM_chr4_DAF
fread('input_files/ESN_chr4_DAF.frq', header=F)->ESN_chr4_DAF
fread('input_files/FIN_chr4_DAF.frq', header=F)->FIN_chr4_DAF
fread('input_files/GBR_chr4_DAF.frq', header=F)->GBR_chr4_DAF
fread('input_files/GIH_chr4_DAF.frq', header=F)->GIH_chr4_DAF
fread('input_files/GWD_chr4_DAF.frq', header=F)->GWD_chr4_DAF
fread('input_files/IBS_chr4_DAF.frq', header=F)->IBS_chr4_DAF
fread('input_files/ITU_chr4_DAF.frq', header=F)->ITU_chr4_DAF
fread('input_files/JPT_chr4_DAF.frq', header=F)->JPT_chr4_DAF
fread('input_files/KHV_chr4_DAF.frq', header=F)->KHV_chr4_DAF
fread('input_files/LWK_chr4_DAF.frq', header=F)->LWK_chr4_DAF
fread('input_files/MSL_chr4_DAF.frq', header=F)->MSL_chr4_DAF
fread('input_files/MXL_chr4_DAF.frq', header=F)->MXL_chr4_DAF
fread('input_files/PEL_chr4_DAF.frq', header=F)->PEL_chr4_DAF
fread('input_files/PJL_chr4_DAF.frq', header=F)->PJL_chr4_DAF
fread('input_files/PUR_chr4_DAF.frq', header=F)->PUR_chr4_DAF
fread('input_files/STU_chr4_DAF.frq', header=F)->STU_chr4_DAF
fread('input_files/TSI_chr4_DAF.frq', header=F)->TSI_chr4_DAF
fread('input_files/YRI_chr4_DAF.frq', header=F)->YRI_chr4_DAF


colnames(ACB_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(ASW_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(BEB_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(CDX_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(CEU_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(CHB_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(CHS_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(CLM_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(ESN_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(FIN_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(GBR_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(GIH_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(GWD_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(IBS_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(ITU_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(JPT_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(KHV_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(MSL_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(MXL_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(PEL_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(PJL_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(PUR_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(STU_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(TSI_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')
colnames(YRI_chr4_DAF)<-c('CHR','POS','nAL', 'N_chr','AAF','DAF')

as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",ACB_chr4_DAF$DAF)))))-> ACB_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",ASW_chr4_DAF$DAF)))))-> ASW_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",BEB_chr4_DAF$DAF)))))-> BEB_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",CDX_chr4_DAF$DAF)))))-> CDX_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",CEU_chr4_DAF$DAF)))))-> CEU_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",CHB_chr4_DAF$DAF)))))-> CHB_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",CHS_chr4_DAF$DAF)))))-> CHS_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",CLM_chr4_DAF$DAF)))))-> CLM_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",ESN_chr4_DAF$DAF)))))-> ESN_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",FIN_chr4_DAF$DAF)))))-> FIN_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",GBR_chr4_DAF$DAF)))))-> GBR_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",GIH_chr4_DAF$DAF)))))-> GIH_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",GWD_chr4_DAF$DAF)))))-> GWD_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",IBS_chr4_DAF$DAF)))))-> IBS_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",ITU_chr4_DAF$DAF)))))-> ITU_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",JPT_chr4_DAF$DAF)))))-> JPT_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",KHV_chr4_DAF$DAF)))))-> KHV_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",LWK_chr4_DAF$DAF)))))-> LWK_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",MSL_chr4_DAF$DAF)))))-> MSL_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",MXL_chr4_DAF$DAF)))))-> MXL_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",PEL_chr4_DAF$DAF)))))-> PEL_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",PJL_chr4_DAF$DAF)))))-> PJL_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",PUR_chr4_DAF$DAF)))))-> PUR_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",STU_chr4_DAF$DAF)))))-> STU_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",TSI_chr4_DAF$DAF)))))-> TSI_chr4_DAF$DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",YRI_chr4_DAF$DAF)))))-> YRI_chr4_DAF$DAF

pdf('figures/test.SFS.all_chr4.pdf')
#par(mfrow=c(4,1))
barplot(table(factor(round(LWK_chr4_DAF$DAF,2),levels=seq(0,1,0.01)))/sum(table(factor(round(LWK_chr4_DAF$DAF,2),levels=seq(0,1,0.01)))), main="LWK all SNPs in chr4", xlab="DAF",ylim=c(0,1), cex.names=0.7) #rounding the DAF makes the plot look nicer
barplot(table(factor(round(YRI_chr4_DAF$DAF,2),levels=seq(0,1,0.01)))/sum(table(factor(round(YRI_chr4_DAF$DAF,2),levels=seq(0,1,0.01)))), main="YRI all SNPs in chr4", xlab="DAF",ylim=c(0,1), cex.names=0.7)
barplot(table(factor(round(GBR_chr4_DAF$DAF,2),levels=seq(0,1,0.01)))/sum(table(factor(round(GBR_chr4_DAF$DAF,2),levels=seq(0,1,0.01)))), main="GBR all SNPs in chr4", xlab="DAF",ylim=c(0,1), cex.names=0.7)
barplot(table(factor(round(TSI_chr4_DAF$DAF,2),levels=seq(0,1,0.01)))/sum(table(factor(round(TSI_chr4_DAF$DAF,2),levels=seq(0,1,0.01)))), main="TSI all SNPs in chr4", xlab="DAF",ylim=c(0,1), cex.names=0.7)
dev.off()



as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(ACB_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> ACB_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(ASW_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> ASW_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(BEB_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> BEB_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(CDX_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> CDX_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(CEU_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> CEU_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(CHB_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> CHB_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(CHS_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> CHS_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(CLM_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> CLM_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(ESN_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> ESN_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(FIN_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> FIN_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(GBR_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> GBR_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(GIH_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> GIH_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(GWD_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> GWD_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(IBS_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> IBS_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(ITU_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> ITU_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(JPT_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> JPT_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(KHV_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> KHV_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(LWK_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> LWK_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(MSL_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> MSL_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(MXL_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> MXL_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(PEL_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> PEL_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(PJL_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> PJL_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(PUR_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> PUR_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(STU_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> STU_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(TSI_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> TSI_inv_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(YRI_chr4_DAF %>% dplyr::filter(POS %in% var_dt$POS))$DAF)))))-> YRI_inv_DAF

#now only PrivateStd SNps in the inversion
dplyr::filter(var_dt, Type=="PrivateStd")$POS-> Std


as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(ACB_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> ACB_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(ASW_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> ASW_invPrivStd_DAF



as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(BEB_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> BEB_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(CDX_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> CDX_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(CEU_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> CEU_invPrivStd_DAF




as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(CHB_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> CHB_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(CHS_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> CHS_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(CLM_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> CLM_invPrivStd_DAF


as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(ESN_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> ESN_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(FIN_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> FIN_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(GBR_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> GBR_invPrivStd_DAF


as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(GIH_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> GIH_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(GWD_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> GWD_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(IBS_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> IBS_invPrivStd_DAF


as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(ITU_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> ITU_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(JPT_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> JPT_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(KHV_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> KHV_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(LWK_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> LWK_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(MSL_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> MSL_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(MXL_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> MXL_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(PEL_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> PEL_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(PJL_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> PJL_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(PUR_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> PUR_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(STU_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> STU_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(TSI_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> TSI_invPrivStd_DAF
as.numeric(gsub("T:","",gsub("G:","",gsub("C:","", gsub("A:","",(YRI_chr4_DAF %>% dplyr::filter(POS %in% Std))$DAF)))))-> YRI_invPrivStd_DAF

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

ev.off()

#Store

Store(LWK_chr4_DAF, YRI_chr4_DAF, GBR_chr4_DAF, TSI_chr4_DAF, LWK_inv_DAF, YRI_inv_DAF, GBR_inv_DAF, TSI_inv_DAF, LWK_invPrivStd_DAF, YRI_invPrivStd_DAF, GBR_invPrivStd_DAF, TSI_invPrivStd_DAF,
ACB_ch4_DAF, ACB_inv_DAF, ACB_invPrivStd_DAF,
ASW_ch4_DAF, ASW_inv_DAF, ASW_invPrivStd_DAF,
BEB_ch4_DAF, BEB_inv_DAF, BEB_invPrivStd_DAF,
CDX_ch4_DAF, CDX_inv_DAF, CDX_invPrivStd_DAF,
CEU_ch4_DAF, CEU_inv_DAF, CEU_invPrivStd_DAF,
CHB_ch4_DAF, CHB_inv_DAF, CHB_invPrivStd_DAF,
CHS_ch4_DAF, CHS_inv_DAF, CHS_invPrivStd_DAF,
CLM_ch4_DAF, CLM_inv_DAF, CLM_invPrivStd_DAF,
ESN_ch4_DAF, ESN_inv_DAF, ESN_invPrivStd_DAF,
FIN_ch4_DAF, FIN_inv_DAF, FIN_invPrivStd_DAF,
GIH_ch4_DAF, GIH_inv_DAF, GIH_invPrivStd_DAF,
GWD_ch4_DAF, GWD_inv_DAF, GWD_invPrivStd_DAF,
IBS_ch4_DAF, IBS_inv_DAF, IBS_invPrivStd_DAF,
ITU_ch4_DAF, ITU_inv_DAF, ITU_invPrivStd_DAF,
JPT_ch4_DAF, JPT_inv_DAF, JPT_invPrivStd_DAF,
KHV_ch4_DAF, KHV_inv_DAF, KHV_invPrivStd_DAF,
MSL_ch4_DAF, MSL_inv_DAF, MSL_invPrivStd_DAF,
MXL_ch4_DAF, MXL_inv_DAF, MXL_invPrivStd_DAF,
PEL_ch4_DAF, PEL_inv_DAF, PEL_invPrivStd_DAF,
PJL_ch4_DAF, PJL_inv_DAF, PJL_invPrivStd_DAF,
PUR_ch4_DAF, PUR_inv_DAF, PUR_invPrivStd_DAF,
STU_ch4_DAF, STU_inv_DAF, STU_invPrivStd_DAF,
)

##################################################################################
### END OF SFS #### END OF SFS #### END OF SFS #### END OF SFS ### END OF SFS ####
##################################################################################
##################################################################################

