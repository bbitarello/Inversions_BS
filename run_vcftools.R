#last modified: 20.0.2017
library(doMC)
registerDoMC(11)
library(parallel)

system.time(mclapply2(1:22, function(i) system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/ACB_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/ACB_chr', i, '_AF'))))

system.time(mclapply2(1:22, function(i) system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/ASW_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/ASW_chr', i, '_AF'))))

system.time(mclapply2(1:22, function(i) system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/BEB_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/BEB_chr', i, '_AF'))))

system.time(mclapply2(1:22, function(i) system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/CDX_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/CDX_chr', i, '_AF'))))

system.time(mclapply2(1:22, function(i) system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/CEU_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/CEU_chr', i, '_AF'))))

system.time(mclapply2(1:22, function(i) system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/CHB_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/CHB_chr', i, '_AF'))))

system.time(mclapply2(1:22, function(i) system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/CHS_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/CHS_chr', i, '_AF'))))

system.time(mclapply2(1:22, function(i) system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/CLM_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/CLM_chr', i, '_AF'))))

system.time(mclapply2(1:22, function(i) system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/ESN_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/ESN_chr', i, '_AF'))))

system.time(mclapply2(1:22, function(i) system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/FIN_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/FIN_chr', i, '_AF'))))


system.time(foreach(i=1:22) %dopar%
system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/GBR_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/GBR_chr', i, '_AF')))


system.time(foreach(i=1:22) %dopar%
system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/GIH_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/GIH_chr', i, '_AF')))

system.time(foreach(i=1:22) %dopar%
system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/GWD_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/GWD_chr', i, '_AF')))


system.time(foreach(i=1:22) %dopar%
system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/IBS_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/IBS_chr', i, '_AF')))

system.time(foreach(i=1:22) %dopar%
system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/ITU_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/ITU_chr', i, '_AF')))


system.time(foreach(i=1:22) %dopar%
system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/JPT_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/JPT_chr', i, '_AF')))


system.time(foreach(i=1:22) %dopar%
system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/KHV_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/KHV_chr', i, '_AF')))

system.time(foreach(i=1:22) %dopar%
system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/LWK_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/LWK_chr', i, '_AF')))

system.time(foreach(i=1:22) %dopar%
system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/MSL_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/MSL_chr', i, '_AF')))


system.time(foreach(i=1:22) %dopar%
system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/MXL_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/MXL_chr', i, '_AF')))

system.time(foreach(i=1:22) %dopar%
system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/PEL_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/PEL_chr', i, '_AF')))


system.time(foreach(i=1:22) %dopar%
system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/PJL_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/PJL_chr', i, '_AF')))

system.time(foreach(i=1:22) %dopar%
system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/PUR_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/PUR_chr', i, '_AF')))

system.time(foreach(i=1:22) %dopar%
system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/STU_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/STU_chr', i, '_AF')))

system.time(foreach(i=1:22) %dopar%
system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/TSI_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/TSI_chr', i, '_AF')))

system.time(foreach(i=1:22) %dopar%
system(paste0('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr',i,'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/YRI_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/YRI_chr', i, '_AF')))

#clean


mclapply2(1:22, function(i)
system(paste0('sed \'s/\t/  /\'  input_files/LWK_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/LWK_chr', i, '_AF_2.frq')))
mclapply2(1:22, function(i)
system(paste0('sed \'s/\t/  /\'  input_files/YRI_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/YRI_chr', i, '_AF_2.frq')))

mclapply2(1:22, function(i)
system(paste0('sed \'s/\t/  /\'  input_files/TSI_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/TSI_chr', i, '_AF_2.frq')))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/ACB_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/ACB_chr', i, '_AF_2.frq'))


foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/ASW_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/ASW_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/BEB_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/BEB_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/CDX_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/CDX_chr', i, '_AF_2.frq'))


foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/CEU_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/CEU_chr', i, '_AF_2.frq'))

ch(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/CHB_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/CHB_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/CHS_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/CHS_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/CLM_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/CLM_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/ESN_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/ESN_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/FIN_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/FIN_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/GBR_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/GBR_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/GIH_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/GIH_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/GWD_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/GWD_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/IBS_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/IBS_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/ITU_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/ITU_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/JPT_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/JPT_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/KHV_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/KHV_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/LWK_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/LWK_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/MSL_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/MSL_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/MXL_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/MXL_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/PEL_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/PEL_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/PJL_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/PJL_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/PUR_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/PUR_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/STU_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/STU_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar% 
system(paste0('sed \'s/\t/  /\'  input_files/TSI_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/TSI_chr', i, '_AF_2.frq'))

foreach(i=1:22) %dopar%
system(paste0('sed \'s/\t/  /\'  input_files/YRI_chr', i, '_AF.frq  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/g\'|sed \'s/  /\t/g\' > input_files/YRI_chr', i, '_AF_2.frq'))
#remove duplicated files
system('rm input_files/*AF.frq')
system('rm input_files/*AF.log')


