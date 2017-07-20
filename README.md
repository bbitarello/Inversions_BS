**********************************************************************************

Inversions Project (participation in Mario Caceres' Project)

Investigating signatures of balancing selection in human chromosomal inversions

Author: Bárbara Domingues Bitarello

Last modified: 19.07.2017

**********************************************************************************

Use Felix's script to recode vcf and make ALT==DERIVED


```
gunzip -c /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz |awk 'length($4)==1 && length($5)==1{print $0}'| /mnt/scratch/josh/1000Genomes/scripts/aa2ref_j_edit_2.py -cwp >  /mnt/sequencedb/PopGen/barbara/collaborations/inversions_BS/chr4_1000gph3_DAF.vcf

gunzip -c /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz|head -250 > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/header.txt

cat input_files/header.txt chr4_1000gph3_DAF.vcf |bgzip -c > chr4_1000gph3_DAF.vcf.gz

rm chr4_1000gph3_DAF.vcf

tabix -p vcf chr4_1000gph3_DAF.vcf.gz 
```


************
Get samples from each population

```

grep ACB /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/ACB_samples.txt
grep ASW /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/ASW_samples.txt
grep BEB /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/BEB_samples.txt
grep CDX /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/CDX_samples.txt
grep CEU /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/CEU_samples.txt
grep CHB /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/CHB_samples.txt
grep CHS /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/CHS_samples.txt
grep CLM /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/CLM_samples.txt
grep ESN /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/ESN_samples.txt
grep FIN /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/FIN_samples.txt
grep GBR /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/GBR_samples.txt
grep GIH /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/GIH_samples.txt
grep GWD /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/GWD_samples.txt
grep IBS /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/IBS_samples.txt
grep ITU /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/ITU_samples.txt
grep JPT /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/JPT_samples.txt
grep KHV /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/KHV_samples.txt
grep LWK /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/LWK_samples.txt
grep MSL /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/MSL_samples.txt
grep MXL /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/MXL_samples.txt
grep PEL /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/PEL_samples.txt
grep PJL /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/PJL_samples.txt
grep PUR /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/PUR_samples.txt
grep STU /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/STU_samples.txt
grep TSI /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/TSI_samples.txt
grep YRI /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/YRI_samples.txt

```


***************************************************************************
Subset VCF files for chr4 for each pop separately and calculate frequencies

see phase3_1000g_inversion_chr4.R


**************************************************************************
For each pop, select SNPs contained in the inversion and plot SFS

see phase3_1000g_inversion_chr4.R



****************************************************************************
For each pop, select PrivateSTD SNPs contained in the inversion and plot SFS

see phase3_1000g_inversions_chr4.R


********************************************************************************************************************
Calculate PtoD for all SNPs and only for PrivateStd SNPs. Compare to distributio for similar-sized windows for chr4.
********************************************************************************************************************

For this one we do not need DAF, so we need a new VCF. Use the sample list and vcftools to calculate frequencies per pop.



********************************************************************************
Use MAF to calculate *NCD1* and *NCD2* for the entire chromosome in windows of 3 kb


Need to update NCD scripts. USe data.table.
