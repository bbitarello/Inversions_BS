**********************************************************************************

Inversions Project (participation in Mario Caceres' Project)

Investigating signatures of balancing selection in human chromosomal inversions

Author: Bárbara Domingues Bitarello

Last modified: 28.03.2017

**********************************************************************************

Use Felix's script to recode vcf and make ALT==DERIVED

```
gunzip -c /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz |awk 'length($4)==1 && length($5)==1{print $0}'| /mnt/scratch/josh/1000Genomes/scripts/aa2ref_j_edit_2.py -cwp >  /mnt/sequencedb/PopGen/barbara/inversions_balsel/chr4_1000gph3_DAF.vcf

cat header.txt chr4_1000gph3_DAF.vcf |bgzip -c > chr4_1000gph3_DAF.vcf.gz

rm chr4_1000gph3_DAF.vcf

tabix -p vcf chr4_1000gph3_DAF.vcf.gz 
```


************
Get samples from each population

```
grep LWK /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/inversions_balsel/input_files/LWK_samples.txt
grep YRI /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/inversions_balsel/input_files/YRI_samples.txt
grep GBR /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/inversions_balsel/input_files/GBR_samples.txt
grep TSI /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/integrated_call_samples_v3.20130502.ALL.panel|awk '{print $1}' > /mnt/sequencedb/PopGen/barbara/inversions_balsel/input_files/TSI_samples.txt
```


***************************************************************************
Subset VCF files for chr4 for each pop separately and calculate frequencies

see phase3_1000g_inversion_chr4.R


**************************************************************************
For each pop, select SNPs contained in the inversion and plot SFS

see phase3_1000g_inversion_chr4.R

****************************************************************************
For each pop, select PrivateSTD SNPs contained in the inversion and plot SFS

********************************************************************************
Use MAF to calculate *NCD1* and *NCD2* for the entire chromosome in windows of 3 kb

