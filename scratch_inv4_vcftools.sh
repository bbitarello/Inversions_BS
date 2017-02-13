cd /mnt/sequencedb/1000Genomes/ftp/phase3/20140910

vcftools --gzvcf ALL.chr4.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --freq --chr 4 --positions /mnt/sequencedb/PopGen/barbara/inversions_balsel/Inversion_chr4.txt --out /mnt/sequencedb/PopGen/barbara/inversions_balsel/OUT_inv_chr4.txt


#vcftools --gzvcf ALL.chr4.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --freq --chr 4 --positions /mnt/sequencedb/PopGen/barbara/inversions_balsel/PrivateInv.txt --out /mnt/scratch/barbara/chr4_inv_private_inv.txt

#vcftools --gzvcf ALL.chr4.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --freq --chr 4 --positions /mnt/sequencedb/PopGen/barbara/inversions_balsel/PrivateHet.txt --out /mnt/scratch/barbara/chr4_inv_private_het.txt

#vcftools --gzvcf ALL.chr4.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --freq --chr 4 --positions /mnt/sequencedb/PopGen/barbara/inversions_balsel/PrivateStd.txt --out /mnt/scratch/barbara/chr4_inv_private_std.txt

#vcftools --gzvcf ALL.chr4.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --freq --chr 4 --positions /mnt/sequencedb/PopGen/barbara/inversions_balsel/MultiAllelic.txt --out /mnt/scratch/barbara/chr4_inv_multiallelic.txt

#vcftools --gzvcf ALL.chr4.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --freq --chr 4 --positions /mnt/sequencedb/PopGen/barbara/inversions_balsel/Monomorphic.txt --out /mnt/scratch/barbara/chr4_inv_monomorphic.txt

#vcftools --gzvcf ALL.chr4.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --freq --chr 4 --positions /mnt/sequencedb/PopGen/barbara/inversions_balsel/Shared.txt --out /mnt/scratch/barbara/chr4_inv_shared.txt

#done
