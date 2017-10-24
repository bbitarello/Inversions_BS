###

#now do the same for chrX

#calculate allele frequencies
#replace space for tab
#add chr
#remove header
#intersect concatenate with the other chromosomes
#run the other stuff again.

for i in CHB CEU LWK GIH JPT TSI YRI; do
vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1a.20130502.genotypes.vcf.gz --keep /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/input_files/${i}_samples.txt --freq --out /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/v2/${i}_chrX_AF
echo ${i}
done

for i in CHB CEU LWK GIH JPT TSI YRI; do
sed 's/\t/  /'  /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/v2/${i}_chrX_AF.frq |sed 's/\t/  /'|sed 's/\t/  /'|sed 's/\t/  /'|sed 's/\t/;/g'|sed 's/  /\t/g' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/v2/${i}_chrX_AF.frq
echo ${i}
done
rm *_chrX_AF.frq

#copy files from input_files that are names <POP>_chr$i_AF_2.frq. Pops are: YRI, LWK, CHB, JPT, CEU, TSI and GIH


rm *frq

for i in {1..22}; do
cp ../input_files/YRI_chr${i}_AF_2.frq .
cp ../input_files/LWK_chr${i}_AF_2.frq .
cp ../input_files/CHB_chr${i}_AF_2.frq .
cp ../input_files/JPT_chr${i}_AF_2.frq .
cp ../input_files/CEU_chr${i}_AF_2.frq .
cp ../input_files/TSI_chr${i}_AF_2.frq .
cp ../input_files/GIH_chr${i}_AF_2.frq .
done

#remove header

for i in {1..22}; do
sed -i '1d' CHB_chr${i}_AF_2.frq
sed -i '1d' CEU_chr${i}_AF_2.frq 
sed -i '1d' LWK_chr${i}_AF_2.frq
sed -i '1d' GIH_chr${i}_AF_2.frq 
sed -i '1d' JPT_chr${i}_AF_2.frq
sed -i '1d' TSI_chr${i}_AF_2.frq
sed -i '1d' YRI_chr${i}_AF_2.frq
echo ${i}
done

#replace spaces by tab
for i in {1..22}; do
sed -i 's/\s+/\t/g' CHB_chr${i}_AF_2.frq
sed -i 's/\s+/\t/g' CEU_chr${i}_AF_2.frq 
sed -i 's/\s+/\t/g' LWK_chr${i}_AF_2.frq 
sed -i 's/\s+/\t/g' GIH_chr${i}_AF_2.frq 
sed -i 's/\s+/\t/g' JPT_chr${i}_AF_2.frq
sed -i 's/\s+/\t/g' TSI_chr${i}_AF_2.frq 
sed -i 's/\s+/\t/g' YRI_chr${i}_AF_2.frq 
echo ${i}
done

#add 'chr' to first column

for i in {1..22}; do
sed -i 's/^/chr/' CHB_chr${i}_AF_2.frq
sed -i 's/^/chr/' CEU_chr${i}_AF_2.frq
sed -i 's/^/chr/' LWK_chr${i}_AF_2.frq 
sed -i 's/^/chr/' GIH_chr${i}_AF_2.frq 
sed -i 's/^/chr/' JPT_chr${i}_AF_2.frq
sed -i 's/^/chr/' TSI_chr${i}_AF_2.frq 
sed -i 's/^/chr/' YRI_chr${i}_AF_2.frq 
echo ${i}
done

#make it look like a bed file

for i in {1..22}; do
awk 'OFS="\t"{print $1,$2,$2,$3,$4,$5}' CHB_chr${i}_AF_2.frq > CHB_chr${i}_AF_3.frq
awk 'OFS="\t"{print $1,$2,$2,$3,$4,$5}' CEU_chr${i}_AF_2.frq > CEU_chr${i}_AF_3.frq
awk 'OFS="\t"{print $1,$2,$2,$3,$4,$5}' LWK_chr${i}_AF_2.frq > LWK_chr${i}_AF_3.frq
awk 'OFS="\t"{print $1,$2,$2,$3,$4,$5}' GIH_chr${i}_AF_2.frq > GIH_chr${i}_AF_3.frq
awk 'OFS="\t"{print $1,$2,$2,$3,$4,$5}' JPT_chr${i}_AF_2.frq > JPT_chr${i}_AF_3.frq
awk 'OFS="\t"{print $1,$2,$2,$3,$4,$5}' TSI_chr${i}_AF_2.frq > TSI_chr${i}_AF_3.frq
awk 'OFS="\t"{print $1,$2,$2,$3,$4,$5}' YRI_chr${i}_AF_2.frq > YRI_chr${i}_AF_3.frq
done

# concatenate all chromsomes

touch CHB_all.frq
touch CEU_all.frq
touch LWK_all.frq
touch GIH_all.frq
touch JPT_all.frq
touch TSI_all.frq
touch YRI_all.frq

for i in {1..22};do
cat CHB_chr${i}_AF_3.frq >> CHB_all.frq
cat CEU_chr${i}_AF_3.frq >> CEU_all.frq
cat LWK_chr${i}_AF_3.frq >> LWK_all.frq
cat GIH_chr${i}_AF_3.frq >> GIH_all.frq
cat JPT_chr${i}_AF_3.frq >> JPT_all.frq
cat TSI_chr${i}_AF_3.frq >> TSI_all.frq
cat YRI_chr${i}_AF_3.frq >> YRI_all.frq
echo ${i}
done
#bgzip and tabix -p bed

for i in CHB CEU LWK GIH JPT TSI YRI; do
bgzip ${i}_all.frq
echo $i
done


for i in CHB CEU LWK GIH JPT TSI YRI; do
tabix -p bed ${i}_all.frq.gz
echo $i
done

#remove duplicates
rm _AF_3*

#bgzip and tabix -p vcf the non strict mask from 1000G
tabix -p bed 20141020.pilot_mask.whole_genome.bed.gz 

#bgzip and tabix -p vcf hg19.pantro

#for each chr and pop, intersectbed twice:with mask and with hg19.pantro
for i in CHB CEU LWK GIH JPT TSI YRI; do
intersectBed -a ${i}_all.frq.gz -b 20141020.pilot_mask.whole_genome.bed.gz  >  ${i}_pilot_mask.bed;
wait
bgzip ${i}_pilot_mask.bed;
wait
tabix -p bed ${i}_pilot_mask.bed.gz;
echo ${i}

done


#add 'chr' to hg19 pantro
 
sed -i 's/^/chr/' hg19.pantro2.bed

bgzip hg19.pantro2.bed 

tabix -p bed hg19.pantro2.bed.gz


#intersect pilot mask and hg19 pantro2
for i in CHB CEU LWK GIH JPT TSI YRI; do
intersectBed -a ${i}_pilot_mask.bed.gz -b hg19.pantro2.bed.gz > ${i}_final.bed
echo $i
done

