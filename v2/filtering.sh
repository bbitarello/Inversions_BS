###
#last updated on 10/12/2017
#pantro2
rm hg19.pantro2.bed
touch hg19.pantro2.bed

for k in {1..22} X; do cat /mnt/sequencedb/PopGen/cesare/hg19/bedfiles/tmp/hg19.pantro2_${k}.bed >> /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/v2/hg19.pantro2.bed; done
sed -i 's/\s+/\t/g' hg19.pantro2.bed

sed -i 's/^/chr/'   hg19.pantro2.bed

sed -i '/-1/d' hg19.pantro2.bed

sort -k1,1V -k2,2n -k3,3n hg19.pantro2.bed > hg19.pantro2.sort.bed

#bgzip hg19.pantro2.bed 
bgzip hg19.pantro2.sort.bed

tabix -p vcf hg19.pantro2.sort.bed.gz

#pantrp 4 


cp /mnt/expressions/Mike/scATAC/2_analyses/ben_ableToLift/panTro4/hg19_panTro4/hg19_sort.bed .
mv hg19_sort.bed hg19_pantro4_sort.bed

bgzip hg19_pantro4_sort.bed

tabix -p vcf hg19_pantro4_sort.bed.gz

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
sed 's/\t/  /'  /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/v2/${i}_chrX_AF.frq |sed 's/\t/  /'|sed 's/\t/  /'|sed 's/\t/  /'|sed 's/\t/;/g'|sed 's/  /\t/g' | sed '1d' | sed 's/^/chr/' > /mnt/sequencedb/PopGen/barbara/collaborations/Inversions_BS/v2/${i}_chrX_AF_2.frq
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

for i in {1..22} X; do
awk 'OFS="\t"{print $1,$2,$2,$3,$4,$5}' CHB_chr${i}_AF_2.frq > CHB_chr${i}_AF_3.frq
awk 'OFS="\t"{print $1,$2,$2,$3,$4,$5}' CEU_chr${i}_AF_2.frq > CEU_chr${i}_AF_3.frq
awk 'OFS="\t"{print $1,$2,$2,$3,$4,$5}' LWK_chr${i}_AF_2.frq > LWK_chr${i}_AF_3.frq
awk 'OFS="\t"{print $1,$2,$2,$3,$4,$5}' GIH_chr${i}_AF_2.frq > GIH_chr${i}_AF_3.frq
awk 'OFS="\t"{print $1,$2,$2,$3,$4,$5}' JPT_chr${i}_AF_2.frq > JPT_chr${i}_AF_3.frq
awk 'OFS="\t"{print $1,$2,$2,$3,$4,$5}' TSI_chr${i}_AF_2.frq > TSI_chr${i}_AF_3.frq
awk 'OFS="\t"{print $1,$2,$2,$3,$4,$5}' YRI_chr${i}_AF_2.frq > YRI_chr${i}_AF_3.frq
echo ${i}
done

#for i in CHB CEU LWK GIH JPT TSI YRI; do
#awk 'OFS="\t"{print $1,$2,$2,$3,$4,$5}' ${i}_chrX_AF_2.frq > ${i}_chrX_AF_3.frq
#echo ${i}
#done

# concatenate all chromsomes

rm *_all.frq*

touch CHB_all.frq
touch CEU_all.frq
touch LWK_all.frq
touch GIH_all.frq
touch JPT_all.frq
touch TSI_all.frq
touch YRI_all.frq

for i in {1..22} X;do
cat CHB_chr${i}_AF_3.frq >> CHB_all.frq
cat CEU_chr${i}_AF_3.frq >> CEU_all.frq
cat LWK_chr${i}_AF_3.frq >> LWK_all.frq
cat GIH_chr${i}_AF_3.frq >> GIH_all.frq
cat JPT_chr${i}_AF_3.frq >> JPT_all.frq
cat TSI_chr${i}_AF_3.frq >> TSI_all.frq
cat YRI_chr${i}_AF_3.frq >> YRI_all.frq
echo ${i}
done


for i in CHB CEU LWK GIH JPT TSI YRI;do
sed -i 's/^chr//' ${i}_all.frq
sort -k1,1V -k2,2n -k3,3n ${i}_all.frq |sed 's/^/chr/' > ${i}_all_2.frq
echo ${i}
done

#bgzip and tabix -p bed

for i in CHB CEU LWK GIH JPT TSI YRI; do
bgzip ${i}_all_2.frq
echo $i
done


for i in CHB CEU LWK GIH JPT TSI YRI; do
tabix -p bed ${i}_all_2.frq.gz
echo $i
done

#remove duplicates
rm _AF_3*

#bgzip and tabix -p vcf the non strict mask from 1000G
bgzip 20141020.pilot_mask.whole_genome.bed
tabix -p bed 20141020.pilot_mask.whole_genome.bed.gz 

#bgzip and tabix -p vcf hg19.pantro

#for each chr and pop, intersectbed twice:with mask and with hg19.pantro
for i in CHB CEU LWK GIH JPT TSI YRI; do
intersectBed -a ${i}_all_2.frq.gz -b 20141020.pilot_mask.whole_genome.bed.gz  >  ${i}_pilot_mask.bed;
wait
bgzip ${i}_pilot_mask.bed;
wait
tabix -p bed ${i}_pilot_mask.bed.gz;
echo ${i}

done

#intersect pilot mask and hg19 pantro2
for i in CHB CEU LWK GIH JPT TSI YRI; do
intersectBed -a ${i}_pilot_mask.bed.gz -b hg19.pantro2.sort.bed.gz > ${i}_final.bed
bgzip ${i}_final.bed
wait
tabix -p bed ${i}_final.bed.gz
echo $i
done

#now with pantro4 10/12/2017
for i in CHB CEU LWK GIH JPT TSI YRI; do
intersectBed -a ${i}_pilot_mask.bed.gz -b hg19_pantro4_sort.bed.gz > ${i}_final.bed
bgzip ${i}_final.bed
wait
tabix -p bed ${i}_final.bed.gz
echo $i
done


###CHIMP AND OTHER PRIMATES INFO

for CHROMOSOME in {1..22} X;  
do
cp /mnt/scratch/kay/ape_table/${CHROMOSOME}.tab.gz .
done

#info about these files

# get ref-human, apes and Altai, Den, Vindija
for i in `seq 1 22` ; do /home/pruefer/src/BamSNPTool/BamSNPJoinGen <(gunzip -c /mnt/solexa/Genomes/hg19_evan/$i.fa.gz | ./faToTable.pl | /home/pruefer/src/BamSNPTool/BamSNPAddMaf /mnt/sequencedb/ucsc/goldenPath/hg19/vsPanTro4/axtNet/mafnet/chr$i.maf panTro4 | /home/pruefer/src/BamSNPTool/BamSNPAddMaf /mnt/454/BonoboGenome/chain_net/bonobo/i7_vs_hg19/mafnet/add_dot/chr$i.maf bonobo| /home/pruefer/src/BamSNPTool/BamSNPAddMaf /mnt/sequencedb/ucsc/goldenPath/hg19/vsGorGor3/axtNet/mafnet/chr$i.hg19.gorGor3.net.maf gorgor3 | /home/pruefer/src/BamSNPTool/BamSNPAddMaf /mnt/sequencedb/ucsc/goldenPath/hg19/vsPonAbe2/axtNet/mafnet/chr$i.hg19.ponAbe2.maf ponabe2 | /home/pruefer/src/BamSNPTool/BamSNPAddMaf /mnt/sequencedb/ucsc/goldenPath/hg19/vsRheMac3/axtNet/mafnet/chr$i.maf rheMac3) <(zcat /mnt/454/Vindija/high_cov/genotypes/Vindija33.19/chr${i}_mq25_mapab100.vcf.gz | vcf2table.pl -A) <(zcat /mnt/454/Vindija/high_cov/genotypes/Altai/chr${i}_mq25_mapab100.vcf.gz | vcf2table.pl -A) <(zcat /mnt/454/Vindija/high_cov/genotypes/Denisova/chr${i}_mq25_mapab100.vcf.gz | vcf2table.pl -A) | perl -lane 'print uc($_)' |  gzip > ape_table/$i.tab.gz & done




# columns:
1 chr
2 pos
3 hg19
4 pantro4
5 panpan1.1
6 gorgor3
7 ponabe2
8 rhemac3
9 Vindija33.19 (with ambiguity codes for hets)
10 Altai (with ambiguity codes for hets)
11 Denisova (with ambiguity codes for hets)

#pruefer@bionc05:/mnt/scratch/kay/ape_table$ for i in `seq 1 22` X ; do bedtools intersect -a <(zcat $i.tab.gz | vcf2bed) -b /mnt/454/HCNDCAM/Hengs_Alignability_Filter/hs37m_filt35_99.bed.gz | bed2vcf | gzip > manifesto/$i.tab.gz & done

