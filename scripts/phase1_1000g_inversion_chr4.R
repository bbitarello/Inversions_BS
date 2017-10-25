##########################################################################################
#	Use Phase 1 1000G data to verify allele frequencies within inversion in chr4
#	Bárbara Bitarello
#	Created: 14.02.2017
#	Last modified: 15.02.2017
##########################################################################################

#preamble
library(pegas)
library(dplyr)
library(plyr)
library(data.table)
library(parallel)
library(lattice)
library(SOAR)
library(ggplot2)
library(plyr)
library(splitstackshape)
biocLite(VariantAnnotation)
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/bedtools_inR.R')
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/SFS_script.r')
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/mclapply2.R')
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/SFS_script.r')

#load('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/read_scan_data/list.SCAN.3.RData')

####
#read in data 

#Inversions file
fread('VariantsClassified_HsInv0102_using434Ph3Samples.txt')-> var_dt

#human-chimp divergence file

'zcat < /mnt/sequencedb/PopGen/barbara/NCV_dir_package/input_data/outgroup_files/fds.chr4.hg19_pantro2.Map50_100.TRF.SDs.bed.gz' %>%
data.table::fread(sep='\t') -> HC_div

colnames(HC_div)<-c('CHROM', 'POS', 'REF', 'Chimp_REF')

#use uppercase

HC_div<-HC_div %>% dplyr::mutate(REF=toupper(REF), Chimp_REF=toupper(Chimp_REF)) %>% dplyr::filter(POS>=40224426 & POS<= 40247234) %>% as.data.table  #select range of inversions +- 10 kb

HC_div2<-HC_div %>% dplyr::mutate(REF=toupper(REF), Chimp_REF=toupper(Chimp_REF)) %>% as.data.table  #select range of inversions +- 10 kb
##################   SNP DATA ###################
#PHASE 3 1000G
###vcf file chr4

#system('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --chr 4 --from-bp 40224426 --to-bp 40247234 --remove-indels --recode --out  /mnt/sequencedb/PopGen/barbara/inversions_balsel/chr4_inversion_vcf.vcf')

#save header
#system('head -n 250 chr4_inversion_vcf.vcf.recode.vcf > header.txt')
#system('tail -n 599 chr4_inversion_vcf.vcf.recode.vcf |sed \'s/#//\' > chr4_inversion_no_header.vcf') #this has 599 lines while their input data has 607. Need to figure out why. THis might be because I am filtering out structural variants --remove-indels, but when I don't, we get 631 lines instead of 599...

#trying an alternative approach

system('awk \'OFS=\"\t\"{print $1, $2}\' VariantsClassified_HsInv0102_using434Ph3Samples.txt |sed \'1d\' > filter_pos.txt') #positions sent by MC


system('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --chr 4 --positions filter_pos.txt --recode --recode-INFO-all --out /mnt/sequencedb/PopGen/barbara/inversions_balsel/test.vcf')#thus we recapitulate the snps they sent us.

#this here is so that we can get ancestral states (obsolete)

#system('vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/functional_annotation/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5_func_anno.20130502.sites.vcf.gz --chr 4 --positions filter_pos.txt --recode --recode-INFO-all --out /mnt/sequencedb/PopGen/barbara/inversions_balsel/func_annot_chr4_phase3_inversions.vcf')

system('head -n 250 test.vcf.recode.vcf > header.txt')
system('tail -n 607 test.vcf.recode.vcf |sed \'s/#//\' > chr4_inversion_no_header.vcf')

#keep only those that have derived allele info

system('grep \\;AA test.vcf.recode.vcf |grep -v insertion|grep -v deletion|sed \'1d\' > test2.vcf.recode.vcf') #583 SNPs

system('cat header.txt test2.vcf.recode.vcf |bgzip -c > chr4_inversion_vcf.gz') #compress vcf
system('tabix -fp vcf chr4_inversion_vcf.gz') #index vcf
#calculate frequencies using vctools

system('vcftools --gzvcf chr4_inversion_vcf.gz --counts --out /mnt/sequencedb/PopGen/barbara/inversions_balsel/freq_chr4_inv.txt')

system('sed \'s/\t/  /\'  freq_chr4_inv.txt.frq.count  |sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/  /\'|sed \'s/\t/;/\'|sed \'s/\t/;/\'|sed \'s/  /\t/g\' > freq_inv.txt')

fread('freq_inv.txt')-> freq_inv #global frequencies of the SNPs in the inversion. If I want pop specific it will require some further filtering.
colnames(freq_inv)[5]<- gsub(":","_",gsub("\\{|\\}","", colnames(freq_inv)[5])) #fix col name
#index vcf

system('cat header2.txt test2.vcf.recode.vcf > real_vcf')

fread('real_vcf', header=T)-> chr4_inv_vcf #read in vcf

#adding a header to this data.table
#colnames(chr4_inv_vcf)<-unlist(strsplit(gsub("#","",system('tail -n 1 header.txt', intern=T)), "\t"))

chr4_simp<- select(chr4_inv_vcf, CHROM:FORMAT)


#join var_dt to frequency info

left_join(freq_inv, var_dt) %>% as.data.table -> comp_dt



#now use chimp info:

nrow(HC_div) #209 FD in the inversion +- 10 kb
nrow(comp_dt) # 606 SNPs

comp_dt %>% dplyr::group_by(Type) %>% dplyr::summarise(N=n(), PtoD=N/nrow(HC_div))  #PtoD without correcting FDs


(na.omit(left_join(comp_dt, HC_div) %>% as.data.table) %>% dplyr::filter(Chimp_REF==ALT))$POS-> temp #FDs that are not FDs.9 FDs actually have the chimp state segreagating as alternate allele

filter(HC_div, !(POS %in% temp)) #180 FDs



comp_dt %>% dplyr::group_by(Type) %>% dplyr::summarise(N=n(), D=209, Dcor=180, PtoD=N/D, PtoD_cor=N/Dcor, PtoD_all=nrow(var_dt)/209, PtoD_all_cor=nrow(var_dt)/180)  #PtoD without correcting FDs


#phase 1 (felix dataa)

#first, use vcf tools to filter the region

#then read in that vcf file


#something new


fread('calculated_fra.40225087-40246984.ALL.txt', sep="\t")-> ens_inv
colnames(ens_inv)[1]<-'CHROM'

#add ancestral info




##


apply(select(chr4_inv_vcf, INFO), 2, function(x) strsplit(x, ";"))[[1]]-> test
ANCESTRAL<-unlist(mclapply(test, function(x) gsub("\\|","",gsub("AA=","",x[[11]]))))

chr4_inv_vcf %>% mutate(ANC=toupper(ANCESTRAL)) %>% as.data.table -> yay
#i love dplyr
ens_inv %>% dplyr::filter(POS %in% chr4_inv_vcf$POS) %>% mutate(ANC=ANCESTRAL) %>% as.data.table -> yay2
#aso far the 57 cases where ALT is NOT derived, we need to redefine frequencies

yay2 %>% dplyr::mutate(RefAnc=ifelse(REF==ANC, 'TRUE', 'FALSE')) %>% as.data.table %>% filter(RefAnc=='FALSE') %>% as.data.table-> yay3





