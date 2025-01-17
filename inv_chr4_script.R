###################################################################################
##        Barbara Bitarello
##
#
#       Signatures of LTBS in a particular region in chr4
#       last modified: 10.02.2017
#       Obs: I am using Phase 3 because the SNs they sent me are from Phase 3
###################################################################################
#preamble
library(dplyr)
library(data.table)
library(parallel)
library(lattice)
library(SOAR)
library(ggplot2)
library(plyr)
library(splitstackshape)
library(doMC)

doMC::registerDoMC(cores=22)


mclapply2 <- function(X, FUN, ..., 
    mc.preschedule = TRUE, mc.set.seed = TRUE,
    mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
    mc.cleanup = TRUE, mc.allow.recursive = TRUE,
    mc.progress=TRUE, mc.style=3) 
{
    if (!is.vector(X) || is.object(X)) X <- as.list(X)

    if (mc.progress) {
        f <- fifo(tempfile(), open="w+b", blocking=T)
        p <- parallel:::mcfork()
        pb <- txtProgressBar(0, length(X), style=mc.style)
        setTxtProgressBar(pb, 0) 
        progress <- 0
        if (inherits(p, "masterProcess")) {
            while (progress < length(X)) {
                readBin(f, "double")
                progress <- progress + 1
                setTxtProgressBar(pb, progress) 
            }
            cat("\n")
            parallel:::mcexit()
        }
    }
    tryCatch({
        result <- mclapply(X, ..., function(...) {
                res <- FUN(...)
                if (mc.progress) writeBin(1, f)
                res
            }, 
            mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed,
            mc.silent = mc.silent, mc.cores = mc.cores,
            mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive
        )

    }, finally = {
        if (mc.progress) close(f)
    })
    result
}

#1.first, read in input data
#Inversions file

fread('VariantsClassified_HsInv0102_using434Ph3Samples.txt')-> var_dt


#human-chimp divergence file

'zcat < /raid/genevol/users/barbara/NCV_dir_package/input_data/outgroup_files/fds.chr4.hg19_pantro2.Map50_100.TRF.SDs.bed.gz' %>%
data.table::fread(sep='\t') -> HC_div

colnames(HC_div)<-c('chr', 'pos', 'human', 'chimp')

#use uppercase

HC_div<-HC_div %>% dplyr::mutate(human=toupper(human), chimp=toupper(chimp)) %>% as.data.table

##########################################################################
########################1000G data Phase 3 ###############################
##########################################################################

data <- 
'bzcat < /raid/genevol/1kg/phase3/annotations/all.chr.cadd.annotation4.tsv.bz2' %>%
data.table::fread(sep='\t') %>%
#dplyr::select(-matches("(all)|(noadm)|(MXL)|(PUR)|(CLM)|(PEL)|(ASW)|(ACB)|(BEB)|(CDX)|(STU)|(CHS)|(PJL)|(KHV)|(JPT)|(ITU)|(CHB)|(GIH)")) %>% #remove admix and asian pops
dplyr::filter(CHR==4) %>%
tidyr::gather(POP, AF, starts_with("AF_ALT")) %>%
#dplyr::filter(!is.na(AF), AF != 0, AF != 1) %>% if I filter in thie way i get 562 snps within the inversion. just checking if this is the reason.
dplyr::mutate(POP = gsub("^AF_ALT_(\\w+)$", "\\1", POP)) %>% as.data.table

##  test whether this is better.
system.time(
data2 <- 
    data.table::fread(
        'bzcat < /raid/genevol/1kg/phase3/annotations/all.chr.cadd.annotation4.tsv.bz2',
         sep='\t') %>%
	dplyr::filter(CHR==4) %>% as.data.table
)


#the first one works better. Now, add a MAF column:

data %>% dplyr::group_by(POP, POS) %>% dplyr::mutate(MAF=min(AF)) %>% as.data.table -> data.maf
# the above thoes not work when there is only two alleles and f>0.5
test<-
bind_rows(
data %>% 
dplyr::group_by(POP, CHR, POS) %>%
dplyr::filter(nAL == 2,
data %>%  
dplyr::group_by(POP, CHR, POS) %>%
dplyr::filter(nAL > 2, AF == min(AF))
 ) %>%
 dplyr::arrange(POP, CHR, POS) %>%
 {.$AF[.$AF > 0.5] <- 1 - .$AF[.$AF > 0.5];.} %>%
 dplyr::rename(MAF = AF) %>% 
 as.data.table
 ) #974.843


#further subset to inlcude only range of snps in inversion

#var_dt %>% summarise(min=min(POS), max=max(POS))
#       min      max
#1 40225085 40246984


data.maf %>% dplyr::group_by(POP) %>% filter(POS>=40225085 & POS<=40246984) %>% as.data.table -> inv_chr4_range  #funny thing, it would appear this has 707 SNPs, whereas var_dt has only 607, which makes no sense because they both come from the 1000G phase 3.
setkey(inv_chr4_range, POS)
unique(inv_chr4_range) -> unique_inv_chr4_range

#But actually this fulldata has repeated lined per SNP, because of conflicting annotation.


#select SNPs per type.


dlply(var_dt, "Type")-> snp_type
names(snp_type)
#one option: write files with these snps and use that to filter Phase 3 vcf and calculate frequencies using vcftools. This is to be done at the MPi to see if results converge.


write.table(select(snp_type[[1]],CHROM, POS), file="Monomorphic.txt", row.names=F, col.names=F, quote=F, sep="\t")
write.table(select(snp_type[[6]],CHROM, POS), file="Shared.txt", row.names=F, col.names=F, quote=F, sep="\t")
write.table(select(snp_type[[4]],CHROM, POS), file="PrivateInv.txt", row.names=F, col.names=F, quote=F, sep="\t")
write.table(select(snp_type[[5]],CHROM, POS), file="PrivateStd.txt", row.names=F, col.names=F, quote=F, sep="\t")
write.table(select(snp_type[[3]],CHROM, POS), file="PrivateHet.txt", row.names=F, col.names=F, quote=F, sep="\t")
write.table(select(snp_type[[2]],CHROM, POS), file="MultiAllelic.txt", row.names=F, col.names=F, quote=F, sep="\t")



#3. do the normal checks: SNP needs to be segregating in the pop and should NOT be a FD, and if it is, it should be excluded as FD, UNLESS Alt is different from Chimp reference, in which case it might be both. use old scripts for this



#4. Take FD info for this genomic interval (from first to last SNP? or absed on the coordiantes they sent before?)


#5. for each category of SNP, calculate PtoD
                                                                        
