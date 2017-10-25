####################################################
#	Barbara Bitarello
#
#	Signatures of LTBS in human inversions
#	last modified: 14.11.2016
##
###################################################
#preamble
library(dplyr)
library(data.table)
library(parallel)
library(lattice)
library(SOAR)
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/bedtools_inR.R')

load('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/read_scan_data/list.SCAN.3.RData')

source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/SFS_script.r')
#

inv_tab<-setDT(read.csv('SummaryInfo_45Invs_Frequencies_v3.1.csv', sep=",", header=T, stringsAsFactors=F))

#this is necessary for bedtools to work.
with(list.SCAN[[2]], paste0("chr",Chr))-> list.SCAN[[2]]$Chr
with(list.SCAN[[3]], paste0("chr",Chr))-> list.SCAN[[3]]$Chr

with(list.SCAN[[6]], paste0("chr",Chr))-> list.SCAN[[6]]$Chr
with(list.SCAN[[7]], paste0("chr",Chr))-> list.SCAN[[7]]$Chr


#convert to data.table

setDT(list.SCAN[[2]])-> LWK
setDT(list.SCAN[[3]])-> YRI
setDT(list.SCAN[[6]])-> GBR
setDT(list.SCAN[[7]])-> TSI

#

inv_tab %>% mutate(lengthLong=EndBP2-StartBP1, lengthShort=StartBP2-EndBP1) -> inv_tab2

setDT(inv_tab2)



pdf('BP.hist.pdf')

hist(inv_tab2$lengthLong, nclass=30, col=rgb(1,0,0,0.3), border='lightgray', main="Length of predicted inversions", xlab="Length BP", freq=F)
par(new=TRUE)
hist(inv_tab2$lengthShort,  nclass=30, col=rgb(0,0,1,0.3), freq=F, add=T,border='lightgray')
legend("topright", c('Longer', 'Shorter'), col=c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)), pch=1, lty=1)
dev.off()

#run the function

LWK.res<-vector('list',nrow(inv_tab2))
YRI.res<-vector('list',nrow(inv_tab2))
GBR.res<-vector('list',nrow(inv_tab2))
TSI.res<-vector('list',nrow(inv_tab2))

#

for(i in 1: nrow(inv_tab)){
try(bedTools.2in(bed1=select(arrange(LWK,Chr),Chr:End.Win, NCDf3:P.val.NCDf0.3, Dist.NCD.f0.5:Dist.NCD.f0.3, Z.f0.5.P.val:Z.f0.3.P.val),bed2=select(inv_tab2[i,], Chr, StartBP1, EndBP2)))-> temp1
try(bedTools.2in(bed1=select(arrange(LWK,Chr),Chr:End.Win, NCDf3:P.val.NCDf0.3, Dist.NCD.f0.5:Dist.NCD.f0.3, Z.f0.5.P.val:Z.f0.3.P.val),bed2=select(inv_tab2[i,], Chr, EndBP1, StartBP2)))-> temp2
list(BP1=temp1, BP2=temp2)-> LWK.res[[i]]
remove(temp1, temp2)
gc()
try(bedTools.2in(bed1=select(arrange(GBR,Chr),Chr:End.Win, NCDf3:P.val.NCDf0.3, Dist.NCD.f0.5:Dist.NCD.f0.3, Z.f0.5.P.val:Z.f0.3.P.val),bed2=select(inv_tab2[i,], Chr, StartBP1, EndBP2)))-> temp1.GBR
try(bedTools.2in(bed1=select(arrange(GBR,Chr),Chr:End.Win, NCDf3:P.val.NCDf0.3, Dist.NCD.f0.5:Dist.NCD.f0.3, Z.f0.5.P.val:Z.f0.3.P.val),bed2=select(inv_tab2[i,], Chr, EndBP1, StartBP2)))-> temp2.GBR
list(BP1=temp1.GBR, BP2=temp2.GBR)-> GBR.res[[i]]
remove(temp1.GBR, temp2.GBR)
gc()
try(bedTools.2in(bed1=select(arrange(YRI,Chr),Chr:End.Win, NCDf3:P.val.NCDf0.3, Dist.NCD.f0.5:Dist.NCD.f0.3, Z.f0.5.P.val:Z.f0.3.P.val),bed2=select(inv_tab2[i,], Chr, StartBP1, EndBP2)))-> temp1.YRI
try(bedTools.2in(bed1=select(arrange(YRI,Chr),Chr:End.Win, NCDf3:P.val.NCDf0.3, Dist.NCD.f0.5:Dist.NCD.f0.3, Z.f0.5.P.val:Z.f0.3.P.val),bed2=select(inv_tab2[i,], Chr, EndBP1, StartBP2)))-> temp2.YRI
list(BP1=temp1.YRI, BP2=temp2.YRI)-> YRI.res[[i]]
remove(temp1.YRI, temp2.YRI)
gc()
try(bedTools.2in(bed1=select(arrange(TSI,Chr),Chr:End.Win, NCDf3:P.val.NCDf0.3, Dist.NCD.f0.5:Dist.NCD.f0.3, Z.f0.5.P.val:Z.f0.3.P.val),bed2=select(inv_tab2[i,], Chr, StartBP1, EndBP2)))-> temp1.TSI
try(bedTools.2in(bed1=select(arrange(TSI,Chr),Chr:End.Win, NCDf3:P.val.NCDf0.3, Dist.NCD.f0.5:Dist.NCD.f0.3, Z.f0.5.P.val:Z.f0.3.P.val),bed2=select(inv_tab2[i,], Chr, EndBP1, StartBP2)))-> temp2.TSI
list(BP1=temp1.TSI, BP2=temp2.TSI)-> TSI.res[[i]]
remove(temp1.TSI, temp2.TSI)
gc()
cat(i,'done\n')
}

LWK.final<-vector('list', 45)
YRI.final<-vector('list', 45)
GBR.final<-vector('list', 45)
TSI.final<-vector('list', 45)
#
for ( i in 1:45){
if(class(LWK.res[[i]]$BP1)!='try-error'){if(class(LWK.res[[i]]$BP2)!='try-error'){
do.call('rbind', LWK.res[[i]])-> t
setDT(t)
setkey(t, V12)
unique(t)-> LWK.final[[i]]
remove(t)
}
else{
setDT(LWK.res[[i]]$BP1)-> LWK.final[[i]]}}
if(class(LWK.res[[i]]$BP2)!='try-error'){if(class(LWK.res[[i]]$BP1)=='try-error'){
setDT(LWK.res[[i]]$BP2)-> LWK.final[[i]]}}
if(class(LWK.res[[i]]$BP2)=='try-error'){if(class(LWK.res[[i]]$BP1)=='try-error'){
LWK.final[[i]]<-NULL
}
}
#
if(class(TSI.res[[i]]$BP1)!='try-error'){if(class(TSI.res[[i]]$BP2)!='try-error'){
do.call('rbind', TSI.res[[i]])-> t2
setDT(t2)
setkey(t2, V12)
unique(t2)-> TSI.final[[i]]
remove(t2)
}
else{
setDT(TSI.res[[i]]$BP1)-> TSI.final[[i]]}}
if(class(TSI.res[[i]]$BP2)!='try-error'){if(class(TSI.res[[i]]$BP1)=='try-error'){
setDT(TSI.res[[i]]$BP2)-> TSI.final[[i]]}}
if(class(TSI.res[[i]]$BP2)=='try-error'){if(class(TSI.res[[i]]$BP1)=='try-error'){
TSI.final[[i]]<-NULL
}
}
#
if(class(GBR.res[[i]]$BP1)!='try-error'){if(class(GBR.res[[i]]$BP2)!='try-error'){
do.call('rbind', GBR.res[[i]])-> t3
setDT(t3)
setkey(t3, V12)
unique(t3)-> GBR.final[[i]]
remove(t3)
}
else{
setDT(GBR.res[[i]]$BP1)-> GBR.final[[i]]}}
if(class(GBR.res[[i]]$BP2)!='try-error'){if(class(GBR.res[[i]]$BP1)=='try-error'){
setDT(GBR.res[[i]]$BP2)-> GBR.final[[i]]}}
if(class(GBR.res[[i]]$BP2)=='try-error'){if(class(GBR.res[[i]]$BP1)=='try-error'){
GBR.final[[i]]<-NULL
}
}
#
if(class(YRI.res[[i]]$BP1)!='try-error'){if(class(YRI.res[[i]]$BP2)!='try-error'){
do.call('rbind', YRI.res[[i]])-> t4
setDT(t4)
setkey(t4, V12)
unique(t4)-> YRI.final[[i]]
remove(t4)
}
else{
setDT(YRI.res[[i]]$BP1)-> YRI.final[[i]]}}
if(class(YRI.res[[i]]$BP2)!='try-error'){if(class(YRI.res[[i]]$BP1)=='try-error'){
setDT(YRI.res[[i]]$BP2)-> YRI.final[[i]]}}
if(class(YRI.res[[i]]$BP2)=='try-error'){if(class(YRI.res[[i]]$BP1)=='try-error'){
YRI.final[[i]]<-NULL
}
}
cat(i, 'done\n')}


do.call('rbind',LWK.final)-> LWK.final
do.call('rbind',TSI.final)-> TSI.final
do.call('rbind',GBR.final)-> GBR.final
do.call('rbind',YRI.final)-> YRI.final

c(colnames(list.SCAN[[2]])[c(1:3, 8:19, 22:24, 27:29)],'inv.chr', 'inv.beg', 'inv.end', 'overlap')-> cols_names

colnames(LWK.final)<-cols_names
colnames(YRI.final)<-cols_names
colnames(GBR.final)<-cols_names
colnames(TSI.final)<-cols_names

#specific query:
data.frame(Chr='chr4', Beg=40234426, End=40237234)-> query
query.res<-bedTools.2in(bed1=select(arrange(TSI,Chr),Chr:End.Win, NCDf3:P.val.NCDf0.3, Dist.NCD.f0.5:Dist.NCD.f0.3, Z.f0.5.P.val:Z.f0.3.P.val),bed2=query)

colnames(query.res)<-cols_names

#####
#SFS for the inversions

inv_tab2 %>% filter(Chr %in% paste0('chr', seq(1:22))) -> inv_tab3

SFS_inv_allpops<-vector('list', 4)

names(SFS_inv_allpops)<-c('LWK','YRI', 'GBR', 'TSI')

SFS_inv_allpops[['LWK']]<-vector('list',37)
SFS_inv_allpops[['YRI']]<-vector('list',37)
SFS_inv_allpops[['GBR']]<-vector('list',37)
SFS_inv_allpops[['TSI']]<-vector('list',37)

mclapply(1:37, function(x) try(SFS.function(CHR=inv_tab2$Chr[x], BEG=inv_tab2$StartBP1[x], END=inv_tab2$EndBP2[x], POP=2)))-> SFS_inv_allpops[['LWK']]
mclapply(1:37, function(x) try(SFS.function(CHR=inv_tab2$Chr[x], BEG=inv_tab2$StartBP1[x], END=inv_tab2$EndBP2[x], POP=3)))-> SFS_inv_allpops[['YRI']]
mclapply(1:37, function(x) try(SFS.function(CHR=inv_tab2$Chr[x], BEG=inv_tab2$StartBP1[x], END=inv_tab2$EndBP2[x], POP=6)))-> SFS_inv_allpops[['GBR']]
mclapply(1:37, function(x) try(SFS.function(CHR=inv_tab2$Chr[x], BEG=inv_tab2$StartBP1[x], END=inv_tab2$EndBP2[x], POP=7)))-> SFS_inv_allpops[['TSI']]

which(unlist(sapply(SFS_inv_allpops[['LWK']], function(x) length(x)))>1)-> A

which(sapply(SFS_inv_allpops[['LWK']][A], function(x) sum(x)==0))->B

if(length(B)>0){
A[-B]->B1
}
if(length(B)==0){A-> B1}


for ( i in A){
pdf(paste0(inv_tab3$Inversion[i], '.LWK.pdf'))
n<-length(SFS_inv_allpops[['LWK']][[i]])
print(histogram(SFS_inv_allpops[['LWK']][[i]][SFS_inv_allpops[['LWK']][[i]]!=0 & SFS_inv_allpops[['LWK']][[i]]!=100], col='lightgray', main=paste0(inv_tab3$Inversion[[i]]," (", n, " SNPs)"), xlab=list(label='DAC', cex=2), ylab=list(label='Relative Frequency (%)',cex=2),breaks=seq(from=0,to=100,by=2), scales=list(x=list(cex=1.5),y=list(cex=1.5)), auto.key = list(lines=F),ylim=c(0, 100), border=F))
dev.off()
cat(inv_tab3$Inversion[[i]], 'done\n')
}

for ( i in A){
pdf(paste0(inv_tab3$Inversion[i], '.YRI.pdf'))
n<-length(SFS_inv_allpops[['YRI']][[i]])
print(histogram(SFS_inv_allpops[['YRI']][[i]][SFS_inv_allpops[['YRI']][[i]]!=0 & SFS_inv_allpops[['YRI']][[i]]!=100], col='lightgray', main=paste0(inv_tab3$Inversion[[i]]," (", n, " SNPs)"), xlab=list(label='DAC', cex=2), ylab=list(label='Relative Frequency (%)',cex=2),breaks=seq(from=0,to=100,by=2), scales=list(x=list(cex=1.5),y=list(cex=1.5)), auto.key = list(lines=F),ylim=c(0, 100), border=F))
dev.off()
cat(inv_tab3$Inversion[[i]], 'done\n')
}

for ( i in A){
pdf(paste0(inv_tab3$Inversion[i], '.GBR.pdf'))
n<-length(SFS_inv_allpops[['GBR']][[i]])
print(histogram(SFS_inv_allpops[['GBR']][[i]][SFS_inv_allpops[['GBR']][[i]]!=0 & SFS_inv_allpops[['GBR']][[i]]!=100], col='lightgray', main=paste0(inv_tab3$Inversion[[i]]," (", n, " SNPs)"), xlab=list(label='DAC', cex=2), ylab=list(label='Relative Frequency (%)',cex=2),breaks=seq(from=0,to=100,by=2), scales=list(x=list(cex=1.5),y=list(cex=1.5)), auto.key = list(lines=F),ylim=c(0, 100), border=F))
dev.off()
cat(inv_tab3$Inversion[[i]], 'done\n')
}

for ( i in A){
pdf(paste0(inv_tab3$Inversion[i], '.TSI.pdf'))
n<-length(SFS_inv_allpops[['TSI']][[i]])
print(histogram(SFS_inv_allpops[['TSI']][[i]][SFS_inv_allpops[['TSI']][[i]]!=0 & SFS_inv_allpops[['TSI']][[i]]!=100], col='lightgray', main=paste0(inv_tab3$Inversion[[i]]," (", n, " SNPs)"), xlab=list(label='DAC', cex=2), ylab=list(label='Relative Frequency (%)',cex=2),breaks=seq(from=0,to=100,by=2), scales=list(x=list(cex=1.5),y=list(cex=1.5)), auto.key = list(lines=F),ylim=c(0, 100), border=F))
dev.off()
cat(inv_tab3$Inversion[[i]], 'done\n')
}

try(SFS.function(CHR=query$Chr, BEG=query$Beg, END=query$End, POP=2))-> query.LWK
try(SFS.function(CHR=query$Chr, BEG=query$Beg, END=query$End, POP=3))-> query.YRI
try(SFS.function(CHR=query$Chr, BEG=query$Beg, END=query$End, POP=6))-> query.GBR
try(SFS.function(CHR=query$Chr, BEG=query$Beg, END=query$End, POP=7))-> query.TSI

n<-length(query.LWK)
pdf(paste0(paste(query[1,], collapse="|"), ".pdf"))

print(histogram(query.LWK[which(query.LWK!=0 & query.LWK!=100)], col='lightgray', main=paste0(paste(query[1,], collapse="|")," ", "LWK", " ", n, " SNPs"), xlab=list(label='DAC', cex=2), ylab=list(label='Relative Frequency (%)',cex=2),breaks=seq(from=0,to=100,by=2), scales=list(x=list(cex=1.5),y=list(cex=1.5)), auto.key = list(lines=F),ylim=c(0, 100), border=F))

print(histogram(query.YRI[which(query.YRI!=0 & query.YRI!=100)], col='lightgray', main=paste0(paste(query[1,], collapse="|")," ", "YRI", " ", n, " SNPs"), xlab=list(label='DAC', cex=2), ylab=list(label='Relative Frequency (%)',cex=2),breaks=seq(from=0,to=100,by=2), scales=list(x=list(cex=1.5),y=list(cex=1.5)), auto.key = list(lines=F),ylim=c(0, 100), border=F))

print(histogram(query.GBR[which(query.GBR!=0 & query.GBR!=100)], col='lightgray', main=paste0(paste(query[1,], collapse="|")," ", "GBR", " ", n, " SNPs"), xlab=list(label='DAC', cex=2), ylab=list(label='Relative Frequency (%)',cex=2),breaks=seq(from=0,to=100,by=2), scales=list(x=list(cex=1.5),y=list(cex=1.5)), auto.key = list(lines=F),ylim=c(0, 100), border=F))

print(histogram(query.TSI[which(query.TSI!=0 & query.TSI!=100)], col='lightgray', main=paste0(paste(query[1,], collapse="|")," ", "TSI", " ", n, " SNPs"), xlab=list(label='DAC', cex=2), ylab=list(label='Relative Frequency (%)',cex=2),breaks=seq(from=0,to=100,by=2), scales=list(x=list(cex=1.5),y=list(cex=1.5)), auto.key = list(lines=F),ylim=c(0, 100), border=F))

dev.off()

