##This needs:

#X-> a data.table with a column called MAF (abolsute frequency). 
#For NCD2: we need a column in X called 'REF' with the reference human (or whatever) allele. And Y, a data.table with CHR, POS, Ref(ingroup), Ref(outgroup). with this we will check the FD status.
#N-> number of chr in the sample 
#CHR-> chr of the data being analysed (it's better to run them separately)
#W-> size of window
#S-> SLIDE szie (e.g. W/2, but could also be 1, to check each SNP...)


NCD1<-function(X,W=100000,S=50000, cores=30){ #: data table, W= window size S: slide window 
b<-seq(X[1,POS],X[nrow(X), POS],S); #beggining positions
e<-b+W; #end positions
CHR<-unique(X$CHR); #make error message if this is longer than one.
id<-paste0(CHR,"|", b, "|", e);
length(b)-1 -> len
tmp<-mclapply2(1:len, 
function(i) dplyr::filter(X, POS>=b[i], POS<e[i]) %>% as.data.table, mc.cores=cores);
print(paste0("Finished setting up SNP for chr", CHR));
ncd1_a<-mclapply2(1:len, function(i)
tmp[[i]] %>% dplyr::summarise(Beg.Win=b[i], End.Win=e[i],N_SNPs=n(), Win.ID=id[i]) %>% mutate(Chr=CHR) %>% as.data.table, mc.cores=cores);
ncd1_b<- mclapply2(1:len, function(j)
tmp[[j]] %>% dplyr::filter(!(MAF==0|MAF==1)) %>% dplyr::summarise(N_SNPs_cor= n(),NCD1_tf0.5=sqrt(sum((MAF-0.5)^2)/N_SNPs_cor),NCD1_tf0.4=sqrt(sum((MAF-0.4)^2)/N_SNPs_cor), NCD1_tf0.3=sqrt(sum((MAF-0.3)^2)/N_SNPs_cor)) %>% mutate(Win.ID=id[j]) %>% as.data.table, mc.cores=cores);
print(paste0("Finished filtering segregating SNPs for chr ", CHR));
ncd1_c<- do.call(rbind, mclapply2(1:len,function(i) left_join(ncd1_a[[i]], ncd1_b[[i]], by="Win.ID"),mc.cores=cores)) %>% as.data.table;
print(paste0("Finished NCD1 calculations for chr", CHR));
return(ncd1_c); }

############################################################################################
####### NCD2 ######### NCD2 ########### NCD2 ########## NCD2 ############ NCD2 #############
############################################################################################
#mclapply2

NCD2<-function(X,Y,W=3000,S=1500, cores=30){
b<-seq(X[1,POS],X[nrow(X), POS],S); #beggining positions
e<-b+W; #end positions
CHR<-unique(X$CHR); #make error message if this is longer than one.
id<-paste0(CHR,"|",b, "|", e);
length(b)-1 -> len
tmp_snp<-mclapply2(1:len,function(i) dplyr::filter(X, POS>=b[i], POS<e[i]) %>% as.data.table, mc.cores=cores); #SNP
print(paste0("Finished setting up SNPs for chr",CHR));
tmp_fd<-mclapply2(1:len,function(i) dplyr::filter(Y, POS>=b[i], POS<e[i]) %>% as.data.table, mc.cores=cores); #FD
#filter MAF==0 and MAF==1
print(paste0("Finished setting up FDs for chr", CHR));
tmp2_snp<-mclapply2(1:len, function(x) tmp_snp[[x]] %>% dplyr::filter(!(MAF==0|MAF==1)) %>% as.data.table, mc.cores=cores);
print(paste0("Finished filtering segregating SNPs for chr", CHR));
tmp_fd_SNP<- mclapply2(1:len, function(x) as.vector(inner_join(tmp_fd[[x]], tmp2_snp[[x]], by="ID") %>% filter(ALT==Chimp_REF) %>% select(ID) %>% as.data.table)[,ID], mc.cores=cores); #Fd positions that are also SNP positions.Save positions to remove from FD

tmp3_snp<-mclapply2(1:len, function(j) tmp_snp[[j]]  %>% dplyr::summarise(Beg.Win=b[j], End.Win=e[j],N_SNPs=n(), Win.ID=id[j]) %>% mutate(Chr=CHR) %>% as.data.table,mc.cores=cores); #initital # SNPs

tmp4_snp<-mclapply2(1:len, function(j) tmp2_snp[[j]]  %>% dplyr::summarise(N_SNPs_cor=n(), Win.ID=id[j]) %>% as.data.table,mc.cores=cores); #actual number of segregating SNPs

tmp2_fd <-mclapply2(1:len, function(j) tmp_fd[[j]]  %>% dplyr::summarise(N_FDs=n(),Win.ID=id[j]) %>%  as.data.table,mc.cores=cores); #initial # FDs

tmp3_fd <- mclapply2(1:len, function(j) filter(tmp_fd[[j]], !(ID %in% tmp_fd_SNP[[j]])) %>% dplyr::summarise(N_FDs_cor=n(), Win.ID=id[j]) %>% as.data.table, mc.cores=cores); #corrected number of FDs
print(paste0("Finished correcting FD number for chr", CHR));

tmp4<-mclapply2(1:len, function(x) data.frame(Win.ID=id[x], NCD2_tf0.5=sqrt(sum((c(rep(0,tmp3_fd[[x]]$N_FDs_cor),tmp4_snp[[x]]$MAF)-0.5)^2)/(tmp4_snp[[x]]$N_SNPs_cor+tmp3_fd[[x]]$N_FDs_cor)), NCD2_tf0.4=sqrt(sum((c(rep(0,tmp3_fd[[x]]$N_FDs_cor),tmp4_snp[[x]]$MAF)-0.4)^2)/(tmp4_snp[[x]]$N_SNPs_cor+tmp3_fd[[x]]$N_FDs_cor)), NCD2_tf0.3=sqrt(sum((c(rep(0,tmp3_fd[[x]]$N_FDs_cor),tmp4_snp[[x]]$MAF)-0.3)^2)/(tmp4_snp[[x]]$N_SNPs_cor+tmp3_fd[[x]]$N_FDs_cor))) %>% as.data.table, mc.cores=cores);
print(paste0("Finished NCD2 calculations for chr", CHR));

tmp5<-mclapply2(1: len, function(x) join_all(list(tmp3_snp[[x]],tmp4_snp[[x]],tmp2_fd[[x]], tmp3_fd[[x]], tmp4[[x]]), by='Win.ID', type='left') %>% as.data.table %>% mutate(PtoD=N_SNPs_cor/(N_FDs_cor+1)) %>% as.data.table,mc.cores=cores);
print(paste0("Combining all info into data table"));
ncd2<- do.call(rbind, tmp5) %>% as.data.table;
return(ncd2);
print("Done!");
}
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

#mclapply wrapper that prints stages of parallel running
#copied from here: http://stackoverflow.com/questions/10984556/is-there-way-to-track-progress-on-a-mclapply
#last mdofied: 14.12.2016
#

mclapply2 <- function(X, FUN, ..., 
    mc.preschedule = FALSE, mc.set.seed = TRUE,
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

