##This needs:

#X-> a data.table with a column called MAF (abolsute frequency). 
#For NCD2: we need a column in X called 'REF' with the reference human (or whatever) allele. And Y, a data.table with CHR, POS, Ref(ingroup), Ref(outgroup). with this we will check the FD status.
#N-> number of chr in the sample 
#CHR-> chr of the data being analysed (it's better to run them separately)
#W-> size of window
#S-> SLIDE szie (e.g. W/2, but could also be 1, to check each SNP...)


NCD1<-function(X,N=198,W=10000,S=5000){
b<-seq(X[1,POS],X[nrow(X), POS],S); #beggining positions
e<-b+W; #end positions
CHR<-unique(X$CHR) #make error message if this is longer than one.
#X$MAF/N-> X$MAF; #make relative frequencies
ncd1<-do.call(
rbind,mclapply2(1:(length(b)-1), function(i) dplyr::filter(X, POS>=b[i], POS<e[i]) %>% dplyr::select(MAF) %>% dplyr::summarise(Chr=CHR, Beg.Win=b[i], End.Win=e[i],tf0.5=sqrt(sum((MAF-0.5)^2)/length(MAF)),tf0.4=sqrt(sum((MAF-0.4)^2)/length(MAF)), tf0.3=sqrt(sum((MAF-0.3)^2)/length(MAF)), N_SNPs=n()) %>% as.data.table));
return(ncd1)
}

############################################################################################
####### NCD2 ######### NCD2 ########### NCD2 ########## NCD2 ############ NCD2 #############
############################################################################################
#mclapply2

#mclapply wrapper that prints stages of parallel running
#copied from here: http://stackoverflow.com/questions/10984556/is-there-way-to-track-progress-on-a-mclapply
#last mdofied: 14.12.2016
#

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

