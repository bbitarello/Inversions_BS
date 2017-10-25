###############################################################################
#NCV scan for 1000g data
#dataAuthor: Barbara Bitarello
#BitarelloCreation: 17.12.2013
#Modified: 28.02.2014 by Cesare de Filippo
#Modified: 21.10.2014 by Barbara Bitarello
#BitarelloNote: attempting to fix two issues:run NCV for all pops at once and
#save window coordinates in output.
#Done
#Last Modified: 6.10.2016 by Barbara Bitarello 
###############################################################################
#NCV is calculated regardless of SNP density. Filtering per snp density should
#happen downstream in case we decide to change the threshold.

#in the future, modify this to have all pops run at once.
NCV.scan4<-function(INPUT.N, pop='YRI',FDs=TRUE, FD.N, WIN, SNP=TRUE) {  
      ##INPUT.N : input data // ##FDs: fixed differences (human vs chimp reference)
        WIN[1,1]->beg; WIN[1,2]->end
	cat('Beg win:', beg, 'End win:', end,'.And pop is', pop,"\n")
        if(pop=='PUR'){n<-88}else{n<-100}  #all pops except PUR have 50 unrelated individuals.
        if(FDs==TRUE){ as.data.frame(FD.N)->z; #if we want to calculate NCD2 and there is an actual FD file for this 3 kb window; do
        nifds<-dim(z)[1]}else{nifds<-0} ; #list of FDs between human and chimp
	if(SNP==TRUE){
        nisnps<-dim(INPUT.N)[1]
	}else{ #number of SNPs in INPUT.N
	nisnps<-0}
	
	cat('I have', nisnps, 'SNPs and', nifds, 'FDs in this input file\n')

       if(SNP==TRUE){
		y2 <- as.data.frame(cbind(counts=as.numeric(INPUT.N[,pop])/n, pos=as.numeric(INPUT.N[,2]), ref=INPUT.N[,4], alt=INPUT.N[,5]), stringsAsFactors=F); #
       		y2[,1]<-as.numeric(y2[,1]);y2[,2]<-as.numeric(y2[,2]); y3<-y2;
             #y3[,1] <- sapply(y2[,1], function(x) if (x>0.5){x<-1-x} else{x<-x})  #use minor allele frequency.
                #up until this point we have all original SNPs in y2. Now we
            #need to filter SNPs and FDs.
	          if(FDs==T){
		#real.snps <- y3[which(y3[,1] > 0),] 
        	tmp<-which(z$pos %in% y2$pos) #positions present in SNP and FD file. Store them and then check what they are

		cat('There are', length(tmp),'positions present in both the SNP and FD\n')

                    tmp.vec<-NA;
                   if(length(tmp)){for (i in 1:length(tmp)){  #for each of these positions
    			temp.pos<-z$pos[tmp[i]]; if(toupper(z$chimp[tmp[i]]) == toupper(y2[y2$pos==temp.pos,4])){#upper and lower case match.  
                                tmp.vec<-c(tmp.vec,tmp[i])}}  #store the positiuons which wuill be eliminated.    
    #if chimp in FD is != Alt in SNP, this is a SNP and also a FD. Keep both. 
    #no need to put any conditional statements for this
			    tmp.vec<-tmp.vec[-1]  #eliminate the NA
	
				    if(length(tmp.vec)){z2<-z[-tmp.vec,]}   #FD without the fake FDs.
				        else{z2<-z} } else{z2<-z; cat('Nothing changed in the FD input\n')}  #check for positions present in SNPs that have f=1
                        #for alt allele and are absent from. FDs. Include them as FDs and exclude from SN                    
					tmp.vec2<-0; tmp.vec3<-NA; tmp2<-which(y2$counts ==1);
					cat('There are', length(tmp2), 'fixed alternate alleles in this pop and window\n');
				            if(length(tmp2)){for (j in 1:length(tmp2)){  #for each of these positions
				                temp.pos2<-y2$pos[tmp2[i]]; tmp.vec3<-c(tmp.vec3, temp.pos2); #exclude from SNPs
					              if(!(temp.pos2 %in% z2$pos)){  #if this SNP position is not present in the FD
					               tmp.vec2<-tmp.vec2+1}} # count number of FDs which should be included in NCV fd . ATTENTION: the FDinput bed will not be changed
						tmp.vec3<-tmp.vec3[-1]  #eliminate NA
					        fxdlen<-dim(z2)[1]+tmp.vec2 #if it is 0 no FDs will be added
			                    if(length(tmp.vec3)){#if this vector has at least one position
                       				real.snps2<-subset(y2, !(pos %in% tmp.vec3))}
		                            if(length(tmp.vec3)==0){ #if there are no positions in tmp.vec3, all stays the same.
		                                real.snps2<-y2}}
					else{real.snps2<-y2;fxdlen<-dim(z2)[1]} #very important detail
                            		real.snps3<-real.snps2[which(real.snps2$counts!=0),] #eliminate remaining snps with f=0 
				cat('So in the end we have', nrow(real.snps3), 'SNPs and', fxdlen, 'FDs left for this pop and window\n');
                                polsites <- dim(real.snps3)[1] ;  #the 'real' number of SNPs used in NCV calculation.
				if(polsites>0){
                                tp<-as.numeric(c(real.snps3$counts,rep(0,fxdlen)));
				tp<-sapply(tp, function(x) if (x>0.5){x<-1-x} else{x<-x})} #use MAF
				if(polsites==0 & fxdlen>0){tp<-rep(0, fxdlen)}
				if(polsites==0 & fxdlen==0){ tp<-NA};cat(tp);}
                ###############################################NCV without
                ###FD############################################
                    if(FDs==F){fxdlen<-0;nifds<-0;real.snps4 <- y2[which(y2[,1] > 0 & y2[,1]<1),];polsites<-dim(real.snps4)[1];
			if(polsites>0){tp<-as.numeric(real.snps4$counts); tp<-sapply(tp, function(x) if (x>0.5){x<-1-x} else{x<-x})}
			if(polsites==0){tp<-NA}
					}#close the if(FD==F)
			}#close the if(SNP==T)
		if(SNP==FALSE){cat('i have no SNPs but I do have Fds in this window\n');
			nrow(z)-> fxdlen;tp<-rep(0,fxdlen);polsites<-0;}
			ncdf1<-sqrt(sum((tp-0.1)^2)/(polsites+fxdlen)); ncdf2<-sqrt(sum((tp-0.2)^2)/(polsites+fxdlen));ncdf3<-sqrt(sum((tp-0.3)^2)/(polsites+fxdlen));
			ncdf4<-sqrt(sum((tp-0.4)^2)/(polsites+fxdlen)); ncdf5<-sqrt(sum((tp-0.5)^2)/(polsites+fxdlen));
                ################################################################fxdlen 
final<- data.frame(Beg.Win=beg, End.Win=end, Initial_seg_sites=nisnps, Initial_fds_sites=nifds, NCDf1=ncdf1,NCDf2=ncdf2, NCDf3=ncdf3,NCDf4=ncdf4, NCDf5=ncdf5,Nr.SNPs=polsites, Nr.FDs=fxdlen);
                
return(final);
}

