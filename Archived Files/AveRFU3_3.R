
#THIS IS A FUNCTION FOR CALCULATING avgRFU (based on Mx estimates)

AveRFU3 = function(popFreq,evidData,refData,hypData,
                  kit,AT,fst,pC,lambda,
                  seed=1,steptol = 1e-6, nDone=3, sig=2 ) {
  library(euroformix) #load euroformix 
  #HELPFUNCTION for obtaining MAC (matching allele counting)
  MAC2 = function(x,y) return(as.numeric(x[1]%in%y)+as.numeric(x[2]%in%y))
  
  #READ DATA:
  #evidData = euroformix::tableReader(evidData_file) #read_excel(evidFn)
  #refData = euroformix::tableReader(refData_file)
  #hypData = euroformix::tableReader(hypData_file)
  #popFreq <- freqImport(popFreq_file)[[1]] #import population freqs.
  locs <- names(popFreq) #loci to consider is restricted to the population database
  
  #LOOP THROUGH EACH HYPOTHESIS
  outcn <- c(colnames(hypData),"loglikhp","loglikhd","LRmle","thetaHp","thetaHd","MAC","nDropout","nAunknown","nRefAlleles","nAlleles","avgRFU")
  bigrestable <- matrix(NA,ncol=length(outcn),nrow=nrow(hypData)) #reset matrix
  colnames(bigrestable) <- outcn
  
  #LOOPING THROUGH EACH HYPOTHESIS SET:
  for(hyprow in 1:nrow(hypData)) { 
    #hyprow=2
    print(paste0("Hypothesis set #",hyprow))
    #hypData[hyprow,]
    evidUse = unlist(hypData[hyprow,1])
    evidUseTable = subset(evidData,subset=`Sample Name`==evidUse) #obtain subset
    evidUseList  = euroformix::sample_tableToList( evidUseTable  )
    
    #Extract hypothesis:
    POI = unlist(hypData[hyprow,2]) #Person of interest (to calc LR for)
    COND = unlist(hypData[hyprow,3]) #conditional references
    if(!is.na(COND) && COND=="") COND = NA
    NOC =  unlist(hypData[hyprow,4]) #number of contributors

    refUse = na.omit(c(POI,COND)) #get refs to consider
    refUseTable = subset(refData,subset=`Sample Name`%in%refUse) #obtain subset
    refUseList  = euroformix::sample_tableToList( refUseTable )[refUse] #ensure same order (POI first if given)
    
    #PREPARE DATA FOR USING euroformix probabilistic functions:
    dat = euroformix::prepareData(evidUseList, refData=refUseList,popFreq=popFreq, threshT=AT )
    
    #BLOCK TO OBTAIN ALLELE STATISTICS
    nAlleles = sum(sapply(dat$samples[[1]] ,function(x) length(x$adata))) #obtain total number of observed alleles
    avgRFU = mean(sapply(dat$samples[[1]] ,function(x) sum(x$hdata))) #obtain average of total PH sum
    #sumRFU = sum(sapply(dat$samples[[1]] ,function(x) sum(x$hdata))) #obtain sum of total PH sum
    #sumRFU/length(popFreq)
    nRefs = length(refUse) #number of references
    
    #compare all references to submix (use only aSTR)
    MACmat <- matrix(0,ncol=1,nrow=nRefs )
    rownames(MACmat) <-refUse 
    for(rr in refUse ) {
      #    rr=refUse[1]
      insind = which(refUse==rr) #index to insert
      #loc=locs[1]
      for(loc in locs) MACmat[insind,1] <- MACmat[insind,1] + MAC2(dat$refData[[loc]][[rr]],dat$samples[[1]][[loc]]$adata)
    }
    totndrops <-  ncol(MACmat)*length(locs)*2 - rowSums(MACmat) #number of dropout per individual (over all samples)
    names(totndrops) = refUse  
    
    # print(totndrops) #number of dropouts
    mactxt = apply(MACmat,1,function(x) paste0(x,collapse="-")) 
    mactxt = paste0(paste0(names(mactxt),"=",mactxt),collapse="/")
    
    dotxt = sapply(totndrops,function(x) paste0(x,collapse="-")) 
    dotxt = paste0(paste0(names(dotxt),"=",dotxt),collapse="/")
    
    #obtain number of tot alleles/heights for unknows
    nAu = NA #number of alleles not belonging to any referencess
    hAu = NA #Height of alleles not belonging to any refs    
    nA_refs = rep(NA,nRefs)  #number of (unique) alleles for references
    nAutxt = paste0("nU=",nAu,"/hU=",hAu)
    nArtxt = paste0(paste0(refUse,"=",nA_refs),collapse="/")
    #END NEW BLOCK
    
    #VIEW DATA:
    # plotEPG2(evidUseList,refUseList ,kit=kit,AT=AT)
    modtype = 1:3 #1="NO STUTTER", 2="ONLY BW STUTTER", 3="BOTH STUTTER"
    for(mod in modtype) { #loop through each model type
      
      #SELECT MODEL FOR PROBABILISTIC CALCULATION:
      kit0 <- NULL
      xiFW <- xiBW <- 0 #no stutter 
      if(mod%in%2:3) xiBW <- NULL #use BW stutter
      if(mod%in%3) xiFW <- NULL #use FW stutter
      if(!is.null(kit)) kit0 <- kit #use degradation
      
      #CALCULATE HYPS FOR DEFINED HYPOTHESES:
      nCOND = length(refUse) #number of references to condition on 
      condHp <- condHd <- nCOND:1 #POI is set last
      knownRef <- NULL #known non-contributors under hd
      loghp <- thhatp <- NA
      if(!is.na(POI)) {
        condHd[1] = 0 #
        knownRef = 1 #first ref is POI
        
        #Under hp:
        tryCatch({
          hpfit <- contLikMLE(nC=NOC,samples=dat$samples,refData=dat$refData,popFreq=dat$popFreq,threshT=AT,verbose=FALSE,xi=xiBW,xiFW=xiFW,kit=kit0,condOrder=condHp,prC=pC,lambda=lambda,fst=fst,nDone=nDone,steptol=steptol,seed=seed)
          loghp <- hpfit$fit$loglik 
          thhatp <- paste0(signif(hpfit$fit$thetahat2,sig),collapse="/")
        },error=function(e) print(e))
      }
      if(is.na(loghp) || is.infinite(loghp)) break #stop loop if model not converging

      #UNRELATED under Hd::
      hdfit <- contLikMLE(nC=NOC,samples=dat$samples,refData=dat$refData,popFreq=dat$popFreq,threshT=AT,verbose=FALSE,xi=xiBW,xiFW=xiFW,kit=kit0,condOrder=condHd,prC=pC,lambda=lambda,fst=fst,knownRef=knownRef,nDone=nDone,steptol=steptol,seed=seed)
      loghd <- hdfit$fit$loglik 
      thhatd <- paste0(signif(hdfit$fit$thetahat2,sig),collapse="/") #SHOW PARAMS UNDER UNRELATEDNESS
      if(is.infinite(loghd)) break #stop loop if model not converging
      
      #Calculate LR
      LRmle <- exp(loghp - loghd)  
    
      #store results:
      bigrestable[hyprow,] = c(unlist(hypData[hyprow,]),loghp,loghd, LRmle, thhatp,thhatd,mactxt,dotxt,nAutxt,nArtxt ,nAlleles,avgRFU)
    } #end foor each model fit type
  } #end for each hypotheses
write.table(bigrestable, file="ResultsFile", sep=";", row.names=FALSE)
  return(bigrestable)
}