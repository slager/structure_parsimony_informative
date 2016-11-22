setwd("/Users/dave/crows/stacks/missing_data/")
library(data.table)

parseStructureFile <- function(){
  missingdatainteger <- 0   ## 0 for stacks, -9 for pyRAD
  ### Parse Stacks STRUCTURE file
  cat("Parsing",filename,"\n")
  cat("Missing data coded as",missingdatainteger,"\n")
  fread(filename,data.table=F,verbose=F) -> s
  
  ## Get locus labels vector
  ## This part is file-specific.
  readLines(filename,n=3)[2] -> a
  unlist(strsplit(a,split="\t"))[-1] -> sitelabels
  cat(length(sitelabels),"SNPs\n")
  
  ## Get sample labels
  s$V1 -> samplelabels  # Save the original sample labels
  # Relabel samples because duplicate rownames not allowed
  paste0(s$V1,c("a","b")) -> rownames(s) #Rownames useful for troubleshooting
  cat(nrow(s)/2,"samples\n")
  
  ## Remove leading columns so that #columns = #sites
  ## This part is file-specific.
  s[,-c(1,2)] -> s
  
  ## Assign locus labels as the column names
  sitelabels -> colnames(s)
  ## s is now a data.frame
  ## with {#samples*2} rows and {#sites} columns.
  ## Data = various integers
  ## Column names of s are the site labels
  
  assign("s", s, envir = .GlobalEnv)
} # End of function parseStructureFile


retainParsimonyInformative <- function(){
  ### Detect parsimony-informative sites in s
  
  # Create L, a list of vectors for indexing:
  # length(L) = #samples
  # L[[1]] = c(1,2)
  # L[[nSamples]] = c(2*nSamples-1,2*nSamples)
  lapply(seq(1,nrow(s)-1,2),function(x){x:(x+1)}) -> L
  
  #For each site (each column of s):
  apply(s,2,function(a){
    #Sort the diploid unlinked SNP data so that (1,2) is same as (2,1)
    lapply(1:length(L),function(x){sort(a[L[[x]]])}) -> g
    #Remove samples with missing data before detecting pi sites
    g[which(sapply(g,function(x){!all(x==c(missingdatainteger,missingdatainteger))}))] -> g
    #Coerce to character & paste (result = e.g. "11", "12", "22", "13", ...)
    lapply(g,function(x){paste(as.character(x),sep="",collapse="")}) -> g
    unlist(g) -> g  # Convert to vector
    #Check if parsimony informative
    #https://en.wikipedia.org/wiki/Informative_site
    rle(sort(g))$lengths -> sitepatterns
    if (length(sitepatterns[which(sitepatterns >= 2)]) >= 2) TRUE else FALSE #This is value returned by apply() for each site
  }) -> p # A vector of length nsites, indicates if site parsimony-informative
  
  ## Get parsimony-informative data.frame
  s[,which(p)] -> pi
  
  cat(ncol(pi),"parsimony-informative SNPs retained\n")
  assign("s", pi, envir = .GlobalEnv)
} # End function retainParsimonyInformative


retainRandomUnlinkedSNP <- function(){
  
  ## Get 1 random unlinked SNP per locus
  sapply(strsplit(colnames(s),"_"),function(x){x[1]}) -> loci
  unique(loci) -> uloci
  cat(ncol(s),"SNPs at",length(uloci),"loci\n")
  rep(as.character(NA),length(uloci)) -> u
  for (i in 1:length(uloci)){
    colnames(s)[which(loci==uloci[i])] -> current_locus_SNPs
    sample(current_locus_SNPs,1)-> u[i]
  }
  s[,u] -> unlinked
  assign("s", unlinked, envir = .GlobalEnv)
  cat(ncol(unlinked),"random unlinked SNPs retained\n")
} # End of retainRandomUnlinkedSNP


writeOutStructure <- function(suffix){
  ## Write New Structure file
  if (nchar(suffix)==0) {warning("Write terminated: Check writeOutSuffix");stop}
  # Write site names header (Optional)
  write(paste(c("",names(s)),collapse="\t",sep=""),file=paste0(filename,suffix))
  # Write Structure data
  write.table(s,file=paste0(filename,suffix),col.names=FALSE,row.names=samplelabels,append=TRUE,quote=FALSE,sep="\t")
  cat("Wrote to",paste0(filename,suffix))
  } #End of writeOutStructure

###  Where the action happens ###

## Optional single file action
#filename <- "str_files/Cb.s64.r0.05.maf0.02.het0.5.str.tsv"

## Optional "for" loop action
filenames <- c("Cb.s62.r0.05.maf0.02.het0.5.str.tsv","Cb.s64.r0.05.maf0.02.het0.5.str.tsv")
directory <- "str_files/" ## Optional
filenames <- paste0(directory,filenames) ## Optional
for (filename in filenames){
  parseStructureFile() #Reads 'filename', writes 's'
  #retainParsimonyInformative()
  retainRandomUnlinkedSNP() #Reads 's', writes 's'
  writeOutStructure(suffix=".unlinked") #Reads 's', writes 'filename'+'suffix'
}

#### OBSOLETE CODE ####
  
## Summary outputs to console
#ncol(s)          # Total number of sites
#length(unique(sapply(strsplit(colnames(s),"_"),function(x){x[1]}))) # Total number of loci
#ncol(pi)         # Total number of parsimony informative sites
#ncol(upi)        # Number of unlinked parsimony informative sites
#ncol(pi)/ncol(s) # Proportion of sites parsimony informative

## Summary outputs to console OLD
#ncol(s)          # Number of sites
#ncol(pi)         # Number of parsimony informative sites
#ncol(pi)/ncol(s) # Proportion of sites parsimony informative

## Write New Structure file PI
# Write site names header (Optional)
#write(paste(c("",names(pi)),collapse="\t",sep=""),file=paste0(filename,".pi"))
# Write Structure data
#write.table(pi,file=paste0(filename,".pi"),col.names=FALSE,row.names=samplelabels,append=TRUE,quote=FALSE,sep="\t")


#}  #Bracket for entire "for" loop


## Older remnant code (FYI) for initial parse of pyRAD STRUCTURE files

#fread("c88m48p9.str",data.table=F,verbose=T) -> s
#paste0(s$V1,c("a","b")) -> rownames(s)
#s[,-c(1:6)] -> s
#str(s)