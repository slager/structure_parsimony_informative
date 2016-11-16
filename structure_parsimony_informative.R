#setwd("/Users/dave/crows/stacks/missing_data/")

library(data.table)

### Parse Stacks STRUCTURE file

filename <- "str_files/b.05.random.structure copy.tsv"
missingdatainteger <- 0   ## 0 for stacks, -9 for pyRAD

## Optional "for" loop
#filenames <- c("b.05.all.structure.tsv","b.05.random.structure.tsv","b.75.all.structure.tsv","b.75.random.structure.tsv","c.05.random.structure.tsv","c.75.all.structure.tsv","c.75.random.structure.tsv")
#filenames <- paste0("str_files/",filenames)
#for (filename in filenames){  #Also a commented bracket at the end

fread(filename,data.table=F,verbose=T) -> s

## Get locus labels vector
## This part is file-specific.
readLines(filename,n=3)[2] -> a
unlist(strsplit(a,split="\t"))[-1] -> sitelabels

## Get sample labels
s$V1 -> samplelabels  # Save the original sample labels
# Relabel samples because duplicate rownames not allowed
paste0(s$V1,c("a","b")) -> rownames(s) #Rownames useful for troubleshooting

## Remove leading columns so that #columns = #sites
## This part is file-specific.
s[,-c(1,2)] -> s

## Assign locus labels as the column names
sitelabels -> colnames(s)

## s is now a data.frame
## with {#samples*2} rows and {#sites} columns.
## Data = various integers
## Column names of s are the site labels

### Detect parsimony-informative sites

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
  if (length(sitepatterns[which(sitepatterns >= 2)]) >= 2) 1 else 0 #This is value returned by apply() for each site
}) -> p # A vector of length nsites, indicates if site parsimony-informative

## Get parsimony-informative data.frame
s[,which(as.logical(p))] -> pi

## Summary outputs to console
ncol(s)          # Number of sites
ncol(pi)         # Number of parsimony informative sites
ncol(pi)/ncol(s) # Proportion of sites parsimony informative

## Write New Structure file
# Write site names header (Optional)
write(paste(c("",names(pi)),collapse="\t",sep=""),file=paste0(filename,".pi"))
# Write Structure data
write.table(pi,file=paste0(filename,".pi"),col.names=FALSE,row.names=samplelabels,append=TRUE,quote=FALSE,sep="\t")

#}  Bracket for entire "for" loop


## Older remnant code (FYI) for initial parse of pyRAD STRUCTURE files

#fread("c88m48p9.str",data.table=F,verbose=T) -> s
#paste0(s$V1,c("a","b")) -> rownames(s)
#s[,-c(1:6)] -> s
#str(s)