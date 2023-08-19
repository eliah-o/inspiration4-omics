library(tidyverse)

### ARGUMENTS
# input directory
# output directory
# database (bacterial/viral)
# thres
# uthres

args = commandArgs(trailingOnly=TRUE)


inputdir = args[[1]]
outputdir = args[[2]]
database = args[[3]]
thres = args[[4]]
uthres = args[[5]]
namemap = args[[6]]
dtype = args[[7]]


#inputdir = "prof_bac_xtree/metag/"#args[[1]]

#outputdir = "xtree_merged/"#args[[2]]
#database = "gtdb"#args[[3]]
#thres = .01#args[[4]]
#uthres = 0.005#args[[5]]
#namemap = "xtree_all_db_mapping"#args[[6]]
#dtype = "metagenomics"#args[[7]]

#inputdir = "prof_vir_xtree/metag/"#args[[1]]
#outputdir = "xtree_merged/"#args[[2]]
#database = "genbank-viral"#args[[3]]
#thres = .01#args[[4]]
#uthres = 0.005#args[[5]]
#namemap = "viral_tax_map.tsv"#args[[6]]
#dtype = "metagenomics"#args[[7]]


covs = data.frame()
covsU = data.frame()
# covsR = data.frame()
files = dir(inputdir,pattern = "*.cov")
for (file in files) {
  t = read.delim(paste(inputdir,file,sep='/'),T,row=1,as.is=T)
  sname = gsub(".cov.*","",file)
  covs[rownames(t),sname] = t[,"Proportion_covered"]
  covsU[rownames(t),sname] = t[,"Unique_proportion_covered"]
  #covsR[rownames(t),sname] = t[,"Reads_covered"]
}
covs[is.na(covs)]=0
covsU[is.na(covsU)]=0
# covsR[is.na(covsR)] = 0

otus = data.frame()
files = dir(inputdir,pattern = "*.ref")
for (file in files) {
  tryCatch({
  t = read.delim(paste(inputdir,file,sep='/'),F,row=1,as.is=T)
  },
  error=function(cond){
  print('No lines in file')
  })
  otus[rownames(t),gsub(".ref.*","",file)] = t
}
otus[is.na(otus)]=0

# Prepare the coverage + ref cross-maps with taxonomy
if(namemap != 'None'){
  tkey= read.delim(namemap,row=1,head=F,quote="")
  covs = covs[rownames(covs) %in% rownames(tkey), ]
  covsU = covsU[rownames(covsU) %in% rownames(tkey), ]
  taxaconv = tkey[rownames(covs),1] 
  otus = otus[rownames(otus) %in% rownames(tkey), ]
  taxaconvR = tkey[rownames(otus),1]
}
if(namemap == 'None'){
  taxaconv = rownames(covs)
  taxaconvR = rownames(otus)
}

orig.abund = colSums(otus)

mask = (covsU <= as.numeric(uthres) | covs <= as.numeric(thres)) #& covs <= hthres
unmasked = rowsum(0+!mask,taxaconv)

tax.tm = rowsum(otus,taxaconvR)
tax.tm[!unmasked[rownames(tax.tm),]]=0
tax.tm = tax.tm[rowSums(tax.tm) > 0,]
tax.tm["Unknown",]=orig.abund - colSums(tax.tm)
tax.cm = rowsum(otus,taxaconvR)[rownames(tax.tm)[-nrow(tax.tm)],]


saveRDS(tax.tm,paste(outputdir,'/',database,'_',as.character(thres),'_',as.character(uthres),'_',dtype,'_','species_counts.rds',sep=''))

write.table(tax.tm,paste(outputdir,'/',database,'_',as.character(thres),'_',as.character(uthres),'_',dtype,'_','species_counts.tsv',sep=''),F,F,'\t')
tax.tmr=sweep(tax.tm,2,colSums(tax.tm),'/')
write.table(tax.tmr,paste(outputdir,'/',database,'_',as.character(thres),'_',as.character(uthres),'_',dtype,'_','species_ra.tsv',sep=''),F,F,'\t')
write.table(tax.cm,paste(outputdir,'/',database,'_',dtype,'_','species_rawcounts.tsv',sep=''),F,F,'\t')

if(database == 'GTDB'){
lv = rev(c("d","p","c","o","f","g","s","t"))
for (l in 1:(length(lv)-1)) {
  t.cm = rowsum(tax.cm,sub(paste0(";",lv[l],"__.*"),"",rownames(tax.cm)))
  write.table(t.cm,paste(outputdir,'/',database,'_',lv[l+1],'_rawcounts.tsv',sep=''),F,F,'\t')
  t.tm = rowsum(tax.tm,sub(paste0(";",lv[l],"__.*"),"",rownames(tax.tm)))
  write.table(t.tm,paste(outputdir,'/',database,'_',as.character(thres),'_',as.character(uthres),'_',lv[l+1],'_',dtype,'_counts.tsv',sep=''),F,F,'\t')
  t.tmr=sweep(t.tm,2,colSums(tax.tm),'/')
  write.table(t.tmr,paste(outputdir,'/',database,'_',as.character(thres),'_',as.character(uthres),'_',lv[l+1],'_',dtype,'_ra.tsv',sep=''),F,F,'\t')
  }
}


if(database == 'genbank-viral'){
  lv = rev(c("p","c","o","f","g","s"))
  for (l in 1:(length(lv)-1)) {
    t.cm = rowsum(tax.cm,sub(paste0(";",lv[l],"__.*"),"",rownames(tax.cm)))
    write.table(t.cm,paste(outputdir,'/',database,'_',lv[l+1],'_rawcounts.tsv',sep=''),F,F,'\t')
    t.tm = rowsum(tax.tm,sub(paste0(";",lv[l],"__.*"),"",rownames(tax.tm)))
    write.table(t.tm,paste(outputdir,'/',database,'_',as.character(thres),'_',as.character(uthres),'_',lv[l+1],'_',dtype,'_counts.tsv',sep=''),F,F,'\t')
    t.tmr=sweep(t.tm,2,colSums(tax.tm),'/')
    write.table(t.tmr,paste(outputdir,'/',database,'_',as.character(thres),'_',as.character(uthres),'_',lv[l+1],'_',dtype,'_ra.tsv',sep=''),F,F,'\t')
  }
}




