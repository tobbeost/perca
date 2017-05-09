library(data.table)
library(fgsea)
gomap<-read.table('wholetranscript/west_coast/go_mf__westcoast.txt',sep='\t',quote="",comment.char="")
loadGSEA<-function(gomap){
  ugo<-unique(gomap[,2])
  npath<-length(ugo)
  go<-list()
  for(i in 1:npath){
    name<-ugo[i]
    index<-which(gomap[,2]==ugo[i])
    go[[i]]<-gomap[index,1]
  }
  names(go)<-ugo
  return(go)
}
gsc<-loadGSEA(gomap)
a<-commandArgs(TRUE)[1]
x1<-read.table(paste('results_',a,'.txt',sep=''),sep='\t',header=TRUE,quote="",comment.char="")
pvals<- -log(x1$PValue)
pvals[x1$logFC<0]<- -pvals[x1$logFC<0]
names(pvals)<-x1$Row.names
pvals<-sort(pvals)
gsares<-fgsea(pathways=gsc,stats=pvals,nperm=100000,nproc=20)
gsares<-gsares[order(gsares$padj),1:(dim(gsares)[2]-1)]
fwrite(gsares,paste('GSEA_results_',a,'_fgsea_GOmf2.txt',sep=''),sep='\t',quote=FALSE)