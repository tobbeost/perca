hygtest2<-function(siggenes,totalNumberofGenes,gomap){
  go<-gomap #list of all genes and GO-terms
  ugo<-unique(go[,2]) #unique GO-terms
  pgo<-NULL
  allgenes<-totalNumberofGenes  # number of unique genes
  for (i in 1:length(ugo)){
    allgenesGO<-go[go[,2] %in% ugo[[i]],1]
    siggenesGO<-length(siggenes[siggenes %in% allgenesGO])
    m<-length(allgenesGO)
    n<-allgenes-length(allgenesGO)
    k<-length(siggenes)
    p<-phyper(siggenesGO-1, k, allgenes-k, m,lower.tail=FALSE) #calculate the hypergeometric probability
    pgo<-rbind(pgo,c(p,siggenesGO,m,n,k))
  }
  rownames(pgo)<-ugo
  colnames(pgo)<-c("p","siggenesGO","allgenesGO","allgenes-allgenesGO","siggenes")
  return(pgo)
}

go.analysis2<-function(output.list,labellist,gomap,fdr.cutoff=0.05,direction="all"){
  gomatrix<-unique(gomap[,2])
  totalNumberofGenes<-nrow(output.list[[length(output.list)]])#background
  for(i in 1:length(output.list)){
    x1<-output.list[[i]]
    if (direction=="all"){
      siggenes<-unique(rownames(x1)[x1$FDR<fdr.cutoff])
    }
    if (direction=="up"){
      siggenes<-unique(rownames(x1)[x1$FDR<fdr.cutoff & x1$logFC>0])
    }
    if (direction=="down"){
      siggenes<-unique(rownames(x1)[x1$FDR<fdr.cutoff & x1$logFC<0])
    }
    go.pvalues<-hygtest2(siggenes,totalNumberofGenes,gomap)
    #go.pvalues<-go.pvalues[order(names(go.pvalues))]
    gomatrix<-merge(gomatrix,go.pvalues)
  }
  colnames(gomatrix)<-labellist
  return(gomatrix)
}

namevec=c("EE2","Mix","Cd","BNF")
output.list<-NULL
for (i in 1:length(namevec)){
  output.list[[i]]<-read.table(paste('results_',namevec[[i]],'vs.control.txt',sep=''),sep='\t',header=TRUE,row.names=1,quote="",comment.char="")
}
names(output.list)<-namevec
gomap<-read.table('wholetranscript/west_coast/go_bp__westcoast.txt',sep='\t',quote="",comment.char="")
goannot<-read.table('all_GOterms_new.txt',sep='\t',quote="",comment.char="",header=TRUE)
goannot<-goannot[,1:3]
fdr.cutoff<-0.05
#up for all
#gomatrix<-go.analysis2(output.list,namevec,gomap,fdr.cutoff=0.05,direction="up")
totalNumberofGenes<-nrow(output.list[[length(output.list)]])#background
for(i in 1:length(output.list)){
  x1<-output.list[[i]]
  siggenes<-unique(rownames(x1)[x1$FDR<fdr.cutoff & x1$logFC>0])
  go.pvalues.up<-hygtest2(siggenes,totalNumberofGenes,gomap)
  go.pvalues.up.sig<-as.data.frame(go.pvalues.up)
  go.pvalues.up.sig<-go.pvalues.up.sig[go.pvalues.up.sig$p<0.05,]
  #go.pvalues.up.sig.ordered<-go.pvalues.up.sig[order(go.pvalues.up.sig$p),]
  goout<-merge(go.pvalues.up.sig,goannot,by.x=0,by.y=1,all.x=TRUE)
  write.table(goout,paste("GOterms/go_results_up_bp_",fdr.cutoff,"_",namevec[[i]],".txt",sep=''),sep='\t',quote=FALSE,row.names=FALSE)
  siggenes<-unique(rownames(x1)[x1$FDR<fdr.cutoff & x1$logFC<0])
  go.pvalues.down<-hygtest2(siggenes,totalNumberofGenes,gomap)
  go.pvalues.down.sig<-as.data.frame(go.pvalues.down)
  go.pvalues.down.sig<-go.pvalues.down.sig[go.pvalues.down.sig$p<0.05,]
  #go.pvalues.down.sig.ordered<-go.pvalues.down.sig[order(go.pvalues.down.sig$p),]
  goout<-merge(go.pvalues.down.sig,goannot,by.x=0,by.y=1,all.x=TRUE)
  write.table(goout,paste("GOterms/go_results_down_bp_",fdr.cutoff,"_",namevec[[i]],".txt",sep=''),sep='\t',quote=FALSE,row.names=FALSE)
}
gomap<-read.table('wholetranscript/west_coast/go_mf__westcoast.txt',sep='\t',quote="",comment.char="")
goannot<-read.table('all_GOterms_new.txt',sep='\t',quote="",comment.char="",header=TRUE)
goannot<-goannot[,1:3]

for(i in 1:length(output.list)){
  x1<-output.list[[i]]
  siggenes<-unique(rownames(x1)[x1$FDR<fdr.cutoff & x1$logFC>0])
  go.pvalues.up<-hygtest2(siggenes,totalNumberofGenes,gomap)
  go.pvalues.up.sig<-as.data.frame(go.pvalues.up)
  go.pvalues.up.sig<-go.pvalues.up.sig[go.pvalues.up.sig$p<0.05,]
  #go.pvalues.up.sig<-go.pvalues.up.sig[order(go.pvalues.up.sig$p),]
  goout<-merge(go.pvalues.up.sig,goannot,by.x=0,by.y=1,all.x=TRUE)
  write.table(goout,paste("GOterms/go_results_up_mf_",fdr.cutoff,"_",namevec[[i]],".txt",sep=''),sep='\t',quote=FALSE,col.names=NA)
  siggenes<-unique(rownames(x1)[x1$FDR<fdr.cutoff & x1$logFC<0])
  go.pvalues.down<-hygtest2(siggenes,totalNumberofGenes,gomap)
  go.pvalues.down.sig<-as.data.frame(go.pvalues.down)
  go.pvalues.down.sig<-go.pvalues.down.sig[go.pvalues.down.sig$p<0.05,]
  #go.pvalues.down.sig<-go.pvalues.down.sig[order(go.pvalues.down.sig$p),]
  goout<-merge(go.pvalues.down.sig,goannot,by.x=0,by.y=1,all.x=TRUE)
  write.table(goout,paste("GOterms/go_results_down_mf_",fdr.cutoff,"_",namevec[[i]],".txt",sep=''),sep='\t',quote=FALSE,col.names=NA)
}
