

namevec=c("EE2","Mix","Cd","BNF")
output.list<-NULL
for (i in 1:length(namevec)){
  output.list[[i]]<-read.table(paste('results_',namevec[[i]],'vs.control.txt',sep=''),sep='\t',header=TRUE,row.names=1,quote="",comment.char="")
}
names(output.list)<-namevec
#gomap<-read.table('wholetranscript/west_coast/go_bp__westcoast.txt',sep='\t',quote="",comment.char="")
#goannot<-read.table('all_GOterms_new.txt',sep='\t',quote="",comment.char="",header=TRUE)
#goannot<-goannot[,1:3]
fdr.cutoff<-0.1
uniprot<-read.table('wholetranscript/west_coast/annotation_transcripts_uniprot_westcoast.txt',sep='\t',header=TRUE)
#up for all
#gomatrix<-go.analysis2(output.list,namevec,gomap,fdr.cutoff=0.05,direction="up")
#totalNumberofGenes<-nrow(output.list[[length(output.list)]])#background
for(i in 1:length(output.list)){
  x1<-output.list[[i]]
  siggenes<-unique(rownames(x1)[x1$FDR<fdr.cutoff])
  upgenes<-merge(as.matrix(siggenes),uniprot,by.x=1,by.y=1,all.x=TRUE)
  upgenes.unique<-unique(upgenes$Annotation_uniprot[!upgenes$Annotation_uniprot %in% ''])
  write.table(as.matrix(upgenes.unique),paste("GOterms/uniprot_",fdr.cutoff,"_",namevec[[i]],".txt",sep=''),sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
}
#gomap<-read.table('wholetranscript/west_coast/go_mf__westcoast.txt',sep='\t',quote="",comment.char="")
#goannot<-read.table('all_GOterms_new.txt',sep='\t',quote="",comment.char="",header=TRUE)
#goannot<-goannot[,1:3]

#for(i in 1:length(output.list)){
#  x1<-output.list[[i]]
#  siggenes<-unique(rownames(x1)[x1$FDR<fdr.cutoff & x1$logFC>0])
#  go.pvalues.up<-hygtest2(siggenes,totalNumberofGenes,gomap)
#  go.pvalues.up.sig<-as.data.frame(go.pvalues.up)
#  go.pvalues.up.sig<-go.pvalues.up.sig[go.pvalues.up.sig$p<0.05,]
  #go.pvalues.up.sig<-go.pvalues.up.sig[order(go.pvalues.up.sig$p),]
#  goout<-merge(go.pvalues.up.sig,goannot,by.x=0,by.y=1,all.x=TRUE)
#  write.table(goout,paste("GOterms/go_results_up_mf_",fdr.cutoff,"_",namevec[[i]],".txt",sep=''),sep='\t',quote=FALSE,col.names=NA)
#  siggenes<-unique(rownames(x1)[x1$FDR<fdr.cutoff & x1$logFC<0])
#  go.pvalues.down<-hygtest2(siggenes,totalNumberofGenes,gomap)
#  go.pvalues.down.sig<-as.data.frame(go.pvalues.down)
#  go.pvalues.down.sig<-go.pvalues.down.sig[go.pvalues.down.sig$p<0.05,]
#  #go.pvalues.down.sig<-go.pvalues.down.sig[order(go.pvalues.down.sig$p),]
#  goout<-merge(go.pvalues.down.sig,goannot,by.x=0,by.y=1,all.x=TRUE)
#  write.table(goout,paste("GOterms/go_results_down_mf_",fdr.cutoff,"_",namevec[[i]],".txt",sep=''),sep='\t',quote=FALSE,col.names=NA)
#}
