library(edgeR)
#set up the experimental design matrix
#metadata<-read.table('metadata2.txt',sep='\t',header=TRUE,row.names=1)

#exp<-metadata$group
#names(exp)<-rownames(metadata)
#exp<-exp[exp %in% c("Kv2010","Kv2014")]
exp<-c(rep("Kv2010",5),rep("Kv2014",5))
names(exp)<-c("X901", "X904", "X910", "X913", "X918", "P408", "P411", "P414", "P422", "P426")
exp<-as.factor(exp)
exp<-relevel(exp,"Kv2010")
#read the data
x<-read.table('/home/tobiaso/abborre/perca/wholetranscript/east_coast/abundance_matrix_eastcoast.txt',sep='\t',header=TRUE,row.names=1)
#x<-read.table('abundance_matrix_eastcoast.txt',sep='\t',header=TRUE,row.names=1)
annotation<-read.table('/home/tobiaso/abborre/perca/wholetranscript/east_coast/annotation_transcripts_eastcoast.txt',sep='\t',header=TRUE,row.names=1,comment.char="",quote="")
#annotation<-read.table('annotation_transcripts_eastcoast.txt',sep='\t',header=TRUE,row.names=1,comment.char="",quote="")


#filter out singletons and low abundant transcripts:
x<-as.matrix(x)
cutoff<-4
n1<-apply(x,1,function(xrow) {return(sum(xrow>0))}) #check if number of nonzero entries is above cutoff
test1<-(n1>=cutoff)
#mean1<-apply(x,1,mean) #test if avg count is above 2.
#test2<-(mean1>=2)
#test <- (test1 & test2)

x<-x[test1,]
ngroup<-length(exp)
head(x)
n<-nrow(x)





#pairwise:

analyze.edgeR<-function(x.sel,exp.sel,name,description){
  n<-nrow(x.sel)
  y <- DGEList(counts=x.sel,group=exp.sel)
  y<-calcNormFactors(y,method="TMM")
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  design.sel<-model.matrix(~exp.sel)
  design.sel<-design.sel[,colnames(design.sel) %in% c("(Intercept)",paste("exp.sel",name,sep=""))]
  rownames(design.sel)<-colnames(x.sel)
  glm<-glmFit(y,design=design.sel)
  glm<-glmLRT(glm,coef=2)#likelihood ratio test
  results<-topTags(glm,n)
  pseudocounts<-y$pseudo.counts
  write.table(pseudocounts,paste('pseudocounts_',name,'vs.control.txt',sep=''),sep='\t',quote=FALSE,col.names=NA)
  outputAll<-results[[1]]
  output<-merge(outputAll,annotation, by.x=0,by.y=0, all.x=TRUE)
  output<-output[order(output$FDR),]
  write.table(output,paste('results_',name,'vs.control.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)
  tmp<-list(results=output,pseudo.counts=pseudocounts)
  return(tmp)
}

out<-analyze.edgeR(x,exp,"Kv2014",annotation)


pdf('pca_kvaddofjarden.pdf',height=5,width=5)
par(mar=c(4,4,4,6))
colors<-c(rep("darkred",5),rep("dodgerblue",5))
par(xpd=NA)
namevec=c("Kv2014")
for (i in 1:length(namevec)){
  pc<-read.table(paste('pseudocounts_',namevec[[i]],'vs.control.txt',sep=''),header=TRUE,row.names=1)
  pca<-prcomp(t(sqrt(x)),scale=TRUE,center=TRUE)
  explained.pc1<-summary(pca)$importance[2,1]
  explained.pc2<-summary(pca)$importance[2,2]
  plot(pca$x[,c(1,2)],bg=colors,pch=21,bty='n',main=paste('PCA ', namevec[[i]], " vs. Kv2010",sep=''),xlab='',ylab='')
  title(xlab=sprintf("PC1(%1.0f%%)",explained.pc1*100),ylab=sprintf("PC2 (%1.0f%%)",explained.pc2*100))
  legend("right",legend=c("Kv2010",namevec[[i]]),pch=21,pt.bg=unique(colors),inset=c(-0.3,-0.05))
  text(pca$x[,c(1,2)],labels=colnames(x),pos=4)
}


dev.off()

