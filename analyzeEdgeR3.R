library(edgeR)
#set up the experimental design matrix
metadata<-read.table('metadata2.txt',sep='\t',header=TRUE,row.names=1)

exp<-metadata$group
names(exp)<-rownames(metadata)
exp<-exp[exp %in% c("control","EE2","Mix","Cd","BNF")]
ee2<-vector(length=length(exp))
ee2[exp %in% c('EE2','Mix')]<-1
cd<-vector(length=length(exp))
cd[exp %in% c('Cd','Mix')]<-1
bnf<-vector(length=length(exp))
bnf[exp %in% c('BNF','Mix')]<-1
design<-model.matrix(~cd+ee2+bnf+cd:ee2:bnf)
rownames(design)<-names(exp)

#read the data
x<-read.table('/home/tobiaso/abborre/perca/wholetranscript/west_coast/abundance_matrix_westcoast.txt',sep='\t',header=TRUE,row.names=1)
annotation<-read.table('/home/tobiaso/abborre/perca/wholetranscript/west_coast/annotation_transcripts_westcoast.txt',sep='\t',header=TRUE,row.names=1,comment.char="",quote="")
x<-x[,c(16:23,1:15)]

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


#DE analysis
y <- DGEList(counts=x,group=exp)
y<-calcNormFactors(y,method="TMM")
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)

glm<-glmFit(y,design=design)
glm1<-glmLRT(glm,coef=2:5)#likelihood ratio test
glm2<-glmLRT(glm,coef=2)#likelihood ratio test
glm3<-glmLRT(glm,coef=3)#likelihood ratio test
glm4<-glmLRT(glm,coef=4)#likelihood ratio test
glm5<-glmLRT(glm,coef=5)#likelihood ratio test
results<-topTags(glm1,n)
results2<-topTags(glm2,n)
results3<-topTags(glm3,n)
results4<-topTags(glm4,n)
results5<-topTags(glm5,n)
pseudocounts<-y$pseudo.counts
write.table(pseudocounts,'pseudocounts_westcoast.txt',sep='\t',quote=F)
outputAll<-results[[1]]
#dom<-split.domNames(rownames(outputAll))
#outputAll$Name<-dom$domains
#outputAll$Subbin<-dom$subbins

output<-merge(outputAll,annotation, by.x=0,by.y=0, all.x=TRUE)
output<-output[order(output$FDR),]
write.table(output,'results_ANOVA_glm.txt',sep='\t',quote=FALSE,row.names=FALSE)
output2<-merge(results2[[1]],annotation, by.x=0,by.y=0, all.x=TRUE)
output2<-output2[order(output2$FDR),]
write.table(output2,'results_cd_glm.txt',sep='\t',quote=FALSE,row.names=FALSE)
output3<-merge(results3[[1]],annotation, by.x=0,by.y=0, all.x=TRUE)
output3<-output3[order(output3$FDR),]
write.table(output3,'results_ee2_glm.txt',sep='\t',quote=FALSE,row.names=FALSE)
output4<-merge(results4[[1]],annotation, by.x=0,by.y=0, all.x=TRUE)
output4<-output4[order(output4$FDR),]
write.table(output4,'results_bnf_glm.txt',sep='\t',quote=FALSE,row.names=FALSE)
output5<-merge(results5[[1]],annotation, by.x=0,by.y=0, all.x=TRUE)
output5<-output5[order(output5$FDR),]
write.table(output5,'results_mix_glm.txt',sep='\t',quote=FALSE,row.names=FALSE)
cutoff<-0.05
sig1<-rownames(outputAll)[outputAll$FDR<cutoff]
sig1.up<-rownames(outputAll)[outputAll$FDR<cutoff & outputAll$logFC>0]
sig1.down<-rownames(outputAll)[outputAll$FDR<cutoff & outputAll$logFC<0]
bins=c(length(rownames(outputAll)),length(sig1),length(sig1.up),length(sig1.down))


#plot pca
plot.pca2<-function(x,main="PCA",groups,invert.y=FALSE,plot.legend=FALSE,plot.text=FALSE,show.percentage=FALSE,...){
  x<-as.matrix(x)
  head(x)
  n<-nrow(x) #number of rows in the matrix
  x.sqrt<-sqrt(x)
  #x.log<-log(x+1)
  #groups<-c("0","316","100","0","0","0.316","3.16","31.6","316","0","316","1000","3.16","1","3.16","10")#specify which group each sample belong to
  groups<-as.factor(groups)
  ugroups<-unique(groups)
  library(RColorBrewer)
  col.tmp<-brewer.pal(9,'Set1')
  pch.tmp<-rep(21:25,length=9)
  #col.tmp<-c("dodgerblue","red3","forestgreen","darkorange1","violetred1","skyblue3","slategray2","thistle1","burlywood4")
  pca<-prcomp(t(x.sqrt),center=TRUE,scale=TRUE)
  #summary(pca)
  explained.pc1<-summary(pca)$importance[2,1]
  explained.pc2<-summary(pca)$importance[2,2]
  par(xpd=NA)
  #postscript('pca_taxa_bacteria.eps',height=6,width=6,font="Helvetica")
  if(invert.y==TRUE){
    pca$x[,2]<- -pca$x[,2]
  }
  groupsfactor<-as.numeric(as.factor(groups))
  #pch.tmp=pch.tmp[c(1:5)]
  pch=pch.tmp[groupsfactor]
  #col.tmp<-col.tmp[c(1:5)]
  col<-col.tmp[groupsfactor]
  
  plot(pca$x[,c(1,2)],col="black",bg=col.tmp[groupsfactor],pch=pch,cex=1.3,bty='n',main=main)
  axis(1,lwd=3,cex.lab=0.9)
  if (show.percentage==TRUE){
    title(xlab=sprintf("PC1(%1.0f%%)",explained.pc1*100),ylab=sprintf("PC2 (%1.0f%%)",explained.pc2*100))
  }else{
    title(xlab='PC1',ylab='PC2')
  }
  axis(2,lwd=3)
  if (plot.text==TRUE){
    text(pca$x,labels=groups,pos=c(1,3,2,3,2,1,2,2,2,4,2,1,3,2,2,2))
  }
  
  if(plot.legend==TRUE){
    plot.new()
    groupsfactor<-as.numeric(as.factor(groups))
    legend2<-ugroups[order(unique(groupsfactor))]
    #legend2<-legend2[c(1:3,7,4,8,5,9,6)]
    #pch2<-pch.tmp[c(1:3,7,4,8,5,9,6)]
    #col2<-col.tmp[c(1:3,7,4,8,5,9,6)]
    legend("topleft",legend=legend2,pch=pch,pt.bg=col)
  }
  return(pca)
}
par(mfrow=c(1,2))
pca<-plot.pca2(x,groups=exp,plot.legend=TRUE,plot.text=TRUE,show.percentage=TRUE)
colors<-c(rep("red",5),rep("green",5),rep("blue",5),rep("yellow",4),rep("black",4))
pdf('pca.pdf',height=5,width=10)
par(mfrow=c(1,2),mar=c(4,4,4,6))
pca<-prcomp(t(x),scale=TRUE,center=TRUE)
plot(pca$x[,c(1,2)],bg=colors,pch=21,bty='n')
pca<-prcomp(t(y$pseudo.counts),scale=TRUE,center=TRUE)
plot(pca$x[,c(1,2)],bg=colors,pch=21,bty='n')
par(xpd=NA)
legend("topright",legend=c("control","EE2","Mix","Cd","BNF"),pch=21,pt.bg=unique(colors),inset=c(-0.3,-0.2))
dev.off()

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
read.results<-function(name){
  output<-read.table(paste('results_',name,'vs.control.txt',sep=''),sep='\t',quote='',comment.char='',header=TRUE)
  pseudocounts<-read.table(paste('pseudocounts_',name,'vs.control.txt',sep=''),sep='\t',quote='',comment.char='',header=TRUE,row.names=1)
  tmp<-list(results=output,pseudo.counts=pseudocounts)
  return(tmp)
}
# x.sel<-x[,exp %in% c("control", "EE2")]
# exp.sel<-exp[exp %in% c("control", "EE2")]
# exp.sel<-relevel(exp.sel,"control")
# out.ee2<-analyze.edgeR(x.sel,exp.sel,"EE2",annotation)
out.ee2<-read.results("EE2")

# x.sel<-x[,exp %in% c("control", "Mix")]
# exp.sel<-exp[exp %in% c("control", "Mix")]
# exp.sel<-relevel(exp.sel,"control")
# out.mix<-analyze.edgeR(x.sel,exp.sel,"Mix",annotation)
out.mix<-read.results("Mix")

# x.sel<-x[,exp %in% c("control", "Cd")]
# exp.sel<-exp[exp %in% c("control", "Cd")]
# exp.sel<-relevel(exp.sel,"control")
# out.cd<-analyze.edgeR(x.sel,exp.sel,"Cd",annotation)
out.cd<-read.results("Cd")

# x.sel<-x[,exp %in% c("control", "BNF")]
# exp.sel<-exp[exp %in% c("control", "BNF")]
# exp.sel<-relevel(exp.sel,"control")
# out.bnf<-analyze.edgeR(x.sel,exp.sel,"BNF",annotation)
out.bnf<-read.results("BNF")

pdf('pca_pairwise.pdf',height=10,width=10)
par(mfrow=c(2,2),mar=c(4,4,4,6))
colors<-c(rep("darkred",5),rep("dodgerblue",5))
par(xpd=NA)
namevec=c("EE2","Mix","Cd","BNF")
for (i in 1:length(namevec)){
  pc<-read.table(paste('pseudocounts_',namevec[[i]],'vs.control.txt',sep=''),header=TRUE,row.names=1)
  pca<-prcomp(t(sqrt(pc)))
  explained.pc1<-summary(pca)$importance[2,1]
  explained.pc2<-summary(pca)$importance[2,2]
  plot(pca$x[,c(1,2)],bg=colors,pch=21,bty='n',main=paste('PCA ', namevec[[i]], " vs. control",sep=''),xlab='',ylab='')
  title(xlab=sprintf("PC1(%1.0f%%)",explained.pc1*100),ylab=sprintf("PC2 (%1.0f%%)",explained.pc2*100))
  legend("topright",legend=c("control",namevec[[i]]),pch=21,pt.bg=unique(colors),inset=c(-0.2,-0.05))
}


dev.off()

pdf('pca_pairwise_nn.pdf',height=10,width=10)
par(mfrow=c(2,2),mar=c(4,4,4,6))
colors<-c(rep("darkred",5),rep("dodgerblue",5))
par(xpd=NA)
namevec=c("EE2","Mix","Cd","BNF")
for (i in 1:length(namevec)){
  x.sel<-x[,exp %in% c("control", namevec[[i]])]
  pca<-prcomp(t(sqrt(x.sel)))
  explained.pc1<-summary(pca)$importance[2,1]
  explained.pc2<-summary(pca)$importance[2,2]
  plot(pca$x[,c(1,2)],bg=colors,pch=21,bty='n',main=paste('PCA ', namevec[[i]], " vs. control",sep=''),xlab='',ylab='')
  title(xlab=sprintf("PC1(%1.0f%%)",explained.pc1*100),ylab=sprintf("PC2 (%1.0f%%)",explained.pc2*100))
  legend("topright",legend=c("control",namevec[[i]]),pch=21,pt.bg=unique(colors),inset=c(-0.2,-0.05))
}


dev.off()
cutoff<-0.05
sig.ee2<-out.ee2$results[out.ee2$results$FDR<cutoff,"Row.names"]
sig.mix<-out.mix$results[out.mix$results$FDR<cutoff,"Row.names"]
sig.cd<-out.cd$results[out.cd$results$FDR<cutoff,"Row.names"]
sig.bnf<-out.bnf$results[out.bnf$results$FDR<cutoff,"Row.names"]
library(gplots)
pdf('venn_0.05.pdf',height=5,width=5)
par(xpd=NA)
v1<-venn(list(EE2=sig.ee2,Mix=sig.mix,Cd=sig.cd,BNF=sig.bnf))
mtext('Number of significant transcripts, FDR<0.05')
dev.off()
#inputvenn<-list(out.ee2)
sig.ee2.up<-out.ee2$results[out.ee2$results$FDR<cutoff & out.ee2$results$logFC>0,"Row.names"]
sig.mix.up<-out.mix$results[out.mix$results$FDR<cutoff & out.mix$results$logFC>0,"Row.names"]
sig.cd.up<-out.cd$results[out.cd$results$FDR<cutoff & out.cd$results$logFC>0,"Row.names"]
sig.bnf.up<-out.bnf$results[out.bnf$results$FDR<cutoff & out.bnf$results$logFC>0,"Row.names"]
pdf('venn_up_0.05.pdf',height=5,width=5)
par(xpd=NA)
v2<-venn(list(EE2=sig.ee2.up,Mix=sig.mix.up,Cd=sig.cd.up,BNF=sig.bnf.up))
mtext("Number of up-regulated transcripts, FDR<0.05")
dev.off()
sig.ee2.down<-out.ee2$results[out.ee2$results$FDR<cutoff & out.ee2$results$logFC<0,"Row.names"]
sig.mix.down<-out.mix$results[out.mix$results$FDR<cutoff & out.mix$results$logFC<0,"Row.names"]
sig.cd.down<-out.cd$results[out.cd$results$FDR<cutoff & out.cd$results$logFC<0,"Row.names"]
sig.bnf.down<-out.bnf$results[out.bnf$results$FDR<cutoff & out.bnf$results$logFC<0,"Row.names"]
pdf('venn_down_0.05.pdf',height=5,width=5)
par(xpd=NA)
v3<-venn(list(EE2=sig.ee2.down,Mix=sig.mix.down,Cd=sig.cd.down,BNF=sig.bnf.down))
mtext("Number of down-regulated transcripts, FDR<0.05")
dev.off()
statstable<-matrix(nrow=4,ncol=3)
statstable[,1]<-c(length(sig.ee2),length(sig.mix),length(sig.cd),length(sig.bnf))
statstable[,2]<-c(length(sig.ee2.up),length(sig.mix.up),length(sig.cd.up),length(sig.bnf.up))
statstable[,3]<-c(length(sig.ee2.down),length(sig.mix.down),length(sig.cd.down),length(sig.bnf.down))
rownames(statstable)<-c("EE2 vs. control","Mix vs. control", "Cd vs. control","BNF vs. control")
colnames(statstable)<-c("Number of significant (FDR<0.05)","Up-regulated (FDR<0.05)","Down-regulated(FDR<0.05)")
statstable<-t(statstable)
statstable<-statstable[,c(1,3,4,2)]
write.table(statstable,"NumberOfSig_westcoast.txt",sep='\t',quote=F,col.names=NA)

#1100
sig1100<-attr(v1,"intersections")$'1100'
sig1100.up<-attr(v2,"intersections")$'1100'
sig1100.down<-attr(v3,"intersections")$'1100'
sig1100.res<-out.ee2$results[out.ee2$results$Row.names %in% sig1100,]
colnames(sig1100.res)[colnames(sig1100.res)=='logFC']<-'logFC EE2'
colnames(sig1100.res)[colnames(sig1100.res)=='FDR']<-'FDR EE2'
tmp<-out.mix$results[,c("Row.names","logFC","FDR")]
colnames(tmp)[colnames(tmp)=='logFC']<-'logFC Mix'
colnames(tmp)[colnames(tmp)=='FDR']<-'FDR Mix'
sig1100.res<-merge(sig1100.res,tmp,by.x=1,by.y=1,all.x=TRUE)
sig1100.res<-sig1100.res[order(sig1100.res$"FDR EE2"),]
write.table(sig1100.res,'results1100.txt',sep='\t',quote=FALSE,row.names=FALSE)
#0110
sig0110<-attr(v1,"intersections")$'0110'
sig0110.up<-attr(v2,"intersections")$'0110'
sig0110.down<-attr(v3,"intersections")$'0110'
sig0110.res<-out.cd$results[out.ee2$results$Row.names %in% sig0110,]
colnames(sig0110.res)[colnames(sig0110.res)=='logFC']<-'logFC Cd'
colnames(sig0110.res)[colnames(sig0110.res)=='FDR']<-'FDR Cd'
tmp<-out.mix$results[,c("Row.names","logFC","FDR")]
colnames(tmp)[colnames(tmp)=='logFC']<-'logFC Mix'
colnames(tmp)[colnames(tmp)=='FDR']<-'FDR Mix'
sig0110.res<-merge(sig0110.res,tmp,by.x=1,by.y=1,all.x=TRUE)
sig0110.res<-sig0110.res[order(sig0110.res$"FDR Cd"),]
write.table(sig0110.res,'results0110.txt',sep='\t',quote=FALSE,row.names=FALSE)
#0101
sig0101<-attr(v1,"intersections")$'0101'
sig0101.up<-attr(v2,"intersections")$'0101'
sig0101.down<-attr(v3,"intersections")$'0101'
sig0101.res<-out.bnf$results[out.ee2$results$Row.names %in% sig0101,]
colnames(sig0101.res)[colnames(sig0101.res)=='logFC']<-'logFC BNF'
colnames(sig0101.res)[colnames(sig0101.res)=='FDR']<-'FDR BNF'
tmp<-out.mix$results[,c("Row.names","logFC","FDR")]
colnames(tmp)[colnames(tmp)=='logFC']<-'logFC Mix'
colnames(tmp)[colnames(tmp)=='FDR']<-'FDR Mix'
sig0101.res<-merge(sig0101.res,tmp,by.x=1,by.y=1,all.x=TRUE)
sig0101.res<-sig0101.res[order(sig0101.res$"FDR BNF"),]
write.table(sig0101.res,'results0101.txt',sep='\t',quote=FALSE,row.names=FALSE)

