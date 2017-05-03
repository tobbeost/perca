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