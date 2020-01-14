setwd("D:\\GEO\\normal")
library(affyPLM)
library(affy)
Data<-ReadAffy()
sampleNames(Data)
N=length(Data)
eset.rma<-rma(Data)
normal_exprs<-exprs(eset.rma)
probeid<-rownames(normal_exprs)
normal_exprs<-cbind(probeid,normal_exprs)
write.table(normal_exprs,file="normal.expres.txt",sep='\t',quote=F,row.names=F)

setwd("D:\\GEO\\tumor")
library(affyPLM)
library(affy)
Data<-ReadAffy()
sampleNames(Data)
N=length(Data)
eset.rma<-rma(Data)
normal_exprs<-exprs(eset.rma)
probeid<-rownames(normal_exprs)
normal_exprs<-cbind(probeid,normal_exprs)
write.table(normal_exprs,file="tumor.expres.txt",sep='\t',quote=F,row.names=F)


setwd("D:\\LC")
normal_exprs<-read.table("normal.expres.txt",header=T,sep="\t")
tumor_exprs<-read.table("tumor.expres.txt",header=T,sep="\t")
probe_exprs<-merge(normal_exprs,tumor_exprs,by="probeid")
write.table(probe_exprs,file="cancer.probeid.exprs.txt",sep='\t',quote=F,row.names=F)

probe_exp<-read.table("cancer.probeid.exprs.txt",header=T,sep="\t",row.names=1)
probeid_geneid<-read.table("GPL.txt",header=T,sep="\t")
probe_name<-rownames(probe_exp)
loc<-match(probeid_geneid[,1],probe_name)
probe_exp<-probe_exp[loc,]
raw_geneid<-as.numeric(as.matrix(probeid_geneid[,3]))
index<-which(!is.na(raw_geneid))
geneid<-raw_geneid[index]
exp_matrix<-probe_exp[index,]
geneidfactor<-factor(geneid)
gene_exp_matrix<-apply(exp_matrix,2,function(x) tapply(x,geneidfactor,mean))
rownames(gene_exp_matrix)<-levels(geneidfactor)
geneid<-rownames(gene_exp_matrix)
gene_exp_matrix2<-cbind(geneid,gene_exp_matrix)
write.table(gene_exp_matrix2,file="Gastric.cancer.geneid.exprs.txt",sep='\t',quote=F,row.names=F)

loc<-match(rownames(gene_exp_matrix),probeid_geneid[,3])
rownames(gene_exp_matrix)=probeid_geneid[loc,2]
genesymbol<-rownames(gene_exp_matrix)
gene_exp_matrix3<-cbind(genesymbol,gene_exp_matrix)
write.table(gene_exp_matrix3,file="Gastric.cancer.genesyb.exprs.txt",sep='\t',quote=F,row.names=F)

library(impute)
gene_exp_matrix<-read.table("Gastric.cancer.genesyb.exprs.txt",header=T,sep="\t",row.names=1)
gene_exp_matrix<-as.matrix(gene_exp_matrix)
imputed_gene_exp<-impute.knn(gene_exp_matrix,k=10,rowmax=0.5,colmax=0.8,maxp=3000,rng.seed=362436069)
GeneExp<-imputed_gene_exp$data
genesymbol<-rownames(GeneExp)
GeneExp<-cbind(genesymbol,GeneExp)
write.table(GeneExp,file="Gastric.cancer.gene.exprs.txt",sep='\t',quote=F,row.names=F)

library(limma)
rt<-read.table("Gastric.cancer.gene.exprs.txt",header=T,sep="\t",row.names="genesymbol")
#differential
class<-c(rep("normal",NUM_normal),rep("tumor",NUM_tumor))
design<-model.matrix(~factor(class))
colnames(design)<-c("normal","tumor")
fit<-lmFit(rt,design)
fit2<-eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=200000)
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)

#write table
diffLab<-allDiff[with(allDiff, ((logFC>1 |logFC<(-1)) & adj.P.Val<0.05)),]
write.table(diffLab,file="diffEXp.xls",sep="\t",quote=F)
diffExpLevel<-rt[rownames(diffLab),]
qvalue=allDiff[rownames(diffLab),]$adj.P.Val
diffExpQvalue=cbind(qvalue,diffExpLevel)
write.table(diffExpQvalue,file="diffExpLevel.xls",sep="\t",quote=F)

hmExp=log10(diffExpLevel+0.00001)
library('gplots')
hmMat=as.matrix(hmExp)
pdf(file="heatmap.pdf",height=120,width=90)
par(oma=c(3,3,3,5))
heatmap.2(hmMat,col='greenred',trace="none",cexCol=1)
dev.off()


##
setwd("D:\\TCGA\\DESeq")
a<-read.table("GPL2.txt",header=T,sep="\t")
b<-read.table("diffSig.xls",header=T,sep="\t")
c<-merge(a,b,by="genesymbol")
x<-unique(c)
write.table(x,"x.txt",sep="\t")

source("http://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")
biocLite("pathview")
library("clusterProfiler")
rt=read.table("X.txt",sep="\t",header=T,check.names=F)
geneFC=rt$logFC
gene<-rt$ID
names(geneFC)=gene

#kegg
kk<-enrichKEGG(gene=gene,organism="human",qvalueCutoff=0.05)
write.table(summary(kk),file="KEGG.xls",sep="\t",quote=F,row.names=F)
pdf(file="KEGG.barplot.pdf")
barplot(kk,drop=TRUE,showCategory=12)
pdf(file="KEGG.cnetplot.pdf")
cnetplot(kk,categorySize="geneNum",foldChange=geneFC)   
library("pathview")
keggxls=read.table("KEGG.xls",sep="\t",header=T)
for(i in keggxls$ID){
pv.out<-pathview(gene.data=geneFC,pathway.id=i,species="hsa",out.suffix="pathview")}
