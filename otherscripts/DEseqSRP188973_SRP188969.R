setwd("/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188973")
gtf_info_genes=read.csv("/nfs/users/mdierssen/idetoma/genomes/gtf_info_genes.csv", header=T, row.names=1)

metadata_SRP188973=read.table("/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188973/SRP188973.csv", sep=",",header=T)
metadata_SRP188969=read.table("/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188969/SRP188969.csv", sep=",",header=T)

countData=data.frame(
  SRR8757123=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188973/SRR8757123/SRR8757123ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR8757124=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188973/SRR8757124/SRR8757124ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR8757125=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188973/SRR8757125/SRR8757125ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR8757127=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188973/SRR8757127/SRR8757127ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR8757128=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188973/SRR8757128/SRR8757128ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR8757129=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188973/SRR8757129/SRR8757129ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR8757130=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188973/SRR8757130/SRR8757130ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR8757131=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188973/SRR8757131/SRR8757131ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR8757132=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188973/SRR8757132/SRR8757132ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR8757133=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188973/SRR8757133/SRR8757133ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR8757134=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188973/SRR8757134/SRR8757134ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR8757135=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188973/SRR8757135/SRR8757135ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR8757136=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188973/SRR8757136/SRR8757136ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR8757126=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188973/SRR8757126/SRR8757126ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR8757082=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188969/SRR8757082/SRR8757082ReadsPerGene.out.tab", header=F, skip=4)[,2],
  SRR8757074=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188969/SRR8757074/SRR8757074ReadsPerGene.out.tab", header=F, skip=4)[,2],
  SRR8757076=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188969/SRR8757076/SRR8757076ReadsPerGene.out.tab", header=F, skip=4)[,2],
  SRR8757078=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188969/SRR8757078/SRR8757078ReadsPerGene.out.tab", header=F, skip=4)[,2],
  SRR8757080=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188969/SRR8757080/SRR8757080ReadsPerGene.out.tab", header=F, skip=4)[,2]
  
)


control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188973/SRR8757123/SRR8757123ReadsPerGene.out.tab", header=F, skip=4)
rownames(control1)=control1$V1

# gtf_neg=gtf_info_genes[which(gtf_info_genes$V7=="-"),]
# gtf_pos=gtf_info_genes[which(gtf_info_genes$V7=="+"),]
# colSums(control1[gtf_neg$id,2:4])
# #data stranded if V3 or V4 have the highest number -strand=reverse -> 4th column


somme=rowSums(countData)
# 
countData$sum=somme


delete_dot=function(x){strsplit(as.character(x), split="[.]")[[1]][1]}

countData$names=unlist(lapply(control1$V1, delete_dot))



countData=countData[base::order(as.numeric(countData$sum), decreasing=T),]

countData=countData[!duplicated(countData$names),]

rownames(countData)=countData$names



countData=countData[,1:(dim(countData)[2]-2)]

metadata_SRP188969=metadata_SRP188969[, intersect(colnames(metadata_SRP188969),colnames(metadata_SRP188973))]
metadata_SRP188973=metadata_SRP188973[, intersect(colnames(metadata_SRP188969),colnames(metadata_SRP188973))]
colData=rbind(metadata_SRP188973, metadata_SRP188969)



colData=as.data.frame(colData)

rownames(colData)=colnames(countData)

colData$batch=c(rep("first", length(rownames(metadata_SRP188973))), rep("second", length(rownames(metadata_SRP188969))))
require(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = as.data.frame(colData),
                              design = ~ sex + Age +batch+ karyotype)

dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 3) >= 3
dds <- dds[filter,]
dds <- DESeq(dds, full=design(dds), reduced=~ sex + Age +batch, test="LRT")
res <- results(dds)


resLFC <- lfcShrink(dds, coef=2)

save.image("SRP188973_SRP188969.RData")
