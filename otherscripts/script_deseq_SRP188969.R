setwd("/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188969")
gtf_info_genes=read.csv("/nfs/users/mdierssen/idetoma/genomes/gtf_info_genes.csv", header=T, row.names=1)

metadata_SRP188969=read.table("/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188969/SRP188969.csv", sep=",",header=T)

countData=data.frame(
  SRR8757082=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188969/SRR8757082/SRR8757082ReadsPerGene.out.tab", header=F, skip=4)[,2],
  SRR8757074=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188969/SRR8757074/SRR8757074ReadsPerGene.out.tab", header=F, skip=4)[,2],
  SRR8757076=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188969/SRR8757076/SRR8757076ReadsPerGene.out.tab", header=F, skip=4)[,2],
  SRR8757078=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188969/SRR8757078/SRR8757078ReadsPerGene.out.tab", header=F, skip=4)[,2],
  SRR8757080=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188969/SRR8757080/SRR8757080ReadsPerGene.out.tab", header=F, skip=4)[,2],
  
)


control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP188969/SRR8757123/SRR8757123ReadsPerGene.out.tab", header=F, skip=4)
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


colData=metadata_SRP188969
colData=as.data.frame(colData)

rownames(colData)=colnames(countData)


require(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = as.data.frame(colData),
                              design = ~   Age + karyotype)

dds <- DESeq(dds, full=design(dds), reduced=~ Age, test="LRT")
res <- results(dds)


resLFC <- lfcShrink(dds, coef=2)

save.image("SRP188969.RData")