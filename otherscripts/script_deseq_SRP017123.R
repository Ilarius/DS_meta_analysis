setwd("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP017123")
gtf_info_genes=read.csv("/nfs/users/mdierssen/idetoma/genomes/gtf_info_genes.csv", header=T, row.names=1)


countData=data.frame(
  control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP017123/SRR611842/SRR611842ReadsPerGene.out.tab", header=F, skip=4)[,2],
  control2=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP017123/SRR611843/SRR611843ReadsPerGene.out.tab", header=F, skip=4)[,2],
  control3=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP017123/SRR611844/SRR611844ReadsPerGene.out.tab", header=F, skip=4)[,2],
  control4=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP017123/SRR611845/SRR611845ReadsPerGene.out.tab", header=F, skip=4)[,2],
  control5=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP017123/SRR611846/SRR611846ReadsPerGene.out.tab", header=F, skip=4)[,2],
  disease1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP017123/SRR611847/SRR611847ReadsPerGene.out.tab", header=F, skip=4)[,2],
  disease2=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP017123/SRR611848/SRR611848ReadsPerGene.out.tab", header=F, skip=4)[,2],
  disease3=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP017123/SRR611849/SRR611849ReadsPerGene.out.tab", header=F, skip=4)[,2],
  disease4=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP017123/SRR611850/SRR611850ReadsPerGene.out.tab", header=F, skip=4)[,2]
)
control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP017123/SRR611842/SRR611842ReadsPerGene.out.tab", header=F, skip=4)
rownames(control1)=control1$V1

# gtf_neg=gtf_info_genes[which(gtf_info_genes$V7=="-"),]
# gtf_pos=gtf_info_genes[which(gtf_info_genes$V7=="+"),]
# colSums(control1[gtf_neg$id,2:4])
# #data unstranded


somme=rowSums(countData)

countData$sum=somme


delete_dot=function(x){strsplit(as.character(x), split="[.]")[[1]][1]}

countData$names=unlist(lapply(control1$V1, delete_dot))



countData=countData[base::order(as.numeric(countData$sum), decreasing=T),]

countData=countData[!duplicated(countData$names),]

rownames(countData)=countData$names



countData=countData[,1:(dim(countData)[2]-2)]



colData=cbind(c(rep("control",5), rep("disease",4)))
rownames(colData)=colnames(countData)
colnames(colData)=c("treatment")
colData=as.data.frame(colData)
colData$sex=c("M", "M", "F", "F", "M", "F", "F", "F", "M")
colData$trimester=c("first", "first", "first", "first", "second", "second", "second", "first", "second")

require(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = as.data.frame(colData),
                              design = ~ sex + trimester + treatment)

dds <- DESeq(dds, full=design(dds), reduced=~ sex + trimester, test="LRT")
res <- results(dds)


resLFC <- lfcShrink(dds, coef=2)

save.image("SRP017123.RData")


