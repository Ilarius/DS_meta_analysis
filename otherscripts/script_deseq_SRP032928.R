setwd("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP032928")

gtf_info_genes=read.csv("/nfs/users/mdierssen/idetoma/genomes/gtf_info_genes.csv", header=T, row.names=1)

countData=data.frame(
control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP032928/SRR1028343/SRR1028343ReadsPerGene.out.tab", header=F, skip=4)[,2],
control2=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP032928/SRR1028344/SRR1028344ReadsPerGene.out.tab", header=F, skip=4)[,2],
control3=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP032928/SRR1028345/SRR1028345ReadsPerGene.out.tab", header=F, skip=4)[,2],
control4=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP032928/SRR1028346/SRR1028346ReadsPerGene.out.tab", header=F, skip=4)[,2],
disease1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP032928/SRR1028347/SRR1028347ReadsPerGene.out.tab", header=F, skip=4)[,2],
disease2=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP032928/SRR1028348/SRR1028348ReadsPerGene.out.tab", header=F, skip=4)[,2],
disease3=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP032928/SRR1028349/SRR1028349ReadsPerGene.out.tab", header=F, skip=4)[,2]
)


somme=rowSums(countData)

countData$sum=somme


control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP032928/SRR1028343/SRR1028343ReadsPerGene.out.tab", header=F, skip=4)

rownames(control1)=control1$V1


delete_dot=function(x){strsplit(as.character(x), split="[.]")[[1]][1]}

countData$names=unlist(lapply(control1$V1, delete_dot))



countData=countData[base::order(as.numeric(countData$sum), decreasing=T),]

countData=countData[!duplicated(countData$names),]

rownames(countData)=countData$names

countData=countData[,1:7]
# gtf_neg=gtf_info_genes[which(gtf_info_genes$V7=="-"),]
# gtf_pos=gtf_info_genes[which(gtf_info_genes$V7=="+"),]
# colSums(control1[gtf_neg$id,2:4])
#data looks unstranded
colData=cbind(c(rep("control",4), rep("disease",3)))
rownames(colData)=colnames(countData)
colnames(colData)=c("treatment")

treatment=as.factor(c(rep("control", 4),rep("disease", 3)))

design=model.matrix(~0+treatment)
colnames(design) <- levels(treatment)

require(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = as.data.frame(colData),
                              design = ~ treatment )

dds <- DESeq(dds)
res <- results(dds)


resLFC <- lfcShrink(dds, coef=2)

save.image("SRP032928.RData")
