setwd("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP164925")

gtf_info_genes=read.csv("/nfs/users/mdierssen/idetoma/genomes/gtf_info_genes.csv", header=T, row.names=1)

countData=data.frame(
control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP164925/SRR7993818/SRR7993818ReadsPerGene.out.tab", header=F, skip=4)[,3],
control2=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP164925/SRR7993819/SRR7993819ReadsPerGene.out.tab", header=F, skip=4)[,3],
control3=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP164925/SRR7993820/SRR7993820ReadsPerGene.out.tab", header=F, skip=4)[,3],
disease1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP164925/SRR7993815/SRR7993815ReadsPerGene.out.tab", header=F, skip=4)[,3],
disease2=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP164925/SRR7993816/SRR7993816ReadsPerGene.out.tab", header=F, skip=4)[,3],
disease3=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP164925/SRR7993817/SRR7993817ReadsPerGene.out.tab", header=F, skip=4)[,3]
)


somme=rowSums(countData)

countData$sum=somme


control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP164925/SRR7993818/SRR7993818ReadsPerGene.out.tab", header=F, skip=4)

rownames(control1)=control1$V1


delete_dot=function(x){strsplit(as.character(x), split="[.]")[[1]][1]}

countData$names=unlist(lapply(control1$V1, delete_dot))



countData=countData[base::order(as.numeric(countData$sum), decreasing=T),]

countData=countData[!duplicated(countData$names),]

rownames(countData)=countData$names

countData=countData[,1:6]

# colSums(control1[,2:4], na.rm=T)
#If the counts from yes and reverse are roughly equal for most genes then the dataset is unstranded 
# If either yes or reverse produces much higher counts than the other then the appropriate setting is the one giving the higher counts.
#V2 unstranded
#V3 yes
#V4 reverse

colData=cbind(c(rep("control",3), rep("disease",3)))
rownames(colData)=colnames(countData)
colnames(colData)=c("treatment")

treatment=as.factor(c(rep("control", 3),rep("disease", 3)))

design=model.matrix(~0+treatment)
colnames(design) <- levels(treatment)

require(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = as.data.frame(colData),
                              design = ~ treatment )

dds <- DESeq(dds)
res <- results(dds)


resLFC <- lfcShrink(dds, coef=2)


save.image("SRP164925.RData")
