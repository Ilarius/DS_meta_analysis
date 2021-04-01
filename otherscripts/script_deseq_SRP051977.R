setwd("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP051977/")

gtf_info_genes=read.csv("/nfs/users/mdierssen/idetoma/genomes/gtf_info_genes_m.csv", header=T, row.names=1)

countData=data.frame(
control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP051977/SRR1748088/SRR1748088ReadsPerGene.out.tab", header=F, skip=4)[,4],
control2=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP051977/SRR1748089/SRR1748089ReadsPerGene.out.tab", header=F, skip=4)[,4],
control3=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP051977/SRR1748090/SRR1748090ReadsPerGene.out.tab", header=F, skip=4)[,4],
disease1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP051977/SRR1748091/SRR1748091ReadsPerGene.out.tab", header=F, skip=4)[,4],
disease2=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP051977/SRR1748092/SRR1748092ReadsPerGene.out.tab", header=F, skip=4)[,4],
disease3=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP051977/SRR1748093/SRR1748093ReadsPerGene.out.tab", header=F, skip=4)[,4]
)


somme=rowSums(countData)

countData$sum=somme


control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP051977/SRR1748088/SRR1748088ReadsPerGene.out.tab", header=F, skip=4)

rownames(control1)=control1$V1


delete_dot=function(x){strsplit(as.character(x), split="[.]")[[1]][1]}

countData$names=unlist(lapply(control1$V1, delete_dot))



countData=countData[base::order(as.numeric(countData$sum), decreasing=T),]

countData=countData[!duplicated(countData$names),]

rownames(countData)=countData$names

countData=countData[,1:6]
# gtf_neg=gtf_info_genes[which(gtf_info_genes$V7=="-"),]
# gtf_pos=gtf_info_genes[which(gtf_info_genes$V7=="+"),]
# colSums(control1[gtf_neg$id,2:4])
#data looks stranded
#V2       V3       V4 
#10623218   423297 10913226 
# htseq-count option -s reverse)
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

save.image("SRP051977.RData")
