setwd("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078911/")

gtf_info_genes=read.csv("/nfs/users/mdierssen/idetoma/genomes/gtf_info_genes_m.csv", header=T, row.names=1)

countData=data.frame(
  control1=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078911/SRR3929997/SRR3929997ReadsPerGene.out.tab", header=F, skip=4)[,2],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078911/SRR3929998/SRR3929998ReadsPerGene.out.tab", header=F, skip=4)[,2]),
  
  control2=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078911/SRR3929999/SRR3929999ReadsPerGene.out.tab", header=F, skip=4)[,2],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078911/SRR3930000/SRR3930000ReadsPerGene.out.tab", header=F, skip=4)[,2]),
  
  control3=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078911/SRR3930001/SRR3930001ReadsPerGene.out.tab", header=F, skip=4)[,2],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078911/SRR3930002/SRR3930002ReadsPerGene.out.tab", header=F, skip=4)[,2]),
  
  disease1=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078911/SRR3930003/SRR3930003ReadsPerGene.out.tab", header=F, skip=4)[,2],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078911/SRR3930004/SRR3930004ReadsPerGene.out.tab", header=F, skip=4)[,2]),
  disease2=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078911/SRR3930005/SRR3930005ReadsPerGene.out.tab", header=F, skip=4)[,2],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078911/SRR3930006/SRR3930006ReadsPerGene.out.tab", header=F, skip=4)[,2]),
  disease3=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078911/SRR3930007/SRR3930007ReadsPerGene.out.tab", header=F, skip=4)[,2],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078911/SRR3930008/SRR3930008ReadsPerGene.out.tab", header=F, skip=4)[,2])
)


somme=rowSums(countData)

countData$sum=somme


control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078911/SRR3929997/SRR3929997ReadsPerGene.out.tab", header=F, skip=4)

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

save.image("SRP078911.RData")
