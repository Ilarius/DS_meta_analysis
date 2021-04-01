setwd("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769")
gtf_info_genes=read.csv("/nfs/users/mdierssen/idetoma/genomes/gtf_info_genes.csv", header=T, row.names=1)

metadata_SRP072769=read.table("/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRP072769_SraRunTable.txt", sep="\t",header=T)

countData=data.frame(
  control1=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322275/SRR3322275ReadsPerGene.out.tab", header=F, skip=4)[,4],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322276/SRR3322276ReadsPerGene.out.tab", header=F, skip=4)[,4]),
  
  control2=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322277/SRR3322277ReadsPerGene.out.tab", header=F, skip=4)[,4],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322278/SRR3322278ReadsPerGene.out.tab", header=F, skip=4)[,4]),
  
  control3=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322279/SRR3322279ReadsPerGene.out.tab", header=F, skip=4)[,4],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322280/SRR3322280ReadsPerGene.out.tab", header=F, skip=4)[,4]),
  
  
  control4=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322281/SRR3322281ReadsPerGene.out.tab", header=F, skip=4)[,4],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322282/SRR3322282ReadsPerGene.out.tab", header=F, skip=4)[,4]),
  
  control5=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322283/SRR3322283ReadsPerGene.out.tab", header=F, skip=4)[,4],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322284/SRR3322284ReadsPerGene.out.tab", header=F, skip=4)[,4]),
  
  control6=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322285/SRR3322285ReadsPerGene.out.tab", header=F, skip=4)[,4],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322286/SRR3322286ReadsPerGene.out.tab", header=F, skip=4)[,4]),
  
  disease1=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322287/SRR3322287ReadsPerGene.out.tab", header=F, skip=4)[,4],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322288/SRR3322288ReadsPerGene.out.tab", header=F, skip=4)[,4]),
  
  
  disease2=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322289/SRR3322289ReadsPerGene.out.tab", header=F, skip=4)[,4],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322290/SRR3322290ReadsPerGene.out.tab", header=F, skip=4)[,4]),
  
  disease3=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322291/SRR3322291ReadsPerGene.out.tab", header=F, skip=4)[,4],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322292/SRR3322292ReadsPerGene.out.tab", header=F, skip=4)[,4]),
  
  disease4=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322293/SRR3322293ReadsPerGene.out.tab", header=F, skip=4)[,4],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322294/SRR3322294ReadsPerGene.out.tab", header=F, skip=4)[,4]),
  
  
  disease5=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322295/SRR3322295ReadsPerGene.out.tab", header=F, skip=4)[,4],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322296/SRR3322296ReadsPerGene.out.tab", header=F, skip=4)[,4]),
  
  disease6=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322297/SRR3322297ReadsPerGene.out.tab", header=F, skip=4)[,4],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322298/SRR3322298ReadsPerGene.out.tab", header=F, skip=4)[,4])
  
)


control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322275/SRR3322275ReadsPerGene.out.tab", header=F, skip=4)
rownames(control1)=control1$V1

# gtf_neg=gtf_info_genes[which(gtf_info_genes$V7=="-"),]
# gtf_pos=gtf_info_genes[which(gtf_info_genes$V7=="+"),]
# colSums(control1[gtf_neg$id,2:4])
# #data stranded if V3 or V4 have the highest number


somme=rowSums(countData)
# 
countData$sum=somme


delete_dot=function(x){strsplit(as.character(x), split="[.]")[[1]][1]}

countData$names=unlist(lapply(control1$V1, delete_dot))



countData=countData[base::order(as.numeric(countData$sum), decreasing=T),]

countData=countData[!duplicated(countData$names),]

rownames(countData)=countData$names



countData=countData[,1:(dim(countData)[2]-2)]



colData=cbind(c(rep("control",6), rep("disease",6)))
rownames(colData)=colnames(countData)
colnames(colData)=c("treatment")
colData=as.data.frame(colData)
colData$sex=metadata_SRP072769$gender[1:24][c(TRUE, FALSE)]
colData$batch=metadata_SRP072769$batch[1:24][c(TRUE, FALSE)]
colData$age=metadata_SRP072769$age[1:24][c(TRUE, FALSE)]

require(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = as.data.frame(colData),
                              design = ~ sex + batch + treatment)

dds <- DESeq(dds, full=design(dds), reduced=~ sex + batch, test="LRT")
res <- results(dds)


resLFC <- lfcShrink(dds, coef=2)

save.image("SRP072769_fibroblasts.RData")



countData=data.frame(
  control1=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322299/SRR3322299ReadsPerGene.out.tab", header=F, skip=4)[,4],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322300/SRR3322300ReadsPerGene.out.tab", header=F, skip=4)[,4]),
  
  control2=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322301/SRR3322301ReadsPerGene.out.tab", header=F, skip=4)[,4],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322302/SRR3322302ReadsPerGene.out.tab", header=F, skip=4)[,4]),
  
  control3=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322303/SRR3322303ReadsPerGene.out.tab", header=F, skip=4)[,4],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322304/SRR3322304ReadsPerGene.out.tab", header=F, skip=4)[,4]),
  
  
 
  
  disease1=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322305/SRR3322305ReadsPerGene.out.tab", header=F, skip=4)[,4],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322306/SRR3322306ReadsPerGene.out.tab", header=F, skip=4)[,4]),
  
  
  disease2=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322307/SRR3322307ReadsPerGene.out.tab", header=F, skip=4)[,4],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322308/SRR3322308ReadsPerGene.out.tab", header=F, skip=4)[,4]),
  
  disease3=c(read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322309/SRR3322309ReadsPerGene.out.tab", header=F, skip=4)[,4],
             read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP072769/SRR3322310/SRR3322310ReadsPerGene.out.tab", header=F, skip=4)[,4])
  
)


somme=rowSums(countData)

countData$sum=somme
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

save.image("SRP072769_lymphoblastoids.RData")
