setwd("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348")
gtf_info_genes=read.csv("/nfs/users/mdierssen/idetoma/genomes/gtf_info_genes.csv", header=T, row.names=1)


countData=data.frame(
  control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182245/SRR1182245ReadsPerGene.out.tab", header=F, skip=4)[,2],
  control2=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182247/SRR1182247ReadsPerGene.out.tab", header=F, skip=4)[,2],
  control3=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182250/SRR1182250ReadsPerGene.out.tab", header=F, skip=4)[,2]+read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182251/SRR1182251ReadsPerGene.out.tab", header=F, skip=4)[,2],
  control4=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182255/SRR1182255ReadsPerGene.out.tab", header=F, skip=4)[,2]+read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182254/SRR1182254ReadsPerGene.out.tab", header=F, skip=4)[,2],
  disease1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182244/SRR1182244ReadsPerGene.out.tab", header=F, skip=4)[,2],
  disease2=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182246/SRR1182246ReadsPerGene.out.tab", header=F, skip=4)[,2],
  disease3=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182248/SRR1182248ReadsPerGene.out.tab", header=F, skip=4)[,2]+read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182249/SRR1182249ReadsPerGene.out.tab", header=F, skip=4)[,2],
  disease4=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182252/SRR1182252ReadsPerGene.out.tab", header=F, skip=4)[,2]+read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182253/SRR1182253ReadsPerGene.out.tab", header=F, skip=4)[,2]
  )
  

control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182245/SRR1182245ReadsPerGene.out.tab", header=F, skip=4)
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



colData=cbind(c(rep("control",4), rep("disease",4)))
rownames(colData)=colnames(countData)
colnames(colData)=c("treatment")

treatment=as.factor(c(rep("control", 4),rep("disease", 4)))

design=model.matrix(~0+treatment)
colnames(design) <- levels(treatment)


require(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = as.data.frame(colData),
                              design = ~ treatment )

dds <- DESeq(dds)
res <- results(dds)


resLFC <- lfcShrink(dds, coef=2)

save.image("SRP039348_twins.RData")


countData=data.frame(
  control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182262/SRR1182262ReadsPerGene.out.tab", header=F, skip=4)[,2],
  control2=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182263/SRR1182263ReadsPerGene.out.tab", header=F, skip=4)[,2],
  control3=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182264/SRR1182264ReadsPerGene.out.tab", header=F, skip=4)[,2],
  control4=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182265/SRR1182265ReadsPerGene.out.tab", header=F, skip=4)[,2],
  control5=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182266/SRR1182266ReadsPerGene.out.tab", header=F, skip=4)[,2],
  control6=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182267/SRR1182267ReadsPerGene.out.tab", header=F, skip=4)[,2],
  control7=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182268/SRR1182268ReadsPerGene.out.tab", header=F, skip=4)[,2],
  control8=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182269/SRR1182269ReadsPerGene.out.tab", header=F, skip=4)[,2],
  disease1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182270/SRR1182270ReadsPerGene.out.tab", header=F, skip=4)[,2],
  disease2=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182271/SRR1182271ReadsPerGene.out.tab", header=F, skip=4)[,2],    
  disease3=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182272/SRR1182272ReadsPerGene.out.tab", header=F, skip=4)[,2],
  disease4=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182273/SRR1182273ReadsPerGene.out.tab", header=F, skip=4)[,2],
  disease5=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182274/SRR1182274ReadsPerGene.out.tab", header=F, skip=4)[,2],
  disease6=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182275/SRR1182275ReadsPerGene.out.tab", header=F, skip=4)[,2],
  disease7=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182276/SRR1182276ReadsPerGene.out.tab", header=F, skip=4)[,2],
  disease8=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182277/SRR1182277ReadsPerGene.out.tab", header=F, skip=4)[,2]
  )
  

control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182262/SRR1182262ReadsPerGene.out.tab", header=F, skip=4)
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



colData=cbind(c(rep("control",8), rep("disease",8)))
rownames(colData)=colnames(countData)
colnames(colData)=c("treatment")

treatment=as.factor(c(rep("control", 8),rep("disease", 8)))

design=model.matrix(~0+treatment)
colnames(design) <- levels(treatment)


require(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = as.data.frame(colData),
                              design = ~ treatment )

dds <- DESeq(dds)
res <- results(dds)


resLFC <- lfcShrink(dds, coef=2)

save.image("SRP039348_fibroblasts.RData")


gtf_info_genes_m=read.csv("/nfs/users/mdierssen/idetoma/genomes/gtf_info_genes_m.csv", header=T, row.names=1)

countData=data.frame(
  control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182259/SRR1182259ReadsPerGene.out.tab", header=F, skip=4)[,2],

  disease1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182258/SRR1182258ReadsPerGene.out.tab", header=F, skip=4)[,2]
   )
  

control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182258/SRR1182258ReadsPerGene.out.tab", header=F, skip=4)
rownames(control1)=control1$V1

# gtf_neg=gtf_info_genes_m[which(gtf_info_genes_m$V7=="-"),]
# gtf_pos=gtf_info_genes_m[which(gtf_info_genes_m$V7=="+"),]
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



colData=cbind(c(rep("control",1), rep("disease",1)))
rownames(colData)=colnames(countData)
colnames(colData)=c("treatment")

treatment=as.factor(c(rep("control", 1),rep("disease", 1)))

design=model.matrix(~0+treatment)
colnames(design) <- levels(treatment)


require(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = as.data.frame(colData),
                              design = ~ treatment )

dds <- DESeq(dds)
res <- results(dds)


resLFC <- lfcShrink(dds, coef=2)

save.image("SRP039348_ipsc.RData")





gtf_info_genes_m=read.csv("/nfs/users/mdierssen/idetoma/genomes/gtf_info_genes_m.csv", header=T, row.names=1)

countData=data.frame(
  control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182261/SRR1182261ReadsPerGene.out.tab", header=F, skip=4)[,2],

  disease1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182260/SRR1182260ReadsPerGene.out.tab", header=F, skip=4)[,2]
   )
  

control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP039348/SRR1182260/SRR1182260ReadsPerGene.out.tab", header=F, skip=4)
rownames(control1)=control1$V1

# gtf_neg=gtf_info_genes_m[which(gtf_info_genes_m$V7=="-"),]
# gtf_pos=gtf_info_genes_m[which(gtf_info_genes_m$V7=="+"),]
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



colData=cbind(c(rep("control",1), rep("disease",1)))
rownames(colData)=colnames(countData)
colnames(colData)=c("treatment")

treatment=as.factor(c(rep("control", 1),rep("disease", 1)))

design=model.matrix(~0+treatment)
colnames(design) <- levels(treatment)


require(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = as.data.frame(colData),
                              design = ~ treatment )

dds <- DESeq(dds)
res <- results(dds)


resLFC <- lfcShrink(dds, coef=2)

save.image("SRP039348_mef.RData")


