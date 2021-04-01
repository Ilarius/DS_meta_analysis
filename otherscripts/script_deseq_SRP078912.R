setwd("/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912")
gtf_info_genes=read.csv("/nfs/users/mdierssen/idetoma/genomes/gtf_info_genes.csv", header=T, row.names=1)

metadata_SRP078912=read.table("/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRP078912_SraRunTable.txt", sep="\t",header=T)

countData=data.frame(
  SRR3930009=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930009/SRR3930009ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930010=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930010/SRR3930010ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930011=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930011/SRR3930011ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930012=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930012/SRR3930012ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930013=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930013/SRR3930013ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930014=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930014/SRR3930014ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930015=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930015/SRR3930015ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930016=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930016/SRR3930016ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930017=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930017/SRR3930017ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930018=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930018/SRR3930018ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930019=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930019/SRR3930019ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930020=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930020/SRR3930020ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930021=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930021/SRR3930021ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930022=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930022/SRR3930022ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930023=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930023/SRR3930023ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930024=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930024/SRR3930024ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930025=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930025/SRR3930025ReadsPerGene.out.tab", header=F, skip=4)[,4]

)


control1=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930009/SRR3930009ReadsPerGene.out.tab", header=F, skip=4)
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

sex=c("female", "female", "male", "male", "female", "female", "female", "female",
       rep("male",4), "female", "female", rep("male",4),rep("female",4), "male", "male",rep("female",10))

colData=metadata_SRP078912$condition[which(metadata_SRP078912$source_name=="Monocyte")]
colData=as.data.frame(colData)

rownames(colData)=colnames(countData)
colnames(colData)=c("treatment")
colData$sex=sex[which(metadata_SRP078912$source_name=="Monocyte")]
colData$batch=metadata_SRP078912$batch[which(metadata_SRP078912$source_name=="Monocyte")]

require(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = as.data.frame(colData),
                              design = ~ sex + batch + treatment)

dds <- DESeq(dds, full=design(dds), reduced=~ sex + batch, test="LRT")
res <- results(dds)


resLFC <- lfcShrink(dds, coef=2)

save.image("SRP078912_monocytes.RData")



countData=data.frame(
  SRR3930026=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930026/SRR3930026ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930027=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930027/SRR3930027ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930028=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930028/SRR3930028ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930029=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930029/SRR3930029ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930030=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930030/SRR3930030ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930031=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930031/SRR3930031ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930032=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930032/SRR3930032ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930033=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930033/SRR3930033ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930034=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930034/SRR3930034ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930035=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930035/SRR3930035ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930036=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930036/SRR3930036ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930037=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930037/SRR3930037ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930038=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930038/SRR3930038ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930039=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930039/SRR3930039ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930040=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930040/SRR3930040ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930041=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930041/SRR3930041ReadsPerGene.out.tab", header=F, skip=4)[,4],
  SRR3930042=read.table("/nfs/users/mdierssen/sequencing_data/DS_transcriptomics/SRP078912/SRR3930042/SRR3930042ReadsPerGene.out.tab", header=F, skip=4)[,4]
  
)


somme=rowSums(countData)

countData$sum=somme
countData$names=unlist(lapply(control1$V1, delete_dot))



countData=countData[base::order(as.numeric(countData$sum), decreasing=T),]

countData=countData[!duplicated(countData$names),]

rownames(countData)=countData$names

countData=countData[,1:(dim(countData)[2]-2)]
# gtf_neg=gtf_info_genes[which(gtf_info_genes$V7=="-"),]
# gtf_pos=gtf_info_genes[which(gtf_info_genes$V7=="+"),]
# colSums(control1[gtf_neg$id,2:4])
#data looks stranded
#V2       V3       V4 
#10623218   423297 10913226 
# htseq-count option -s reverse)
colData=metadata_SRP078912$condition[which(metadata_SRP078912$source_name=="T cells")]
colData=as.data.frame(colData)

rownames(colData)=colnames(countData)
colnames(colData)=c("treatment")

colData$sex=sex[which(metadata_SRP078912$source_name=="T cells")]

require(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = as.data.frame(colData),
                              design = ~ sex  + treatment)

dds <- DESeq(dds, full=design(dds), reduced=~ sex , test="LRT")
res <- results(dds)


resLFC <- lfcShrink(dds, coef=2)

save.image("SRP078912_tcells.RData")
