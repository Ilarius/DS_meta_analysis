setwd("/users/mdierssen/sequencing_data/DS_transcriptomics")
load("201911_de_list_ds.RData")


subset_metanalysis=function(list_array, adjpval=0.1, max_n_genes=20000, abslog2FC=0){
  #list_array is a list of list with names "list_ensembl", "list_adjpval", "list_log2FC"
  
  indexes1=lapply(list_array$adjpval, FUN=function(x) which(x < adjpval))
  nomi=names(list_array$names)
  list_array$names= lapply(1:length(list_array$names), FUN=function(x) (list_array$names[[x]][(indexes1[[x]])]))
  list_array$adjpval= lapply(1:length(list_array$adjpval), FUN=function(x) (list_array$adjpval[[x]][(indexes1[[x]])]))
  list_array$log2FC= lapply(1:length(list_array$log2FC), FUN=function(x) (list_array$log2FC[[x]][(indexes1[[x]])]))
  names(list_array$names)=nomi
  names(list_array$adjpval)=nomi
  names(list_array$log2FC)=nomi
  
  
  subset=which(lapply(list_array$names, length)>0)
  list_array$names=list_array$names[subset]
  list_array$adjpval=list_array$adjpval[subset]
  list_array$log2FC=list_array$log2FC[subset]
  nomi=names(list_array$names)
  
  indexes1=lapply(list_array$log2FC, FUN=function(x) which((x > abslog2FC )| (x < -abslog2FC )))
  
  list_array$names= lapply(1:length(list_array$names), FUN=function(x) (list_array$names[[x]][(indexes1[[x]])]))
  list_array$adjpval= lapply(1:length(list_array$adjpval), FUN=function(x) (list_array$adjpval[[x]][(indexes1[[x]])]))
  list_array$log2FC= lapply(1:length(list_array$log2FC), FUN=function(x) (list_array$log2FC[[x]][(indexes1[[x]])]))
  names(list_array$names)=nomi
  names(list_array$adjpval)=nomi
  names(list_array$log2FC)=nomi
  
  subset=which(lapply(list_array$names, length)>0)
  list_array$names=list_array$names[subset]
  list_array$adjpval=list_array$adjpval[subset]
  list_array$log2FC=list_array$log2FC[subset]
  
  
  for(x in 1:length(list_array$adjpval)){
    if(!is.na(list_array$adjpval[[x]][max_n_genes])){
      max_genes=max_n_genes
      a=TRUE
      while(a){
        if(list_array$adjpval[[x]][max_genes] == list_array$adjpval[[x]][max_genes+1]){
          max_genes=max_genes+1
          
        } else{a=FALSE}}
      
      list_array$names[[x]]= list_array$names[[x]][1:max_genes][!is.na(list_array$names[[x]][1:max_genes])]
      list_array$adjpval[[x]]=list_array$adjpval[[x]][1:max_genes][!is.na(list_array$adjpval[[x]][1:max_genes])]
      list_array$log2FC[[x]]=list_array$log2FC[[x]][1:max_genes][!is.na(list_array$log2FC[[x]][1:max_genes])]
      
    }
  }
  return(list_array)
}


list_array=list(list_ensembl, list_adjpval, list_log2FC)

names(list_array)=c("names", "adjpval", "log2FC")

## correct names
list_array$names$EMTAB312[grep(" /// ",list_array$names$EMTAB312)]=unlist(lapply(list_array$names$EMTAB312[grep(" /// ",list_array$names$EMTAB312)], FUN=function(x)(strsplit(x, " /// ")[[1]][1])))


list_array.1=subset_metanalysis(list_array)
list_array.05=subset_metanalysis(list_array, adjpval = 0.05)
list_array.01=subset_metanalysis(list_array, adjpval = 0.01)
list_array.001=subset_metanalysis(list_array, ƒadjpval = 0.001)


require(poweRlaw)
m4 = dislnorm$new(table(unlist(list_array.001$names)))
m4$setPars(estimate_pars(m4))
pvalues=rev(dist_cdf(m4, lower_tail = FALSE))

m4 = dislnorm$new(table(unlist(list_array.1$names)))
m4$setPars(estimate_pars(m4))
pvalues.1=rev(dist_cdf(m4, lower_tail = FALSE))


m4 = dislnorm$new(table(unlist(list_array.05$names)))
m4$setPars(estimate_pars(m4))
pvalues.05=rev(dist_cdf(m4, lower_tail = FALSE))

m4 = dislnorm$new(table(unlist(list_array.01$names)))
m4$setPars(estimate_pars(m4))
pvalues.01=rev(dist_cdf(m4, lower_tail = FALSE))

hubs.1=unique(names(table(unlist(list_array.1$names))[names(pvalues.1[which(pvalues.1<0.05)])]))
hubs.05=unique(names(table(unlist(list_array.05$names))[names(pvalues.05[which(pvalues.05<0.05)])]))
hubs.01=unique(names(table(unlist(list_array.01$names))[names(pvalues.1[which(pvalues.01<0.05)])]))
hubs.001=unique(names(table(unlist(list_array.001$names))[names(pvalues[which(pvalues<0.05)])]))

find_fold_changes=function(ensgene, dataset=list_array){
  dataset_true=dataset
  dataset_true$names=dataset$names[unlist(lapply(dataset$names, FUN=function(x) (ensgene %in% x)))]
  dataset_true$adjpval=dataset$adjpval[unlist(lapply(dataset$names, FUN=function(x) (ensgene %in% x)))]
  dataset_true$log2FC=dataset$log2FC[unlist(lapply(dataset$names, FUN=function(x) (ensgene %in% x)))]
  
  
  output=unlist(lapply(1:length(dataset_true$log2FC), FUN=function(x)(dataset_true$log2FC[[x]][which(dataset_true$names[[x]]==ensgene)])))
  names(output)=names(dataset_true$names)
  return(output)
}

list_array.1_fc1.5=subset_metanalysis(list_array, adjpval = 0.1, abslog2FC = log2(1.5))
list_array.05_fc1.5=subset_metanalysis(list_array, adjpval = 0.05, abslog2FC = log2(1.5))
list_array.01_fc1.5=subset_metanalysis(list_array, adjpval = 0.01, abslog2FC = log2(1.5))
list_array.001_fc1.5=subset_metanalysis(list_array, adjpval = 0.001, abslog2FC = log2(1.5))


list_array.1_fc1.5_max500=subset_metanalysis(list_array, adjpval = 0.1, abslog2FC = log2(1.5), max_n_genes = 500 )
list_array.05_fc1.5_max500=subset_metanalysis(list_array, adjpval = 0.05, abslog2FC = log2(1.5), max_n_genes = 500 )
list_array.01_fc1.5_max500=subset_metanalysis(list_array, adjpval = 0.01, abslog2FC = log2(1.5), max_n_genes = 500 )
list_array.001_fc1.5_max500=subset_metanalysis(list_array, adjpval = 0.001, abslog2FC = log2(1.5), max_n_genes = 500 )


m4 = dislnorm$new(table(unlist(list_array.1_fc1.5$names)))
m4$setPars(estimate_pars(m4))
pvalues.1_fc1.5=rev(dist_cdf(m4, lower_tail = FALSE))

m4 = dislnorm$new(table(unlist(list_array.001_fc1.5$names)))
m4$setPars(estimate_pars(m4))
pvalues.001_fc1.5=rev(dist_cdf(m4, lower_tail = FALSE))


m4 = dislnorm$new(table(unlist(list_array.05_fc1.5$names)))
m4$setPars(estimate_pars(m4))
pvalues.05_fc1.5=rev(dist_cdf(m4, lower_tail = FALSE))

m4 = dislnorm$new(table(unlist(list_array.01_fc1.5$names)))
m4$setPars(estimate_pars(m4))
pvalues.01_fc1.5=rev(dist_cdf(m4, lower_tail = FALSE))

hubs.1_fc1.5=unique(names(table(unlist(list_array.1_fc1.5$names))[names(pvalues.1_fc1.5[which(pvalues.1_fc1.5<0.05)])]))
hubs.05_fc1.5=unique(names(table(unlist(list_array.05_fc1.5$names))[names(pvalues.05_fc1.5[which(pvalues.05_fc1.5<0.05)])]))
hubs.01_fc1.5=unique(names(table(unlist(list_array.01_fc1.5$names))[names(pvalues.1_fc1.5[which(pvalues.01_fc1.5<0.05)])]))
hubs.001_fc1.5=unique(names(table(unlist(list_array.001_fc1.5$names))[names(pvalues.001_fc1.5[which(pvalues.001_fc1.5<0.05)])]))

m4 = dislnorm$new(table(unlist(list_array.1_fc1.5_max500$names)))
m4$setPars(estimate_pars(m4))
pvalues.1_fc1.5_max500=rev(dist_cdf(m4, lower_tail = FALSE))


m4 = dislnorm$new(table(unlist(list_array.001_fc1.5_max500$names)))
m4$setPars(estimate_pars(m4))
pvalues.001_fc1.5_max500=rev(dist_cdf(m4, lower_tail = FALSE))


m4 = dislnorm$new(table(unlist(list_array.05_fc1.5_max500$names)))
m4$setPars(estimate_pars(m4))
pvalues.05_fc1.5_max500=rev(dist_cdf(m4, lower_tail = FALSE))

m4 = dislnorm$new(table(unlist(list_array.01_fc1.5_max500$names)))
m4$setPars(estimate_pars(m4))
pvalues.01_fc1.5_max500=rev(dist_cdf(m4, lower_tail = FALSE))


hubs.1_fc1.5_max500=unique(names(table(unlist(list_array.1_fc1.5_max500$names))[names(pvalues.1_fc1.5_max500[which(pvalues.1_fc1.5_max500<0.05)])]))
hubs.05_fc1.5_max500=unique(names(table(unlist(list_array.05_fc1.5_max500$names))[names(pvalues.05_fc1.5_max500[which(pvalues.05_fc1.5_max500<0.05)])]))
hubs.01_fc1.5_max500=unique(names(table(unlist(list_array.01_fc1.5_max500$names))[names(pvalues.1_fc1.5_max500[which(pvalues.01_fc1.5_max500<0.05)])]))
hubs.001_fc1.5_max500=unique(names(table(unlist(list_array.001_fc1.5_max500$names))[names(pvalues.001_fc1.5_max500[which(pvalues.001_fc1.5_max500<0.05)])]))




indeces=t(combn(1:length(list_array.05_fc1.5_max500$names),2))
mat_int=matrix(rep(list(),length(list_array.05_fc1.5_max500$names)*length(list_array.05_fc1.5_max500$names)), length(list_array.05_fc1.5_max500$names),length(list_array.05_fc1.5_max500$names)) #number

for(i in 1:dim(indeces)[1]){
  mat_int[indeces[i,1],indeces[i,2]][[1]]=intersect(list_array.05_fc1.5_max500$names[[(indeces[i,1])]], list_array.05_fc1.5_max500$names[[(indeces[i,2])]])
}


# mat_int_num=matrix(rep(NA,length(list_array.05_fc1.5_max500$names)*length(list_array.05_fc1.5_max500$names)), length(list_array.05_fc1.5_max500$names),length(list_array.05_fc1.5_max500$names)) #number
#
# for(i in 1:dim(indeces)[1]){
#   mat_int_num[indeces[i,1],indeces[i,2]][[1]]=length(intersect(list_array.05_fc1.5_max500$names[[(indeces[i,1])]], list_array.05_fc1.5_max500$names[[(indeces[i,2])]]))
# }

gene_frequency=sort(table(unlist(list_array.05_fc1.5_max500$names)), decreasing=T)


gene_frequency_top=gene_frequency[which(gene_frequency>=gene_frequency[round(length(gene_frequency)*0.05)])]
gene_frequency_top_genes=unique(unlist(lapply(names(gene_frequency_top), function(x)(strsplit(x, split=".", fixed=T)))))

require(arrangements)

require(parallel)

# Calculate the number of cores
no_cores <- 22
mat_int_comb2=matrix(rep(list(),length(list_array.05_fc1.5_max500$names)*length(list_array.05_fc1.5_max500$names)), length(list_array.05_fc1.5_max500$names),length(list_array.05_fc1.5_max500$names)) #number


# combinazioni=parLapply(cl,list_array_pair, fun = function(x)(if(length(x)>1)(t(combn(x,2)))))
#
# combina=function(i, k){
#
#   if(isTRUE((length(mat_int[indeces[i,1],indeces[i,2]][[1]])>(k-1)))){
#
#     toadd=combinations(mat_int[indeces[i,1],indeces[i,2]][[1]],k, layout="list")
#
#     lapply(toadd,sort)
#
#   }}

combina=function(i, k, subset){
  #number of dataset where it is found
  #level of combination
  # subset of genes to restrict the permutation
  if(isTRUE((length(mat_int[indeces[i,1],indeces[i,2]][[1]])>(k-1)))){
    
    toadd=combinations(intersect(mat_int[indeces[i,1],indeces[i,2]][[1]],subset),k, layout="list")
    
    lapply(toadd,sort)
    
  }}



cl <- makeForkCluster(no_cores)
# clusterEvalQ(cl, library(arrangements))

res=parLapply(cl, 1:dim(indeces)[1], function(x)(combina(x,2, gene_frequency_top_genes)))

mat_int_comb2=matrix(rep(list(),length(list_array.05_fc1.5_max500$names)*length(list_array.05_fc1.5_max500$names)), length(list_array.05_fc1.5_max500$names),length(list_array.05_fc1.5_max500$names)) #number

for(i in 1:dim(indeces)[1]){
  if(isTRUE((length(mat_int[indeces[i,1],indeces[i,2]][[1]])>1))){
    mat_int_comb2[indeces[i,1],indeces[i,2]][[1]]=res[i][[1]]
  }
  
}

combination_common_k2=lapply(1:dim(mat_int_comb2)[1], function(i)(unique(c(do.call(c, mat_int_comb2[i,]),do.call(c, mat_int_comb2[,i])))))


all_comb2=do.call(c,combination_common_k2)

names(combination_common_k2)=names(list_array.05_fc1.5_max500$names)


lung2=parLapply(cl, combination_common_k2, length)

# save(lung2, file="lung2.RData")

all_comb2dot=parSapply(cl, all_comb2, function(x)(paste(x[1], x[2],sep=".")))



lunghezze=unlist(lung2)
names_comb=rep(NA,sum(lunghezze))
index=0

for(i in 1:length(lunghezze)){
  if(lunghezze[i]>0){
    names_comb[(index+1):(index+lunghezze[i])]=names(lunghezze[i])
    index=index+lunghezze[i]
  }
}

names(all_comb2dot)=names_comb




all_comb2_rank=sort(table(all_comb2dot), decreasing = T)



all_comb2_rank_top=all_comb2_rank[which(all_comb2_rank>=gene_frequency[round(length(gene_frequency)*0.05)])]

all_comb2_rank_top_genes=unique(unlist(lapply(names(all_comb2_rank_top), function(x)(strsplit(x, split=".", fixed=T)))))


save(all_comb2_rank_top, all_comb2_rank_top_genes, file="all_comb2_rank_top_genes.RData")

save(all_comb2_rank,all_comb2dot, file="all_comb2_rank.RData")

### combn 3
stopCluster(cl)

cl <- makeForkCluster(no_cores)
# clusterEvalQ(cl, library(arrangements))

res=parLapply(cl, 1:dim(indeces)[1], function(x)(combina(x,3, all_comb2_rank_top_genes)))


mat_int_comb3=matrix(rep(list(),length(list_array.05_fc1.5_max500$names)*length(list_array.05_fc1.5_max500$names)), length(list_array.05_fc1.5_max500$names),length(list_array.05_fc1.5_max500$names)) #number

for(i in 1:dim(indeces)[1]){
  if(isTRUE((length(mat_int[indeces[i,1],indeces[i,2]][[1]])>2))){
    mat_int_comb3[indeces[i,1],indeces[i,2]][[1]]=res[i][[1]]
  }
  
}

combination_common_k3=lapply(1:dim(mat_int_comb3)[1], function(i)(unique(c(do.call(c, mat_int_comb3[i,]),do.call(c, mat_int_comb3[,i])))))


all_comb3=do.call(c,combination_common_k3)

names(combination_common_k3)=names(list_array.05_fc1.5_max500$names)


lung3=parLapply(cl, combination_common_k3, length)

save(lung3, file="lung3.RData")

all_comb3dot=parSapply(cl, all_comb3, function(x)(paste(x[1], x[2],x[3],sep=".")))



lunghezze=unlist(lung3)
names_comb=rep(NA,sum(lunghezze))
index=0

for(i in 1:length(lunghezze)){
  if(lunghezze[i]>0){
    names_comb[(index+1):(index+lunghezze[i])]=names(lunghezze[i])
    index=index+lunghezze[i]
  }
}

names(all_comb3dot)=names_comb




all_comb3_rank=sort(table(all_comb3dot), decreasing = T)



all_comb3_rank_top=all_comb3_rank[which(all_comb3_rank>=gene_frequency[round(length(gene_frequency)*0.05)])]

all_comb3_rank_top_genes=unique(unlist(lapply(names(all_comb3_rank_top), function(x)(strsplit(x, split=".", fixed=T)))))


save(all_comb3_rank_top, all_comb3_rank_top_genes, file="all_comb3_rank_top_genes.RData")

save(all_comb3_rank,all_comb3dot, file="all_comb3_rank.RData")


###combn 4

stopCluster(cl)

cl <- makeForkCluster(no_cores)
# clusterEvalQ(cl, library(arrangements))

res=parLapply(cl, 1:dim(indeces)[1], function(x)(combina(x,4, all_comb3_rank_top_genes)))


mat_int_comb4=matrix(rep(list(),length(list_array.05_fc1.5_max500$names)*length(list_array.05_fc1.5_max500$names)), length(list_array.05_fc1.5_max500$names),length(list_array.05_fc1.5_max500$names)) #number

for(i in 1:dim(indeces)[1]){
  if(isTRUE((length(mat_int[indeces[i,1],indeces[i,2]][[1]])>3))){
    mat_int_comb4[indeces[i,1],indeces[i,2]][[1]]=res[i][[1]]
  }
  
}

combination_common_k4=lapply(1:dim(mat_int_comb4)[1], function(i)(unique(c(do.call(c, mat_int_comb4[i,]),do.call(c, mat_int_comb4[,i])))))


all_comb4=do.call(c,combination_common_k4)

names(combination_common_k4)=names(list_array.05_fc1.5_max500$names)


lung4=parLapply(cl, combination_common_k4, length)

save(lung4, file="lung4.RData")

all_comb4dot=parSapply(cl, all_comb4, function(x)(paste(x[1], x[2],x[3],x[4],sep=".")))



lunghezze=unlist(lung4)
names_comb=rep(NA,sum(lunghezze))
index=0

for(i in 1:length(lunghezze)){
  if(lunghezze[i]>0){
    names_comb[(index+1):(index+lunghezze[i])]=names(lunghezze[i])
    index=index+lunghezze[i]
  }
}

names(all_comb4dot)=names_comb




all_comb4_rank=sort(table(all_comb4dot), decreasing = T)



all_comb4_rank_top=all_comb4_rank[which(all_comb4_rank>=gene_frequency[round(length(gene_frequency)*0.05)])]

all_comb4_rank_top_genes=unique(unlist(lapply(names(all_comb4_rank_top), function(x)(strsplit(x, split=".", fixed=T)))))


save(all_comb4_rank_top, all_comb4_rank_top_genes, file="all_comb4_rank_top_genes.RData")

save(all_comb4_rank,all_comb4dot, file="all_comb4_rank.RData")