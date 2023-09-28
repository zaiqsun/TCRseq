rm(list=ls())
setwd("./")
#####convert######
name <- read.table("name.txt")
n <- 3 
TCR_split <- function(matrix,output_TRA,output_TRB,output_TRD,output_TRG){
  trust4_results <- matrix
  TRA_results <- trust4_results[which(substr(trust4_results$v,1,4)=="TRAV"),]
  TRB_results <- trust4_results[which(substr(trust4_results$v,1,4)=="TRBV"),]
  TRD_results <- trust4_results[which(substr(trust4_results$v,1,4)=="TRDV"),]
  TRG_results <- trust4_results[which(substr(trust4_results$v,1,4)=="TRGV"),]
  TRA_results$freq <- unlist(lapply(TRA_results$count,function(x) x/sum(TRA_results$count)))
  TRB_results$freq <- unlist(lapply(TRB_results$count,function(x) x/sum(TRB_results$count)))
  TRD_results$freq <- unlist(lapply(TRD_results$count,function(x) x/sum(TRD_results$count)))
  TRG_results$freq <- unlist(lapply(TRG_results$count,function(x) x/sum(TRG_results$count)))
  write.table(TRA_results,file = output_TRA,quote = F,sep = "\t",row.names = F)
  write.table(TRB_results,file = output_TRB,quote = F,sep = "\t",row.names = F)
  write.table(TRD_results,file = output_TRD,quote = F,sep = "\t",row.names = F)
  write.table(TRG_results,file = output_TRG,quote = F,sep = "\t",row.names = F)
}
metadata_all <- data.frame(file_name="",sample_id="",
                           sample_type=c(rep("Control",n),rep("Treatment",n)),
                           tcr_type="")
metadata <- data.frame(file_name="",sample_id="",
                       sample_type=c(rep("Control",n*4),rep("Treatment",n*4)),
                       tcr_type="")
count_sum <- data.frame(name=name$V1,sum=rep(0,length(name)))
for (i in 1:nrow(name)){
  input <- paste0("./",name[i,"V1"],"_report.tsv")
  output <- paste0("./",name[i,"V1"],".txt")
  trust4_results <- read.csv(input,sep="\t")
  trust4_results <- trust4_results[,1:8]
  colnames(trust4_results) <- c("count","freq","cdr3nt","cdr3aa","v","d","j","c")
  trust4_results <- trust4_results[-which(trust4_results$cdr3aa=="out_of_frame"),]
  row.names(trust4_results) <- 1:nrow(trust4_results)
  count_sum[i,"sum"] <- sum(trust4_results$count)
  write.table(trust4_results,file = output,quote = F,sep = "\t",row.names = F)
  output_TRA <- paste0("./convert_results/",name[i,"V1"],"_TRA.txt")
  output_TRB <- paste0("./convert_results/",name[i,"V1"],"_TRB.txt")
  output_TRD <- paste0("./convert_results/",name[i,"V1"],"_TRD.txt")
  output_TRG <- paste0("./convert_results/",name[i,"V1"],"_TRG.txt")
  TCR_split(trust4_results,output_TRA,output_TRB,output_TRD,output_TRG)
  metadata_all[i,"file_name"] <- output
  metadata_all[i,"sample_id"] <- strsplit(name[i,"V1"],"_")[[1]][1]
  metadata_all[i,"tcr_type"] <- "TCR"
  metadata[4*i-3,"file_name"] <- output_TRA
  metadata[4*i-3,"tcr_type"] <- "TRA"
  metadata[4*i-2,"file_name"] <- output_TRB
  metadata[4*i-2,"tcr_type"] <- "TRB"
  metadata[4*i-1,"file_name"] <- output_TRD
  metadata[4*i-1,"tcr_type"] <- "TRD"
  metadata[4*i,"file_name"] <- output_TRG
  metadata[4*i,"tcr_type"] <- "TRG"
  for (k in 3:0){
    metadata[4*i-k,"sample_id"] <- strsplit(name[i,"V1"],"_")[[1]][1]
  }
}
metadata_all$sample_id <- paste0(metadata_all$sample_id,"_",metadata_all$tcr_type)
metadata$sample_id <- paste0(metadata$sample_id,"_",metadata$tcr_type)
write.table(metadata_all,file = "./metadata_all.txt",sep="\t",row.names = F,quote = F)
write.table(metadata,file = "./metadata.txt",sep="\t",row.names = F,quote = F)
write.table(count_sum,file = "./count_sum.txt",sep="\t",row.names = F,quote = F)


#####normalizaiton#####
size=5500
name <- read.table("./name.txt")
TCR_split <- function(matrix,output_TRA,output_TRB,output_TRD,output_TRG){
  trsut4_results <- matrix
  TRA_results <- trsut4_results[which(substr(trsut4_results$v,1,4)=="TRAV"),]
  TRB_results <- trsut4_results[which(substr(trsut4_results$v,1,4)=="TRBV"),]
  TRD_results <- trsut4_results[which(substr(trsut4_results$v,1,4)=="TRDV"),]
  TRG_results <- trsut4_results[which(substr(trsut4_results$v,1,4)=="TRGV"),]
  write.table(TRA_results,file = output_TRA,quote = F,sep = "\t",row.names = F)
  write.table(TRB_results,file = output_TRB,quote = F,sep = "\t",row.names = F)
  write.table(TRD_results,file = output_TRD,quote = F,sep = "\t",row.names = F)
  write.table(TRG_results,file = output_TRG,quote = F,sep = "\t",row.names = F)
}
tcr_normalization <- function(trust4_results,size){
  pool <- c()
  for (i in 1:nrow(trust4_results)){
    pool <- append(pool,rep(i,trust4_results[i,"count"]))
  }
  all <- as.numeric(sample(pool,size = size,replace = F))
  count <- table(all)
  trust4_results <- trust4_results[unique(all),]
  trust4_results$count <- count
  trust4_results <- trust4_results[order(trust4_results$count,decreasing = T),]
  rownames(trust4_results) <- 1:nrow(trust4_results)
  return(trust4_results)
}
metadata_all_normal <- data.frame(file_name="",sample_id="",
                                  sample_type=c(rep("Con",n),rep("Bgal",n),rep("AH1",n),rep("A5",n)),
                                  tcr_type="")
metadata_normal <- data.frame(file_name="",sample_id="",
                              sample_type=c(rep("Con",n*4),rep("Bgal",n*4),rep("AH1",n*4),rep("A5",n*4)),
                              tcr_type="")
for (i in 1:nrow(name)){
  input <- paste0("./trust4/",name[i,"V1"],"_report.tsv")
  output <- paste0("./convert_results_normal/",name[i,"V1"],".txt")
  trust4_results <- read.csv(input,sep="\t")
  trust4_results <- trust4_results[,1:8]
  colnames(trust4_results) <- c("count","freq","cdr3nt","cdr3aa","v","d","j","c")
  trust4_results <- trust4_results[-which(trust4_results$cdr3aa=="out_of_frame"),]
  row.names(trust4_results) <- 1:nrow(trust4_results)
  trust4_results <- tcr_normalization(trust4_results,size)
  write.table(trust4_results,file = output,quote = F,sep = "\t",row.names = F)
  output_TRA <- paste0("./convert_results_normal/",name[i,"V1"],"_TRA.txt")
  output_TRB <- paste0("./convert_results_normal/",name[i,"V1"],"_TRB.txt")
  output_TRD <- paste0("./convert_results_normal/",name[i,"V1"],"_TRD.txt")
  output_TRG <- paste0("./convert_results_normal/",name[i,"V1"],"_TRG.txt")
  TCR_split(trust4_results,output_TRA,output_TRB,output_TRD,output_TRG)
  metadata_all_normal[i,"file_name"] <- output
  metadata_all_normal[i,"sample_id"] <- strsplit(name[i,"V1"],"_")[[1]][1]
  metadata_all_normal[i,"tcr_type"] <- "TCR"
  metadata_normal[4*i-3,"file_name"] <- output_TRA
  metadata_normal[4*i-3,"tcr_type"] <- "TRA"
  metadata_normal[4*i-2,"file_name"] <- output_TRB
  metadata_normal[4*i-2,"tcr_type"] <- "TRB"
  metadata_normal[4*i-1,"file_name"] <- output_TRD
  metadata_normal[4*i-1,"tcr_type"] <- "TRD"
  metadata_normal[4*i,"file_name"] <- output_TRG
  metadata_normal[4*i,"tcr_type"] <- "TRG"
  for (k in 3:0){
    metadata_normal[4*i-k,"sample_id"] <- strsplit(name[i,"V1"],"_")[[1]][1]
  }
}
metadata_all_normal$sample_id <- paste0(metadata_all_normal$sample_id,"_",metadata_all_normal$tcr_type)
metadata_normal$sample_id <- paste0(metadata_normal$sample_id,"_",metadata_normal$tcr_type)
write.table(metadata_all_normal,file = "./metadata_all_normal.txt",sep="\t",row.names = F,quote = F)
write.table(metadata_normal,file = "./metadata_normal.txt",sep="\t",row.names = F,quote = F)

#####Basic#####
inputpath <- "./results/"
outputpath <- "./results/"
df <- read.csv(paste0(inputpath,"diversity.strict.exact.txt"), sep = "\t",stringsAsFactors = FALSE)
df$real_diversity <- log(df$shannonWienerIndex_mean)
df$clonality <- 1 - df$normalizedShannonWienerIndex_mean
for (k in c("TRA","TRB","TRD","TRG")){
  work_df <- df[which(df$tcr_type==k),]
  sample_name <- work_df$sample_id
  richness<- work_df$observedDiversity_mean
  diversity <- work_df$real_diversity
  clonality <- work_df$clonality
  work_df <- data.frame(sample_name,richness,diversity,clonality)
  write.table(work_df,file = paste0(outputpath,k,"-Basic.txt"),row.names = F,quote = F,sep="\t")
}

#?ֶ?????
metadata <- read.csv("./metadata_all.txt",sep = "\t")
Shannon_index <- function(list_p){
  sum=0
  for (p in list_p){
    sum <- p*log(p) + sum
  }
  H <- -sum
  return(H)
}
clonality <- function(diversity,richness){
  E<- diversity/log(richness)
  C <- 1-E
  return(C)
}
file_path <- metadata$file_name
results <- data.frame(sample=metadata$sample_id)
for (i in 1:length(file_path)){
  df <- read.csv(file = file_path[i],
                 stringsAsFactors = F,sep = "\t")
  temp <- try(list_p <- df[,"freq"])
  if ('try-error' %in% class(temp)){next}
  results[i,"richness"] <- nrow(df)
  results[i,"diversity"] <- Shannon_index(list_p)
  results[i,"clonality"] <- 1-clonality(Shannon_index(list_p),nrow(df))
}
write.table(results,"./results/basic.txt",sep = "\t",row.names = F,quote = F)

#####Basic-plot#####
rm(list = ls())
inputpath <- "./results/"
outputpath <- "./output/"
library(ggplot2)
library(reshape2)
library(ggpubr)
df <- read.csv(paste0(inputpath,"diversity.strict.exact.txt"), sep = "\t",stringsAsFactors = FALSE)
df$real_diversity <- log(df$shannonWienerIndex_mean)
df$clonality <- 1 - df$normalizedShannonWienerIndex_mean
bcr_type <- "TRB"
work_df <- df[which(df$tcr_type==bcr_type),]
sample_name <- work_df$sample_id
Group <- work_df$sample_type
richness<- work_df$observedDiversity_mean
diversity <- work_df$real_diversity
clonality <- work_df$clonality
work_df <- data.frame(Group,richness,diversity,clonality)
work_df$Group <- factor(work_df$Group)  
work_res<- melt(work_df, id.vars = "Group")
work_richness <- work_res[which(work_res$variable=="richness"),]
work_diversity <- work_res[which(work_res$variable=="diversity"),]
work_clonality <- work_res[which(work_res$variable=="clonality"),]

p1 <- ggboxplot(work_diversity, x = "Group", y = "value",
                legend = "right",color = "Group", palette = "nejm",
                add = "jitter",add.params=list(color = "Group",size=3),
                title = "Diversity",  xlab = "",size = 0.5,width =0.5,
                ylab = paste0("Shannon Diversity Index of \n ",bcr_type," repertoire")) +
  ylim(c(0,max(work_diversity$value)+min(work_diversity$value)))+
  theme(axis.text.x =element_text(size=13), axis.title.y=element_text(size=16,face="bold"))
p1 + stat_compare_means(aes(group = work_df$Group),
                        label.x = 1.5,
                        label.y = max(work_diversity$value)+min(work_diversity$value),
                        label = "p.signif",
                        method = 'wilcox.test',hide.ns = TRUE,size=8)  #p.format
ggsave(paste0(outputpath,paste0(bcr_type,"-Diversity.pdf")),dpi = 600,width = 8,height = 8)
dev.off()
p2 <- ggboxplot(work_richness, x = "Group", y = "value",
                legend = "right",color = "Group", palette = "nejm",
                add = "jitter",add.params=list(color = "Group",size=3),
                title = "Richness",  xlab = "",size = 0.5,width =0.5,
                ylab = paste0("Richness of \n ",bcr_type," repertoire"))+
  ylim(c(0,max(work_richness$value)+min(work_richness$value)))+
  theme(axis.text.x =element_text(size=13), axis.title.y=element_text(size=16,face="bold"))
  
p2 + stat_compare_means(aes(group =work_df$Group),
                        label = "p.signif",label.x = 1.5,
                        label.y = max(work_richness$value)+min(work_richness$value),
                        method = 'wilcox.test',hide.ns = TRUE,size=8) #p.format
ggsave(paste0(outputpath,paste0(bcr_type,"-Richness.pdf")),dpi = 600,width = 8,height = 8)
dev.off()
p3 <- ggboxplot(work_clonality, x = "Group", y = "value",
                legend = "right",color = "Group", palette = "nejm",
                add = "jitter",add.params=list(color = "Group",size=3),
                title = "Clonality",  xlab = "",size = 0.5,width =0.5,
                ylab = paste0("Clonality Index of \n ",bcr_type," repertoire"))+
  ylim(c(0,max(work_clonality$value)+min(work_clonality$value)))+
  theme(axis.text.x =element_text(size=13), axis.title.y=element_text(size=16,face="bold"))
p3 + stat_compare_means(aes(group =work_df$Group),  
                        label = "p.signif",label.x = 1.5,
                        label.y = max(work_clonality$value)+min(work_clonality$value),
                        method = 'wilcox.test',hide.ns = TRUE,size=8) #p.format
ggsave(paste0(outputpath,paste0(bcr_type,"-Clonality.pdf")),dpi = 600,width = 8,height = 8)
dev.off()



#####CDR3_length_distribution#####
rm(list = ls())
library(ggplot2)
library(reshape2)
library(ggpubr)
lenth_distribution <- function(filepath,all_length,id){
  df <- read.csv(filepath,sep="\t",stringsAsFactors = F)
  for (i in 1:nrow(df)){
    cdr3 <- df[i,"cdr3aa"]
    all_length[as.character(nchar(cdr3)),id] <- all_length[as.character(nchar(cdr3)),id]+1
  }
  return(all_length)
}
all_lenth_distribution <- function(file_list,id_list,group_list){
  all_length <- data.frame()
  for (i in 1:length(file_list)){
    df <- read.csv(file_list[i],sep="\t",stringsAsFactors = F)
    for (k in 1:nrow(df)){
      cdr3 <- df[k,"cdr3aa"]
      all_length[as.character(nchar(cdr3)),id_list[i]] <- 0
    }
  }
  all_length[is.na(all_length)] <- 0
  for (i in 1:length(file_list)){
    all_length <- lenth_distribution(file_list[i],all_length,id_list[i])
    all_length["sample_type",id_list[i]]=group_list[i]
  }
  return(all_length)
}
name <- read.table("./metadata.txt",sep = "\t",header = 1)
work_tcr_type <- "TRB"
file_list <- name$file_name[which(name$tcr_type==work_tcr_type)]
id_list <- name$sample_id[which(name$tcr_type==work_tcr_type)]
group_list <- c(rep("Control",4),rep("Treatment",4))
all_length <- all_lenth_distribution(file_list,id_list,group_list)
#---------absolute----------------#
if(F){
temp <- all_length[-nrow(all_length),]
temp <- temp[order(as.numeric(rownames(temp))),]
all_length <- rbind(temp,all_length[nrow(all_length),])}
####-------------percent----------###
if(T){
temp <- all_length[-nrow(all_length),]
temp=as.data.frame(lapply(temp,as.numeric),row.names = rownames(temp))
colSums(temp)
temp <- data.frame(apply(temp,2,function(x) x/sum(x)))
temp <- temp[order(as.numeric(rownames(temp))),]
all_length <- rbind(temp,all_length[nrow(all_length),])}
####------------------------------###

####---check this and change--###
all_length <- all_length[c(6:14,nrow(all_length)),]
####--------------------------###
##draw picture##
df <- all_length
work_df <- data.frame(t(df))
df <- melt(work_df, id="sample_type", variable.name="length",value.name = "count")
df$length <- substring(df$length,2)
df$count <- as.numeric(df$count)
df$length <- as.numeric(df$length)
df <- df[order(df$length),]
df$length <- as.factor(df$length)
mean <- aggregate(df$count, by=list(df$length,df$sample_type), FUN=mean)
sd <- aggregate(df$count, by=list(df$length, df$sample_type), FUN=sd)
len <- aggregate(df$count, by=list(df$length, df$sample_type), FUN=length)
df_res <- data.frame(mean, sd=sd$x, len=len$x)
colnames(df_res) = c("CDR3_Length", "Group", "Mean", "Sd", "Count")
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
ggplot(df_res, mapping = aes(x=CDR3_Length, y=Mean, fill=Group,group=Group,color=Group)) +
  labs(xlab = "CDR3 Length",ylab="Mean",title =paste0("Distribution of ",work_tcr_type," CDR3 Length"))+
  geom_bar(stat="identity", position=position_dodge(),color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean-Se, ymax=Mean +Se),
                position=position_dodge(.8), width=.2,color="black") +
 # scale_y_continuous(expand=c(0,0))+
 # coord_cartesian(ylim = c(0, 31000))+
  theme(axis.text.x = element_text(size = 14, color = "black"))+##????x????????С
  theme(axis.text.y = element_text(size = 14, color = "black"))+##????y????????С
  theme(title=element_text(size=13))+#???ñ?????????С
  theme_bw()+
  geom_line(size=0.5,position = position_dodge(0.8))+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.5)) #?޸?Y?᷶Χ
ggsave(paste0("./output/CDR3_lenth",work_tcr_type,".pdf"),width = 10,height = 8)

#?ѵ?ͼ
library(ggsci)
ggplot(df_res,mapping = aes(x=Group,y=Mean,fill=CDR3_Length))+
  geom_col(position = 'stack', width = 0.6)+
  scale_color_aaas()+
  theme_bw()+scale_y_continuous(expand = c(0,0))
ggsave(paste0("./output/CDR3_lenth",work_tcr_type,"_stack.pdf"),width = 7,height = 7) 




#####VJ-Usage#####
rm(list=ls())
library(RColorBrewer)
library(pheatmap)
for (region in c("V","J")){
  for (type in c("TRA","TRB","TRD","TRG")){
    file_name <- paste0("./results/segments.wt.",region,".txt")
    df <- read.csv(file_name,sep="\t",row.names = 1)
    work_tcr_type <- type
    work_df <- df[which(df$tcr_type==work_tcr_type),]
    work_df <- work_df[,-1:-2]
    work_df <- t(work_df)
    work_df <- work_df[which(rowSums(work_df)>0),]
    work_df <- work_df*100
    #---------------------------------------------------------#
    annotation_col <- data.frame(c(rep("Control",4),rep("Treatment",4)))
    #---------------------------------------------------------#
    rownames(annotation_col) <- colnames(work_df)
    colnames(annotation_col) <- "Type"
    #coul <- colorRampPalette(brewer.pal(9, "OrRd"))(50)
    coul <-  colorRampPalette(colors = c("blue","white","red"))(100)
    pheatmap(work_df,
             filename = paste0("./output/",work_tcr_type,"_GeneUsage_",region,".pdf"),
             scale = "row",
             annotation_col = annotation_col,cluster_cols = FALSE,
             height = 15,width = 10,cutree_col =2,display_numbers = FALSE,labels_col = "",color = coul)
  }
}

#####aa-fre#####
metadata <- read.table("./metadata_normal.txt",sep="\t",header = 1)
chr_freq <- function(str,chr){
  df <- data.frame()
  for (i in 1:nchar(str)){
    df[i,"cdr3"] <- substr(str,i,i)
  }
  p <- table(df)[chr]/nchar(str)
  return(p)
}
aa_freq <- function(df,chr){
  p <- c()
  for (i in 1:nrow(df)){
    cdr3 <- df[i,"cdr3aa"]
    p[i] <- chr_freq(cdr3,chr)
  }
  p_mean <- mean(p,na.rm = T)
  return(p_mean)
}
aa_name <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
aa_name <- "K"
file_list <- metadata$file_name
sample_id <- metadata$sample_id
p=data.frame()
for (i in 1:length(file_list)){
  print(paste0(file_list[i],"---running!"))
  df <- read.csv(file_list[i],sep="\t",stringsAsFactors = F)
  for (aa in aa_name){
    print(paste0(aa,"-------runnning!----------"))
    p[sample_id[i],aa] <- aa_freq(df,aa)
  }
}
p1 <- p*100
write.table(p1,file="./results/aa_frequency.txt",sep="\t",row.names = T,quote = F)

library(ggplot2)
library(reshape2)
library(ggpubr)
aa_fre <- read.csv("./results/aa_frequency.txt",sep="\t",row.names = 1)
aa_fre$Sample_type <- c(rep("WT",24),rep("KO",24))
aa_fre$tcr_type <- rep( c("TRA","TRB","TRD","TRG"),4)
aa_fre$Sample_type <- factor(aa_fre$Sample_type)
aa_fre$tcr_type <- factor(aa_fre$tcr_type)

for (tcr in c("TRA","TRB","TRD","TRG")){
  aa_fre_work <- aa_fre[which(aa_fre$tcr_type==tcr),]
  aa_fre_work <- aa_fre_work[,-ncol(aa_fre_work)]
  work_res <- melt(aa_fre_work,id.vars = "Sample_type")
  for (i in 1:(ncol(aa_fre_work)-1)){
    temp <- work_res[which(work_res$variable==colnames(aa_fre)[i]),]
    p1 <- ggpaired(temp, x = "Sample_type", y = "value",
                   line.color = "gray", line.size = 0.5,
                    legend = "right",color = "Sample_type", palette = "nejm",
                    add = "jitter",add.params=list(color = "Sample_type",size=3),
                    title = colnames(aa_fre)[i],  xlab = "",size = 0.5,width =0.5,
                    ylab = paste0("Percentage of ",tcr," aa frequency")) +
      theme(axis.text.x =element_text(size=13), axis.title.y=element_text(size=16,face="bold"))+
      theme(legend.position = 'none')
    p1 + stat_compare_means(aes(group = Sample_type),
                            label.x = 1.5,label = "p.signif",
                           # comparisons = list(c("WT","KO")),
                            method = 'wilcox.test',hide.ns = TRUE,size=8)
    ggsave(paste0("./output/aa_fre/",tcr,"_",colnames(aa_fre_work)[i],".png"))
  }
}






#####VJ-Paired#####
rm(list=ls())
library(circlize);library(RColorBrewer);library(pheatmap)
inputdir <- "./results/VJPAIR/"
outputdir <- "./output/VJPAIR/"
#*********************************individual*********************************#
VJpair_circle <- function(input,output){
  args <- c(input,output)
  file_in  <- args[1]
  file_out <- args[2]
  
  
  
  # load data and preproc to fit formats
  
  temp <- read.table(file_in, sep="\t", comment="")
  n <- nrow(temp)
  m <- ncol(temp)
  rn = as.character(temp[2:n,1])
  cn = apply(temp[1,2:m], 2 , as.character)
  mat <- matrix(apply(temp[2:n, 2:m], 1:2, as.numeric), n - 1, m-1) * 100
  
  n <- nrow(temp)
  m <- ncol(temp)
  
  # Here columns and rows correspond to V and J segments respectively
  # Also replace possible duplicates (undef, '.', ...)
  
  duplicates <- intersect(rn, cn)
  
  rownames(mat) <- replace(rn, rn==duplicates, paste("V", duplicates, sep=""))
  colnames(mat) <- replace(cn, cn==duplicates, paste("J", duplicates, sep=""))
  
  # sort
  
  col_sum = apply(mat, 2, sum)
  row_sum = apply(mat, 1, sum)
  
  mat <- mat[order(row_sum), order(col_sum)]
  
  # equal number of characters for visualizaiton
  
  rn <- rownames(mat)
  cn <- colnames(mat)
  
  maxrn <- max(nchar(rn))
  maxcn <- max(nchar(cn))
  
  for(i in seq_len(length(rn))) {
    rn[i] <- paste(rn[i], paste(rep(" ", maxrn - nchar(rn[i])), collapse = ''))
  }
  
  for(i in seq_len(length(cn))) {
    cn[i] <- paste(cn[i], paste(rep(" ", maxcn - nchar(cn[i])), collapse = ''))
  }
  
  rownames(mat) <- rn
  colnames(mat) <- cn
  
  # viz using circlize
  
  if (grepl("\\.pdf$",file_out)){
    pdf(file_out,width=10,height=10)
  } else if (grepl("\\.png$",file_out)) {
    png(file_out, width     = 3.25,
        height    = 3.25,
        units     = "in",
        res       = 1200,
        pointsize = 4)
  } else {
    stop('Unknown plotting format')
  }
  
  circos.par(gap.degree = c(rep(1, nrow(mat)-1), 10, rep(1, ncol(mat)-1), 15), start.degree = 5)
  
  rcols <- rep(brewer.pal(12, "Paired"), nrow(mat)/12 + 1)[1:nrow(mat)]
  ccols <- rep(brewer.pal(12, "Paired"), ncol(mat)/12 + 1)[1:ncol(mat)]
  
  names(rcols) <- rev(rownames(mat))
  names(ccols) <- rev(colnames(mat))
  
  chordDiagram(mat, annotationTrack = "grid",reduce=0,
               grid.col = c(rcols, ccols),
               preAllocateTracks = list(track.height = 0.2), transparency = 0.5)
  
  circos.trackPlotRegion(track.index = 1, bg.border = NA,
                         panel.fun = function(x, y) {
                           sector.name = get.cell.meta.data("sector.index")
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           circos.text(mean(xlim), ylim[1], cex = 0.5, sector.name, facing = "clockwise", adj = c(0, 0.5))
                         }
  )
  
  circos.clear()
  dev.off()
  print(paste0(args[1],"---have finished"))
}
VJpair_heat <- function(input,output){
  args <- c(input,output)
  file_in  <- args[1]
  file_out <- args[2]
  temp <- read.csv(file_in, sep="\t", comment="",row.names = 1)
  work_df <- temp
  work_df <- work_df[rowSums(work_df)!=0,]
  pheatmap(work_df,filename = file_out,
           cluster_cols = T,  cluster_rows = T,height = 15,width = 15,
           display_numbers = FALSE)
}
#name.txt ls *txt>name.txt
name <- read.table("./results/VJPAIR/name.txt")
file_list <- name[,"V1"][1:nrow(name)-1]
metadata<- read.table("./metadata.txt",sep = "\t",header = 1)
metadata <- metadata[which(metadata$tcr_type=="TRA" | metadata$tcr_type=="TRB"),]
out_list <- paste0(metadata$sample_id,"_VJpaired_.pdf")
length(file_list)==length(out_list)
for (i in 1:length(file_list)){
  VJpair_circle(paste0(inputdir,file_list[i]),
                paste0(outputdir,out_list[i]))
  VJpair_heat(paste0(inputdir,file_list[i]),
              paste0(outputdir,"heat_",out_list[i]))
}

#****************************************************************************************************#
#*********************************average pair*********************************#
ave_VJ <- function(file_list,output,heat_output){
  file_out <- output
  all_V <- c()
  all_J <- c()
  for (i in 1:length(file_list)){
    temp <- read.table(paste0(inputdir,file_list[i]), sep="\t", comment="")
    #temp <- rbind(temp[1,],temp[grepl("TRBJ2",temp[,1]),])
    #temp <- cbind(temp[,1],temp[,grepl("TRBV13",temp[1,])])
    V <- as.character(temp[1,2:ncol(temp)])
    J <- as.character(temp[2:nrow(temp),1])
    all_V <- append(all_V,V)
    all_J <- append(all_J,J)
  }
  pair_matrix <- data.frame(matrix(0,nrow = length(unique(all_J))+1, ncol = length(unique(all_V))+1))
  pair_matrix[1,2:ncol(pair_matrix)] <- unique(all_V)
  pair_matrix[2:nrow(pair_matrix),1] <- unique(all_J)
  rownames(pair_matrix) <- pair_matrix[,1]
  colnames(pair_matrix) <- pair_matrix[1,]
  pair_matrix <- pair_matrix[-1,-1]
  pair_matrix=as.data.frame(lapply(pair_matrix,as.numeric),row.names = rownames(pair_matrix))
  for (i in 1:length(file_list)){
    temp <- read.table(paste0(inputdir,file_list[i]), sep="\t", comment="")
    #temp <- rbind(temp[1,],temp[grepl("TRBJ2",temp[,1]),])
    #temp <- cbind(temp[,1],temp[,grepl("TRBV13",temp[1,])])
    rownames(temp) <- temp[,1]
    colnames(temp) <- temp[1,]
    temp <- temp[-1,-1]
    temp=as.data.frame(lapply(temp,as.numeric),row.names = rownames(temp))
    for (i in 1:nrow(temp)){
      for (k in 1:ncol(temp)){
        pair_matrix[rownames(temp)[i],colnames(temp)[k]] <- temp[i,k]+pair_matrix[rownames(temp)[i],colnames(temp)[k]]
      }
    }
  }
  pair_matrix <- pair_matrix[which(rowSums(pair_matrix)!=0),]
  pair_matrix <- pair_matrix[,which(colSums(pair_matrix)!=0)]
  pair_matrix <- pair_matrix/length(file_list)
  work_heat <- pair_matrix
  pair_matrix <- cbind(rownames(pair_matrix),pair_matrix)
  pair_matrix <- rbind(colnames(pair_matrix),pair_matrix)
  write.table(pair_matrix,file = paste0(output,".txt"),sep = "\t",quote = F,row.names = F,col.names = F)
  temp <- pair_matrix
  n <- nrow(temp)
  m <- ncol(temp)
  rn = as.character(temp[2:n,1])
  cn = apply(temp[1,2:m], 2 , as.character)
  mat <- matrix(apply(temp[2:n, 2:m], 1:2, as.numeric), n - 1, m-1) * 100
  
  n <- nrow(temp)
  m <- ncol(temp)
  
  # Here columns and rows correspond to V and J segments respectively
  # Also replace possible duplicates (undef, '.', ...)
  
  duplicates <- intersect(rn, cn)
  
  rownames(mat) <- replace(rn, rn==duplicates, paste("V", duplicates, sep=""))
  colnames(mat) <- replace(cn, cn==duplicates, paste("J", duplicates, sep=""))
  
  # sort
  
  col_sum = apply(mat, 2, sum)
  row_sum = apply(mat, 1, sum)
  # subset
  mat <- mat[order(row_sum,decreasing = T), order(col_sum,decreasing = T)]
  #mat <- mat[1:12,1:24]
  
  #mat <- mat[order(row_sum), order(col_sum)]
  
  # equal number of characters for visualizaiton
  
  rn <- rownames(mat)
  cn <- colnames(mat)
  
  maxrn <- max(nchar(rn))
  maxcn <- max(nchar(cn))
  
  for(i in seq_len(length(rn))) {
    rn[i] <- paste(rn[i], paste(rep(" ", maxrn - nchar(rn[i])), collapse = ''))
  }
  
  for(i in seq_len(length(cn))) {
    cn[i] <- paste(cn[i], paste(rep(" ", maxcn - nchar(cn[i])), collapse = ''))
  }
  
  rownames(mat) <- rn
  colnames(mat) <- cn
  
  # viz using circlize
  
  if (grepl("\\.pdf$",file_out)){
    pdf(file_out,width=10,height=10)
  } else if (grepl("\\.png$",file_out)) {
    png(file_out, width     = 3.25,
        height    = 3.25,
        units     = "in",
        res       = 1200,
        pointsize = 4)
  } else {
    stop('Unknown plotting format')
  }
  
  circos.par(gap.degree = c(rep(1, nrow(mat)-1), 10, rep(1, ncol(mat)-1), 15), start.degree = 5)
  
  rcols <- rep(brewer.pal(12, "Paired"), nrow(mat)/12 + 1)[1:nrow(mat)]
  ccols <- rep(brewer.pal(12, "Paired"), ncol(mat)/12 + 1)[1:ncol(mat)]
  
  names(rcols) <- rev(rownames(mat))
  names(ccols) <- rev(colnames(mat))
  
  chordDiagram(mat, annotationTrack = "grid",reduce=0,
               grid.col = c(rcols, ccols),
               preAllocateTracks = list(track.height = 0.2), transparency = 0.5)
  
  circos.trackPlotRegion(track.index = 1, bg.border = NA,
                         panel.fun = function(x, y) {
                           sector.name = get.cell.meta.data("sector.index")
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           circos.text(mean(xlim), ylim[1], cex = 0.5, sector.name, facing = "clockwise", adj = c(0, 0.5))
                         }
  )
  
  circos.clear()
  dev.off()
  pheatmap(work_heat,filename = heat_output,
           color = colorRampPalette(colors = c("blue","white","red"))(100),
           cluster_cols = F,  cluster_rows = F,height = 15,width = 15,
           display_numbers = FALSE)
  
}
metadata<- read.table("./metadata.txt",sep = "\t",header = 1)
metadata$file_name <- lapply(metadata$file_name, function(x) paste0(strsplit(x,"/")[[1]][3],".fancyvj.wt.txt"))
#*******************************#
tcr_type="TRB";sample_type="KO"
#*******************************#
file_list <- metadata[which(metadata$tcr_type==tcr_type & metadata$sample_type==sample_type),]$file_name
output <- paste0(sample_type,"_",tcr_type,"_VJpaired_.pdf")
ave_VJ(file_list,paste0(outputdir,output),paste0(outputdir,"heat_",output))



#****************************************************************************************************#
#*********************************special pair*********************************#
metadata<- read.table("./metadata.txt",sep = "\t",header = 1)
metadata$file_name <- lapply(metadata$file_name, function(x) paste0(strsplit(x,"/")[[1]][3],".fancyvj.wt.txt"))
#*******************************#
tcr_type="TRB";sample_type="MTCR2"
#*******************************#
file_list <- metadata[which(metadata$tcr_type==tcr_type & metadata$sample_type==sample_type),]$file_name
output <- paste0(sample_type,"_",tcr_type,"V13_VJpaired_.pdf")
res <- data.frame(pair=c("V17-J2.1","V29-J2.1","V14-J2.1","V19-J2.1","V13.1-J2.7"))
for (i in 1:length(file_list)){
  temp <- read.table(paste0(inputdir,file_list[i]), sep="\t", comment="")
  v <- c("TRBV17","TRBV29","TRBV14","TRBV19","TRBV13-1")
  j <- c("TRBJ2-1","TRBJ2-1","TRBJ2-1","TRBJ2-1","TRBJ2-7")
  sample_id <- metadata$sample_id[i]
  fancy <- c()
  for (i in 1:5){
    fancy[i] <- temp[which(temp$V1==j[i]),which(temp[1,]==v[i])]
  }
  res <- cbind(res,data.frame(M = fancy))
}
assign(sample_type,res)
work_df <- cbind(MTCR1,MTCR2)
rownames(work_df) <- work_df$pair
work_df <- work_df[,c(-1,-6)]
colnames(work_df) <- metadata$sample_id
work_df=as.data.frame(lapply(work_df,as.numeric),row.names = rownames(work_df))


pheatmap(work_df,filename = "./output/sepcial_heatmap.pdf",
        # color = colorRampPalette(colors = c("blue","white","red"))(100),
         cluster_cols = F,  cluster_rows = F,height = 15,width = 15,
         display_numbers = FALSE)


#####CDR3aa_motif#####
rm(list=ls())
library(ggplot2)
library(ggseqlogo)
library(dplyr)
draw_aa_fre <- function(input,output,n){
  df <- read.csv(input,sep="\t",stringsAsFactors = F)
  df <- df %>% group_by(cdr3aa) %>% dplyr::summarise(count=sum(count)) %>% data.frame()
  df$length <- nchar(df[,"cdr3aa"])
  cdr3_length <- df %>% group_by(length) %>% dplyr::summarise(n=sum(count),richness=n())
  cdr3_length <- cdr3_length[order(cdr3_length$richness,decreasing = T),]
  list_name <- paste0(cdr3_length$richness,",",cdr3_length$n)
  cdr3_length <- cdr3_length[1:n,"length"]
  cdr3_length <- unlist(lapply(cdr3_length[1],function(x) as.numeric(x)))
  aa_list <- c()
  for (n in 1:length(cdr3_length)){
    cdr3aa <- df[which(nchar(df[,"cdr3aa"])==cdr3_length[n]),"cdr3aa"]
    aa_list[[n]] <- unlist(lapply(cdr3aa, function(x) rep(x,df[which(df$cdr3aa==x),]$count)))
    names(aa_list)[n] <- list_name[n]
  }
  ggseqlogo(aa_list,method="prob")
  ggsave(output,height = 7.98,width = 20)
}
metadata <- read.table("./metadata.txt",sep = "\t",header = 1)
file_list   <- metadata$file_name
name <- metadata$sample_id
n <- 9 #number of CDR3aa length ???????ֳ???
for (i in 1:length(file_list)){
  draw_aa_fre(file_list[i],paste0("./output/aa_motif/",name[i],".pdf"),n)
}
#****************************************************************************************************#
#average motif
ave_motif <- function(file_list,output,form){
  all_matrix <- data.frame()
  for (i in 1:length(file_list)){
    df <- read.csv(file_list[i],sep="\t",stringsAsFactors = F)
    all_matrix <- rbind(all_matrix,df)
  }
  
  cdr3_length <- sort(table(nchar(all_matrix[,"cdr3aa"])),decreasing = T)
  cdr3_length <- cdr3_length[1:9]
  cdr3_length <- as.numeric(names(cdr3_length))
  cdr3_length <- sort(cdr3_length)
  aa_list <- c()
  for (n in 1:length(cdr3_length)){
    aa_list[[n]] <- all_matrix[which(nchar(all_matrix[,"cdr3aa"])==cdr3_length[n]),"cdr3aa"]
  }
  if (form=="up"){
    for (i in 1:length(aa_list)){
      ggseqlogo(aa_list[[i]],method = "prob",rev_stack_order=F)+
        theme(axis.text.y = element_blank())+ylab("")+theme(legend.position = "None")+
        theme(axis.text.x = element_text(size=15,face = "bold"))
      ggsave(paste0(output,nchar(aa_list[[i]][1]),".pdf"),height = 7.98,width = 11.88)
    }
  }else if (form == "down"){
    for (i in 1:length(aa_list)){
      ggseqlogo(aa_list[[i]],method = "prob",rev_stack_order=T)+
        theme(axis.text.y = element_blank())+ylab("")+theme(axis.text.x = element_blank())
      ggsave(paste0(output,nchar(aa_list[[i]][1]),".pdf"),height = 7.98,width = 11.88)
    }
  }
  
}
metadata<- read.table("./metadata.txt",sep = "\t",header = 1)
#****************************************#
sample_type="Severe";
tcr_type="TRB";form="up"
#****************************************#
file_list <- metadata[which(metadata$tcr_type==tcr_type & metadata$sample_type==sample_type),]$file_name
output <- paste0(sample_type,"_",tcr_type,"_aamotif_")
ave_motif(file_list,paste0("./output/aa_motif/",output),form)
#****************************************************************************************************#



ggseqlogo(aa_list[[i]],method = "prob",rev_stack_order=T)+
  theme(axis.text.y = element_blank())+ylab("")+theme(axis.text.x = element_blank())




#####aa_character#####
rm(list=ls())
library(ggplot2)
library(reshape2)
name <- read.table("./immu/results/metadata_phy.txt",sep = "\t",header = 1)
aa_character <- read.table("../IMGT/aa_feature.txt",sep = "\t",header = 1)
get_all_cdr3aa <- function(file_path,most_frq_length){
  temp <- read.table(file_path,sep = "\t",header = 1)
  all_cdr3 <- c()
  for (i in 1:nrow(temp)){
    if (nchar(temp$cdr3aa[i])==most_frq_length & !("_" %in% strsplit(temp$cdr3aa[i],"")[[1]]) & !("?" %in% strsplit(temp$cdr3aa[i],"")[[1]]) ){
      all_cdr3 <- append(all_cdr3,temp$cdr3aa[i])
    }
  }
  return(all_cdr3)
}
#---------------------------------------------------------#
tcr.type <- "Middle";most_frq_length <- 14
#---------------------------------------------------------#
for (sample.type in c("spc","com")){
  file_list <- name[which(name$tcr_type==tcr.type & name$sample_type==sample.type),]$file_name
  all_cdr3 <- c()
  #get all cdr3
  for (i in 1:length(file_list)){
    file_path <- file_list[i]
    all_cdr3 <- append(all_cdr3,get_all_cdr3aa(file_path,most_frq_length))
  }
  #average character index
  aa_feature_hy <- data.frame();aa_feature_vol <- data.frame();aa_feature_pI <- data.frame()
  for (i in 1:length(all_cdr3)){
    for (k in 1:nchar(all_cdr3[i])) {
      aa_feature_hy[i,k] <- aa_character[which(aa_character$Abbreviations==substr(all_cdr3[i],k,k)),]$Hydropathyindex
      aa_feature_vol[i,k] <- aa_character[which(aa_character$Abbreviations==substr(all_cdr3[i],k,k)),]$Volume
      aa_feature_pI[i,k] <- aa_character[which(aa_character$Abbreviations==substr(all_cdr3[i],k,k)),]$pI
    }
  }
  assign(paste0(sample.type,"_Hydrophobicity"),colSums(aa_feature_hy)/nrow(aa_feature_hy))
  assign(paste0(sample.type,"_vol"),colSums(aa_feature_vol)/nrow(aa_feature_vol))
  assign(paste0(sample.type,"_pI"),colSums(aa_feature_pI)/nrow(aa_feature_pI))
}
Hydrophobicity <- com_Hydrophobicity/spc_Hydrophobicity
vol <- com_vol/spc_vol
pI <- com_pI/spc_pI
cdr3_character <- data.frame(Hydropathyindex=Hydrophobicity,Volume=vol,pI=pI,position = 1:most_frq_length)
df <- melt(cdr3_character, id="position", variable.name="Feature",value.name = "value")
ggplot(data = df, mapping = aes(x = position, y = value, group = Feature,color=Feature)) + geom_line(size=1)+
  labs(title = "Middle Physiochemical feature")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+##????x????????С
  theme(axis.text.y = element_text(size = 14, color = "black"))+##????y????????С
  theme(title=element_text(size=13))+#???ñ?????????С
  theme_bw()+xlab("TCR Position")+ylab("com vs spc")+
  scale_x_continuous(breaks=df$position) 
ggsave(paste0("./output/Physiochemical_feature_",tcr.type,".png"),dpi = 600)



#####VJPAIR_FAMILY#####
rm(list=ls())
library(circlize);library(RColorBrewer);library(pheatmap);library(dplyr)
inputdir <- "./results/VJPAIR/"
outputdir <- "./results/VJPAIR/family/"
v_j_family <- function(file,outfile){
df <- read.csv(paste0(file),sep="\t")
rownames(df) <- df$.
df <- df[,-1]
TRBJ_name <- rownames(df);TRBV_name <- colnames(df)
TRBJ_name <- unlist(unique(lapply(TRBJ_name,function(x) 
  strsplit(x,"\\-")[[1]][1])))
TRBV_name <- unlist(unique(lapply(TRBV_name,function(x) 
  strsplit(x,"\\.")[[1]][1])))
df2 <- data.frame()
for (TRBV in TRBV_name){
  for (TRBJ in TRBJ_name){
    temp <- df[which(sapply(strsplit(rownames(df),"\\-"),function(x) x[1])==TRBJ),
               which(sapply(strsplit(colnames(df),"\\."),function(x) x[1])==TRBV)]
    df2[TRBJ,TRBV] <- sum(temp)
  }
}
df2 <- cbind(data.frame(.=rownames(df2)),df2)
df2 <- rbind(colnames(df2),df2)
write.table(df2,outfile,row.names = F,sep="\t",col.names = F,quote = F)
}

name <- read.table("./results/VJPAIR/name.txt")
file_list <- name[,"V1"][1:nrow(name)-1]
for (i in 1:length(file_list)){
  v_j_family(paste0(inputdir,file_list[i]),
             paste0(outputdir,file_list[i]))
}


#*********************************plot*********************************#
library(ggplot2)
library(ggsci)
library(gtools)
rm(list=ls())
inputdir <- "./results/VJPAIR/family/"
outputdir <- "./output/VJPAIR/family/"
name <- read.table("./results/VJPAIR/name.txt")
file_list <- name[,"V1"][1:nrow(name)-1]
all_Vgenes <- function(file_list){
  all_V <- c()
  all_J <- c()
  for (i in 1:length(file_list)){
    temp <- read.table(paste0(inputdir,file_list[i]), sep="\t", comment="")
    V <- as.character(temp[1,2:ncol(temp)])
    J <- as.character(temp[2:nrow(temp),1])
    all_V <- append(all_V,V)
    all_J <- append(all_J,J)
  }
  return(all_V)
}
all_v <- unique(all_Vgenes(file_list))
ind_df <- function(file,all_v,group){
  df <- data.frame(vGene=all_v,fancy=rep(0,length(all_v)),group=group)
  rownames(df) <- df$vGene
  temp <- read.csv(file,sep="\t")
  temp <- temp[which(temp$.=="TRBJ2"),]
  rownames(temp) <- temp$.
  temp <- temp[,-1]
  for (j in 1:ncol(temp)){
    df[colnames(temp)[j],"fancy"] <- temp[1,j]
  }
  return(df)
}
all <- data.frame()
for (i in 1:length(file_list)){
  if (i < 5){
    group <- "MTCR1"
  }else{
    group <- "MTCR2"
  }
  df <- ind_df(paste0(inputdir,file_list[i]),all_v,group)
  all <- rbind(all,df)
}
#x order
v <- all$vGene
v_o <- mixedorder(v)
out = all[v_o,]
level <- as.list(unique(out$vGene))
out$vGene <- factor(out$vGene,levels = level)
write.table(out,file = paste0(inputdir,"Vb2-Jb2_GeneUsage.txt"),sep="\t",quote = F,row.names = F)
#plot
p <- ggplot(out,aes(x=vGene,y=fancy))
p+geom_jitter(
  aes(shape = group, color = group), 
  position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
  size = 1.2)+
  stat_summary(
    aes(shape = group, color = group),
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.2,
    position = position_dodge(0.8)
  )+
  theme_classic()+scale_color_aaas()+
  xlab("")+
  ylab("Mean frequency")+labs(title = "V??-J??2 gene usage")+scale_y_continuous(limits = c(0,0.5))
ggsave(paste0(outputdir,"Vb2-Jb2_geneUsage.pdf"),width = 15)



#####CDR3 calculate#####
library(dplyr)
library(plyr)
rm(list=ls())
metadata <- read.table("./metadata.txt",header = 1)
n=3 #number of goups
all_cdr3aa <- data.frame()
for (i in 1:nrow(metadata)){
  df <- read.csv(metadata[i,"file_name"],sep="\t")
  #df <- df[which(df$v=="TRBV13-1*01" & df$j=="TRBJ2-7*01"),]
  unq_cdr3 <- unique(df$cdr3aa)
  a <- unlist(lapply(unq_cdr3, function(x)
    sum(df$count[which(df$cdr3aa==x)])
  ))
  temp <- data.frame(cdr3aa=unq_cdr3,count=a,group=metadata[i,"sample_type"],sample=metadata[i,"sample_id"])
  all_cdr3aa <- rbind(all_cdr3aa,temp)
}
group_cdr3 <- all_cdr3aa %>% group_by(cdr3aa) %>% dplyr::summarise(count = n())
group_cdr3$n <- unlist(lapply(group_cdr3$cdr3aa,function(x) length(unique(all_cdr3aa[which(all_cdr3aa$cdr3aa==x),]$group))))

#------------------specific cdr3---------------------------#
delay <- group_cdr3[which(group_cdr3$n == 1),]
delay <- delay[order(delay$count,decreasing = T),]
cdr3 <- delay$cdr3aa
for (i in 1:n){
  df <- all_cdr3aa[which(all_cdr3aa$cdr3aa %in% cdr3),]
  type <- unique(metadata$sample_type)[i]
  df <- df[which(df$group==type),]
  temp <- aggregate(count ~ cdr3aa ,data=df,sum)
  temp <- temp[order(temp$count,decreasing = T),]
  res <- df[unlist(lapply(temp$cdr3aa,function(x) which(df$cdr3aa==x))),]
  write.table(res,file = paste0("./results/specific_",type,"_CDR3.txt"),sep = "\t",quote = F,row.names = F)
}
#------------------shared cdr3---------------------------#
delay <- group_cdr3[which(group_cdr3$n > 1),]
for (i in n:2){
  cdr3 <- delay[which(delay$n==i),]$cdr3aa
  shared_cdr3aa <- all_cdr3aa[unlist(lapply(cdr3,function(x) which(all_cdr3aa$cdr3aa==x))),]
  temp <- try(aggregate(count ~ cdr3aa ,data=shared_cdr3aa,sum),silent=T)
  if ('try-error' %in% class(temp)){
    next
  }
  temp <- temp[order(temp$count,decreasing = T),]
  cdr3_aa <- temp$cdr3aa
  temp$group <- unlist(lapply(cdr3_aa,function(x) paste(all_cdr3aa[which(all_cdr3aa$cdr3aa==x),]$group,collapse = ",")))
  temp$sn <- unlist(lapply(cdr3_aa,function(x) paste(all_cdr3aa[which(all_cdr3aa$cdr3aa==x),]$count,collapse = ",")))
  temp$ns <- unlist(lapply(cdr3_aa,function(x) paste(all_cdr3aa[which(all_cdr3aa$cdr3aa==x),]$sample,collapse = ",")))
  write.table(temp,file = paste0("./results/",i," Groups_shared_cdr3aa.txt"),sep = "\t",quote = F,row.names = F)
}



#shared cdr3 analysis#
df <- read.csv("E:/TCR/spatial2/bluk_immu/results/bluk_2 Groups_shared_cdr3aa.txt",sep="\t")
 
#------------------ĳһ??CDR3??ÿ?????еĺ?��---------------------------#
A5 <- read.csv("./results/6 Times_shared_A5_cdr3aa.txt",sep = "\t")
cdr3 <- unique(A5$cdr3)
cdr3Content <- function(cdr3,df){
  sum <- sum(df$count)
  df <- df[which(df$cdr3aa==cdr3[1]),]
  count <- sum(df$count)
  ratio <- (count/sum)*100
  return(ratio)
}



res <- data.frame(row.names = metadata$sample_id)

for (i in 1:nrow(metadata)){
  df <- read.csv(metadata[i,"file_name"],sep="\t")
  sample <- metadata$sample_id[i]
  ratio <- c()
  for (k in 1:length(cdr3)){
    res[i,k] <- cdr3Content(cdr3[k],df)
  }
}
colnames(res) <- cdr3
res <- cbind(rownames(res),res)
write.table(res,"./results/specificA5_in_others.txt",sep = "\t",row.names = F,col.names = T,quote = F)
#####PI calculate#####
library(ggplot2)
library(ggsci)
library(ggpubr)
aa_character <- read.table("../IMGT/aa_feature.txt",sep = "\t",header = 1)
name <- read.table("./metadata.txt",sep = "\t",header = 1)
name <- name[which(name$tcr_type=="TRA"|name$tcr_type=="TRB"),]
pI <- function(df){
  cdr3_pI <- function(cdr3){
    sum <- 0
    n <- 0
    for (i in 1:nchar(cdr3)){
      if(substr(cdr3,i,i)=="_"|substr(cdr3,i,i)=="?"){
        sum <- sum
      }else{
        sum <- sum+aa_character[which(aa_character$Abbreviations==substr(cdr3,i,i)),"pI"]
        n <- n+1
      }
      
    }
    pI <- sum/n
    return(pI)
  }
  for (n in 1:nrow(df)){
    cdr3 <- df[n,"cdr3aa"]
    df[n,"pI"] <- cdr3_pI(cdr3)
  }
  return(df)
}
pool <- c()
for (i in 1:nrow(name)){
  df <- read.csv(name$file_name[i],sep = "\t")
  temp <- pI(df)
  res <- data.frame(pI=rep(temp$pI,temp$count),sample=name$sample_id[i],Group=name$sample_type[i],tcr_type=name$tcr_type[i])
  pool <- rbind(pool,res)
}
pool$sample <- unlist(lapply(pool$sample, function(x)
  strsplit(x,"_")[[1]][1]))
write.table(pool,file = "./0913/results/pI_results.txt",sep="\t",row.names = F,quote = F)


pool <- read.csv("./0913/results/pI_results.txt",sep = "\t")
group <- "CON"
work_df <- pool[which(pool$Group==group),]
p <- ggplot(work_df,aes(x=sample,y=pI))
p+geom_jitter(
  aes(shape = tcr_type, color = tcr_type), 
  position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
  size = 1.2)+
  stat_summary(
    aes(shape = tcr_type, color = tcr_type),
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.2,
    position = position_dodge(0.8)
  )+
  theme_classic()+scale_color_aaas()+
  xlab("")+
  ylab("pI")+labs(group)+
  #geom_boxplot(aes(color=tcr_type))+
  stat_compare_means(aes(group = tcr_type),
                     label.x = 1.5,
                    # label.y = max(work_diversity$value)+min(work_diversity$value),
                     label = "p.signif",
                     method = 'wilcox.test',hide.ns = TRUE,size=8)


ggsave(paste0(outputdir,group,"_pI.pdf"))




######clonoVSnon-clono#####
library(ggplot2)
control <- read.csv("./convert_results/control.txt",sep = "\t")
middle <- read.csv("./convert_results/middle.txt",sep = "\t")
severe <- read.csv("./convert_results/severe.txt",sep = "\t")
control_cdr3 <- data.frame(table(control$count),condition="control")
middle_cdr3 <- data.frame(table(middle$count),condition="middle")
severe_cdr3 <- data.frame(table(severe$count),condition="severe")
#---------------overlap------------------------------------#
df1 <- read.csv("./results/2 Groups_shared_cdr3aa.txt",sep = "\t") 
df1 <- data.frame(table(df1$count),condition="sp_2groups")
df2 <- read.csv("../bluk_immu/results/bluk_2 Groups_shared_cdr3aa.txt",sep = "\t")
df2 <- data.frame(table(df2$count),condition="bluk_2groups")
df3 <- read.csv("../bluk_immu/results/bluk_3 Groups_shared_cdr3aa.txt",sep = "\t")
df3 <- data.frame(table(df3$count),condition="bluk_3groups")
df4 <- read.csv("../bluk_immu/results/sp_bluk_mid_2 Groups_shared_cdr3aa.txt",sep = "\t")
df4 <- data.frame(table(df4$count),condition="mid_sp_bluk")
df5 <- read.csv("../bluk_immu/results/sp_blku_sev_2 Groups_shared_cdr3aa.txt",sep = "\t")
df5 <- data.frame(table(df5$count),condition="sev_sp_bluk")
df <- rbind(df1,df2,df3,df4,df5)
#-----------------------------------------------------------#




df <- rbind(control_cdr3,middle_cdr3,severe_cdr3)
df[1] <- lapply(df[1],as.character)
df[1] <- lapply(df[1],as.numeric)
df[,1:2] <- log2(df[,1:2])
colnames(df) <- c("counts","num","condition")
ggplot(data = df, mapping = aes(x = counts, y = num, group = condition,color=condition))+
  geom_point()+
  geom_smooth(method = loess,se = F)+
  labs(title = "CDR3_clonal")+
  theme(axis.text.x = element_text(size = 14, color = "black"))+##????x????????С
  theme(axis.text.y = element_text(size = 14, color = "black"))+##????y????????С
  theme(title=element_text(size=13))+#???ñ?????????С
  theme_bw()+xlab("log2(number of UMI per clonotype)")+ylab("log2(number of clonotypes)")

#####Gliph2#####
rm(list=ls())
c <- read.csv("./control.txt",sep = "\t")
m <- read.csv("./middle.txt",sep = "\t")
s <- read.csv("./severe.txt",sep = "\t")
data <- data.frame(CDR3baa=c(c$cdr3aa,m$cdr3aa,s$cdr3aa),
                   TRBV=c(c$v,m$v,s$v),
                   TRBJ=c(c$j,m$j,s$j),
                   CDR3aaa=NA,
                   subject=c(rep("1:Control",379),rep("2:Middle",6160),rep("3:Severe",6369)),
                   count=c(c$count,m$count,s$count))
write.table(data,file = "./gliph2_test.txt",sep = "\t",row.names = F,col.names = F,quote = F)
