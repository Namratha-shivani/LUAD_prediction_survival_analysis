# import libraries

library(tidyverse)
library(TCGAbiolinks)
library(GEOquery)
library(DESeq2)
library(SummarizedExperiment)
library(dplyr)
library(WGCNA)
library(ggplot2)
library(gridExtra)
library(flashClust)
library(survival)
library(survminer)
library(gProfileR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(KEGGREST)
library(KEGGgraph)
library(edgeR)

#get clinical data from TCGA 
clinical_TCGA = GDCquery_clinic('TCGA-LUAD')


#checking the number of cases Alive and Dead
table(clinical_TCGA$vital_status)
clinical_TCGA$deceased <- ifelse(clinical_TCGA$vital_status=='Alive', 1, 0)
clinical_TCGA$overallSurvival <- ifelse(clinical_TCGA$vital_status=='Alive',
                                        clinical_TCGA$days_to_last_follow_up,
                                        clinical_TCGA$days_to_death)


#get expression data

LUAD_all <- GDCquery(
  project = 'TCGA-LUAD',
  data.category = 'Transcriptome Profiling',
  experimental.strategy = 'RNA-Seq',
  workflow.type = 'STAR - Counts',
  data.type = 'Gene Expression Quantification',
  sample.type = c('Primary Tumor', 'Solid Tissue Normal'),
  access = 'open'
)

output_LUAD = getResults(LUAD_all)

GDCdownload(LUAD_all)

tcga_luad_data <- GDCprepare(LUAD_all, summarizedExperiment = TRUE)
luad_matrix <- assay(tcga_luad_data,'unstranded')

# get the column and row metadata 
gene_metadata <- as.data.frame(rowData(tcga_luad_data))
coldata <- as.data.frame(colData(tcga_luad_data))

#add overall survival data to the coldata and subset it to include only the wanted columns
coldata$deceased <- ifelse(coldata$vital_status =='Alive',1,0)
coldata$overallSurvival <- ifelse(coldata$vital_status == 'Alive',
                                  coldata$days_to_last_follow_up,
                                  coldata$days_to_death)

TCGA_clinical <- coldata[,c(5,24,51:52,55,88:89)]
TCGA_clinical <- TCGA_clinical[which(! is.na(TCGA_clinical$age_at_index)),]


# fill in the missing stage information based on the TMN staging

Overall_stage = c('Stage 0', 'Stage I', 'Stage IIA', 'Stage IIA', 'Stage IIA',
                  'Stage IIB','Stage IIB','Stage IIIA','Stage IIIA','Stage IIIA',
                  'Stage IIIA','Stage IIIA','Stage IIIB','Stage IIIC','Stage IV') 
T_category = c('Tis', 'T1','T0','T1','T2','T2','T3','T0','T1','T2',
               'T3','T3','T4','Any T','Any T')
N_category = c('N0','N0','N1','N1','N0','N1','N0','N2','N2',
               'N2','N1','N2','Any N','N3','Any N')
M_category = c('M0','M0','M0','M0','M0','M0','M0','M0','M0','M0',
               'M0','M0','M0','M0','M1')

cancer_stage <- data.frame(Overall_stage,T_category, N_category, M_category)

subset <- rownames(TCGA_clinical)[which(is.na(TCGA_clinical$ajcc_pathologic_stage))]
data <- coldata[subset,][,c(36,38:39)]
data <- data[which(!is.na(data$ajcc_pathologic_m)),]

x <- c()
for (i in rownames(data)){
  for (j in c(1:nrow(cancer_stage))){
    if (substring(data[i,1],1,2) == cancer_stage[j,]$T_category){
      if (cancer_stage[j,]$N_category == data[i,2]){
        if (cancer_stage[j,]$M_category == data[i,3]){
          TCGA_clinical[i,]$ajcc_pathologic_stage <- cancer_stage$Overall_stage[j]
          x = c(x,cancer_stage$Overall_stage[j])
        }else{
          if (cancer_stage[j,]$T_category == substring(data[i,1], 1,2) & cancer_stage[j,]$N_category == 'Any N'){
            TCGA_clinical[i,]$ajcc_pathologic_stage <- cancer_stage$Overall_stage[j]
            x = c(x,cancer_stage$Overall_stage[j]) 
          }
        }
      }
      
    }
  }
}

# filter out the clinical data to remove out the na's 
# in the matrix remove the filtered out files
TCGA_clinical <- TCGA_clinical[which(! is.na(TCGA_clinical$ajcc_pathologic_stage)),]
luad_matrix <- luad_matrix[,which(rownames(TCGA_clinical) %in% colnames(luad_matrix))]

# Differential Expression
#create a design matrix to model for DESeq2
design <- data.frame(row.names=colnames(luad_matrix))
design$age_at_index <- as.numeric(TCGA_clinical$age_at_index)
design$gender <- factor(TCGA_clinical$gender)
design$race <- factor(TCGA_clinical$race)
design$type_of_sample <- factor(TCGA_clinical$definition)

dds <- DESeqDataSetFromMatrix(countData = luad_matrix,
                              colData   = design,
                              design = ~ age_at_index+gender+race+type_of_sample
)

# keep only those genes with greater than 20 reads expressed
dds75 <- dds[rowSums(counts(dds) >=10) >= 433,]

dds75$type_of_sample <- relevel(dds75$type_of_sample, ref = 'Solid Tissue Normal')

dds <- DESeq(dds75)

# perform variant stabilization tranformation to normalize the RNA-seq data
# it can improve the accuracy and power of differential expression.


vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
luad_matrix_vst <- assay(vsd)

res <- results(dds, alpha = 0.05)

result_de <- subset(res,padj<0.05)
result_de <- data.frame(result_de[order(result_de$padj),])

luad_matrix1 <- counts(dds)

# WGCNA
tumor_samples <- rownames(TCGA_clinical)[which(TCGA_clinical$definition == 'Primary solid Tumor')]
luad_data <- as.data.frame(luad_matrix1) %>% t()
luad_data <- luad_data[which(rownames(luad_data) %in% tumor_samples),]

#SOFT POWER
power <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function

sft <- pickSoftThreshold(luad_data,
                         powerVector = power,
                         verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.8;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=power,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=power, cex=cex1,col="red")


adjacencyluad = adjacency((luad_data),power = 6, type= "unsigned")
#diag(adjacencyluad=0)
dissTOM = 1-TOMsimilarity(adjacencyluad, TOMType = "unsigned")
geneTree = flashClust(as.dist(dissTOM), method = 'average')

plot(geneTree, main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)



mcolorh = NULL
for (ds in 0:3){
  tree = cutreeHybrid(dendro = geneTree, pamStage = FALSE, minClusterSize = 30, cutHeight = 0.99, deepSplit = ds, distM = dissTOM)
  mcolorh = cbind(mcolorh, labels2colors(tree$labels))
}

plotDendroAndColors(geneTree, mcolorh, paste("dsplt =",0:3), main = "",dendroLabels = FALSE)

module = mcolorh[,2]
plotDendroAndColors(geneTree, module, "Modules", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendogram and module colors")


PCs = moduleEigengenes(luad_data2, colors = module)
MEs = PCs$eigengenes
distPC = 1-abs(cor(ME, use = 'p'))
distPC = ifelse(is.na(distPC),0, distPC)
pcTree = hclust(as.dist(distPC),method = 'a')
MDS = cmdscale(as.dist(distPC),2)
colors = names(table(module))



MM = cor(luad_data2,MEs,method = 'pearson')
geneModuleMembership = signedKME(luad_data2, MEs)
colnames(geneModuleMembership)=paste("PC",colors,".cor",sep="");
MMPvalue=corPvalueStudent(as.matrix(geneModuleMembership),dim(luad_data2)[[2]]);
colnames(MMPvalue)=c("blue","brown","grey","turquoise");
Gene = rownames(t(luad_data2))
kMEtable = cbind(Gene,Gene,module)
for (i in 1:length(colors))
  kMEtable = cbind(kMEtable, geneModuleMembership[,i], MMPvalue[,i]) 
colnames(kMEtable)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership), 
                                                   colnames(MMPvalue))))

kMEtable <- data.frame(kMEtable)
blue_genes_module <-kMEtable[which(kMEtable$Module == 'blue'),]
brown_genes_module <- kMEtable[which(kMEtable$Module == 'brown'),]
turquoise_genes_module <- kMEtable[which(kMEtable$Module == 'turquoise'),]
grey_genes_module <-kMEtable[which(kMEtable$Module == 'grey'),]

survival_data <- survival_data[which(row.names(survival_data) %in% rownames(luad_data)),]
survival_data <- survival_data[,c(6,7)]


blue_module_expression <- luad_data2[,which(colnames(luad_data2) %in% blue_genes_module$Gene)]
r_blue <- rowMeans(blue_module_expression)
lung_cox <- coxph(Surv(overallSurvival, deceased) ~ ., data = cbind(r_blue, survival_data))
lung_risk <- predict(lung_cox, type = "risk")
median_risk <- median(lung_risk)
high_risk <- lung_risk >= median_risk
low_risk <- lung_risk < median_risk

ggsurvplot(survfit(Surv(survival_data$overallSurvival, survival_data$deceased) ~ high_risk),
           data = survival_data,
           risk.table = TRUE,
           conf.int = TRUE,
           pval = TRUE,
           legend.title = "Risk group",
           legend.labs = c("Low risk", "High risk"),
           xlab = "Time (days)",
           ylab = "Survival probability",
           title = "Blue Module Kaplan-Meier curves for high and low risk groups in LUAD")



brown_module_expression <- luad_data2[,which(colnames(luad_data2) %in% brown_genes_module$Gene)]

r <- rowMeans(brown_module_expression)
lung_cox <- coxph(Surv(overallSurvival, deceased) ~ ., data = cbind(r, survival_data))
lung_risk_brown <- predict(lung_cox, type = "risk")
median_risk_brown <- median(lung_risk_brown)
high_risk_brown <- lung_risk_brown >= median_risk_brown
low_risk_brown <- lung_risk_brown < median_risk_brown

ggsurvplot(survfit(Surv(survival_data$overallSurvival, survival_data$deceased) ~ high_risk_brown),
           data = survival_data,
           risk.table = TRUE,
           conf.int = TRUE,
           pval = TRUE,
           legend.title = "Risk group",
           legend.labs = c("Low risk", "High risk"),
           xlab = "Time (days)",
           ylab = "Survival probability",
           title = "Brown Module Kaplan-Meier curves for high low risk groups in LUAD")


turquoise_module_expression <- luad_data[,which(colnames(luad_data) %in% turquoise_genes_module$Gene)]

r_turquoise <- rowMeans(turquoise_module_expression)
lung_cox <- coxph(Surv(overallSurvival, deceased) ~ ., data = cbind(r_turquoise, survival_data))
lung_risk_turquoise <- predict(lung_cox, type = "risk")
median_risk_turquoise <- median(lung_risk_turquoise)
high_risk_turquoise <- lung_risk_turquoise >= median_risk_turquoise
low_risk_turquoise <- lung_risk_turquoise < median_risk_turquoise

ggsurvplot(survfit(Surv(survival_data$overallSurvival, survival_data$deceased) ~ high_risk_turquoise),
           data = survival_data,
           risk.table = TRUE,
           conf.int = TRUE,
           pval = TRUE,
           legend.title = "Risk group",
           legend.labs = c("Low risk", "High risk"),
           xlab = "Time (days)",
           ylab = "Survival probability",
           title = "Turquoise Module Kaplan-Meier curves for high and low risk groups in LUAD")


#clustering samples to see the expression difference

d <- dist(t(luad_data), method="euclidean")
hc <- hclust(d)

# Use the elbow method to determine the number of clusters
wss <- (nrow(luad_data)-1)*sum(apply(luad_data,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(luad_data, centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of clusters", ylab="Within groups sum of squares")

plot(hc, main = "Dendrogram of Hierarchical Clustering", xlab = "Sample Index", sub = NULL)

# Use k-means clustering to assign samples to clusters
num_clusters <- 5
clusters <- kmeans(luad_data, centers=num_clusters)$cluster


group <- clusters ==3
survival_clu3 <- survival_data[which(rownames(survival_data) %in% gsub("\\.", "-", colnames(DEGs_cluster3))),]
DEGs_cluster3 <- data.frame(DEGs_cluster3[,which(colnames(DEGs_cluster3) %in% rownames((survival_clu3)))])
survival_clu3$strata <- ifelse(means_col3 >= median_clu3, 'HIGH', 'LOW')
km <- survfit(Surv(overallSurvival, deceased) ~ strata, data = survival_clu3)
plot(km, main=paste("Cluster", 3), xlab="Time (days)", ylab="Survival probability")
ggsurvplot(km, data = survival_clu3, 
           pval = TRUE, conf.int = TRUE, 
           risk.table = TRUE, risk.table.col = "strata",
           xlim = c(0, 8000), ylim = c(0, 1), 
           legend.title = "Group", legend.labs =  c("Low Expression", "High Expression"),
           ggtheme = theme_minimal(),
           palette = c("#0072B2", "#D55E00"),
           main = paste("Cluster", 3),
           xlab = "Time (days)", ylab = "Survival probability")


group <- clusters ==1
DEGs_cluster1 <- data.frame(DEGs[,which(group == TRUE)])
colnames(DEGs_cluster1) <- gsub("\\.", "-", colnames(DEGs_cluster1))
survival_clu1 <- survival_data[which(rownames(survival_data) %in% gsub("\\.", "-", colnames(DEGs_cluster1))),]
DEGs_cluster1 <- data.frame(DEGs_cluster1[,which(colnames(DEGs_cluster1) %in% rownames((survival_clu1)))])
means_col1 <- colMeans(DEGs_cluster1) 
median_clu1 <- median(means_col1) 
survival_clu1$strata <- ifelse(means_col1 >= median_clu1, 'HIGH', 'LOW')
km <- survfit(Surv(overallSurvival, deceased) ~ strata, data = survival_clu1)
ggsurvplot(km, data = survival_clu1, 
           pval = TRUE, conf.int = TRUE, 
           risk.table = TRUE, risk.table.col = "strata",
           xlim = c(0, 8000), ylim = c(0, 1), 
           legend.title = "Group", legend.labs =  c("Low Expression", "High Expression"),
           ggtheme = theme_minimal(),
           palette = c("#0072B2", "#D55E00"),
           main = paste("Cluster", 1),
           xlab = "Time (days)", ylab = "Survival probability")


group <- clusters ==2
DEGs_cluster2 <- data.frame(DEGs[,which(group == TRUE)])
colnames(DEGs_cluster2) <- gsub("\\.", "-", colnames(DEGs_cluster2))
survival_clu2 <- survival_data[which(rownames(survival_data) %in% gsub("\\.", "-", colnames(DEGs_cluster2))),]
DEGs_cluster2 <- data.frame(DEGs_cluster2[,which(colnames(DEGs_cluster2) %in% rownames((survival_clu2)))])
means_col2 <- colMeans(DEGs_cluster2) 
median_clu2 <- median(means_col2) 
survival_clu2$strata <- ifelse(means_col2 >= median_clu2, 'HIGH', 'LOW')
km <- survfit(Surv(overallSurvival, deceased) ~ strata, data = survival_clu2)
ggsurvplot(km, data = survival_clu2, 
           pval = TRUE, conf.int = TRUE, 
           risk.table = TRUE, risk.table.col = "strata",
           xlim = c(0, 8000), ylim = c(0, 1), 
           legend.title = "Group", legend.labs =  c("Low Expression", "High Expression"),
           ggtheme = theme_minimal(),
           palette = c("#0072B2", "#D55E00"),
           main = paste("Cluster", 2),
           xlab = "Time (days)", ylab = "Survival probability")


group <- clusters ==4
DEGs_cluster4 <- data.frame(DEGs[,which(group == TRUE)])
colnames(DEGs_cluster4) <- gsub("\\.", "-", colnames(DEGs_cluster4))
survival_clu4 <- survival_data[which(rownames(survival_data) %in% gsub("\\.", "-", colnames(DEGs_cluster4))),]
DEGs_cluster4 <- data.frame(DEGs_cluster4[,which(colnames(DEGs_cluster4) %in% rownames((survival_clu4)))])
means_col4 <- colMeans(DEGs_cluster4) 
median_clu4 <- median(means_col4) 
survival_clu4$strata <- ifelse(means_col4 >= median_clu4, 'HIGH', 'LOW')
km <- survfit(Surv(overallSurvival, deceased) ~ strata, data = survival_clu4)
ggsurvplot(km, data = survival_clu4, 
           pval = TRUE, conf.int = TRUE, 
           risk.table = TRUE, risk.table.col = "strata",
           xlim = c(0, 8000), ylim = c(0, 1), 
           legend.title = "Group", legend.labs =  c("Low Expression", "High Expression"),
           ggtheme = theme_minimal(),
           palette = c("#0072B2", "#D55E00"),
           main = paste("Cluster", 4),
           xlab = "Time (days)", ylab = "Survival probability")


group <- clusters ==5
DEGs_cluster5 <- data.frame(DEGs[,which(group == TRUE)])
colnames(DEGs_cluster5) <- gsub("\\.", "-", colnames(DEGs_cluster5))
survival_clu5 <- survival_data[which(rownames(survival_data) %in% gsub("\\.", "-", colnames(DEGs_cluster5))),]
DEGs_cluster5 <- data.frame(DEGs_cluster5[,which(colnames(DEGs_cluster5) %in% rownames((survival_clu5)))])
means_col5 <- colMeans(DEGs_cluster5) 
median_clu5 <- median(means_col5) 
survival_clu5$strata <- ifelse(means_col5 >= median_clu5, 'HIGH', 'LOW')
km <- survfit(Surv(overallSurvival, deceased) ~ strata, data = survival_clu5)
ggsurvplot(km, data = survival_clu5, 
           pval = TRUE, conf.int = TRUE, 
           risk.table = TRUE, risk.table.col = "strata",
           xlim = c(0, 2000), ylim = c(0, 1), 
           legend.title = "Group", legend.labs =  c("Low Expression", "High Expression"),
           ggtheme = theme_minimal(),
           palette = c("#0072B2", "#D55E00"),
           main = paste("Cluster", 5),
           xlab = "Time (days)", ylab = "Survival probability")


# GO 

genes <- gene_metadata[which(gene_metadata$gene_id %in% rownames(result_de)), ]$gene_name
result <- enrichGO(gene     = genes, 
                   OrgDb    = org.Hs.eg.db, 
                   keyType  = "SYMBOL",
                   ont      = "BP", 
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.1)

dotplot(result)




























