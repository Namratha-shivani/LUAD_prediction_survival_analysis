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
