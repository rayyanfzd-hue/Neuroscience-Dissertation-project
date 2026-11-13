# Install packages
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install()
}
BiocManager::install("IsoformSwitchAnalyzeR")
install.packages("dplyr")
BiocManager::install("edgeR")

#Load packages
library(IsoformSwitchAnalyzeR)
library(dplyr)
library(edgeR)

#Saves
dir.create("results", showWarnings = FALSE)

saveRDS(salmonQuant,          "results/adult_salmonQuant.rds")
saveRDS(aSwitchList,          "results/adult_aSwitchList.rds")
saveRDS(aSwitchListFiltered,  "results/adult_aSwitchListFiltered.rds")
saveRDS(aSwitchListAnalyzed2,  "results/adult_aSwitchListAnalyzed.rds")
saveRDS(aSwitchAnalyzed_0.025,  "results/adult_aSwitchAnalyzed_0.025.rds")
saveRDS(aSwitchList_Adult_DG,  "results/adult_aSwitchList_Adult_DG.rds")
saveRDS(aSwitchListFiltered_Adult_DG,  "results/adult_aSwitchListFiltered_Adult_DG.rds")
saveRDS(aSwitchListAnalyzed_Adult_DG,  "results/adult_aSwitchListAnalyzed_Adult_DG.rds")
saveRDS(aSwitchListAnalyzed_Adult_DG_0.025, "results/adult_aSwitchListAnalyzed_Adult_DG_0.025")
saveRDS(aSwitchList_Adult_PFC, "results/adult_aSwitchList_Adult_PFC")
saveRDS(aSwitchListFiltered_Adult_PFC, "results/adult_aSwitchListFiltered_Adult_PFC")


aSwitchList_Adult_DG <- readRDS("results/adult_aSwitchList_Adult_DG.rds")

#Misc
getwd()
setwd("~/Dissertation try 2")

#ext
res <- aSwitchListAnalyzed2$isoformSwitchAnalysis
head(res)
sum(res$isoform_switch_q_value < 0.05, na.rm = TRUE)
summary(abs(res$dIF))

res <- aSwitchListAnalyzed2$isoformSwitchAnalysis

colnames(res)
sum(res$pval < 0.05, na.rm = TRUE)
sum(res$padj < 0.05, na.rm = TRUE)
min(res$padj, na.rm = TRUE)

isoFeat <- aSwitchListAnalyzed2$isoformFeatures
sum(isoFeat$isoform_switch_q_value < 0.05, na.rm = TRUE)
summary(abs(isoFeat$dIF))

range(isoFeat$isoform_switch_q_value, na.rm = TRUE)
hist(res$pval, breaks = 40, main = "SatuRn p values, het vs wt")

#GitHub


#Imporint into R
salmonQuant <- importIsoformExpression(
  parentDir = "Adult_data"
)

metaData <- read.table("Adult_data/samples_Cyfip1.txt", sep = "\t", header = T)
metaAdult <- metaData[metaData$age == "adult",]

myDesign <- data.frame(
  sampleID = metaAdult$SAMPLE_ID,
  condition = metaAdult$genotype,
  region = metaAdult$region
)

colnames(salmonQuant$counts)
myDesign$sampleID
setdiff(myDesign$sampleID, colnames(salmonQuant$counts))

myDesign$sampleID <- paste0(myDesign$sampleID, "_salmon_transcripts")
myDesign$sampleID[myDesign$sampleID == "CYFIP1_PFC_12_salmon_transcripts"] <-
  "CYFIP1_PFC_12_PFC_salmon_transcripts"
myDesign$sampleID[myDesign$sampleID == "CYFIP1_PFC_14_salmon_transcripts"] <-
  "CYFIP1_PFC_14_PFC_salmon_transcripts"

aSwitchList <- importRdata(
  isoformCountMatrix = salmonQuant$counts,
  isoformRepExpression = salmonQuant$abundance,
  designMatrix = myDesign,
  isoformExonAnnoation = "/Users/rayyanahamad/Dissertation try 2/gencode.vM37.primary_assembly.annotation.gtf",
  isoformNtFasta = "/Users/rayyanahamad/Dissertation try 2/gencode.vM37.transcripts.fa",
  estimateDifferentialGeneRange = FALSE,
  showProgress = FALSE
)


#Filtering

aSwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  isoCount = 10,
  min.Count.prop = 0.7,
  IFcutoff = 0.1,
  min.IF.prop = 0.5,
  removeSingleIsoformGenes = TRUE
)

#Identifying Isofrom Switches
aSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = aSwitchListFiltered,
  reduceToSwitchingGenes = TRUE
)
#With reduceToSwitchingGenes = TRUE results came back with 0, so retried with it FALSE

aSwitchListAnalyzed2 <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = aSwitchListFiltered,
  reduceToSwitchingGenes = FALSE
)

#Results still = null, so attempted to adjust for dIFcutoff

extractSwitchSummary(aSwitchListAnalyzed2)

aSwitchAnalyzed_0.025 <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = aSwitchListFiltered,
  dIFcutoff = 0.025,
  reduceToSwitchingGenes = FALSE
)

extractSwitchSummary(aSwitchAnalyzed_0.025)

#Still no results, try to analyse brain regions seperatley
#Adult DG

salmonQuant <- importIsoformExpression(
  parentDir = "Adult_data"
)

metaData <- read.table("Adult_data/samples_Cyfip1.txt", sep = "\t", header = T)
metaAdult <- metaData[metaData$age == "adult",]

metaAdult_DG <- metaAdult[grepl(pattern = "^CYFIP1_DG", x = metaAdult$SAMPLE_ID),]
myDesign_DG <- data.frame(
  sampleID = metaAdult_DG$SAMPLE_ID,
  condition = metaAdult_DG$genotype
)

colnames(salmonQuant$counts)
myDesign_DG$sampleID
setdiff(myDesign_DG$sampleID, colnames(salmonQuant$counts))

myDesign_DG$sampleID <- paste0(myDesign_DG$sampleID, "_salmon_transcripts")

aSwitchList_Adult_DG <- importRdata(
  isoformCountMatrix = salmonQuant$counts,
  isoformRepExpression = salmonQuant$abundance,
  designMatrix = myDesign_DG,
  isoformExonAnnoation = "/Users/rayyanahamad/Dissertation try 2/gencode.vM37.primary_assembly.annotation.gtf",
  isoformNtFasta = "/Users/rayyanahamad/Dissertation try 2/gencode.vM37.transcripts.fa",
  estimateDifferentialGeneRange = FALSE,
  showProgress = FALSE
)

aSwitchListFiltered_Adult_DG <- preFilter(
  switchAnalyzeRlist = aSwitchList_Adult_DG,
  isoCount = 10,
  min.Count.prop = 0.7,
  IFcutoff = 0.1,
  min.IF.prop = 0.5,
  removeSingleIsoformGenes = TRUE
)

summary(aSwitchListFiltered_Adult_DG)

#With reduceToSwitchingGenes = TRUE results came back with 0, so retried with it FALSE

aSwitchListAnalyzed_Adult_DG <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = aSwitchListFiltered_Adult_DG,
  reduceToSwitchingGenes = FALSE
)

extractSwitchSummary(aSwitchListAnalyzed_Adult_DG)

#Try SatuRN with dIFcutoff = 0.025 for Adult DG also, because null

aSwitchListAnalyzed_Adult_DG_0.025 <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = aSwitchListFiltered_Adult_DG,
  dIFcutoff = 0.025,
  reduceToSwitchingGenes = FALSE
)

extractSwitchSummary(aSwitchListAnalyzed_Adult_DG_0.025)


#Adult PFC

salmonQuant <- importIsoformExpression(
  parentDir = "Adult_data"
)

metaData <- read.table("Adult_data/samples_Cyfip1.txt", sep = "\t", header = T)
metaAdult <- metaData[metaData$age == "adult",]

metaAdult_PFC <- metaAdult[grepl(pattern = "^CYFIP1_PFC", x = metaAdult$SAMPLE_ID),]
myDesign_PFC <- data.frame(
  sampleID = metaAdult_PFC$SAMPLE_ID,
  condition = metaAdult_PFC$genotype
)

colnames(salmonQuant$counts)
myDesign_PFC$sampleID
setdiff(myDesign_PFC$sampleID, colnames(salmonQuant$counts))

myDesign_PFC$sampleID <- paste0(myDesign_PFC$sampleID, "_salmon_transcripts")
myDesign_PFC$sampleID[myDesign_PFC$sampleID == "CYFIP1_PFC_12_salmon_transcripts"] <-
  "CYFIP1_PFC_12_PFC_salmon_transcripts"
myDesign_PFC$sampleID[myDesign_PFC$sampleID == "CYFIP1_PFC_14_salmon_transcripts"] <-
  "CYFIP1_PFC_14_PFC_salmon_transcripts"

aSwitchList_Adult_PFC <- importRdata(
  isoformCountMatrix = salmonQuant$counts,
  isoformRepExpression = salmonQuant$abundance,
  designMatrix = myDesign_PFC,
  isoformExonAnnoation = "/Users/rayyanahamad/Dissertation try 2/gencode.vM37.primary_assembly.annotation.gtf",
  isoformNtFasta = "/Users/rayyanahamad/Dissertation try 2/gencode.vM37.transcripts.fa",
  estimateDifferentialGeneRange = FALSE,
  showProgress = FALSE
)

aSwitchListFiltered_Adult_PFC <- preFilter(
  switchAnalyzeRlist = aSwitchList_Adult_PFC,
  isoCount = 10,
  min.Count.prop = 0.7,
  IFcutoff = 0.1,
  min.IF.prop = 0.5,
  removeSingleIsoformGenes = TRUE
)

summary(aSwitchListFiltered_Adult_PFC)

#With reduceToSwitchingGenes = TRUE results came back with 0, so retried with it FALSE

aSwitchListAnalyzed_Adult_DG <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = aSwitchListFiltered_Adult_DG,
  reduceToSwitchingGenes = FALSE
)

extractSwitchSummary(aSwitchListAnalyzed_Adult_DG)

#Try SatuRN with dIFcutoff = 0.025 for Adult DG also, because null

aSwitchListAnalyzed_Adult_DG_0.025 <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = aSwitchListFiltered_Adult_DG,
  dIFcutoff = 0.025,
  reduceToSwitchingGenes = FALSE
)

extractSwitchSummary(aSwitchListAnalyzed_Adult_DG)







