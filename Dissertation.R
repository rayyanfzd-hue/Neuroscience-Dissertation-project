# Install packages
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install()
}
BiocManager::install("IsoformSwitchAnalyzeR")
install.packages("dplyr")
BiocManager::install("edgeR")

# Load packages
library(IsoformSwitchAnalyzeR)
library(dplyr)
library(edgeR)

#Misc
getwd()
file.exists("/Users/rayyanahamad/Dissertation Project/gencode.vM37.primary_assembly.annotation.gtf")

setwd("~/Dissertation Project")

#Github
install.packages("usethis")
usethis::use_git_config(user.name = "rayyanfzd-hue", user.email = "rayyanfzd@gmail.com")

usethis::create_github_token()

gitcreds::gitcreds_set()

usethis::use_git_ignore(c("Adult_data", "P14_data", "*.gtf", "*.fa"))

usethis::use_git()

#Importing into R

salmonQuant <- importIsoformExpression(
  parentDir = "Adult_data")

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
  isoformExonAnnoation = "/Users/rayyanahamad/Dissertation Project/gencode.vM37.primary_assembly.annotation.gtf",
  isoformNtFasta = "/Users/rayyanahamad/Dissertation Project/gencode.vM37.transcripts.fa",
  showProgress = FALSE
)

# Filtering

data("aSwitchList")

aSwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  isoCount = 10,
  min.Count.prop = 0.7,
  IFcutoff = 0.1,
  min.IF.prop = 0.5,
  removeSingleIsoformGenes = TRUE
)

#Identifying isoform switches

