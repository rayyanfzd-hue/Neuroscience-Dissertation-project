#P14

salmonQuant <- importIsoformExpression(
  parentDir = "P_14 data"
)

metaData2 <- read.table("P_14 data/samples_Cyfip1.txt", sep = "\t", header = T)
metaP14 <- metaData[metaData$age == "juvenile",]

myDesign_P14 <- data.frame(
  sampleID = metaP14$SAMPLE_ID,
  condition = metaP14$genotype
)

colnames(salmonQuant$counts)
myDesign_P14$sampleID
setdiff(myDesign_P14$sampleID, colnames(salmonQuant$counts))

myDesign_P14$sampleID <- paste0(myDesign_P14$sampleID, "_salmon_transcripts")

aSwitchList_P14 <- importRdata(
  isoformCountMatrix = salmonQuant$counts,
  isoformRepExpression = salmonQuant$abundance,
  designMatrix = myDesign_P14,
  isoformExonAnnoation = "/Users/rayyanahamad/Dissertation try 2/gencode.vM37.primary_assembly.annotation.gtf",
  isoformNtFasta = "/Users/rayyanahamad/Dissertation try 2/gencode.vM37.transcripts.fa",
  estimateDifferentialGeneRange = FALSE,
  showProgress = FALSE
)

summary(aSwitchList_P14)

aSwitchListFiltered_P14 <- preFilter(
  switchAnalyzeRlist = aSwitchList_P14,
  isoCount = 10,
  min.Count.prop = 0.7,
  IFcutoff = 0.1,
  min.IF.prop = 0.5,
  removeSingleIsoformGenes = TRUE
)

summary(aSwitchListFiltered_P14)

#With reduceToSwitchingGenes = TRUE results came back with 0, so retried with it FALSE

aSwitchListAnalyzed_P14 <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = aSwitchListFiltered_P14,
  reduceToSwitchingGenes = TRUE
)

extractSwitchSummary(aSwitchListAnalyzed_P14)

#Try SatuRN with dIFcutoff = 0.025 for P14 also, because null

aSwitchListAnalyzed_P14_0.025 <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = aSwitchListFiltered_P14,
  dIFcutoff = 0.025,
  reduceToSwitchingGenes = FALSE
)

extractSwitchSummary(aSwitchListAnalyzed_P14_0.025)

