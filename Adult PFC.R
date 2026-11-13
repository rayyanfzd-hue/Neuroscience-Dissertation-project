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

aSwitchListAnalyzed_Adult_PFC <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = aSwitchListFiltered_Adult_PFC,
  reduceToSwitchingGenes = FALSE
)

extractSwitchSummary(aSwitchListAnalyzed_Adult_PFC)

#Try SatuRN with dIFcutoff = 0.025 for Adult PFC also, because null

aSwitchListAnalyzed_Adult_PFC_0.025 <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = aSwitchListFiltered_Adult_PFC,
  dIFcutoff = 0.025,
  reduceToSwitchingGenes = FALSE
)

extractSwitchSummary(aSwitchListAnalyzed_Adult_PFC_0.025)



