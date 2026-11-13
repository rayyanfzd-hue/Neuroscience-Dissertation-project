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




