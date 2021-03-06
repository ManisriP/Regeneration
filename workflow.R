options(stringsAsFactors = F)
library(xlsx)
library(openxlsx)
library(dplyr)
#library(IRanges)
library(matrixStats)
#library(RColorBrewer)
#library(pheatmap)
#library(EnrichR)

dd=read.csv(file="GSE95135_Rib_et_al.RPKM_log2.csv",header = T,sep="\t")

dd[,-c(1:3)]=log2(2^dd[,-c(1:3)] + 1)

exp=dd[-c(10:17,29:32)]  
mydata1 <- exp[,-c(2,3)]
rownames(mydata1) <- mydata1$geneid
mydata1$geneid <- NULL
mydata1$CV <- rowVars(as.matrix(mydata1))
mydata1 <- mydata1[order(mydata1$CV, decreasing = TRUE),]

### DREM input for PH samples(top 5000 varying genes)
mydata <- mydata1[1:5000,c(1:6,27:72)] ### Extracting 0hr and post-PH samples
# colnames(mydata)=substr(colnames(mydata),1,  nchar(colnames(mydata))-1)


s=unique(substr(colnames(mydata),1,  nchar(colnames(mydata))-1))
drem_inp_5000=matrix(,ncol=length(s)-1,nrow=0)
for (j in 1:nrow(mydata)){
t=lapply(sapply(s, grep, x=colnames(mydata)), 
function(m,n)  mean(as.numeric(m[n])),  m=mydata[j,])
drem_inp_5000=rbind(drem_inp_5000,(unlist(t)-unlist(t[1]))[-1])
}

drem_inp_5000=cbind(exp$GeneName[match(rownames(mydata), exp$geneid)],drem_inp_5000)
colnames(drem_inp_5000)[1]="Probe"

write.table(drem_inp_5000, file="gse95135_phdrem_5000.txt", row.names=F, quote=F, sep="\t")  

#############################################################################################################
#############################################################################################################

############ WGCNA 
options(java.parameters = "-Xmx8000m")
library(WGCNA)                      
library(flashClust)
library(matrixStats)
library(knitr)
library(dplyr)
library(openxlsx)
library(stringi)



options(stringsAsFactors = FALSE);
enableWGCNAThreads()

traitData = read.csv("gse95135_traits.tsv", sep = "\t", stringsAsFactors = FALSE);
exp <- read.csv("gse95135_expression_wgcna.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)[-c(10:17,29:32)]
mydata1 <- exp[,-c(2,3)]
rownames(mydata1) <- mydata1$geneid
mydata1$geneid <- NULL
mydata1$CV <- rowVars(as.matrix(mydata1))
mydata1 <- mydata1[order(mydata1$CV, decreasing = TRUE),]



mydata <- mydata1[1:5000,c(1:6,27:72)]
n=nrow(mydata);
gene.names = rownames(mydata)
mydata.trans = t(mydata)
datExpr=mydata.trans[,gene.names[1:n]];
SubGeneNames=gene.names[1:n];

#Checking data for missing values
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}
dim(datExpr)
n=ncol(datExpr);
SubGeneNames=colnames(datExpr);

powers = c(c(1:9), seq(from = 10, to=30, by=1))
sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,powerVector = powers,corFnc = "cor", corOptions=list(use ='p'), networkType = "signed")
cex1 = 1.0;
mm=cbind(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],sft$fitIndices[,5])
colnames(mm)=c("Power","R-sqd","Mean connectivity")


softPower = 15;
####
adj1 = adjacency(datExpr, power = softPower,type="signed")
TOM = TOMsimilarityFromExpr(datExpr, networkType = "signed",
                            corType = "pearson",
                            nThreads = 5, power = softPower)
colnames(TOM) =rownames(TOM) = SubGeneNames
dissTOM = 1-TOM


geneTree = flashClust(as.dist(dissTOM), method = "average");

minModuleSize = 100;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method = "hybrid",
                            deepSplit = 2, pamStage = TRUE, pamRespectsDendro = TRUE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)



MEList = moduleEigengenes(datExpr, colors = dynamicColors, softPower =softPower)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = flashClust(as.dist(MEDiss), method = "average");
MEDissThres = 0.15
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
unique(moduleColors)
table(moduleColors)

module_colors = unique(moduleColors)
merged_modules=vector("list",length(module_colors))
names(merged_modules)=module_colors

for (color in module_colors){
  module=SubGeneNames[which(moduleColors==color)]
  merged_modules[[color]] = exp %>% subset(geneid %in% module) %>% dplyr::select(geneid,GeneName)
  }


dim(traitData)
names(traitData)
allTraits = traitData
dim(allTraits)
names(allTraits)
samples = rownames(datExpr);
traitRows = match(samples, allTraits$Sample);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
MEs0 = moduleEigengenes(datExpr, moduleColors, softPower = softPower)$eigengenes
MEs = orderMEs(MEs0)

moduleTraitCor = cor(MEs, datTraits, use = "p");            #### Module Eigen gene-trait relationship
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
