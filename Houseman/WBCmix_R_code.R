## Objects to save
setwd("/users/waccomando/Desktop")
save(file="WBCmixLMEobjects.RData", compress=TRUE, list=c("Betas_WBC","pheno0", "annot"))

### R code for Quantitative Reconstruction of Leukocyte Subsets Using DNA Methylation

# Part 1: Linear Mixed Effects Model for Purified WBC Subset Samples

load(WBCmixLMEobjects.RData)

library(nlme)

theModel = y~PanTcell+NK+Bcell+Mono+Gran
sizeModel = 6 #Must match "theModel"

Luse = diag(sizeModel)[-1,]

sigmaResid = sigmaIcept = nObserved = nClusters = Fstat = rep(NA, M)
coefEsts = matrix(NA, M, sizeModel)
coefVcovs =list()
for(j in 1:M){
  ii = !is.na(Betas_WBC[j,])
  nObserved[j] = sum(ii)
  pheno0$y = Betas_WBC[j,]

  if(j%%100==0) cat(j,"\n")

  try({
    fit = try(lme(theModel, random=~1|PLATE, data=pheno0[ii,]))

    if(inherits(fit,"try-error")){
      fit = lm(theModel, data=pheno0[ii,])
      fitCoef = fit$coef
      sigmaResid[j] = summary(fit)$sigma
      sigmaIcept[j] = 0
      nClusters[j] = 0
    }
    else{
      fitCoef = fit$coef$fixed
      sigmaResid[j] = fit$sigma
      sigmaIcept[j] = sqrt(getVarCov(fit)[1])
      nClusters[j] = length(fit$coef$random[[1]])
    }
    coefEsts[j,] = fitCoef
    coefVcovs[[j]] = vcov(fit)

    useCoef = Luse %*% fitCoef
    useV = Luse %*% coefVcovs[[j]] %*% t(Luse)
    Fstat[j] = (t(useCoef) %*% solve(useV, useCoef))/sizeModel
  })
}

rownames(coefEsts) = rownames(Betas_WBC)[1:M]
colnames(coefEsts) = names(fitCoef)

results = data.frame (Fstat, sigmaIcept, sigmaResid)
rownames(results) = rownames(Betas_WBC)
annot.results = cbind(annot, results)

#save(file="WBC-Analysis.RData", compress=TRUE, list=c("coefEsts", "coefVcovs", "sigmaIcept", "sigmaResid" , "Luse", "Fstat", "nClusters", "nObserved"))

#############################
#############################
# Part 2: Stochastic Search Algorithm to identify candidate panel of CpG loci

# load(WBC-Analysis)	# only if workspace cleared after part 1
# load (WBCmixLMEobjects.RData) # only if workspace cleared after part 1

autos = !(annot$CHR %in% c("X","Y"))
names(autos) = annot$NAME
autosIx = which(autos)
autoBetas = Betas_WBC[autosIx,] # restrict to autosomal loci

getOptimalSet = function(
  Helems,                                # Hessian elements (coef ^x2 for each locus)
  initialSet,                            # Set of indices for initial locus set
  coreSet = numeric(0),                  # Set of indices for core locus set (always included)
  sampleProb = rep(1,dim(Helems)[3]),    # Sampling probabilities for all loci
  Contrast = NULL,                       # Contrast for determining objective (see below)
  M = 1000,                              # Number of iterations
  seed = 1,                              # Random number seed
  showPlot = TRUE                        # Show progress plots?
){
  nCandidates = dim(Helems)[3]
  if(is.null(Contrast)){               
#DEFAULT CONTRAST is to pick the variance of the first coefficient
    Contrast = matrix(0, 1, dim(Helems)[1])
    Contrast[1,1] = 1
    }
    CpGs = union(initialSet, coreSet)
    
  HH = apply(Helems[,,CpGs], 1:2, sum)
  VV = solve(HH)
  SS = setdiff(1:nCandidates, CpGs)
  PS = sampleProb[SS]/sum(sampleProb[SS])
  
  Vlist = matrix(NA,M+1,8)
  Vlist[1,] = diag(VV)
  ObjHistory = AcceptRate = rep(NA,M+1)
  AcceptHistory = rep(0,M+1)
  ObjHistory[1] = Objective0 = sum(diag(Contrast %*% diag(diag(VV)) %*% t(Contrast)))
  
  set.seed(seed)
  for(m in 1:M){
  	mp1 = m+1
  	i0 = sample(initialSet, 1)
  	i1 = sample(SS, 1, prob=PS)
  	HHc = HH - Helems[,,i0] + Helems[,,i1]
  	#if(any(eigen(HHc,only.values = TRUE)$val<0)) browser()
  	
  	VVc = solve(HHc)
  	Objective1 = sum(diag(Contrast %*% diag(diag(VVc)) %*% t(Contrast)))
  	if (Objective1 < Objective0) {
  		HH = HHc
  		VV = VVc
  		initialSet = union(setdiff(initialSet,i0), i1)
  		SS = union(setdiff(SS,i1), i0)
  		PS = sampleProb[SS]/sum(sampleProb[SS])
  		Objective0 = Objective1
  		AcceptHistory[mp1] = 1
  		}
  		ObjHistory[mp1] = Objective0
  		Vlist[mp1,] = diag(VV)
  		if(m>100) AcceptRate[mp1] = mean(AcceptHistory[mp1-(1:100)])
  		if(showPlot & (m%%10==0)) {
  			plot(c(0,m), c(0,ObjHistory[1]), type="n", xlab="Iteration", ylab="Objective")
  			lines(0:m, ObjHistory[1:mp1], col="blue")
  			lines(0:m, AcceptRate[1:mp1], col="red")
  			}
  			}
  			out = list(loci=CpGs, history=cbind(AcceptHistory=AcceptHistory, Rate=AcceptRate, Objective=ObjHistory),
  			seed=seed)
  			class(out) = "OptimalSet"
  			out
  			}

print.OptimalSet = function(x){
  M = dim(x$history)[1]-1
  cat("Optimal set starting from seed =", x$seed,"with", M, "iterations:\n")
  cat("Final acceptance rate = ", 100*x$history[M+1,2],"%\n", sep="")
  cat("Final objective = ", x$history[M+1,3],"\n\n", sep="")
}

combineOptimalSets = function(setlist){
  allloci = unlist(lapply(setlist, function(u)u$loci))
  allObj = unlist(lapply(setlist, function(u)rep(u$history[dim(u$history)[1],3], length(u$loci))))
  tab = table(allloci)
  best = tapply(allObj, allloci, min)
  out = cbind(locus=as.numeric(names(tab)), TimesChosen=tab, BestObjective=best)
  out[order(-out[,2], out[,3]),]
}


prepareHessianElements = function(coefficients){
  Q=dim(coefficients)[1]
  Helems = array(NA, dim=c(8,8,Q))
  for(m in 1:Q) {
    Helems[,,m] = coefficients[m,] %*% t(coefficients[m,])
  }
  Helems
}

Lwbc = diag(ncol(coefEsts)-1)
rownames(Lwbc) = colnames(Lwbc) = colnames(coefEsts)[2:6]


# Get the contrasted coefficients
LGamma = coefEsts %*% t(Lwbc)
head(LGamma) # view first few

# Subset to autosomal only
LGamma = LGamma[autosIx,]

# Loci unsuitable for GGMA must be removed as informed by Illumina
# LGamma = LGamma[-FAILED,]

Helements = prepareHessianElements(LGamma)

# Set core as those that passed Illumina scoring

# Create order index of F statistics (autosomal only)

Fstat = Fstat[autosIx]
Fstat = Fstat[-FAILED]
FstatOrd = order(-Fstat)

initRuns1 = list()
initSeeds = 1:10

for(r in 1:10){
  initRuns1[[r]] = getOptimalSet(
    Helements,
    initialSet = FstatOrd[1:500],   # 500 largest F stats
    coreSet = CORE,
    Contrast = matrix(1,1,8),
    M = 50000,
    seed = initSeeds[r],
    showPlot = FALSE 
  )
  print(initRuns1[[r]])
}

# Now do it again but start with absolute effect size
EFFSIZEFUN = function(L) order(-apply(abs(L),1,sum))
EffSizeOrd = EFFSIZEFUN(LGamma)


initRuns2 = list()
initSeeds = 11:21

for(r in 1:10){
  initRuns2[[r]] = getOptimalSet(
    Helements,
    initialSet = EffSizeOrd[1:500],   # 500 largest abs effect sizes
    coreSet = CORE,
    Contrast = matrix(1,1,8),
    M = 50000,
    seed = initSeeds[r],
    showPlot = FALSE 
  )
  print(initRuns2[[r]])
}

# Combine the inital sets
run1 = combineOptimalSets(c(initRuns1,initRuns2))
run1[1:100,]  # See the loci most frequently chosen

# Now do it AGAIN but start with top 500 from the combined

initRuns3 = list()
initSeeds = 22:27

for(r in 1:5){
  initRuns3[[r]] = getOptimalSet(
    Helements,
    initialSet = run1[1:500,"locus"],   # 500 most frequently chosen from above
    coreSet = CORE,
    Contrast = matrix(1,1,8),
    M = 50000,
    seed = initSeeds[r],
    showPlot = FALSE 
  )
  print(initRuns3[[r]])
}

run2 = combineOptimalSets(initRuns3)
run2[1:100,]

# Now do the final run for a very long time
NUM_CPGS_VERACODE = 96

runFinal = getOptimalSet(
  Helements,
  initialSet = run2[1:NUM_CPGS_VERACODE,"locus"],
  coreSet = CORE,
  Contrast = matrix(1,1,8),
  M = 500000,
  seed = 35,
  showPlot = TRUE
)

runFinal

lociNUM = runFinal$loci # this is adjusted for the removal of FAILED CpGs from Illumina scoring, and w.r.t *autosomal only*

allNames = auto.annot$NAME[-FAILED]
lociNAME = allNames[runFinal$loci] # Selected CpGs

allSymbols = auto.annot$SYMBOL[-FAILED]
lociSYMBOL = allSymbols[runFinal$loci]

Hess = apply(Helements[,,runFinal$loci], 1:2, sum)
finalSE = sqrt(diag(solve(Hess)))
names(finalSE) = colnames(LGamma)
finalSE

### Plot selected loci
useBetas = autoBetas [-FAILED,]

DMRs = useBetas[lociNUM, which(pheno0$CONTROL != 1)]
samplePheno = subset(pheno0, CONTROL != 1)

tmp = DMRs
rownames(tmp) = c(lociSYMBOL)
colnames(tmp) = samplePheno$SP
heatmap(tmp, scale="n", cexRow = 0.3, cexCol = 0.75,
  col=colorRampPalette(c("yellow","black","blue"),space="Lab")(64))

#save(file="FinalRun.RData", compress=TRUE, list=c("runFinal", "lociNUM","lociNAME", "lociSYMBOL", "finalSE","run1","run2", "initRuns1","initRuns2", "initRuns3"))

# NOTES: this stochastic search was run man times as selected CpG loci needed to be checked with Illumina to determine viability on the GGMA. The final list of loci selected appears in Supplementary Table S1 and also includes CpGs selected for published biological significance


################################################
################################################
# Step 3: perform immune profiling in VeraCode GGMA data

assayDat = read.csv("Veracode_WBC_DMR_betas.csv", head=TRUE, stringsAsFactors=FALSE)
phenoDat = read.csv("Veracode key.csv", head=TRUE, stringsAsFactors=FALSE)
phenoDat$Plate = as.numeric(substring(phenoDat$PlateWell,2,2))

sampIDassay = gsub("^X","", names(assayDat)[-1])
sampIDassay2 = gsub("[.]1$","",sampIDassay)
all(phenoDat$SampleID==sampIDassay2)

S0pheno = phenoDat[phenoDat$SampleType=="Purified WBC subtype",]
S0assay = as.matrix(assayDat[,-1][,phenoDat$SampleType=="Purified WBC subtype"])
colnames(S0assay)= S0pheno$Details
rownames(S0assay) = assayDat[,1]

CellType <- S0pheno$Details
CellType <- gsub("Granuloctyes","Granulocytes",CellType)

# Replace CpG IDs with gene symbols
cgsAtEnd = merge(assayDat, annot[,1:2], sort = FALSE)
all(assayDat$X == cgsAtEnd$X)
assayDat2 = cgsAtEnd[,c(194,3:193)]

# Make S0assay rownames gene symbols
rownames(S0assay) = assayDat2[,1]

tmp = as.matrix(S0assay)
heatmap(tmp, scale="n", cexRow = 0.4, cexCol = 0.55, xlab = "Purified Leukocyte Type", ylab = "Gene", col=colorRampPalette(c("yellow","black","blue"),space="Lab")(64))


thePool <- c(1:96) # Pool of CpGs to consider (by index)

CellType2 <- CellType 

# Forget NK specifics 
CellType2[CellType2=="CD16- NK cells"] <- "NK cells"
CellType2[CellType2=="CD16+ NK cells"] <- "NK cells"
CellType2[CellType2=="CD8- NK cells"] <- "NK cells"
CellType2[CellType2=="CD8+ NK cells"] <- "NK cells"
CellType2[CellType2=="Pan-NK cells"] <- "NK cells"

# Forget T-cell specifics
CellType2[CellType2=="Pan-T cells"] <- "T-cells"
CellType2[CellType2=="CD4+ T-cells"] <- "T-cells"
CellType2[CellType2=="CD8+ T-cells"] <- "T-cells"
CellType2[CellType2=="NKT cells"] <- "T-cells"
CellType2[CellType2=="Tregs"] <- "T-cells"

# Utility function for distinguishing types

distinguish <- function(type, S0, types, stat="p"){
  n = dim(S0)[2]
  z = rep(0,n)
  z[types==type] = 1
  if(stat=="p"){
    out <- apply(S0, 1, function(u)wilcox.test(u~z, exact=FALSE)$p.value)
  }
  else{
    st <- apply(S0, 1, function(u)wilcox.test(u~z, exact=FALSE)$statistic)/prod(table(z))
    out <- pmax(st,1-st)
  }
  attr(out,"sd") <- cbind(apply(S0[,z==0],1,sd), apply(S0[,z==1],1,sd))

  rg0 = apply(S0[,z==0],1,range)
  rg1 = apply(S0[,z==1],1,range)
  md0 = apply(S0[,z==0],1,median)
  md1 = apply(S0[,z==1],1,median)

  attr(out,"sepExtreme") <- ifelse(md1<md0, rg0[1,]-rg1[2,], rg1[1,]-rg0[2,])
  #      pmax(, rg1[1,]-rg0[2,]), pmax(rg0[2,]-rg1[1,], rg1[2,]-rg0[1,]))

  attr(out,"sepMedian") <- abs(md1-md0)  
  out
}

# Utility function for picking the best one by p-value and separation of extremes
pickBest <- function(pObj,n=1){
  pmin = min(pObj)
  ii = which(pObj==pmin)
  exts = attr(pObj,"sepExtreme")[ii]
  ii[order(-exts)[1:n]]
}

# Get p-values for distinguishing each type

pNeut <- distinguish("Neutrophils", S0assay[thePool,], CellType2)
pEos <- distinguish("Eosinophils", S0assay[thePool,], CellType2)
pBaso <- distinguish("Basophils", S0assay[thePool,], CellType2)
pGran <- distinguish("Granulocytes", S0assay[thePool,], CellType2)
pPanT <- distinguish("T-cells", S0assay[thePool,], CellType2)
pBCell <- distinguish("B-cells", S0assay[thePool,], CellType2)
pMono <- distinguish("Monocytes", S0assay[thePool,], CellType2)
pNK <- distinguish("NK cells", S0assay[thePool,], CellType2)

# Order CpGs by p-value

oNeut <- order(pNeut)
oEos <- order(pEos)
oBaso <- order(pBaso)
oGran <- order(pGran)
oPanT <- order(pPanT)
oBCell <- order(pBCell)
oMono <- order(pMono)
oNK <- order(pNK)


cglist = rownames(S0assay[thePool,])

# Plot all VeraCode CpG loci beta values by cell type

par(mfrow = c(3,3))

for(i in 1:96){boxplot(S0assay[i,]~CellType2, las=2, cex.axis = 0.8, ylab = "Beta value", ylim = c(0,1), main = paste(i,":",cglist[i]))}


# get gene numbers in S0assay

for(i in 1:96){paste(i,":",cglist[i])}

# Select the best loci to use based on WBC specific methylation

theInitialSet <- thePool[c(31,32,22,42,29,5,46,43,70,56,67,39,94,7,96,10,35,33,92,95)]

length(theInitialSet)

which(is.na(S0assay[theInitialSet])) == TRUE # no NAs in the reference set

S0assay2 = S0assay
#for(i in 1:96){rownames(S0assay2)[i] = paste(rownames(S0assay)[i], "(",i,")")} # Adds index numbers to gene names

currS0 = S0assay2[theInitialSet,]

##Get CpG IDs to select same loci in Infinium data

annotReOrd = annot[c(94:96,1:93) , ]
all(annotReOrd$Symbol == rownames(S0assay)) # Check correspondence

cgFinalList = annotReOrd$X[theInitialSet]

# Run the projection
source("wbcUtilities.R") 

myParameterization <- wbcParameterize("param.txt",
  assay= currS0, celltypes=CellType2, rowID=S0pheno$SP,
  checkANOVA=TRUE
)

xchk <- xChk <- wbcCrossCheck(currS0, myParameterization)
rownames(xchk) <- CellType

heatmap(xchk, scale="n", col=gray(seq(1,0,-.01)), cexRow = 0.6, cexCol = 0.8, xlab = "Predicted Cell Type", ylab = "Sample")


createCrossCheckMap <- function(filename){
   txt <- read.delim(filename, sep="\t", head=FALSE, stringsAsFactors=FALSE)
   split(txt[[2]],txt[[1]])
}

evalCrossCheck <- function(x, map, types=rownames(x)){
  n = length(map)
  mapvals = names(map)
  s = rep(0, dim(x)[1])
  for(i in 1:n){
    ii = which(types==mapvals[i])
    s[ii] <- apply(x[ii,map[[i]],drop=FALSE],1,sum)
  }
  names(s) = rownames(x)
  s
}


# Identify outiers (those that fail the cross-check)

mapXCheck <- createCrossCheckMap("map-CrossCheck.txt")
xChkSum <- evalCrossCheck(xChk,mapXCheck,CellType2)
specOutliers <- which(xChkSum < 0.85)
heatmap(xChk[-specOutliers,], scale="n", col=gray(seq(1,0,-.01)), cexRow = 0.4, cexCol = 0.7, xlab = "Predicted Cell Type", ylab = "True Cell Type")
length(specOutliers)

# Remove outliers in reference

myParameterization2 <- wbcParameterize("param.txt",
  assay=currS0[,-specOutliers], 
  celltypes=CellType2[-specOutliers], rowID=S0pheno$SampleID[-specOutliers],
  checkANOVA=TRUE
)

xChk2 <- wbcEstimateTarget(currS0[,-specOutliers], currS0, 
  myParameterization2, normalize=TRUE)
annotRow = rep("white", dim(xchk)[1])
annotRow[specOutliers] <- "red"

rownames(xChk2) <- CellType 

heatmap(xChk2, scale="n", col=gray(seq(1,0,-.01)), RowSide=annotRow, cexRow = 0.5, cexCol = 0.6, xlab = "Predicted Cell Type", ylab = "Sample")

## Examine validation data
# CelL DNA Mixtures
S1phenoMix =  phenoDat[phenoDat$SampleType=="Reconstruction mixture",]
S1assayMix = S1assayProfWB = as.matrix(assayDat[theInitialSet,-1][,phenoDat$SampleType=="Reconstruction mixture"])
# Make rownames gene symbols
rownames(S1assayMix) = assayDat2[theInitialSet,1]
all(rownames(S1assayMix) == rownames(currS0))

properOrder = c(5,7,3,4,2,1,6)
S1assayMix = S1assayMix[,properOrder]
S1phenoMix = S1phenoMix[properOrder,]

expMix = read.csv("veracode reconstruction mixture expected.csv", row.names=1)
colnames(expMix) = c("B-cells","Granulocytes","T-cells","NK cells","Monocytes")

estMix <- wbcEstimateTarget(currS0[,-specOutliers], S1assayMix[,], myParameterization2, pretty=TRUE)


S1phenoWB = phenoDat[phenoDat$SampleType=="Profiled whole blood",]
S1assayWB = S1assayProfWB = as.matrix(assayDat[theInitialSet,-1][,phenoDat$SampleType=="Profiled whole blood"])
# Make rownames gene symbols
rownames(S1assayWB) = assayDat2[theInitialSet,1]
all(rownames(S1assayWB) == rownames(currS0))

estWB <- wbcEstimateTarget(currS0[,-specOutliers], S1assayWB[,], myParameterization2, pretty=TRUE)

# Determine expected cell number estimates

setwd(PROFILE)
profWBexp = read.csv("Profiled_Bloods_differential.csv", head=TRUE)

# FACS

leukAfacs = profWBexp$Exp.A.Absolute.No.of.Total.Leukocytes.CD45pos.per.ul
leukBfacs = profWBexp$Exp.B.Absolute.No.of.Total.Leukocytes.CD45pos.per.ul
granBfacs = profWBexp$Exp.B.Absolute.No.of.Granulocytes.per.ul
neutBfacs = profWBexp$Exp.B.Absolute.No.of.Neutrophils.CD45pos.CD15pos.CD16pos.per.ul
eosBfacs = profWBexp$Exp.B.Absolute.No.of.Eosinophils.CD45pos.CD15pos.CD16neg.per.ul
basoBfacs = profWBexp$Exp.B.Absolute.No.of.Basophils.CD45pos.CD123pos.per.ul
monoBfacs = profWBexp$Exp.B.Absolute.No.of.Monocytes.per.ul
TcellAfacs = profWBexp$Exp.A.Absolute.No.of.T.cells.CD45pos.CD3pos.per.ul
CD8TcellAfacs = profWBexp$Exp.A.Absolute.No.of.CD8pos.T.cells.CD45pos.CD3pos.CD8pos.per.ul
CD4TcellAfacs = profWBexp$Exp.A.Absolute.No.of.CD4pos.T.cells.CD45pos.CD3pos.CD4pos.per.ul
BcellBfacs = profWBexp$Exp.B.Absolute.No.of.B.cells.CD45pos.CD19pos.per.ul
NKcellAfacs = profWBexp$Exp.A.Absolute.No.of.Total.NK.Cells.CD45pos.CD56pos.per.ul


granPerFacs = 100*granBfacs/leukBfacs
neutPerFacs = 100*neutBfacs/leukBfacs
eosPerFacs = 100*eosBfacs/leukBfacs
basoPerFacs = 100*basoBfacs/leukBfacs
monoPerFacs = 100*monoBfacs/leukBfacs
TcellPerFacs = 100*TcellAfacs/leukAfacs
CD8TcellPerFacs = 100*CD8TcellAfacs/leukAfacs
CD4TcellPerFacs = 100*CD8TcellAfacs/leukAfacs
BcellPerFacs = 100*BcellBfacs/leukBfacs
NKcellPerFacs = 100*NKcellAfacs/leukAfacs

lymphPerFacs = TcellPerFacs + BcellPerFacs + NKcellPerFacs

allFacs = data.frame(BcellPerFacs,basoPerFacs,eosPerFacs,monoPerFacs,neutPerFacs,NKcellPerFacs,TcellPerFacs)

rownames(allFacs) = c("X3021facs","X3045facs","X3082facs","X3093facs","X3159facs","X3177facs")
colnames(allFacs) = colnames(estWBCitfr)

# Manual differential

neutMan = 100*(profWBexp$Manual.Diff.Percent.Neutrophils.)
eosMan = 100*(profWBexp$Manual.Diff.Percent.Eosinophils)
basoMan = 100*(profWBexp$Manual.Diff.Percent.Basophils)
granMan = 100*(profWBexp$Manual.Diff.Percent.Neutrophils.+ profWBexp$Manual.Diff.Percent.Eosinophils + profWBexp$Manual.Diff.Percent.Basophils)
monoMan = 100*(profWBexp$Manual.Diff.Percent.Monocytes)
lymphMan = 100*(profWBexp$Manual.Diff.Percent.Lymphocytes)

allManDiff = data.frame(basoMan, eosMan, monoMan, neutMan, lymphMan)
rownames(allManDiff) = c("X3021man","X3045man","X3082man","X3093man","X3159man","X3177man")
colnames(allManDiff) = c("Basophils", "Eosinophils", "Monocytes", "Neutrophils", "Lymphocytes")

neutManAb = profWBexp$Manual.Absolute.No.of.Neutrophils.cell.ul
eosManAb = profWBexp$Manual.Absolute.No.of.Eosinophils.cell.ul
basoManAb = profWBexp$Manual.Absolute.No.of.Basophils.cell.ul
granManAb = profWBexp$Manual.Absolute.No.of.Neutrophils.cell.ul+ profWBexp$Manual.Absolute.No.of.Eosinophils.cell.ul + profWBexp$Manual.Absolute.No.of.Basophils.cell.ul
monoManAb = profWBexp$Manual.Absolute.No.of.Monocytes.cell.ul
lymphManAb = profWBexp$Manual.Absolute.No.of.Lymphocytes.cell.ul


# Automated differential

neutCBC = 100*(profWBexp$CBC.Diff.Percent.Neutrophils.cell.ul)
eosCBC = 100*(profWBexp$CBC.Diff.Percent.Eosinophils)
basoCBC = 100*(profWBexp$CBC.Diff.Percent.Basophils)
granCBC = 100*(profWBexp$CBC.Diff.Percent.Neutrophils.cell.ul + profWBexp$CBC.Diff.Percent.Eosinophils + profWBexp$CBC.Diff.Percent.Basophils)
monoCBC = 100*(profWBexp$CBC.Diff.Percent.Monocytes)
lymphCBC = 100*(profWBexp$CBC.Diff.Percent.Lymphocytes)

allAutoDiff = data.frame(basoCBC,eosCBC,monoCBC,neutCBC,lymphCBC)
rownames(allAutoDiff) = c("X3021aut","X3045aut","X3082aut","X3093aut","X3159aut","X3177aut")
colnames(allAutoDiff) = c("Basophils", "Eosinophils", "Monocytes", "Neutrophils", "Lymphocytes")

#################################################
# Subset whole bloods by storage conditions
# HEPARIN

toShowHep <- intersect(grep("Heparin",S1phenoWB$Details),grep("Fresh$",S1phenoWB$Details))
estWBhepfr = estWB[toShowHep,]
estWBhepfr = estWBhepfr[c(5,4,6,1,3,2),] # put the estimates in the proper order

toShowHep4C <- intersect(grep("Heparin",S1phenoWB$Details),grep("o/n at 4C",S1phenoWB$Details))
estWBhep4C = estWB[toShowHep4C,]
estWBhep4C = estWBhep4C[c(5,4,2,6,1,3),] # put the estimates in the proper order

toShowHepRT <- intersect(grep("Heparin",S1phenoWB$Details),grep("o/n at room temp",S1phenoWB$Details))
estWBhepRT = estWB[toShowHepRT,]
estWBhepRT = estWBhepRT[c(2,1,6,3,5,4),] # put the estimates in the proper order

toShowHep80 <- intersect(grep("Heparin",S1phenoWB$Details),grep("o/n at -80C",S1phenoWB$Details))
estWBhep80 = estWB[toShowHep80,]
estWBhep80 = estWBhep80[c(4,6,3,2,1,5),] # put the estimates in the proper order

### CITRATE

toShowCit <- intersect(grep("Citrate",S1phenoWB$Details),grep("Fresh$",S1phenoWB$Details))
estWBCitfr = estWB[toShowCit,]
estWBCitfr = estWBCitfr[c(1,6,3,2,4,5),] # put the estimates in the proper order

toShowCit4C <- intersect(grep("Citrate",S1phenoWB$Details),grep("o/n at 4C",S1phenoWB$Details))
estWBCit4C = estWB[toShowCit4C,]
estWBCit4C = estWBCit4C[c(2,3,1,5,4,6),] # put the estimates in the proper order

toShowCitRT <- intersect(grep("Citrate",S1phenoWB$Details),grep("o/n at room temp",S1phenoWB$Details))
estWBCitRT = estWB[toShowCitRT,]
estWBCitRT = estWBCitRT[c(2,1,6,5,4,3),] # put the estimates in the proper order

toShowCit80 <- intersect(grep("Citrate",S1phenoWB$Details),grep("o/n at -80C",S1phenoWB$Details))
estWBCit80 = estWB[toShowCit80,]
estWBCit80 = estWBCit80[c(5,4,1,3,6,2),] # put the estimates in the proper order

# EDTA

toShowEDTA <- intersect(grep("EDTA",S1phenoWB$Details),grep("Fresh$",S1phenoWB$Details))
estWBEDTAfr = estWB[toShowEDTA,]
estWBEDTAfr = estWBEDTAfr[c(3,1,6,5,2,4),] # put the estimates in the proper order

toShowEDTA4C <- intersect(grep("EDTA",S1phenoWB$Details),grep("o/n at 4C",S1phenoWB$Details))
estWBEDTA4C = estWB[toShowEDTA4C,]
estWBEDTA4C = estWBEDTA4C[c(3,1,2,6,5,4),] # put the estimates in the proper order

toShowEDTART <- intersect(grep("EDTA",S1phenoWB$Details),grep("o/n at room temp",S1phenoWB$Details))
estWBEDTART = estWB[toShowEDTART,]
estWBEDTART = estWBEDTART[c(6,2,1,4,3,5),] # put the estimates in the proper order

toShowEDTA80 <- intersect(grep("EDA",S1phenoWB$Details),grep("o/n at -80C",S1phenoWB$Details))
estWBEDTA80 = estWB[toShowEDTA80,]
estWBEDTA80 = estWBEDTA80[c(6,4,3,2,5,1),] # put the estimates in the proper order

### Four panel figure of X-Y plots

par(mfrow=c(2,2),mar = c(2.5,2,1,1),mgp = c(1.1,0.5,0))

#Mixtures

plot(estMix[1,1], expMix[1,1], xlim = c(0,100), ylim = c(0,100), xlab = "DNA methylation estimate (% of leukocytes)", ylab = "Expected (% of leukocytes)", pch = 0, col = "firebrick1", cex.axis = 0.7, cex.lab = 0.7, cex = 0.7)
abline(0,1, col = "gray46")
for(i in 1:7){
points(estMix[i,1], expMix[i,1], pch = i-1, col = "firebrick1", lwd=1, cex=0.7)
points((estMix[i,2]+estMix[i,3]+estMix[i,5]), expMix[i,2], pch = i-1, col = "darkorchid4", lwd=1, cex=0.7)
points(estMix[i,7], expMix[i,3], pch = i-1, col = "gray0", lwd=1, cex=0.7)
points(estMix[i,6], expMix[i,4], pch = i-1, col = "goldenrod4", lwd=1, cex=0.7)
points(estMix[i,4], expMix[i,5], pch = i-1, col = "green", lwd=1, cex=0.7)
}

legend("topleft", title = "Color Legend",pch = 15, col = c("firebrick1", "darkorchid4", "gray0", "goldenrod4", "green"), legend = c("B-cells", "Granulocytes", "T-cells", "NK cells", "Monocytes"), cex =0.7)
legend("bottomright", title = "Symbol Legend",pch = c(0,1,2,3,4,5,6), col = "black", legend = c("Normal", "T-cell lymphopenia (1)", "T-cell lymphopenia (2)", "Granulocytosis", "Granulocytopenia", "B-cell lymphopenia", "Monocytosis"), pt.lwd = 1, cex = 0.7)

# Manual Differential
plot(estWBhepfr[1,2], basoMan[1], ylab = "Manual differential estimate (% of leukocytes)", xlab = "DNA methylation estimate (% of leukocytes)", xlim = c(0,100), ylim = c(0,100), pch = 0, col = "turquoise2", cex.axis = 0.7, cex.lab = 0.7, cex = 0.7)
abline(0,1, col = "gray46")
for(i in 1:6){
points(estWBhepfr[i,2], basoMan[i], pch = i-1, col = "turquoise2", lwd=1, cex=0.7)
points(estWBhepfr[i,3], eosMan[i], pch = i-1, col = "blue", lwd=1, cex=0.7)
points(estWBhepfr[i,5], neutMan[i], pch = i-1, col = "darkorchid", lwd=1, cex=0.7)
points(estWBhepfr[i,4], monoMan[i], pch=i-1, col = "green", lwd=1, cex=0.7)
points(lymphest[i], lymphMan[i], pch = i-1, col = "orange2", lwd=1, cex=0.7)
}

legend("topleft", title = "Color Legend", pch = 15, col = c("turquoise2", "blue", "darkorchid", "green", "orange2"), legend = c("Basophils", "Eosinophils", "Neutrophils", "Monocytes", "Lymphocytes"), cex =0.7)
legend("bottomright", title = "Symbol Legend",pch = c(0,1,2,3,4,5), col = "black", legend = c("Donor #1", "Donor #2", "Donor #3", "Donor #4", "Donor #5", "Donor #6"), pt.lwd = 1, cex = 0.7)

# Auto Differential
plot(estWBhepfr[2,2], basoCBC[2], ylab = "Automated differential estimate (% of leukocytes)", xlab = "DNA methylation estimate (% of leukocytes)", xlim = c(0,100), ylim = c(0,100), pch = 1, col = "turquoise2", cex.axis = 0.7, cex.lab = 0.7, cex = 0.7)
abline(0,1, col = "gray46")
for (i in 2:6){
points(estWBhepfr[i,2], basoCBC[i], pch = i-1, col = "turquoise2", lwd=1, cex=0.7)
points(estWBhepfr[i,3], eosCBC[i], pch = i-1, col = "blue", lwd=1, cex=0.7)
points(estWBhepfr[i,5], neutCBC[i], pch = i-1, col = "darkorchid", lwd=1, cex=0.7)
points(estWBhepfr[i,4], monoCBC[i], pch=i-1, col = "green", lwd=1, cex=0.7)
points(lymphest[i], lymphCBC[i], pch = i-1, col = "orange2", lwd=1, cex=0.7)
}

legend("topleft", title = "Color Legend", pch = 15, col = c("turquoise2", "blue", "darkorchid", "green", "orange2"), legend = c("Basophils", "Eosinophils", "Neutrophils", "Monocytes", "Lymphocytes"), cex =0.7)
legend("bottomright", title = "Symbol Legend",pch = c(1,2,3,4,5), col = "black", legend = c("Donor #2", "Donor #3", "Donor #4", "Donor #5", "Donor #6"), pt.lwd = 1, cex = 0.7)

#FACS
plot(estWBhepfr[1,1], BcellPerFacs[1],  ylab = "FACS estimate (% of leukocytes)", xlab = "DNA methylation estimate (% of leukocytes)", xlim = c(0,100), ylim = c(0,100), pch = 0, col = "firebrick1", cex.axis = 0.7, cex.lab = 0.7, cex = 0.7)
abline(0,1, col = "gray46")
for(i in 1:6){
points(estWBhepfr[i,1], BcellPerFacs[i], pch = i-1, col = "firebrick1", lwd=1, cex=0.7)
points(estWBhepfr[i,2], basoPerFacs[i], pch = i-1, col = "turquoise2", lwd=1, cex=0.7)
points(estWBhepfr[i,3], eosPerFacs[i], pch = i-1, col = "blue", lwd=1, cex=0.7)
points(estWBhepfr[i,5], neutPerFacs[i], pch = i-1, col = "darkorchid", lwd=1, cex=0.7)
points(estWBhepfr[i,7], TcellPerFacs[i], pch = i-1, col = "gray0", lwd=1, cex=0.7)
points(estWBhepfr[i,6], NKcellPerFacs[i], pch = i-1, col = "goldenrod4", lwd=1, cex=0.7)
points(estWBhepfr[i,4], monoPerFacs[i], pch = i-1, col = "green", lwd=1, cex=0.7)
}

legend("topleft", title = "Color Legend", pch = 15, col = c("firebrick1", "turquoise2", "blue", "darkorchid", "gray0", "goldenrod4", "green"), legend = c("B-cells", "Basophils", "Eosinophils", "Neutrophils", "T-cells", "NK cells", "Monocytes"), cex =0.7)

legend("bottomright", title = "Symbol Legend",pch = c(0,1,2,3,4,5), col = "black", legend = c("Donor #1", "Donor #2", "Donor #3", "Donor #4", "Donor #5", "Donor #6"), pt.lwd = 1, cex = 0.7)

### Three panel traditional methods

par(mfrow=c(1,3),mar = c(2.5,2,1,1),mgp = c(1.1,0.5,0))

# Manual and auto differential
plot(basoMan[1], basoCBC[1], ylab = "Automated differential estimate (% of leukocytes)", xlab = "Manual differential estimate (% of leukocytes)", xlim = c(0,100), ylim = c(0,100), pch = 0, col = "turquoise2", cex.axis = 0.7, cex.lab = 0.7, cex = 0.7)
abline(0,1, col = "gray46")
for (i in 2:6){
points(basoMan[i], basoCBC[i], pch = i-1, col = "turquoise2", lwd=1, cex=0.7)
points(eosMan[i], eosCBC[i], pch = i-1, col = "blue", lwd=1, cex=0.7)
points(neutMan[i], neutCBC[i], pch = i-1, col = "darkorchid", lwd=1, cex=0.7)
points(monoMan[i], monoCBC[i], pch= i-1, col = "green", lwd=1, cex=0.7)
points(lymphMan[i], lymphCBC[i], pch = i-1, col = "orange2", lwd=1, cex=0.7)
}

### FACS and manual differentials
plot(basoMan[1], basoPerFacs[1], ylab = "FACS estimates (% of leukocytes)", xlab = "Manual differential estimates (% of leukocytes)", xlim = c(0,100), ylim = c(0,100), pch = 0, col = "turquoise2", cex.axis = 0.7, cex.lab = 0.7, cex = 0.7)
abline(0,1, col = "gray46")
for (i in 1:6){
points(basoMan[i], basoPerFacs[i], pch = i-1, col = "turquoise2",  lwd=1, cex=0.7)
points(eosMan[i], eosPerFacs[i], pch = i-1, col = "blue",  lwd=1, cex=0.7)
points(neutMan[i], neutPerFacs[i], pch = i-1, col = "darkorchid4",  lwd=1, cex=0.7)
points(monoMan[i], monoPerFacs[i], pch=i-1, col = "green",  lwd=1, cex=0.7)
points(lymphMan[i], lymphPerFacs[i], pch = i-1, col = "orange",  lwd=1, cex=0.7)
}

# CBC automated differential and FACS

plot(basoCBC[2], basoPerFacs[2], xlab = "Automated differential estimate (% of leukocytes)", ylab = "FACS estimate (% of leukocytes)", xlim = c(0,100), ylim = c(0,100), pch = 1, col = "turquoise2", cex.axis = 0.7, cex.lab = 0.7, cex = 0.7)
abline(0,1, col = "gray46")
for (i in 2:6){
points(basoCBC[i], basoPerFacs[i], pch = i-1, col = "turquoise2",  lwd=1, cex=0.7)
points(eosCBC[i], eosPerFacs[i], pch = i-1, col = "blue",  lwd=1, cex=0.7)
points(neutCBC[i], neutPerFacs[i], pch = i-1, col = "darkorchid",  lwd=1, cex=0.7)
points(monoCBC[i], monoPerFacs[i], pch= i-1, col = "green",  lwd=1, cex=0.7)
points(lymphCBC[i], lymphPerFacs[i], pch = i-1, col = "orange2", lwd=1, cex=0.7)
}

# Legend
plot.new()
legend("left", title = "Cell Types", pch = 15, col = c("turquoise2", "blue", "darkorchid", "green", "orange2"), legend = c("Basophils", "Eosinophils", "Neutrophils", "Monocytes", "Lymphocytes"), cex = 0.7)
legend("right", title = "Blood Samples",pch = c(0,1,2,3,4,5), col = "black", legend = c("Donor #1", "Donor #2", "Donor #3", "Donor #4", "Donor #5", "Donor #6"), pt.lwd = 1, cex = 0.7)

# Bland-Altman agreement and RMSE

### Agreement by Bland-Altman

estMixVect = c(estMix[,1],(estMix[,2]+estMix[,3]+estMix[,5]),estMix[,7],estMix[,6],estMix[,4])
expMixVect = c(expMix[,1],expMix[,2],expMix[,3],expMix[,4],expMix[,5])

estVect2 = c(estWBhepfr[,2], estWBhepfr[,3], estWBhepfr[,5], estWBhepfr[,4], lymphest)
manVect = c(basoMan, eosMan, neutMan, monoMan, lymphMan)
autoVect = c(basoCBC, eosCBC, neutCBC, monoCBC, lymphCBC)
estVect = c(estWBhepfr[,1],estWBhepfr[,2],estWBhepfr[,3],estWBhepfr[,4],estWBhepfr[,5],estWBhepfr[,6],estWBhepfr[,7])
facsVect = c(BcellPerFacs, basoPerFacs, eosPerFacs, monoPerFacs, neutPerFacs,NKcellPerFacs,TcellPerFacs)

library(epade)
library(hydroGOF)

par(mfrow=c(2,2))

MixRmse = rmse(estMixVect, expMixVect)
bland.altman.ade(estMixVect, expMixVect, data=NULL, ltext=TRUE, main="WBC DNA Mixtures", xlab="Mean of paired estimates (expected, observed)", ylab="Difference in estimates (expected-observed)", xlim=NULL, ylim=c(-10,10), lwd=2, cex=1, pch=16, lty=c(1,2,2), xticks=NULL, yticks=NULL, col=NULL, tcol=NULL, bgcol=NULL, lcol=c(4,2,2), alpha=NULL, fitline=0, wall=0, v=NULL, h=NULL, span=0.75)
legend("topleft", legend = paste("RMSE =", signif(MixRmse, digits = 3)), bty = "n")

manMethRmse = rmse(manVect, estVect2)
bland.altman.ade(manVect, estVect2, data=NULL, ltext=TRUE, main="Manual diff. vs. DNAm", xlab="Mean of paired estimates (manual diff, DNAm)", ylab="Difference in estimates (Manual diff - DNAm)", xlim=NULL, ylim=c(-10,10), lwd=2, cex=1, pch=16, lty=c(1,2,2), xticks=NULL, yticks=NULL, col=NULL, tcol=NULL, bgcol=NULL, lcol=c(4,2,2), alpha=NULL, fitline=0, wall=0, v=NULL, h=NULL, span=0.75)
legend("topleft", legend = paste("RMSE =", signif(manMethRmse, digits = 3)), bty = "n")

autoMethRmse = rmse(autoVect, estVect2)
bland.altman.ade(autoVect, estVect2, data=NULL, ltext=TRUE, main="Automated diff. vs. DNAm", xlab="Mean of paired estimates (auto diff, DNAm)", ylab="Difference in estimates (auto diff - DNAm)", xlim=NULL, ylim=c(-10,10), lwd=2, cex=1, pch=16, lty=c(1,2,2), xticks=NULL, yticks=NULL, col=NULL, tcol=NULL, bgcol=NULL, lcol=c(4,2,2), alpha=NULL, fitline=0, wall=0, v=NULL, h=NULL, span=0.75)
legend("topleft", legend = paste("RMSE =", signif(autoMethRmse, digits = 3)), bty = "n")

facsMethRmse = rmse(facsVect, estVect)
bland.altman.ade(facsVect, estVect, data=NULL, ltext=TRUE, main="FACS vs. DNAm", xlab="Mean of paired estimates (FACS, DNAm)", ylab="Difference in estimates (FACS - DNAm)", xlim=NULL, ylim=c(-10,10), lwd=2, cex=1, pch=16, lty=c(1,2,2), xticks=NULL, yticks=NULL, tcol=NULL, bgcol=NULL, lcol=c(4,2,2), alpha=NULL, fitline=0, wall=0, v=NULL, h=NULL, span=0.75)
legend("topleft", legend = paste("RMSE =", signif(facsMethRmse, digits = 3)), bty = "n")

# Bland-Altman plot of traditional methods

facsVect2 = c(basoPerFacs, eosPerFacs, neutPerFacs, monoPerFacs, lymphPerFacs)

bland.altman.ade(autoVect, facsVect2, data=NULL, ltext=TRUE, main="Automated diff. vs. FACS", xlab="Mean of paired estimates (auto diff, FACS)", ylab="Difference in estimates (auto diff - FACS)", xlim=NULL, ylim=c(-10,10), lwd=2, cex=1, pch=16, lty=c(1,2,2), xticks=NULL, yticks=NULL, col=NULL, tcol=NULL, bgcol=NULL, lcol=c(4,2,2), alpha=NULL, fitline=0, wall=0, v=NULL, h=NULL, span=0.75)
autoFacsRmse = rmse(autoVect, facsVect2)
legend("topleft", legend = paste("RMSE =", signif(autoFacsRmse, digits = 3)), bty = "n")


bland.altman.ade(manVect, facsVect2, data=NULL, ltext=TRUE, main="Manual diff. vs. FACS", xlab="Mean of paired estimates (manual diff, FACS)", ylab="Difference in estimates (manual diff - FACS)", xlim=NULL, ylim=c(-10,10), lwd=2, cex=1, pch=16, lty=c(1,2,2), xticks=NULL, yticks=NULL, col=NULL, tcol=NULL, bgcol=NULL, lcol=c(4,2,2), alpha=NULL, fitline=0, wall=0, v=NULL, h=NULL, span=0.75)
manFacsRmse = rmse(manVect, facsVect2)
legend("topleft", legend = paste("RMSE =", signif(manFacsRmse, digits = 3)), bty = "n")


bland.altman.ade(manVect, autoVect, data=NULL, ltext=TRUE, main="Manual diff vs. automated diff", xlab="Mean of paired estimates (manual diff, auto diff)", ylab="Difference in estimates (manual diff - auto diff)", xlim=NULL, ylim=c(-10,10), lwd=2, cex=1, pch=16, lty=c(1,2,2), xticks=NULL, yticks=NULL, col=NULL, tcol=NULL, bgcol=NULL, lcol=c(4,2,2), alpha=NULL, fitline=0, wall=0, v=NULL, h=NULL, span=0.75)
manAutoRmse = rmse(manVect, autoVect)
legend("topleft", legend = paste("RMSE =", signif(manAutoRmse, digits = 3)), bty = "n")


# Consider different storage conditions
coagTest1 = cor.test(estWBhepfr, estWBCitfr)
coagTest2 = cor.test(estWBhepfr, estWBEDTAfr)
coagTest3 = cor.test(estWBCitfr, estWBEDTAfr)
hepTest1 = cor.test(estWBhepfr, estWBhepRT)
hepTest2 = cor.test(estWBhepfr, estWBhep4C)
hepTest3 = cor.test(estWBhepfr, estWBhep80)
citTest1 = cor.test(estWBCitfr, estWBCitRT)
citTest2 = cor.test(estWBCitfr, estWBCit4C)
citTest3 = cor.test(estWBCitfr, estWBCit80)
eTest1 = cor.test(estWBEDTAfr, estWBEDTART)
eTest2 = cor.test(estWBEDTAfr, estWBEDTA4C)
eTest3 = cor.test(estWBEDTAfr, estWBEDTA80)

# make 4 panel

par(mfrow=c(2,2),mar = c(2.5,2,1,1),mgp = c(1.1,0.5,0))


# Plot fresh samples for all anticoagulants

plot(estWBhepfr[1,1], estWBCitfr[1],  ylab = "Estimate for sample with citrate or EDTA (% of leukocytes)", xlab = "Estimate for sample with Heparin (% of leukocytes)", main = "Fresh samples", xlim = c(0,100), ylim = c(0,100), pch = 0, col = "firebrick1", cex.axis = 0.7, cex.lab = 0.7, cex = 0.7)
abline(0,1, col = "gray46")
for(i in 1:6){
points(estWBhepfr[i,1], estWBCitfr[i,1], pch = 0, col = "firebrick1", lwd=1, cex=0.7)
points(estWBhepfr[i,2], estWBCitfr[i,2], pch = 0, col = "turquoise2", lwd=1, cex=0.7)
points(estWBhepfr[i,3], estWBCitfr[i,3], pch = 0, col = "blue", lwd=1, cex=0.7)
points(estWBhepfr[i,5], estWBCitfr[i,5], pch = 0, col = "darkorchid", lwd=1, cex=0.7)
points(estWBhepfr[i,7], estWBCitfr[i,7], pch = 0, col = "gray0", lwd=1, cex=0.7)
points(estWBhepfr[i,6], estWBCitfr[i,6], pch = 0, col = "goldenrod4", lwd=1, cex=0.7)
points(estWBhepfr[i,4], estWBCitfr[i,4], pch = 0, col = "green", lwd=1, cex=0.7)
}

for(i in 1:6){
points(estWBhepfr[i,1], estWBEDTAfr[i,1], pch = 1, col = "firebrick1", lwd=1, cex=0.7)
points(estWBhepfr[i,2], estWBEDTAfr[i,2], pch = 1, col = "turquoise2", lwd=1, cex=0.7)
points(estWBhepfr[i,3], estWBEDTAfr[i,3], pch = 1, col = "blue", lwd=1, cex=0.7)
points(estWBhepfr[i,5], estWBEDTAfr[i,5], pch = 1, col = "darkorchid", lwd=1, cex=0.7)
points(estWBhepfr[i,7], estWBEDTAfr[i,7], pch = 1, col = "gray0", lwd=1, cex=0.7)
points(estWBhepfr[i,6], estWBEDTAfr[i,6], pch = 1, col = "goldenrod4", lwd=1, cex=0.7)
points(estWBhepfr[i,4], estWBEDTAfr[i,4], pch = 1, col = "green", lwd=1, cex=0.7)
}

legend("topleft", title = "Cell Type", pch = 15, col = c("firebrick1", "turquoise2", "blue", "darkorchid", "gray0", "goldenrod4", "green"), legend = c("B-cells", "Basophils", "Eosinophils", "Neutrophils", "T-cells", "NK cells", "Monocytes"), cex =0.7)

legend("bottomright", title = "Anticoagulant",pch = c(0,1), col = "black", legend = c("Citrate","EDTA"), pt.lwd = 1, cex = 0.7)


# HEPARIN

plot(estWBhepfr[1,1], estWBhep4C[1],  ylab = "Estimate for sample stored overnight (% of leukocytes)", xlab = "Estimate for fresh sample (% of leukocytes)", main = "Heparin anticoagulant", xlim = c(0,100), ylim = c(0,100), pch = 0, col = "firebrick1", cex.axis = 0.7, cex.lab = 0.7, cex = 0.7)
abline(0,1, col = "gray46")
for(i in 1:6){
points(estWBhepfr[i,1], estWBhep4C[i,1], pch = 0, col = "firebrick1", lwd=1, cex=0.7)
points(estWBhepfr[i,2], estWBhep4C[i,2], pch = 0, col = "turquoise2", lwd=1, cex=0.7)
points(estWBhepfr[i,3], estWBhep4C[i,3], pch = 0, col = "blue", lwd=1, cex=0.7)
points(estWBhepfr[i,5], estWBhep4C[i,5], pch = 0, col = "darkorchid", lwd=1, cex=0.7)
points(estWBhepfr[i,7], estWBhep4C[i,7], pch = 0, col = "gray0", lwd=1, cex=0.7)
points(estWBhepfr[i,6], estWBhep4C[i,6], pch = 0, col = "goldenrod4", lwd=1, cex=0.7)
points(estWBhepfr[i,4], estWBhep4C[i,4], pch = 0, col = "green", lwd=1, cex=0.7)
}

for(i in 1:6){
points(estWBhepfr[i,1], estWBhepRT[i,1], pch = 1, col = "firebrick1", lwd=1, cex=0.7)
points(estWBhepfr[i,2], estWBhepRT[i,2], pch = 1, col = "turquoise2", lwd=1, cex=0.7)
points(estWBhepfr[i,3], estWBhepRT[i,3], pch = 1, col = "blue", lwd=1, cex=0.7)
points(estWBhepfr[i,5], estWBhepRT[i,5], pch = 1, col = "darkorchid", lwd=1, cex=0.7)
points(estWBhepfr[i,7], estWBhepRT[i,7], pch = 1, col = "gray0", lwd=1, cex=0.7)
points(estWBhepfr[i,6], estWBhepRT[i,6], pch = 1, col = "goldenrod4", lwd=1, cex=0.7)
points(estWBhepfr[i,4], estWBhepRT[i,4], pch = 1, col = "green", lwd=1, cex=0.7)
}

for(i in 1:6){
points(estWBhepfr[i,1], estWBhep80[i,1], pch = 2, col = "firebrick1", lwd=1, cex=0.7)
points(estWBhepfr[i,2], estWBhep80[i,2], pch = 2, col = "turquoise2", lwd=1, cex=0.7)
points(estWBhepfr[i,3], estWBhep80[i,3], pch = 2, col = "blue", lwd=1, cex=0.7)
points(estWBhepfr[i,5], estWBhep80[i,5], pch = 2, col = "darkorchid", lwd=1, cex=0.7)
points(estWBhepfr[i,7], estWBhep80[i,7], pch = 2, col = "gray0", lwd=1, cex=0.7)
points(estWBhepfr[i,6], estWBhep80[i,6], pch = 2, col = "goldenrod4", lwd=1, cex=0.7)
points(estWBhepfr[i,4], estWBhep80[i,4], pch = 2, col = "green", lwd=1, cex=0.7)
}

legend("topleft", title = "Cell Type", pch = 15, col = c("firebrick1", "turquoise2", "blue", "darkorchid", "gray0", "goldenrod4", "green"), legend = c("B-cells", "Basophils", "Eosinophils", "Neutrophils", "T-cells", "NK cells", "Monocytes"), cex =0.7)

legend("bottomright", title = "Sample Storage",pch = c(1,0,2), col = "black", legend = c("o/n at room temp","o/n at 4C", "o/n at -80C"), pt.lwd = 1, cex = 0.7)



# EDTA

plot(estWBEDTAfr[1,1], estWBEDTA4C[1],  ylab = "Estimate for sample stored overnight (% of leukocytes)", xlab = "Estimate for fresh sample (% of leukocytes)", main = "EDTA anticoagulant", xlim = c(0,100), ylim = c(0,100), pch = 0, col = "firebrick1", cex.axis = 0.7, cex.lab = 0.7, cex = 0.7)
abline(0,1, col = "gray46")
for(i in 1:6){
points(estWBEDTAfr[i,1], estWBEDTA4C[i,1], pch = 0, col = "firebrick1", lwd=1, cex=0.7)
points(estWBEDTAfr[i,2], estWBEDTA4C[i,2], pch = 0, col = "turquoise2", lwd=1, cex=0.7)
points(estWBEDTAfr[i,3], estWBEDTA4C[i,3], pch = 0, col = "blue", lwd=1, cex=0.7)
points(estWBEDTAfr[i,5], estWBEDTA4C[i,5], pch = 0, col = "darkorchid", lwd=1, cex=0.7)
points(estWBEDTAfr[i,7], estWBEDTA4C[i,7], pch = 0, col = "gray0", lwd=1, cex=0.7)
points(estWBEDTAfr[i,6], estWBEDTA4C[i,6], pch = 0, col = "goldenrod4", lwd=1, cex=0.7)
points(estWBEDTAfr[i,4], estWBEDTA4C[i,4], pch = 0, col = "green", lwd=1, cex=0.7)
}

for(i in 1:6){
points(estWBEDTAfr[i,1], estWBEDTART[i,1], pch = 1, col = "firebrick1", lwd=1, cex=0.7)
points(estWBEDTAfr[i,2], estWBEDTART[i,2], pch = 1, col = "turquoise2", lwd=1, cex=0.7)
points(estWBEDTAfr[i,3], estWBEDTART[i,3], pch = 1, col = "blue", lwd=1, cex=0.7)
points(estWBEDTAfr[i,5], estWBEDTART[i,5], pch = 1, col = "darkorchid", lwd=1, cex=0.7)
points(estWBEDTAfr[i,7], estWBEDTART[i,7], pch = 1, col = "gray0", lwd=1, cex=0.7)
points(estWBEDTAfr[i,6], estWBEDTART[i,6], pch = 1, col = "goldenrod4", lwd=1, cex=0.7)
points(estWBEDTAfr[i,4], estWBEDTART[i,4], pch = 1, col = "green", lwd=1, cex=0.7)
}

for(i in 1:6){
points(estWBEDTAfr[i,1], estWBEDTA80[i,1], pch = 2, col = "firebrick1", lwd=1, cex=0.7)
points(estWBEDTAfr[i,2], estWBEDTA80[i,2], pch = 2, col = "turquoise2", lwd=1, cex=0.7)
points(estWBEDTAfr[i,3], estWBEDTA80[i,3], pch = 2, col = "blue", lwd=1, cex=0.7)
points(estWBEDTAfr[i,5], estWBEDTA80[i,5], pch = 2, col = "darkorchid", lwd=1, cex=0.7)
points(estWBEDTAfr[i,7], estWBEDTA80[i,7], pch = 2, col = "gray0", lwd=1, cex=0.7)
points(estWBEDTAfr[i,6], estWBEDTA80[i,6], pch = 2, col = "goldenrod4", lwd=1, cex=0.7)
points(estWBEDTAfr[i,4], estWBEDTA80[i,4], pch = 2, col = "green", lwd=1, cex=0.7)
}

legend("topleft", title = "Cell Type", pch = 15, col = c("firebrick1", "turquoise2", "blue", "darkorchid", "gray0", "goldenrod4", "green"), legend = c("B-cells", "Basophils", "Eosinophils", "Neutrophils", "T-cells", "NK cells", "Monocytes"), cex =0.7)

legend("bottomright", title = "Sample Storage",pch = c(1,0,2), col = "black", legend = c("o/n at room temp","o/n at 4C", "o/n at -80C"), pt.lwd = 1, cex = 0.7)

# Citrate

plot(estWBCitfr[1,1], estWBCit4C[1],  ylab = "Estimate for sample stored overnight (% of leukocytes)", xlab = "Estimate for fresh sample (% of leukocytes)", main = "Citrate anticoagulant", xlim = c(0,100), ylim = c(0,100), pch = 0, col = "firebrick1", cex.axis = 0.7, cex.lab = 0.7, cex = 0.7)
abline(0,1, col = "gray46")
for(i in 1:6){
points(estWBCitfr[i,1], estWBCit4C[i,1], pch = 0, col = "firebrick1", lwd=1, cex=0.7)
points(estWBCitfr[i,2], estWBCit4C[i,2], pch = 0, col = "turquoise2", lwd=1, cex=0.7)
points(estWBCitfr[i,3], estWBCit4C[i,3], pch = 0, col = "blue", lwd=1, cex=0.7)
points(estWBCitfr[i,5], estWBCit4C[i,5], pch = 0, col = "darkorchid", lwd=1, cex=0.7)
points(estWBCitfr[i,7], estWBCit4C[i,7], pch = 0, col = "gray0", lwd=1, cex=0.7)
points(estWBCitfr[i,6], estWBCit4C[i,6], pch = 0, col = "goldenrod4", lwd=1, cex=0.7)
points(estWBCitfr[i,4], estWBCit4C[i,4], pch = 0, col = "green", lwd=1, cex=0.7)
}

for(i in 1:6){
points(estWBCitfr[i,1], estWBCitRT[i,1], pch = 1, col = "firebrick1", lwd=1, cex=0.7)
points(estWBCitfr[i,2], estWBCitRT[i,2], pch = 1, col = "turquoise2", lwd=1, cex=0.7)
points(estWBCitfr[i,3], estWBCitRT[i,3], pch = 1, col = "blue", lwd=1, cex=0.7)
points(estWBCitfr[i,5], estWBCitRT[i,5], pch = 1, col = "darkorchid", lwd=1, cex=0.7)
points(estWBCitfr[i,7], estWBCitRT[i,7], pch = 1, col = "gray0", lwd=1, cex=0.7)
points(estWBCitfr[i,6], estWBCitRT[i,6], pch = 1, col = "goldenrod4", lwd=1, cex=0.7)
points(estWBCitfr[i,4], estWBCitRT[i,4], pch = 1, col = "green", lwd=1, cex=0.7)
}

for(i in 1:6){
points(estWBCitfr[i,1], estWBCit80[i,1], pch = 2, col = "firebrick1", lwd=1, cex=0.7)
points(estWBCitfr[i,2], estWBCit80[i,2], pch = 2, col = "turquoise2", lwd=1, cex=0.7)
points(estWBCitfr[i,3], estWBCit80[i,3], pch = 2, col = "blue", lwd=1, cex=0.7)
points(estWBCitfr[i,5], estWBCit80[i,5], pch = 2, col = "darkorchid", lwd=1, cex=0.7)
points(estWBCitfr[i,7], estWBCit80[i,7], pch = 2, col = "gray0", lwd=1, cex=0.7)
points(estWBCitfr[i,6], estWBCit80[i,6], pch = 2, col = "goldenrod4", lwd=1, cex=0.7)
points(estWBCitfr[i,4], estWBCit80[i,4], pch = 2, col = "green", lwd=1, cex=0.7)
}

legend("topleft", title = "Cell type", pch = 15, col = c("firebrick1", "turquoise2", "blue", "darkorchid", "gray0", "goldenrod4", "green"), legend = c("B-cells", "Basophils", "Eosinophils", "Neutrophils", "T-cells", "NK cells", "Monocytes"), cex =0.7)

legend("bottomright", title = "Sample Storage",pch = c(1,0,2), col = "black", legend = c("o/n at room temp","o/n at 4C", "o/n at -80C"), pt.lwd = 1, cex = 0.7)

# Include additional samples on crosscheck heatmap
####

WBtoXchk = estWBhepfr
rownames(WBtoXchk) = c("WB3021", "WB3045", "WB3082", "WB3093", "WB3159", "WB3177")

mixToXchk = estMix
rownames(mixToXchk) = c("Normal mix", "T-cell lymphopenia 1", "T-cell lymphopenia 2", "Granulocytosis", "Granulocytopenia", "B-cell lymphopenia", "Monocytosis")

FullXchk = rbind(xChk2*100, WBtoXchk, mixToXchk)

heatmap(FullXchk, scale="n", col=gray(seq(1,0,-.01)), cexRow = 0.4, cexCol = 0.6, xlab = "Predicted Cell Type", ylab = "Sample")

################################################
################################################
# Step 3 Perform immune profiling in Infinium data

S0pheno = pheno0[which(pheno0$CONTROL != 1),]
S0assay = as.matrix(Betas_WBC)

# Define cell types
CellType = ifelse(S0pheno$PanTcell == 1, "T-cells", ifelse(S0pheno$Bcell == 1, "B-cells", ifelse(S0pheno$NK == 1, "NK cells", ifelse(S0pheno$Mono == 1, "Monocytes", ifelse(S0pheno$Neut == 1, "Neutrophils", ifelse(S0pheno$Eos == 1, "Eosinophils", ifelse(S0pheno$Baso == 1, "Basophils", ifelse(S0pheno$Gran == 1, "Granulocytes", NA))))))))
colnames(S0assay) = CellType

# Define the pool of CpGs to work with
thePool <- c(1:500)

# Utility function for distinguishing types
distinguish <- function(type, S0, types, stat="p"){
  n = dim(S0)[2]
  z = rep(0,n)
  z[types==type] = 1
  if(stat=="p"){
    out <- apply(S0, 1, function(u)wilcox.test(u~z, exact=FALSE)$p.value)
  }
  else{
    st <- apply(S0, 1, function(u)wilcox.test(u~z, exact=FALSE)$statistic)/prod(table(z))
    out <- pmax(st,1-st)
  }
  attr(out,"sd") <- cbind(apply(S0[,z==0],1,sd), apply(S0[,z==1],1,sd))

  rg0 = apply(S0[,z==0],1,range)
  rg1 = apply(S0[,z==1],1,range)
  md0 = apply(S0[,z==0],1,median)
  md1 = apply(S0[,z==1],1,median)

  attr(out,"sepExtreme") <- ifelse(md1<md0, rg0[1,]-rg1[2,], rg1[1,]-rg0[2,])
  #      pmax(, rg1[1,]-rg0[2,]), pmax(rg0[2,]-rg1[1,], rg1[2,]-rg0[1,]))

  attr(out,"sepMedian") <- abs(md1-md0)  
  out
}

# Utility function for picking the best one by p-value and separation of extremes
pickBest <- function(pObj,n=1){
  pmin = min(pObj)
  ii = which(pObj==pmin)
  exts = attr(pObj,"sepExtreme")[ii]
  ii[order(-exts)[1:n]]
}

# Get p-values for distinguishing each type
pNeut <- distinguish("Neutrophils", S0assay[thePool,], CellType2)
pEos <- distinguish("Eosinophils", S0assay[thePool,], CellType2)
pBaso <- distinguish("Basophils", S0assay[thePool,], CellType2)
pGran <- distinguish("Granulocytes", S0assay[thePool,], CellType2)
pPanT <- distinguish("T-cells", S0assay[thePool,], CellType2)
pBCell <- distinguish("B-cells", S0assay[thePool,], CellType2)
pMono <- distinguish("Monocytes", S0assay[thePool,], CellType2)
pNK <- distinguish("NK cells", S0assay[thePool,], CellType2)

# Order CpGs by p-value
oNeut <- order(pNeut)
oEos <- order(pEos)
oBaso <- order(pBaso)
oGran <- order(pGran)
oPanT <- order(pPanT)
oBCell <- order(pBCell)
oMono <- order(pMono)
oNK <- order(pNK)

#### Create a list of either CpG names or gene symbols
#cglist = rownames(S0assay[thePool,]) # CpG names
cglist = annotS0$SYMBOL[thePool] # gene symbols

##### SELECT THE BEST SET OF CPG LOCI
##### Identified from the above ranked lists

theInitialSet <- thePool[unique(c(oBaso[c(6,9)], oEos[c(3,18)], oNeut[c(19,15,13,4,10)], oPanT[c(1,2,12,24,25,22,22,83)], oBCell[c(2,3,4,5,1,8,15,16)], oNK[c(1,3)], oMono[c(1,2,28,39,6,9,25,3)]))]

currS0 = S0assay[theInitialSet,]

# change NA values to cell type mean for that locus

anyMissing = apply(is.na(currS0),1,any) #per CpG
if(any(anyMissing>0)) {
	for(i in 1:dim(currS0)[1]) if(anyMissing[i]>0){
		mu =tapply(currS0[i,], CellType2,mean,na.rm=TRUE)
		for (j in which(is.na(currS0[i,]))) {
			currS0[i,j] = mu[CellType2[j]]
			}
			}
			}

source("wbcUtilities.R")

myParameterization <- wbcParameterize("param.txt",
  assay= currS0, celltypes=CellType2, rowID=S0pheno$SP,
  checkANOVA=TRUE
)
xchk <- xChk <- wbcCrossCheck(currS0, myParameterization)
rownames(xchk) <- CellTypeDetail
heatmap(xchk, scale="n", col=gray(seq(1,0,-.01)), cexRow = 0.6, cexCol = 0.8, xlab = "Predicted Cell Type", ylab = "Sample")
createCrossCheckMap <- function(filename){
   txt <- read.delim(filename, sep="\t", head=FALSE, stringsAsFactors=FALSE)
   split(txt[[2]],txt[[1]])
}
evalCrossCheck <- function(x, map, types=rownames(x)){
  n = length(map)
  mapvals = names(map)
  s = rep(0, dim(x)[1])
  for(i in 1:n){
    ii = which(types==mapvals[i])
    s[ii] <- apply(x[ii,map[[i]],drop=FALSE],1,sum)
  }
  names(s) = rownames(x)
  s
}


# Identify outiers (those that fail the cross-check)
mapXCheck <- createCrossCheckMap("map-CrossCheck.txt")
xChkSum <- evalCrossCheck(xChk,mapXCheck,CellType2)
specOutliers <- which(xChkSum < 0.85)
heatmap(xChk[-specOutliers,], scale="n", col=gray(seq(1,0,-.01)), cexRow = 0.6, cexCol = 0.7, xlab = "Predicted Cell Type", ylab = "True Cell Type")
length(specOutliers)


# Remove outliers in reference
myParameterization2 <- wbcParameterize("param.txt",
  assay=currS0[,-specOutliers], 
  celltypes=CellType2[-specOutliers], rowID=S0pheno$SampleID[-specOutliers],
  checkANOVA=TRUE
)
xChk2 <- wbcEstimateTarget(currS0[,-specOutliers], currS0, 
  myParameterization2, normalize=TRUE)
annotRow = rep("white", dim(xchk)[1])
annotRow[specOutliers] <- "red"
rownames(xChk2) <- CellTypeDetail

heatmap(xChk2, scale="n", col=gray(seq(1,0,-.01)), RowSide=annotRow, cexRow = 0.7, cexCol = 0.7, xlab = "Predicted Cell Type", ylab = "Sample")

### Run projection in target data sets

S1assay # Infinium betas for the target samples
S1pheno # covariate data for the target samples
all(S1pheno$AR == colnames(S1assay)) # check correspondance

# Must subset S1assay to only loci in currS0 (the reference set)
S1assay2 = subset(S1assay, rownames(S1assay) %in% rownames(currS0))
# Now reorder S1assay to be exactly the same as S0assay CpG loci
cglist = rownames(currS0)
S1order = rep(NA,34)
for (i in 1:34){
	S1order[i] = which(rownames(S1assay2) %in% cglist[i])
}
S1assay3 = as.matrix(S1assay2[S1order,])
all(rownames(currS0) == rownames(S1assay3)) # check correspondence
colnames(S1assay3) = S1pheno$Details
# Generate estimates
S1est <- wbcEstimateTarget(currS0[,-specOutliers], S1assay3, myParameterization2, pretty=TRUE)
estSum <- apply(as.matrix(S1est),1,sum)/100
