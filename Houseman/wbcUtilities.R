#######################################
# Includes new wbcProject function from 12/13/2012

wbcProject <- function(Y, coefWBC, contrastWBC=NULL, nonnegative=TRUE){
	
	if(is.null(contrastWBC)) Xmat = coefWBC
	else Xmat = coefWBC %*% t(contrastWBC)
	
	nCol = dim(Xmat)[2]
	nSubj = dim(Y)[2]
	
	mixCoef = matrix(0, nSubj, nCol)
	rownames(mixCoef) = colnames(Y)
	colnames(mixCoef) = colnames(Xmat)
	if(nonnegative){
		library(quadprog)
		Amat = cbind(rep(-1,nCol), diag(nCol))
		b0vec = c(-1,rep(0,nCol))
		for(i in 1:nSubj){
			obs = which(!is.na(Y[,i]))
			Dmat = t(Xmat[obs,])%*%Xmat[obs,]
			mixCoef[i,] = solve.QP(Dmat, t(Xmat[obs,])%*%Y[obs,i], Amat, b0vec, meq = 0)$sol
		}
	}
	else{
		for(i in 1:nSubj){
			obs = which(!is.na(Y[,i]))
			Dmat = t(Xmat[obs,])%*%Xmat[obs,]
			mixCoef[i,] = solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
		}
	}
	return(mixCoef)
}


wbcParameterize <- function(paramFile, assay, celltypes, rowID=celltypes, 
    adjust=NULL, checkANOVA=FALSE){
  param <- read.delim(paramFile, sep="\t", head=TRUE, stringsAsFactors=FALSE)
  if(is.null(param$S0types)) stop("Cannot find variable S0types.")
  if(is.null(param$S1report)) stop("Cannot find variable S1report.")
  if(length(missType <- setdiff(celltypes, param$S0types))>0) {
    print(missType)
    stop("Some types are not represented in this parameterization.")
  }
  header <- strsplit(readLines(paramFile,1),"\t")[[1]]
  rmv <- c(grep("^S0types$", names(param)),grep("^S1report$", names(param)))
  designTab <- as.matrix(param[,-rmv])
  colnames(designTab) <- header[-rmv]
  rownames(designTab) <- param$S0types
  designMatrix <- designTab[celltypes,]
  reportMatrix <- designTab[param$S1report==1,]
  rownames(designMatrix) <- rowID

  if(checkANOVA){
    if(all(apply(designMatrix,1,sum)==1)) cat("All rows sum to 1\n")
    else warning("Not all rows sum to 1!\n")
  }

  if(!is.null(adjust)){
    if(!is.matrix(adjust)){
      adjust = matrix(adjust)
    }
    designMatrix <- cbind(designMatrix, adjust)
    reportMatrix <- cbind(reportMatrix, matrix(0, dim(reportMatrix)[1], dim(adjust)[2]))
  }
  list(design=designMatrix, report=reportMatrix)  
}


wbcEstimateValidationMatrix <- function(assay, param){
  x = param$design
  t(solve(t(x)%*%x, t(x)%*%t(assay)))
}

wbcCrossCheck <- function(assay, param, normalize=TRUE, rowID=NULL, pretty=FALSE, digits=0){
  B = wbcEstimateValidationMatrix(assay, param)
  omega = wbcProject(assay, B, param$report)
  if(normalize) omega = (1/apply(omega,1,sum))*omega
  if(!is.null(rowID)) rownames(omega) = rowID
  if(pretty){
    omega = round(100*omega, digits)
    omega = omega[order(rownames(omega)),]
  }
  omega
}

wbcEstimateTarget <- function(assayS0, assayS1, param, normalize=TRUE, pretty=FALSE, digits=1){
  B = wbcEstimateValidationMatrix(assayS0, param)
  omega = wbcProject(assayS1, B, param$report)
  if(normalize) omega = (1/apply(omega,1,sum))*omega 
  if(pretty){
    omega = round(100*omega, digits)
  }
  omega
}
