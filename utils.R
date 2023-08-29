get_toptables <- function(comparison, designlevels, fit, pval, lgfch) {
  cont.matrix<-makeContrasts(contrasts =  comparison , levels=designlevels)
  fit2 <- contrasts.fit(fit,cont.matrix)
  fit2 <- eBayes(fit2, trend=TRUE)
  top <- topTable(fit2, number=Inf)
  
  top$"log10q" <- -log10(top$adj.P.Val)
  #top$"Sig" <- top$adj.P.Val <= 0.05 & abs(top$logFC) >= 1.5
  top$"Sig" <- top$adj.P.Val <= pval & abs(top$logFC) >= lgfch
  top$"Up" <- top$logFC > 0
  top$"color" <- 'Not.Significant'
  
  top<-mutate(top, color=ifelse(Sig==T & Up==T, 'Up.Significant','Not.Significant'))
  top<-mutate(top, color=ifelse(Sig==T & Up==F, 'Down.Significant',top$color))
  top<-mutate(top, color=ifelse(Sig==F, 'Not.Significant', top$color))
  
  # ct.means<-makeContrasts(Normal.means=normal,
  #                         SCC.means=SCC,
  #                         levels=designlevels)
  # group.means<-eBayes(contrasts.fit(fit,ct.means),trend=TRUE)
  # top2.means <- as.data.frame(group.means$coefficients)[rownames(top),]
  # top2.means$"logFC" <- top2.means$SCC.means - top2.means$Normal.means
  # identical(top2$logFC, top2.means$logFC)
  # top2 <- cbind(top, top2.means)
  
  return(top)
}

calc_toptables <- function(contrastCombo, values) {
  target <- values$target
  design <- values$design
  
  # Filter 0 Sum rows
  target2 <- target[rowSums(target)>0,] 
  # Filter genes where 80% or more samples have counts < 10 
  target.final <- target2[!apply(target2, 1, function(x){ length(x[x < 10]) >= 0.8 * length(x)}),] 
  genenames<-rownames(target.final)
  
  if (values$reference == "Human") {
    load("reference/grch37.genemap.Rdata")
    values$genemap <- genemap.1
    idx <- match(genenames, genemap.1$ensembl_gene_id )
    genesymbol <- genemap.1$external_gene_name[idx]
  }
  if (values$reference == "Mouse") {
    load("reference/mmusculus.genemap.Rdata")
    values$genemap <- genemap.1
    idx <- match(genenames, genemap.1$ensembl_gene_id )
    genesymbol <- genemap.1$mgi_symbol[idx]
    
  }
  if (values$reference == "Other") {
    genesymbol <- rep("No Information",nrow(target.final))
  }
  
  if (values$batch != "No Batch Correction needed.") {
    adjusted_target <- ComBat_seq(as.matrix(target.final), batch=design[,values$batch], group=NULL)  
    counts <- DGEList(counts=adjusted_target,genes=genesymbol) 
  } else {
    counts <- DGEList(counts=target.final,genes=genesymbol) # Create Limma data object
  }
  isexpr <- rowSums(cpm(counts)>1) >= (ncol(counts)*0.7) # Change this to 70% of nColumns
  counts <- counts[isexpr,keep.lib.sizes=FALSE]
  counts <- calcNormFactors(counts)
  
  # Design matrix
  f <-factor(design$condition)
  design.limma <-model.matrix(~0+f)
  colnames(design.limma)<-levels(f)
  design.limma <- as.data.frame(design.limma)
  rownames(design.limma) <- colnames(target.final)
  
  # Voom
  v <- voom(counts, design.limma)
  
  # Linear Modeling after Voom
  fit <- lmFit(v, design.limma)
  
  # TopTables
  tops <- lapply(contrastCombo, get_toptables, design.limma, fit, values$pval, values$lgfch)
  names(tops) <- names(contrastCombo)
  
  
  # Update ReactiveVals
  values$tops <- tops
  values$v <- v
  
  # PCA data (AllSamples)
  pcaData <- prcomp(t(v$E))
  pcaData <- data.frame(cbind(pcaData$x[,c(1,2)], var.explained=pcaData$sdev^2/sum(pcaData$sdev^2)))
  pcaData <- cbind(pcaData, "condition"=design[rownames(pcaData), "condition"])
  values$pcaData <- pcaData
  
  # PCA Pairwise
  pairPCA <- lapply(contrastCombo, getPairPCA, v$E, design, values$pval, values$lgfch)
  names(pairPCA) <- names(contrastCombo)
  
  values$pairPCA <- pairPCA
}

getPairPCA <- function(comparison, exprs, design, pval, lgfch) {
  sel <- unlist(strsplit(comparison, split = "-"))
  
  selSmp <- rownames(design[design$condition %in% sel,])
  
  pcaData <- prcomp(t(exprs[,selSmp]))
  pcaData <- data.frame(cbind(pcaData$x[,c(1,2)], var.explained=pcaData$sdev^2/sum(pcaData$sdev^2)))
  pcaData <- cbind(pcaData, "condition"=design[rownames(pcaData), "condition"])
  
  return(pcaData)
}

pop_Error <- function(title, msg) {
  showModal(modalDialog(
    title = title,
    h4(msg)
  ))
}

tpmNorm <- function(target, len, genemap) {
  target <- target[rowSums(target) > 0,]
  targetTPM <- target

  len <- len[rownames(target),]
  len$KB <- len$Length/1000
  
  
  # Divide each gene by transcript length
  # (Q::Length Reference comparable? GRCh37 to CIBERSORTx Knowledge Base)
  expPKB <- apply(targetTPM, 2, function(x){ x / len$KB } )
  # Divide by the transcript length
  targetTPM <- apply(expPKB, 2, function(x) { x / sum(x) * 1E6})
  
  genemap.1 <- genemap
  targetTPM.m <- as.data.frame(cbind("Gene"=genemap.1[match(rownames(targetTPM), genemap.1$ensembl_gene_id),"external_gene_name"],
                                     targetTPM))
  targetTPM.m <- targetTPM.m[!duplicated(targetTPM.m$Gene),]
  rownames(targetTPM.m) <- targetTPM.m$Gene
  targetTPM.m <- targetTPM.m[,-1]
  return(targetTPM.m)
}

CoreAlg <- function(X, y, absolute, abs_method){
  
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  
  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  if(!absolute || abs_method == 'sig.score') w <- (q/sum(q)) #relative space (returns fractions)
  if(absolute && abs_method == 'no.sumto1') w <- q #absolute space (returns scores)
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}

doPerm <- function(perm, X, Y, absolute, abs_method){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  
  while(itor <= perm){
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    
    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)
    
    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr, absolute, abs_method)
    
    mix_r <- result$mix_r
    
    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}
    
    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}

CIBERSORT <- function(sig_matrix, mixture, perm=100, QN=TRUE, absolute=FALSE, abs_method='sig.score', interruptCallback=function(){}, progressSet=function(){}, progressStart=0.0, progressScale=1.0, currStep=1, totalSteps=1){
  
  if(absolute && abs_method != 'no.sumto1' && abs_method != 'sig.score') stop("abs_method must be set to either 'sig.score' or 'no.sumto1'")
  
  X <- sig_matrix
  Y <- mixture
  
  #to prevent crashing on duplicated gene symbols, add unique numbers to identical names
  dups <- dim(Y)[1] - length(unique(rownames(Y)))
  if(dups > 0) {
    warning(paste(dups," duplicated gene symbol(s) found in mixture file!",sep=""))
    rownames(Y) <- make.names(rownames(Y), unique=TRUE)
  }
  
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  
  #order
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]
  
  P <- perm #number of permutations
  
  #anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}
  
  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  #store original mixtures
  Yorig <- Y
  Ymedian <- max(median(Y),1)
  
  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  table(YintX)
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]
  
  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))
  
  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y, absolute, abs_method)$dist)}
  
  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  if(absolute) header <- c(header, paste('Absolute score (',abs_method,')',sep=""))
  
  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  
  #iterate through mixtures
  time.eta <- "This step can take a while"
  while(itor <= mixtures){
    start.time <- Sys.time()
    
    interruptCallback()
    
    # progressSet(value=progressStart + progressScale * (itor-1) / mixtures, message=sprintf("Step %d/%d: Running CIBERSORT", currStep, totalSteps), detail=sprintf("%s, please DO NOT refresh or close the page", time.eta))
    
    y <- Y[,itor]
    
    #standardize mixture
    y <- (y - mean(y)) / sd(y)
    
    #run SVR core algorithm
    result <- CoreAlg(X, y, absolute, abs_method)
    
    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    
    if(absolute && abs_method == 'sig.score') {
      w <- w * median(Y[,itor]) / Ymedian
    }
    
    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
    
    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(absolute) {out <- c(out, sum(w))}
    
    if (itor == 1) {
      output <- out
    } else {
      output <- rbind(output, out)
    }
    
    itor <- itor + 1
    
    end.time <- Sys.time()
    time.eta <- sprintf("ETA: %s", format((end.time - start.time) * (mixtures - itor + 1), format="%H:%M:%S"))
  }
  
  #save results
  #write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)
  
  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  if(!absolute) {
    colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
  } else {
    colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE",paste('Absolute score (',abs_method,')',sep=""))
  }
  obj
}

