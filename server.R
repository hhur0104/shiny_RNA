# Shiny & DE
library(BiocManager)
library(edgeR)
library(limma)
library(genefilter)
library(shiny)
library(shinyjqui)
library(shinycssloaders)
library(plotly)
library(htmlwidgets)
library(DT)
library(heatmaply)
library(RColorBrewer)
library(ggVennDiagram)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(sva)
# Deconvolution
library(preprocessCore)
library(parallel)
library(e1071)
library(progress)
library(Metrics)
library(DTWBI)

options(shiny.maxRequestSize=30*1024^2)       # 30MB
options(repos = BiocManager::repositories())  # to use Bioconductor packages for shinyapps.io 

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

tpmNorm <- function(raw, len, genemap) {
  target <- raw 
  rownames(target) <- target$Geneid
  target <- target[rowSums(target[,-c(1:6)]) > 0,-c(1:6)]
  targetTPM <- target
  rownames(len) <- len$Geneid
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


# Server main
function(input, output, session) {
  values <- reactiveValues(target=NULL, 
                           design=NULL, 
                           tops=NULL, 
                           pcaData=NULL,
                           v=NULL
                           )
  
  ## counts
  observeEvent(input$file1, {
    # Values should flush
    target <- read.table(input$file1$datapath, header=TRUE)
    if (!identical(colnames(target)[1:6],c("Geneid", "Chr","Start","End","Strand","Length"))) {
      showModal(modalDialog(
        title = "Input Error",
          h4("First 6 columns must be gene metadata from featureCounts."),
          h4("ex) Geneid, Chr, Start, End, Strand, Length")
      ))
      return()
    }
    if (all(sapply(target[,-c(1:6)], is.numeric))!=TRUE) {
      pop_Error("Input Error", "Some count matrix columns include non-numeric values.")
      return()
    }
    values$raw <- target
    values$len <- target[,c("Geneid","Length")]
    
    
    rownames(target) <- target$Geneid
    target <- target[,-c(1:6)]
    target <- target[,c(order(colnames(target), decreasing=F))]
    
    values$target <- target
  })
  observeEvent(input$file1tip, {
    showModal(modalDialog(
      title = "Example Format of Counts.tsv",
      h4("- Count file must be in featureCount output format."),
      h4("- Expected coloumn order: Gene metadata followed by counts of each samples."),
      h4("- For each column names (sample) in 'counts', there must be one-to-one match to 'samplename' in condition.csv."),
      p(),
      img(src='count.png'),
      downloadButton("example", "Download Example Count Data", style='color: #444;')
    ))
  })

  ## design
  observeEvent(input$file2, {
    # Values should flush
    design <- read.csv(input$file2$datapath)
    rownames(design) <- design$samplename
    design <- design[c(order(rownames(design), decreasing=F)),]
  
    values$design <- design
   })
  observeEvent(input$file2tip, {
    showModal(modalDialog(
      title = "Example Format of Design file",
      h4("- This file consists sample metadata."),
      h4("- Column names must match the example below."),
      h4("- If there are more than one phenotype, concatenate conditions with a dot. "),
      h4("    ex) trt1.pheno1 or trt2.pheno2"),
      h4("- Optionally, provide Batch information for batch correction among samples."),
      p(),
      img(src='design.png'),
      downloadButton("example", "Download Example DataSet", style='color: #444;')
    ))
  })
  
  ## Run Differential Analysis
  observeEvent(input$runAnalysis, {
    req(input$file1)
    req(input$file2)
    req(input$pval)
    req(input$lgfch)
    req(input$comparisons)
    
    choices <- input$comparisons
    values$reference <- input$reference
    values$choices <- choices
    values$pval <- input$pval
    values$lgfch <- input$lgfch
    values$batch <- input$batch
    
    # design <- read.csv(input$file2$datapath)
    # rownames(design) <- design$samplename
    # design <- design[c(order(rownames(design), decreasing=F)),]
    # 
    # values$design <- design
    # 
    # target <- read.table(input$file1$datapath, header=TRUE)
    # rownames(target) <- target$Geneid
    # target <- target[,-c(1:6)]
    # target <- target[,c(order(colnames(target), decreasing=F))]
    # 
    # values$target <- target
    # top <- values$combo[rownames(values$combo) %in% choices,]
    # contrast.combo<-gsub(".vs.","-",choices)
    # Make All Possible Combo?
    contrastCombo <- gsub("\\.vs\\.","-",choices) 
    names(contrastCombo) <- choices
    target <- values$target
    design <- values$design
    
    # Validate
    if (!identical(colnames(target), rownames(design))) {
      pop_Error("Input Erorr","Column names in counts.tsv do not match with condition.csv")
      return() # Exit Function Here. Do nothing else.
    } 
    
    #### Tab container ####
    output$multabs <- renderUI({ 
      # Comparison Tab
      tabs <- lapply(choices, function(i) {
        ### Include more output here..
        tabPanel(
          i,
          fluidRow(
            box(
              title="Volcano Plot", width = 6, status="primary", solidHeader = TRUE,
              plotlyOutput(outputId=paste0("VC_",i)) ),#%>% withSpinner(color="#0dc5a1")),
            box(
              title="Gene Expression", width=6, status="primary", solidHeader = TRUE,
              plotlyOutput(outputId=paste0("GeneExp_",i)) )
          ),
          fluidRow(
            box(
              title="Top Genes", width=12, status="primary", solidHeader = TRUE,
              dataTableOutput(outputId=paste0("Top_",i)) %>% withSpinner(color="#0dc5a1"))
          ),
          fluidRow(
            box(
              title= "PCA (Pairwise)", status="primary", solidHeader = TRUE,
              plotlyOutput(outputId=paste0("PCA_pair_",i)) %>% withSpinner(color="#0dc5a1")),
            box(
              title="Heatmap (Pairwise)", status="primary", solidHeader = TRUE,
              plotlyOutput(outputId=paste0("HM_pair_",i)) %>% withSpinner(color="#0dc5a1"))
          ),
          fluidRow(
            box(
              title= "PCA (All Samples)", status="primary", solidHeader = TRUE,
              plotlyOutput(outputId=paste0("PCA_",i)) %>% withSpinner(color="#0dc5a1")),
            box(
              title="Heatmap (All Samples)", status="primary", solidHeader = TRUE,
              plotlyOutput(outputId=paste0("HM_",i)) %>% withSpinner(color="#0dc5a1"))
          )
        )
      })
      
      # VennDiagram Tab
      if (length(choices) > 1) {
        join <- data.frame(t(combn(unique(choices),2)))
        join <- paste(join$X1, " & ", join$X2, sep="")
        #print(join)
        tabs <- append(tabs, list(tabPanel(
            "Venn Diagram",
            fluidRow(
              box(
                title = "Select Genelists", width="6", status="primary", solidHeader=TRUE,
                # InputSelections
                radioButtons("vennComparison", "Choose Genelists to be joined:", choices = join)
              )
            ),
            fluidRow(
              box(
                title="VennDiagram of Differentially Expressed Genes ", width="12", status="primary", solidHeader=TRUE,
                plotOutput("venn.all")
              )
            ),
            fluidRow(
              box(
                title="Left Circle Only", width="4", status="primary", solidHeader=TRUE,
                dataTableOutput(outputId="vennList_L")
              ), 
              box(
                title="Intersection", width="4", status="primary", solidHeader=TRUE,
                dataTableOutput(outputId="vennList_C")
              ), 
              box(
                title="Right Circle Only", width="4", status="primary", solidHeader=TRUE,
                dataTableOutput(outputId="vennList_R")
              )
            )
        )))
      }
      
      # Insert Deconvolution Tab
      tabs <- append(tabs, list(tabPanel(
          "Cell Fractions",
          # Input reference
          fluidRow(
            box(
              title = "CIBERSORT Settings", width="4", status="primary", solidHeader=TRUE,
              selectInput("sigmtrx", "Choose a Reference Matrix:", choices = c("LM22","Derm22")),
              actionButton("runDeconv", "Run CIBERSORT"),
            )
          ),
          # Output CFR
          fluidRow(
            box(
              title="Cell Fraction of All Samples", width="12", status="primary", solidHeader=TRUE,
              dataTableOutput(outputId="cfr_table")
            )
          )
      )))
      
      do.call(tabsetPanel, c(tabs, id="tabs"))
    })
  
    isolate(calc_toptables(contrastCombo, values))
    
    #### Render All Graphs Here ####
    lapply(choices, function(i) {
      
      tb <- as.data.frame(values$tops[i])
      colnames(tb) <- gsub(paste0(i,"\\."), "", colnames(tb))
      tb$ids <- rownames(tb)
      
      ## TopTable
      output[[paste0("Top_",i)]] <- renderDataTable({
          datatable(
            tb[,c("genes","logFC","AveExpr","P.Value","adj.P.Val")],
            selection=list(mode = "single", selected = c(1)),
            options=list(
              rowCallback = JS(c(
                  "function(row, data, displayNum, index){",
                  "  var x2 = data[2];", # logFC
                  "  var x3 = data[3];", # AveExpr
                  "  var x4 = data[4];", # P.Value
                  "  var x5 = data[5];", # adj.P.Val
                  "  $('td:eq(2)', row).html(x2.toFixed(3));",
                  "  $('td:eq(3)', row).html(x3.toExponential(3));",
                  "  $('td:eq(4)', row).html(x4.toExponential(3));",
                  "  $('td:eq(5)', row).html(x5.toExponential(3));",
                  "}"))
              )
            ) 
      }) 
        
      observeEvent(input[[paste0("Top_",i,"_rows_selected")]], {
        info <- input[[paste0("Top_",i,"_rows_selected")]]
        vE <- t(t(values$v$E[rownames(tb[info,]),]))
        # Important line. Make Sure it works. # Need to assert identical(colnames(target),design$samplename)
        dt <- cbind(vE[order(rownames(vE)),], values$design[order(values$design$samplename),]) 
        rownames(dt) <- dt$samplename
        colnames(dt)[1] <- "Value"
        dt <- dt[dt$condition %in% unlist(strsplit(i, ".vs.")),]
        
        m <- tb[info,]
        a <- list (
          x = m$logFC,
          y = m$log10q,
          text = m$ids,
          xref = "x",
          yref = "y",
          showarrow = TRUE,
          arrowhead = 3,
          ax = 20,
          ay = -40
        )
        
        output[[paste0("VC_",i)]] <- renderPlotly({
          p <- plot_ly(tb, 
                  x=~logFC, 
                  y=~log10q, 
                  color=~color,
                  colors=c("Up.Significant"="red","Down.Significant"="blue","Not.Significant"="grey"),
                  type='scatter',
                  mode = 'markers',
                  marker=list(size=5),
                  text= ~ids,
                  hovertemplate = ~paste(
                    "<b>%{text}</b><br>",
                    "GeneSymbol:",genes,"<br>",
                    "logFC: %{x}<br>",
                    "Adj.P: ",adj.P.Val,"<br>",
                    "<extra></extra>")
          )
          p <- p %>% layout(annotations = a,
                            legend=list(orientation='h', y = -0.2))
          p
        })
        
        output[[paste0("GeneExp_",i)]] <- renderPlotly({
          # Gene Expression Per Group
          plot_ly(dt, 
                  x=~condition, 
                  y=~Value, 
                  color=~condition, 
                  type='scatter',
                  mode = 'markers',
                  marker=list(size=15),
                  text= ~samplename,
                  hovertemplate = ~paste(
                    "<b>%{text}</b><br>",
                    "NormalizedCount: %{y}<br>",
                    "<extra></extra>")
          ) %>% layout(title = m$ids)
        }) 
      })
      
      ## PCA
      output[[paste0("PCA_",i)]] <- renderPlotly({
        # Axis 
        lx <- list(title=paste0("PC1 ",round(values$pcaData[1,"var.explained"]*100)," % variance"), 
                   zeroline = FALSE, showticklabels = TRUE, showgrid = TRUE)
        ly <- list(title=paste0("PC2 ",round(values$pcaData[2,"var.explained"]*100)," % variance"), 
                   zeroline = FALSE, showticklabels = TRUE, showgrid = TRUE)
        # Plot
        p1 <- plot_ly(as.data.frame(values$pcaData), x = ~PC1, y=~PC2, 
                      color=~condition, type='scatter', marker=list(size=10),
                      text=~rownames(values$pcaData),
                      hovertemplate = ~paste(
                        "<b>%{text}</b><br>",
                        "PC1: %{x}<br>",
                        "PC2: %{y}<br>",
                        "<extra></extra>"))
      
        p1 <- p1 %>% layout(xaxis=lx,
                            yaxis=ly)
        
        # Return
        p1
      })
      
      ## PCA Pairwise
      output[[paste0("PCA_pair_",i)]] <- renderPlotly({
        pcatbl <- values$pairPCA[i][[1]]
        # Axis 
        lx <- list(title=paste0("PC1 ",round(pcatbl[1,"var.explained"]*100)," % variance"), 
                   zeroline = FALSE, showticklabels = TRUE, showgrid = TRUE)
        ly <- list(title=paste0("PC2 ",round(pcatbl[2,"var.explained"]*100)," % variance"), 
                   zeroline = FALSE, showticklabels = TRUE, showgrid = TRUE)
        # Plot
        p1 <- plot_ly(as.data.frame(pcatbl), x = ~PC1, y=~PC2, 
                      color=~condition, type='scatter', marker=list(size=10),
                      text=~rownames(pcatbl),
                      hovertemplate = ~paste(
                        "<b>%{text}</b><br>",
                        "PC1: %{x}<br>",
                        "PC2: %{y}<br>",
                        "<extra></extra>"))
        
        p1 <- p1 %>% layout(xaxis=lx,
                            yaxis=ly)
        
        # Return
        p1
      })
      
      ## Heatmap
      output[[paste0("HM_",i)]] <- renderPlotly({
        DEGs <- rownames(tb[tb$adj.P.Val <= values$pval & abs(tb$logFC) >= values$lgfch,])
        
        p1 <- heatmaply(
          values$v$E[DEGs,],
          showticklabels=c(TRUE, FALSE),
          scale="row",
          color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
          Rowv = TRUE,
          Colv = TRUE
        )
        
        p1
      })
    
      ## Heatmap Pairwise
      output[[paste0("HM_pair_",i)]] <- renderPlotly({
        DEGs <- rownames(tb[tb$adj.P.Val <= values$pval & abs(tb$logFC) >= values$lgfch,])
        sel <- unlist(strsplit(i, split = ".vs."))
        selSmp <- rownames(design[design$condition %in% sel,])
        
        p1 <- heatmaply(
          values$v$E[DEGs,selSmp],
          showticklabels=c(TRUE, FALSE),
          scale="row",
          color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
          Rowv = TRUE,
          Colv = TRUE
        )
        
        p1
      })
    })
    # Render VennDiagram
    if (length(choices) > 1 ) {
      output$venn.all <- renderPlot({
        sel <- unlist(strsplit(input$vennComparison, split = " & "))
        
        s1 <- values$tops[sel[1]][[1]]
        s2 <- values$tops[sel[2]][[1]]
        
        s1.deg <- s1[abs(s1$logFC) >= values$lgfch & s1$adj.P.Val <=values$pval, ]
        s2.deg <- s2[abs(s2$logFC) >= values$lgfch & s2$adj.P.Val <=values$pval, ]
        
        values$l <- s1.deg[!rownames(s1.deg) %in% rownames(s2.deg),]
        values$c <- s1.deg[ rownames(s1.deg) %in% rownames(s2.deg),]
        values$c2 <- s2.deg[rownames(s2.deg) %in% rownames(s1.deg),]
        values$r <- s2.deg[!rownames(s2.deg) %in% rownames(s1.deg),]
        
        x <- list(rownames(s1.deg), rownames(s2.deg))
        names(x) <- sel
        
        ggVennDiagram(x)
      })
      
      # Make Tables for Genelist
      output$vennList_L <- renderDataTable({
        datatable(
          values$l[,c("genes","logFC","adj.P.Val")],
          options=list(
            scrollX = TRUE,
            rowCallback = JS(c(
              "function(row, data, displayNum, index){",
              "  var x2 = data[2];", # logFC
              "  var x3 = data[3];", # AveExpr
              "  $('td:eq(2)', row).html(x2.toFixed(3));",
              "  $('td:eq(3)', row).html(x3.toExponential(3));",
              "}"))
          )
        ) 
      }) 
      output$vennList_C <- renderDataTable({
        dt <- merge(values$c[,c("genes","logFC","adj.P.Val")], values$c2[,c("logFC","adj.P.Val")], by="row.names")
     
        colnames(dt) <- gsub(".x",".Left",colnames(dt))
        colnames(dt) <- gsub(".y",".Right",colnames(dt))
        
        datatable(
          dt,
          options=list(
            scrollX = TRUE,
            rowCallback = JS(c(
              "function(row, data, displayNum, index){",
              "  var x3 = data[3];", # logFC
              "  var x4 = data[4];", # adjPVal
              "  var x5 = data[5];", # logFC
              "  var x6 = data[6];", # adjPVal
              "  $('td:eq(3)', row).html(x3.toFixed(3));",
              "  $('td:eq(4)', row).html(x4.toExponential(3));",
              "  $('td:eq(5)', row).html(x5.toFixed(3));",
              "  $('td:eq(6)', row).html(x6.toExponential(3));",
              "}"))
          )
        ) 
      })
      output$vennList_R <- renderDataTable({
        datatable(
          values$r[,c("genes","logFC","adj.P.Val")],
          options=list(
            scrollX = TRUE,
            rowCallback = JS(c(
              "function(row, data, displayNum, index){",
              "  var x2 = data[2];", # logFC
              "  var x3 = data[3];", # AveExpr
              "  $('td:eq(2)', row).html(x2.toFixed(3));",
              "  $('td:eq(3)', row).html(x3.toExponential(3));",
              "}"))
          )
        ) 
      }) 
    }
    
    # Render Export Button
    output$export <- renderUI({
      downloadButton("download", "Download Plots & Tables", style='color: #444;')
    })
  })
  
  ## Run Deconvolution
  observeEvent(input$runDeconv, {
    # TPM Normalize on Input file1
    mixture <- tpmNorm(values$raw, values$len, values$genemap)
    
    # Run CIBERSORT
    if (input$sigmtrx == "LM22") {
      sig_matrix <- read.table("reference/sig_LM22.txt",header=T,sep="\t",row.names=1,check.names=F)
    } else if (input$sigmtrx == "Derm22") {
      sig_matrix <- read.table("reference/sig_derm22.txt",header=T,sep="\t",row.names=1,check.names=F)
    }
    
    result <- CIBERSORT(sig_matrix, mixture, QN=FALSE)
    values$cfr_table <- result
    result <- round(result, 3)
    
    output$cfr_table <- renderDataTable({
      datatable(
        result,
        options=list(
          scrollX = TRUE #,
          # rowCallback = JS(c(
          #   "function(row, data, displayNum, index){",
          #   "  var x2 = data[2];", # logFC
          #   "  var x3 = data[3];", # AveExpr
          #   "  $('td:eq(2)', row).html(x2.toFixed(3));",
          #   "  $('td:eq(3)', row).html(x3.toExponential(3));",
          #   "}"))
        ))
    })
    
    
  })
  
  ## Download Handlers
  output$example <- downloadHandler(
    filename = "RNA-DifferentialAnalysis-ExampleDataSet.zip",
    content = function(file) {
      
      zip(zipfile=file, files=c("ex_condition.csv","ex_counts.tsv"))
      
    },
    contentType = "application/zip"
  )
  
  
  ## This DownloadHandler has potential to disrupt the server process when wd is not correctly returned to prior.
  ## Needs an exception handling. For Gracious Failure.
  output$download <- downloadHandler(
    filename = function () {
      paste0("RNA-DifferentialAnalysis-",Sys.Date(),".zip")
    },
    content = function(file) {
      fs <- c()
      deg <- c()
      wd <- getwd()
      setwd(tempdir())
      
      # CSV & Heatmap & Volcano & PCA _ PairWise
      for (i in values$choices) {
        ## BigTable_Pairwise
        bt <- values$tops[i][[1]][,c("genes"	,"logFC",	"AveExpr",	"t",	"P.Value",	"adj.P.Val",	"B")]
            csv <- paste0("BigTable_",i,".csv")
            write.csv(bt, csv)
        
        ## Heatmap_Pairwise
        htm <- paste0("Heatmap_",i,".png")
            DEGs <- rownames(bt[bt$adj.P.Val <= values$pval & abs(bt$logFC) >= values$lgfch,])
            sel <- unlist(strsplit(i, split = ".vs."))
            selSmp <- rownames(values$design[values$design$condition %in% sel,])
            
            df <- as.data.frame(values$design[selSmp,"condition"])
            colnames(df) <- "Condition"
            rownames(df) <- rownames(values$design[selSmp,])
            color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)
            
            hm <- pheatmap(values$v$E[DEGs,selSmp], show_rownames = F, scale = "row", 
                           col=color, annotation_col = df, main=paste0("DEGs in ",i))
            ggsave(htm, hm, width=4, height=4, dpi=200, units="in", device="png")
        
        ## VolcanoPlot
        vcp <- paste0("VolcanoPlot_",i,".png")
        
        ## PCA_Pairwise
        pca_name <- paste0("PCA_Plot_",i,".png")
            pcatbl <- values$pairPCA[i][[1]]
            p <- ggplot(pcatbl, aes(PC1, PC2, color = condition)) + 
              geom_label(data=pcatbl, aes(x=PC1, y=PC2, label=rownames(pcatbl)),  size=4) +
              labs(title = "PCA") +
              xlab(paste("PC1",round(pcatbl[1,"var.explained"] * 100),"% variance")) +
              ylab(paste("PC2",round(pcatbl[2,"var.explained"] * 100),"% variance")) +
              theme(plot.title=element_text(size=14),
                    axis.title.x=element_text(size=10),
                    axis.title.y=element_text(size=10))
            
          ggsave(pca_name, p, width=6, height=4, dpi=200, units="in", device="png")
        
        fs <- c(fs, csv, htm, pca_name)
        
        deg <- unique(c(deg,rownames(bt[abs(bt$logFC) >= values$lgfch & bt$adj.P.Val <= values$pval, ])))
      }
      
      
      # For all vennComparison choices 
      if (length(values$choices) > 1 ) {
        vennFs <- c()
        join <- data.frame(t(combn(unique(values$choices),2)))
        join <- paste(join$X1, ".&.", join$X2, sep="")
        for (i in join) {
          sel <- unlist(strsplit(i, split = ".&."))
          
          s1 <- values$tops[sel[1]][[1]]
          s2 <- values$tops[sel[2]][[1]]
          
          s1.deg <- s1[abs(s1$logFC) >= values$lgfch & s1$adj.P.Val <=values$pval, ]
          s2.deg <- s2[abs(s2$logFC) >= values$lgfch & s2$adj.P.Val <=values$pval, ]
          
          values$l <- s1.deg[!rownames(s1.deg) %in% rownames(s2.deg),]
          values$c <- s1.deg[ rownames(s1.deg) %in% rownames(s2.deg),]
          values$c2 <- s2.deg[rownames(s2.deg) %in% rownames(s1.deg),]
          values$r <- s2.deg[!rownames(s2.deg) %in% rownames(s1.deg),]
          
          cc2 <- merge(values$c[,c("genes","logFC","adj.P.Val")], values$c2[,c("logFC","adj.P.Val")], by="row.names")
          colnames(cc2) <- gsub(".x",".Left",colnames(cc2))
          colnames(cc2) <- gsub(".y",".Right",colnames(cc2))
          
          x <- list(rownames(s1.deg), rownames(s2.deg))
          names(x) <- sel
          
          csv_l <- paste0("VennDiagram_",join,"_Genelist_LeftCircle.csv")
          csv_r <- paste0("VennDiagram_",join,"_Genelist_RightCircle.csv")
          csv_c <- paste0("VennDiagram_",join,"_Genelist_CenterCircle.csv")
          write.csv(values$l[,c("genes","logFC","adj.P.Val")], csv_l)
          write.csv(values$r[,c("genes","logFC","adj.P.Val")], csv_r)
          write.csv(cc2, csv_c)
          
          vennP <- ggVennDiagram(x)
          vennName <- paste0("VennDiagram_",join,".png")
          ggsave(vennName,vennP,width=6, height=4, dpi=200, units="in", device="png")
          vennFs <- c(vennFs, vennName, csv_l, csv_r, csv_c)
        }
      }
      
      # PCA All Samples
      pca <- "PCA.png"
      pcaData.l <- values$pcaData
      p <- ggplot(pcaData.l, aes(PC1, PC2, color = condition)) + 
        geom_label(data=pcaData.l, aes(x=PC1, y=PC2, label=rownames(pcaData.l)),  size=4) +
        labs(title = "PCA") +
        xlab(paste("PC1",round(pcaData.l[1,"var.explained"] * 100),"% variance")) +
        ylab(paste("PC2",round(pcaData.l[2,"var.explained"] * 100),"% variance")) +
        theme(plot.title=element_text(size=14),
              axis.title.x=element_text(size=10),
              axis.title.y=element_text(size=10))
      ggsave(pca, p, width=6, height=4, dpi=200, units="in", device="png")
      
      # Heatmap All Samples
      hfile <- "Heatmap_All.DEGs.png"
      df <- as.data.frame(values$design[,"condition"])
      colnames(df) <- "Condition"
      rownames(df) <- rownames(values$design)
      color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)
      
      hm <- pheatmap(values$v$E[deg,], show_rownames = F, scale = "row", 
               col=color, annotation_col = df, main="Union of DEGs in each comparisons.")
      ggsave(hfile, hm, width=4, height=4, dpi=200, units="in", device="png")
      
      
      # Cell Fraction Result
      if (!is.null(values$cfr_table)) {
        cfr <- paste0("CellFractionResult.",input$sigmtrx,".csv")
        write.csv(values$cfr_table, cfr)
        fs <- c(fs, cfr)
      }
      
      # Venn Diagram
      if (length(values$choices) > 1) {
        zip(zipfile=file, files=c(fs, pca, hfile, vennFs))  
      } else {
        zip(zipfile=file, files=c(fs, pca, hfile))
      }
      
      setwd(wd)
    },
    contentType = "application/zip"
  )
  # Get current tab w/ following:
  # print(input$tabs)
  
  output$comparisonBoxes <- renderUI({
    req(input$file2)
    
    design <- read.csv(input$file2$datapath)
    rownames(design) <- design$samplename
    design <- design[c(order(rownames(design), decreasing=F)),]
    
    values$design <- design
    
    combo <- data.frame(t(combn(unique(values$design$condition),2)))
    combo$pair <- paste(combo$X1, ".vs.", combo$X2, sep="")
    rownames(combo) <- combo$pair
    #values$combo <- combo
    
    checkboxGroupInput("comparisons", "Choose Comparisons to be made:", combo$pair)
  })
  
  output$batchCrrct <- renderUI({
    req(input$file2)
    
    design <- read.csv(input$file2$datapath)
    rownames(design) <- design$samplename
    design <- design[c(order(rownames(design), decreasing=F)),]
    
    cs <- c("No Batch Correction needed.", colnames(design)[-c(1:3)])
    radioButtons("batch","Column for Batch Correction (optional):",
                 choices=cs,
                 selected="No Batch Correction needed.")
  })
}