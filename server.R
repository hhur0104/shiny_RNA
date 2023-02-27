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
    load("./grch37.genemap.Rdata")
    idx <- match(genenames, genemap.1$ensembl_gene_id )
    genesymbol <- genemap.1$external_gene_name[idx]
  }
  if (values$reference == "Mouse") {
    load("mmusculus.genemap.Rdata")
    idx <- match(genenames, genemap.1$ensembl_gene_id )
    genesymbol <- genemap.1$mgi_symbol[idx]
   
  }
  if (values$reference == "Other") {
    genesymbol <- rep("No Information",nrow(target.final))
  }
  
  counts <- DGEList(counts=target.final,genes=genesymbol) # Create Limma data object
  isexpr<-rowSums(cpm(counts)>1)>=3
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
  
  # PCA data (EveryThing)
  pcaData <- prcomp(t(v$E))
  pcaData <- data.frame(cbind(pcaData$x[,c(1,2)], var.explained=pcaData$sdev^2/sum(pcaData$sdev^2)))
  pcaData <- cbind(pcaData, "condition"=design[rownames(pcaData), "condition"])
  values$pcaData <- pcaData
}

pop_Error <- function(title, msg) {
    showModal(modalDialog(
      title = title,
      msg
    ))
}

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
    
    if (colnames(target)[1] != "Geneid") {
      pop_Error("Input Error", "First column of the count matrix must be 'Geneid'.")
      return()
    }
    if (all(sapply(target[,-1], is.numeric))!=TRUE) {
      pop_Error("Input Error", "Some count matrix columns include non-numeric values.")
      return()
    }
    
    rownames(target) <- target$Geneid
    target <- target[,-1]
    target <- target[,c(order(colnames(target), decreasing=F))]
  
    values$target <- target
  })
  observeEvent(input$file1tip, {
    showModal(modalDialog(
      title = "Example Format of Counts.tsv",
      h4("- Do not include any gene metadata."),
      h4("- First column must be GeneID."),
      h4("- Subsequent columns must be counts of each samples."),
      p(),
      img(src='count.png'),
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
      p(),
      img(src='design.png'),
    ))
  })
  
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
      pop_Erorr("Input Erorr","Column names in counts.tsv do not match with condition.csv")
      return() # Exit Function Here. Do nothing else.
    } 
    
    #### Tab container ####
    output$multabs <- renderUI({ 
      
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
              title= "PCA", status="primary", solidHeader = TRUE,
              plotlyOutput(outputId=paste0("PCA_",i)) %>% withSpinner(color="#0dc5a1")),
            box(
              title="Heatmap", status="primary", solidHeader = TRUE,
              plotlyOutput(outputId=paste0("HM_",i)) %>% withSpinner(color="#0dc5a1"))
          )
        )
      })
      do.call(tabsetPanel, tabs)
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
    })
  })
  
  
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
  
}