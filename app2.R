# RNA-Seq Analysis App
# written for CDSBF 591-R
# code written by: Bella Pfeiffer
# date: 17 December 2024

# loading all the packages we need for data analysis and visualization
library(shiny)
library(shinythemes)
library(DT)               # for nice interactive tables
library(tidyverse)        # for data wrangling
library(pheatmap)         # for making heatmaps
library(plotly)           # for interactive plots
library(reshape2)         # for reshaping data
library(RColorBrewer)     # for better color palettes
library(matrixStats)      # for matrix operations
library(Matrix)           # for handling sparse matrices
library(fgsea)           # for gene set enrichment analysis
library(igraph)          # for network analysis
library(ggraph)          # for network visualization

# helper functions section 

# this function reads in our data files and makes sure they have the columns we need
# i used this a lot for inputs + error checking
read_input_file <- function(file_path, req_cols = NULL) {
  data <- read.csv(file_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  
  tryCatch({
    # check if we're missing any columns we absolutely need
    if (!is.null(req_cols)) {
      missing_cols <- setdiff(req_cols, colnames(data))
      if (length(missing_cols) > 0) {
        return(list(
          success = FALSE,
          message = paste("looks like we're missing these columns:", paste(missing_cols, collapse = ", "))
        ))
      }
    }
    
    return(list(success = TRUE, data = data))
  }, error = function(e) {
    return(list(
      success = FALSE,
      message = paste("oops, had trouble reading the file:", e$message)
    ))
  })
}

# this function gives us a nice summary of our dataset
# tells us things like how many rows/columns we have and what kind of data is in each column
generate_summary <- function(data) {
  if (is.null(data) || nrow(data) == 0) {
    return("no data to show yet")
  }
  
  summary_text <- paste("here's what your dataset looks like:\n",
                        "number of rows:", nrow(data), "\n",
                        "number of columns:", ncol(data), "\n\n",
                        "let's look at each column:\n")
  
  for (col in names(data[-1])) {
    summary_text <- paste0(summary_text, "\n", col, ":\n")
    summary_text <- paste0(summary_text, "type: ", class(data[[col]])[1], "\n")
    
    # for numeric columns, show the range and missing values
    if (is.numeric(data[[col]])) {
      summary_text <- paste0(summary_text, 
                             "range: [", min(data[[col]], na.rm = TRUE), ", ",
                             max(data[[col]], na.rm = TRUE), "]\n",
                             "missing values: ", sum(is.na(data[[col]])), "\n")
    } else {
      # for non-numeric columns, show unique values and missing values
      summary_text <- paste0(summary_text,
                             "unique values: ", length(unique(data[[col]])), "\n",
                             "missing values: ", sum(is.na(data[[col]])), "\n")
    }
  }
  
  return(summary_text)
}

# this function reads in our GMT files for pathway analysis
# GMT files contain our gene sets - each line is a pathway and its genes
read_gmt <- function(gmt_file) {
  con <- file(gmt_file, "r")
  lines <- readLines(con)
  close(con)
  
  # create a list where each element is a pathway and its genes
  pathways <- list()
  for (line in lines) {
    parts <- unlist(strsplit(line, "\t"))
    pathway_name <- parts[1]
    # we skip the description (it's in parts[2]) and just get the genes
    genes <- parts[3:length(parts)]
    pathways[[pathway_name]] <- genes
  }
  return(pathways)
}

# now let's build our UI
ui <- fluidPage(
  theme = shinytheme("flatly"),
  
  # main title and subtitle
  titlePanel("RNA-Seq Data Analysis Dashboard"),
  h4("interactive analysis and visualization of RNA sequencing data", style = "color: #666; margin-bottom: 20px;"),
  
  # organizing everything into tabs for better navigation
  tabsetPanel(
    # first tab: looking at our sample information
    tabPanel("Sample Information",
             sidebarLayout(
               sidebarPanel(
                 # let users upload their sample info
                 fileInput("sampleFile", "upload your sample information (CSV)",
                           accept = c("text/csv", 
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 actionButton("analyzeSamples", "analyze samples",
                              class = "btn-primary"),
                 width = 3
               ),
               mainPanel(
                 tabsetPanel(
                   # show summary stats
                   tabPanel("Summary", 
                            verbatimTextOutput("sampleSummary")),
                   # show the actual data
                   tabPanel("Data Table", 
                            DTOutput("sampleTable")),
                   # visualize the data
                   tabPanel("Visualization",
                            fluidRow(
                              column(4,
                                     selectInput("plotColumn", "which variable should we plot?",
                                                 choices = NULL)),
                              column(4,
                                     selectInput("plotType", "how should we visualize it?",
                                                 choices = c("Histogram", "Density", "Box"),
                                                 selected = "Histogram"))
                            ),
                            plotlyOutput("samplePlot", height = "600px"))
                 )
               )
             )
    ),
    # second tab: analyzing expression data
    tabPanel("Expression Analysis",
             sidebarLayout(
               sidebarPanel(
                 # let users upload their normalized count data
                 fileInput("countsFile", "upload your normalized counts (CSV)",
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 # sliders to filter our data
                 sliderInput("variancePercentile", 
                             "how variable should genes be to include them?",
                             min = 0, max = 100, value = 50),
                 sliderInput("nonZeroSamples",
                             "what % of samples should have non-zero expression?",
                             min = 0, max = 100, value = 50),
                 actionButton("analyzeExpression", "analyze expression",
                              class = "btn-primary"),
                 width = 3
               ),
               mainPanel(
                 # organizing our analysis outputs
                 tabsetPanel(
                   tabPanel("Summary",
                            verbatimTextOutput("expressionSummary")),
                   tabPanel("Diagnostic Plots",
                            fluidRow(
                              column(12,
                                     plotlyOutput("diagExpression", height = "400px"),
                                     br(),
                                     plotlyOutput("diagDetection", height = "400px")
                              )
                            )
                   ),
                   # heatmap to visualize expression patterns
                   tabPanel("Heatmap",
                            plotOutput("expressionHeatmap", height = "800px")),
                   # pca plot to look at sample relationships
                   tabPanel("PCA",
                            fluidRow(
                              column(3,
                                     selectInput("pcX", "which PC on x-axis?", choices = NULL),
                                     selectInput("pcY", "which PC on y-axis?", choices = NULL)),
                              column(9,
                                     plotlyOutput("pcaPlot", height = "600px"))
                            ))
                 )
               )
             )
    ),
    
    # third tab: differential expression analysis
    tabPanel("Differential Expression",
             sidebarLayout(
               sidebarPanel(
                 # upload results from differential expression analysis
                 fileInput("deFile", "upload your DE results file",
                           accept = c(
                             "text/csv",
                             "text/comma-separated-values",
                             "text/tab-separated-values",
                             "text/plain",
                             ".csv",
                             ".tsv",
                             ".txt"
                           )
                 ),
                 # search for specific genes you're interested in
                 textInput("geneSearch", "looking for specific genes?", ""),
                 # set cutoffs for significance
                 numericInput("pvalCutoff", 
                              "how significant should genes be? (p-value cutoff)",
                              value = 0.05,
                              min = 0,
                              max = 1,
                              step = 0.01
                 ),
                 numericInput("fcCutoff",
                              "how much should expression change? (log2 fold change cutoff)",
                              value = 1,
                              min = 0,
                              step = 0.1
                 ),
                 width = 3
               ),
               mainPanel(
                 tabsetPanel(
                   # interactive table of results
                   tabPanel("Results Table", 
                            DTOutput("deTable")
                   ),
                   # volcano plot to visualize changes
                   tabPanel("Volcano Plot",
                            plotlyOutput("volcanoPlot", height = "600px")
                   )
                 )
               )
             )
    ),
    
    # fourth tab: gene set enrichment analysis
    tabPanel("Gene Set Enrichment",
             sidebarLayout(
               sidebarPanel(
                 # file inputs for GSEA
                 fileInput("rankedGenesFile", "upload your ranked genes (CSV)",
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 fileInput("gmtFile", "upload your gene sets (GMT)",
                           accept = c("text/plain",
                                      ".gmt",
                                      ".txt")),
                 
                 # parameters for the analysis
                 numericInput("minSize", 
                              "what's the smallest gene set to consider?",
                              value = 15,
                              min = 5,
                              max = 100),
                 numericInput("maxSize",
                              "what's the largest gene set to consider?",
                              value = 500,
                              min = 50,
                              max = 1000),
                 
                 # visualization controls
                 sliderInput("pvalFilterSlider",
                             "how significant should pathways be?",
                             min = -20,
                             max = 0,
                             value = -2,
                             step = 1),
                 
                 # buttons to run analysis and download results
                 actionButton("runGSEA", "run GSEA analysis",
                              class = "btn-primary"),
                 downloadButton("downloadGSEA", "download results"),
                 width = 3
               ),
               
               mainPanel(
                 tabsetPanel(
                   # table of pathway results
                   tabPanel("Results Table",
                            DTOutput("gseaResultsTable")),
                   # overview of enrichment scores
                   tabPanel("Pathway Overview",
                            plotOutput("gseaNESPlot", height = "600px")),
                   # detailed enrichment plots
                   tabPanel("Enrichment Plots",
                            selectInput("pathwaySelect", "which pathway should we look at?",
                                        choices = NULL),
                            plotOutput("enrichmentPlot", height = "400px"),
                            verbatimTextOutput("pathwayStats"))
                 )
               )
             ))
  )
)

# now for the server logic - this is where all the analysis happens
server <- function(input, output, session) {
  # let's handle bigger files (up to 50MB)
  options(shiny.maxRequestSize = 50 * 1024^2)
  
  # setting up reactive values to store our data
  values <- reactiveValues(
    sample_data = NULL,     # for sample information
    counts_data = NULL,     # for expression counts
    filtered_counts = NULL, # for filtered expression data
    pca_data = NULL,       # for PCA results
    de_data = NULL,        # for differential expression results
    gsea = reactiveValues(  # for GSEA analysis
      ranked_genes = NULL,
      pathways = NULL,
      results = NULL
    )
  )
  
  # handling sample data uploads and processing
  observeEvent(input$sampleFile, {
    req(input$sampleFile)
    result <- read_input_file(input$sampleFile$datapath)
    
    if (result$success) {
      # store the data if upload worked
      values$sample_data <- result$data
      
      # update our plot controls to show available numeric columns
      updateSelectInput(session, "plotColumn",
                        choices = names(values$sample_data)[sapply(values$sample_data, is.numeric)])
    } else {
      # let the user know if something went wrong
      showNotification(result$message, type = "error")
    }
  })
  
  # create our sample data summary
  output$sampleSummary <- renderPrint({
    req(values$sample_data)
    cat(generate_summary(values$sample_data))
  })
  
  # create an interactive table of sample data
  output$sampleTable <- renderDT({
    req(values$sample_data)
    datatable(values$sample_data,
              options = list(pageLength = 25,
                             scrollX = TRUE,
                             dom = 'Bfrtip',
                             buttons = c('copy', 'csv', 'excel')))
  })
  
  # make plots of sample data
  output$samplePlot <- renderPlotly({
    req(values$sample_data, input$plotColumn)
    
    data <- values$sample_data
    
    # create different types of plots based on user selection
    p <- switch(input$plotType,
                "Histogram" = ggplot(data, aes_string(x = input$plotColumn)) +
                  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7),
                "Density" = ggplot(data, aes_string(x = input$plotColumn)) +
                  geom_density(fill = "steelblue", alpha = 0.7),
                "Box" = ggplot(data, aes_string(y = input$plotColumn)) +
                  geom_boxplot(fill = "steelblue", alpha = 0.7))
    
    # add some nice formatting
    p <- p + theme_minimal() +
      labs(title = paste("looking at the distribution of", input$plotColumn))
    
    ggplotly(p)
  })
  
  # handle expression data uploads
  observeEvent(input$countsFile, {
    req(input$countsFile)
    result <- read_input_file(input$countsFile$datapath)
    
    if (result$success) {
      values$counts_data <- result$data
      
      # update PCA plot controls - we'll show up to 10 PCs
      updateSelectInput(session, "pcX", 
                        choices = paste0("PC", 1:min(10, ncol(values$counts_data)-1)),
                        selected = "PC1")
      updateSelectInput(session, "pcY",
                        choices = paste0("PC", 1:min(10, ncol(values$counts_data)-1)),
                        selected = "PC2")
    } else {
      showNotification(result$message, type = "error")
    }
  })
  
  # filter our expression data based on user settings
  observe({
    req(values$counts_data, input$variancePercentile, input$nonZeroSamples)
    
    withProgress(message = 'filtering your expression data...', value = 0, {
      # get our numeric data (assuming first column is gene names)
      counts_matrix <- as.matrix(values$counts_data[,-1])
      gene_names <- values$counts_data[[1]]
      
      # calculate variance and expression frequency for each gene
      var_vec <- rowVars(counts_matrix)
      nonzero_prop <- rowMeans(counts_matrix > 0) * 100
      
      # apply our filters
      var_cutoff <- quantile(var_vec, input$variancePercentile/100)
      keep_genes <- var_vec >= var_cutoff & nonzero_prop >= input$nonZeroSamples
      
      # store filtered data and stats
      values$filtered_counts <- list(
        matrix = counts_matrix[keep_genes,],
        gene_names = gene_names[keep_genes],
        stats = list(
          total_genes = length(gene_names),
          passing_genes = sum(keep_genes),
          variance_cutoff = var_cutoff
        )
      )
      
      # calculate PCA if we have enough genes
      if (sum(keep_genes) > 0) {
        scaled_counts <- t(scale(t(values$filtered_counts$matrix)))
        values$pca_data <- prcomp(t(scaled_counts))
      } else {
        showNotification("no genes passed the current filtering criteria. 
                       try adjusting the variance percentile or non-zero sample threshold.",
                         type = "warning")
      }
    })
  })
  
  # show expression analysis summary
  output$expressionSummary <- renderPrint({
    req(values$filtered_counts)
    
    stats <- values$filtered_counts$stats
    cat("here's what we found in your expression data:\n\n")
    cat("total number of genes:", stats$total_genes, "\n")
    cat("genes that passed filtering:", stats$passing_genes, 
        sprintf("(%.1f%%)\n", stats$passing_genes/stats$total_genes*100))
    cat("\nfiltering criteria used:\n")
    cat("- variance percentile cutoff:", input$variancePercentile, "\n")
    cat("- minimum % non-zero samples:", input$nonZeroSamples, "\n")
  })
  
  # diagnostic plot for expression variance vs median
  output$diagExpression <- renderPlotly({
    req(values$counts_data)
    
    withProgress(message = 'generating expression diagnostic plot...', value = 0, {
      # Calculate metrics using the raw counts data
      counts_matrix <- as.matrix(values$counts_data[,-1]) # Remove gene names
      var_vec <- rowVars(counts_matrix, na.rm = TRUE)
      medians <- rowMedians(counts_matrix, na.rm = TRUE)
      
      # Calculate whether each gene passes our filters
      nonzero_prop <- rowMeans(counts_matrix > 0, na.rm = TRUE) * 100
      var_cutoff <- quantile(var_vec, input$variancePercentile/100, na.rm = TRUE)
      passed_filter <- var_vec >= var_cutoff & nonzero_prop >= input$nonZeroSamples
      
      data <- data.frame(
        variance = log2(var_vec + 1),
        median = log2(medians + 1),
        passed_filter = passed_filter
      )
      
      p <- ggplot(data, aes(x = median, y = variance, color = passed_filter)) +
        geom_point(alpha = 0.6) +
        scale_color_manual(values = c("grey70", "#0072B2")) +
        theme_minimal() +
        labs(x = "log2 median expression",
             y = "log2 variance",
             color = "passed filtering",
             title = "expression variance vs median") +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10)
        )
      
      ggplotly(p, tooltip = c("x", "y"))
    })
  })
  
  # diagnostic plot for detection rate vs expression
  output$diagDetection <- renderPlotly({
    req(values$counts_data)
    
    withProgress(message = 'generating detection rate plot...', value = 0, {
      # calculate metrics using the raw counts data
      counts_matrix <- as.matrix(values$counts_data[,-1]) # remove gene names
      var_vec <- rowVars(counts_matrix, na.rm = TRUE)
      medians <- rowMedians(counts_matrix, na.rm = TRUE)
      
      # calculate whether each gene passes our filters
      nonzero_prop <- rowMeans(counts_matrix > 0, na.rm = TRUE) * 100
      var_cutoff <- quantile(var_vec, input$variancePercentile/100, na.rm = TRUE)
      passed_filter <- var_vec >= var_cutoff & nonzero_prop >= input$nonZeroSamples
      
      data <- data.frame(
        nonzero = nonzero_prop,
        median = log2(medians + 1),
        passed_filter = passed_filter
      )
      
      p <- ggplot(data, aes(x = median, y = nonzero, color = passed_filter)) +
        geom_point(alpha = 0.6) +
        scale_color_manual(values = c("grey70", "#0072B2")) +
        theme_minimal() +
        labs(x = "log2 median expression",
             y = "% non-zero samples",
             color = "passed filtering",
             title = "detection rate vs expression level") +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10)
        )
      
      ggplotly(p, tooltip = c("x", "y"))
    })
  })
  
  # create our expression heatmap with better error handling
  output$expressionHeatmap <- renderPlot({
    req(values$filtered_counts)
    
    withProgress(message = 'generating heatmap...', value = 0, {
      tryCatch({
        # Validate we have data
        if (is.null(values$filtered_counts$matrix) || nrow(values$filtered_counts$matrix) == 0) {
          stop("no genes passed filtering criteria")
        }
        
        # get the most variable genes for visualization
        n_genes <- min(50, nrow(values$filtered_counts$matrix))
        var_vec <- rowVars(values$filtered_counts$matrix)
        top_idx <- order(var_vec, decreasing = TRUE)[1:n_genes]
        
        # make a nice heatmap
        pheatmap(log2(values$filtered_counts$matrix[top_idx,] + 1),
                 scale = "row",
                 show_rownames = FALSE,
                 main = paste("expression patterns of the top", n_genes, "most variable genes"),
                 color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
      }, error = function(e) {
        # If there's an error, show it to the user
        showNotification(paste("oops, had trouble making the heatmap:", e$message), 
                         type = "error")
        return(NULL)
      })
    })
  })
  
  # create our PCA plot
  output$pcaPlot <- renderPlotly({
    req(values$pca_data, input$pcX, input$pcY)
    
    # get the PCA data ready
    pc_data <- as.data.frame(values$pca_data$x)
    var_explained <- summary(values$pca_data)$importance[2,] * 100
    
    # create an interactive scatter plot
    p <- ggplot(pc_data, aes_string(x = input$pcX, y = input$pcY)) +
      geom_point(size = 3) +
      theme_minimal() +
      labs(x = sprintf("%s (explains %.1f%% of variance)", input$pcX, 
                       var_explained[as.numeric(gsub("PC", "", input$pcX))]),
           y = sprintf("%s (explains %.1f%% of variance)", input$pcY, 
                       var_explained[as.numeric(gsub("PC", "", input$pcY))]))
    
    ggplotly(p)
  })
  
  # handle differential expression file uploads
  observeEvent(input$deFile, {
    req(input$deFile)
    result <- read_input_file(input$deFile$datapath,
                              req_cols = c("log2FoldChange", "padj"))
    
    if (result$success) {
      values$de_data <- result$data
    } else {
      showNotification(result$message, type = "error")
    }
  })
  
  # create our interactive DE results table
  output$deTable <- renderDT({
    req(values$de_data)
    
    # start with all our data
    filtered_data <- values$de_data
    
    # filter for specific genes if requested
    if (nchar(input$geneSearch) > 0) {
      filtered_data <- filtered_data[grep(input$geneSearch, 
                                          filtered_data[[1]], 
                                          ignore.case = TRUE),]
    }
    
    # apply our significance filters
    filtered_data <- filtered_data[!is.na(filtered_data$padj) & 
                                     filtered_data$padj <= input$pvalCutoff &
                                     abs(filtered_data$log2FoldChange) >= input$fcCutoff,]
    
    # make an interactive table
    datatable(filtered_data,
              options = list(pageLength = 25,
                             scrollX = TRUE,
                             dom = 'Bfrtip',
                             buttons = c('copy', 'csv', 'excel'))) %>%
      formatSignif(columns = c("padj", "log2FoldChange"), digits = 3)
  })
  
  # create our volcano plot
  output$volcanoPlot <- renderPlotly({
    req(values$de_data)
    
    plot_data <- values$de_data
    
    # categorize our genes based on significance and fold change
    plot_data$sig <- "not significant"
    plot_data$sig[plot_data$padj <= input$pvalCutoff & 
                    plot_data$log2FoldChange >= input$fcCutoff] <- "upregulated"
    plot_data$sig[plot_data$padj <= input$pvalCutoff & 
                    plot_data$log2FoldChange <= -input$fcCutoff] <- "downregulated"
    
    # create an interactive volcano plot
    p <- ggplot(plot_data, aes(x = log2FoldChange, 
                               y = -log10(padj),
                               color = sig,
                               text = plot_data[,1])) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = c(
        "upregulated" = "red",
        "downregulated" = "blue",
        "not significant" = "grey50")) +
      theme_minimal() +
      labs(x = "log2 fold change",
           y = "-log10 adjusted p-value",
           color = "regulation status",
           title = "volcano plot of expression changes")
    
    ggplotly(p, tooltip = c("text", "x", "y"))
  })
  
  # handle GSEA data uploads and processing
  observeEvent(input$rankedGenesFile, {
    req(input$rankedGenesFile)
    result <- read_input_file(input$rankedGenesFile$datapath)
    
    if (result$success) {
      tryCatch({
        # prepare our ranked gene list
        ranked_genes <- result$data
        values$gsea$ranked_genes <- sort(setNames(ranked_genes$rank, ranked_genes$gene),
                                         decreasing = TRUE)
        
        # check if our ranking looks okay
        rank_range <- range(ranked_genes$rank)
        if (all(ranked_genes$rank >= 0) || all(ranked_genes$rank <= 0)) {
          showNotification("heads up: your ranking seems to be one-sided. 
                         you might want to use a metric that goes both positive and negative,
                         like log2 fold change or signal-to-noise ratio.", 
                           type = "warning")
        }
        
        showNotification("successfully loaded your ranked genes", type = "message")
      }, error = function(e) {
        showNotification(paste("oops, had trouble with the ranked genes:", e$message),
                         type = "error")
      })
    } else {
      showNotification(result$message, type = "error")
    }
  })
  
  # handle GMT file uploads
  observeEvent(input$gmtFile, {
    req(input$gmtFile)
    tryCatch({
      values$gsea$pathways <- read_gmt(input$gmtFile$datapath)
      showNotification("great! your gene sets are loaded", type = "message")
    }, error = function(e) {
      showNotification(paste("oops, had trouble with the GMT file:", e$message), 
                       type = "error")
    })
  })
  
  # run the GSEA analysis when the button is clicked
  observeEvent(input$runGSEA, {
    req(values$gsea$ranked_genes, values$gsea$pathways)
    
    withProgress(message = 'running GSEA analysis...', value = 0, {
      tryCatch({
        # run fgsea with some optimized parameters
        results <- fgsea::fgsea(
          pathways = values$gsea$pathways,
          stats = values$gsea$ranked_genes,
          minSize = input$minSize,
          maxSize = input$maxSize,
          nperm = 1000,    # more permutations = more accurate p-values
          scoreType = "std" # standardized scoring type
        )
        
        # store our results
        values$gsea$results <- as.data.frame(results)
        
        # check if our enrichment scores look balanced
        nes_summary <- summary(results$NES)
        if (all(results$NES > 0) || all(results$NES < 0)) {
          showNotification("hey, all your enrichment scores are going one way. 
                         might want to double-check your gene ranking method.", 
                           type = "warning")
        }
        
        # update the pathway selector with our results
        updateSelectInput(session, "pathwaySelect",
                          choices = values$gsea$results$pathway)
        
        showNotification("GSEA analysis completed successfully", type = "message")
      }, error = function(e) {
        showNotification(paste("oops, had trouble running GSEA:", e$message), 
                         type = "error")
      })
    })
  })
  
  # create our enrichment score overview plot
  output$gseaNESPlot <- renderPlot({
    req(values$gsea$results)
    
    withProgress(message = 'creating pathway overview plot...', value = 0, {
      tryCatch({
        # prepare our data for plotting
        plot_data <- values$gsea$results %>%
          filter(padj <= 10^input$pvalFilterSlider) %>%
          mutate(
            direction = ifelse(NES > 0, "upregulated", "downregulated"),
            pathway = factor(pathway, levels = pathway[order(NES)]),
            logPadj = -log10(padj)
          )
        
        if (nrow(plot_data) == 0) {
          stop("no pathways meet the current significance threshold")
        }
        
        # create a scatter plot of enrichment scores
        ggplot(plot_data, aes(x = NES, y = logPadj)) +
          geom_point(aes(color = direction), size = 3, alpha = 0.7) +
          scale_color_manual(values = c("downregulated" = "#0072B2", 
                                        "upregulated" = "#D55E00")) +
          theme_minimal() +
          labs(
            title = "overview of GSEA results",
            x = "normalized enrichment score (NES)",
            y = "-log10 adjusted p-value",
            color = "direction"
          )
      }, error = function(e) {
        showNotification(paste("oops, had trouble making the plot:", e$message), 
                         type = "error")
        return(NULL)
      })
    })
  })
  
  # create our results table
  output$gseaResultsTable <- renderDT({
    req(values$gsea$results)
    
    # filter by significance
    filtered_results <- values$gsea$results %>%
      filter(padj <= 10^input$pvalFilterSlider)
    
    if (nrow(filtered_results) == 0) {
      showNotification("no pathways meet the current significance threshold",
                       type = "warning")
      return(NULL)
    }
    
    # make an interactive table
    datatable(filtered_results,
              options = list(
                pageLength = 15,
                scrollX = TRUE,
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel')
              )
    ) %>%
      formatSignif(columns = c("pval", "padj", "NES"), digits = 3)
  })
  
  # create enrichment plots for individual pathways
  output$enrichmentPlot <- renderPlot({
    req(values$gsea$results, input$pathwaySelect)
    
    withProgress(message = 'creating enrichment plot...', value = 0, {
      tryCatch({
        # use fgsea's built-in plotting function
        plotEnrichment(values$gsea$pathways[[input$pathwaySelect]],
                       values$gsea$ranked_genes) +
          labs(title = paste("enrichment plot for", input$pathwaySelect))
      }, error = function(e) {
        showNotification(paste("oops, had trouble making the enrichment plot:", e$message), 
                         type = "error")
        return(NULL)
      })
    })
  })
  
  # show detailed stats for selected pathway
  output$pathwayStats <- renderPrint({
    req(values$gsea$results, input$pathwaySelect)
    
    pathway_data <- values$gsea$results %>%
      filter(pathway == input$pathwaySelect)
    
    cat("let's look at the details for this pathway:\n")
    cat("name:", pathway_data$pathway, "\n")
    cat("enrichment score (NES):", round(pathway_data$NES, 3), "\n")
    cat("p-value:", format.pval(pathway_data$pval, digits = 3), "\n")
    cat("adjusted p-value:", format.pval(pathway_data$padj, digits = 3), "\n")
    cat("number of genes in set:", pathway_data$size, "\n")
  })
  
  # set up the download handler for GSEA results
  output$downloadGSEA <- downloadHandler(
    filename = function() {
      paste0("gsea_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(values$gsea$results, file, row.names = FALSE)
    }
  )
}
        
# run the application!
shinyApp(ui = ui, server = server)