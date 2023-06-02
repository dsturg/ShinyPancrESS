#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# I found a discussion of my python problem here:
# https://stackoverflow.com/questions/64221199/problems-with-deployment-of-shiny-app-due-to-reticulate-python

library('shiny')
library('dplyr')
library('Seurat')
library('dittoSeq')
#library('Rtsne')
#library('RColorBrewer')
# For color coding tSNE by gene
#library('mltools')
#library('BiocManager')
options(repos = BiocManager::repositories())

#~~~~~~~~~~~~
# Load data
# Note: This folder includes full size and reduced objects
#~~~~~~~~~~~~
#pan <- readRDS("pan_reduced.rds")
pan <- readRDS("tinypan.rds")

#load(file="ESS_cellXstudy_DESeqnorm_reduced.Rdata")
#ess <- ESS_cellXstudy_DESeqnorm_reduced
#rm(ESS_cellXstudy_DESeqnorm_reduced)

ess <- readRDS("ESS_cellXstudy_DESeqnorm_reduced.rds")
#ess <- ess_all_kallisto_counts
#ess$geneSymbol <- as.character(row.names(ess))

# Pull out gene names
gene_names <- ess$geneSymbol
gene_names <- sort(gene_names)

gene_names <- c("INS","PRSS1","ACTB")
#pan <- subset(x = pan,features = gene_names)
#ess <- ess[ess$geneSymbol %in% gene_names,]



#~~~~~~~~~~~~
#~~~~~~~~~~~~


# Use a fluid Bootstrap layout
ui <- fluidPage(    
  
  # Give the page a title
  titlePanel("Expression Specificity in Human Pancreas"),
  
  # Generate a row with a sidebar
  sidebarLayout(      
    
    # Define the sidebar with two input
    sidebarPanel(
      #selectInput("geneSymbol", "geneSymbol:", selected=NULL, multiple=FALSE, choices = character(0)),
      #selectInput("geneSymbol", "geneSymbol:", selected=NULL, multiple=FALSE, choices = character(0)),
      selectizeInput("geneSymbol", h4("geneSymbol"),
                     choices = "", selected = "", width = 230,
                     options = list(placeholder = "select a gene", maxItems = 1, maxOptions = 30, selectOnTab = T)),
    
      
      # See https://shiny.rstudio.com/reference/shiny/latest/radioButtons.html
      #radioButtons("rb", "Choose an ESS method:",
       #            choiceNames = list(
        #             "TPM",
        #             "CPM",
        #             "DESeq"
        #           ),
        #           choiceValues = list(
        #             "TPM", "CPM", "DES"
        #           )),
      #textOutput("txt"),
      
      
      hr(),
      helpText("Data from an integrated analysis of 4k\n
                     scRNA-Seq samples from published studies"),
      img(src = "logo.png", height = 200, width = 200)
    ),
    
# Consider a tab layout instead, looks pretty good:
# https://shiny.rstudio.com/articles/layout-guide.html

mainPanel("ESS and expression data",
          fluidRow(
            splitLayout(cellWidths = c("50%", "50%"), plotOutput("UMAP_celltype"),plotOutput("UMAP_gene"))
          ),
          plotOutput("ViolinPlot"),
          tableOutput("obs.ess"),
          tableOutput("obs.expr")
          
)



  ),
  # Try putting other plors on another row
  #plotOutput("ViolinPlot"),  
  #plotOutput("BarPlot"),  
  #tableOutput("obs")
  
  
)


# Define a server for the Shiny app
server <- function(input, output, session) {
  updateSelectizeInput(session, "geneSymbol", selected=NULL, choices = gene_names, server = TRUE)
  

  output$obs.ess <- renderTable({ 
  mytable.ess <- ess[ess$geneSymbol %in% input$geneSymbol,1:6]
  mytable.ess
  })
  output$obs.expr <- renderTable({ 
   mytable.expr <- ess[ess$geneSymbol %in% input$geneSymbol,7:12]
    mytable.expr
  })
  
  # Umap color coded by cell type (Not changed by selection)
#  output$UMAP_celltype <- plot(rnorm(21))
  output$UMAP_celltype <- renderPlot({   
    px <- dittoDimPlot(
      pan, "celltype", do.label = TRUE, labels.size = 3,
      main = "UMAP plot, by celltype",
      sub = "From Seurat normalized counts",
      legend.show = TRUE) +
      coord_fixed()
    px        
  })
# dittoplot Violins
  output$ViolinPlot <- renderPlot({
    if (is.null(input$geneSymbol)) 
      return(NULL)
    p3 <- dittoPlot(
      pan, input$geneSymbol, group.by = "celltype",
      plots = c("jitter", "vlnplot", "boxplot"),
      #boxplot.color = "white", boxplot.fill = FALSE,
      vlnplot.lineweight = 0.5,
      legend.show = FALSE,
      main = paste0(input$geneSymbol, " expression"),
      sub = "Units: Normalized counts",
      xlab = NULL)
    
    p3

  })
#Regular barplot
  output$BarPlot <- renderPlot({
    if (is.null(input$geneSymbol)) 
      return(NULL)
    # Render a barplot
    tobar <- unlist(ess[ess$geneSymbol %in% input$geneSymbol,1:5])
    barplot(tobar,
            main = paste0(input$geneSymbol, " ESS"),
            ylab="ESS",
            xlab="Cell type",las=1,ylim=c(0,1))
  })
  
# UMAP colored by gene
  output$UMAP_gene <- renderPlot({
    p1 <- dittoDimPlot(pan, input$geneSymbol,
                       max.color = "red", 
                       min.color = "gray70",
                       do.label = TRUE, labels.size = 3,
                       main = paste0("UMAP plot, ",input$geneSymbol),
                       sub = "Color by selected gene",
                       legend.show = TRUE) +
                        coord_fixed()
    p1
  })
  
}
# Run the application 
shinyApp(ui = ui, server = server)

