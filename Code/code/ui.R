shinyUI(fluidPage(
  titlePanel("My Shiny App"),
  sidebarLayout(
    sidebarPanel(
      selectInput("Function", 
                  label = "Choose an analysis function",
                  choices = list("Phenotype geration", "GWAs","Locus zoom",
                                 "TWAs", "Omic QTL","two trait MR","Omic SMR"),
                  selected = "GWAs"),
      
      
      conditionalPanel(condition="input.Function == 'Phenotype geration'",
                       fileInput(
                         inputId = "phenodata",
                         label = "Upload Phenotype Data File",
                         accept = c(".txt", ".csv"))
      ),
        
      conditionalPanel(condition="input.Function == 'GWAs'",
                       fileInput(
                         inputId = "phenodata2",
                         label = "Upload Phenotype Data File",
                         accept = c(".txt", ".csv")),
                       uiOutput("dynamic_trait_select"),
                       fileInput(
                         inputId = "vcf",
                         label = "Upload VCF Data File",
                         accept = c(".vcf.gz", ".vcf")),
                       textInput("threshold", "Threshold", "5e-8"),
                       selectInput("showtop", 
                                   label = "Show Top SNPs",
                                   choices = list("Yes","No"),
                                   selected = "Yes"),
                       textInput("phenum", "Trait name", "Null"),
                       
      ),
      
    actionButton(inputId = "runAnalysis", label = "Run Analysis")
    ),
    
    mainPanel(
      h3("This is a shiny app for complex trait and cross environment genetic architecture analysis "),
      h1(textOutput("text")),
      
    )
  )
))