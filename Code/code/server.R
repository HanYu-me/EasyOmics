library(data.table)
options(shiny.maxRequestSize = 3000*1024^2)
shinyServer(
  function(input, output) {
    
    observeEvent(input$phenodata2,{
        data=data.frame(fread(input$phenodata2$datapath))
        names=names(data)
        output$dynamic_trait_select <- renderUI({
          selectInput("selected_trait", "Select a trait:", choices = names)
        })
    })
    observeEvent(input$runAnalysis,{
      outpath=getwd()
      out=paste0(outpath,"/Analysis_Result/")
      dir.create(out)
      
      if(input$Function=="Phenotype geration"){
        out=paste0(out,"Phenotype_geration/")
        dir.create(out)
        
        output$text <- renderText({ 
          paste("You have selected", input$phenodata$datapath)
        })
        shell(paste("code/Rscript 230724_phenotype_process.R",input$phenodata$datapath,out))
      }
      if(input$Function=="GWAs"){
        shell(paste("code/Rscript 230720_GCTA_singletrait_GWAs.R",input$vcf$datapath,input$phenodata2$datapath,out,1,"aa",input$threshold))
      }
    })
  }
)