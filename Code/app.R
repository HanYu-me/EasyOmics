library(shiny)
library(shinybusy)
library(png)
library(grid)
library(gridExtra)
library(data.table)
library(shinyFiles)
library(stringr)
options(shiny.maxRequestSize = 3000*1024^2)
load("code/session.RData")
ui=fluidPage(
  add_busy_spinner(spin = "cube-grid",position = "full-page",height = "100px",width="100px"),
  titlePanel("EasyOmics"),
  sidebarLayout(
    sidebarPanel(
      selectInput("Function", 
                  label = "Choose an analysis function",
                  choices = list("Data Matching","Phenotype Analysis", "GWAs","COJO","Locus zoom",
                                 "OmicWAs", "Omic QTL","Two Trait MR","Omic SMR"),
                  selected = "GWAs"),
      
      # Data Matching --------------------------------
      conditionalPanel(condition="input.Function == 'Data Matching'",
                       fileInput(
                         inputId = "phe_match",
                         label = "Upload Phenotype Data File",
                         accept = c(".txt", ".csv")),
                       fileInput(
                         inputId = "vcf_match",
                         label = "Upload VCF Data File",
                         accept = c(".vcf.gz", ".vcf")),                  
      ),
      # phenotype analysis --------------------------------
      conditionalPanel(condition="input.Function == 'Phenotype Analysis'",
                       fileInput(
                         inputId = "phenodata",
                         label = "Upload Phenotype Data File",
                         accept = c(".txt", ".csv")),
                       fileInput(
                         inputId = "vcf_pp",
                         label = "Upload VCF Data File",
                         accept = c(".vcf.gz", ".vcf")),
                       selectInput("analysis_type", 
                                   label = "Analysis type",
                                   choices = list("Single_Trait","Multi_Traits"),
                                   selected = "Single_Trait")
      ),
      # GWAs --------------------------------
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
                       textInput("color", "Color", "EAA4A4_96DCF8"),
                       selectInput("showtop", 
                                   label = "Show Top SNPs",
                                   choices = list("Yes","No"),
                                   selected = "Yes")
                       ),
      # locuszoom --------------------------------
      conditionalPanel(condition="input.Function == 'Locus zoom'",
                        fileInput(
                          inputId = "vcf_locus",
                          label = "Upload VCF Data File",
                          accept = c(".vcf.gz", ".vcf")),
                        fileInput(
                          inputId = "gff_locus",
                          label = "Upload GFF Data File",
                          accept = c(".gff3.gz", ".txt",".gff3")),
                        fileInput(
                          inputId = "gwa_locus",
                          label = "Upload GWA Data File",
                          accept = c(".mlma", ".txt")),
                        uiOutput("dynamic_trait_select2"),
                        fileInput(
                         inputId = "phenodata3",
                         label = "Upload Phenotype Data File",
                         accept = c(".txt", ".csv")),
                        textInput("topsnp", "SNP ID", "Eg. 1:1111"),
                        textInput("region", "Region (bp)", "10000"),
      ),
      # cojo --------------------------------          
      conditionalPanel(condition="input.Function == 'COJO'",
                        fileInput(
                          inputId = "mlm_cojo",
                          label = "Upload GWAs Data File",
                          accept = c(".txt", ".mlma")),
                        fileInput(
                          inputId = "vcf_cojo",
                          label = "Upload VCF Data File",
                          accept = c(".vcf.gz", ".vcf")),
                        textInput("threshold_cojo", "Threshold", "5e-8"),
                        textInput("color_cojo", "Color", "EAA4A4_96DCF8"),
                        selectInput("showtop_cojo", 
                                   label = "Show Top SNPs",
                                   choices = list("Yes","No"),
                                   selected = "Yes")
                        ),
      # OmicWAs --------------------------------
      conditionalPanel(condition="input.Function == 'OmicWAs'",
                       fileInput(
                         inputId = "phenodata_twas",
                         label = "Upload Phenotype Data File",
                         accept = c(".txt", ".csv")),
                       fileInput(
                         inputId = "exp_twas",
                         label = "Upload Omic Data File",
                         accept = c(".txt", ".csv")),
                       fileInput(
                         inputId = "gtf",
                         label = "Upload GFF Data File",
                         accept = c(".gff3.gz", ".gff3")),
                       textInput("threshold_twas", "Threshold", "5e-8"),
                       textInput("color_twa", "Color", "EAA4A4_96DCF8"),
                       selectInput("showtop_twas", 
                                   label = "Show Top Genes",
                                   choices = list("Yes","No"),
                                   selected = "Yes")
      ),      
      # Omic QTL --------------------------------
      conditionalPanel(condition="input.Function == 'Omic QTL'",
                       fileInput(
                         inputId = "omic_phe",
                         label = "Upload Omic Phenotype Data File",
                         accept = c(".txt", ".csv")),
                       selectInput("omic_type", 
                                   label = "Select Omic Type",
                                   choices = list("Transcriptom","Methylomics","Others"),
                                   selected = "Transcriptom"),
                       fileInput(
                         inputId = "gtf_qtl",
                         label = "Upload GFF Data File",
                         accept = c(".gff3.gz", ".gff3")),
                       fileInput(
                         inputId = "vcf_qtl",
                         label = "Upload VCF Data File",
                         accept = c(".vcf.gz", ".vcf")),
                       selectInput("norm_qtl", 
                                   label = "Normalization Data",
                                   choices = list("Yes","No"),
                                   selected = "Yes"),
                       textInput("threshold_qtl", "Threshold", "5e-8"),
      ),
      # MR --------------------------------                
      conditionalPanel(condition="input.Function == 'Two Trait MR'",
                       fileInput(
                         inputId = "exposure_single",
                         label = "Upload GWAs Result of Exposure",
                         accept = c(".mlma")),
                       fileInput(
                         inputId = "outcome_single",
                         label = "Upload GWAs Result of Outcome",
                         accept = c(".mlma")),
                       fileInput(
                         inputId = "vcf_single",
                         label = "Upload VCF Data File",
                         accept = c(".vcf",".vcf.gz")),
                       textInput("threshold_single", "Threshold", "5e-8"),
      ),
      # SMR --------------------------------
      conditionalPanel(condition="input.Function == 'Omic SMR'",
                       fileInput(
                         inputId = "exposure_smr",
                         label = "Upload Omic QTL Result of Exposure",
                         accept = c(".txt")),
                       fileInput(
                         inputId = "outcome_smr",
                         label = "Upload GWAs Result of Outcome",
                         accept = c(".mlma")),
                       fileInput(
                         inputId = "vcf_smr",
                         label = "Upload VCF Data File",
                         accept = c(".vcf",".vcf.gz")),
                       fileInput(
                         inputId = "gtf_smr",
                         label = "Upload GFF Data File",
                         accept = c(".gff3",".gff3.gz")),
                       textInput("threshold_smr", "Threshold", "5e-8"),
      ),
      
    actionButton(inputId = "runAnalysis", label = "Run Analysis")
    ),
    # main panel --------------------------------
    mainPanel(
      add_busy_spinner(spin = "fading-circle"),
      h3("A graphical interface for population-scale omics data association, integration and visualization"),
      uiOutput("dynamic_result_select"),
      uiOutput("dynamic_output"),
    )
  )
)


#dir.create("Analysis_Result/")
write.table("",file="info.txt",append=T,row.names = F,col.names = F,quote=F)
tryCatch({dir.create("Analysis_Result")}, warning = function(w){message("")})
tryCatch({dir.create("Analysis_Result/Data_Matching")}, warning = function(w){message("")})
tryCatch({dir.create("Analysis_Result/Phenotype_Analysis")}, warning = function(w){message("")})
tryCatch({dir.create("Analysis_Result/GWAs")}, warning = function(w){message("")})
tryCatch({dir.create("Analysis_Result/locuszoom")}, warning = function(w){message("")})
tryCatch({dir.create("Analysis_Result/COJO")}, warning = function(w){message("")})
tryCatch({dir.create("Analysis_Result/OmicWAs")}, warning = function(w){message("")})
tryCatch({dir.create("Analysis_Result/Omic_QTL")}, warning = function(w){message("")})
tryCatch({dir.create("Analysis_Result/MR")}, warning = function(w){message("")})
tryCatch({dir.create("Analysis_Result/SMR")}, warning = function(w){message("")})
message("-------------------Analysis Begin-------------------")
message("All files used to store the analysis results have been created in your binded file !")
message("")
message("----------------------------------------------------")
message("")
message("")


server=function(input, output) {
    outpath=getwd()
    #out=paste0(outpath,"/inter_result/")
    #dir.create(out)
    #system(paste("chmod 777",out))
    ##add an info file and set 777
    
    #system(paste("chmod 777","info.txt"))
    names=""
    observeEvent(input$phenodata2,{
        data=data.frame(fread(input$phenodata2$datapath))
        names<<-names(data)
        names<<-names[3:length(names)]
        write.table(names,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
        output$dynamic_trait_select <- renderUI({
          selectInput("selected_trait", "Select a trait:", choices=names)
        })
    })
    observeEvent(input$phenodata3,{
        data=data.frame(fread(input$phenodata3$datapath))
        names<<-names(data)
        names<<-names[3:length(names)]
        write.table(names,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
        output$dynamic_trait_select2 <- renderUI({
          selectInput("selected_trait2", "Select a trait:", choices=names)
        })
    })
    
    #set the layout of result --------------------------------
    observeEvent(input$Function,{
      output$dynamic_output=renderUI({
          tabsetPanel(
            tabPanel("Panel and Workflow", 
                     fluidRow(
                       column(6,div(style="height:0px;background-color: yellow;") ,plotOutput("Main")),
                     )
            )
            
          )
          
        })
        output$Main <- renderImage({
              list(src = "code/Main_fig.png", contentType = "image/png",height = 700)
            }, deleteFile = FALSE)
      # "Phenotype Analysis", "GWAs","COJO","Locus zoom","TWAs", "Omic QTL","Two Trait MR","Omic SMR"
     #if(input$Function=="Data Matching"){
     #  
     #}else if(input$Function=="Phenotype Analysis"){
     #  
     #}else if(input$Function=="GWAs"){
     # 
     #} else if(input$Function=="COJO"){
     #  
     #} else if(input$Function=="Locus zoom"){
     #  
     #}else if(input$Function=="OmicWAs"){
     # 
     #}else if(input$Function=="Omic QTL"){
     # 
     #}else if(input$Function=="Two Trait MR"){
     # 
     #}else if(input$Function=="Omic SMR"){
     # 
     #}
    })
    
    
    observeEvent(input$runAnalysis,{
      
      outpath=getwd()
      out=paste0(outpath,"/Analysis_Result/")
      #dir.create(out)
      #system(paste("chmod 777",out))
      #cat("1")
      # Data Matching --------------------------------
      if(input$Function=="Data Matching"){
        output$dynamic_output=renderUI({
          tabsetPanel(
            tabPanel("Results", 
                     fluidRow(
                       column(6,div(style="height:0px;background-color: yellow;") ,plotOutput("match_plot1")),
                     )
            )
          )
        })
        out=paste0(out,"Data_Matching/")
        #dir.create(out)
        ###filter vcf based on phe data
        message("\n-------------------Data Matching-------------------")
        message("VCF file format checking......")
        if(length(grep(".gz$",input$vcf_match$datapath))==1){
          system(paste0('zcat ', input$vcf_match$datapath,' | grep "^#[a-zA-Z]">inter_result/ids.txt'))
        }else{
          system(paste0('cat ', input$vcf_match$datapath,' | grep "^#[a-zA-Z]">inter_result/ids.txt'))
        }
        vcf_id=t.data.frame(data.frame(fread("inter_result/ids.txt",header = F)))
        vcf_id=data.frame(vcf_id[10:nrow(vcf_id),],stringsAsFactors=F)
        if(sum(grepl("_",vcf_id[,1]))==nrow(vcf_id)){
          message("The VCF file contains family information")
          vcf_id$family=str_split_fixed(vcf_id[,1],"_",2)[,1]
          vcf_id$iid=str_split_fixed(vcf_id[,1],"_",2)[,2]
          phe_id=data.frame(fread(input$phe_match$datapath))
          names(phe_id)[1:2]=c("family","id")
          ids=paste0(phe_id[,1],"_",phe_id[,2])
          all_id=intersect(ids,vcf_id[,1])
          phe_id=phe_id[ids%in%all_id,]
          message("Extracting the common individuals from VCF file and phenotype file")
          write.table(phe_id,file=paste0(out,"Matched_phe.txt"),row.names=F,quote=F)
          plink="../../home/software/plink"
          system(paste0("awk '{print $1,$2}' ",paste0(out,"Matched_phe.txt")," > inter_result/","keepid.txt"))
          message("Filtering VCF based on MAF and LD")
          system(paste0(plink," --vcf ",input$vcf_match$datapath," --keep inter_result/","keepid.txt"," --maf 0.03 --indep-pairwise 20 10 0.99 --allow-extra-chr --keep-allele-order --recode vcf --out ",paste0(out,"Matched_before_vcf")))
          system(paste("/usr/bin/python3 code/vcf_format.py",paste0(out,"Matched_before_vcf.vcf"),paste0(out,"Matched_vcf.vcf")))
        }else{
          message("The VCF file don't contain family information")
          vcf_id2=paste0(vcf_id[,1],"_",vcf_id[,1])
          phe_id=data.frame(fread(input$phe_match$datapath))
          names(phe_id)[1:2]=c("family","id")
          ids=paste0(phe_id[,1],"_",phe_id[,2])
          all_id=intersect(ids,vcf_id2)
          phe_id=phe_id[ids%in%all_id,]
          message("Extracting the common individuals from VCF file and phenotype file")
          write.table(phe_id,file=paste0(out,"Matched_phe.txt"),row.names=F,quote=F)
          plink="../../home/software/plink"
          system(paste("/usr/bin/python3 code/vcf_format.py",input$vcf_match$datapath,paste0("inter_result/","Matched_vcf.vcf")))
          system(paste0("awk '{print $1,$2}' ",paste0(out,"Matched_phe.txt")," > inter_result/","keepid.txt"))
          message("Filtering VCF based on MAF and LD")
          system(paste0(plink," --vcf ",paste0("inter_result/","Matched_vcf.vcf")," --keep inter_result/","keepid.txt"," --maf 0.03 --indep-pairwise 20 10 0.99 --allow-extra-chr --keep-allele-order --recode vcf --out ",paste0(out,"Matched_vcf")))
        }
        if(file.exists(paste0(out,"Matched_vcf.vcf.gz"))){
          system(paste("rm ",paste0(out,"Matched_vcf.vcf.gz")))
          message("Removing previous files")
        }
        system(paste("gzip",paste0(out,"Matched_vcf.vcf")))
        system(paste("rm ",paste0(out,"Matched_before_vcf.vcf")))
        write.table("aa",file="inter_result/info.txt",append=T,row.names = F,col.names = F,quote=F)
        system("rm inter_result/*")
        message("---------------------Finished----------------------\n")
      }
      
      # Phenotype Analysis --------------------------------
      if(input$Function=="Phenotype Analysis"){
        output$dynamic_output=renderUI({
          tabsetPanel(
            tabPanel("Results", 
                     fluidRow(
                       column(6,div(style="height:0px;background-color: yellow;") ,plotOutput("phe_plot1")),
                     )
            )
          )
        })
        message("\n-------------------Phenotype Analysis-------------------")
        out=paste0(out,"Phenotype_Analysis/")
        #dir.create(out)
        system(paste("Rscript code/230724_phenotype_process.R",input$phenodata$datapath,out,input$vcf_pp$datapath,input$analysis_type))
        if(input$analysis_type=="Multil_env"&file.exists((paste0(out,"Rawenv_phenotype_summary.png")))){
          output$phe_plot1 <- renderImage({
            list(src = paste0(out,"Rawenv_phenotype_summary.png"), contentType = "image/png",height = 800)
          }, deleteFile = FALSE)  
        }else if(input$analysis_type=="Multi_Traits"&file.exists((paste0(out,"Rawtra_phenotype_summary.png")))){
          output$phe_plot1 <- renderImage({
            list(src = paste0(out,"Rawtra_phenotype_summary.png"), contentType = "image/png",height = 800)
          }, deleteFile = FALSE)  
        }else if(input$analysis_type=="Single_Trait"&file.exists((paste0(out,"Single_Trait_result.png")))){
            output$phe_plot1 <- renderImage({
            list(src = paste0(out,"Single_Trait_result.png"), contentType = "image/png",height = 800)
          }, deleteFile = FALSE)  
        }else{
          output$phe_plot1 <- renderImage({
            list(src ="code/Feedback.png", contentType = "image/png",height = 500)
          }, deleteFile = FALSE)
        }
        write.table("aa",file="inter_result/info.txt",append=T,row.names = F,col.names = F,quote=F)
        system("rm inter_result/*")
        message("---------------------Finished----------------------\n")
      }
      
      # GWAs --------------------------------
      if(input$Function=="GWAs"){
         output$dynamic_output=renderUI({
          tabsetPanel(
            tabPanel("Results", 
                     fluidRow(
                       column(6,div(style="height:0px;background-color: yellow;") ,plotOutput("GWAs_plot1")),
                     )
            )
          )
        })
        out=paste0(out,"GWAs/")
        #dir.create(out)
        #system(paste("chmod 777",out))
        #write.table(input$color,file="info.txt",append=T,row.names = F,col.names = F,quote=F) 
        message("\n-------------------GWAS-------------------")
        system(paste("Rscript code/230720_GCTA_singletrait_GWAs.R",
                    input$vcf$datapath,
                    input$phenodata2$datapath,
                    out,
                    which(names==input$selected_trait),
                    input$selected_trait,
                    input$threshold,
                    ifelse(input$showtop=="Yes",T,F),as.character(input$color)))
        
        #check the existance of corrected result
        if(file.exists(paste0(out,input$selected_trait,"manhattan.png"))){
          if(file.exists(paste0(out,input$selected_trait,"manhattan_corrected.png"))){
            output$GWAs_plot1 <- renderImage({
              list(src = paste0(out,input$selected_trait,"manhattan_corrected.png"), contentType = "image/png",height = 800)
            }, deleteFile = FALSE)
          }else{
            output$GWAs_plot1 <- renderImage({
              list(src = paste0(out,input$selected_trait,"manhattan.png"), contentType = "image/png",height = 800)
            }, deleteFile = FALSE)
          }
        }else{
          output$GWAs_plot1 <- renderImage({
            list(src = "code/Feedback.png", contentType = "image/png",height = 800)
          }, deleteFile = FALSE)
        }
        #system("rm inter_result/*")
        message("---------------------Finished----------------------\n")
      }

      # Locuszoom --------------------------------
      if(input$Function=="Locus zoom"){
        message("\n-------------------Locus zoom-------------------")
        output$dynamic_output=renderUI({
          tabsetPanel(
            tabPanel("Results", 
                     fluidRow(
                       column(6,div(style="height:0px;background-color: yellow;") ,plotOutput("locus_plot1")),
                     )
            )
          )
        })
        out=paste0(out,"locuszoom/")
        #dir.create(out)
        #system(paste("chmod 777",out))
        system(paste("Rscript code/230920_locus.R",input$vcf_locus$datapath,input$gff_locus$datapath,input$gwa_locus$datapath,input$region,out,input$topsnp,input$phenodata3$datapath))
        ###draw phenotypic boxplot
        write.table("over",file="info.txt",append=T,row.names = F,col.names = F,quote=F) 
        write.table(paste0('cat ', input$vcf_locus$datapath,' | grep "^#[a-zA-Z]">inter_result/ids.txt'),file="info.txt",append=T,row.names = F,col.names = F,quote=F) 
        if(length(grep(".gz",input$vcf_locus$datapath))==1){
          system(paste0('zcat ', input$vcf_locus$datapath,' | grep "^#[a-zA-Z]">inter_result/ids.txt'))
          system(paste0('zcat ', input$vcf_locus$datapath,' | grep "',input$topsnp,'">inter_result/genotype.txt'))
        }else{
          system(paste0('cat ', input$vcf_locus$datapath,' | grep "^#[a-zA-Z]">inter_result/ids.txt'))
          system(paste0('cat ', input$vcf_locus$datapath,' | grep "',input$topsnp,'">inter_result/genotype.txt'))
        }
        tryCatch({
          locus_phe=data.frame(fread(input$phenodata3$datapath))
          locus_phe=locus_phe[!is.na(locus_phe[,3]),]
          locus_id=data.frame(fread("inter_result/ids.txt",header = F))
          locus_gen=data.frame(fread("inter_result/genotype.txt",header = F))
          names(locus_phe)[1:2]=c("id","family")
          locus_phe$id=as.character(locus_phe$id)
          locus_phe=locus_phe[,c("family","id",input$selected_trait2)]
          locus_gen=t.data.frame(rbind(locus_id,locus_gen))[-1:-10,]
          row.names(locus_phe)=locus_phe$id
          row.names(locus_gen)=str_split_fixed(locus_gen[,1],"_",2)[,2]
          all_id=intersect(row.names(locus_gen),locus_phe$id)
          locus_phe=locus_phe[all_id,]
          locus_gen=locus_gen[all_id,]
          locus_phe=locus_phe[!is.na(locus_phe[,3]),]
          locus_phe$gen=locus_gen[locus_phe$id,2]
          ###plot phenotype boxplot and write ids
          locus_phe_dis=data.frame(matrix(nrow=nrow(locus_phe),ncol=2*length(names(table(locus_phe$gen)))))
          for(i in 1:length(names(table(locus_phe$gen)))){ # nolint
            sub=locus_phe[locus_phe$gen==names(table(locus_phe$gen))[i],]
            locus_phe_dis[1:nrow(sub),(i*2-1)]=sub[,3]
            locus_phe_dis[1:nrow(sub),(i*2)]=sub[,2]
          }
          pdf(paste0(out,input$selected_trait2,sub(":","_",input$topsnp),"_boxplot.pdf"),width=8,height=6)
          boxplot(locus_phe_dis[,c(1:length(names(table(locus_phe$gen)))*2-1)],
                  names=paste("Gen",names(table(locus_phe$gen)),"num",apply(locus_phe_dis[,c(1:length(names(table(locus_phe$gen)))*2)],2,function(x)sum(!is.na(x)))),
                  frame.plot=F,col="skyblue",ylab="Phenotype")        
          dev.off()
          names(locus_phe_dis)=rep(names(table(locus_phe$gen)),each=2)
          for(i in 1:length(names(table(locus_phe$gen)))){
            order_data=order(locus_phe_dis[!is.na(locus_phe_dis[,(i*2-1)]),(i*2-1)],decreasing=T)
            locus_phe_dis[!is.na(locus_phe_dis[,(i*2-1)]),(i*2-1)]=locus_phe_dis[!is.na(locus_phe_dis[,(i*2-1)]),(i*2-1)][order_data]
            locus_phe_dis[!is.na(locus_phe_dis[,(i*2-1)]),(i*2)]=locus_phe_dis[!is.na(locus_phe_dis[,(i*2-1)]),(i*2)][order_data]
          }
          write.csv(locus_phe_dis,file=paste0(out,input$selected_trait2,sub(":","_",input$topsnp),"_boxplot_data.csv"),
                    row.names=F,quote=F)
        }, warning = function(w){message("")},
        error = function(e){message("")})
        

        if(file.exists((paste0(out,"result_",sub(":","_",input$topsnp),".png")))){
          output$locus_plot1 <- renderImage({
            list(src = paste0(out,"result_",sub(":","_",input$topsnp),".png"), contentType = "image/png",height = 800)
          }, deleteFile = FALSE)  
        }else{
          output$locus_plot1 <- renderImage({
            list(src ="code/Feedback.png", contentType = "image/png",height = 500)
          }, deleteFile = FALSE)
        }
      write.table("aa",file="inter_result/info.txt",append=T,row.names = F,col.names = F,quote=F)
      system("rm inter_result/*")
      message("---------------------Finished----------------------\n")
      }
      
      # cojo --------------------------------
      if(input$Function=="COJO"){
        output$dynamic_output=renderUI({
          tabsetPanel(
            tabPanel("Results", 
                     fluidRow(
                       column(6,div(style="height:0px;background-color: yellow;") ,plotOutput("COJO_plot1")),
                     )
            )
          )
        })
        message("\n-------------------COJO-------------------")
        out=paste0(out,"COJO/")
        system("rm Analysis_Result/COJO/*")
        #dir.create(out)
        #system(paste("chmod 777",out))
        system(paste("Rscript code/230720_cojo.R",input$mlm_cojo$datapath,input$vcf_cojo$datapath,input$threshold_cojo,out,input$color_cojo,ifelse(input$showtop_cojo=="Yes",T,F)))

        if(file.exists(paste0(out,"GWAs","_cojoed","manhattan.png"))){
          if(file.exists(paste0(out,"GWAs","_cojoed","manhattan_corrected.png"))){
            output$COJO_plot1 <- renderImage({
              list(src = paste0(out,"GWAs","_cojoed","manhattan_corrected.png"), contentType = "image/png",height = 800)
            }, deleteFile = FALSE)
          }else{
            output$COJO_plot1 <- renderImage({
              list(src = paste0(out,"GWAs","_cojoed","manhattan.png"), contentType = "image/png",height = 800)
            }, deleteFile = FALSE)
          }
        }else{
          output$COJO_plot1 <- renderImage({
            list(src = "code/Feedback.png", contentType = "image/png",height = 800)
          }, deleteFile = FALSE)
        }
      write.table("aa",file="inter_result/info.txt",append=T,row.names = F,col.names = F,quote=F)
      system("rm inter_result/*")
      message("---------------------Finished----------------------\n")
      }
      # TWAs --------------------------------
      if(input$Function=="OmicWAs"){
         output$dynamic_output=renderUI({
          tabsetPanel(
            tabPanel("Results", 
                     fluidRow(
                       column(6,div(style="height:0px;background-color: yellow;") ,plotOutput("TWAs_plot1")),
                     )
            )
          )
        })
        message("\n-------------------OmicWAs-------------------")
        out=paste0(out,"OmicWAs/")
        
        #dir.create(out)
        #system(paste("chmod 777",out))
        write.table(paste("Rscript code/230721_TWAs.R",input$exp_twas$datapath,input$phenodata_twas$datapath,input$gtf$datapath,out,input$threshold_twas,ifelse(input$showtop_twas=="Yes",T,F)),file="info.txt",append=T,row.names = F,col.names = F,quote=F)
        system(paste("Rscript code/230721_TWAs.R",input$exp_twas$datapath,input$phenodata_twas$datapath,input$gtf$datapath,out,input$threshold_twas,ifelse(input$showtop_twas=="Yes",T,F),input$color_twa))
        if(file.exists(paste0(out,"OmicWAs_manhattan.png"))){
          output$TWAs_plot1 <- renderImage({
            list(src = paste0(out,"OmicWAs_manhattan.png"), contentType = "image/png",height = 800)
          }, deleteFile = FALSE)  
        }else{
          output$TWAs_plot1 <- renderImage({
            list(src ="code/Feedback.png", contentType = "image/png",height = 500)
          }, deleteFile = FALSE)
        }
        write.table("aa",file="inter_result/info.txt",append=T,row.names = F,col.names = F,quote=F)
        system("rm inter_result/*")
        message("---------------------Finished----------------------\n")
      }
      
      # Omic QTL --------------------------------
      if(input$Function=="Omic QTL"){
         output$dynamic_output=renderUI({
          tabsetPanel(
            tabPanel("Results", 
                     fluidRow(
                       column(6,div(style="height:0px;background-color: yellow;") ,plotOutput("Omic_plot1")),
                     )
            )
          )
        })
        out=paste0(out,"Omic_QTL/")
        message("\n-------------------Omic QTL-------------------")
        #dir.create(out)
        #system("rm Analysis_Result/Omic_QTL/*")
        
        write.table(paste("Rscript code/230727_omicQTL.R",input$omic_phe$datapath,input$omic_type,input$gtf_qtl$datapath,input$vcf_qtl$datapath,ifelse(input$norm_qtl=="Yes","zscore",F),input$threshold_qtl,out),file="info.txt",append=T,row.names = F,col.names = F,quote=F)
        system(paste("Rscript code/230727_omicQTL.R",input$omic_phe$datapath,input$omic_type,input$gtf_qtl$datapath,input$vcf_qtl$datapath,ifelse(input$norm_qtl=="Yes","zscore",F),input$threshold_qtl,out))
        if(file.exists((paste0(out,"cis-trans_plot.png")))){
          output$Omic_plot1 <- renderImage({
            list(src = paste0(out,"cis-trans_plot.png"), contentType = "image/png",height = 800)
          }, deleteFile = FALSE)  
        }else{
          output$Omic_plot1 <- renderImage({
            list(src ="code/Feedback.png", contentType = "image/png",height = 500)
          }, deleteFile = FALSE)
        }
        write.table("aa",file="inter_result/info.txt",append=T,row.names = F,col.names = F,quote=F)
        system("rm inter_result/*")
        message("---------------------Finished----------------------\n")
      }
      
      # MR --------------------------------
      if(input$Function=="Two Trait MR"){
         output$dynamic_output=renderUI({
          tabsetPanel(
            tabPanel("Results", 
                     fluidRow(
                       column(6,div(style="height:0px;background-color: yellow;") ,plotOutput("MR_plot1")),
                     )
            )
          )
        })
        out=paste0(out,"MR/")
        message("\n-------------------Two Trait MR-------------------")

        #dir.create(out)
        #system(paste("chmod 777",out))
        
        system(paste("Rscript code/230726_single_molecular_SMR.R",input$exposure_single$datapath,input$outcome_single$datapath,input$vcf_single$datapath,input$threshold_single,out))
        if(file.exists(paste0(out,"gsmr_result.eff_plot.gz"))){
          output$MR_plot1 <- renderImage({
            list(src = paste0(out,"gsmr_result.png"), contentType = "image/png",height = 800)
          }, deleteFile = FALSE)  
        }else{
          output$MR_plot1 <- renderImage({
            list(src ="code/Feedback.png", contentType = "image/png",height = 500)
          }, deleteFile = FALSE)
        }
        write.table("aa",file="inter_result/info.txt",append=T,row.names = F,col.names = F,quote=F)
        system("rm inter_result/*")
        message("---------------------Finished----------------------\n")
      }
      
      # SMR --------------------------------
      if(input$Function=="Omic SMR"){
         output$dynamic_output=renderUI({
          tabsetPanel(
            tabPanel("Results", 
                     fluidRow(
                       column(6,div(style="height:0px;background-color: yellow;") ,plotOutput("SMR_plot1")),
                     )
            )
          )
        })
        message("\n-------------------Omic SMR-------------------")
        out=paste0(out,"SMR/")
        #dir.create(out)
        #system(paste("chmod 777",out))
        write.table(paste("Rscript code/230817_SMR_withmartrixeqtl.R",input$exposure_smr$datapath,input$outcome_smr$datapath,input$vcf_smr$datapath,input$gtf_smr$datapath,input$threshold_smr,out),file="info.txt",append=T,row.names = F,col.names = F,quote=F)    
        system(paste("Rscript code/230817_SMR_withmartrixeqtl.R",input$exposure_smr$datapath,input$outcome_smr$datapath,input$vcf_smr$datapath,input$gtf_smr$datapath,input$threshold_smr,out))
        if(dir.exists(paste0(out,"plot"))){
          output$dynamic_result_select <- renderUI({
            selectInput("selected_result", "Select a result:", choices = sub("smr.png","",list.files(out)[grep("\\.png$",list.files(out))]))
          })
        }else{
          output$SMR_plot1 <- renderImage({
            list(src ="code/Feedback.png", contentType = "image/png",height = 500)
          }, deleteFile = FALSE)
        }
        write.table("aa",file="inter_result/info.txt",append=T,row.names = F,col.names = F,quote=F)
        system("rm inter_result/*")
        message("---------------------Finished----------------------\n")
      }
    })
    
    observeEvent(input$selected_result,{
      outpath=getwd()
      out=paste0(outpath,"/Analysis_Result/")
      out=paste0(out,"SMR/")
      output$SMR_plot1 <- renderImage({
        list(src = paste0(out,input$selected_result,"smr.png"), contentType = "image/png",height = 1000)
      }, deleteFile = FALSE)
      print(paste0(out,input$selected_result,"smr.png"))
    })
  }


shinyApp(ui = ui, server = server)


