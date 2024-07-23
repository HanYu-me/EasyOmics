### Used the GSMR from GCTA to do SMR analysis and output the result

# Function ----------------------------------------------------------------
SMR_single_fun=function(exposure,outcome,out,threshold){
  suppressMessages(library(data.table))
  inter="inter_result/"
  bfile=paste0(inter,"bfile")
  
  #gcta=" ../softwares/gcta/exe/gcta-win-1.94.1.exe"
  gcta=" ../../home/software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"
  
  ##data format check and manipulation
  message("Converting mlm file format to GSMR")
  data1=data.table(fread(exposure))
  data2=data.table(fread(outcome))
  data1=data1[,c(2,4,5,6,7,8,9,10)]
  data2=data2[,c(2,4,5,6,7,8,9,10)]
  names(data1)=c("SNP","A1","A2","freq","b","se","p","N")
  names(data2)=c("SNP","A1","A2","freq","b","se","p","N")
  ## generate the gwa data
  write.table(data1,file=paste0(inter,"exp",".raw"),
              quote = F,col.names = T,row.names = F)
  write.table(data2,file=paste0(inter,"out",".raw"),
              quote = F,col.names = T,row.names = F)
  ## generate the summary data
  write.table(paste("exp",paste0(inter,"exp",".raw")),file=paste0(inter,"exposure.txt"),
              quote = F,col.names = F,row.names = F)
  write.table(paste("out",paste0(inter,"out",".raw")),file=paste0(inter,"outcome.txt"),
              quote = F,col.names = F,row.names = F)
  
  if(threshold=="Bonferroni"){
    thershold=0.05/nrow(data1)
  }else{
    threshold=as.numeric(threshold)
  }
  system("Performing two trait MR")
  cmd_gsmr=paste(gcta,"--bfile",bfile,"--gsmr-file",paste0(inter,"exposure.txt"),paste0(inter,"outcome.txt"),
                 "--gsmr-direction",0,"--effect-plot","--gwas-thresh",threshold,"--out",paste0(out,"gsmr_result"))
  system(cmd_gsmr)
  message("The MR result was saved in gsmr_result.gsmr")
}


##for raw mlm result

# Code ------------------------------------------------------------------
source("code/All_basic_function.R")
source("code/gsmr_plot.r")
args <- commandArgs(TRUE)
exposure=args[1]
outcome=args[2]
vcf=args[3]
threshold=args[4]
out=args[5]

data_convert_fun(vcf)
SMR_single_fun(exposure=exposure,outcome=outcome,out=out,threshold=threshold)


if(file.exists(paste0(out,"gsmr_result.eff_plot.gz"))){  
  message("There is a siginicant result")
  message("Plotting the effect point plot")
  gsmr_data=read_gsmr_data(paste0(out,"gsmr_result.eff_plot.gz"))
  #gsmr_summary(gsmr_data)
  png(paste0(out,"gsmr_result.png"),width=1200,height = 900)
  par(mar = c(5, 5, 5, 5))
  plot_gsmr_effect(gsmr_data, "exp", "out", colors()[75])          
  dev.off()
}else{
  message("No significant MR result")
  message("Maybe the over setted threshold")
  message("Please Check the .log file")
}



