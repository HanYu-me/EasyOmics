
suppressMessages(library(stringr))
# Functions ---------------------------------------------------------------

###Using the GCTA to perform the gmlm GWAs
#phenotype data

GWA_fun=function(phe,out,phenum=1,name=""){
  #the name parameter should pre load the phe file and let user select the phe to analysis,
  #the name of the phe need to take to hear
  inter="inter_result/"
  #gcta address
  gcta="../../home/software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"
  cat("1")
  #make grm matrix
  
  cmd_grm=paste(gcta,"--bfile",paste0(inter,"bfile"),"--make-grm  --thread-num 10 --out",paste0(inter,"grm"))
  write.table(cmd_grm,file="info.txt",append=T,row.names = F,col.names = F,quote=F) 
  system(cmd_grm)
  #make sparse grm matrix
  # cmd_grm=paste(gcta,"--grm",paste0(inter,"grm"),"--make-bK-sparse",0.05,"--out",paste0(inter,"grm_sparse"))
  # system(cmd_grm)
  
  #select phe
  
  #mlm use this model
  message("Performing GWAS analysis with mixed linear model")
  cmd_mlm=paste(gcta,"--mlma","--bfile",paste0(inter,"bfile"), "--grm",paste0(inter,"grm"),"--pheno",phe,"--mpheno",phenum,"--thread-num 10","--out", paste0(out,name,"mlm"))
  system(cmd_mlm)
  
  #add a N col into mlma file
  data=data.frame(fread(paste0(out,name,"mlm.mlma")))
  phedata=data.frame(fread(phe))
  data$N=sum(!is.na(phedata[,3]))
  write.table(data,file=paste0(out,name,"mlm.mlma"),
              quote = F,col.names = T,row.names = F)
  message(paste0("The GWAS result was saved in ",name,"mlm.mlma"))
}


# code --------------------------------------------------------------------
#Rscript code/230720_GCTA_singletrait_GWAs.R "data/input/vcf.vcf" "data/input/single_trait/phe.txt" "result/single_trait/" 
#
write.table("begin",file="info.txt",append=T,row.names = F,col.names = F,quote=F) 
source("code/All_basic_function.R")
args <- commandArgs(TRUE)
vcf=args[1]
phe=args[2]
out=args[3]
phenum=args[4]
name=args[5]
threshold=args[6]
showtop=args[7]
color=args[8]
write.table(args,file="info.txt",append=T,row.names = F,col.names = F,quote=F) 

data_convert_fun(vcf)
GWA_fun(phe,out=out,phenum=phenum,name)

plot_fun(result=paste0(out,name,"mlm.mlma"),out=paste0(out,name),threshold=threshold,show_peakloci=showtop,color_manh=paste0("#",str_split(color,"_")[[1]]))



