library(stringr)

# Function ----------------------------------------------------------------
###This cojo was just for single trait
cojo_fun=function(mlm,cojop,out){
  inter="inter_result/"
  gcta=" ../../home/software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"
  data=data.frame(fread(mlm))
  data=data[,c(2,4,5,6,7,8,9,10)]
  names(data)=c("SNP","A1","A2","freq","b","se","p","N")
  write.table(data,file=paste0(inter,"mlm_for_cojo.ma"),
              quote = F,col.names = T,row.names = F)
  if(cojop=="Bonferroni"){
    cojop=0.05/nrow(data)
  }
  cmd_cojo=paste(gcta,"--bfile",paste0(inter,"bfile"),"--cojo-file",paste0(inter,"mlm_for_cojo.ma"),
                 "--cojo-slct","--cojo-p",cojop,"--out",paste0(out,"GWAs_cojoed"))
  write.table(cmd_cojo,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  system(cmd_cojo)
}


# code --------------------------------------------------------------------
source("code/All_basic_function.R")

args <- commandArgs(TRUE)
mlm=args[1]
vcf=args[2]
threshold=args[3]
out=args[4]
color=args[5]
showtop=args[6]
print(args)


data_convert_fun(vcf)
cojo_fun(mlm,threshold,out)
plot_fun(result=paste0(out,"GWAs","_cojoed.cma.cojo"),out=paste0(out,"GWAs","_cojoed"),threshold=threshold,show_peakloci=showtop,color_manh=paste0("#",str_split(color,"_")[[1]]))




