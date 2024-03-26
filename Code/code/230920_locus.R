
# function ----------------------------------------------------------------
locus_fun=function(vcf,gff,gwa,region,out,snp){
  library(data.table)
  ldbs="../../home/software/LDBlockShow-1.40/bin/LDBlockShow"
  inter="inter_result/"
  data=data.frame(fread(gwa))
  if(ncol(data)==10){
    data=data[,c(1,2,3,9)]
  }
  region=as.numeric(region)
  write.table(data[,c(1,3,4)],file=paste0(inter,"ingwa.txt"),row.names = F,col.names = F,quote = F)
  reg=as.numeric(as.vector(data[which(data[,2]==snp),c(1,3)]))
  reg=paste0(reg[1],":",reg[2]-region,":",reg[2]+region)
  #process gff
  system(paste("/usr/bin/python3 code/gff_format.py",gff,paste0(inter,"newgff.gff3")))
  cmd_locus=paste(ldbs,"-InVCF",vcf,"-OutPut",paste0(out,"result_",snp),"-InGWAS", paste0(inter,"ingwa.txt"),
                  "-InGFF",paste0(inter,"newgff.gff3"),"-Region",reg,"-OutPng -SeleVar 2 -TopSite")
  write.table(cmd_locus,file="info.txt",append=T,row.names = F,col.names = F,quote=F) 
  system(cmd_locus)
}

# code --------------------------------------------------------------------
args <- commandArgs(TRUE)
vcf=args[1]
gff=args[2]
gwa=args[3]
region=args[4]
out=args[5]
snp=args[6]
phe=args[7]

#print(args)
locus_fun(vcf=vcf,gff=gff,gwa=gwa,region=region,out=out,snp=snp)