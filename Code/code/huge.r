###This function was used to convert the VCF file to bfile format by running in plink software
#This input should be provided by users
vcf="data/input/vcf.vcf"
a=function(){
  print(100)
}

a()

  
library(stringr)
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
  cmd_mlm=paste(gcta,"--mlma","--bfile",paste0(inter,"bfile"), "--grm",paste0(inter,"grm"),"--pheno",phe,"--mpheno",phenum,"--thread-num 10","--out", paste0(out,name,"mlm"))
  system(cmd_mlm)
  
  #add a N col into mlma file
  data=data.frame(fread(paste0(out,name,"mlm.mlma")))
  phedata=data.frame(fread(phe))
  data$N=sum(!is.na(phedata[,3]))
  write.table(data,file=paste0(out,name,"mlm.mlma"),
              quote = F,col.names = T,row.names = F)
  print("mlm finshed")
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
print("ploting")
plot_fun(result=paste0(out,name,"mlm.mlma"),out=paste0(out,name),threshold=threshold,show_peakloci=showtop,color_manh=paste0("#",str_split(color,"_")[[1]]))
print("plot finished")


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




load("C:\\Users\\MSI\\OneDrive\\桌面\\花期可塑性\\data230605.RData")


export.plink(data_af)

#phe
phe=data_af@phdata$THP
idname=data_af@phdata$id
phe=data.frame("family"=idname,"id"=idname,"phe"=phe)
phe$family=1:nrow(phe)
write.table(phe,file = "data/input/single_trait/phe.txt",col.names = F,row.names = F,quote = F)
#for multi-phe
phe=phdata(data_af)[,c(1,1,3:12)]
phe[,1]=1:nrow(phe)
names(phe)[1:2]=c("family","id")
write.table(phe,file = "phe.txt",col.names = T,row.names = F,quote = F)

#for single molecualr SMR
ph=ft1016[,c(1,1,3)]
ph[,1]=1:nrow(ph)
names(ph)=c("family","id","FT10")
write.table(ph,file = "C:\\Users\\MSI\\OneDrive\\桌面\\GWA_shiny_app\\data\\input\\SMR_type1\\phe_FT10.txt",col.names = T,row.names = F,quote = F)


ph=data.frame(family=1:nrow(ft1016),id=ft1016[,1],AT5G55080=p10[,"AT5G55080"])
write.table(ph,file = "C:\\Users\\MSI\\OneDrive\\桌面\\GWA_shiny_app\\data\\input\\SMR_type1\\phe_AT5G55080.txt",col.names = T,row.names = F,quote = F)

#for eQTL analysis
ph=p10
ph=ph[,-2]
ph=cbind(1:nrow(ph),ph)
names(ph)[1]="family"
write.table(ph[,c(1,2,sample(c(3:ncol(ph)),5000,replace=F))],file="exp10C.txt",quote=F,col.names = T,row.names = F)

#for docker test
ids=sample(data_af@phdata$id,200,replace=F)
export.plink(data_af[ids,])

phe=phdata(data_af)
phe=phe[phe$id%in%ids,c(1,3:12)]
phe$id=paste0(1:200,"_",phe$id)
write.table(phe,file = "phe.txt",col.names = F,row.names = F,quote = F)

d=data.frame(fread(file.choose()))
row.names(d)=d$V2
phe=phdata(data_af)
phe=phe[phe$id%in%d$V2,c(1,3:12)]
phe$id=paste0(d[phe$id,1],"_",d[phe$id,2])
write.table(phe,file = "phe2.txt",col.names = T,row.names = F,quote = F)



ph=data.frame(family=1:nrow(ft1016),id=ft1016[,1],AT2G33600=p10[,"AT2G33600"])
write.table(ph,file = "C:\\Users\\MSI\\OneDrive\\桌面\\GWA_shiny_app\\data\\input\\SMR_type1\\phe_AT2G33600.txt",col.names = T,row.names = F,quote = F)



samp_data=cbind.data.frame(data[,1:2],data[,sample(ncol(data)-2,5000,replace=F)+2])
dim(samp_data)
genes=c("AT1G35180","AT2G20440","AT2G45660","AT3G57230","AT4G33040","AT5G56380")
genes%in%names(samp_data)
samp_data=cbind.data.frame(samp_data,data[,genes[c(3,5,6)]])
write.table(samp_data,file="../input_V2/omicQTL/Matched_phe_sub5000.txt",row.names=F,quote=F)

snps
gene

exp_lw=data.frame(fread("../input_V2/omicQTL/exp_lw.txt"),stringAsFactor=F)
exp_lw[1:10,1:10]
exp=t(exp_lw)
exp[1:10,1:10]
exp=as.data.frame(exp,stringAsFactor=F)
names(exp)=as.character(exp[1,])
row.names(exp)=gsub("X","",row.names(exp))
exp=exp[c("V1",snps$columnNames),]
dim(exp)
exp=cbind(row.names(exp),exp)
exp=cbind(row.names(exp),exp)
names(exp)[1:2]=c("family","id")
write.table(exp,file="../input_V2/omicQTL/exp_lw2.txt",row.names=F,col.names=F,quote=F)

all_exp=data.frame(fread("../input_V2/omicQTL/Matched_phe.txt"))
sub=all_exp[,c("family","id","AT4G02260")]
write.table(sub,file="../input_V2/omicQTL/exp_AT4G02260.txt",row.names=F,col.names=T,quote=F)
library(stringr)
# Function ----------------------------------------------------------------
TWA_fun=function(exp,phe,out){
  inter="inter_result/"
  osca="../../home/software/osca-0.46.1-linux-x86_64/osca-0.46.1"
  #make bod expression file
  cmd_b=paste(osca,"--efile",exp,"--gene-expression","--make-bod","--no-fid","--out",paste0(inter,"bod"))
  write.table(cmd_b,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  system(cmd_b)
  
  #moa TWA
  cmd_moa=paste(osca,"--moa-exact","--befile",paste0(inter,"bod"),"--reml-maxit",200,"--thread-num 10 --pheno",phe,"--out",paste0(out,"moa"))
  write.table(cmd_moa,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  system(cmd_moa)
}
plot_twa_fun=function(result,out,gtf,color_manh=c("tomato","skyblue"),threshold=5e-8,corrected=T,show_peakloci=T){
  inter="inter_result/"
  write.table(threshold,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  if(threshold!="Bonferroni"){
    threshold=as.numeric(threshold)
    write.table(threshold,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  }
  write.table(threshold,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  library(data.table)
  ##get gene position
  ##Must set the input data type of chr as numeic(1,2,3...)
  system(paste("/usr/bin/python3 code/gff_format.py",gtf,paste0(inter,"twasgff.gff3")))
  gtf=data.frame(fread(paste0(inter,"twasgff.gff3"),fill=T))
  gene=data.frame(str_split_fixed(str_split_fixed(gtf$V9[which(gtf$V3=="gene")],";",2)[,1],":",2)[,2])
  gene$chr=gtf$V1[which(gtf$V3=="gene")]
  gene$BP=(gtf$V4[which(gtf$V3=="gene")]+gtf$V5[which(gtf$V3=="gene")])/2
  names(gene)[1]="Probe"
  row.names(gene)=gene$Probe
  chrlen=aggregate(BP~chr,data=gene,FUN=max)
  chrlen[,1]=as.character(chrlen[,1])
  chrlen[,2]=as.numeric(chrlen[,2])
  
  ##add a extend length to chr len, making the plot more beautiful
  chrlen[,2]=chrlen[,2]+sum(chrlen[,2])*0.02
  ##convert chr position to genomic position
  #get chr postion
  for(i in seq(nrow(chrlen),2)){
    chrlen$BP[i]=sum(chrlen$BP[1:(i-1)])
  }
  chrlen$BP[1]=0
  for(i in chrlen$chr){
    gene$BP[gene$chr==i]=gene$BP[gene$chr==i]+chrlen[chrlen[,1]==i,2]
  }
  
  
  data=data.frame(fread(result))
  print("this row may make some error in using new data, please check")
  data[,2]=sub("_10C","",data[,2])
  data[,1]=gene[data[,2],2]
  data[,3]=gene[data[,2],3]
  data=data[!is.na(data[,3]),]
  names(data)=c("CHR","GENE","POS","no1","no1","BETA","SE","P")
  
  #set color
  if(class(data[,1])=="integer"){
    cols=ifelse(data[,1]%%2,color_manh[1],color_manh[2])
  }else if(class(data[,1])!="integer"){
    #if the chr format are chr1 chr2 ....
    #I used the mean of chr position to order the chr, and give different colors to  adjacent chr.
    chr_mean=aggregate(POS~CHR,data=data,mean)
    cols=ifelse(order(chr_mean$POS)%%2,color_manh[1],color_manh[2])
    names(cols)=chr_mean$CHR[order(chr_mean$POS)]
    cols=cols[data$CHR]
  }

  png(file=paste0(out,"TWAs_manhattan.png"),width=15,height=10,units = "in",res=600)
  plot(x=data[,3],y=-log10(data[,8]),
       cex=1,pch=20,col=cols,frame.plot=F,xaxt="n",yaxt="n",
       xlab="Chromosome",ylab=expression("-Log"["10"]*"(P)"))
  chr_mean=aggregate(POS~CHR,data=data,mean)[,2]
  axis(1,chr_mean,c(unique(data$CHR)),las=1,cex.axis=1)
  axis(2,seq(0,max(-log10(data$P),na.rm=T),2),seq(0,max(-log10(data$P),na.rm=T),2),las=1,cex.axis=1)
  write.table(threshold,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  ##show the threshold
  if(threshold==5e-8){
    abline(h=-log10(threshold),col="red",lty="dashed")
  }else if((threshold!=5e-8)&is.numeric(threshold)){
    #threshold can be setted
    abline(h=-log10(threshold),col="red",lty="dashed")
    
  }else if(threshold=="Bonferroni"){
    #threshold can calculated
    threshold=0.05/nrow(data)
    abline(h=-log10(threshold),col="red",lty="dashed")
    write.table(threshold,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  }
  ##show top genes
  top_data=data[data$P<threshold,c(1:3,6:8)]
  write.table(threshold,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  if(nrow(top_data)>0){
    write.table(top_data,file=paste0(out,"top_genes.txt"),
                quote = F,col.names = T,row.names = F)
    if(show_peakloci==T){
      points(top_data$POS,-log10(top_data$P),cex=1.2,pch=17,col="red")
      text(top_data$POS-0.05*max(data$POS,na.rm=T),-log10(top_data$P),top_data$GENE,cex=0.8)
    }
  }else{
    write.table("No siginificant genes",file=paste0(out,"top_genes.txt"),
                quote = F,col.names = T,row.names = F)
  }
  dev.off()
}

# code --------------------------------------------------------------------
args <- commandArgs(TRUE)
exp=args[1]
phe=args[2]
gtf=args[3]
out=args[4]
thresh=args[5]
show_peak=args[6]
color=args[7]

write.table(args,file="info.txt",append=T,row.names = F,col.names = F,quote=F)

TWA_fun(exp=exp,phe=phe,out=out)
plot_twa_fun(result=paste0(out,"moa.moa"),out=out,gtf=gtf,threshold=thresh,show_peakloci=show_peak,color_manh=paste0("#",str_split(color,"_")[[1]]))


###the phe file should contain the family id phe info

# Function ----------------------------------------------------------------
##specific phenotypic plasticity of pairs phe difference
SP_fun=function(phe){
  data=phe
  envs=names(data)[3:ncol(data)]
  sp_data=data.frame(matrix(ncol=(length(envs)*(length(envs)-1))/2,nrow=nrow(data)))
  mark=1
  for(i in 1:(length(envs)-1)){
    for(j in (i+1):length(envs)){
      sp_data[,mark]=data[,envs[i]]-data[,envs[j]]
      names(sp_data)[mark]=paste0(envs[i],"_",envs[j])
      mark=mark+1
    }
  }
  return(sp_data)
}
OP_fun=function(phe){
  data=phe
  op_data=data.frame(matrix(ncol=6,nrow=nrow(data)))
  names(op_data)=c("blup","fwr","pc1","pc2","var","cv")
  ##BLUP 
  env_phe=tidyr::gather(data[,2:ncol(data)],key="envs",value="phe",names(data)[3:ncol(data)],na.rm = F)
  lm1=lmer(phe~envs+(1|id),data=env_phe)
  op_data[,"blup"]=lm1@u
  ##fwr
  env_phe=tidyr::gather(data[,2:ncol(data)],key="envs",value="phe",names(data)[3:ncol(data)],na.rm = F)
  lm2=lmFWh0(y=env_phe$phe,VAR=env_phe$id,ENV=env_phe$envs)
  op_data["fwr"]=as.numeric(lm2$b)
  ##pc
  mat=data[,3:ncol(data)]
  mat[is.na(mat)]=mean(as.numeric(unlist(mat)),na.rm=T)
  pca=prcomp(x=mat,center=T,scale=T)
  op_data[,"pc1"]=pca$x[,1]
  op_data[,"pc2"]=pca$x[,2]
  ##var
  env_phe=tidyr::gather(data[,2:ncol(data)],key="envs",value="phe",names(data)[3:ncol(data)],na.rm = F)
  envs=as.factor(env_phe$envs)
  phe=env_phe$phe
  bx=boxplot(phe~envs,plot=F)
  m_bx=bx$stats[3,]
  envs=factor(envs,levels=bx$names[order(m_bx,decreasing = F)])
  lm3=lm(phe~envs)
  beta=summary(lm3)$coefficients[,1]
  names(beta)=bx$names[order(m_bx,decreasing = F)]
  beta[1]=0
  phe_new=phe-beta[as.character(env_phe$envs)] #remove the environment effect
  env_phe$phe=phe_new
  op_data[,"var"]=aggregate(phe~id,data=env_phe,function(x){var(x,na.rm=T)})[,2]
  ##cv
  env_phe=tidyr::gather(data[,2:ncol(data)],key="envs",value="phe",names(data)[3:ncol(data)],na.rm = F)
  op_data[,"cv"]=aggregate(phe~id,data=env_phe,function(x){sd(as.numeric(x),na.rm=T)/mean(as.numeric(x),na.rm=T)})[,2]
  
  return(op_data)
}

##duplicated from https://github.com/lian0090/FW/blob/master/R/lmFWh0.R
lmFWh0=function(y,VAR,ENV){
  #if genotype or environment is completely missing for a GxE combination, the predicted value of  y is still NA.
  VAR = factor(VAR)
  ENV = factor(ENV)
  IDL = as.numeric(VAR)
  IDE = as.numeric(ENV)
  VARlevels = levels(VAR)
  ENVlevels = levels(ENV)
  n.var=length(VARlevels)
  n.env=length(ENVlevels)
  ##step 1 obtain environment effect (with sum contrast)
  # lm0=lm(y~-1+fENV+fVAR)
  # h=coef(lm0)[paste("fENV",ENVlevels,sep="")]
  #  h=h-mean(h,na.rm=T)
  fVARc=VAR; attr(fVARc, "contrasts") <- contr.sum(n.var) 
  fENVc=ENV; attr(fENVc,"contrasts")<-contr.sum(n.env)
  mf=model.frame(y~fENVc+fVARc)
  lm0=lm(mf)
  h=coef(lm0)[c(2:n.env)]
  h=c(h,-sum(h,na.rm=T))
  # h=tapply(y,INDEX=IDE,function(a)mean(a,na.rm=T))-mean(y,na.rm=T) 
  g=rep(0,n.var)
  b=rep(0,n.var)
  var_e=rep(0,n.var)
  df=rep(0,n.var)
  for(i in 1:n.var){
    whVar=which(IDL==i)
    lmi=lm(y[whVar]~h[IDE[whVar]])
    sum_lmi=summary(lmi)
    g[i]=coef(lmi)[1]
    b[i]=coef(lmi)[2]-1
    df[i]=sum_lmi$df[2]
    var_e[i]=(sum_lmi$sigma)^2
  }
  g=matrix(g)
  b=matrix(b)  
  h=matrix(h)
  rownames(g)=VARlevels
  rownames(b)=VARlevels
  rownames(h)=ENVlevels
  yhat=matrix(g[IDL,]+(1+b[IDL,])*h[IDE,]) #do not use VAR or ENV for index, because they might be factor
  ##this var_e will be exactly the same as if fitting all observations in a single linear model
  var_e_weighted=sum(var_e*df)/sum(df)
  
  LSvalue=list(y=y, whichNa = which(is.na(y)), VAR = VAR, ENV = ENV, VARlevels=VARlevels,ENVlevels=ENVlevels, mu = 0, g = g, b = b, h = h, yhat = yhat,var_e=var_e, var_e_weighted=var_e_weighted) 
  class(LSvalue)=c("FW","list")
  return(LSvalue)
}

genetic_data_fun=function(phe,vcf,name="",out){
  
  phe_data=data.frame(fread(phe))
  
  inter="inter_result/"
  #gcta="../softwares/gcta/exe/gcta-win-1.94.1.exe"
  gcta="../../home/software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"

  ##estimate h2 vp va
  cmd_grm=paste(gcta,"--bfile",paste0(inter,"bfile"),"--make-grm  --thread-num 10 --out",paste0(inter,"grm"))
  system(cmd_grm)
  summ=data.frame(va=1,vp=1,h2=1)[-1,]
  for(i in 1:(ncol(phe_data)-2)){
    print(i)
    cmd_h2=paste(gcta,"--reml","--grm",paste0(inter,"grm"),"--pheno",phe,"--mpheno",i,"--out",paste0(inter,"h2_",i))
    system(cmd_h2)
    h2=data.frame(fread(paste0(inter,"h2_",i,".hsq")))
    h2=h2[c(1,3,4),2]
    summ[i,]=h2
  }
  write.table(summ,file="info.txt",append=T,row.names = F,col.names = F,quote=F) 
  cols=colorRampPalette(brewer.pal(12, 'Set3'))(length(3:ncol(phe_data)))
  write.table(paste0(out,name,"phenotype_summary.png"),file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  kk=length(3:ncol(phe_data))
  if(kk>20){kk=20}
  png(file=paste0(out,name,"phenotype_summary.png"),width=kk,height=13,units = "in",res=600)
  par(mfrow=c(4,1))
  par(mar=c(3,5,3,1))
  boxplot(phe_data[,3:ncol(phe_data)],col=cols,frame.plot=F,ylab="Phenotype value")
  barplot(summ$vp,names=names(phe_data)[3:ncol(phe_data)],col=cols,ylab="Phenotypic variance")
  barplot(summ$va,names=names(phe_data)[3:ncol(phe_data)],col=cols,ylab="Additive variance")
  barplot(summ$h2,names=names(phe_data)[3:ncol(phe_data)],col=cols,ylab="Kinship heritability",ylim=c(0,1))
  write.table("what",file="info.txt",append=T,row.names = F,col.names = F,quote=F) 
  dev.off()
  row.names(summ)=names(phe_data)[3:ncol(phe_data)]
  if(name=="Plasticity_"){
    write.table(summ,file=paste0(out,"Phe_Pla_info.txt"),quote = F,col.names = T,row.names = T)  
  }else{
    write.table(summ,file=paste0(out,"Phe_info.txt"),quote = F,col.names = T,row.names = T)  
  }
  
  ##phe cor
  phe_cor=data.frame(matrix(nrow=length(3:ncol(phe_data)),ncol=length(3:ncol(phe_data))))
  for(i in 3:ncol(phe_data)){
    print(i)
    for(j in 3:ncol(phe_data)){
      phe_cor[(i-2),(j-2)]=cor.test(phe_data[,i],phe_data[,j])$estimate    
    }
  }
  phe_cor[upper.tri(phe_cor)]=NA
  row.names(phe_cor)=names(phe_data)[3:ncol(phe_data)]
  names(phe_cor)=names(phe_data)[3:ncol(phe_data)]
  png(file=paste0(out,name,"phenotype_cor.png"),width=11.5,height=10,units = "in",res=600)
  pheatmap(phe_cor,cluster_rows = F,cluster_cols = F,na_col = "white",
           border_color = NA)
  dev.off()
  if(name=="Plasticity_"){
    write.table(phe_cor,file=paste0(out,"Phe_Pla_cor.txt"),quote = F,col.names = T,row.names = T)  
  }else{
    write.table(phe_cor,file=paste0(out,"Phe_cor.txt"),quote = F,col.names = T,row.names = T)  
  }
  
  # ##genetic cor
  # for(i in 3:ncol(phe_data)){
  #   for(j in i:ncol(phe_data)){
  #     data=data.frame(1:nrow(phe_data),id=phe_data$id,phe1=phe_data[,i],phe2=phe_data[,j])
  #     write.table(data,file = paste0(inter,"p.txt"),row.names = F,col.names = F,quote = F)  
  #     cmd=paste(gcta,"--reml-bivar --grm",grm,"--pheno", paste0(inter,"p.txt")," --reml-maxit 200 --out",paste0(inter,"gc_",i,"_",j))
  #     system(cmd)    
  #   }
  # }
  # 
}
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}
c.z.hglm <- function(kin){
  relmat <- kin
  relmat[upper.tri(relmat)] <- t(relmat)[upper.tri(relmat)]#上三角变成下三角
  svd <- svd(relmat)
  Z <- svd$u %*% diag(sqrt(svd$d)) #左边的矩阵乘奇异值
  return(Z)
}

# Code --------------------------------------------------------------------
library(lme4)
library(tidyr)
library(data.table)
library(pheatmap)
library(png)
library(RColorBrewer)
#library(devtools)
#install_github("lian0090/FW",force=T)
#library(FW)
args <- commandArgs(TRUE)
print(args)
phe=args[1]
out=args[2]
vcf=args[3]
type=args[4]
source("code/All_basic_function.R")
data_convert_fun(vcf) 
if(type=="Mutli_Environments"){
  phe_data=data.frame(fread(phe))
  genetic_data_fun(phe=phe,vcf=vcf,name="Rawenv_",out=out)
  sp=SP_fun(phe_data)
  op=OP_fun(phe_data)
  phe_data=cbind.data.frame(phe_data,sp)
  phe_data=cbind.data.frame(phe_data,op)
  write.table(phe_data,file=paste0(out,"precessed_phe.txt"),quote = F,col.names = T,row.names = F)
  genetic_data_fun(phe=paste0(out,"precessed_phe.txt"),vcf=vcf,name="Plasticity_",out=out)
}else if(type=="Multli_Traits"){
  print(12)
  genetic_data_fun(phe=phe,vcf=vcf,name="Rawtra_",out=out)
}else if( type=="Single_Trait"){
  phe_data=data.frame(fread(phe))
  inter="inter_result/"
  gcta="../../home/software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"
  ##estimate h2
  cmd_grm=paste(gcta,"--bfile",paste0(inter,"bfile"),"--make-grm-gz  --thread-num 10 --out",paste0(inter,"grm"))
  system(cmd_grm)
  cmd_h2=paste(gcta,"--reml","--grm",paste0(inter,"grm"),"--pheno",phe,"--mpheno",1,"--out",paste0(out,"h2_single_",1))
  system(cmd_h2)
  ##load grm
  if(file.exists(paste0(inter,"grm",".grm"))){
          system(paste("rm ",paste0(inter,"grm",".grm")))
  }
  system(paste("gunzip ",paste0(inter,"grm",".grm.gz")))
  grmfile=data.frame(fread(paste0(inter,"grm",".grm")),stringsAsFactors = F)
  grm_data=matrix(ncol=max(grmfile[,1]),nrow=max(grmfile[,1]))
  for(i in 1:nrow(grmfile)){
    a=as.numeric(grmfile[i,1:4])
    grm_data[a[1],a[2]]=a[4] 
  }
  grm_data[upper.tri(grm_data)] <- t(grm_data)[upper.tri(grm_data)]
  Z=c.z.hglm(grm_data)
  pca <- prcomp(Z)
  pc1=pca$x[,1]
  pc2=pca$x[,2]
  explainability <- round(summary(pca)$importance[2, 1:2],2)
  
  if(sum(phe_data[,3]%in%c(1,0,-9))!=nrow(phe_data)){
    png(file=paste0(out,"Single_Trait_result.png"),width=15,height=7,units = "in",res=600)
    par(mfrow=c(1,2))
    plot(density(phe_data[,3],na.rm=T),frame.plot=F,main="Phenotype distribution")
    text(x=phe_data[,3],y=rep(0,nrow(phe_data)),phe_data[,2], srt=90,cex=0.8)
  }else{
    png(file=paste0(out,"Single_Trait_result.png"),width=8,height=7,units = "in",res=600)
  }
  plot(pc1,pc2,col="white",frame=F,
    xlab=paste("PC1 ",explainability[1]),ylab=paste("PC2 ",explainability[2]),
    main="Genetic structure")
  text(x=pc1,y=pc2,phe_data[,2],cex=0.8)
  dev.off()
}


#phe="data/input/phenotype/phe.txt"### Used the GSMR from GCTA to do SMR analysis and output the result

# Function ----------------------------------------------------------------
SMR_single_fun=function(exposure,outcome,out,threshold){
  library(data.table)
  inter="inter_result/"
  bfile=paste0(inter,"bfile")
  
  #gcta=" ../softwares/gcta/exe/gcta-win-1.94.1.exe"
  gcta=" ../../home/software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"
  
  ##data format check and manipulation
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
  
  cmd_gsmr=paste(gcta,"--bfile",bfile,"--gsmr-file",paste0(inter,"exposure.txt"),paste0(inter,"outcome.txt"),
                 "--gsmr-direction",0,"--effect-plot","--gwas-thresh",threshold,"--out",paste0(out,"gsmr_result"))
  system(cmd_gsmr)
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
  gsmr_data=read_gsmr_data(paste0(out,"gsmr_result.eff_plot.gz"))
  #gsmr_summary(gsmr_data)
  png(paste0(out,"gsmr_result.png"),width=1200,height = 900)
  par(mar = c(5, 5, 5, 5))
  plot_gsmr_effect(gsmr_data, "exp", "out", colors()[75])          
  dev.off()
}



# Function ----------------------------------------------------------------

omicQTL_fun=function(omic,type,gtf,vcf,norm,out,threshold){
  source("code/All_basic_function.R")
  osca="../../home/software/osca-0.46.1-linux-x86_64/osca-0.46.1"
  library(data.table)
  library(MatrixEQTL)
  inter="inter_result/"
  data=data.frame(fread(omic))
  #data[1:10,1:10]
  ##filter data
  index=as.logical(apply(data,2,function(x){
    return(
      sum(is.na(x))<=(length(x)/2)
      )
    }))
  index2=as.logical(apply(data,2,function(x){
    return(
      sum(na.omit(x==0))<=(length(x)/2)
    )
  }))
  index=index&index2
  names=names(data)[index]
  names=names[names%in%names(data)]
  data=data[,names]
  ##normalization
  if(norm=="zscore"){
    data=cbind.data.frame(data[,1:2],
                          apply(data[,3:ncol(data)],2,zscore))
  }
  
  ##matrix eQTL analysis
  library(MatrixEQTL)
  ##make exp file
  data=data.frame(t(data))
  data=cbind.data.frame(row.names(data),data)
  data=data[-1,]
  names(data)=as.character(as.numeric(data[1,]))
  names(data)[1]="gene"
  data=data[-1,]
  ###separate the data into five file to parallel analysis
  mark=split(1:nrow(data),cut(seq_along(1:nrow(data)),breaks=2,labels = F))
  for(i in 1:2){
    sub_data=data[mark[[i]],]
    write.table(sub_data,file=paste0(inter,"omic_",i,".txt"),
                row.names = F,col.names = T,quote=F,sep="\t",eol = "\n")
  }
  
  ##make snp file
  library(dplyr)
  #BiocManager::install("VariantAnnotation",type="binary")
  #BiocManager::install('snpStats')
  library(VariantAnnotation)
  library(snpStats)
  #convert.vcf(vcf.file=vcf,genotype_file_name = paste0(inter,"snp_for_eqtl.txt"))
  write.table("maked",file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  ##make cvrt(PC1-40) file 
  plink="softwares\\plink\\plink.exe"
  cmd_cvrt=paste("../../home/software/plink","--vcf",vcf,"--pca 20","--out",paste0(inter,"pca_40"))
  write.table(cmd_cvrt,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  system(cmd_cvrt)
  pca=read.table(file=paste0(inter,"pca_40",".eigenvec"))
  row.names(pca)=pca[,2]
  pca=t(pca[,-1:-2])
  row.names(pca)=paste(rep("pc_",20),1:20)
  write.table(pca,file=paste0(inter,"pca_40_eqtl.txt"),quote = F,sep="\t")
  
  #snp_text=data.frame(fread(paste0(inter,'snpmatrix2.txt')))
  #snp_text[1:10,1:10]
  #write.table(snp_text,row.names = F, col.names = F, quote = F, sep = "\t",file =paste0(inter,'snpmatrix3.txt') )
  ##load snp and cvrt
  snps = SlicedData$new();
  snps$fileDelimiter = " ";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 10000;      # read file in slices of 2,000 rows
  snps$LoadFile(paste0(inter,'snpmatrix2.txt'))

  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  cvrt$LoadFile(paste0(inter,"pca_40_eqtl.txt"))
  
  vcf_data=data.frame(fread(vcf))
  vcf_save=vcf_data[,1:10]
  row.names(vcf_save)=vcf_save[,3]
  if(threshold=="Bonferroni"){
    threshold=0.05/nrow(snps)
  }else{
    threshold=as.numeric(threshold)
  }
  ##do eQTL scan
  write.table("dataload",file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  library(parallel)
  clnum<-2
  cl <- makeCluster(getOption("cl.cores", clnum));
  clusterExport(cl,deparse(substitute(eQTL_para)))
  clusterExport(cl,c("inter","snps","cvrt"),envir=environment())
  clusterExport(cl,c("threshold"),envir=environment())
  clusterEvalQ(cl,library(MatrixEQTL))
  #5000 genes runtime 5min
  result_me=parLapply(cl,1:2,eQTL_para)
  stopCluster(cl)
  write.table("ok",file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  ##QQplot

  all_qtl=data.frame()
  for(i in 1:2){
    qtl=data.table(fread(paste0(inter,"omic_qtl_",i,".txt")))
    all_qtl=rbind.data.frame(all_qtl,qtl)
  }
  write.table(all_qtl,file=paste0(out,"qtls.txt"),quote = F,row.names = F,col.names = T)
  all_qtl=as.data.frame(all_qtl)
  ## remove linked QTL
  genes=unique(all_qtl$gene)
  vcf_data=vcf_data[vcf_data[,3]%in%all_qtl$SNP,]
  print("Filted the QTL based on distance")

  #all_qtl_filte=data.frame()
  #for(i in seq(1,length(genes))){
    #sub_qtl=all_qtl[all_qtl$gene==genes[i],]
    #sub_qtl=as.data.frame(sub_qtl)
    #if(nrow(sub_qtl)==1){
    #  all_qtl_filte=rbind.data.frame(all_qtl_filte,sub_qtl)
    #}else{
    #  sub_vcf=vcf_data[vcf_data[,3]%in%sub_qtl$SNP,]
    #  row.names(sub_vcf)=sub_vcf[,3]
    #  row.names(sub_qtl)=as.character(sub_qtl$SNP)
    #  sub_vcf=sub_vcf[order(sub_vcf[,2]),]
    #  pos=sub_vcf[sub_qtl[,1],2]
    #  sub_qtl=sub_qtl[order(pos),]
    #  sub_qtl=try(qtl_finder(sub_qtl,sub_vcf),silent=T)
    #  all_qtl_filte=rbind.data.frame(all_qtl_filte,sub_qtl)
    #}
  #}
  # clnum<-5
  # cl <- makeCluster(getOption("cl.cores", clnum));
  # clusterExport(cl,deparse(substitute(remove_linksnp)))
  # clusterExport(cl,c("all_qtl","vcf_data"),envir=environment())
  # clusterExport(cl,c("genes"),envir=environment())
  # #500 genes runtime 5mins
  # result=parLapply(cl,1:5,remove_linksnp)
  # stopCluster(cl)
  #write.table(all_qtl_filte,file=paste0(out,"qtls_filted.txt"),quote = F,row.names = F,col.names = T)
  ##cis-trans analysis
  if(type=="Transcriptom"){
    library(stringr)
    system(paste("/usr/bin/python3 code/gff_format.py",gtf,paste0(inter,"qtlgff.gff3")))
    gtf=data.frame(fread(paste0(inter,"qtlgff.gff3"),fill=T))
    gene=data.frame(str_split_fixed(str_split_fixed(gtf$V9[which(gtf$V3=="gene")],";",2)[,1],":",2)[,2])
    gene$chr=gtf$V1[which(gtf$V3=="gene")]
    gene$BP=(gtf$V4[which(gtf$V3=="gene")]+gtf$V5[which(gtf$V3=="gene")])/2
    names(gene)[1]="Probe"
    row.names(gene)=gene$Probe
    
    chrlen=aggregate(BP~chr,data=gene,FUN=max)
    chrlen=chrlen[grep("[0-9]+",chrlen[,1]),]
    chrlen[,1]=as.character(chrlen[,1])
    chrlen[,2]=as.numeric(chrlen[,2])
    if(sum(diff(chrlen[,2])>0)!=(nrow(chrlen)-1)){
      for(i in chrlen[,1]){
        if(i!=1){
          gene$BP[gene$chr==i]=gene$BP[gene$chr==i]+sum(chrlen[1:(as.numeric(i)-1),2])
        }      
      }
    } 
  }else{
    gene=read.csv(gtf,sep="\t",header=F)
    chrlen=aggregate(BP~chr,data=gene,FUN=max)
    chrlen=aggregate(BP~chr,data=gene,FUN=max)
    chrlen[,1]=as.character(chrlen[,1])
    chrlen[,2]=as.numeric(chrlen[,2])
    
    ##add a extend length to chr len, making the plot more beautiful
    chrlen[,2]=chrlen[,2]+sum(chrlen[,2])*0.02
    ##convert chr position to genomic position
    #get chr postion
    for(i in seq(nrow(chrlen),2)){
      chrlen$BP[i]=sum(chrlen$BP[1:(i-1)])
    }
    chrlen$BP[1]=0
    for(i in chrlen$chr){
      gene$BP[gene$chr==i]=gene$BP[gene$chr==i]+chrlen[chrlen[,1]==i,2]
    }  
  }
  chr_mean=aggregate(BP~chr,data=gene,mean)
  chr_mean=chr_mean[chr_mean[,1]%in%unique(vcf_data[,1]),2]
  chr_mean=chr_mean+max(chrlen$BP)
  
  print("Plot result")
  all_qtl=as.data.frame(all_qtl)
  gene_pos=gene[as.character(all_qtl[,2]),3]
  row.names(vcf_save)=vcf_save[,3]
  names(vcf_save)
  snp_name=gsub("_[A-Z]","",all_qtl[,1])
  snp_pos=vcf_save[snp_name,c(1,2)]
  
  snp_location=vcf_save[,1:3]
  names(snp_location)=c("Chr","POS","SNP")
  chr_info=aggregate(POS~Chr,data=snp_location,max)  
  if(sum(diff(chr_info[,2])>0)!=(nrow(chr_info)-1)){
    for(i in chr_info[,1]){
      if(i!=1){
        snp_pos[snp_pos[,1]==i,2]=snp_pos[snp_pos[,1]==i,2]+sum(chr_info[1:(i-1),2])
      }      
    }
  } 
  snp_pos=snp_pos[,2]
  ##plot result
  png(paste0(out,"cis-trans_plot.png"),width=10,height=10,unit="in",res=600)
  plot(snp_pos,gene_pos,
       frame.plot=F,pch=20,xaxt="n",yaxt="n",ylab="Gene position",xlab="QTL position")
  axis(1,chr_mean,unique(vcf_data[,1]),las=1,cex.axis=1)
  axis(2,chr_mean,unique(vcf_data[,1]),las=1,cex.axis=1)
  dev.off()
  
  histdata=hist(snp_pos,breaks=seq(0,max(snp_pos,na.rm=T)+5e4,5e4))
  png(paste0(out,"QTL_times.png"),width=10,height=10,unit="in",res=600)
  par(mar = c(5, 5, 5, 5))
  x=histdata$breaks;y=c(0,histdata$counts)
  plot(x,log10(y)+1,cex=log10(y)+1,
       frame.plot=F,pch=20,xaxt="n",ylab=expression("-Log"["10"]*"(QTL times)"),xlab="QTL position")
  axis(1,chr_mean,unique(vcf_data[,1]),las=1,cex.axis=1)
  dev.off()
  
  #all_qtl_filte=as.data.frame(all_qtl_filte)
  #gene_pos=gene[as.character(all_qtl_filte[,2]),3]
  #snp_pos=vcf_save[as.character(all_qtl_filte[,1]),2]
  #pdf(paste0(out,"cis-trans_plot_filted.pdf"),width=10,height=10)
  #plot(gene_pos,snp_pos,
  #     frame.plot=F,pch=20,xaxt="n",yaxt="n",ylab="QTL position",xlab="Gene position")
  #axis(1,chr_mean,unique(vcf_data[,1]),las=1,cex.axis=1)
  #axis(2,chr_mean,unique(vcf_data[,1]),las=1,cex.axis=1)
  #dev.off()
  
  #histdata=hist(snp_pos,breaks=seq(0,max(snp_pos,na.rm=T)+5e4,5e4))
  #pdf(paste0(out,"QTL_times_filted.pdf"),width=10,height=10)
  #par(mar = c(5, 5, 5, 5))
  #x=histdata$breaks;y=c(0,histdata$counts)
  #plot(x,log10(y)+1,cex=log10(y)+1,
  #     frame.plot=F,pch=20,xaxt="n",ylab=expression("-Log"["10"]*"(QTL times)"),xlab="QTL position")
  #axis(1,chr_mean,unique(vcf_data[,1]),las=1,cex.axis=1)
  #dev.off()

  ##plot merged result
  #png(file=paste0(out,"QTL_plot_merge.png"),width=1200,height=1000)
  #par(mfcol=c(2,2))
  #par(mar = c(2, 2, 2, 2))
  #plot(gene_pos,snp_pos,
  #     frame.plot=F,pch=20,xaxt="n",yaxt="n",ylab="QTL position",xlab="Gene position",main="Raw QTLs")
  #axis(1,chr_mean,unique(vcf_data[,1]),las=1,cex.axis=1)
  #axis(2,chr_mean,unique(vcf_data[,1]),las=1,cex.axis=1)
  
  
  #histdata=hist(snp_pos,breaks=seq(0,max(snp_pos,na.rm=T)+5e4,5e4),plot=F)
  #par(mar = c(5, 5, 5, 5))
  #x=histdata$breaks;y=c(0,histdata$counts)
  #plot(x,log10(y)+1,cex=log10(y)+1,
  #     frame.plot=F,pch=20,xaxt="n",ylab=expression("-Log"["10"]*"(QTL times)"),xlab="QTL position")
  #axis(1,chr_mean,unique(vcf_data[,1]),las=1,cex.axis=1)
  
  
  #all_qtl_filte=as.data.frame(all_qtl_filte)
  #gene_pos=gene[as.character(all_qtl_filte[,2]),3]
  #snp_pos=vcf_save[as.character(all_qtl_filte[,1]),2]
  #plot(gene_pos,snp_pos,
   #    frame.plot=F,pch=20,xaxt="n",yaxt="n",ylab="QTL position",xlab="Gene position",main="Filtered QTLs")
  #axis(1,chr_mean,unique(vcf_data[,1]),las=1,cex.axis=1)
  #axis(2,chr_mean,unique(vcf_data[,1]),las=1,cex.axis=1)
  
  
  #histdata=hist(snp_pos,breaks=seq(0,max(snp_pos,na.rm=T)+5e4,5e4),plot=F)
  #par(mar = c(5, 5, 5, 5))
  #x=histdata$breaks;y=c(0,histdata$counts)
  #plot(x,log10(y)+1,cex=log10(y)+1,
  #     frame.plot=F,pch=20,xaxt="n",ylab=expression("-Log"["10"]*"(QTL times)"),xlab="QTL position")
  #axis(1,chr_mean,unique(vcf_data[,1]),las=1,cex.axis=1)
  #dev.off()
}


# Code --------------------------------------------------------------------
source("code/All_basic_function.R")
library(parallel)
inter="inter_result/"
args <- commandArgs(TRUE)
omic=args[1]
type=args[2]
#if type=Expression gtf file is the simple gtf file
#if type!=Expression gtf file must be a tab delimtated txt file with 3 colums 
#as probe name(same as the omic head name), chr, pos in chr or in genome. Can also be a empty file
gtf=args[3]
vcf=args[4]
inter="inter_result/"
if(length(grep(".gz",vcf))==1){
  system(paste("gunzip -c",vcf, paste0("> ",inter,"vcf_qtl.vcf")))
  vcf=paste0(inter,"vcf_qtl.vcf")
}
plink="../../home/software/plink"
system(paste0(plink," --vcf ",vcf," --recodeA --out ",paste0(inter,"vcf_qtl2")))
system(paste0('cat ',paste0(inter,'vcf_qtl2.raw'),' | cut -d" " -f2,7- |',"sed 's/_[A-Z]//g' > ",paste0(inter,'snpmatrix.txt')))
system("bash code/convert.sh")

#add normalization function,now only support Zsocre
norm=args[5]
threshold=args[6]
out=args[7]
write.table(args,file="info.txt",append=T,row.names = F,col.names = F,quote=F)

omicQTL_fun(omic=omic,type=type,gtf=gtf,vcf=vcf,norm=norm,out=out,threshold=threshold)
#data_convert_fun(vcf)
###the input data of SMR was generated by matrixeqtl package
###smr need chromomose is a numeric data
#Function -------------------------------------------------------------------------
smr_fun=function(exposure,outcome,out,threshold,gtf){
  inter="inter_result/"
  smr="../../home/software/smr-1.3.1-linux-x86_64/smr-1.3.1"
  #exposure="../input_V2/SMR/qtls.txt"
  #outcome="../input_V2/SMR/FT16mlm.mlma"


  ##load mlm result of output and delect some colums
  library(data.table)
  library(stringr)
  print("make file")
  data=data.frame(fread(outcome))
  gwa=data
  data=data[,c(2,4,5,6,7,8,9,10)]
  names(data)=c("SNP","A1","A2","freq","b","se","p","n")
  write.table(data,paste0(inter,"smr_output.ma"),row.names = F,col.names = T,quote = F)

  ##make esi file
  system(paste("/usr/bin/python3 code/gff_format.py",gtf,paste0(inter,"smrgff.gff3")))
  gtf=data.frame(fread(paste0(inter,"smrgff.gff3"),fill=T))
  gene=data.frame(str_split_fixed(str_split_fixed(gtf$V9[which(gtf$V3=="gene")],";",2)[,1],":",2)[,2])
  gene$chr=gtf$V1[which(gtf$V3=="gene")]
  gene$BP=(gtf$V4[which(gtf$V3=="gene")]+gtf$V5[which(gtf$V3=="gene")])/2
  gene$start=gtf$V4[which(gtf$V3=="gene")]
  gene$end=gtf$V5[which(gtf$V3=="gene")]
  names(gene)[1]="Probe"
  row.names(gene)=gene$Probe
  gene=gene[grep("[0-9]+",gene[,2]),]
  
  qtls=data.frame(fread(exposure))
  qtls=qtls[qtls[,1]%in%gwa[,2],]
  qtls=qtls[qtls[,2]%in%gene[,1],]
  exposure=paste0(inter,"qtls_V2.txt")
  write.table(qtls,file=exposure,row.names = F,col.names =T,quote=F)  
  
  row.names(gwa)=gwa$SNP
  esi=data.frame(gwa[qtls$SNP,"Chr"],qtls$SNP,0,gwa[qtls$SNP,"bp"],gwa[qtls$SNP,"A1"],gwa[qtls$SNP,"A2"],gwa[qtls$SNP,"Freq"])


  ##make epi file
  print("make epi file")

  
  epi=data.frame(gene[qtls$gene,"chr"],qtls$gene,0,gene[qtls$gene,"BP"],qtls$gene,"+")
  nanum=is.na(esi[,4])|is.na(epi[,4])
  esi=esi[!nanum,]
  epi=epi[!nanum,]
  write.table(esi,paste0(inter,"new.esi"),row.names = F,col.names = F,quote = F)
  write.table(epi,paste0(inter,"new.epi"),row.names = F,col.names = F,quote = F)
  ##make gene list file
  glist=data.frame(gene$chr,gene$start,gene$end,gene$Probe)
  write.table(glist,paste0(inter,"glist.txt"),row.names = F,col.names = F,quote = F)
  ##make a BESD file from matrixeqtl file
  cmd_besd=paste(smr,"--eqtl-summary",exposure,"--matrix-eqtl-format","--make-besd","--out",
                 paste0(inter,"BESD"))
  write.table(cmd_besd,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  system(cmd_besd)

  
  ##update of BESD
  print("update file")
  cmd_updata=paste(smr,"--beqtl-summary",paste0(inter,"BESD"),"--update-esi", paste0(inter,"new.esi"))
  write.table(cmd_updata,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  system(cmd_updata)
  cmd_updata=paste(smr,"--beqtl-summary",paste0(inter,"BESD"),"--update-epi", paste0(inter,"new.epi"))
  write.table(cmd_updata,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  system(cmd_updata)

  esi=data.frame(fread(paste0(inter,"BESD.esi")))
  nanum2=is.na(esi[,1])
  esi=esi[!nanum2,]
  write.table(esi,file=paste0(inter,"BESD.esi"),row.names = F,col.names = F,quote=F)

  if(threshold=="Bonferroni"){
    #threshold can calculated
    threshold=0.05/nrow(gwa)
  }else{
    threshold=as.numeric(threshold)
  }
  ##SMR
  print("SMR......")
  cmd_smr=paste(smr,"--bfile",paste0(inter,"bfile"),"--gwas-summary", paste0(inter,"smr_output.ma"),"--beqtl-summary",paste0(inter,"BESD"),"--peqtl-smr",threshold,"--thread-num 10",
                "--out",paste0(out,"smr"))
  write.table(cmd_smr,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  system(cmd_smr)
  #generate plot file
  data=data.frame(fread(paste0(out,"smr",".smr")))
  write.table(sum(data$p_GWAS<threshold),file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  if(sum(data$p_GWAS<threshold)>1){
    data=data[data$p_GWAS<threshold,]
    genename=data$Gene
    source("code/plot_SMR.r")
    write.table(genename,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
    if(length(genename)>0){
      for(i in genename){
        cmd_smrplot=paste(smr,"--bfile",paste0(inter,"bfile"),"--gwas-summary", paste0(inter,"smr_output.ma"),"--beqtl-summary",paste0(inter,"BESD"),"--peqtl-smr",threshold,"--thread-num 10",
                          "--out",paste0(out,i,"_smr"),"--plot","--probe",i,"--probe-wind 500","--gene-list",paste0(inter,"glist.txt"))
        system(cmd_smrplot)
      }
      for(i in genename){
        print(i)
        SMRData = ReadSMRData(paste0(out,"plot/",i,"_smr.",i,".txt"))
        png(paste0(out,i,"smr.png"),width=1000,height=1000)
        par(mfrow=c(2,2),mar=c(10, 10, 10, 10))
        layout(matrix(c(1,1,2,3),ncol=2,byrow = T),heights=c(1,1))
        try(SMRLocusPlot(data=SMRData, smr_thresh=threshold, heidi_thresh=0.05, plotWindow=50, max_anno_probe=10),silent = T)
        try(SMREffectPlot(data=SMRData, trait_name=""),silent = T)
        dev.off()
      }
    }
  }else{
    print("No significant smr result")
  }
}
# Code --------------------------------                   ------------------------------------

source("code/All_basic_function.R")
args <- commandArgs(TRUE)
print(args)
exposure=args[1]
outcome=args[2]
vcf=args[3]
gtf=args[4]
threshold=args[5]
out=args[6]
data_convert_fun(vcf)
smr_fun(exposure,outcome,out,threshold,gtf)



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
locus_fun(vcf=vcf,gff=gff,gwa=gwa,region=region,out=out,snp=snp)vcf="../vcf.vcf"
phe="../phe.txt"
snp="4:23943698"
region=1000000
library(data.table)
library(stringr)
library(pheatmap)
library(ggplot2)
library(grid)
vcf_hap=fread(vcf)
phe_hap=data.frame(fread(phe))
row.names(phe_hap)=phe_hap$id
vcf_hap=vcf_hap[vcf_hap$`#CHROM`==str_split(snp,":")[[1]][1],]
vcf_hap=vcf_hap[vcf_hap$POS<as.numeric(str_split(snp,":")[[1]][2])+as.numeric(region)&vcf_hap$POS>as.numeric(str_split(snp,":")[[1]][2])-as.numeric(region),]
nams_eco=c(names(vcf_hap)[1:5],paste0(phe_hap[,1],"_",phe_hap[,2]))
sub=vcf_hap[,..nams_eco]
hap_data=as.data.frame(sub[,6:ncol(sub)],stringsAsFactors = F)
hap_data=as.matrix(hap_data)
hap_data[hap_data=="./."]=NA
hap_data[hap_data=="0/0"]=1
hap_data[hap_data=="0/1"]=2
hap_data[hap_data=="1/1"]=3
hap_data=apply(hap_data,2,as.numeric)

region2=50000
gene_region=c(23944625,23945365)
index=sub$POS<as.numeric(str_split(snp,":")[[1]][2])+as.numeric(region2)&sub$POS>as.numeric(str_split(snp,":")[[1]][2])-as.numeric(region2)
sub2=as.data.frame(sub)
gap=findInterval(as.numeric(sub2[index,"POS"]),gene_region)
gap=as.numeric(table(gap))
gap[2]=gap[1]+gap[2]
gap[3]=gap[2]+gap[3]
sub_hap=hap_data[index,]
pheatmap(sub_hap,cluster_rows = F,gaps_row=gap,labels_col = "")
mapdata=pheatmap(sub_hap,cluster_rows = F,labels_col = "")

phes=phe_hap[str_split_fixed(colnames(sub_hap)[mapdata$tree_col$order],"_",2)[,2],3]
plot(phes,cols="",pch=20,xaxt="n",ylim=c(-0.2,1.2),yaxt="n")
library(data.table)
#library(GenABEL)


# Data convert ------------------------------------------------------------
data_convert_fun=function(vcf){
  inter="inter_result/"
  #plink software address
  #plink="../softwares/plink/plink.exe"
  plink="../../home/software/plink"
  cmd_plink=paste(plink,"--vcf",vcf,"--maf", 0.05,"--r2 --ld-window-r2", 0.99,"--make-bed","--allow-extra-chr","--out",paste0(inter,"bfile"))
  write.table(cmd_plink,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  system(cmd_plink)
}
zscore=function(x) qnorm((rank(x, na.last = "keep") - 0.5)/sum(!is.na(x)))

convert.vcf <- function(vcf.file, which="GT", map=FALSE, snp.pos=FALSE, genotype_file_name=NULL){
  require(VariantAnnotation)
  require(dplyr)
  
  vcf <- readVcf(vcf.file)
  if(map | snp.pos){
    res <- list()
    snpMat <- genotypeToSnpMatrix(vcf)
    res[["gt.mat"]] <- t(as(snpMat$genotype, "numeric")) %>% as.data.frame
    if(snp.pos){
      snp.pos <- rowRanges(vcf)
      res[["snp.pos"]] <- data.frame(snps=names(snp.pos),
                                     chr=as.numeric(as.character(seqnames(snp.pos))),
                                     pos=start(snp.pos),
                                     stringsAsFactors = F)
    }
    if(map){
      res[["map"]] <- snpMat$map
    }
    if(!is.null(genotype_file_name)){
      write.table(cbind(snpid=rownames(res[["gt.mat"]]), res[["gt.mat"]]),
                  row.names = F, col.names = T, quote = F, sep = "\t",
                  file = genotype_file_name)
      cat("Created genotype file from .vcf file at:\n", genotype_file_name, "\n")
    }else{return(res)}
  } else {
    res <- t(as(genotypeToSnpMatrix(vcf)$genotype, "numeric")) %>% as.data.frame
    if(!is.null(genotype_file_name)){
      write.table(cbind(snpid=rownames(res), res),
                  row.names = F, col.names = T, quote = F, sep = "\t",
                  file = genotype_file_name)
      cat("Created genotype file from .vcf file at:\n", genotype_file_name, "\n")
    }else{return(res)}
  }
}


# Display the manhattan result --------------------------------------------
##This function was used to illustrate the single trait, cojo result, and multi-trait result
plot_fun=function(result,out,show_peakloci=T,color_manh=c("tomato","skyblue"),threshold=5e-8,corrected=T){
  data=data.frame(fread(result))
  
  ##check the file format to continue right downstream analysis
  if(ncol(data)==10){
    #For single trait and multi-trait
    names(data)=c("CHR","SNP","POS","A1","A2","AF1","BETA","SE","P","N")  
  }else if(ncol(data)==13){
    #For cojo
    data=data[,c(1,2,3,4,4,5,11,12,13,9)]
    #Note, the A2 col isn't the real A2, because the cojo result file don't contain this cols
    names(data)=c("CHR","SNP","POS","A1","A2","AF1","BETA","SE","P","N")
  }
  chr_info=aggregate(POS~CHR,data,max)  
  if(sum(diff(chr_info[,2])>0)!=(nrow(chr_info)-1)){
    for(i in chr_info[,1]){
      if(i!=1){
        data[data$CHR==i,"POS"]=data[data$CHR==i,"POS"]+sum(chr_info[1:(i-1),2])
      }      
    }
  }  
  ###manhattan plot
  
  #set the point color
  if(class(data[,1])=="integer"){
    cols=ifelse(data[,1]%%2,color_manh[1],color_manh[2])
  }else if(class(data[,1])!="integer"){
    #if the chr format are chr1 chr2 ....
    #I used the mean of chr position to order the chr, and give different colors to  adjacent chr.
    chr_mean=aggregate(POS~CHR,data=data,mean)
    cols=ifelse(order(chr_mean$POS)%%2,color_manh[1],color_manh[2])
    names(cols)=chr_mean$CHR[order(chr_mean$POS)]
    cols=cols[data$CHR]
  }
  
  png(file=paste0(out,"manhattan.png"),width=15,height=10,units = "in",res=600)
  plot(x=data[,3],y=-log10(as.numeric(data[,9])),
       cex=1,pch=20,col=cols,frame.plot=F,xaxt="n",yaxt="n",
       xlab="Chromosome",ylab=expression("-Log"["10"]*"(P)"))
  chr_mean=aggregate(POS~CHR,data=data,mean)[,2]
  ##show the threshold
  if(threshold=="Bonferroni"){
     #threshold can calculated
    threshold=0.05/nrow(data)
    abline(h=-log10(threshold),col="red",lty="dashed")
  }else{
    #threshold can be setted
    threshold=as.numeric(threshold)
    abline(h=-log10(as.numeric(threshold)),col="red",lty="dashed")
  }
  print("threshold")
  print(threshold)
  ##show top loci(grouped by distance) and save the loci passed the threshold.
  snps=snp_finder(data,threshold)
  print("threshold")
  print(threshold)
  #save the snp for locuszoom analysis and SMR
  if(length(snps)==0){
    write.table("No siginificant SNPs",file=paste0(out,"top_snps.txt"),
                quote = F,col.names = T,row.names = F)
  }else{
    top_data=data[data$SNP%in%snps,]
    print(max(top_data$P,na.rm=T))

    write.table(top_data,file=paste0(out,"top_snps.txt"),
                quote = F,col.names = T,row.names = F)
    #show the SNPs
    if(show_peakloci==T){
      points(top_data$POS,-log10(top_data$P),cex=1,pch=17,col="red")
      text(top_data$POS-0.04*max(data$POS,na.rm=T),-log10(top_data$P),top_data$SNP,cex=0.8)
    }
  }
  axis(1,chr_mean,c(unique(data$CHR)),las=1,cex.axis=1)
  axis(2,seq(0,max(-log10(data[,9]),na.rm=T),2),seq(0,max(-log10(data[,9]),na.rm=T),2),las=1,cex.axis=1)
  dev.off()
  
  ###QQplot
  png(file=paste0(out,"qqplot.png"),width=10,height = 10,units="in",res=600)
  lambda=qqplot_fun(data$P,plot=T,frame.plot=F,pch=20,col="skyblue")
  dev.off()
  
  ### corrected manhattan
  if(lambda$estimate>1.1 & corrected==T){
    data$P=pchisq(qchisq(data$P,df=1,lower.tail=F)/lambda$estimate,df=1,lower.tail = F)
    write.table(data,file=paste0(out,"mlm_corrected.mlma"),
                quote = F,col.names = T,row.names = F)
    ##The below code are same as above manhattan. 
    ##If changed the above code,can directely copy to here, except the filename.
    png(file=paste0(out,"manhattan_corrected.png"),width=15,height=10,units = "in",res=600)
    plot(x=data[,3],y=-log10(data[,9]),
         cex=1,pch=20,col=cols,frame.plot=F,xaxt="n",yaxt="n",
         xlab="Chromosome",ylab=expression("-Log"["10"]*"(P)"))
    chr_mean=aggregate(POS~CHR,data=data,mean)[,2]
    ##show the threshold
    if(threshold=="Bonferroni"){
      #threshold can calculated
      threshold=0.05/nrow(data)
      abline(h=-log10(threshold),col="red",lty="dashed")
    }else{
      #threshold can be setted
      threshold=as.numeric(threshold)
      abline(h=-log10(as.numeric(threshold)),col="red",lty="dashed")
    }
    ##show top loci(grouped by distance) and save the loci passed the threshold.
    snps=snp_finder(data,threshold)
    
    #save the snp for locuszoom analysis and SMR
    if(length(snps)==0){
      write.table("No siginificant SNPs",file=paste0(out,"top_snps_corrected.txt"),
                  quote = F,col.names = T,row.names = F)
    }else{
      top_data=data[data$SNP%in%snps,]
      write.table(top_data,file=paste0(out,"top_snps_corrected.txt"),
                  quote = F,col.names = T,row.names = F)
      #show the SNPs
      if(show_peakloci==T){
        points(top_data$POS,-log10(top_data$P),cex=1,pch=17,col="red")
        text(top_data$POS-0.06*max(data$POS,na.rm=T),-log10(top_data$P),top_data$SNP,cex=0.8)
      }
    }
    axis(1,chr_mean,c("1","2","3","4","5"),las=1,cex.axis=1)
    axis(2,seq(0,max(-log10(data[,9]),na.rm=T),2),seq(0,max(-log10(data[,9]),na.rm=T),2),las=1,cex.axis=1)
    dev.off()
  }
}



snp_finder=function(gctadata,threshold){
  snps=c()
  gctadata=na.omit(gctadata)
  window <- seq(0,max(gctadata$POS),5e04)
  gctadata=gctadata[gctadata$P<threshold,]
  if(nrow(gctadata)==0){
    return(c())
  }else if(nrow(gctadata)==1){
    return(gctadata$SNP)
  }else{
    snp=gctadata$SNP
    chr=gctadata$CHR
    pos=gctadata$POS
    p=gctadata$P
    for(i in unique(chr)){
      snp_in_chr=snp[chr==i]
      pos_in_chr=pos[chr==i]
      p_in_chr=p[chr==i]
      inter_id=findInterval(pos_in_chr,window)
      inter_num=unique(inter_id)
      for(j in inter_num){
        #analysis within window
        snp_in_win=snp_in_chr[inter_id==j][p_in_chr[inter_id==j]==min(p_in_chr[inter_id==j])]
        snps=c(snps,snp_in_win)
      }
    }
    if(length(snps)==1){
      return(snps)
    }
    row.names(gctadata)=gctadata$SNP
    pos=gctadata[snps,"POS"]
    p=gctadata[snps,"P"]
    snps2=c()
    snps2=c(snps2,snps[1])
    diff=diff(pos)
    ##filter again to eliminated the situation where the adjacent SNPs were in closed two windows
    for(i in 1:length(diff)){
      if(diff[i]<5e4){
        if(p[i]>p[i+1]){
          snps2[i]=snps[i+1]
        }else{
          
        }
      }else{
        snps2=c(snps2,snps[i+1])
      }
    }
    return(snps2)
  }
}


qqplot_fun=function (data, plot = FALSE, proportion = 1, method = "regression", 
                     filter = TRUE, df = 1, ...) 
{
  data <- data[which(!is.na(data))]
  if (proportion > 1 || proportion <= 0) 
    stop("proportion argument should be greater then zero and less than or equal to one")
  ntp <- round(proportion * length(data))
  if (ntp < 1) 
    stop("no valid measurements")
  if (ntp == 1) {
    warning(paste("One measurement, lambda = 1 returned"))
    return(list(estimate = 1, se = 999.99))
  }
  if (ntp < 10) 
    warning(paste("number of points is too small:", 
                  ntp))
  if (min(data) < 0) 
    stop("data argument has values <0")
  if (max(data) <= 1) {
    data <- qchisq(data, 1, lower.tail = FALSE)
  }
  if (filter) {
    data[which(abs(data) < 1e-08)] <- NA
  }
  data <- sort(data)
  ppoi <- ppoints(data)
  ppoi <- sort(qchisq(ppoi, df = df, lower.tail = FALSE))
  data <- data[1:ntp]
  ppoi <- ppoi[1:ntp]
  out <- list()
  if (method == "regression") {
    s <- summary(lm(data ~ 0 + ppoi))$coeff
    out$estimate <- s[1, 1]
    out$se <- s[1, 2]
  }
  else if (method == "median") {
    out$estimate <- median(data, na.rm = TRUE)/qchisq(0.5, 
                                                      df)
    out$se <- NA
  }
  else if (method == "KS") {
    limits <- c(0.5, 100)
    out$estimate <- estLambdaKS(data, limits = limits, df = df)
    if (abs(out$estimate - limits[1]) < 1e-04 || abs(out$estimate - 
                                                     limits[2]) < 1e-04) 
      warning("using method='KS' lambda too close to limits, use other method")
    out$se <- NA
  }
  else {
    stop("'method' should be either 'regression' or 'median'!")
  }
  if (plot) {
    lim <- c(0, max(data, ppoi, na.rm = TRUE))
    oldmargins <- par()$mar
    par(mar = oldmargins + 0.2)
    plot(ppoi, data, xlab = expression("Expected " ~ chi^2), ylab = expression("Observed " ~ chi^2), 
         main=paste("Lambda =",round(out$estimate,3)),...)
    abline(a = 0, b = 1)
    abline(a = 0, b = out$estimate, col = "red")
    par(mar = oldmargins)
  }
  out
}
estLambdaKS <- function(chi2values,limits=c(0.5,100),df=1) {
	iniLambda <- 1
	optRes <- optimize(lossFunctionLambdaKS, interval=limits, chi2values=chi2values, "pchisq", 1,df=df)
	return(optRes$minimum)
}

# QTL analysis ------------------------------------------------------------
qtl_finder=function(sub_qtl,sub_vcf){
  sub_qtl=na.omit(sub_qtl)
  snps=c()
  window <- seq(0,max(sub_vcf[,2]),5e04)
  snp=sub_qtl$SNP
  chr=sub_vcf[snp,1]
  pos=sub_vcf[snp,2]
  p=sub_qtl[,5]
  for(i in unique(chr)){
    snp_in_chr=snp[chr==i]
    pos_in_chr=pos[chr==i]
    p_in_chr=p[chr==i]
    inter_id=findInterval(pos_in_chr,window)
    inter_num=unique(inter_id)
    for(j in inter_num){
      #analysis within window
      snp_in_win=snp_in_chr[inter_id==j][p_in_chr[inter_id==j]==min(p_in_chr[inter_id==j])]
      if(length(snp_in_win)>1){
        snp_in_win=sample(snp_in_win,1)
      }
      snps=c(snps,snp_in_win)
    }
  }
  mark=10000000
  if(length(snps)!=1){
    while(mark!=length(snps)){
      mark=length(snps)
      pos=sub_vcf[snps,2]
      p=sub_qtl[snps,5]
      diff=abs(diff(pos))
      ##filter again to eliminated the situation where the adjacent SNPs were in closed two windows
      for(j in 1:length(diff)){
        if(diff[j]<5e4){
          if(p[j]<p[j+1]){
            snps[j+1]=NA
          }else{
            snps[j]=NA
          }
        }
      }
      p=p[!is.na(snps)]
      pos=pos[!is.na(snps)]
      snps=snps[!is.na(snps)]
      if(length(snps)==1){
        break
      }
    }  
  }
  return(sub_qtl[snps,])
}
eQTL_para=function(x){
  
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 100;      # read file in slices of 1,000 rows
  gene$LoadFile(paste0(inter,"omic_",x,".txt"));
  output_file_name=paste0(inter,"omic_qtl_",x,".txt")
  all_names=intersect(intersect(snps$columnNames,gene$columnNames),cvrt$columnNames)
  snps$ColumnSubsample(c(which(snps$columnNames%in%all_names)))
  gene$ColumnSubsample(c(which(gene$columnNames%in%all_names)))
  cvrt$ColumnSubsample(c(which(cvrt$columnNames%in%all_names)))
  useModel = modelLINEAR
  pvOutputThreshold=threshold
  write.table("run eqtl",file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    verbose = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE,
    pvalue.hist = "qqplot")
  write.table("finished eqtl",file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  return(me)
}
remove_linksnp=function(x){
  source("code\\All_basic_function.R")
  all_qtl_filte=data.frame()
  for(i in seq(x,length(genes),5)){
    sub_qtl=all_qtl[all_qtl$gene==genes[i],]
    sub_qtl=as.data.frame(sub_qtl)
    if(nrow(sub_qtl)==1){
      all_qtl_filte=rbind.data.frame(all_qtl_filte,sub_qtl)
    }else{
      sub_vcf=vcf_data[vcf_data[,3]%in%sub_qtl$SNP,]
      row.names(sub_vcf)=sub_vcf[,3]
      row.names(sub_qtl)=as.character(sub_qtl$SNP)
      sub_vcf=sub_vcf[order(sub_vcf[,2]),]
      pos=sub_vcf[sub_qtl[,1],2]
      sub_qtl=sub_qtl[order(pos),]
      sub_qtl=try(qtl_finder(sub_qtl,sub_vcf),silent=T)
      all_qtl_filte=rbind.data.frame(all_qtl_filte,sub_qtl)
    }
  }
  return(all_qtl_filte)
}



�PNG

   
���DTA�$�� �ePQ�I� ����b�5Q4&&�t��M�����t:�ｫKR�|��[poթ:�^�Y����N�󟡬���1�fV��^�p|X.��;�=�������O���/�߇��5t         �![�l��=�)[�l��M�F)[�l��]ʆ)[�l��m��鸐�S�Oۇl���2�c���ö�pR�$�
~�%T�#         `.�6*�l����v*�l���ʶ*+c�1��L�+��¾!����s���_C�C        ���U6W�^e��-V6Y�f�K��c�fV;������`�u�~x        ��d���V�[�pe˕M�1�3g��a�pZ�-�$T?        �=�ve㕭W6_�~c�1#7����ῇ��2��ʫtk��N���[t[m��n�����|���Ao>�{�q���N��;y�5��>tCw��>՝y�g�s�������u\����]y�w�k��Aw����x�O�[������        �W�F�e{�
������V�v���n�%�-?�        `��v.�l���Qd���^6{�36�s�)T��M6߮{���n��g�0         ��l鲩˶�j�z���3Ƙy;��B�!�x�*O�������+�)?p        `�ec��]6wU��(��˦�c��~����[��{��ݭ_�y��
         ��.ۻl�6�Qdۗ��1���Y3��.Tr�_����
��#�������I�         ��5��j�d��
�i��?����/@          �lI�)�ZӐ
�)���]���O          �W��fK;�׆ln��5ӜCsA�y�G�'           �li��6d{k�1ۅ�B�˞�         �q�Mm�چlp��]a�>�I�u����         ��Ȧ6����6d�k�`�<��//8          ��lk��6d�kc~.ܺ�o\^h          W��w�![\�3�����+ˋ          �*۪�
��l�G~�8�����          �#���7d<������^�@�          �ّ�o���l�Gv�
wx�w*           0�6�f��wR6�#9��w��S�/          0����C6������WyJw�^>8          `ve�ʪOh�'e<r��0pG�����          �!���8d<R�Ch��W�S>(          `4d�[��!ᑙK���d���          ��M7�~�����H��῅�;��?T>          `�d�;��l����92ܹ�+���t���          ��lWXq��&xR�³>_	wl�}-          0���C�³:��-���           F��K�h��I��ڼ7ܡWZ�|           �h[q�'�������N�C�ujy�         �і-�p���yAh��E��jy�         �і-p��l�g|��;��f۔w          ��	�C��3>��a�xy�         ����
          �M�w�![�e>?���w���          �ܔ��p7�%^��Uh��}���          �ܔ�p��l��ټ'|�-�٥��          �ܖ��p?�)^f��0�
          ��l��;���f�0p#�-�|w�^~S          `<dS�m�po�A^�9<���|q�
          ��9��������}ٝ�          �l��;�mr3�5|�U��~y�           )���9d�<0O_��J+�7          0U���=r�F���6|�:�o\�          �T���!�G�0���kyc           Se{<�#�l��������          L���p��Q~d.	_p�ы�          �*���9d�������3�(o          `�l��{䐍�#�0���Kyc           Se{<�#���#�'a�λ����           �:��/�ȓ�Q~d~���O��1          ���=�C6ʏ̿��/�����7          0U���=r�FybV��+��JyC           �l���䐭�
+��hi�|���>tMy��s�1�/ߣ����������+ L������Y3��   0����C��  ��u�g��m��v��s�j�?���?*������Y���uz�}ۉ��f��{u����u�<m�gt�`�	�l���[l۽p��-�~��_��z��v����n��v�^�˫���l����}���y�˻-�٩�l���
�_���~@��ޯ�y2��9X����e   FA�;kh��a   X�s�E��s�>o8�||�����[~���5z4/��U�m�k��^T��G�(���o���k�C��<
�{Uy�`�AozW��d��k���c�_^��,  0
��YC��  ����7��n:���*߸�a�����hNXt~y;��z��/�ϣ`T���:�m'��{�_���3 L��S�-�̴u�{~y����~��
+����   ���w��.��   �,��ʽ��O�󮼣||���������������Q�nG=XN�~���>�VZ�������>�o->�����L�`����S^���   ���w��.��Gq    f�����z�ۨ�u�Y��Ϝx,�cW���嵩<�Y�)omǞrN���rq�G��zKy�g�3�����ۿ^�g H�_}g�3d�M�����N�����0��|�<   3���5���0   ,k���@�;�\���w)�8Y~�����a�����h;u���{z̕`9m�­��0�6��E���t��[���i�<w����G7���n�'�Z^�a�e   FA�;kh��a   X�>��O��Q�}zk���s7ؤ�.�|ӱ�y`�-�p4�S�,_~Ӄ����c����� �}��ώ�6N��Kwۧ��2   ����5���  `��>����{�~�|��`�W���&�������_�[R�i��K���.��;����Kpn��=�k�@��'��gf��~_�����C],��^U�����9n֧�k�{w�kպ>��v�#��?�]�ׁ���]u�E7�}t�m�30m����-��
1�xs�%��������   mJ����b  �����r�^�{�8��x=��k��@�=4��g����\�y���^m��׍y��u���Ba���Y�6�q�_Ga  �.H����b  ������|,\��ٿ���gg_��x=��;���@�=�����.���#��8�?����iף/�4��hےˬ�'�o�K�{/
�   tA�g-��0   4m<�+�l�E��xv�=��k�Uw=�v����W?���`˕-��=��|�͘�i��r;�;����o��  @�{�R}1
�   tA�g-��0   4m<�+g|�Ƹ����I?��૾��o�<�m��o��.��re��=uAUNK���t���Բ㳰����T,���qϣ��  @�{�R}1
��� L[������r+M��+�e   � ݳ��i   �6-�+g^ts���̳������7�ݧ�ܜ�_�4�+�|�͘�i��r�m��vq�SBa  �.H����b  �������˟��atϓ��k1(,����ܬ��_(�N8��iGW
��,�r�7lv��฿)��  @�{�R}1
˕=86��a���[��u/f�y���/),��RXn�ѧ\������r�q�Mq_SKa  �.H����b  �����\Ye�:�5�f��}�3�4S���ax),7�{_�{���y?�`ڠ�<��{�qO���  @�{�R}1
�0���7����u�r+M���6(,O���x��k�y�^�Ia  �.H�����{�  �m/����f+�^|��Y�m���M����k�8t�+/,���fɥW����X8f`��XG
�ˮ�j��E<�a1����}�ۣ/�$f   �6�{�R}1
��;�s����+��RX^~��1� ]v�C1k��  �t�Z�/�a   hZ����K.�3l.��޸��J{���2/����s�	q]q�M��� �
��yW�s�Ma  �.H����b  ��5UX����g��U�_����m�9�c��k�7����(ξ��⤳�*?ἑ�r��a�>_y��i\W\|����w���ɷ��^�E�����	�Da�y��󐸿��~���|��>�]q��o���pqΥ�+N8���c�.�<�ؑ�����Ï?wd��y^���?��x,F�O�4����G^+.������-�9��⠣��?���#O+�=���Knyo���V~?�|���q�_:��R�nN��.RX�;弫c�APX  ��=k����  �iM��x�x�a��bK��uQ�?9�屩
�U�j��O)6����2˯R�8�L�ڎք	�+��V����wM�^)^����ar
���z�=�����G^�����${��7;�z@1q�
J�<��S�L��t�C�_?kұ��������A{�ş������r�ش;{�Xb��c�A�>�{�ï�"]   ��{�R}1
   *�>�T_L�   д6˕��|.�l��q���NnZ-,��؛�ĵ6��dXU_������e|RXn�7=����N=�ڡ���ꉲ]s_��0��gF�`��8����fq�����K��v]�~��7�B��]���[[������ES�i]�^   PI����b  ���]X�~���J��85�o�R��UO!Ls��)�]�E���~��ބ�6�j��������ܬ�N�N�[W<���1�x���,6�z�x
�)��z�ǜzq���u�y��*,���V�=3�m[UXN����7�l��y��t�   ���#K��4   MTa���B��s
��7���<��Gq_]q�����ǃ�CU�5�{Ju҅�\���/��Xw�-�q�T���ɺ�|]Q}�'e��.�Jq�C���Ʃ�];�w�t�~�y�Y��T} *e���{�X���?�l��?o�����J���M���v:Y��   PI����b  ��
�i�������ҹ��Ig]��/�wn</�Ia�9]|��:���y��|�#��Ҟ��Z�m�5%��@UWK�;�~P�;��'f��L��<-�oj�s��z�o�6(]�^a���߇�Ja�7���K�   *�>�T_L�   д.�+���o<ﰹ����ڒ2Mn����x�o��9�<Ť�~�=�6�j�x�~�
^���ft��V���1�0���3�<K�0a���W~�7��(QV�����u<� �u��1�Za�5��?��{&w��O�,�d<��J��V^=�������?��f�#f��ny<���\X�u�����5  �J��,��0   4�+���2˯�=l��ޱHy&7�,�ƹav�%�Ľ��W���O��yߟ2�U{t|<7�Ea���wM�Sl��>1�0k��\i�[V[c�x�����[�s��g?,��0!f�RՓ�ӹ���Y�����A��~�_�<m�Ra�k��ԥ'�+,���jk�|Sc�]�/�g�)2��?ǦT�֔t
��Xl�e➺��g>���USe�y�[�x�g��ZW��X<����Г��ڴښ��lSj��7��i�AG�s��+�>�m�\zk�7�(,W����U�<3�x.,O�Z4�l[���|   Цt�Z�/�a   hZ�
�O���x�a��V;��5-e��x*,?����9�'n�����6������Ӟ��pPX�s/�^�O��A1�z�{��)�^��/U7�wj]p��|mh��]}�A:W�N:�ʘe,���/��V=�:�����/flR�t�e�{SXΞx�Ә   ڔ�YK���Y   @�^�`���0۔�n)ð9攋����rLn����[���Xv�U��i��&��7ᤳ��t5�]q{<?���:\X>踘��[q���.x����荟��Xa���>���˭��O���`<�Ԛ0a����ٴ���4f�Ͻ��x��s��1�h=��O�q���O�-��K��?fl�#���w��D�8h��ԑ���ub�A�La��Oc>   hS�g-��0   4����JU~K9��m��ה�ar�W�a�����뷑�v��7������m�9c�~�~����0f�����Ko�{��ϸ<fV[��{�g?\t�=��6ǜs��O����1��I���l�25���Î;'��;{#s��{�h����'����ro
���|�  @��=k����  �i]-,WV_{Øeؤ�5%�r㡰\}�~�[��L�Д��;2��UV_7��nSX�ꉪ3�0c�ˠ���&1�:��3�>��zBv:gv�爘�μ�xΦl��61��Xd��ڴ�.��l������7�/�l�ضK� �k��s����۞��Ma9SX  ��H����b  ��u��\�{��b�a��R�Ž5!�r�^X~�՟��6{�[����f1C�n����	��{D�@w),��Z�m��w?�V�<�ο⎸�~9䘳�y�p�/��P�^{��_����-��3�C�3*��MS�a��nx(k��a�A��]a���m
˽u�o]),W�))   �)ݳ��i   �����+�>�6�m�S�_��sOn��l�m�W�<����i+��F�ӄ�/�%f����㠣N�{��.�.fF���m1�\��}�˝���݄7~�o1C�l��y�m�o���U�;��MM�Q��˹���5[l�K�ضKox0�kZ����۟��Ma9�Ja���~�  @��=k����  �i]/,W���N��͉g]��O鼓�i����08霫➚�^�9����'�<M�a�����u�A�(,O��.�)�<6fVU7��_Yl�x�&����1K���������7WT���9�V�}+��:'���x�Ah���h
˕���!f6io���7�a-,{�C�~������mx��_�LM���7c�Eay�]tͽ1w�������c�=�&��7Ŝs���O7�=)��I�^������|��Ǟ��O{||<w�����/����1a����6����;��Ma�����9ۦ�  @�{�R}���   ж�FYXN�������|�d�eV�{�t���8�Lq��n�����&�{ȉ1K[�~R\U�N9��:^XN���ĳ�����*+���?���j���{��9�7��i�=�J��O�~k�x�~Xy�:��V]�t�����9�'g��A�e��bƶU����5m��<a��;�����S�AYh�%bζU��   ڔ�YK��4   M��r%�6[�G���J��0�'��A�K�N;�ژ�-�l����IW��H�Bw(,�M���M��)���XV~h��q��������7���~���3��?5&���x�&,���ū�>�hۊ��3N�ɒ�X���Q1c��w��<aZa�7��La  �.H����b  ��
�]y���w>�
˫��^�7(/��   _J����b  ��
+��=5θ��x�&���Fŋ��}�Ӧ�6�2����l�AaYa���La   �"ݳ��i   �6���t&��A����Ƽ�&�mJ��On����<��+W���b��TOxN����2+�<t��rv�3�ǟqY��bK�l]Q�n�{R��x��ӕ+��PI[����7���T�T�t��-��
�}O�3��{��پԥoPXVX�Ea9����e   � ݳ��i   �6���	���.Xy��c�a�ϒh:�䆩�<��+Ͽ�"1O��rW�ֆ��x&fb����o�����]�:�Xh�n<U�l��.#�o�~ƃ�_�i�wf�aƘ�-�o�u��o�||<����c�y�0�,�W��X�Ֆͷ�9f�,��Rqf��{�Lay��c�A����mSX  ��=k����  �i�),O?�q�+f�}��{���~qoc��=�a),���+O\'fjӤ�~��a�v���i������P<<�G#O?�+�o�yH���s�=o<WUO�����8�r̙q�MXa��c��l��1W���D�V�^N>窘�
�U"�o��[�3����3fk�so�*�b��\X��ħ����o(ξ���/�n�x|�)�nq�Q��rb���G;�~`��f�+����Ϋ���s�������>,\x�x������ڲ�އ�\Mx��Sb�
�CPr����b�aS}�~��h�cN������ϼ��e�����6�'�u�1����26��G�x~����������ע)Ǟzq�і�8%�jB���aJ,���x�AXc퍊G��$�l�m��T��AaYa����'��A�7(
�   �鞵T_L�   д��g�i�8�5])yL����JǛ�0��?�䘽-qj�նZ4�kCU�H�,���Zg�-��\}O|�ƻM��V�&M�>��r���S�s5a��w����l�1(3�<K����6i�ͷ�y�kAaYa���l�%��9ۦ�  @�{�R}1
˳�g�j����&Ua#�훤cMn��g�s]QdR�6��1[�V[c���-;�~P���(,wC�D����۟�)�OӒSν:^�&���k1K[��>�j�=O�s�����+��V��9�þ�����+��m
�
˽t������|���   ��YK��4   MMay�Yg��]��ǿ��6S�U��8��zay������T=}1ek۠�^�S.Gayp6�|��«��˴��.�z5�I�,mi��5�_p]�1VG�x~<~�T����κ�ڹ��;�׶MaYa���l����9ۦ�  @�{�R}1
�
˽(,gK-�b��6�e   � ݳ���/  ���<�������ar���Ľ
�W��H�״];TXN��Ɏ��\w�oP�RX~���c>   hS�g-��0   4m4���/���/�.�o�L7aB�ۗ���i~PN:����m���S�A8��cƶm���1��ray�CO��{9���㱺���{��{����4m��y�4���[}+f�������t�MW�q��q�²�r/
���ˮs�Ma  �.H����b  �����<ϼ���a��=�{&���n�[%�~r]-,���1o�Zd�o�9�;1c۪�֥|��Ӆ��N����F�m��eU�3�eZ�����kӴ.|�d��V^=f�ՇP�9��6�ڳx��{f
�
˽(,gK/�r��6�e   � ݳ��i   �6���|��aU�G���a��a'ǽ��N���媌��K�ܮ�+��1#��ra��~&}����m1�L3�cv��l�2-Yt��i��3�����>�ٚ4����,Sⶇ^��F�/�Hq�-��}+�e��^�SX�$��e   ��t�Z�/�a   h�h
��ο`�fU	;�u�\s��}��M��*���ZU��r��m[y��c�A�Ra���?�iW��>e��������e'�qY�˴���>^�6T�+�_`�Kٚ��O�5�Sb�-w��VS���.RXVX�Ea9[fy�e   �R�g-��0   4m4����_(��g^�E��y����}��L������B�
]*,���_�ra��#N��Gk�������^�6��w?�f�4��~ߏ)q��/�s�U&�S<��{q��DaYa���La   �"ݳ��i   �6��r���iv�]r�}q��d��jO�5��ba��;��Ya�Ͷ��K��*K�H�:]X>�Ԙy,[b�x�Zh���>ƻ+nz8^�u�C/��cJm���<�l�fyb{��PXVX�+����4��eWX%�l��2   ]��YK��4   MMay����������<L�\w������'����^���3B�
�;�_�H��\X>��b��ރ/�cwٞ�2��|�U�ZЬ�o|(�S�ƻ��3z�Yq��@aYa���l�W�9ۦ�  @�{�R}1
�K����E�~K��U&�s��3B�
ˋ.�T�H�:]X>�̘yJ<�ʧ�]�ɖ;ƽ�G���5�Y�zb|?��5�=�7^�1���uw<��E
�
˽t�������|���   ��YK��4   MMay�Ŗ����M�L��o���*æ����a'Ō�Х�����3Ү.�=���yJ�s���<]VeN{o����Ӭw�?�����G�s�'�^vk�{�(,+,����-��Ę�m
�   tA�g-��0   4m4��i��G�|a��I�� Ϳ��1� t�i1� t��<�L3Ō��˅�Î;;f����w<W�=����^ƓA��^o㭊K��o(\z��}S=	9���֧(�Zf�x�ǓcO�8�K��{QX�VXy���m
�   tA�g-�?���   �{�����Q�j�%����Վ����x��<H�׵������e�ٶcN�(f�w��1'�y����M~�91��x�g(��x��Zn�U�^Ɠ9�'�-�n�S���y���-�x����~�����),?�5�;���b�A{��Oc޶����1ߠt���E�   mJ����b  �����������l��׍�b<H����)� u�1� t����G��9iO��ǟ3O���x&���=�����b�Yf��nˊ��s1���b����}<9�3���@aYa���l�Uֈ9ۦ�  @�{�R}1
�7��\�7h]*,�|�ҕ5=��1   �)ݳ��i   �6���Rˬg�{xl�.�*�q�R�A�Ra��cΌ�w�6�=].,}�1s?�z����]vɵ�Ž��>�}�g�|���~��b�
˽),g
�   tA�g-��0   4m4�套])�NK6�d�xm�Q�� �<�,1� t��|�Q�ǌ�����)�=].,s�wb�~�l�����lқ��{F�|���Ƕ���k1��i���0^����q߃����܋�rV�N9ۦ�  @�{�R}�_�   вW>����2˯g�5�,�d�>�&�m��{ޘsF�2���3ʻ?��������.,�򝘹�^����>G��Uk��I�˰�n���>�t���l4��Ύ��xqБ��}�m��;TX����^*,�|���+*,�|�ҙ�����   �M鞵T_L�   д���]a�8;���#t�>�&�m�Xhјs�?���q�=䄘q&L�3�.��?��Ƈb�.;���^��,���ئ�ξ2f�Y\yG|?Ƌ�^�q�w���{QX�V��N��6�e   � ݳ��i   �6���r+*,�����h��}
˙�2   �E�g-��0   4mt��qvZV=�/]�a��3H��W�9g|���qv����q��Aa�����c�.�����2L�]aո�6�8�L����T���3�V]c��׶(,+,����UnSζ),  �鞵T_L�   д��WXy�8;�����zu]�� p��1� �s�-1� t����=�i��r�1�^�t٠�����z��}����_��h�M�L*Yl��
�
˽(,g
�   �鞵T_L�   �4��)W�g�y�xͺ,�e�N=��sν�֘q6���1� q¹1#�RX���h˘��^z���^������ԶÎ;;�}�~��N}�fj����q�m�����R��),�ֵ��jk�s�Ma  �.H����b  ��),O����Y��}ҕ7=s��\3�F�o3B����2���7~YL�0!��Ͷ�)�e��qOm[}�
ˣs�1�(��ŵ�?��U��>g�� <��ocF���Ύ�_�UO�N{iґ'���MaYa���La  �.H����b  ��),��y���cW�̃��cֶ��ہ1� �<�,1c۪2y�G��Go�5׏9�j�9�.^����^�h�5֋�����2�=���B��r+�����ǝvq��6�e��^�VX^c��bζ),  �鞵T_��|1   ���Q��,_�#N�ײR�.�w�c�6m��=c��}���|����ǌ��S�Ay�Ob�.�~�}�^���pG�� ,��b1#���/���s���{�5o�}�C�N>�ʘ�mUa9�kZW
�7��\�7hOu����
������   ڔ�YK��4   MSXn��[��砥�]�÷��y۴��;�lm{�G���ڶ�2+�|���؜q�u1k��}�Mq/]s�%7���r�5�Ĝ����:q�
˙�2   ]��YK��4   MSXnV��'L��려�]���~�i��6�����~�m�E6����S1�����]5,�ۮ]�{�x#��/�{���u[o�-�k�t�Y��ru�����6�ު\y�#1_��{SXΪ����mSX  ��=k����  �i
�ͫ�R��려�]��Jc涬2q���mw=�z�׶�������PX�2W��X��e�{V�K�l���1����!1'���?�����~_ߦ;~�Xf��c�Az�?ƼM����b��]v��1_��{SX���x똳m
�   tA�g-��0   4Ma�ǝzq�����u��ǝ3�e�E����v����|m��S7�k
�Sn�������{&Žtű�^s���C�:��=������6��)U/]�#N87��m�~���V���O�%fm�u�?s������5Ma�7��l�-w�9ۦ�  @�{�R}1
��z�'���k�y�unK��%_sO�ݖ��a�զ�O� fk�zms1X
�S�3.���lPO:�}9>f��6�.f[l�K�W�A����JǮT�M��+n|(fm�K-s5��W?�Y�6�oYPX�Ma9�ʿ7���   Ƈt�Z�/~��r    Z6��r�e�����:�%e�A��/�ၘ�M{xL�֖3�s}��`���C�.�p�m��j	�S�~�2h7�7)f��Ϲ2���s/�%�g�%������_���7�� ]{��1kV]}ݘ�I]�9_}�#�k��{��ƻ�����]��m���
˽),g�zB�ٶ�{?�  �6�{�R}1
?��3c�A[{��b�.����Y�<q����kmg��s��󥫿�x�����+�m�i�_�4m�E��y�T@S�&=��Ͽ���6��޺VX>��n����7c>   hS�g-��0   4�|\�G����+�2�6�l�x͛�rt�^�7m�mw�yڴ�2+�lM������S���u��\��R�.�e�C�^�춇^�{��d��v�1�|'f�w�?�\uˣqfJUBK��RU�Ns]��ST�|����i�L\'�iӮ{�5i�u6�YAa��������m���Wb>   hS�g-��0   4���X�GM>���y��bK,�ySR�.T)�e�g�)fkځG���
˽u��|��WĜm���I1   �)ݳ��i   �v�=߯ݣ&����<S�gߋ׼))CW���AqM�w�c��Lz�W1W�f�a��w3�
+��߆���)S�VZu͘�	���1� ),�ֵ��Wu�ߡ��x*�뇻{=�  �W�{�R}1
+O����E]r䃍i�   ����R}1
�u�����U�A����c�A���{b�A��W�5�a��y��x�eq   �K��,��0   �K����L��I�F���t>������c4��ƃmv�=�_��d�'S�sO����5��_6�l۾g�]�sf|o�l�V-^~���~�����u]*-�{�-1�b�Yg����J��/��uH��7����������^��E;�_�ð�qƙ��~���e�[9�~P�x��Z�Ax��b�A[i�5c�.x�城|@$�{xL�   |�t�Y�/�a   ��+�ӽh�����?�}<7�x��O�駟!��$o�h�鋓[m���V�j�|T=�;������������{;�_i�≗~�6̺��ϱ8�ԋ�~��\3�w�η@�O ����9��⭟�K<�XUŽt��q�a'�su�+�}�Բ+ƽ�9朻����>�[q�83(�7x��m��W�]l��1_TO�O���۟���n����   ���5K��4   S��{'����[`����9hƔ~Ex:�xr��'�}�˂/�
8ZU��zZs:v�T��tn��-�Mj��m�{���ˮ�?�qXU�Ӵ�aP}=~W�z�yWǌ��K-W<������f�9�췹枷8��Ӌ�y7��&���t�����cO�-��5������H�7�iX̿���}O��WYq�5�� m���1k&���b�Z�\]r�7��]p�9W��]������   L���f������?  `����-N���b��&��=�P��N=��۟�[�G�p�e�}�%g�9��+���i�=*^|�o�����mq��&L���g���������{�|����c�{;̶�n���۟,>��Ӱ�as�)�}�mwګ��*.�������<��3�'���a1�_�k0(�ο`�ڤ��՟���y�x�*�\�G���8��u6�<gj���������&�}
�|��0:U�/]ׯJ�ӊ�ox�Xq���u&����p�3q��h?<�����K�:w��v�?fvi�m9����I�����x*�k����q_�w��=�;b������{ō[�{1a8�EG�t~m����/���O̮�A%��&c)�w�V����K�w�t��zu��P��^b�ecΦ-�Т�EW�s  ��J����b   �e��Kn��vrinZs��K.�|�>]��*k���iOt�����?�L��|�ł/s���6���.<��b��f���d�g)N>������y��^*���W�Ə�)�U���ߊ��{g�����.�h��F����W�sF�m�
�=�~�8l�u<n�m��.q?��2q�x��K��ʇ����s�*�Zf����_r���.*���c   �tOZ�/�a   �����#�#����������R��V]��#_㝮S��˾�m��@�(,�W�˕[�s���Ax�*N���b��W�9i�u6*ξ�Ƒru��eUi0��*���������.V_{�x��m�շ:��I���[a�xm�V�����1�Xl�yw��LIa��O�nZ�
˓���m�߿ͷީ���'�y  ��ҽi����   �D�'��O�^��Ͼ�Xa��
��=����U���mͨ�.R���{ꭑ"�y��Z�s�M#����.�>�f��+?k��   @3R7�T_L�           ��nr����          zI��R}1
tR�����:ڄ��L���N�S��c��JI-5'Ҝ�i�QF'T��$&�����"��.3�;{���<�|�}����򾟳�g��������������(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o��          J���d����a          �Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'���t
           h�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7�x�0          @��69�7F�           %Qmr2o��          J���d�          �D��ɼ1
    (���|����7��eO��7��ٟ=�۳��(�T�\��}���t�C��O6��=  ���  ���69�7F�  \��ru�������q��/z�{��}�r�e�3-����._4߭�:G��|C�i��5��~��F����.�m����Ԭ�1�� Q^ n��gA�(?�+Z�V~���j~���$  l/�MN�Q0  �m�b��mY@�o�	�[�iQn���p���n����������֚�V~��ޚB�)�p7D?zE��\�����ݳ��  @��69�7F�  \��26�����s�}Y[�4�+�_Ëg�.kϱ('��~�lkF����Z��??��3��k��!&��������˲�'�kZ�n����?׳��1"Q^.���h탒��  `��69�7F�  \�5/��Emι�K�qB�(�_L�u�=o�\\?������_����F�G�?F����w�ﮭ��E9��=r;K�xtԞ��L �ڟSQ.�G��K�>w���]�E�  �vQmr2o�� �n�/h��/�{s�����,�e�-��rq�|�x�z�F�F��,5��߷��n)��1~�<���!���vͬ����Ż��s,���3�  �E��ɼ1
 ��~y�E����k�(��/&�zE���(�C��K�<�����ڟ�Kz�rL���Q��W��k~�5�g��#ϣ��(��\��}�g���%��,���z�}�v��o���=��\  @��69�7��+ �[(z[��� _�˞��a����+Ys/༢gzI���!Z�%Q~�������%��A..G�n�h-{
JZM�$�Q���uZ�@}�2�l�����%��,���D�Z⼺l[�m��wi�  �E��ɼ1
 ��E/`k�� n)D��� �U�	�^CAD�(�C��K�<>[[���>�`���|��eo�@��=�D9J��N4���lk�����%�����k���m��i�� �  l'�MN�Q0  �/z	[��qϗDYO�R���Kg�^=gK���!Z�%Q��n���e����h-[�z�O�ލ��=�_�F9��E����K��Ż)Z������o�q��. ��D��ɼ1
 ��E/ak�� n)H���D9J|1��ЌE�$�s�E�[})�Z���b�r��籵�>;5���r����\%~����3<�Q�7�V�p���v>r	Z�dQ�K��%Ϋ�����t���  �E��ɼ1
 ��E/bk�� ���h-V�r�x���ЌE�$�s�Est�/�
�/Wm����Z.���i5�h;���=�_4�%Q�l�w���0��g��K��Ż)Z���e��j���j���7  `{Qmr2o�� �~���=/�{�(�Z_G9J�|���s�Dy���^幫�=O[})\�����z���2�i��ߥ���Mcr���:�W#�!�CY�y���eǞg�#�`���F����k����{v���.��8�G�w��H  �GT���_�> ��3~!�"����L_ ���ŗ��Q����$������(�]Tz�z~�G�4����\�_��{-k����}z΃(p^�g���%�����k��겔�ۭ�6  �'�MN�Q0  �/z�[����"�(��p�|�X��K�<w�R��V_
/�'�`�r���
��9�<�y�~8��^���u-q^]�S�m
  �OT����`  �_��F�b���|�X��K�<wIM��V_
��kJ���j=,3�sDy��Y�9�|��{�n�ֵ�yuN��)  p>Qmr2o�� �~���=/�{
-_D�%�����hƢ�^�+j.o��p��
�/W���`���� ��G�� �#��ߋwS��%Ϋ�;�ߦ  ��D��ɼ1
 ��E/k�� �)L�|yŗ�b��/���{I��.h���՗�=�(X�\�篂e�z΃(pz���K��Ż)Z���y��oS  �|���d� p����5z^�&Z�\��_L���4c�z/���v�gj4�l��p��
�/W���`���� ��N~n[�~�#��ߋwS��%Ϋ�9�ߦ  ��D��ɼ1
 ��E/k�� �)L(X"��f,Z�%Q�۬���K��,_�ֽ�`���� ��F�3�|��{�n�ֵ�yu=��V�  ��&'��( ���������k�X�����=���c��=��ÿ��{�oV�c6��9�:�_�_Ql6����|���~�;��eü��D�k�gI�glX�i�0�s��0�Y4���Y*C����^s�kI��6�s}l=�l�N9O��di/�D4�ܾ���Ǟ��߷~N�u?�{-jE}/����°��=�{��(����Q�=����>���ݱa��]ҸJ������}����8��}͹��.����]�KY��
 ��E/k��8���h�W��d�/&�^8O�Q���^Ζ�_(��ޫ����/���������u:�9ǥ��)�F���Ҳ���Ai~�Y?����Z���-}~�~.���s���s=�k�[�9&�#�c��1k�y��/{Zc��(�d�o{�#�>m%���y����|�J�~��7[��%=�^��d���r_����~X�|������oK�c��[^�i�����ﱴ�%����<K�q�����Y�Ec������������2������Y���Q?N��b,ǮY�;����c֌瘚~
 ��E/kl�x/QKj_���Pm}�?�,�}ˋ�u�yQݲ?Z�^�:��\�:�u�{���`:ǹ�w�{�ŗԞ'Y�|G9J��R�[�e�<��"�i}��<��<Y�>���$�����|N����e��[��кWZ���z��2Y�1���DyJj����(GI�#�v�ҳ�={�u�fQ����ғ?�r�)�CI�^왇=��A뾬ݓ��Ӻ�[�},�W��^�p�k����:�,�s�)�j��)�[�<��s*Һ�%Q�H�^�bKj��=�ִ��������5�io�2��9�sV;�(��g�Z�8k]�=�'k�� �zE��ɼ1
 ��E/kl�x/QKj_�F�Ǵ�pnyI_3�5/�kװ�Eu���r�h}	��s�u�ZK��s�qj�/:����=�:�Q��(G������l����kI��t�]¾_�u�e[��4�Ǵ��{�-�?���YK�/e���kψ��� �SR;����(�rD.��zֳg�Gy�DyJ�ޛ=�ϱ�Q?JZ�x�3��9�z�e{�k�27K����(_ɚ{�ʵ����C�X�(�1Q|ɚ��{,{�9���9��Z��(�rDj�B[�G�=��3��1k��-���{9���O�>��>׮�g  �/�MN�Q0  �/z	Xc��{��[��{iy	�?;}��j�/��%rֺ?z�s̚��;��/-�����r�%�{�ŗԞ'Y�|�)�_}�����%�גi�����
���ܖ/&Z�^@��k��z���yq��v=n�:���Xr	���/Q�c����9��KN}�>٩�8�iI�\�p�u��y����c��?"����9�{��<%��bK�}��g=[��L�٧=cɢ\�S���~���u�<��{�iU:c��*}.ҳ�Q��5��$��{?�c�j����l�\��%�sn=k^�vݢؒ���ZQ�Jj�X��]�o{�<_�\�s���y�9C��%=��9�j�?�k~=�?�> `{Qmr2o�� �~�K����0�o�_LdQ�K�׮�rx����E������z�][���Jj�to={%�sL_�2'Q|I�@[�y���!���.��a�����e��9���Q�ۤ��]�3c��3��Ԏ+�-���5�Y�ֳ���٧=cɢ\�S���~���u�<�>G={f�.�!���~�-��1={|�+Ys�Kp�k����:���Y�y���U�X�<��g�kԬ]WҲ�����f����(Ϟz�X����\��y�j�t��g<=�Z�3����� �u�j��yc ���^ָ��QKn�K�-��ֵ�E����K�^�߆u^#ߒ(�)�|��9&�/�9OQ|I�@[�yN�ǩ�O�z�;��r�zƹř����-{jϱ����&��o���{��<%��bKz�~Mzֳ�,�G�>��>���-D�(i��yh}���9�d����+�G�?�g��|%k�uNװ{?�c�j����l�\��%�sIzƴ$��XSҲ�����f����(Ϟz�xε���?w�:C��%=��9Z����  p���d� p����5.�Ea�ߒ�/&z^�^����X�9��g�/���9�K�o�9���z�G��(���<���%�z�<̢\�������٩�f��u���ߖ3fϱ��}�y�T�>���?�(OI��N=_���T�=z�}���b�����u�{��8��ȥ�4WK�S��>�+Ys�s�����yhKV�,��=�r��K�����W%K���{Q�Jj�۳w��Ѣw_D�N�\?w{ֲ������>{���T� ��j��yc ���^ָ��QK��b��e�j^�ֳg�[��ǔ���y=�\��쩧�Q�c�����ػ�k�<�����?K�9z�%��� �3��Vk��������5j}VZ��ȩ� �SR;�K?�O�g=k��A�=z����9�\m!�GI�\��C���ϑH͘���l���(_ɚ{�ڵ����C�X��g���Y��F�X�<��w.���1��y��SR�<���m�wO�>�[;����^ky����<�5?�=�  p���d� p����5.�Ea�ߒ���Ȣ\�R�����}��ѻS��^�:o��K�,�u*=}��ŗԜ'���������}����$���,�q�m�T_�Fz�T�~�sl��C��h}NZ�����g�)�ץ���ֳ����X���g���=-k{��Z+�GI�\���Y�޹�=��-={|�+Ys�S�Ƶ��yhKV�,��=�r��K����1RZÞ{���[��SR�<�z�������?w�z�f�N5�(ϒ(�1{��`ϟ1  \��69�7�4} ��'z	X#�(��]���%��_晊b��k�Z^.G�c=/��w^��-�����K.�yh���J��|�г�<�D�%��I�3�-����F����,����{٩קGQ�V=��r��=�(~�m��i=�>{�A���v\�p�R�z�<O�ا=k۲7O5WkE�(i��S�=k�羬��)�� �W��^�r�k����:���Y����k=c��\��=1U�#=��=�Vԟ�ھF�Kz�n��5�\{��W[�g��Rs����O����Y�k~{�� ��E��ɼ1
 ��E/k\�¨�%�/n{^�����Rԗ�%�p��UO�Ƣ��kZ���k�[�Y�(�1Q|I�\��w�k/=󛵜[���$ʓ��:���K?zƵ�z������{lw��C�u>Z�����g�)�WO��(�m�3-���\�٧=gO��S��ZQ?JZ�gZπ��l��{D�S�<��|%k�u
׼{?�c�Z��(~I��3�(�9}��혞�J{�g�eQ�S��RR�<��o˳���=�����r��g�~�]�(��g<Q���=���d={�� ��D��ɼ1
 ��E/k\�¨�%�_L��H϶���'�u)_ˋߚ���D�%-c�׮�X�K�,�5��u�K4��_�l�gDy���KZ�q�cIϳrL�UZ���!���-�ǒ(Ϡge�z~k����[����l�{l��Ö���V��j}>֞]{�A���e\Q����6�s���zֳg���'�S�z6����j��%������sη�u�j�'.�?�u�Ǣ|%kk_����ֱd-�Bϳ��<z��9���-s���li���%k��T�U�������g�f�kxLΓ�Z����,ʷ���n1�=�H��kSmޱ(OI�=���l8�Zl��  \��69�7F�  \��%`�K~Q��������Ym�%���K��>F9�z�ݲ?��[�=/��<���f�\�=��� ʷ���*�/i�����L��z��Y��cG�X���<�C^���%[췬g>[���c;�~΋s�|�j=������g�)iW�
 �zm�R}���[��Y�w,�/qk��i��L�6��X��/�g�5�%K�[cm�Na�/��j�����K͜�9#�|����\�;�c9�%>{�k��la�5�Ϛy��ǚ�D�jm5��ֽyJ���Ӛ��r&G�5�1�����}Ԩ�ݦ��5c��[�7�]����p�=�;���k�05��yk�(�j�k��k_�5��ɿ�Y���Zs��0֭�����j�z�_K9֌��k��-�U��\�z���i*��%�s����5��wO�޳��f���/{�Ϛ��� p����d���W�*  W�_=�Y�����E� �̹�{�JM���</y�cɟ���/���-��ƿf�5c�ݓ�s�Yx���.�m�\�]�:���Ϣ�V�?���c��4���<���Z;�H�\/٣�9gt�c���jQ�^[�l����ڱԬ�V�T���ܾ՚D�kl��j����8�k[�3����Y��V��(�G�{x��{�f˹���"��E��C�����Y�7r�-��-����[��<l1��f�׎�ؾZ�����c�\-zk\���k��۫��T��1[�?�������%���ط�Vs<�ӷ-��<�|�����G���}���k�U?��m��r�M�����O���3tl~jǱ�X~  �KT����`  .�/{-�(�Jt�-��?�,���k��-�G�ܗ(����Kߞ�T;�[���ڹ���u>�s��g��}ﳳ��\�l5ޥg{�s�q��ploF�]+�;��)�>���i,�[k<g�}/�]����O�ܞR�ǭM����7��-�죭Ɯ�י�cϵ�yf�|v��kM�(�V����k]�{w�u�#��-�#��Z�WQ̒5gE��r.��3��"�Y+�ȿ��6�ߢwo
          (�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7�īS      ��B�  <yIDAT     @��69�7F�           %Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�Q0          @IT���2}          �ET����`          ���69�7F�           %Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'�Ɨ��          4�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o��          J���d�����5          �&Qmr2o��          J���d����a          �Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�/��           � �MN�Q0          @IT����`          ���69�7F�           %Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7�"}          �ET����`          ���69�7F�           %Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o|�ϧ           �Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�?�>          �"�MN�Q0          @IT����`          ���69�7F�           %Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�?��_          h�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'��W�          ��j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o|�kS           @��69�7F�           %Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o��          J���d��3��           -���d�          �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o��          J���d�          �D��ɼ�g^�           D��ɼ1
          (�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o|U�0          @��69�7��u�          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o��          J���d�          �D��ɼ1
          (�j��y��7           �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o��          J���d�����           -���d�          �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o��ק           �Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o���a          �Qmr2o��          J���d�          �D��ɼ1
          (�j��y�Ͻ��          4�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o��          J���d�          �D��ɼ��~1           4�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7�|�0          @��69�7F�           %Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'��(          �$�MN�Q0          @IT����`          ���69�7F�           %Qmr2o��          J���d�          �D��ɼ1
          (�j��yc          P�&'�?���W����           �\�<�KNr���>j���~C�           �k��u�I�U~����{^��׆I           "�yZ���Z�^;j��/���פ           ��yZ���Z�^5j��y/|��kސ          *��i]r�k�x�����ȋ�$           �\�<�KNr��?>j����y~�           �k��u�I�U~�GG
;          �N��xZW��Z��׃���$�3�v          �]r��8�5ƹ�x�뻒�|��$�          p����i=q�k�7�>;9�����y}�9          p����i=q�k�7�~2��O���          �C��j��\c������&���o�y��RG          �[)�O눓\[���5���>�>9�          p;��iq�k�7��Trp���{~o�)          �v�5��:�$�o~=4���Y?��a�          ��k���$��r����fO|ғ��          �-�
O뇓\S������
_���䠃��~�          �,���8��/�_���ɬ��x��          .C���j��\#|Q�s��N�����          �ːk~�u�I�
          ��\s�ko���I��uݿ�?�M�_���'          �M��6��Mr��kt}m2���}�s         ��.��F5�I��uM�wH~=9��?���'          �\k;��MrMn��u��%�I�_����	�O          �}��6��MrM��p�(9��G>���I         ��*��N�n�\��Z��L2��/zғÉ         ��&��F5�I��uU\ߑL�o�-�������	         ��"�����i�m�kp]���$�I���'nޘ&          �\S��&���p��d6�O��'          n�\K��&����x=(yEr0����������υ           �U��͵����$����[W����٤~��?�捿�&          �\C��&��ֵ�zj2��o��o          n�\;��&�$q���>yCr0�����r��7�F�           p[��\;;��Mr�m��ump}j2����K�.
          ��f6��Mr��k��_%���?���Å         �k�ke��$�ֺ6��-����~��~�pq          ���Z�i�l�kjsm�k��s�٤�����         ���5�Q�l�kj];^�Nf������B         ��ɵ�Q�l�ki];_�Mf������n^���.          \�\�kc���$�ҺNp}Q2[����O
         ���5��V4�!Mrmi�1u�6�>0��d���F|�/�����         �ȵ��b�\S�kK]�ͯO��d����p�S�zM�a         ��&4׆F5�I�%�5�.�n�G'��y�{���K~�Í         ��˵��~�{����}T�r�~��$ڀ7�|�n~��^r�+i�         p=r
����.�Ο�B��������         V�lزe˦�j�z �l鲩33��M÷P�eZ�
_
��m��P}h�۰a��&�l�v�e�v�!_o'O��.��v����3s�+?�         �#[�l����+[�l����F�j�>�lɲ)˶,�l����zm�	G�υ�c�kF��F��v�}l�;0�{�u���ƛni�g�ߞ������ʏ0        �@�mT6R3�{��Le;�



�� �   } !1AQa"q2���#B��R��$3br�	
%&'()*456789:CDEFGHIJSTUVWXYZcdefghijstuvwxyz���������������������������������������������������������������������������        	
�� �  w !1AQaq"2�B����	#3R�br�
$4�%�&'()*56789:CDEFGHIJSTUVWXYZcdefghijstuvwxyz��������������������������������������������������������������������������   ? �S��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ������M�����+[XT��L�u%� }hz+�o�_�S/����s��L5X�>��t�0�fȄs��}��o���ZOx���>�L�#����56}B��
��Eob}h�
����)��G��-�����?�iq�Fn}�
�!���K_���%��S����Z���,v����8�!��� �է���_����V���Zr}�����k�Z�>�ʵ8�϶�u� ���Y%�Gox�E�d���&"}��`=�/�pod�#�»x?�s�j�&~�G���5���'g��
��N|@f_
!��x��6;�r�3�ҸgF�?�6:#R2�Y��E~q�� ��x�FO��Լr�5������2��(�
� Z�{�o�χ_�ρ�g���m���ܩ�1��''�	Ebhw�QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE Q^e�����q�:��^<�E��][���o.�� ��r�21�A�E zmy����~�ڞ<�^���KE
�W��s7�w�]{F�2E��;w+ͼ��N�ɭ?��*���I{��g����j���&�bG�kjtjVv���+ɞ� ���TŮ4߄~���x���ܖѷ�?�S�|ըx'����n�� �ږ��,��D����N�'��(����2k�?�����o�����T���mKd��a�� ���5�R�z�q�<_H#�� ���L�Uǌ����Z��}���H&F��_�}G�_���ej �t�cK[u�����k���Zt)R�"qʤ��0��+��(�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� *+�X/�䷹�;�$^)T2�� �jZ(�>%~�� 	>$y�?�Wú��?l�XZ���<�>��z�gǟ�M��"�u��&�Y�ݼ�#y
3�I~���:�;+��G���[��(�LWc������߳� �3��ƹ-��^�O�~!��~ u[i���|��?�OaSx��:G���t�sL����~���4m� `Fk䯌_�M���Ͼ�-���[-�)�\Y9�;��а����l�7s���$�~��<wP�42,�ȡ�H�2���AA%~x��;� � ��`��� �R9x��'P����rX��ҿ@f��+7��� d��t�᷉��y����'n�-���N%
M�kɔei+3�II]tQQZ�C}mż��[̂H扃#���9ȩjFQE QE QE QE QE QE QE QE QE QE QE QE QE QE QE ��xwᯆo<C�j�@�,�t׷��_A�Տ@�$� &�f��� �~ ����B�)�/�W��`#�ld�vu ˟Eu~b�A���
3�5ּC�I� ���%V�I�-��|pq���u��%Q��]�)(���ڃ�	�I?�>�[�[�� ��R��<��>�lG�sѤ���y�	�c��G��I��:桤%���>��}V�>����y#�������^�m�u���Dk����'��t���I�i�v�\��[_#Ω�oH��/��	��f��̻f��y�S� �+|�g�Fv��+ٌTU���ܝ�QEB
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��"�����[k�c���JI�OPA���7�w���^=�F��G����}�4/a+z���  �t���V5)B����#9A�,�������� �wx��9n�&�i�k�*�X�����������,� �B>~�K�gv|-�6\��uiI!��yxYǰ�㒀W=�/
uwlp�$ H ���{M�y�k��^�gMsyy*�1����@ z���������5+� |V�I�V�E�o���'m#����|y��N�-������� �����4���0��A�|C�q}7Fa�g�^���o��f���� �݊_�ߌ&�mƱ2qG̐)��؟���8�|,�Ms
��My��6� �=e����7�<�.#�� �Y��ݺ��ROT���k�
��T�*��
(����(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(�&���2�/���LZ������]j�B��q�c� -���0y��t�C��k�"�ZE�\�n�a�`<�N���L��ܘ��q��#9�C�� xM�V�y�k0jZe�f)�nP<r)�A�  ��#
��$tү*zn���C����Zx}N�t/�Ǻ� �7S4x�d��<��@?�.F~������4�W�C��
�~"�
�aZ<�F��*n�>S���w:-Տ¿����"qg���P$I�[ޖ� ��S�<?0�XGYYX2��e9z��/���|d���_������&�b�"���,�߀�`�1��D�O��շ�/��\E��&�
��S�Ǎ�M
9�������YZ���,�f��.>!� �K?h��o[�M;�����k}"�q�Y�ny�Y�1 ��D�m�]?Þ��N�lc��:�Vc՘��ǒI5�`�������:������_J�_��t=�7J��Cok�QG�$�I�I�5�E��M�QLAEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEP_/~�_���
�h���1QJ1��m݅QT ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ������Wo�M�k��Ih�0��-l��ȼ�o+�c~:?��~�Tw�^[�oqO�RH�P��F
�x ���^�kÒF��:r�FO���ڷ��>/��Iv��+���/��� !V��A ��p��'��᷊?a���/���Oa��v�钮Y-���Q��6R�����I�ٿٯ���� �1��G�Ǉ��[��o�X���6�?Ĭ�޾B�9R���G�)�dz�QY�QE QE QE QE �G�o��_��	�ЮTx� ĩ%�����(�ėd����u
��Ϗ<q�|5�n���%���Z�ywp���v=�I M~Aq�/�(��u�x�ZY��8�$� �\lDV�Í�2:������.��c�&RQWg�� �=f�n�� ��⸞y�g�ņ�fbH{��RNB���������n�m�i��6P%���KF0�Ơ*�� ¬W�P�PG�R��.fQEtQ@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@��_�z?�� ���"�e}ՕT�A�J���pG�x&����8x��	��Nj��|�ΖZ�k�(�O�5���rP���1�O_-�ޟ���-χ�����g��;H��溴夋�yu��ޯ/���������+ٟ���P�Z�sm4w� �9�`��FC)A�EK_��H��������^n��;���+ nt�a�	<��Gd&�G+�\(�� (�� (�� (�������~ο|Q��Kd�M�?c�s��]?�]s�r��!C���� ���Zմ߁~���W� �9�d���Z`rJ�HGr�wS^��%|��/�[:xTx�Q{��9>qD�c\/�w⯏� a��Z�ǟ���ş���i׭z�?j��%��<��bc�_�5�u/m.�f*��"
(����
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��?1i��?���1��K�:������WO��
�9�I�;X�Y�����QĐ�?�2�u���� h��V�������hӮ\��@|���?�W���H��"���]o���e{K=fi$ӡ�l}�S�,<�1� ��P9c_+���j]l�gS�6{��
�(�8�
(��
(��
�t� ��|~����B�+�v���@�3wo�ڵY�T���h�G�SHJ�J���M��o�~<�6�5�M.q��p��=�FE� �W�g��=K�'�_|U�3�����x�3s��J����زڍ7Z���EI(Eɟs|�Oi�O�N��;m�5�;�.���W�n$E
;W�QE}�b��V��m݅QT ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( �����Ꭳ�����k��&���u�uo����!�A�\ a��nOZ�"�<��>����|) Qsu�e+�;���l��ؑ޹1T}�'�
(����(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(���a�����5�	���7�`��E~s~ž$��N� ��'�/fh4�J�o��6�;+Z9ϫ�s�ٍ~�W����sxS◄�y����R��x�)sl�Q�ٶ:� \���i�RS�v�ei��?sh��oĈ~/��_�!+�sI���W�r�`ȟU}�� ��h���muC�T� �IuM}�ڴ��彍�����m���a ��=~�~ҟ��T~�� �Z�yW^�u5�g�2�����?��� �Z�7����W$|[[C�B�u21�@>�T��uaa�+F&5��M��:�(��< ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��?ࠞ� ���s�/b�}�s�;1�>�$f?�
�N��e�<e�
��'o�W��m��M�k��7�H����1G����� �񧙪|+�o�*�V�<���E�6M�׺|���!��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��?6�3x~� �R�{���ַ*�3tQ��d� ��t�M����_�_�R-2��?��<Ub|��t�fI1� /�??����#¾ ��o��}rӛ]N��y��"_х|f"<��3ޥ.h&~0� �U��۫H�ڱ��^���Y^\}�~�*�j��T`( zW�/�K��#�
������]�x=1ek����Ѻ�2����Ìz�QE{g�QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE |/� N������k۫B�� ����A?�~����%� ���@�K|_�1��E��m��?�*�� ��h�T���շ'M�-ns�I�կ�?��~":��_��Bۿ�5=B�{fs>?�=|�aWo��c
�L���M7�%��S��~/�I�Av?�� J	�n_ʿJ��+�s�~ٞ#���ũ��� i�PO�>k�ֽL�~������Q^��QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE x7�ӧ� i~�8P2�%���e�,@j��C�ŧ�f�i��K�Us"n=�Zq����_�WZ�����6��'�p�E~{~�/<)�;�Ky��}NIH�1D?�Z���~�/��p�GM� �o�~� ��� K�s������u_��L�/�� ��7�Ui_����� ������Q^��QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE y��%��>|H_O_7���_�>	�'�>�,;��:��W�7�	� $�W��:����_����k�s?�>����Y����K�� b����ZW��~a� �0�k� �,��U�~�W~_��sb��QEzg QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE�~П�@�%س�� �,���_��'��?�_�,��K%~"W�f}S	����&����ś��*������� �`� �|�� �Y�� ҫJ�<������(����@��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��(���?��J� �gS� �Y+��n� hO� ��Y�� ��J�D����8���g���L�/�� ��7�Ui_����� �� ���� سq� �V��}]�{�����e�+ҹ�2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��y��	� $�W��:����_���{�B���� bΧ� ��W�|�g���0�>�� �_� �}�� �Y�� ҫJ�>���%� ����ś��*������� ��9�_�
(��#�(�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� �� ��HĿ�u?�%�������� ��/�ŝO� Id��*�3/�'���Y����K��� b����ZW�~_� �/�k� �,��U�~�Wn_��͊� QE��EPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEP��п�@~%� س�� �,��_��/����,��K%~ׁ�|q=L'�ϰ?����_�� �n?��ҿP+��	y� %�_� �b�� J�+��p���د�
�g�����I_����_OO	��/�%�����ŋ��+���:��� �]�~�� �X�� һJ�C���NlW��(�@�
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��<���� ��3�ōO� I%�ú�������L� �cS� �Ik������>}�� �� �� ��.?��ҿP����	s� '��.?��ҿQk����A�S�@��)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO����)�PE>��(��@\e�(���Eq�S�.2�}�QO���矴7�����,j�I-~W�?�� &� �3�ōO� I%�Ê�s/�'���Y����N��,\�]�~�������N��,\�]�~��5ۀ�	�'��EMMz'(QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4 QF�F�@m4m4�ߴ?��� ���5?�$��+����� �o�����K_���?OK��O�%����?�X�� һJ�F�˟�%����?�X�� һJ�F���N|O��(�B�(QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(���ߴG��� ���5?�$��
�������~'د��$��_^a�D���>�� �Y��x���{��+��Ժ��� �X� ��x���{��+���ݢ�0?�0�m������sXm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�Xm��������N�(�(�X��+�M��w������K_�����E� ��;��}S� I%����3�'���Y�G���N��+��]�~������N��+��]�~��^�&��Q]�0QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE矴W�������O�$��k�K��� �}��� b��� ����mx������,�'�	_� '���?��ҿS����	_� '	���?��ҿS��G�����9���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� ����_�!��h�ߎ�
{�Ok���H����!����c��b x�� u����e)�{z�ʄ��}�EU�u{{M��t��}F����Zʲ�*��VRA��i�%�PQKE %�PQKE %�PQKE %�PQKE %�PQKE %�PQKE %�PQKE %�PQKE %�PQKE %�PQKE %�PQKE %�PW⯏_|�O�x�������g{z�J��2���A�k����h��'��������K�t=R�(�0��'�Lq�jq�[���{�U��We�E�y��� 
;(�%c�?~ӟm%]蚴wX�����ۚ��AY͚���)�W�p���E 6�u �)�Ph�Q@
��� &5����)�WҞA����υ<iG����0���S� �Jk�?�� �M� ��<��_�q�}X�Rm=̶Żo��Odd�_XQYΝ:��q���?�>~�~3�����%�˧\�-^Ж���=U���9�*���
��uȄ�.�o"�jk��0J�g��z4QEy'hW���� &����u� ����%~�~���17����� J�K�G�rb~{��+�<ѴS���N��(kzՏ��{�[S���N���nne8H�E,�}������Z���`xc�R]�[iS��O����HHmV�!Pq�	f� iڽ��
S�M	�7���U+6�q�OX�r=8w� ��
���O�����1i��Sxت��8�B�>X�>����+�<�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S��:��?���'� د��$��U_�� �g����� ��O�$��
+��>(���g�_�J��8O� دq� �v���_�� �J��8O� دq� �v�����V�F��J)v�6��9���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���o���N������ ��~�m5���}��x���?�"��;�5�u�7�x+�]Gc
�mb*Iƫh�zo�.���+�i��-�hg��B�E"�O*��<�*
�/� ������]�S�}��Nm#ݮi��}£�>P�h�|���g���x��:2�g�Nj��
(��4
�~�V�~
�@Ҽ[�����>h\�*�#�� U���0@#���n.�M_F~���š|v�w���@������-'\o�� �RG�a�
��[�G���ٿ�$wr���%�2ìi񜒟�4c���$xn\��~�x^��U��kE�z��}
�[]Br���*��5�Xz�|��U��"�Q]F_�U?���=� cM���]��u��S�M����4�� �%�s�?�#Z_>8� �{��x#��� �!��d���	�� 's��o� ��d���7��F؟�(���N@��( ��( ��a���oO�	���tȏ�<]*�V����oI���('�o���[�����jY�����c��3B&�?�J��?�c\uqT�nΈQ������?�x%�<C�Cu)��0�ߓ�5�]����|G������?��k�IiI�ؖfc�I�I���� �_�x̞𦷯Ƨi}/N���S\^�~�N��En��}+���A�H�m�'�;��V}J(I>�9ף�zņ�j�Zu����}٭eY�I��O����$> �5M
W���Vr[��QQxgŚ߃5H�-W��o�9[�>�����Mc���<2���"����/� 4��'���6�x�HV�
��яP�?ї'xW�g���� tW�����R��Ӯ�wm��dy���ʞ��},E:���)���z(����(��(��(��(��(��(��?�?�h�7� z�� K������� ��� ɣ�����.��k��z��~�EW�u�}�� �� ����G� ��|/_t�'� �x���� ��]X_�DƷ���kEWў@QE W�~������A/���R?��W��{m|�� ��C���W��2�YE�M�n���G��;S��.�Ə�*(��\���(�v��e�ڻ��C�����6҆�	��z��� �n�RO~�:M��|:%�����`_�zq�W��{����y���Q^��QE QT��f�ú]֧�^���u�g���GJ:�18�@*;��l���T���F
�=I=|�@�O���.�_�z|w�ψ�$>Q�xa��g�O�G5��O��8��|�^.�F��;��&"� r%�'�EpT�BGS�yKW��C�~��	�7!�P��ᨥcMN)}U�X���7�[ɼ����+g�V�� W���;��d>𦩮���Z۟%���A�'5s�����]��|W��WF��#7�ú Ǡ2.T۞k��ṃ����1���xwƶ�������ێL�m�w
>���=�V����}���iװ������H��2�G�_^|� ��x��w6�o�Cx�@�V�l.����K�~O��kO;MX�a������Es
U5��ׇ>���9�Oqp�p~�c?�����c�_��q�jϩ��[��oܜ�j3����
���ah�u�)k-�}K���5��Y�$�u�� �k՜~q�4oڻ��2Ei�#âG�V��m���ͼ���0�q2Em,�B� %������K�z�����^أr���>�
���� �����@Zn�e�Y�w��[��I�Om*�}I�W�?��%x�ᦦ����A�h7JۋY\2+�:�k�fWٿ
� �Z}�����"�5-*E��� �d$���?�k��6�z�%��~����o�4K�cX��M��i�.�,q �I5��� ⧄�*��^xG^���md�-��
x�����?nMS���-�7���o�q�&�	'M} �S& T�g�1< >�� �L^y���i��W�r�����%\qJ�U�L��Ù�}�EWq�W�_���?h/��J��R���di*n� �,c��K.�T�>vy}+���c��u-�ߎ�+[��l��en��ЃA�*��ӗ,�:#FSWL����?Ṿ:�B�� �[o�5G�77�O�(W��m� ƫ?�S��������W���77�O�(W��m� ƨ� ����� E
�� m���^�هէ�������o�F��o���u���կ��g��eU29�
 ��zEz|�3��;�����Þ;�.~���J��˸�(ϝ�_A�־,���>>i�R�]��R������V��}1d㩊�9r�3xє��?hh�ŏ�no���P���� �Q� 
  p �?�T��nuҢ嬶5�,|[�7Ư]���W�{7ʑ�V+x�%b��zu<�I$�E����E+h����� �\>�ŧ�Lu�U�X���a��ӛ]B��;�i�c�20p@<�Qf�=��F��Y� ����� E�� -��������%�� �o���3�������W���7G�_�(7��	m� ƫ���4�ֽ�+�m���<K�>������r"�I��xׅ p���ڶ���ir�ԣ*j���(����(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(�;�����(ح��$��M_�_�o������O�$��	k���Q=/�ϳ?��?�p�!� �Z�� J�+�R�*� ��?�p�!� �Z�� J�+�R��_�1�?|(���9�QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.���~��x� ��?�!�����a��?��<������ם���:��>�� �T� ɼ���ƛ��$���+�?�%O��ψ�i�� �KJ�2�p� �Wﰯʯ���C?	�ɼ�
��Y�!������ F�,��4��Z���a��6*
u����b�Q�w?�
+�k�ُQ��~!=�^u��"�h������C!y��	�!��2TxU|ܢ�'nzђ��
(�����?��?���+Mr��ǋ����S��ˆ�4ۗ9�I��HǞ�������]�:��%(�(��3� ��"� �z����m��ŷa�O�@F�{3|���Lm��ƣ�VE��Y�۵���F�T�y3���
���
�� &����m� ���̯�� ���o>� ���� I.�<G�U'��o�'���ׁ�ݿ� �����Ɵ�'���ׁ~�� �Aq_��́��� #\O�QEz'-(�.���u���=������ �u`�5�g��n�o�:3�AʌH�_������ �GЮ|���+I��p$��nv��D��"� ]�Ggv,�rY�I>����?���N����3HŘ�f9,NI4�S���&H�F�G`��2X�� ��>���� c};��r�� ���j�Ŵ�;Lc���O��c�j~]����
C~����ivp�Y��ii�	j8
�8 z
�>x��
�+�8#TV�
(���QE(�.QEs�?�(W��7�?ް� ��z�k��O�(W��7�?ް� ��z�k��U�zXo���Q^q��O�w�Jg�?����}�� �� ����E� ��ua���g�Q_Fy7
(���W��o�.��{ipUw�o�V�����c����_g�r��῀^��7�� U;K�y׳���1�{��2OJ�I���#Y����gź��v���e`�v%	�*�UÞk��UQ��[��'��9j(��3�
(�j��f�_�����Im�I�}KlQg�Z��>Qۖ親1sj1�Rj*��g�	w�v_
�9�|y���W~#�`��ai ��!o��Oz�v��z��t['K���M��-��aH�E
�� U��jT�8(�ǜ���QE�
(�f
���I4�x�C�i�=Gľ#��M��#�&���
��1$ ��H� ?j���ߴv�-�<�?�m��e��`ɏ�,�}����,zO۳������|7�^���$�m����NָoU2��9?�k���X�7��J�.U�-½��G������-�h���xv�>۫�Fp�  ��wbz
#��m�2k��2���ê۰�G/>�\�q�ҿ+�qT�:�[3գ>xjz����F���s�z��+]is�]GG�Ȇ�?���;_��h>|Iо.xK�W�.�ޗ�D�th�vu9z�Nk�"��� �t���/����K���D���m|Fؤ��#>�C�5����%�"�>e̷?Y����O6�_��R_ښH^o�^�h�M~�� ��j�A� U_�
���w�-��:�'����:>�5��c"��f��5�3�j0���공֧�\�ws3uy�1��y�ʮ�]N�<9�3�fV� �� �_<m�x_G�I�j�Ikﺥ�7���`k���,�v���{�z��$PX��KV�c����HQ��H���Z񠔤���i6�կ��7��g����:l3�X[�r� n�_1��D�DS�=NI�K�ho-��$�	׎E�==EKE}Db������ۻ>f���|5��euu�i��/Ĭ����#�߷� �0�+���OJ�����;� �a'��Uc�JA���Z��<���#�x W�
���o�;�M������^�_UO�^��'�0��*ɹSV� �]��q����_ж��  ����� �	�禼�����v�QE��_������o������ ҩ��j�i�a?�4����q� �SW���#�91?
=������Q�_��-Ҽn��,����-���.q!e���"������G��gs�����>
���'X�� ����V��>'|��:���:�N�5�٤�m>eef�9�F��q���*o����� c%����וWN��z��JRI��(� �
��������i���Q�nc����]���rH,@�8��|��á� S�� �Q�E]�'����?��$�� ���=G�;����	1� ���� �W���F�w�|�;�3�~6ü�=� Г�,��~�~ǿ5� �?����-�{!�k�k�E��%	��i�F*r����ֽ��ږe�g:Ҩ�(���Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEs�� h��7���+j��I-~����� &��C��mS� I%��:�� OC
&U~6QEtG�o�:�/����|E�K��
(��.h����5�-WK���Q��.-��b��AVS؂�ٟ�� �����~��;O�j���+�f������?��/ oŪ��
�� &����m� ��CR�&�
zTI�� �=� ���	��� �����Ə�'������� �Aq_����o���1
(���*���Z�]�}:������'ݎ4R��� O�V���(���>��w�m��]� ��SJ\�BA�c�*��t��5N.O�Q�3H����>1^�v���x���\K��[�� Qj����O�3���Q_/)96��%ed�|Ї��8�?ҙwGu�X� � c�M���k��p����~�_�#;o�����T�旘��[?l(�����
��Ɔ�(�W�4g]������D��c���(�@?��+C�V�� ��q�"�X����+>�H�¾�� �w|p���_^\l�� ���L��*\��y=K/�'����[;��.��������H�C�F �؂+Jst䤺(�&��
��'�_���z|�w�,�e�)�G��H>�Bdt2�ھ����)ħ��{��7��xf�=>%S����Y�����+�S��vݛя4ϕ袊���%�����x�4�#DY���@���>𦋡ۅi�P�F`m�5A��k��}��>4�Ȯ�s�
�� @�����W�O�Z� �^�_UO�^��/��QVIWU� �]��q����_ж��  ����� �	�禼�����v�QE䝡_�� ���������� ҩ��b�j?a_�4� ����?��M^������|(��(����
���
�� &���m� ��5�%|{� J�:g��v��"���/��4QE|��w � �|9� ��M� Ҩ넮����O�?�2i��U\>$)l��Z(����B�(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(�;�����~(ح��$��G_���w������O�$������Q;��3���%G��7��V�� һJ�U�ʯ�%G��7��V�� һJ�U����1QE�s�Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@~0~���w>?� ~�� Hm��~�?o��;�� �e� �6��c���N�?�϶?��������k��&������%����{���� I���Z��� 
&U~6QEtQ@~�߲	������ ��/m _�T�A� �� ���F�HA_��W�o�;�A������4[�k�/�1� ���nO�?��=�b�� ��?3��O�σ袊�Nࢊ(�c�ڊ� �m���\�.���ëخIQ�n#��3���r�v���� ?�,�A�1�GT�n��ӯ|GiqmuIc{+�WS�AZ������N��>������Y�Y�y�a�$J?��>�q�y�q�
�4�M�c*i�M�� �� ������ ��fk��	�� '}�/�� �}�~�W����� #�񠢊+�9B�2��^,7�<ᥐ����d��˰g�� ��~������V}K�������}����/��1ڗ�Ӈ^��
� &����?��޿k�[�
� &����?��޿k��� z��aEW�u}�� �������qs��1`J�����������$�� ���� P��utΉ�?I&�ޓq�K��K�Q�{t*�F �A���j~̳��� ݴ�d�,�\�!9��1��dc=T������y� �o�:7ǯ����u�ؗ͵�Q���P|�W�IwVa޽��UX�npR��/#�v��!�Z�[�M[��!�6z��1�h���*�{�)pA�v�}����QHf��h�� ��╮��[�^���lP�ͷ-��t�#9e��K���v��Xx�E�������������7P���E=u���g����ԣ�K�{�i����ĭ�3�,�ğᐝ����9.1��+r�g-��^�2�G�MQ^��Q@Q@
(�� +���� $������ �/]=s�������%�e�����+�h(�� ��������� ��� A�Z�%������� ��� @�Z����/CƗ�¹o���K�c� `k�����-�[�Io��
�� E��z��<i|L(���J�����?�� �&��z��5_���� �M<�������u
(����+���� �O�{� ^�� �L���_����������� Je�G�G�rb>{�Q^��|{� J� �p��b�� ��|� H� �o�� �b�� �7ψ��i|h����+�X+��� %������o��G\%w � ��9� ��M� Ҩ��!Kf~��E�'�QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE y��� &��C��mS� I%���������~(ح��$��C^>;�߆ٟf� �)� ��C� b�����W����� �S��x���k��+���z��� �񉶍��Wu�a6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6������ 'u�� ���޿h��{��� ���� ��_�Co^v;�k����l�w�	f��������������� �K_�6�C��+��o_`�N��&U~6&�6��]2mih��&ڂ�O��,�,�-⻳�����t���VS�<j�\���W�O����h5&9'�&�3:rKI9ck!<��'�(=J�>l���_
��+�[1{���z<m�dF�S���ҿhO��������´��2�ڈM��ےv�=f\���`�Cپh��J�Nug��tQEp�!EP�?�O�;� � �C� M�5�7��?��� �w�� �����k�r���߯�v#�Bm�m-���m~2� �Bch� k���0�+������7_�� �S-�?�
�� �V�S��?��%G�?o���/6�p��5����Χ���M�m�����&�6��E��~�"�������R�� �V�z��_?ڼY���̾�� 9�M|�����ER�� ����✙;4�#ܛ�#_�k���$��m���cf�S��
� &����?��޿���z��aEW�u}�� �\�B��� �\�8������_�J����5Ն�4Lk|�3�F�Z+�yG�� �P��g���O��9h��E����,��ut嗹��W�U���� �;�Z?	�h|w��2���cq)���%0:$���2�6���Q� �����*}�|uEW�v�Ioq-��S�+�<L9cb��C9��-��ӑ~��
�mW���K�a��J����>��z7��Ǖ���#��0*�� ��_��� �NO5��������qᬄ��=I�'��� p���^�ȯ��T�㺝~�?4(��慨�gX�ҵ{�7S���qiu�X�uVS�5F���(�� +���	��^]����]�;�?U���P����4kw'����{��yQ��
��u��L�t����4넺��)"0`~�;��j�2�Fu ����Gmk��S��O��
(��?rd�� �����v��
����o�G�/������V�$����J�?�� �&��Z���o���� �M<��������u
(����+����� ?��o�)������� �O�y� ^��L�����LGw�F�Z+۹牶�?� ��/�cm���6��&z��@� ��ɶ���k� ���|C�ԍi|h����+��X+��� %������o��G\%w� ��|8� ��M� Ҩ��!Kf~�m�m-�W<Q6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.����� ��G��mS� I%��
����� ��~(� ح��$��?^6;�߆ٟg�)y���E� b������������ �R� ������k��+��՚��� ���F�u�s
(��p����� �l����zgu���o���]�~��u��/���3�O=ţ_6�T��e��j:T��>~�m�m:��S���<�o���F���Ե�|DԆ��� ���v�].�r}6���)=����f���W��1s�'4�(��=������u��4�<1��e�uk�������I=�4���/�'o�_�� �߇�|�����U���F�C����k3¾��w�4�N�����H����8�"��
կ��y"���.f�ݴm�QZ7mi�PvѶ�E 7mi�PvѶ�E |�� _��?�a� �����_�_�P��4?� �a� �����^&7����p� 
(���+�/�$�?�{� `�?�i��k�/�$���/� �.�k�
{���1�
~<C�jW�D�LK���p�q�6��E���5�JΝX���!O�
(��?r� d�� �f�i� `;��omx��w8����09H�?�W�W���������cv�+�a���������!묮C��V	<m4�0��zZI*���OJ�l�[���QE|��Q@�߲j� �4|4� ����+�6ה~ɿ�l� � �k� �
���z<i|Ln�6Ө�	)��� ����� �	�玿��[�Aw���� �_�y8��;��B�(�(�
����W�1?���zM� �2��~��¿�i�?��o�)����Ɉ�Q��h�N������A� ����c]�����*z��I� ��C�~�l�� U�Z?�"� Z¿�kO�G�mQ_6z�]��/�.�?�d���:����?�_��_��� �TuQ����?{6Ѷ�E}I��F�u ݴm�Q@
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
�@��5�u� ړ�=��#�
�0�`x �ⱫMU�+4�����Q_[~�_�#|
���	3r��S�����h� �:� �
�^��� �g���¿�ͽ���״��62�
d�L�>�� z�z�ûU��]`�ת(����B�(��8?��ɳ����,ɆfX�Akp�¥�f�_�kwc_��SV�l��.�M�m�������-�]�,n�YXw?ƭ5V<��p����Q^��^��ן����I|-���o�-F�� ~2B�PU���󒋃qg����
(��f��k��9�v��m��j�]��[�@a؎�FA$w�����㶍�B�2ӼQ��Et@�Q�
[�z�����D�+뒽5ن��g�0�O�i��yEfxg��_��?��%�:��
�mue$B8#��AEi׼y�_(�F���O�-�
��3�_W�d�&���U�`U��A�"��H���FN-4;4W�� �'�7����-&|�J�i�v�I��jޛz����j���nqt��#֌���QEY�ÿ�i���|/⦚�����h�0]�َ�W����Q��9
��9���\�zW�@��r=����+ns{	�?E�#|H�������:�ZV�h>id9gc�G.�(����1�Hx�����g���h���MA��l��p�9ۜp8Q�'��(�f�ƍhj�2�ֵp�G! ��%} �z���	��$˧��ſZl�H��#��̪F
� �k�袾�1�J+��wp��*�~|X���O◌4Yci��ݮ1��3��8�R��� ���#��?��$2��+�oc�/ȷ(N���Y�u��k�*E�n,�!.h�QEfY�1� ��D>%���'H�K��O�P~�,�����3��������Z�����n���|%��k��,�Ƞ��"x݃���9��_�)?��Ӗ��z��FM�Ҧ�G㍟��{�1pJN�mJR�m#�J��� ��~���k���3Px~�Lk7�핽�ND
ìhFI��e����?च������ �n|3�܆��U�p/�c<P��T��%����?�>#�����v��-6�VԯeȆ��g�8��aT^��(f^lF#�~*\���.�(�,�
(��?r�dY�����ՇA�[����z�|������ ���t���Yʽգ�� �;OЊ�����"��zI�QZ�S�G��1�H��_�~�|b�~�M�=nG����܆?�XX��'��z�1�X��
(����+����M�%�=e�-��S���N�`� ���$�\���k$�<�=�唋�KJg � �g�z'���r�>�Z(����
���
9����'��e��uepq����� ڕ��r>��T�i�_���2ډ��r�$�Up�� ��G�=ˋ�g����
�~�6���fլ�]6��n`q�]N=Gp{�
�����B�EA(�%��+��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��<���� �x��� b��� �����~����o�V�?��Z���D��lϳ� ��|��"� �Z�� J���km~S�(� ��<E� b��������u`� �c����F�Z+��mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� �/~�_����C�����^ɂ\Y�3H�P��{��� �O�xj�q��ռ_z�D����?�<�e��[�����>0|Bּc��~)���YVYⱻ�HT�TC[��X�%� ��I� C�?�:�� �k�Ov���e���?���g�?�ST��a�K����-4+6&�y�Ǚ�Ǡ' d����k�G���/�u�� "�� ��I� C�?�:�� �k�XZ�w��J�N*�����?g�)f���K����
om>�����.�J�j`䵆�t1	�G�_���%�ٷ����V7Ac�t��|������2v��d���G��)��iqO����/v�2��O�p�ׁ\0�B@>¾?�g���$��_�~%�}�g}�u�P� ?Z����K㽫� /O�s���SQNX�:(�T�*��}��k�
��isE�;+�k��k���3���#c��Fq��Z�����]���
�:*�   �dӿ��� �)G��V<�q��F܃.!^���%Očbdo��B����%���� U?��)�z������>(��?foأƿ�=����j�3kW��/qn�i� k��<��
��@� ;
ޖ
��:}��RK�:��!�l֖۱�w�>�j��:��'�^4� ��O�E��G���/�u�� "�x��H�J����O���O������7)%� ����J�i$��RI$��I$�k��կ�u�O��i� �֟��G�:��'�^4� ��O�E���'vt��Ҳ?)h�կ�u�O��i� �֟��G�:��'�^4� ��O�E���Q�b�-Y�5+�R��,.e����g��*�H�2��ЂڿT��	?�b��Z�-��?���1x�� �?���ꇷ���
��(DQ�U �
}�R�Q�N	��݉����V�bm��r�լ�'�o�������+��������
���c��2��_Zȧ�+)��G���ѯ�9h�Tλ�����E��I��&�����
�K�o�W��V��v�lx�1΄����� ~� W�|h�࿏�:7����Tɶ���Z1�(��Gl��8��u[U��J����~Q_q|R� �U��C������JZ߷��Ǣ����������~;Y�d�y|�{�um"� �K)�y�R.�'z��<��{�� �N�^"�ׇ�|?yڦ� �+������������P���M�˄!��,U�l��vϙ �l��j��beZ�|I�6��~1���A��07����u�!�?礞�=�Ts_���/�+�\�ǅ�-g�[j
Ί�V er@ ���O�a��N�o�xS��k��.�5����[�P =	�����2x�� -?��߂�?�~|B�񎁬���R��XRJ���"D(�	���vэzv���j��+���������O¿�:�ޕ��iڵ�+˧�2�8���D=A�Y��~5~��;|D�<e��>&��/R(�M��H@�5�pn�3�<�eNt�w��yo�!~؟�P(�3xn�������3t���N�主U`�4E'dd�NI,00� �Q_�_��?���1x�� �?���u�O��i� �֟��^]L=z��Gljӂ�?)h�կ�u�O��i� �֟��G�:��'�^4� ��O�E���T��@����?c�����b�5k�>Mg�Z�+�Y��f�E,�d�'N7 ��W����?���1x�� �?���u�O��i� �֟��W5hK�;�*���P�o���O���;o%��ìkg�ZM�� .v�''�^� ��;�?���
t]Z�P��1i�'I�d��*YX0Ì�Z��Zn��x�)�}�6Ѷ����������[㍩񧃣�?Y��{3�]V%*�<	Tp��F��G�6���ծ��V��M�md1OiwG,N:�+����v���7�����տ�U���PT��jL��H�x�������ª��:3��n]%��OE~����	#g%���緃������� #u� �+��� �Ix�IH�����xhl%���Y�y�U���_��s�z�������?h/C��cOv�0�^�20��N�#��E3v����	_���7��x�]�<b��mB�g?�*��A_^�O����+}úM���[�Gie����Ory=뢞Mަ�S�/�q��� �?��?�|3���/���JD-���F��s� NI����W���#���؛h�KEP��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���<���_���?�+j��I-~
��_�$��� �U���Z��'�D���ϴ� ����q^"� �V�� J���~�(?����q^"� �V�� J����mt�?�a[��]�m��(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh�?i/�7_��*��I-~	��o�$��� �U���Z��#�#��>�� �N��x���[��,������ �N��x���[��,������ Ʒ�QEu��EPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEP�~��n�� �U��Z�����O�M���������_�5�c>$wa�gڟ�I��8�ثs� ��u��_���I��8�ثs� ��u��]X_���(���(�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� ��O�M���������_�5���I� ɺ�S� �WU� �9k�B��gĎ�>��W�	7� '�?�.��οX��s�	7� '�?�.��οX��3��
(���.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.y��)� &��O��MW� H�������?��>)� ة�� ���^V3�Gnf}�� �� ������?�Yg_���7� �� ������?�Yg_��Յ��o�(���1
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��<���� �s��� b��� �r��e~��ҟ�n?�T��Z��+�#����_�	/� '�?�.��οY��g�	/� '�?�.��οY���3�QE�bQE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE y��+� &��S��MW� H����������*ة�� ���^^/�Ge
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��<���� �q��� b��� �r��-~�~���n?�T��Z��/�#����g�	'���?�.��οZ����K�N;��)�� �e�~��F�f5�16Ѷ��걈�h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m��� �h�KE6Ѷ��,m�m-X�F�Z(�	����Q`mih��&�6��E�M�m���濴�� �8�U� �SU� �9k��?i�7���)��G-~כ����Cf}�� �� �����w?�Yg_�u�+� �� ���'��w?�Yg_��ц��o�e�+��(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.���� �n?�S��Z���� �c�M����z����_��b�$vPٟm� �$���I� b����Y��m~I�$?���I� b����Y��mta��gW�
(���B�(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(�6��?��>+ا����_�� ������o��o�#�� +��|H룳>�� �G� ��x���;��,������� �G� ��x���;��,������3*�QEթ���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE���(�P�QEj
(��B�EQ�X(��5QF�`��(�,QE��6��?���+ا����� W��Lɷ�V� �OV� �9k���7�:��Ϸ?����r$� �N�� K,���m~H� �#��<I� b����Y��ta��̪�CvѶ�Et��n�6Ө��7mi�Qp����(�
(��L�QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
*���ei�W�f���#i>�i�d���^������%�h;F��iŸ_��$Y�n�TJ�c��J����tW�� �ZxC�������� ��� 
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`����!|T��m�6����v�o��۫�8��̠j�IE^CQrvGiEx��5��?��� ߈��u_�7h_5y��.�P�xa3�]����`mv�ڦ5!'h��NQWh�*(��"�EP
(���QE��(�,QE`��d�y0�&Ɠb�ڃ,p:�Fڇ-�Ex���I�.��tO�\/ފ{H���sU� �<!� @�o��C� �kmO����c�(���t[
(���QE��J��W����{%����K��sZF����TJ�a��Q�)l{
I��|N��l�aA��Hȏf~���s�IV�NI��M;X�h��,� k�>]Ic���Ž�l��h����'U��4�k������J�MT'��d�.̷Ey_��h��������$�F�dkDX��q�H7)����kO�;[� �� �ڏmO�^�}�l�������:,���
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(��k������ g��YQ�]6:v����3��^]}�SxkL�{{���	����$u��.Eg*���e�r���Ex��5��?��� ߈��{��E�iZ�
��):,�
��ϭTe�,�E�F\��*�`��(Q@X(���EP
(���QX>2�t~�F�.��j���=23��;�0�}�y��5��M�� ��?�v�uaf�T�%t�l���A��<-�-j�K������e�6��c�I�}+�j�%%x����(�+Q@X(���EP
(���QE��(�,W	�C��������j �l� 1������<���x/���<I�ç�\�'��#��@�=Ǡs�v�~k5V\�Կg.^kh{�QZ`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP<���� �n��� b��� �r���_�7�5� &��_��=[� H����V�ꣳ>�� �E� ��x���+��,���=��� �� ���'��W?�Yg_�u��2��	����WI��h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�����o�x/[���e���;\g�?P@?�o�o�� �[տ��o� �U��/FiN���?;�Կf��+�g�q�� E5ymz���� %wL� �3� 覯7�TwV�>��g����r[��l��̇�\r�5��O��j,�eX1JF7��*�\}��������8�!i�.ɧ!^$p���l\Tyd���'�Y�|(���|5բY.tY�&ś+����5��{o�X��Z�&���e�E��� �F�9k�����[��:LK2��������nQq}�QjK����:��M5����`F�_?W�_��~�g����Ҿu�:��/VwC������� �/@� �}���Z��X�� �/@� �}���Zگ~_<h�(M�m����6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6�}F�-CO���C�4M���EY���������9z1���?8dP�0"���I� ��S� �s�į��t���^��$� ���� �9��bW���*��L�q_�_����6��^�扶���Pm�m- �h�KE &�6��@	����Pm�m-��c����Uo�"����c�ξj��l��k���� ��|�^%]*�y���i� ?�Sￆ� ��� �
 �z�k�� ��/~(k�B�I�$ZZg���z��:{���e��%i� _���F��RO�����R���]o�3�v׀~׶Q�{��O|s��H��W��ߵ���:���b���_�?���>V��>����� �o����z���w��
�
	�t�-��-N���냕�>���f����3Ge+!#��	���|h�Ok��$�ϓ� m�+�U� ��� ��/�J�p�r���w��KZn+���i�{��*\�4-2[��HG(\|6G�[��^��I�
<qJvnS�~���_�V�#��Ү�J�;69�����<C� a�����)�M7�,:��|�G��1/��6k[�H5�5e�a&Cǘ��<ؑ�^�_��W
�%9�(�oJ
1��>����/�/E���[����|�#YC8C�Y��g=p1��
�o��:������B�Z��ԕJ��.�e8(��H�¿���6�Z����H��yЫI��_3~K0Hc8�������V�:�Յ՜���adb*xa�Њ��cEUUF �|_�B�P�׏u-O��6������UI@��)=�l�{b�g%*z/�ã'Q8�S���`_��\����ݚ�2��J��"���Gc�Z��?շ���� ���E�-��I���n���&ve8_P6�~������k�Rr�yn�9$�f���ߛ�t���^��O���]gV�.�{���ooYd.��O�z�ٿ�I���u���,5� ��%���6�u
7�P���#��a�sJ.��ǥ�i&����$?<]�M�kV7��<�k�$�1݀?���^m�� j��K�'P���NwB��ãc�~�k��u����[Ȅ7B
�_��6�����T�Y4ߴ�����
�N� ���Ss:,Z���n�>�8�p;�A�@��c�]s�V�;�)1�G��_OM2[�$�6��R�ǰ$׫�
o�<�G�n(�>.|]��]�!d[�Z��Y���
��+	+K�X��c���'8�����ό5-bvb�JD*O܈p�>��k��MH|��`y��Y=��)���'��s�kc���0�w}N/�����1��ºƢ����,y8�c����_�q�_A��J�y/�i��cvG�=�~��kվ0[���� ��#�27�p�E|�Z��9�X�:2�A ���k/��({jW{����|M���|3����}̄����g��;��Rj���x7FՏ߻�I��a�Pk�_��8��5	ą�m�ڮx����'���&����vcB��{#���E�O��k�ˬͧ��Y��%!� ��r@��՟��`�c��kZ>�wqJ!�\��W ���;q�u?��{|MxFeA�S���?*���^���ýOM�_2�v�g�v�H��kR^Ǚ|[� _#X�~ח����������/P�J�mcEf$7^H��ѱ�c����־����,m�-d[\F��"�ea�#_�K��M�
�JTj4��M�?H!�'�$����Yzy��"�bTd����1�k� �~��,�|'� �W�Nx�_�
E�DP��� `
�h��jt՗�j0�����^/��[�&RK�kv����A�O�H��Fk��t� ��T��Ԁ�w�r��*}Gp{�ʾ{����?i��O�	Ψ���'�I���Qs�mRR���@]O�~4����/;|��������]�9W��P�k��+��Þ*���$E��Ѹ��v���{=�⦬ˌ�v?8�ym5�'2h�dw��	������/�5�_�� �ױۙ$��2��7(ۺ0I�:q\ω��c���_�׮~�2m���&>����?�+i(���=D�"��������'�Ǉ�?km�]�)h��Y�4*H9<W��p�|Q�l/4� x��O�v�mơ/����aҾ��
�����r5�� l޻q:S�[���A�;=�Ϝ� �fx����o� 3�U�g����SᏇ�.&����t�N����������_���$�� �o���p��W��G�V��� �#�H���� G%|I_n��_�Gu� �?�::���k� _�Cz�������O� �ow� a?���k�� d��&����� @J���<i'��\�YW���8�N
�	v�G�+ӕOgIK�~H��\�⻿�󿌟���w�'��T1��e��8����T�|������|_��nיkh�`������ӿ5��Q�eDY����cG��K�l�����@�T����}�s���Z^ͨ�C�z�ş�ߍ�N��&��4+&L�FN���x���_
�&�|N���`_���Dwv�����=T�?_J��
�WW���|��~������5��O�N���E���k"�'��w ?KQƣ���O�u��=���:�}�s�ey�o�^>�W����ω� ���n����2�-�5#��⽆��qSVe�N;��Ν������p�I҆-��͓���&�?�^���N�^���!�7h����'��=�3���<x����1�пey6�R��e0� �MyXX���}�CѮ������� ���:�^*�� ��Uh�r�yjf@1ϥvZ��e��T�f�vɾI��� =MiW�?�����S��/�`����;d ?A�� �F�Oe
�'��|#��{��������X�q1H=:_a�}�������dH��B$\�n��z5-m�B�G���A���kP'���"H��.@p�����K*A�#�q�,�� ԓ^+�$��j�'Pl{��k��޷�|)�������#���g�9��Z�9i�ۥ� τy���o��o�_�桨jRi�t���� iH��/8ʆE�#=���f���O��Zյ\�pbRG@1#��|�_m������|+�$ԚGti"�I3��S��9�
�j���VtV�;rh��~ x?U�?����e
ri��Q�^N)[���� 	>���z��׊<K��e%��8~�<�_r��dP �ך��}���I#����QO������~ƶq%��2oX\�.=���e��<Ï���Ɵu�x����r�`~��,����ʒ}k��4-��R�'��q1�b���� ,���mw����C��+v���]��67����kͣS�O��ӱ�R��'�o��? 药hڍ�ݜ,�q�S*n �2@�^��x�g��	{K��8$������P�������� �C�����J�v����CCHZ�?|���e��-cObַI���z>��?
��o
���Y[�Z��2�$[�6c�6�:�9�G W�~�wo�2�I>U��>�P� 2k�kՕ�~�� 3�M����O�^,�|I�i��Z�[)�"� �<n�9���,kz�t��o����T2��k�dFP�����Vw���*�&� ���B��xr/x�F��v���Iz��@��5�P���g����QQ�և�����1��^�L����{�Gm ��Qԓ���W t��^k}y�^	��&����$��:��9�Ell`�lവ�`��qă
� �G�R��?���B��3�GA�W^>�j��3�.x�Տ�?gO�W�4�}[���-c�a�o�4`�Cz���>�>�_~ͷ&�������I�>� �c�+�j��Mʚos��yf�7�wċ�~mB�|���kj�Wǯ`:��"�A�?��|KՒ���\�#�O�c �¨Py9=['޷?h�I�o�W��!k=/�%�����ݑ� G��WT�������|{�3�
���"�M����ʨ�r�{O�� f7���֟�KM�*�o>݂ۆ��*⽾�>���>7��S�zN��O���E�xn�#�0�d��
O�q��
�{���m�u�T�6���,2�>���R}�r?Zbb�˵�2�.i���~�|<� �ß��� �k_���	���D�:�� E�^���'�G���S�ۜ�+�>:�\?�H`����a q��GY����ҽJ�<���v��]�R��"�M�l���rbj8�%՝"�'~���z�Ǎ_P���{%���1���6p��*�`�����M�Gğ�$���-A��o�O��vP9 �O�J�3ҼO�ĭSឹ��<���\�H~I��؎Ƿ�"���w�[A��b$�n<��\�a1X�N?��}:TS�*B��EԔ�+�X�'�4x�S�W��v�W�����O�H��"` 0 �+��� �;��5�+)F����J%���в��$C�^��]r�Q�g�+8b=v�q��?ޜ�(�.����f�����܏7������cM��;�~T��R�OF|u'��'����W�_��{y�-z��"Հ���3���>���8���oX���y�^�e�������?ʾ���t԰�K�Ƞo��i��ʏ�EsQ�y9T�t7�j1J���*|4��\Ykֵӧ(���H��@b�۱�G�_�C�����["�cB�̀*�(��fHȭO�j�g�K�3b�_��?�5���ˏ
��O�G�>?x�-!�Y,l���
����<Z��
����N3��=[.�G��d|S�|R��\�7���e1�[+�2��t�T��G�OÏ��gĭu<��B���^�C�c���� �z�m�,^8]�}�2s�n������i'�>#i�d�e����Rx�����MF������|˭Mr�H�u>���� �n�����e ]����9V��Bk���������E>3�.�EYb��̥��FF����W�����S� %K��a�mO�����Ň��� [n��7�K�Y��U��d2o�&+�p�~n��`�o�|^����eu�lJ=�.!�s���z(9�^�;�ڛ��{`��,{�2@�~��O�^������#���.q˷v>��Tb�A�=VJ������������y��}F
�$$H�9=s_@~ο�+�j/��r_%�"��f�"���y �q�A�'�/
C�sS(�w��;�,���|6�6_�76��u2}���Ya�(�Q��cJ�R�������}���QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E̿i��6ߊ� �)j��G->��~ӟ�m�?�Rտ�Z�}k��[3��$G����?�R�� ��:�u�ȯ�$?����?�R�� ��:�u��� �O�(���2
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
��?��z��zM� �Ӭ�� ȷ�פ����� ^�ҟƽO���?ٟ�J��� \g� �M^Y^��3�^�?�� �)���� ]���g�TQEzǜQE |��a� ��?�� �+�:�3��� ���9� �W�u�O㗫=h��?B�� "^�� `����Ǜ�����6|���cw�������� ȗ��>�� E�|��I诤�Xԥe";�㹌��S� �)��oğ�Ʉٯ$s�	t[����M���<=up�*_Z�޸B��YO8�Z��'�?�k�����_0x;�RxO�Z^��l�R����~##��-'T��4�]B�U���1,R/B�f��!:w{�V��=6>~���nu+�n�<o-�̧t�O��w>���&���� ��N�y�����P��_d�v��;�=�w�T�=+�Zn��=WE�y&��
�6�h� ?�S������f�G���Q�@��wLe��1�o!��G�	�=�>���π�{o��֢�mq���e4g����~�����Z��߳����	�ߜ����S������W�Ys��3�u���������z��i&0H��&)�ї����O��� ��MJ�fC��ڱ�C&:{��� \��CR�n�{�쯭䵻��I��R;^��#j��}V�X��k#�m�"�?�7�Z�j6��x�$���_�Ӷ��~v,I��?�|�l�n""F��H�W���W�߶��X�{VU&8����8�e� �Z�e�I{���1�
ݏ�c���#W_�Z�VM� ���~�:��d�MW��$cj=��TzӜ
��~6�����@ז����g�u�� x`�5�W�*4ߚ<�՚�s���.7��U�kK⏴�)��O�~��u�ӌ�c^�Xz�����N�;����� ��B�$�z/��ܫ��ch�L����
(����O������%�?��ĸ���^!��
�Ϩ5��;t�����?շ��4��:��gb����ޖO�m�5�O�^�-�����I���u�߲?��ڟ���� F%x���� x� :����G�O���� ���/�����b�z����R� �}��ro�k�O������K�A�_�ɿ���9?�7֞'���X����� � �կג� ��_B�H�k?��$�3�4����>[W�_���Z��y/���������7������%�}Y����x[.��g>���u�������7����_�u/
(����)��?��� ���0W�W�~ӟ�H�����͈��~f�?����}}��,� �+O���� e��k�o�g�IZ����-a������lO½Fz�x'�}� "��� _�� E�{�x/�}� "~�� _�� E�k���~h��E����}��7�I/����� gj�R���� $��_��� ��Fi|�R���?l��=� _r���}K�`ȷ������ @������e�tGM?����Ͻ�� �0���xE��¿l�oxzߝ�k#�L����+�~� �0���xE���l
���>��?|5,�.�m�'�٭x��� ^� ���+�|?�� b�Zv�9����]�����?�������%z߅����^���5~�|7�[���@t�:{�'����~��?���_��:ߧ�sZ�/�/��i�����
��H.�:����0;ׄW�_���A����J�v�k%�#z-�}� ����M� �oy� a?��k������FO�:���#���~G��~��>��� %[������S�� �X���}�L��� %[������S��Y���}�^N���W�zu��/F}�_�X� ����M� ������� %3�� ��� C5�+���1��/��n~�� �X|?��� D�}�_��� �X|?��� D�}�]t?��g=o�E����(�k�jӾwIw3���ɯ@��� �i�������L�_�� <A`��.�B��۔�DV��
�MB��u�1��ƾC��ýg�ήlu{}��a����_U?��+̔%I����g�C鏁?��>Ɖ�m�]�7G2�%ҁ��gH�ҷ� h����C��dGǡ�3�¾@��>���݊��8�|�~ ����_4i5"��KI�c2����� �k�Ru0�O�9TU:�Km���O�?|E���w�MWC��$b��HR"��G=z���})�%x��Q�ā&g�]�c���{���8uI�]�k����K�o��f�����R���-�2H����i��X���� �7� �3� �W�u��i�
��kq�Y����;dN��"�eе�GN�J�kq$,�,Ez���h|�
��Ee�G�I�1 �m��*��r�I=��G2�kt{$� |]"2?ōiц
��A��"���c��� �� ۫�:l�$<�2���p��J�{�q*��3��S���]�c�� i��<�o#����noO^�������Ě�Ą�%��ؓ�K�_|�GŶ6����+���B�I]���Oc_|P�$���
�Z��
6��Ԍ�ǿ�O��/�q���}�_��_�T��xC�ں�}Fsaw���-�?�=�� �����j���?�=�� �����j�
���k���%���-���+�l����`��0���^�������3�ÿ���v�� �Q� ���������� ��<<�n��U��F�¿��"��%��S4VP�,h�n���#�r�7©~!xv+�6/3Zӷ4q��F~�}x��Gz�|3�g��
غ��I9lx��E�����(N�473Aj�]�	L����� �ڽ;�Y�us����Z�-Ƥ��n0�s����1��]W�g?xV�.��mR�3�R�H�� )�A�N����Kٷ9n� QT�ΔV��~ �2�����nN�=3�ϗ �[�?�k�� x_R�~�q���mw	�V�ó)��:��g��/Y����C~���WO�\��5�s�h�U*܋��|3�?�z�Ýc�CH�T��5��tS/���W���ajSX���3n��0%k�x��0��������4���X�Rs��qQ��FO�]o���C��	uc����9[���dS�§� iB�[r�Y)�o��f���uu�O;M�=I6���k�q�$��`=I�[���IE$�9[m݅QLG;��v�<){����e�R�&)*��BE|!��	�~��4�V٭�?�^̧��Z��Ox+C񵐵���u�����/�k��w�ΚU�+��G�^8���I���L�mx�n�U��/q��W�/�}����[��c�m���f�������<<��u�[�9�qQ��FO�]_�~x?�7���j���n�[�u>��O� �ӧU.[�9�o��g/�"�V��;�w�����*���_g��p3�9=���I���7���V�E6�\.�A���
ԮW��1��{�{L�t��c���,{T���a�}kyE�r@�I9sL�
F�#7��{췯��&K��V�����`�_r�\��^�� ���������� ��8�;� �{� 'n��sң:rRV:*U�D�:�x����/�����`}�G �m8U�$���NI5����s��$��m�#�i�[U�L���:/��=:���� �q�w� B��N�� ��uhΤ���Ռ#fxg������1Z���ˉ�p��:����g�3�G"���������;s� �+����MK�Ӭ"�,�cE��Tt$���馜`�.�=F�'%��G���Is�~�Ěu�I�_?�?��[�O9�V<��H�ω�#�"�2:���pA5�E41������U�uX���W����f�>�p�Gis����|�T��X2�� W��.���tզy�� k-sE���4�uƍv����+{����|
��u����v3�4oZȲ^-�'���B��1� Ԏ��i?���4ہ,��:��|���� !��ֽ_L���[��-a���a!�"� ��Ӝ�����f�W�տ�����Y[Eok �8�`*���_9~�_�.���5���P�1��ߎ�Z�J����Ug����z�߇?����]6-6��=oN�m�f��,k���<{q_Ax��{�W�.�]9��9y���}Ws�5�c�'�*��d�mV���P)��E?�rX?u�2�NK�G�<k��x�zRh�x��k)��N���@A���������,�k2]!\�3}�d<����0*߇|/��OO[-O�O���q�����ܒkFhR�"�n�E*˜d�]��*ov� ��唹�[#�o���v�h:]��7V�4��1d(P��N��}�� ���������� ��8�;� �{� 'n��qK9Iɵ��ЌT{����^�k����2�kiB���6�Q���9���cs����;}Bݖ�Ub�Bp1��L�W��g�σ��E�h�?����/ڦ|<3���U��)�/���t��\�����c��a��Uhʴ|�zRT�~����W����U���ֽx��˶������2�v�J�?�~� н� ��?�r��g��/���� ��S�I���J���I���n�o����V��V��QG���� d�x��v�\GƏ���᮳�Y)��Yb�uvF
��_
���yy��?C���\��?<�a��O����V�n�\F~V��"�u=�� <������Z]KG�K��E ������;���
������o��
�ԇ.�=9F��}w�*�������E�CP����ܰ��C�u�zWG� ���������� �cM��='P���м��iVh��e�mu ��&#�Es�&��B�ԍH�t3� h�i�O�Z���5��D�}�{��I� �������x�X��5M�W�e�T鸀 �\�g�8�;� �{� 'n��eR��K�ظU�#ʏ=��|Aa����E�\����H�0v��s�Q]��E�]_�,�푝m���`���d��jE�:�=�T�<?�у)�m��J�Y���!�5�)���e#z��%R�#�O��MF�2��گf�U�C�������]
�l�e)"�9ـ��d�8<��`տej��sԴ�c�oi:�c���k��O|?��1�5���1���O��V4�ԋ�Ƶ*�Kk�,|u��`��V��Վ)��Ek �mNDg���Sׯ���!�Dq�R�+�7�~�|m����٭��C#��9�^ao�'�*���է���K�}���C��N�صZ<���A��
|N�h������n/nYB,�b4;���p:s[��7q��h3u<���C�J���A��3�C��v�X�D0��0>����MhWl�����k���42[��J���Y`�G�{'�� �_t;�/S����IM�k�!ʀU���9��{� �~�O�׍{{i-������� �*O����~̞	���]<7z��ܩ�J�������9�N�ZrѝS�NkT_�;�x��K�x�Y�e��SN�3��%��9��^3���zM5UcP��� 
uw%d�8ۻ�R�ih�} ���F�u���S��F���En��y�z�5����6����CT��.��=����	�r�����W���R���� �l�q���l���_
�=�^M?P��<�1�t������ |��6�rNo5+�7׌I3L�NO$�?^���_h^	�6�&�
�P��B���u$~S���dʒ0�GP+��ʐUci����ϯ����'T���X܃4da�����5�|>,Ѯ41���jt������{n���/����>��!�a��۱�S�+���@�d��7}�Z��h��+5�#V�y��p�.g���/'��!j���w��`���J(�{��{��� �����n�"֠0��1췶�|�񞥽�;���o���d��t�k��wt�Y�'�?�]�h�Z�X����-�e�K1
�d�p|U�Ik�> ��s%��wP�oXX2�$:�v?
����U��iz��ڬ. Ž�pʐG t5��8�;� �{� 'n��*��V�l�J�iݽ��V�r�z��Ì�S$�����-7Z��4��;�n,�M�x��3��\�3�ÿ���v�� �V���w�N��h��G���E
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��<���� �l��� b��� �r���_�o�9� &��c��-[� H��ʸ1;�������	� ''�_����ʿ^6��G��NSĿ�)\� �e�~�V��Χ�&�6��]&Bm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����P��Ӌ� ��c��-[� H��ƿ�� �s�M�����Z����_υpbwGM-����NSĿ�)\� �e�~�W�?���SĿ�)\� �e�~�m��|u>!�S�Ѷ���N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m��N�F� m��m�2��?���,� إ�����W�!�N��c_ş��o�#������n��[3�?�$�����R�� ��*�z��_�$�����R�� ��*�z��|u> ��+s0��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��(�� i��6��?�)j��G-=������m�Rտ�Z�{���:)lϸ� ��_�r�%� �F�� K,���m~B� ��|K� b����YW��kC�"��&�6��]Z�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�e�O/�c_ş�uo�#������?i��6��?�(���G-=uň��gܟ�H�9_� أs� ��U��_���H�9o� أs� ��U��Z����
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE ��(�QE�(��,EQ`
(�� QEX�(��QE̿j�6��(�� �E-<U���?�l� � �Q�� �Z�x����ޞ���G��9�ءu� ��U�	��?��?�s&� �B�� Kl�����=��F�Z+b�F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih�0����1��������K_�
孺6���O���Nc���(]�m�~�m�ǯ�#���ω��P�� ��*���)|$�q6Ѷ����M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� �ډ��.ء�� ����W�E�Qɳ|\� �CW� �)k�ݮJۣj{u�� ��|M� b�����W�5~=���M� b�����W�=kK�&{���Elf6�u �)�Ph�Q@
���������/ʿS<�H|g�����5k{����,j��UR�)�^�u�ga�S��,6�u��N����)�PE:��h�Q@Xm�(
C�x)� K�Þ�l��H>�ܾ~��o~%~����F�J��m*��+[�S����壷�^x1���jᓻ�ж3�I��x�ú��~�e��k-��z9�ɍ~.�0�믃?�_��k/��j��$�7�>�i!q���?گ�J�I��4�� f�(����bx����<U�8�o��>��ˑ�{�"��T~��P����Y�wk*�mq��r��2{EM]�EPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEP^9�b|F_�?��_��M�zL����q?�"#�"��ױ��_�S�w�O�_x�m�k:���ҩ�ɷM����0#�*�;+�nQ� �!|8d������d��C����\� ����6�k� �t�-� �O� x�X|��Z��{����S#���_JW!�W�������/
|D��uφ����˭���YQ�S_����|�[m�࿌�rޙ5�N�#�n�d� ����4�������-��O\M�j�
/c�-:f�K��m�1F��ᅋ^�1�X�� �տ�⥷���1����t}:I��C�����FO�R�/�
�9���_
��k�Ɵ�u{���K{�y2^�Ҹ=�cW���ƣ�i��ݍ�m����K��kQ ¢(T{ OE��QE��� �~��&��Z_���b��nb�=�ߕ�_������q_���V��	��#�&�N��Y�}k/����������
��?�Ưٚ�Wӭ��!�c��je�������fS
���_�H��C�� �?�ڕ�}[�r��A�-%��H�˔��F;V�޶"K����h�]BQK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R�����i��w����Z����,q��f'� O�@�_�Y��b��%���^���DSȉ	���?����M}���
��#��~^|ѯ�� �P��A�A$����P�	Xt�fU��=�r"F���5��\�ww7J�(����Q@�wP�Ѭ�ȥ92���=A�a�Cw� �� ��I"�4~fX�~Ѣ]���.Gl�o�_�u����i� ���-|u��y�&�H{���/>���Bq����3�$�婁�����coyi4w6����I�*����jj�c�	?�F����}�\��!�hX����m9�����c�Q�}ӴWR����h�mm�!�S��6�.h�mm\�N�(�(�
�X��u/�)g��^U����d�-���UiK�%��+9�Ȩ�����$���~��m�MR�����[�2(����� ��ב"zW������g��1�Z�ƱEJ#EUP:  ���(�� (�� *+�Xo�����;�y��$2�du#H<G5-�[�c�:��C��ִ8f��SJo�#�E�3m�Ա�z|�3Фnz��������x_I��v��F�k��g!�u����;+�� ������]��#��@/�hnG�#����}%P ��c'�_� �&� jyt�����-�x$W����tv�n	i�z�G��?�+X=lȒ�~��K�Ѵ��bQK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK����և�/ᖽ�oM��U����Hx������4�W����� ��v�	t��߉#��
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,��� ���o��㦋�o�ϛ��n�,��շ&���(҂81�(n��8 �Կ�S��~
�
�r�\WS��_�w�� f��'�4�.{�W����oo\6^y�T��;W��EdhQE QE QE QE QE |��q~�zO�i�m9R+O�+%Ɓ��% n�C� <���![��?�?�N� ��V��>#�|���i%Ơ���@�mP-�G$�E#u=�n�ƿf����
u�����>���*���
�ڽ�
��x�q��m�)�gqq�_���L�ۨ������-O�*;T�R��k�qk#��)?}F>���G6�JwW3�]�m�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�QK��� �R��m %�h�@	E.�6�W�~�ߴ��+�'��>����7�h�I|5��g��>��8��s��׆��>�|e�����4���O!��ƿ��x�8 �����G� ��~��^�.��k4�^�D���G٦�ݛ,p��ĥa�s����e�{���٫�^���j~��77�\.�{�[e����p*�o�Tq�1�"�EUT` : +��s��A�O�m��];B�-����NHQ�f?��Ifc�,I�]%s�Q@Q@Q@Q@Q@Q@Q@���S��A�/~7�)����/۵�'N$�����/#��@�������O/���G���7����=�8�����£�t(uye�p_�Ddu�0U�A��� ���O��$�_a���`�څ����D�6����y �?�}ʌ�D���*+��7�
�~�]��|_-���.�,m�!��G2�:	02яv^2�������R�4m4�J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J��}�����Ox�R�I���L�7S�����@�_x��� �#�x��:�6���e���l*�����$� &��>||��� (�ե|>�}�]A�h�-��l�Wh8k��
 >�ڻ��h��F��>(�I��� C����	�N�kZJ�O�����@����Ǯ22I-��7�ξ��~�x7°e#����"�=�� <���z*�J���d� �ɿ
�6��C�K�Jp��`0e����*��z�a��QE�(�� (�� (�� (�� (�� (�� (�� (�� (�� (�� ����� �q�O�����'���¿�o���}H��]*��`J��0 �>~�� ���>��k�7D�5/Z>�4+��h��Oe)ʴm�۝��
���ï-��?f� ��>o��G[�L���o��,�� 2���S�l
�Kx�O�$�w�� [ ��������E�C>'~�~!?~
�<A�
W�~�~Ͽ���_�cG^ף}E=·{�o��ދ'p� i/�z��>���Y<=�Y��N;�hLKn�c|2<m�Ґ{t�˟��	!�υ�����f�u���9��Mk��j��d��e*��w���RhV?Y��?��UO��5������6O�\O<c֭1�Dp�)8F=K��'�!�S�0��,Vox��Q�	�]*c�_C��`s� WК�I2OU��EP�6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\mꫪj�:�q�jW��}��y�]]J�E���� =�r�VԵ;=O���.౲�C$�72,qƣ�31 �k��#�
����w:g�"?|B�A=�t�[զ#2��#N1�WŶz?�O� 0��YZ�|&� [.�?@���932�O6A�J�$�>������7�������[\����e�`lu�pӐ�N�3+毄߱����� �x��N��xj���
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��<���������E��<-i�Kl��Pom���tÁ�v��
�P��Ѽ]��3ĺ[�6���N�� yw��`�ƾ��� i�q�&�� �ι}�=A�aas��z�a*}w�aj���}�e}o�ZEuiqմ�9�p��zÂ>�5~,j�_���M�{�Z� H����
�|�V?]诐>� �T�|D�a�u��߿G�-
ǻ�&��`{�_¾����/�A��g��I������Pr�j�#~�Z)��R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E-6GX�gv�2�� S@Ex�ď�+��x�� ���F��y�ecq�ې}Poq��_$|P� ��x3H����W�S�o5���� x"�����4��~�מ�W��>���O�|e�xt�2%��ᮥQ���G��T��1{�Y~׿�5��w��k.F(��6ɭ`�<bK�%���Uһ��F���6�]_ⷍ-<<'2��Ŏ�!�RDj����Q��,v�� ೚e�����BMJnUu�� ����=:9Z�='���g� � ԼE>�o�g�̊�^'O�aRG����8�r#����~��
��C
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
�$|
�-��&���k6���P�I���pE;�X�������w�>j�y�6�uX�8�,ש�w�
��ZU�o�$�=~ߥ��� ~d���Ŀ���o�w�_D_��셉��)^U�/�$���3O�{i�����2S�adli?�SٻW��DV���J���C�ֺ�/ۋ��x�,�aA� �ע#�>
import sys
file1=sys.argv[1]
file2=sys.argv[2]
f_in=open(file1)
f_out=open(file2,"w")
a=f_in.readline()
while a:
    a=a.split("\t")
    a[0]=a[0].replace("Chr","")
    a[0]=a[0].replace("chr","")
    a[0]=a[0].replace("Chromosome","")
    a[0]=a[0].replace("chromosome","")
    a[0]=a[0].replace("Ay","")
    a[0]=a[0].replace("_RagTag","")
    a='\t'.join(a)          
    f_out.writelines(a)
    a=f_in.readline()
f_in.close()
f_out.close()
# ************************************************** #
#              Read exposure and outcome             #
# ************************************************** #
read_gsmr_trait = function(file_con) {
    expo_str = scan(file_con, nlines=1, quiet=TRUE, what="");
    outcome_str = scan(file_con, nlines=1, quiet=TRUE, what="");
    strbuf = scan(file_con, nlines=1, quiet=TRUE, what="");
    return(list(expo_str=expo_str, outcome_str=outcome_str))
}

# ************************************************** #
#                  Read GSMR result                  #
# ************************************************** #
read_gsmr_result = function(file_con) {
    expo_str = outcome_str = bxy = bxy_se = bxy_pval = bxy_m = c()
    while(1) {
        strbuf = scan(file_con, nlines=1, quiet=TRUE, what="");
        if(strbuf[1] == "#gsmr_end") break;
        if(strbuf[1] == "Exposure") next;
        expo_str = c(expo_str, strbuf[1]);
        outcome_str = c(outcome_str, strbuf[2]);
        bxy = c(bxy, as.numeric(strbuf[3]));
        bxy_se = c(bxy_se, as.numeric(strbuf[4]));
        bxy_pval = c(bxy_pval, as.numeric(strbuf[5]));
        bxy_m = c(bxy_m, as.numeric(strbuf[6]));
    }
    return(cbind(expo_str, outcome_str, bxy, bxy_se, bxy_pval, bxy_m))
}

# ************************************************** #
#                  Read SNP effects                  #
# ************************************************** #
read_snp_effect = function(file_con) {
    snp_effect = c()
    while(1) {
        strbuf = scan(file_con, nlines=1, quiet=TRUE, what="");
        if(strbuf[1] == "#effect_end") break;
        snp_effect = rbind(snp_effect, strbuf);
        print(length(strbuf))
        if(length(strbuf)<14) print(strbuf)
    }
    return(snp_effect)
}

# ************************************************** #
#                  Read SNP instruments              #
# ************************************************** #
read_snp_instru = function(file_con, snplist, nexpo, noutcome) {
    nrow = length(snplist); ncol = nexpo+noutcome
    snp_instru = matrix(NA, nrow, ncol)
    while(1) {
        strbuf = scan(file_con, nlines=1, quiet=TRUE, what="");
        if(strbuf[1] == "#marker_end") break;
        expo_indx = as.numeric(strbuf[1]); outcome_indx = as.numeric(strbuf[2]);
        forward_flag = TRUE;
        if(expo_indx < outcome_indx) {
            outcome_indx = outcome_indx - nexpo
        } else {
            expo_indx = expo_indx - nexpo
            forward_flag = FALSE;
        }
        snpbuf = scan(file_con, nlines=1, quiet=TRUE, what="");
        snp_indx = match(snpbuf, snplist)
        posbuf = rep(0, nrow); posbuf[snp_indx] = 1;
        indxbuf = expo_indx
        if(!forward_flag) indxbuf = indxbuf + nexpo
        if(length(which(!is.na(snp_instru[,indxbuf])))==0) {
            snp_instru[,indxbuf] = posbuf;
        } else {
            snp_instru[,indxbuf] = paste(snp_instru[,indxbuf], posbuf, sep="")
        }
    }
    return(snp_instru)
}

# ************************************************** #
#          Read output by GCTA-GSMR for plot         #
# ************************************************** #
read_gsmr_data = function(gsmr_effect_file) {
    trait_flag = gsmr_flag = marker_flag = effect_flag = FALSE;
    file_con = file(gsmr_effect_file, "r")
    while(1) {
        strbuf = scan(file_con, nlines=1, quiet=TRUE, what="");
        if(strbuf == "#trait_begin") {
            # Read the exposures and outcomes
            resbuf = read_gsmr_trait(file_con);
            expo_str = resbuf$expo_str; 
            outcome_str = resbuf$outcome_str;
            pheno_str = c(expo_str, outcome_str);
            nexpo = length(expo_str); noutcome = length(outcome_str)
            trait_flag = TRUE;
        } else if(strbuf == "#gsmr_begin") {
            # Read the GSMR result
            bxy_result = read_gsmr_result(file_con);
            colnames(bxy_result) = c("Exposure", "Outcome", "bxy", "se", "p", "n_snps")
            gsmr_flag = TRUE;
        } else if(strbuf == "#effect_begin") {
            # Read the summary statistics
            snp_effect = read_snp_effect(file_con);
            snplist = as.character(snp_effect[,1])
            effect_flag = TRUE;
        } else if(strbuf == "#marker_begin") {
            # Read the SNPs
            snp_instru = read_snp_instru(file_con, snplist, nexpo, noutcome);
            snp_effect = cbind(snp_effect[,1], snp_instru, snp_effect[,-1])
            marker_flag = TRUE;
        }
        if(trait_flag==T & gsmr_flag==T & marker_flag==T & effect_flag==T) break;
    }
    return(list(pheno=c(nexpo, noutcome, pheno_str), bxy_result=bxy_result, snp_effect = snp_effect))
}

# ************************************************** #
#         Display summary of the gsmr data           #
# ************************************************** #
gsmr_summary = function(gsmr_data) {
    message("\n## Exposure and outcome")
    pheno_str = gsmr_data$pheno[c(-1,-2)]
    # exposure
    nexpo = as.numeric(gsmr_data$pheno[1]);
    noutcome = as.numeric(gsmr_data$pheno[2]);
    logger_m = paste(nexpo, "exposure(s):");
    logger_m = paste(logger_m, gsmr_data$pheno[3])
    if(nexpo > 1) {
        for(i in 2 : nexpo) 
            logger_m = paste(logger_m, gsmr_data$pheno[i+2], sep=", ")
    } 
    message(logger_m)
    # outcome
    logger_m = paste(noutcome, "outcome(s):");
    logger_m = paste(logger_m, gsmr_data$pheno[3+nexpo])
    if(noutcome > 1) {
        for(i in 2 : noutcome) 
            logger_m = paste(logger_m, gsmr_data$pheno[i+2+nexpo], sep=", ")
    } 
    message(logger_m)

    message("\n## GSMR result")
    m_bxy_rst = data.frame(gsmr_data$bxy_result)
    print(m_bxy_rst)
}


# ************************************************** #
#               Retrieve SNP effects                 #
# ************************************************** #
gsmr_snp_effect = function(gsmr_data, expo_str, outcome_str) {
   # index of SNP instruments
    pheno_str = as.character(gsmr_data$pheno[c(-1,-2)])
    nexpo = as.numeric(gsmr_data$pheno[1])
    noutcome = as.numeric(gsmr_data$pheno[2])
    expo_indx = match(expo_str, pheno_str)
    if(is.na(expo_indx)) stop("\"", expo_str, "\" is not found.")
    outcome_indx = match(outcome_str, pheno_str)
    if(is.na(outcome_indx)) stop("\"", outcome_str, "\" is not found.")
    forward_flag = TRUE;
    if(expo_indx > outcome_indx) forward_flag = FALSE;
    if(forward_flag) {
        outcome_indx = outcome_indx - nexpo;
    } else {
        expo_indx = expo_indx - nexpo;
    }
    indxbuf = expo_indx + 1
    if(!forward_flag) indxbuf = indxbuf + nexpo
    strbuf = as.character(substr(gsmr_data$snp_effect[,indxbuf], outcome_indx, outcome_indx))
    snpindx = which(strbuf=="1")
    if(length(snpindx) < 3) stop("Not enough SNPs retained.")
    # bxy
    indxbuf = which(gsmr_data$bxy_result[,1]==expo_str & gsmr_data$bxy_result[,2]==outcome_str)
    bxy = as.numeric(gsmr_data$bxy_result[indxbuf, 3])
    # SNP effects
    if(forward_flag) {
        indxbuf1 = 1 + nexpo + noutcome + 3 + (expo_indx-1)*2 + 1
        indxbuf2 = 1 + nexpo + noutcome + 3 + nexpo*2 + (outcome_indx-1)*2 + 1
    } else {
        indxbuf1 = 1 + nexpo + noutcome + 3 + nexpo*2 + (expo_indx-1)*2 + 1
        indxbuf2 = 1 + nexpo + noutcome + 3 + (outcome_indx-1)*2 + 1
    }
    snpid = as.character(gsmr_data$snp_effect[snpindx,1])
    bzx = as.numeric(gsmr_data$snp_effect[snpindx,indxbuf1]); indxbuf1 = indxbuf1 + 1;
    bzx_se = as.numeric(gsmr_data$snp_effect[snpindx,indxbuf1]);
    bzx_pval = pchisq((bzx/bzx_se)^2, 1, lower.tail=F);
    bzy = as.numeric(gsmr_data$snp_effect[snpindx,indxbuf2]); indxbuf2 = indxbuf2 + 1;
    bzy_se = as.numeric(gsmr_data$snp_effect[snpindx,indxbuf2]);
    bzy_pval = pchisq((bzy/bzy_se)^2, 1, lower.tail=F);
    return(list(snp=snpid, bxy=bxy, bzx=bzx, bzx_se=bzx_se, bzx_pval=bzx_pval, bzy=bzy, bzy_se=bzy_se, bzy_pval=bzy_pval))
}

# ************************************************** #
#                  Plot bzy vs bzx                   #
# ************************************************** #
plot_snp_effect = function(expo_str, outcome_str, bxy, bzx, bzx_se, bzy, bzy_se, effect_col=colors()[75]) {
    vals = c(bzx-bzx_se, bzx+bzx_se)
    xmin = min(vals); xmax = max(vals)
    vals = c(bzy-bzy_se, bzy+bzy_se)
    ymin = min(vals); ymax = max(vals)
    plot(bzx, bzy, pch=20, cex=0.8, bty="n", cex.axis=1.1, cex.lab=1.2,
         col=effect_col, xlim=c(xmin, xmax), ylim=c(ymin, ymax),
         xlab=substitute(paste(trait, " (", italic(b[zx]), ")", sep=""), list(trait=expo_str)),
         ylab=substitute(paste(trait, " (", italic(b[zy]), ")", sep=""), list(trait=outcome_str)))
    if(!is.na(bxy)) abline(0, bxy, lwd=1.5, lty=2, col="dim grey")
    ## Standard errors
    nsnps = length(bzx)
    for( i in 1:nsnps ) {
        # x axis
        xstart = bzx[i] - bzx_se[i]; xend = bzx[i] + bzx_se[i]
        ystart = bzy[i]; yend = bzy[i]
        segments(xstart, ystart, xend, yend, lwd=1.5, col=effect_col)
        # y axis
        xstart = bzx[i]; xend = bzx[i] 
        ystart = bzy[i] - bzy_se[i]; yend = bzy[i] + bzy_se[i]
        segments(xstart, ystart, xend, yend, lwd=1.5, col=effect_col)
    }
}

# ************************************************** #
#             Plot bzy_pval vs bzx_pval              #
# ************************************************** #
plot_snp_pval = function(expo_str, outcome_str, bzx_pval, bzy_pval, gwas_thresh, truncation, effect_col) {
    eps = 1e-300; truncation = -log10(truncation);
    if(truncation > 300) {
        warning("The minimal truncated p-value would be 1e-300.")
        truncation = 300
    }
    bzx_pval = -log10(bzx_pval + eps);
    bzy_pval = -log10(bzy_pval + eps);
    pval = c(bzx_pval, bzy_pval)
    min_val = 0; max_val = max(pval);
    max_val = ifelse(max_val > truncation, truncation, max_val)
    gwas_thresh = -log10(gwas_thresh);
    plot(bzx_pval, bzy_pval, pch=20, cex=0.8, bty="n", cex.axis=1.1, cex.lab=1.2,
         col=effect_col, xlim=c(min_val, max_val), ylim=c(min_val, max_val),
         xlab=substitute(paste(trait, " (", -log[10], italic(P)[zx], ")", sep=""), list(trait=expo_str)),
         ylab=substitute(paste(trait, " (", -log[10], italic(P[zy]), ")", sep=""), list(trait=outcome_str)))
    abline(h=gwas_thresh, lty=2, lwd=1.5, col="maroon")
}

# ************************************************** #
#                Plot bxy vs bzx_pval                #
# ************************************************** #
plot_snp_bxy = function(expo_str, outcome_str, bxy, bzx_pval, effect_col) {
    eps = 1e-300;
    bzx_pval = -log10(bzx_pval + eps);
    xmin = min(bxy, na.rm=T); xmax = max(bxy, na.rm=T)
    ymin = min(bzx_pval); ymax = max(bzx_pval);
    plot(bxy, bzx_pval, pch=20, cex=0.8, bty="n", cex.axis=1.1, cex.lab=1.2,
         col=effect_col, xlim=c(xmin, xmax), ylim=c(ymin, ymax),
         xlab=substitute(paste(italic(hat(b)[xy]), " (", trait1, " -> ", trait2, ")", sep=""), list(trait1=expo_str, trait2=outcome_str)),
         ylab=substitute(paste(trait, " (", -log[10], italic(P[zx]), ")", sep=""), list(trait=expo_str)))
}

# ************************************************** #
#                  Effect size plot                  #
# ************************************************** #
# expo_str, exposure
# outcome_str, outcome
# effect_col, plotting colour
plot_gsmr_effect = function(gsmr_data, expo_str, outcome_str, effect_col=colors()[75]) {
    resbuf = gsmr_snp_effect(gsmr_data, expo_str, outcome_str);
    bxy = resbuf$bxy
    bzx = resbuf$bzx; bzx_se = resbuf$bzx_se;
    bzy = resbuf$bzy; bzy_se = resbuf$bzy_se;
    # plot
    plot_snp_effect(expo_str, outcome_str, bxy, bzx, bzx_se, bzy, bzy_se, effect_col)
}

# ************************************************** #
#                    P-value plot                    #
# ************************************************** #
# expo_str, exposure
# outcome_str, outcome
# effect_col, plotting colour
plot_gsmr_pvalue = function(gsmr_data, expo_str, outcome_str, gwas_thresh=5e-8, truncation=1e-50, effect_col=colors()[75]) {
    resbuf = gsmr_snp_effect(gsmr_data, expo_str, outcome_str);
    bzx_pval = resbuf$bzx_pval; bzy_pval = resbuf$bzy_pval;
    # plot
    plot_snp_pval(expo_str, outcome_str, bzx_pval, bzy_pval, gwas_thresh, truncation, effect_col)
}

# ************************************************** #
#                     bxy distribution plot                         #
# ************************************************** #

# expo_str, exposure
# outcome_str, outcome
# effect_col, plotting colour
plot_bxy_distribution = function(gsmr_data, expo_str, outcome_str, effect_col=colors()[75]) {
    resbuf = gsmr_snp_effect(gsmr_data, expo_str, outcome_str);
    bzx = resbuf$bzx; bzx_pval = resbuf$bzx_pval;
    bzy = resbuf$bzy; 
    bxy = bzy/bzx
    # plot
    plot_snp_bxy(expo_str, outcome_str, bxy, bzx_pval, effect_col)
}
# ===============================================================================
# Author: Futao Zhang, Zhihong Zhu
# Date started: 18/03/2016
# Date last updated: 24/03/2017
# R script to draw regional plot and effect size plot for SMR analysis
#
#   For regional plot, users should specify the parameters such as SMR threshold,
# Heidi threshold, plot window size et al. We also provide parameters smr_thresh_plot
# and probeNEARBY to draw eQTL plots of probes of interest. We set these 2 parameters
# to NULL as default.
#
# AS DEFAULT ONLY THE TARGET PROBE AND THE PROBES THAT PASSED THE PSMR THRESHOLD WOULD BE SHOWN IN THE EQTL LAYER.
# IF YOU WANT TO PLOT SELECTED PROBES IN EQTL LAYER, PLEASE USE PARAMETER smr_thresh_plot OR probeNEARBY
#
#   In order to get your satisfactory plot, it is necessary to include information in
# the plot file as many as possible. so when generating plot files, we recommend to
# use these two parameters: --psmr 1 and --phet 0.
#
# Amendment:
#  1. In eQTL layers, we use maroon color to indicate this probe passed SMR threshold,
#                         navy color to indicate this probe did not pass SMR threshold.
#  2. In GWAS layer, we use maroon color to indicate this probe passed SMR threshold,
#                         navy color to indicate this probe did not pass SMR threshold.
#                 we use solid rhombus to indicate this probe passed HEIDI threshold,
#                        hollow rhombus to indicate this probe did not pass HEIDI threshold.
#  3. shifted label of SMR threshold a little right in case of shading GWAS signals.
#  4. enlarged vetical axis of eQTL layer in order to prevent the label of probe name shading
#       eQTL signals.
#  5. fixed a bug that probe names could overlaps sometimes.
#  example:
#        source("plot_SMR.r")
#        smrdata=ReadSMRData("ILMN_1719097.ILMN_1719097.txt")
#        # target probe passed pSMR
#        SMRLocusPlot(data=smrdata,smr_thresh = 8.4e-6,heidi_thresh = 0.05,plotWindow = 1000)
#        # target probe did not pass pSMR
#        SMRLocusPlot(data=smrdata,smr_thresh = 8.4e-8,heidi_thresh = 0.05,plotWindow = 1000)
#        # use "smr_thresh_plot" to show ILMN_2404135
#        SMRLocusPlot(data=smrdata,smr_thresh = 8.4e-6,heidi_thresh = 0.05,plotWindow = 1000, smr_thresh_plot=1e-1)
#        # use "probeNEARBY" to show ILMN_1724700
#        SMRLocusPlot(data=smrdata,smr_thresh = 8.4e-6,heidi_thresh = 0.05,plotWindow = 1000, probeNEARBY=c("ILMN_1724700","ILMN_does_not_exist"))
#        # use library of the third part to arrange the probe names
#        SMRLocusPlot(data=smrdata,smr_thresh = 8.4e-6,heidi_thresh = 0.05,plotWindow = 1000,anno_selfdef=FALSE) # default anno_selfdef=TRUE
# ===============================================================================
is.installed <- function(mypkg){
    is.element(mypkg, installed.packages()[,1])
}
# check if package "TeachingDemos" is installed
if (!is.installed("TeachingDemos")){
    install.packages("TeachingDemos");
}
library("TeachingDemos")

# parameters for plot
genemove = 0.01; txt=1.1;  cex =1.3; lab=1.1; axis=1; top_cex=1.2;


GeneRowNum = function(GENELIST) {
    BP_THRESH = 0.03; MAX_ROW = 5
    # get the start and end position
    GENELIST = GENELIST[!duplicated(GENELIST$GENE),]
    START1 = as.numeric(GENELIST$GENESTART); END1 = as.numeric(GENELIST$GENEEND)
    STRLENGTH = nchar(as.character(GENELIST$GENE))
    MIDPOINT = (START1 + END1)/2
    START2 = MIDPOINT-STRLENGTH/250; END2 = MIDPOINT+STRLENGTH/250
    START = cbind(START1, START2); END = cbind(END1, END2);
    START = apply(START, 1, min); END = apply(END, 1, max)
    GENELIST = data.frame(GENELIST, START, END)
    GENELIST = GENELIST[order(as.numeric(GENELIST$END)),]
    START = as.numeric(GENELIST$START); END = as.numeric(GENELIST$END)
    # get the row index for each gene
    NBUF = dim(GENELIST)[1]
    ROWINDX = rep(1, NBUF)
    ROWEND = as.numeric(rep(0, MAX_ROW))
    MOVEFLAG = as.numeric(rep(0, NBUF))
    if(NBUF>1) {
        for( k in 2 : NBUF ) {
            ITERFLAG=FALSE
            if(START[k] < END[k-1]) {
                INDXBUF=ROWINDX[k-1]+1
            } else INDXBUF = 1
            if(INDXBUF>MAX_ROW) INDXBUF=1;
            REPTIME=0
            repeat{
                if( ROWEND[INDXBUF] > START[k] ) {
                    ITERFLAG=FALSE
                    INDXBUF=INDXBUF+1
                    if(INDXBUF>MAX_ROW) INDXBUF = 1
                } else {
                    ITERFLAG=TRUE
                }
                if(ITERFLAG) break;
                REPTIME = REPTIME+1
                if(REPTIME==MAX_ROW) break;
            }
            ROWINDX[k]=INDXBUF;
            
            if( (abs(ROWEND[ROWINDX[k]]-START[k]) < BP_THRESH)
            | ((ROWEND[ROWINDX[k]]-START[k])>0) ) {
                MOVEFLAG[k] = 1
                SNBUF = tail(which(ROWINDX[c(1:k)]==ROWINDX[k]), n=2)[1]
                MOVEFLAG[SNBUF] = MOVEFLAG[SNBUF] - 1
            }
            if(ROWEND[ROWINDX[k]]<END[k]) {
                ROWEND[ROWINDX[k]] = END[k]  }
        }
    }
    GENEROW = data.frame(as.character(GENELIST$GENE),
    as.character(GENELIST$ORIENTATION),
    as.numeric(GENELIST$GENESTART),
    as.numeric(GENELIST$GENEEND),
    ROWINDX, MOVEFLAG)
    colnames(GENEROW) = c("GENE", "ORIENTATION", "START", "END", "ROW", "MOVEFLAG")
    return(GENEROW)
}
plot_probe = function(probeinfobuf, k, colplot, x.min, x.max, y.min, y.max,pchbuf,heidi) {
    xcenter = as.numeric(probeinfobuf[k,3])
    pvalbuf = as.numeric(probeinfobuf[k,8])
    strbuf = probeinfobuf[k,1]
    par(new=TRUE)
    if(heidi==TRUE) {
        plot(xcenter, pvalbuf, ylim=c(y.min,y.max),  xlim=c(x.min,x.max),cex.axis=axis,
        xlab="", ylab="", col=colplot, bg=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
    } else {
        plot(xcenter, pvalbuf, ylim=c(y.min,y.max),  xlim=c(x.min,x.max),cex.axis=axis,
        xlab="", ylab="", col=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
    }
}
ReadSMRData = function(plotfile)
{
    SMRData = list();
    key=c("$probe","$SNP","$GWAS","$eQTL");
    skiplines=0;
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    if(keywords[1]!=key[1])
    {
        print("ERROR: plot file is not correct!");
        quit();
    }
    nprobes=as.numeric(keywords[2]);
    SMRData$probeID=keywords[3];
  
    
    skiplines=skiplines+1;
    SMRData$SMR=read.table(plotfile, header=F, nrows=nprobes, skip=skiplines);
    skiplines=skiplines+nprobes;
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    if(keywords[1]!=key[2])
    {
        print("ERROR: plot file is not correct!");
        quit();
    }
    nrs=as.numeric(keywords[2]);
    skiplines=skiplines+1;
    SMRData$SNP=read.table(plotfile, header=F, nrows=nrs, skip=skiplines);
    skiplines=skiplines+nrs;
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    if(keywords[1]!=key[3])
    {
        print("ERROR: plot file is not correct!");
        quit();
    }
    ngwas=as.numeric(keywords[2]);
    skiplines=skiplines+1;
    SMRData$GWAS=read.table(plotfile, header=F, nrows=ngwas, skip=skiplines);
    skiplines=skiplines+ngwas;
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    if(keywords[1]!=key[4])
    {
        print("ERROR: plot file is not correct!");
        quit();
    }
    neqtl=as.numeric(keywords[2]);
    skiplines=skiplines+1;
    
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    prbname=keywords[1];
    neqtlsnp=as.numeric(keywords[2]);
    skiplines=skiplines+1;
    SMRData$eQTL=read.table(plotfile, header=F, nrows=neqtlsnp, skip=skiplines);
    SMRData$eQTL=cbind(prbname,SMRData$eQTL)
    skiplines=skiplines+neqtlsnp;
    if(neqtl>1)
    {
        for(i in 2:neqtl)
        {
            keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
            prbname=keywords[1];
            neqtlsnp=as.numeric(keywords[2]);
            skiplines=skiplines+1;
            raweQTLtmp=read.table(plotfile, header=F, nrows=neqtlsnp, skip=skiplines);
            raweQTLtmp=cbind(prbname,raweQTLtmp);
            SMRData$eQTL=rbind(SMRData$eQTL,raweQTLtmp);
            skiplines=skiplines+neqtlsnp;
        }
    }
    
    keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
    if(length(keywords)>0)
    {
        if(keywords[1]!="$Gene")
        {
            print("ERROR: plot file is not correct!");
            quit();
        }
        ngenes=as.numeric(keywords[2]);
        skiplines=skiplines+1;
        SMRData$Gene=read.table(plotfile, header=F, nrows=ngenes, skip=skiplines);
    }
    return(SMRData)
}
SMRLocusPlot = function(data=SMRData, probeNEARBY=NULL,smr_thresh=NULL, smr_thresh_plot=NULL, heidi_thresh=NULL, plotWindow=NULL,pointsize=20,max_anno_probe=16,anno_selfdef=TRUE)
{
    
    cex_coeff=3/4 * pointsize/15;
    if(length(smr_thresh)==0){
        print("ERROR: please specify the threshold of SMR test!");
        quit();
    }
    if(length(heidi_thresh)==0){
        print("ERROR: please specify the threshold of HEIDI test!");
        quit();
    }
    if(length(plotWindow)==0){
        print("ERROR: please specify the plot window size!");
        quit();
    }
    if(length(which(is.na(data$SMR[,3])))>0)
    {
        print("ERROR: Some probes' physical positon is missing!");
        quit();
    }
    idx=match(data$probeID,data$SMR[,1]);
    if(length(idx)==0){
        print("ERROR: Plot file is not generated correctly, can't find target probe!");
        quit();
    }
    if(length(smr_thresh_plot)==0){
        smr_thresh_plot=smr_thresh;
    }
    cis_start=data$SMR[idx,3]-plotWindow*1000;
    if(cis_start<0) cis_start=0
    cis_end=data$SMR[idx,3]+plotWindow*1000;
    idx=which(data$SMR[,3]>=cis_start & data$SMR[,3]<=cis_end)
    data$SMR=data$SMR[idx,]
    idx=match(data$GWAS[,1],data$SNP[,1])
    tmpsnpbp=data$SNP[idx,3]
    idx=which(tmpsnpbp>=cis_start &tmpsnpbp<=cis_end)
    data$GWAS=data$GWAS[idx,]
    idx=match(data$eQTL[,2],data$SNP[,1])
    tmpsnpbp=data$SNP[idx,3]
    idx=which(tmpsnpbp>=cis_start &tmpsnpbp<=cis_end)
    data$eQTL=data$eQTL[idx,]
    
    if(!is.null(data$Gene))
    {
        idx=which(data$Gene[,2]>=cis_start & data$Gene[,3]<=cis_end )
        data$Gene=data$Gene[idx,]
    }
    
    #start to plot
    smrindx = which(data$SMR[,8] <= smr_thresh_plot)
    #heidiindx = which((data$SMR[,8] <= smr_thresh_plot) & (data$SMR[,9] >= heidi_thresh_plot))
    smrprobes = NULL; heidiprobes = NULL;
    if(length(smrindx)>0) { smrprobes =  as.character(data$SMR[smrindx,1]) }
    #if(length(heidiindx)>0) { heidiprobes = as.character(data$SMR[heidiindx,1]) }
    
    smrindx_bonferr = which(data$SMR[,8] <= smr_thresh)
    heidiindx_strengent = which((data$SMR[,9] >= heidi_thresh))
    smrprobes_red = NA; heidiprobes_solid = NA;
    if(length(smrindx_bonferr)>0) { smrprobes_red =  as.character(data$SMR[smrindx_bonferr,1]) }
    if(length(heidiindx_strengent)>0) { heidiprobes_solid = as.character(data$SMR[heidiindx_strengent,1]) }
    
    if(length(probeNEARBY)>0)
    {
        idx=match(probeNEARBY,data$SMR[,1])
        idxx=which(is.na(idx))
        if(length(idxx)>0)
        {
            for(ii in 1:length(idxx)) {
                print(paste("WARNING: cann't find probe ",probeNEARBY[idxx[ii]], " in plot region.",sep=""))
            }
            probeNEARBY=probeNEARBY[-idxx]
        }
        
    }
    probePLOT=smrprobes #draw the eQTL of all the probes that passed smr_thresh_plot
    probePLOT=unique(c(data$probeID,probePLOT,probeNEARBY)) # draw the target probe anyway
    nprobePLOT = length(probePLOT)
    
	idx=which(is.na(data$GWAS[,2]) | is.na(data$GWAS[,3]))
    if(length(idx)>0) data$GWAS=data$GWAS[-idx,]
	pZY=-log10(pchisq((data$GWAS[,2]/data$GWAS[,3])^2,1,lower.tail=F))
    
    idx=match(data$probeID,data$SMR[,1]);
    if(length(idx)>0){
        chrPLOT = data$SMR[idx,2]
    }else{
        print("ERROR: Plot file is not generated correctly, please report this bug!");
        quit();
    }
    idx=which(is.na(data$SMR[,8]) )
    if(length(idx)>0) {
        probeINFO=data$SMR[-idx,];
    }else{
        probeINFO=data$SMR;
    }
    idx=which(is.na(probeINFO[,5]) | is.na(probeINFO[,6]));
    idx2=which(is.na(probeINFO[,3]));
    if(length(intersect(idx,idx2))>0)
    {
        print("ERROR: Some probes' physical positon is missing!");
        quit();
    }
    probeINFO[idx,5]=probeINFO[idx,3]-7500;
    probeINFO[idx,6]=probeINFO[idx,3]+7500;
    probeINFO[,8]=-log10(probeINFO[,8]);
    probeINFO[,3]=probeINFO[,3]/1e6;
    probeINFO[,5]=probeINFO[,5]/1e6;
    probeINFO[,6]=probeINFO[,6]/1e6;
    pXY=probeINFO[,8];
    yMAX = ceiling(max(c(pZY, pXY), na.rm=T)) + 1;
    if(is.null(data$Gene))
    {
        glist=cbind(probeINFO[,2],probeINFO[,5:6],as.character(probeINFO[,4]),probeINFO[,7]);
    } else {
        glist=data$Gene;
        glist[,2]=glist[,2]/1e6;
        glist[,3]=glist[,3]/1e6;
    }
    colnames(glist)=c("CHR", "GENESTART",  "GENEEND",   "GENE", "ORIENTATION");
    idx=which(is.na(glist[,2]) | is.na(glist[,3]));
    if(length(idx>0)) glist=glist[-idx,];
    generow = GeneRowNum(glist);
    num_row = max(as.numeric(generow$ROW));
    offset_map = ceiling(yMAX);
    offset_probe = yMAX / 2.5;
    num_probe = nprobePLOT
    offset_eqtl = ceiling(yMAX / 2.5) + 0.5;
    dev_axis = 0.1*yMAX;
    if(dev_axis<1.5) dev_axis = 1.5;
    yaxis.min = -offset_map - offset_eqtl*num_probe - dev_axis*(num_probe+1);
    yaxis.max = yMAX + ceiling(offset_probe) + 1;
    # scales of x-axis
    idx=match(data$GWAS[,1],data$SNP[,1]);
    gwasBP = as.numeric(data$SNP[idx,3])/1e6;
    #min.pos = min(gwasBP);
    #max.pos = max(gwasBP);
    min.pos = cis_start/1e6
    max.pos = cis_end/1e6
    start = min(as.numeric(glist[,2]));
    end = max(as.numeric(glist[,3]));
    bp = c(min.pos, max.pos, start, end);
    xmin = min(bp, na.rm=T) - 0.001;  xmax = max(bp, na.rm=T) +0.001;
     xmax=xmax+(xmax-xmin)*0.1 #extend
    ylab = expression(paste("-", log[10], "(", italic(P), " GWAS or SMR)", sep=""));
    xlab = paste("Chromosome", chrPLOT, "Mb");
    # plot GWAS p value
    par(mar=c(5,5,3,2), xpd=TRUE)
    plot(gwasBP, pZY, yaxt="n", bty="n", ylim=c(yaxis.min,yaxis.max),
    ylab="", xlab=xlab, cex.lab=lab, cex.axis=axis,cex=0.6,
    xlim=c(xmin, xmax), pch=20, col="gray68");
    
    # y1 axis
    devbuf1 = yMAX/4
    axis(2, at=seq(0,yMAX,devbuf1), labels=round(seq(0,yMAX,devbuf1),0), las=1, cex.axis=axis);
    mtext(ylab, side=2, line=3, at=(yMAX*2/3), cex=cex_coeff);
    eqtl.lab = expression(paste("-", log[10], "(", italic(P), " eQTL)", sep=""));
    axis.start = 0; axis.down = offset_eqtl + dev_axis;
    for( k in 1 : nprobePLOT ) {
        axis.start = axis.start - axis.down
        eqtlinfobuf = data$eQTL[which(data$eQTL[,1]==probePLOT[k]),]
        if(dim(eqtlinfobuf)[1]==0) next;
        pvalbuf=-log10(pchisq((eqtlinfobuf[,3]/eqtlinfobuf[,4])^2,1,lower.tail=F));
        pvalbuf[which(is.infinite(pvalbuf))]=1e-300;
        if(length(which(smrprobes_red==probePLOT[k]))==0) {
            col_eqtl = "navy"
        } else col_eqtl = "maroon"
        eqtl.min = 0; eqtl.max = ceiling(max(pvalbuf))
        eqtl.max =ceiling(eqtl.max *1.25) #extend
        pvalbuf = pvalbuf/eqtl.max * offset_eqtl + axis.start
        idx=match(eqtlinfobuf[,2],data$SNP[,1]);
        eqtlbp = as.numeric(data$SNP[idx,3])/1e6;
        probegene = unique(as.character(data$SMR[which(data$SMR[,1]==probePLOT[k]),4]))
        par(new=TRUE)
        pchbuf = 4;
        #if(k%%2==0) pchbuf = 20;
        plot(eqtlbp, pvalbuf, yaxt="n", bty="n", ylim=c(yaxis.min,yaxis.max), xaxt="n",
        ylab="", xlab="", cex=0.8, pch=pchbuf, col=col_eqtl, xlim=c(xmin, xmax))
        # annotate the eQTLs
        text(xmin, axis.start+offset_eqtl-dev_axis/2 , label=substitute(paste(probeid, " (",italic(geneid), ")", sep=""),list(probeid=probePLOT[k], geneid=probegene)),col="black", cex=1, adj=0)
        # axis
        devbuf1 = offset_eqtl/3; devbuf2 = eqtl.max/3
        axis(2, at=seq(axis.start,(axis.start+offset_eqtl),devbuf1),
        labels=round(seq(0,eqtl.max,devbuf2),0),
        las=1, cex.axis=axis)
        # add separator line
        segments(xmin, axis.start+offset_eqtl+dev_axis/2, xmax, axis.start+offset_eqtl+dev_axis/2,
        col="dim grey", lty="24", lwd=1)
    }
    #ypos = (axis.start - dev_axis)/2
    ypos = (axis.start - dev_axis)*2/3
    mtext(eqtl.lab, side=2, at=ypos, line=3, cex=cex_coeff)
    
    # plot p value of bTG
    # all the probes
    num_gene = dim(generow)[1]
    dist = offset_map/num_row
    for( k in 1 : num_row ) {
        generowbuf = generow[which(as.numeric(generow[,5])==k),]
        xstart = as.numeric(generowbuf[,3])
        xend = as.numeric(generowbuf[,4])
        snbuf = which(xend-xstart< 1e-3)
        if(length(snbuf)>0) {
            xstart[snbuf] = xstart[snbuf] - 0.0025
            xend[snbuf] = xend[snbuf] + 0.0025
        }
        xcenter = (xstart+xend)/2
        xcenter = spread.labs(xcenter, mindiff=0.01, maxiter=1000, min = xmin, max = xmax)
        num_genebuf = dim(generowbuf)[1]
        for( l in 1 : num_genebuf ) {
            ofs=0.3
            if(l%%2==0) ofs=-0.8
            m = num_row - k
            ypos = m*dist + yaxis.min
            code = 1
            if(generowbuf[l,2]=="+") code = 2;
            arrows(xstart[l], ypos, xend[l], ypos, code=code, length=0.07, ylim=c(yaxis.min,yaxis.max),
            col=colors()[75], lwd=1)
            movebuf = as.numeric(generowbuf[l,6])*genemove
            text(xcenter[l]+movebuf, ypos,label=substitute(italic(genename), list(genename=as.character(generowbuf[l,1]))), pos=3, offset=ofs, col="black", cex=0.9)
        }
    }
    
    # plot the probes
    probeINFO=probeINFO[order(probeINFO[,8],decreasing = TRUE),];
    nprobeINFO=dim(probeINFO)[1];
    if(nprobeINFO>max_anno_probe){
        probeINFO=probeINFO[c(1:max_anno_probe),]
        nprobeINFO=dim(probeINFO)[1];
    }
    if(anno_selfdef) probeINFO=probeINFO[order(probeINFO[2],probeINFO[3]),] ####20170217
    xcenter = as.numeric(probeINFO[,3])
    xcbuf = xcenter
    ####20170217####
    if(anno_selfdef)
    {
        reginlength=(xmax-(xmax-xmin)*0.15)-xmin
        leftspot=xmin+reginlength/20
        rightspot=(xmax-(xmax-xmin)*0.15)-reginlength/20
        itvl=(rightspot-leftspot)/dim(probeINFO)[1]
        if(dim(probeINFO)[1]==1) {
            xcenter=as.numeric(probeINFO[,3])
        } else {
            xcenter=leftspot+itvl/2
            for( k in 2:dim(probeINFO)[1]) xcenter=c(xcenter,leftspot+k*itvl)
        }
        
    } else {
        xcenter = spread.labs(xcenter[1:nprobeINFO], mindiff=0.08, maxiter=1000, min = xmin, max = xmax-1)
    }
    # adjust the line position
    
    adjflag = rep(0, nprobeINFO)
    if(nprobeINFO>1) {
        dbuf = c(0, xcbuf[1:(nprobeINFO-1)])
        mflag = as.numeric(abs(xcbuf[1:(nprobeINFO)] - dbuf) < 0.01)
        adjflag = as.numeric( mflag | c(mflag[2:nprobeINFO],0) )
    }
    
    for( k in 1 : nprobeINFO)  {
         hitflag=FALSE
        if(length(which(heidiprobes_solid==probeINFO[k,1]))>0 & length(which(smrprobes_red==probeINFO[k,1]))>0) {
             hitflag=TRUE
             colplot = "maroon"; colfont=2; pchbuf=23;
        } else if(length(which(smrprobes_red==probeINFO[k,1]))>0) {
            colplot = "maroon"; colfont=2; pchbuf=5
            #} else if (length(which(heidiprobes_solid==probeINFO[k,1]))>0) {
            #hitflag=TRUE
            # colplot = "navy"; colfont=1; pchbuf=23
        } else {
            colplot = "navy"; colfont=1; pchbuf=5
        }
        if( as.numeric(probeINFO[k,8]) < 0 ) {
            colplot = "black"; colfont=1;
        }
        # plot p value of bxy
        plot_probe(probeINFO, k, colplot, xmin, xmax, yaxis.min, yaxis.max,pchbuf,hitflag)
        # annotate the probes
        if(k<=max_anno_probe)
        {
            ypos = 1.02*yMAX
            strbuf =
            text(xcenter[k], ypos,
            labels=substitute(paste(probeid, " (", italic(genename), ")", sep=""),
            list(probeid=as.character(probeINFO[k,1]),
            genename=as.character(probeINFO[k,4]))),
            ylim=c(yaxis.min, yaxis.max),
            srt=30, col=colplot, font=colfont, cex=1, adj=0)
            # plot the lines
            # 1st step
            xstart = xcbuf[k]
            ystart = as.numeric(probeINFO[k,8]); yend = yMAX*(1-1/20);
            if( nprobeINFO > 1 ) {
                if(adjflag[k]==1) {
                    xstart = (xcbuf[k] + xcenter[k])/2
                    segments(xcbuf[k], ystart, xstart, ystart, col=colplot, lwd=axis, lty=3)
                }
            }
            segments(xstart, ystart, xstart, yend, col=colplot, lwd=axis, lty=3)
            # 2nd step
            xend = xcenter[k]; ystart = yMAX*(1-1/20); yend = yMAX*1.01;
            segments(xstart, ystart, xend, yend, col=colplot, lwd=axis, lty=3)
        }
    }
    # plot the threshold
    # SMR threshold
    ybuf = -log10(as.numeric(smr_thresh)); dev_anno = yMAX/9;
    strbuf = paste("pSMR = ",smr_thresh, sep="")
    segments(xmin, ybuf, xmax, ybuf, col="maroon", lty=2, lwd=1);
    text(xmax, ybuf+dev_anno, labels=strbuf, adj=1, col="maroon", cex=axis,font=3);
}


SMREffectPlot = function(data=SMRData, trait_name="",cisWindow=2000, transWindow=5000, pointsize=20)
{    
    # parameters for plot
    pch_top = 24; pch_cis = 21; pch_trans = 22
    col_top = "red"; col_cis = "Navy"; col_trans = "green"
    cex_coeff=3/4 * pointsize/15;
    
    # Extract the probe for plot
    snbuf = which(as.character(data$eQTL[,1])==data$probeID)
    if(length(snbuf)==0) {
        print(paste("ERROR: no eQTL infomation found for probe",data$probeID,"!",sep=""));
        quit();
    }
    plotData = data$eQTL[snbuf,]
    idx=which(is.na(plotData[,5]))
    if(length(idx)>0) plotData=plotData[-idx,]
    
    # SNPs in common
    snpbuf = Reduce(intersect, list(as.character(plotData[,2]), data$GWAS[,1]))
    plotData = plotData[match(snpbuf, as.character(plotData[,2])),]
    plotGWAS = data$GWAS[match(snpbuf, as.character(data$GWAS[,1])),]
    # Effect size
    snplist = as.character(plotData[,2])
    bZX = as.numeric(as.character(plotData[,3]));
    seZX = as.numeric(as.character(plotData[,4]));
    snpCorr=as.numeric(as.character(plotData[,5]));
    bZY = as.numeric(as.character(plotGWAS[,2]));
    seZY = as.numeric(as.character(plotGWAS[,3]));
    # Limit
    xmin =  min(bZX - seZX, na.rm=T)
    xmax =  max(bZX + seZX, na.rm=T)
    ymin =  min(bZY - seZY, na.rm=T)
    ymax =  max(bZY + seZY, na.rm=T)
    
    if(xmin>0) xmin = -xmax/2
    if(xmax<0) xmax = -xmin/2
    if(ymin>0) ymin = -ymax/2
    if(ymax<0) ymax = -ymin/2
    
    # Plots
    #par(mar=c(5,6.5,5,1), xpd=FALSE)
    #layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(4.5,1), heights=c(1,1))

    # Start to plot
    nsnps = dim(plotData)[1]
    # Split data to cis- and trans-
    idx=which(data$SMR[,1]==data$probeID);
    if(length(idx)!=1)
    {
        print("ERROR: plot file is not correct!");
        quit();
    }
    if(is.na(data$SMR[idx,8]))
    {
        print(paste("ERROR: no SMR reslult for probe",data$probeID,"!",sep=""));
        quit();
    }
    probeChr = as.numeric(as.character(data$SMR[idx,2]))
    probeBP = as.numeric(as.character(data$SMR[idx,3]))
    GeneID =as.character(data$SMR[idx,4])
    idx=match(snplist,data$SNP[,1]);
    snpChr = as.numeric(as.character(data$SNP[idx,2]))
    snpBP = as.numeric(as.character(data$SNP[idx,3]))
    cisIndx = which(probeChr==snpChr & abs(snpBP-probeBP)<cisWindow*1000);
    ncis = length(cisIndx)
    transIndx = which(probeChr!=snpChr | (probeChr==snpChr & abs(snpBP-probeBP)>transWindow*1000));
    ntrans = length(transIndx)
    # Plot the cis-eQTL
    snplist_tmp = snplist[cisIndx]
    maxsnpCorr = snpCorr[cisIndx]
    bZX_tmp = bZX[cisIndx]; seZX_tmp = seZX[cisIndx]; zZX_tmp = bZX_tmp/seZX_tmp;
    bZY_tmp = bZY[cisIndx]; seZY_tmp = seZY[cisIndx]; zZY_tmp = bZY_tmp/seZY_tmp;
    maxid = which.max(zZX_tmp^2)
    maxsnp = snplist[maxid]
    maxsnpCorr = maxsnpCorr^2;
    for( k in 1 : ncis ) {
        # effect sizes
        colbuf = rgb(0, 0, 128/255, maxsnpCorr[k])        
        colcir = rgb(0, 0, 1-maxsnpCorr[k]);
        cex = 1
        plot(bZX_tmp[k], bZY_tmp[k], pch=pch_cis, col=colcir, bg=colbuf,
        bty="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax),
        cex=cex, xlab="", ylab="", xaxt="n", yaxt="n")
        par(new=TRUE)
        plot(bZX_tmp[k], bZY_tmp[k], pch=20, col=colcir, bg=colbuf,
        bty="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax),
        cex=0.1, xlab="", ylab="", xaxt="n", yaxt="n")
        par(new=TRUE)
    }
    
    # standard error
    # cis eQTL
    for( k in 1 : ncis ) {
        colcir = rgb(0, 0, 1-maxsnpCorr[k]);
        segments(bZX_tmp[k]-seZX_tmp[k], bZY_tmp[k], bZX_tmp[k]+seZX_tmp[k], bZY_tmp[k],
        col=colcir, lwd=0.5+maxsnpCorr[k])
        segments(bZX_tmp[k], bZY_tmp[k]-seZY_tmp[k], bZX_tmp[k], bZY_tmp[k]+seZY_tmp[k],
        col=colcir, lwd=0.5+maxsnpCorr[k])
    }
    
    # line
    colline = rgb(244/255,164/255,96/255,1)
    bXY = bZY_tmp[maxid]/bZX_tmp[maxid]
    abline(0, bXY, col=colline, lwd=2, lty=2)
    
    # plot effect size of the top SNP
    colbuf = "white"
    colcir = col_top
    cex=2.3
    par(new=TRUE)
    plot(bZX_tmp[maxid], bZY_tmp[maxid], pch=pch_top, col=colcir, bg=colbuf,
    bty="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax),
    cex=cex, xlab="", ylab="", xaxt="n", yaxt="n")
    colbuf = col_top
    colcir = col_top
    cex = 1
    par(new=TRUE)
    plot(bZX_tmp[maxid], bZY_tmp[maxid], pch=pch_top, col=colcir, bg=colbuf,
    bty="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax),
    cex=cex, xlab="", ylab="", xaxt="n", yaxt="n")
    
    # se of the top SNP
    colcir = rgb(1,0,0)
    segments(bZX_tmp[maxid]-seZX_tmp[maxid], bZY_tmp[maxid],
    bZX_tmp[maxid]+seZX_tmp[maxid], bZY_tmp[maxid],
    col=colcir, lwd=1.5)
    segments(bZX_tmp[maxid], bZY_tmp[maxid]-seZY_tmp[maxid],
    bZX_tmp[maxid], bZY_tmp[maxid]+seZY_tmp[maxid],
    col=colcir, lwd=1.5)
    
    # Plot the trans-eQTLs
    if(ntrans>0) {
        snplist_tmp = snplist[transIndx]
        bZX_tmp = bZX[transIndx]; seZX_tmp = seZX[transIndx]; zZX_tmp = bZX_tmp/seZX_tmp;
        bZY_tmp = bZY[transIndx]; seZY_tmp = seZY[transIndx]; zZY_tmp = bZY_tmp/seZY_tmp;
        par(new=TRUE)
        for( k in 1 : ntrans ) {
            # effect sizes
            colbuf = col_trans;
            colcir = col_trans;
            cex = 1
            plot(bZX_tmp[k], bZY_tmp[k], pch=pch_cis, col=colcir, bg=colbuf,
            bty="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax),
            cex=cex, xlab="", ylab="", xaxt="n", yaxt="n")
            par(new=TRUE)
            plot(bZX_tmp[k], bZY_tmp[k], pch=20, col=colcir, bg=colbuf,
            bty="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax),
            cex=0.1, xlab="", ylab="", xaxt="n", yaxt="n")
            par(new=TRUE)
        }
        # standard error
        # trans-eQTL
        for( k in 1 : ntrans ) {
            colcir = col_trans
            segments(bZX_tmp[k]-seZX_tmp[k], bZY_tmp[k], bZX_tmp[k]+seZX_tmp[k], bZY_tmp[k],
            col=colcir, lwd=1)
            segments(bZX_tmp[k], bZY_tmp[k]-seZY_tmp[k], bZX_tmp[k], bZY_tmp[k]+seZY_tmp[k],
            col=colcir, lwd=1)
        }
    }
    
    # plot the axis
    # x axis
    devbuf = (xmax - xmin)/5
    if(xmax!=0 & xmin!=0) {
        numbuf = min(abs(xmin), abs(xmax))
        if( devbuf > numbuf ) devbuf = numbuf
    }
    numbuf = as.numeric()
    if( xmin < 0 ) numbuf = c(numbuf, -seq(0, abs(xmin), devbuf))
    if( xmax > 0 ) numbuf = c(numbuf, seq(0, xmax, devbuf))
    axis(1, at=numbuf, labels=round(numbuf,2), las=1, cex.axis=axis)
    xmid = (xmax+xmin)/2
    mtext("eQTL effect sizes", side=1, at=xmid, line=3, cex=cex_coeff)
    # y axis
    devbuf = (ymax - ymin)/5
    if(ymax!=0 & ymin!=0) {
        numbuf = min(abs(ymin), abs(ymax))
        if( devbuf > numbuf ) devbuf = numbuf
    }
    numbuf = as.numeric()
    if( ymin < 0 ) numbuf = c(numbuf, -seq(0, abs(ymin), devbuf))
    if( ymax > 0 ) numbuf = c(numbuf, seq(0,ymax,devbuf))
    axis(2, at=numbuf, labels=round(numbuf,3), las=1, cex.axis=axis)
    ymid = (ymax + ymin)/2
    mtext("GWAS effect sizes", side=2, at=ymid, line=4.5, cex=cex_coeff)
    
    mainstr1 = trait_name
    mainstr2 = substitute(paste(probeid, " (", italic(gene), ")", sep=""),
    list(probeid=as.character(data$probeID),
    gene=as.character(GeneID)))
    mtext(mainstr1, side=3, at=xmin, adj=0, line=2.5, cex=cex_coeff)
    mtext(mainstr2, side=3, at=xmin, adj=0, line=0.5, cex=cex_coeff)
    # Plot legend
    lstr = c("top cis-eQTL", "cis-eQTL")
    col_led = c(col_top, col_cis); pch_led = c(pch_top, pch_cis)
    if(ntrans>0) {
        lstr=c(lstr, "trans-eQTL"); col_led = c(col_led, col_trans)
        pch_led = c(pch_led, pch_trans)
    }
    
    if(bXY>0) {
        legend("topleft", lstr, bty="n", border="white", pch=pch_led, col=col_led, pt.bg=col_led, cex=axis)
    } else {
        legend("topright", lstr, bty="n", border="white", pch=pch_led, col=col_led, pt.bg=col_led, cex=axis)
    }

   # add the scale bar
   par(mar=c(5,1,5,4.5))
   pal=colorRampPalette(c(rgb(0, 0, 1),  rgb(0, 0, 0)))
   breaks <- seq(min(snpCorr^2), max(snpCorr^2),length.out=2000)
   #image.scale(snpCorr^2, col=pal(length(breaks)-1), breaks=breaks, horiz=FALSE, xlab="", yaxt="n")
   dvd = (max(snpCorr^2) - min(snpCorr^2))/5
   pos = seq(min(snpCorr^2), max(snpCorr^2), dvd)
   axis(4, at=pos, label=sprintf("%.2f", pos), las=2)
   mtext(expression(italic(r)^2), side=4, line=3.5, las=2 )
}

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
)�      r�0��b```f`afd`f2XCC�t-X��FN�$�@�0� �'pA   shinyUI(fluidPage(
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
))##blup
ph <- readxl::read_excel("~/Documents/databaseMB/database/Phenotype/Generations and Phneotype.xlsx")
ph2 <- ph[!is.na(ph$Generation...6),]

dim(ph2)
View(ph2)

ph2 <- ph2[,!(colnames(ph2) %in% "19_JM_SG")]

ts <- c("IL","SG","wilt","TMV","CMV", "LN","LL","LW","PH") 

blup <- data.frame("ID"= unique(ph2$`20_GeID`))
require(lme4)


for(i in 1:length(ts)){
  
  t_now <- ts[i]
  mg1_now <- ph2[,c("20_GeID",colnames(ph2)[grepl(t_now,colnames(ph2))])]
  #mg1_now <- mg1_now[,-2]
  #cat(i,"\n")
  colnames(mg1_now)[1] <- "ID"
  #mg1_now <- mg1_now[-sample(x = mg1_now$ID,size = nrow()-50,]
  
  inData_flatten <- tidyr::gather(mg1_now,key = "Envs",value = "FT_times",colnames(mg1_now)[2:ncol(mg1_now)],na.rm = F)
  inData_flatten$FT_times <- as.numeric(inData_flatten$FT_times)
  inData_flatten <- inData_flatten[!is.na(inData_flatten$FT_times),]
  
  if(ncol(mg1_now) ==2){
    blup2 <- mg1_now
    colnames(blup2) <- c("20_GeID",t_now)
  }else{
    lmer1 <- lmer(FT_times ~ Envs + (1|ID)  ,data = inData_flatten)
    #str(lmer1)
    #lmer1@u
    blup2 <- data.frame("ID"= rownames(lmer1@pp$Ut), lmer1@u)
    colnames(blup2) <- c("ID",t_now)
  }
  
  blup <- merge(blup,blup2,by = "ID")
  cat(i,"\n")
}


##------------------------------------------------------------------
# Kneral weight blup and cv and lp
##------------------------------------------------------------------

ph_KW <- c()
kw_n <- c()
for( i in 1:5){
  t1 <- paste0("EW_",loc[i])
  t2 <- paste0("CW_",loc[i])
  kw_n <- c(kw_n,paste0("KW_",loc[i]))
  ph_now <- ph_all[,t1] - ph_all[,t2] 
  ph_KW <- cbind(ph_KW,ph_now)
}
colnames(ph_KW) <- kw_n
inData <- data.frame(ph_KW) %>% mutate("id"=idnames(data_cubic))
kw_flatten <- tidyr::gather(inData,key = "Envs",value = "FT_times",colnames(inData)[1:(ncol(inData)-1)],na.rm = F)
kw_flatten <- kw_flatten[!is.na(kw_flatten$FT_times),]
##------------------------------------------------------------------
# BLUP
##------------------------------------------------------------------
require(hglm)
hg <- hglm2( FT_times ~ 1 + (1|id) + (1|Envs),data = kw_flatten)
raf <- hg$ranef
mgs <- gsub(".*:(.*)","\\1",names(raf))
names(raf) <- mgs
KW <- raf[idnames(data_cubic)]
##------------------------------------------------------------------
# var
##------------------------------------------------------------------
out_kw  <- remove_e_effects(inData)

KW_Var <- -log10(out_kw$v)

##------------------------------------------------------------------
# CV LP 
##------------------------------------------------------------------
require(dplyr)
require(goeveg)
library("factoextra")

#phe_now_adj <- data.frame(ph_KW)

get_lp <- function(mat){
  mat[is.na(mat)] <- mean(as.numeric(unlist(mat)),na.rm = T)
  res.pca <- prcomp(x = mat,center = T, scale = T)
  res.ind <- get_pca_ind(res.pca)
  if(identical(names(res.ind$coord[,2]),idnames(data_cubic)))
    return(res.ind$coord[,2])
}

KW_CV <- apply(data.frame(ph_KW),1,FUN = cv,na.rm = FALSE)
#CV_up <- apply(out_kw$update_mat[,2:6],1,FUN = cv,na.rm = FALSE)
ph_KW <- data.frame(ph_KW)
rownames(ph_KW) <- idnames(data_cubic)
KW_LP <- get_lp(mat=ph_KW)

hy <- polygenic(ph_out_all[,2],data = data_cubic,kinship.matrix = ibs_cubic)
hy$esth2

##------------------------------------------------------------------
# diff
##------------------------------------------------------------------
pheno_name <- c()
pheno <- c()

t_now <- "KW"
ts_now <- paste0(t_now,"_",loc)

# Dot and diff
# dot
all_ts_now <- combn(ts_now,2)
for(j in 1:ncol(all_ts_now)){
  t1_now <-  all_ts_now[1,j]
  t2_now <-  all_ts_now[2,j]
  #dot_now <- scale(ph_all[,t1_now])*scale(ph_all[,t2_now])
  #dot_name_now <- paste0("dot_",t1_now,"_",t2_now)
  diff_now <- ph_KW[,t1_now] - ph_KW[,t2_now]
  diff_name_now <- paste0("diff_",t1_now,"_",t2_now)
  pheno <- cbind(pheno,diff_now)
  pheno_name <- c(pheno_name,diff_name_now)
}
colnames(pheno) <- pheno_name
ph_out_all <- cbind(ph_KW,pheno,KW,KW_Var,KW_CV,KW_LP)

write.table(ph_out_all,file = "./data/KW_traits.txt",row.names = F,col.names = T,quote = F,sep = "\t")


#-------------------------------------------------------------------------------
# Demo for estimating the plasticity measurement
# Date:200910
# Version:v1
#-------------------------------------------------------------------------------

# 1.Demo data: A tab delimtated txt file with 1404 row and 6 columns
#              Each row is a individual and the first 5 columns are DTT 
#              from different sites, the last column is the individual id 
pheno_dtt <- read.table(file = "./200910_demo_plasticity_measurment.txt",sep = "\t",stringsAsFactors = F)

#-------------------------------------------------------------------------------
# Load required packages
#-------------------------------------------------------------------------------
require(FW)
require(tidyr)
require(dplyr)
require(goeveg)
require(factoextra)

#-------------------------------------------------------------------------------
# Define a few functions
#-------------------------------------------------------------------------------
remove_e_effects <- function(inData){
  # remove E 
  inData_flatten <- tidyr::gather(inData,key = "Envs",value = "FT_times",colnames(inData)[1:(ncol(inData)-1)],na.rm = F)
   y <- inData_flatten$FT_times
  x1 <- as.factor(inData_flatten$id)
  x2 <- as.factor(inData_flatten$Envs)
  bxp1 <- boxplot(y ~ x2,plot=F)
  m <- bxp1$stats[3,]#get the mean value by env
  x3 <- factor(x2,levels = bxp1$names[order(m,decreasing = F)]) #c("FT16","FT10","MHI","MHP","MLI","MLP","THI","THP","TLI","TLP"))
  lm3 <- lm(y~x3)
  slm3 <- summary(lm3)
  beta <- slm3$coefficients[,1]#effect size of env
  names(beta) <- bxp1$names[order(m,decreasing = F)] #c("FT16","FT10","MHI","MHP","MLI","MLP","THI","THP","TLI","TLP")
  beta[1] <- 0#intercept
  y_update <- inData_flatten$FT_times - beta[as.character(inData_flatten$Envs)]#remove the ebv
  bxp2 <- boxplot(y_update ~ x2)
  inData_flatten$FT_times <- y_update
  inData_update <- tidyr::spread(data =inData_flatten,key = "Envs",value = "FT_times")
  rownames(inData_update) <- inData_update$id
  inData_update <- inData_update[inData$id,]
  v <- apply(inData_update[,2:6], 1,FUN = function(x) var(x,na.rm = T))
  return(list("v"=v,"flatten"=inData_flatten,"update_mat"=inData_update,"update"=data.frame("id"=inData_flatten$id,"Envs"=inData_flatten$Envs,"FT_times"=y_update),"beta"=beta))
}

get_lp <- function(mat){
  # pc analysis
  mat[is.na(mat)] <- mean(as.numeric(unlist(mat)),na.rm = T)
  res.pca <- prcomp(x = mat,center = T, scale = T)
  res.ind <- get_pca_ind(res.pca)
  if(identical(names(res.ind$coord[,2]),idnames(data_cubic)))
    return(res.ind$coord[,2])
}


#-------------------------------------------------------------------------------
# measuring linear plasticity
#-------------------------------------------------------------------------------
inData_flatten <- tidyr::gather(pheno_dtt,key = "Envs",value = "GR",colnames(pheno_dtt)[1:(ncol(pheno_dtt)-1)],na.rm = F)
lm1 = FW(y=inData_flatten$GR,VAR=inData_flatten$id,ENV=inData_flatten$Envs)
lp <- lm1$b # This is the linear plasticity if I understand it correctly 

#-------------------------------------------------------------------------------
# measuring across environment variance after removing the main environmental effect
#-------------------------------------------------------------------------------
out <- remove_e_effects(inData=pheno_dtt)
Var <- -log10(out_kw$v) # this is the variance after remove E

#-------------------------------------------------------------------------------
# measuring cv and pc
#-------------------------------------------------------------------------------

CV <- apply(data.frame(pheno_dtt[,1:5]),1,FUN = cv,na.rm = FALSE) # this is the cv
rownames(pheno_dtt) <-pheno_dtt$id
PC <- get_lp(mat=pheno_dtt[,1:5])# this is the second Pc
#-------------------------------------------------------------------------------
# measuring pairwise diff
#-------------------------------------------------------------------------------
pheno_name <- c()
pheno <- c()
ts_now <- colnames(pheno_dtt)[1:5]
all_ts_now <- combn(ts_now,2)
for(j in 1:ncol(all_ts_now)){
  t1_now <-  all_ts_now[1,j]
  t2_now <-  all_ts_now[2,j]
  diff_now <- pheno_dtt[,t1_now] - pheno_dtt[,t2_now]
  diff_name_now <- paste0("diff_",t1_now,"_",t2_now)
  pheno <- cbind(pheno,diff_now)
  pheno_name <- c(pheno_name,diff_name_now)
}
colnames(pheno) <- pheno_name
pheno # this is the paireise differenceimport sys
file1=sys.argv[1]
file2=sys.argv[2]
f_in=open(file1)
f_out=open(file2,"w")
a=f_in.readline()
while a:
    if a[0]=="#" and a[1]!="#":
        a=a.replace("|","_")
        a=a.split("\t")
        for i in range(len(a)):
            if i>=9:
                if "_" not in a[i]:
                    a[i]=a[i]+"_"+a[i]
        a='\t'.join(a)
        #a=a+"\n"
    
    if a[0]!= "#":
        if "scaffold" in a:
            a=f_in.readline()
            continue
        if "*" in a:
            a=f_in.readline()
            continue
        
        a=a.replace("Chr","")
        a=a.replace("chr","")
        a=a.replace("Chromosome","")
        a=a.replace("chromosome","")
        a=a.split("\t")
        if "rs" not in a[2] and ":" not in a[2]:
            a[2]=a[0]+":"+a[1]
        a='\t'.join(a)  
        #a=a+"\n"
        
    f_out.writelines(a)
    a=f_in.readline()
f_in.close()
f_out.close()