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
   IHDR  ,  ,   ^u�   sRGB ���   gAMA  ���a   	pHYs  2�  2�(dZ�  ��IDATx^��	��eu���2���" "�� �8
���DTA�$�� �ePQ�I� ����b�5Q4&&�t��M�����t:�ｫKR�|��[poթ:�^�Y����N�󟡬���1�fV��^�p|X.��;�=�������O���/�߇��5t         �![�l��=�)[�l��M�F)[�l��]ʆ)[�l��m��鸐�S�Oۇl���2�c���ö�pR�$�
~�%T�#         `.�6*�l����v*�l���ʶ*+c�1��L�+��¾!����s���_C�C        ���U6W�^e��-V6Y�f�K��c�fV;������`�u�~x        ��d���V�[�pe˕M�1�3g��a�pZ�-�$T?        �=�ve㕭W6_�~c�1#7����ῇ��2��ʫtk��N���[t[m��n�����|���Ao>�{�q���N��;y�5��>tCw��>՝y�g�s�������u\����]y�w�k��Aw����x�O�[������        �W�F�e{�R�H�$e���R�J�,e��S�L�4e�t��O�h��y�r��M4P�BeU�R3 �l��	�6,1c�1fF��������oC�k�Yq���u�߸�f�]�=�;�;��E�	g\ѝ~�-�y�~������n��g��         �˲��F�k�h���ʆ*[�l����ƪj���lŲ�v,2c�1f����5���B���q[~�����n0�W�w���N���Cg_vgw�g�_�        ��ek��U�W�`e��MV�Y�hU���-Y6eٖec���1��D�b�'\���m�����y��Ï���>~�+�O         ��l����v+�l����	��,۳lвE3�c����oB��dZVY����|��v���_zGw�=?)        0��ʶ+�l������j��lѲI�6-5c�1c6O����OC��b�-���?��=����.��7�n         ��l���&,۰lĪvl��U�f-۵l،1����*�'|1T?��j����˞�a�x����v�~���-         �lĲ�f,۱lȪ�l��e˦-�6c�1sx^>~��%��n�����;O>���_)         0>�%˦,۲l̪�l�q��-�7c�1s`�g��C������̵�]_����/﮹���         �V�fٜe{�Zզ-�l�N��c��yAX�8T��i�w�}���_uo��         �T�h٤e�V5kK�B�q��c����8�7|'TԽ�|ֺ��=�;�̏w7���         <Q٨e���Z�kUӶ���f.�9c�1�p�G���������~���k�X�P         �e-�lٲi�Z�%�]�t��c�YJ�C�$��P}�>�Ͷڱ;�]gvW����         fK�mG����[��=�l겭���c��U��ᛡ��}T�n�}w�1�w�}��<         ���o��D��\��=�l����c�c�K�U�C��Z��[w�����%7}��         ��"[�l�6�l�����]6x��c��2���������J+w�ujw���         0�]�ɯN�r��U-]�l���F�c�v6
������V�v���n�%�-?�        `��v.�l���Qd���^6{�36�s�)T��M6߮{���n��g�0         ��l鲩˶�j�z���3Ƙy;��B�!�x�*O�������+�)?p        `�ec��]6wU��(��˦�c��~����[��{��ݭ_�y��
         ��.ۻl�6�Qdۗ��1���Y3��.Tr�_����o�>z��)         �d���&/ۼ��+d��_��3�V8/�&Tj��̵�C�:����?+?4        ��'ۼl��ի�B�� fh�1#7��s¿��Cl�oٽ����H         `��f/۽��+d�M`���3�j�@������͎�v��{}�a         ,[��e�W5~�l��V�cf|V��?��Cj�+�<�;����?         `feӗm_����f0�Ac�Y��;����P}(=b����=����[�]~�         �+�l����Z�!�fC�-�1�,�9>�2TB^��Gw�����         -��e�W5��l	�)4Ƙ�6���Շ΀����۾[~�         �-�l�F��ma6����%�P}�x��o�.����         0�d�m`��5���c�x��Շʀ]�>������+         `n�F0[��!,d{��1���)�B�A����c�=X~8         �K6��VM�l�E4Ƙf?Շ�#vx�݇����0         �l�%��!�$f�h�1�w�����o���䳮)?|         ��Ma��Us8$�l�1c8+�C����'���c�(?l         ��a��U�8$��l�1c2��
��#�������I�         ��5��j�d���1f�Z���x�6;��}����*          �l�}��e�8���M�1f������:�oܝ��+�         �%�-b6�U�8E6��6c��l�՛}��\w�;�_~h          <�&f�X��Sd㘭�1f��{¿��>a�]��.��k�         ���b��U�8E���<c�����7���<u�|��p          X��Y�v�j���1Hc̈�B�&~�+���������@          X�]|Ş�m�l!�1#8��B�Ɲ��ߝ����C          `&��WO4�U�8Ń!�Hc̈���?��;����)��          3)��l��q�l#��4����՛t�n՝q�m�         `6e㘭c�@N��d6�Ƙ�=��B�ƜpЛO,��          �$�Ǫ��"��l'�134g���8�9�m؝~�-�         `e��d�FNqv0�,��<<�7����>�����P��         FY6��BV��_�Tc��3�&To�n��W�=����         0�d�md�LNʦ2�Jc�R��B�f��^�]z�7�7,         �\�md6�U;9E6�Ƙ'0���6�ෞR�I         �l%��r�l-��4�Ls^�!To�n��mҝy���         `>�f2�ɪ�����>����P��&��u�w�~���         `>�v2ʪ����`���'C�����.߄          � [ʪ��"[�l2�1C�u�~��8�Zk�7��̫7         �8ɦ2�ʪ���Mf��Ƙ�9,��P�a��vzUwݝ?*�p          �(��l,��rR���h3�sn��$^��7          ���;�ز��"[Mc�v>�7F��ݻ]\��          �w�\f{Y5����4f���P�!��m�y����+�P          �����j3'e���1�~�?��{٫��>��_�o$          ]6��bV��?�r3og��_B��~�{�7          K.�̪՜�-���y7G��E?��~�|�          0}�fV���v3ofQ�^��*OY�[t�M�         ��/�l5��sR6�����,T/�n��7�>v݃�         �'.[�l6��sR�����Y1��v��v/�����o          ��l6�ݬ��I�|f�i̜����B���^���7          �N6�U�9)��l@��yV�N�^��~�[�          X����Iـfj�����C����{��>          3'�Ϊ���-h6�ƌ�l�2T/��3?^��         �y'}����Mh��ƌ�l~��+t�>rc�B         `�d㙭gՀ�lC�	�����!4/�UV}jw�E��/p          f_���|V-h�F4[Qcfmv�#4/Ч��fw�w�/l          FǹW�ݭ��g6=�lE�5f�g��ya>���u����         ����3Ъٌ�63cs`�^��:�o�]|���2          �+�lA�FtR6��,��:�zvl�Ew姿W��         }قfZ����	�,��=��n^|�n�}wݝ?*_�          �لfZ5�![�lJ�Y�s��y�m��K�����         �ܓmh6�U;�)ͶԘ�6ۆ�._��~�o�*          sW6�=�r��٘�g���м��O}���          �W��ٌV-i��4[Sc�</�eh^`l�Ewݝ?*_�          �ٌf;Z5�![��1Ӟg��慵��wW~�{�         ��'��lH��4ds��1K<O�	���^���Ư�/D          �lH�%�Ӑ�i6��<�
�i��?����/@          �lI�)�ZӐj������y���S�s���|�         0>ι⮉��jNC���<�\��+tg^t{��         `�d[�p���tR6��4�(T/�n�Gn,_h          ��lL��tR���<2G��ҝ����          dkZ5���Q5f��z�to����          ~+�ӪE������2��м8�;���          ò=��Ԑ�j6�f��ǡya�|�7�/$          x4�����K���j��f��м ���e�          K��U��]5c4��a��7�G�          K��٤V�jȆՌ����*OY���u�/          XR٤f�Z5�![V3��P=�ݢ�n*_0          0]٦V��lZ�<����	͓~�{?V�P          ���F�jWC6�ٶ�y4+���	?����/          x��U�֐mk6�f��'C�D��U��/          XZ�Y�Z֐���sjh���m�y����|Q          �Ғ�j��U����}B��.�pa����+_          ��e��kն��3���B�{���          ��lX��5d��c�Ph���9�|          ����C������}5sh.���N�*�x          �)����s�������<�k��^wݝ?*�t          �)ٴ���s��uR��f�g���<yg^t{��         �L˶�j^C���Ě����{��g�O4          ̖l\��5<��١y�^����'          f۫�}SӿN:+��=C�D���M�[�����         �ٖ��zl�t����5#0�����I:���'          FE6�U���V���\�'�ෞR>�          0j�}��ؐ����y[h������         �Q���_�t����5�0��O�J+��]z�7�'          FU6�����!��lg�σ�yB�=���	         �Q�-l�Ȇ����@h����>�|�          `��&�jeC6�ff��<�Yo����C��         �\�Ml��U3��5�x�������O          �5��V�lȖ�,�yOh.�Ao>�|�          `��F�jgC6�f�&�����pӭ�'          �le��ِMm��f)ϝ���g\p[��          �\��l�Іlk�R����B�����'          �lf��6dck�¬�s���z���-_���I         ��"��lg�{ڐ�m���	�����'/��|B          `��v�jj�'�ys@h.�+����          ��*ڪ��ܚ�1+��
�)���]���O          �W��fK;�׆ln��5ӜCsA�y�G�'           �li��6d{k�1ۅ�B�˞�         �q�Mm�چlp��]a�>�I�u����         ��Ȧ6����6d�k�`�<��//8          ��lk��6d�kc~.ܺ�o\^h          W��w�![\�3�����+ˋ          �*۪���b��)\�m_�{y�         `�ek;�߆lr��5CsEh.؇>~oyq         `�ek[5�!�\3e^����G�          �����Zܐ�������z�*O�����          ������v���蚘CBs��8���          �����ܐ����_�����/,/$          P�w������s�Yה          �e�[��!�ݱ��	�d���Y^@          �_���}n�f7�ݱ��CsA>|ս��          �e�[5��0V�b�u�/{��          �L6�Ýn�O!ޱ�E������          K&�ܪ����Ŭ�1\�]�>��`          ��d�;��lx�����\�ox��X          ��d�[5���0�g��?����o*/          ����u�4���)d�;o��<�Ko�fy�          ��'ݪ�����Y+�Kx��9�-�          ��lu����C���n�����          xb�խ���0�f��0�@_��G�          X:���xÿ�l|�͜�r�/�]s�ʋ          ,��f�;��l|���]x��|LyA          ��+���7�}�sth��~��          �ҕ�n��l}���0<�W�yPy!          �e#��7d�;����y`�^qwy          �e#ު���9;���͎��           X����{C6�srv�:����          ,[��V�o��w�͍a��l���          f��7�r��tS�S�Qh�;O>�|�          ��Ȧ�j}C6�sf����\�|�          ��ʶw��� ωY-��0� 9����          3+����7d�-���qa��/�paw��V>X          `feۻp�
��l�G~�8�����          �#���7d<������^�@�          �ّ�o���l�Gv�
wx�w*           0�6�f��wR6�#9��w��S�/          0����C6������WyJw�^>8          `ve�ʪOh�'e<r��0pG�����          �!���8d<R�Ch��W�S>(          `4d�[��!ᑙK���d���          ��M7�~�����H��῅�;��?T>          `�d�;��l����92ܹ�+���t���          ��lWXq��&xR�³>_	wl�}-          0���C�³:��-���           F��K�h��I��ڼ7ܡWZ�|           �h[q�'�������N�C�ujy�         �і-�p���yAh��E��jy�         �і-p��l�g|��;��f۔w          ��	�C��3>��a�xy�         ����4�8��l�;r�M_/�4          07d\��!������t���;          �-�n��@+<)�����8���;          �-���!����\~˷�;          �-�W�pȖx����7�l��;
          �M�w�![�e>?���w���          �ܔ��p7�%^��Uh��}���          �ܔ�p��l��ټ'|�-�٥��          �ܖ��p?�)^f��0�?zQy�          ��-[��~8dS�L�ɡ���]����          s[��UC�-^�_�Fk>k��          �C6��qȶx���a��ϡ�          ��_{�@C<)��>?��3?^�)          `~�fx�#�/�yah������S          ����p��l��ڜ�����\�!          `~�|�Z�I�/��B����}�          �l��{␍�R��o��7��U_(�          0�|��{Z�I�gk��g�0p�O���          �l���␭���������           �l���␭������p���           �l���␭���Cs����py'          ��)�-�?�yM��u�߸�          ���-�p_�9~�sV��]�>���          ������œ�9~��0p��<���          ��;O��@[<)���=�+�����J��         ��-[��8ds����a��V[}��          �!����8d{<�9>����Y~S          `<dS<��l��=���:��/�)          0?z�@c<)��i�߆�Z|�g�o
          ��l��;���f�0p#�-�|w�^~S          `<dS�m�po�A^�9<���|q�         ��m�po�A^�"���S~3          `��w�1��l��x�n����-�          0^�-�C6�K<�=��e��f��          ��m�po�A^�y~8�ʪO-�          0��1�C�ȏ9����/����7          �S6���q��1�0pp���(�	          0��1�C�ȏ9����G�xn�M          ���pw�E~��I8���;�o          ��l����-r���7���           �)�=�$?�����z�7           �[����q�&�Q�-a���;�Q�8          0޲5�C6ɏ:��~\y�          �x��x�?�$?�<�p��          �-[���8d�����0p��O|��q          `�ek<��l��Y-|���7          ��9�C���l����nP�(          @��x�C�&7�o�­�yy�           )���9d���qa�w����F          R6��r�6����rԩ�          �l��;�mr3�_x�i��7
          ��9��������}ٝ�          �l��;�mr3�5|�U��~y�           )���9d�<0O_��J+�7          0U���=r�F���6|�:�o\�          �T���!�G�0���kyc           Se{<�#�l��������          L���p��Q~d.	_p�ы�          �*���9d�������3�(o          `�l��{䐍�#�0���Kyc           Se{<�#���#�'a�λ����           �:��/�ȓ�Q~d~���O��1          ���=�C6ʏ̿��/�����7          0U���=r�FybV��+��JyC           �l���䐭��,&���:�           T�A�C���������[�7          P�y�K;�{MYL�r���7          P�j��4ɓ�U^pؔń]vۯ�          �J6��]r�Vy�qS�|���           T�A���a��)�	����F           *��ā&yR��.�������7          P�y�K�*/�e�b�q�]Z�          @%��.9d����)�	'/���          �J6��]r�Vy�=S����          �d�<�%�l�<0e1ᴏ}��          �J6��]r�Vy�NYL8�ϔ7          P�y�K�*/���ńs����          �J6��]r�Vy��,&|����          �d�<�%�l���ń��ry#           �l���䐭�MYL�䦯�7          P�y�K�*/�Ŕń+o�ny#           �l���䐭򂿟��p�?(o          ���p��U^��S�����          �d�<�%�l��Ӕń��iy#           �l���䐭����p�^�          @%��.9d��,�           �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]�_          0U��eu          �O�&�vY          �S�ɡ]���8           0U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh����          ���M�:          Чj�C��          �����.��           }�69���0          @��M�:          Чj�C��          �����.?�q           `�69���0          @��M�:          Чj�C��          �����.��           }�69���0          @��M�:          Чj�C��          �����.��           }�69���0          @��M����           �Q�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]~�kq           `�69���0          @��M�:          Чj�C��#�          `:�69���0          @��M�:          Чj�C��          �����.��           }�69���0          @��M�:          Чj�C��          �����.��           }�69���0          @��M�:          Чj�C��          �����.��           }�69���0          @��M�:          Чj�C��          �����.��           }�69���0          @��M�:          Чj�C��\|1          �tTmrh�����           ��j�C��          �����.��           }�69���0          @��M�:          Чj�C��          �����.��           }�69���0          @��M�:          Чj�C��          �����.��           }�69���0          @��M�:          Чj�C��          �����.��           }�69���0          @��M����b          �����.��           }�69���0          @��M�:          Чj�C��          �����.��           }�69���0          @��M�:          Чj�C��          �����.��F           ���M�:          Чj�C��          �����.��           }�69���0          @��M�:          Чj�C��          �����.��           }�69����          LG�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh���f           ���M�:          Чj�C��          �����.�/          ���M�:          Чj�C��          �����.��           }�69���0          @��M��o�          ��Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]�_          0U��eu          �O�&�vY          �S�ɡ]���8           0U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vyO|1          �tTmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��|;           LC�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vyo|1          �tTmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vy��          `Z�69���0          @��M�:          Чj�C��          �����.��           }�69���0          @��M�:          Чj�C���;q           `�69���0          @��M�:          Чj�C��B|1          �tTmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh���          LG�&�vY          �S�ɡ]V�          �Tmrh��}7           LC�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U������           �Q�ɡ]V�          �Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a          �>U������5          ��Tmrh��a          �>U��eu          �O�&�vY          �S�ɡ]V�          �Tmrh��a   H���ο�	\{����gҽ.��/�i�c�o>v՝���7}�|�0��?�?�/���	�|����W~�\�����o�zw�g��]���v�����M�>��r���O?���g�����?����o�m������� 0T�{z.�   0��69���0   |����!�׼����'�Yw���u{�wh�xa=u����	Kד��\������^�y�ƛ���v�Wt/վ���|�q�ۏ?�{�wg��ɉ`���|����(�3 �In�e��m.x��땏	   �}���]V�  `���t�Ig��a�Xg�痏{�����`y��񌵺M7ߦ{����p�1ݻ���n����'�zs���Ҵ��*F��e   ]���]V�  `܂����>]^��`��nX>�.��i�xa�_���'���w��ݱ���]z��� ���/T?���2   ���w��.��   �������|���1׭��F�c�����V[������i�'���˫���������|�`:6y���Ϝ����yn��   ��W�.���  `����S_�u��!o=���Gv����_�]n���/�-�~Iy=��O>���ǟ�|�q���xT��~�v����ճ��h���������ؓ��>r�g����o��{Uw�itǼ���'�5���7�^���޽t�}�m�K��5�Q��mk��v��G�����/��? <���ɉ��p���v{��n�-�/挢���   ̾�w��.��   ���o��;��[������9�n���9׽��#��=��|����{u���o�6�|��̴���n���wG������?����_�����?���}g]92���s��ڡ{�{?�}�KQ>� ��n�ϻE�\���ʟ9�B�   ���]>���0   L�i���j۝��9���c�c�<�/0W�fY�p�-�7���O�_�/`�8❧N����%�-�R�U�|�|^`:.����4X���mk��~y�  ��W�.�eu   ���w^����h~���>|����'���wvL��
+��hi�|���>tMy��s�1�/ߣ����������+ L������Y3��   0����C��  ��u�g��m��v��s�j�?���?*������Y���uz�}ۉ��f��{u����u�<m�gt�`�	�l���[l۽p��-�~��_��z��v����n��v�^�˫���l����}���y�˻-�٩�l��6�lⶪ�1��񌵺�����y����{�?kf�`   FW��|h��a   x�v́��sQ�l��G7��p���۔���:���̮�>|m��5�����?^�}�W����������N�p�/Ln��˺��X����j���% �'��X���-��_ww~��&��U�C��GY[m�sw���(�   fR��|h��a   Xv����йh�7��||�讯�������N�u�Ig��G���{w�,�`�O���w�rn�����V{���5o}� �Į{�_�|�k����>���pDs]�   ����W'���0   ,-/�m��wѹ�ؓ�)�8��[���\k��:-����F�Z�^�|���������&���~��W�}��_�� �9���?[føˏ���   ���w��.��   �4�����t����)�8:��k�$VY���_�iy��h9�'���Q1���o�p���7����'=i��>���^���_���? <�[��q�se6�{��͋_^^�2   ����5���0   ,M7��p�ֳ�i~'�k���3�����1����ݩ�N��o>��=`�\r�}��xT�B��[�|�G��������`�W������ �f�-�/�̴笻Ay����κ��&I�  �(�~g�:   K��W���N:e�[=�qt�ş*��c�������U��Z�^�,��9�ҭ����;
�_���~@��ޯ�y2��9X����e   FA�;kh��a   X�s�E��s�>o8�||�����[~���5z4/��U�m�k��^T��G�(���o���k�C��<
�{Uy�`�AozW��d��k���c�_^��,  0
��YC��  ����7��n:���*߸�a�����hNXt~y;��z��/�ϣ`T���:�m'��{�_���3 L��S�-�̴u�{~y����~��
+����   ���w��.��   �,��ʽ��O�󮼣||���������������Q�nG=XN�~���>�VZ�������>�o->�����L�`����S^���   ���w��.��Gq    f�����z�ۨ�u�Y��Ϝx,�cW���嵩<�Y�)omǞrN���rq�G��zKy�g�3�����ۿ^�g H�_}g�3d�M�����N�����0��|�<   3���5���0   ,k���@�;�\���w)�8Y~�����a�����h;u���{z̕`9m�­��0�6��E���t��[���i�<w����G7���n�'�Z^�a�e   FA�;kh��a   X�>��O��Q�}zk���s7ؤ�.�|ӱ�y`�-�p4�S�,_~Ӄ����c����� �}��ώ�6N��Kwۧ��2   ����5���  `��>����{�~�|��`�W���&�������_�[R�i��K���.��;����Kpn��=�k�@��'��gf��~_�����C],��^U�����9n֧�k�{w�kպ>��v�#��?�]�ׁ���]u�E7�}t�m�30m����-��
1�xs�%��������   mJ����b  �����r�^�{�8��x=��k��@�=4��g����\�y���^m��׍y��u���Ba���Y�6�q�_Ga  �.H����b  ������|,\��ٿ���gg_��x=��;���@�=�����.���#��8�?����iף/�4��hےˬ�'�o�K�{/
�   tA�g-��0   4m<�+�l�E��xv�=��k�Uw=�v���W?���`˕-��=��|�͘�i��r;�;����o��  @�{�R}1  @��ka�����=�WU�-]��z��O�<�mϼ��g����|ӽ/��t����30mRXnޭ��<
�   tA�g-��0   4m<�+g|�Ƹ����I?��૾��o�<�m��o��.��re��=uAUNK���t���Բ㳰����T,���qϣ��  @�{�R}1  @��{a�r�#�ǽ�7O��E��W���?�y�ۺ��z��\uW�S�{�I13 ���f���q����  @�{�R}1  @Ӧ����.�>ތ������n�>l��Lw�0�+/�d�נ�6���?�c���Ea�9g]|K��X(,  �鞵T_L�   дi��\Yw�-��Ǔ��q~��,�}/�������\�O��SϿ6f`�ҕ���˭����h�)�u,�  �t�Z�/�a   hڴRX��u�p�F��c{���4t��rs.��=q_]��
��� L[������r+M��+�e   � ݳ��i   �6-�+g^ts���̳������7�ݧ�ܜ�_�4�+�|�͘�i��r�m��vq�SBa  �.H����b  �������˟��atϓ��k1(,����ܬ��_(�N8��iGW
��,�r�7lv��฿)��  @�{�R}1  @�F[X�^����Z�vƋo*,�9�<q�>��fm��Nqo]PeK��v(,�ϡǞ�65�  �t�Z�/�a   h�X
˕=86��a���[��u/f�y���/),��RXn�ѧ\������r�q�Mq_SKa  �.H����b  �����\Ye�:�5�f��}�3�4S���ax),7�{_�{���y?�`ڠ�<��{�qO���  @�{�R}1  @Ӧ��\�n��as�ŷ��6�f�aƸ�/͡�CKa�Y�}�/qo]q����� L:SX^~����&�0c�{����  @�{�R}1  @Ӧ�������FN�Qm�j������Ka�y3�<K�_l��13 ���)w�ï��6{�O�(,  �鞵T_L�   д)-,W������a�����馛.��K
�0���7����u�r+M���6(,O���x��k�y�^�Ia  �.H�����{�  �m/����f+�^|��Y�m���M����k�8t�+/,���fɥW����X8f`��XG
�ˮ�j��E<�a1����}�ۣ/�$f   �6�{�R}1  @Ӧ��\�`�m�̰������I����s�8t��r�V��N�_L�0}���Aayl��"�-���e   � ݳ��i   �֏�r�zro�6�^v[�߰H{���2/��歷�Vq]Q��N���G��W~���&(,  �鞵T_|�r    Z����,,�ٯJs���?��i?�),�9��^����y�l��nq]��>����RX^n��b��x��ŖX&foң/�$�  �6�{�R}1  @��YX������Yh�����A����ax),7o�=���[x)�`�SX�f�?�Z1�|��MSX  ��=k����  �i�,,WN=�8?l6�l����K{���s���SXn޾�����^�@��������w=W�<ˬ1s�  �t�Z�/�a   hZ�˕�v�'c��w�Iq]��19�e^
��;�s����+��RX^~��1� ]v�C1k��  �t�Z�/�a   hZ����K.�3l.��޸��J{���2/����s�	q]q�M��� �
��yW�s�Ma  �.H����b  ��5UX����g��U�_����m�9�c��k�7����(ξ��⤳�*?ἑ�r��a�>_y��i\W\|����w���ɷ��^�E�����	�Da�y��󐸿��~���|��>�]q��o���pqΥ�+N8���c�.�<�ؑ�����Ï?wd��y^���?��x,F�O�4����G^+.������-�9��⠣��?���#O+�=���Knyo���V~?�|���q�_:��R�nN��.RX�;弫c�APX  ��=k����  �iM��x�x�a��bK��uQ�?9�屩
�U�j��O)6����2˯R�8�L�ڎք	�+��V����wM�^)^����ar
���z�=�����G^�����${��7;�z@1q����?^��Xr���vأ8�;7*�}��_�Yq��w��{�H�1]�Ѩ~?��ʚŷv;p�?����7����W�Y��k����Z���oP��   ��=k����  �iM�+]so<ް�h�����&e����7��Y=-y�u6�װ)U9�*}���1(,7��Y���<�a�=-��:��a'˭�Z�f����K�r�'0�Ӛ��љ�T���f�z��J��]v�9S���ާߍ�?�����Cu�����1_[���7�.f$�e   � ݳ���}�?   к��_�{�$͎Vמ~5��BT�_��ܓ�sa9�N�y��#_U����ֶ��uǞvi��뿌y�6���?�_��υ�{���A����O�-�ּ��?ǝ~Y��2+��Դ�kmX�}�Ż�W� g��ߗ��wT�.M��|^t�}1W/Ͽ��#
J�<��S�L��t�C�_?kұ��������A{�ş������r�ش;{�Xb��c�A�>�{�ï�"]   ��{�R}1  @��(,WV��n<��G���"e��l���EU���+�,��h�x���z��W����Ea�y3�<K�_T�H��%�?��b�]���g\x����YǛG_�i��N���ж��x�Wo��?b֯Z��q*�S��L�,�Т�]w�w���՜s���u���P�ϠM���{O1݄	1Ӵ&]   ����R}1  @��*,Wf�u�x�aS=�1�R��),��i\?R4Kר��Zv����n��aڠ�ܼf�)���`��yZp��v�,�T?�������	g]�ɒ����Wu҅�k��K��҇�Xl��q�K��o�헮|s�X���<�����w��ro+��F�פÏ?/f�V�k   �tY�/�a   hZ���?�}<��Yd���� �ܴ^X��̶��+�k3,�]a�q]���),7��~������`���?+����zt���|�x���Ľ�����ަq�]2�l��ڏ9��⠣�y��>��X�0Ì����d�oŽ��"�-��u��{������u���p�ϠM����[�G�1-K�
   *�>�T_L�   д6˕��|.�l��q���NnZ-,��؛�ĵ6��dXU_������e|RXn�7=����N=�ڡ���ꉲ]s_��0��gF�`��8����fq�����K��v]�~��7�B��]���[[������ES�i]�^   PI����b  ���]X�~���J��85�o�R��UO!Ls��)�]�E���~��ބ�6�j��������ܬ�N�N�[W<���1�x���,6�z�x���3�R����'�?�A��y�8�F~�o�ݮ��*�<�ظ�ap�9ߍ{�V\e�x��[�(=���ϖ_ib��E\yW��X]y���!ǜUl��.#����E���H�Ϡu������ʱ�]�͟�k   �tY�/�_�   ���(�ivjl��v�<��۞����qrUa9͍G��/�v;��C?�2q�����+n}���V����:r��?�]���?..������(�\z�8�3�4sqΥ�֮�˫.,�}�q1�0�n�}�޺`�%���Ǜs/�m� �5�UA����?��x�/��韊Î=��n���q��V��^����sv�߹1����*���#�{�CO.���^��l\,��2q�im��z�O��O��Xm��b���
�)��z�ǜzq���u�y��*,���V�=3�m[UXN����7�l��y��t�   ���#K��4   MTa���B��s��~��q���MnZ),������K-�A�l�ͷG�#��_�[/�Zv�x�~�q����>����?��fU�Ŵ�.�y��c��d�CO�{�R��w������u�{�W�ƛ��7��2�&��kλ�����U��o���x��=��/��.���r�]G�����o�η`�҆�����O]��c�~�α��+���7=3�x/,_{�3�<�-��_K�   *�>�T_L�   дA�+�\�f�ŗ�{��or��:[�O���M7aB�?�9�<�w�{O<�h����,���x�~Xd�F
��7���<��Gq_]q�����ǃ�CU�5�{Ju҅�\���/��Xw�-�q�T���ɺ�|]Q}�'e��.�Jq�C���Ʃ�];�w�t�~�y�Y��T} *e���{�X���?�l��?o�����J���M���v:Y��   PI����b  ����|�o���M��V�_�R�ɍ���Y���/UA����=V�{v<G�s�E�/�����q_]��"K������|P,�ªq�Sj�����_�U&��?��\Z��D���WO�~�����w�Rˬ���m[k�Mc�Ai��\�J�/�D�1(W��X�:h㹰\�="��,]C   ����R}1  @�]X��v�u���cω�kS�5��\X�
�i�������ҹ��Ig]��/�wn</�Ia�9]|��:���y��|�#��Ҟ��Z�m�5%��@UWK�;�~P�;��'f��L��<-�oj�s��z�o�6(]�^a���߇�Ja�7���K�   *�>�T_L�   д.�+��o<ﰹ����ڒ2Mn����x�o��9�<Ť�~�=�6�j�x�~�
^���ft��V���1�0���3�<K�0a���W~�7��(QV�����u<� �u��1�Za�5��?��{&w��O�,�d<��J��V^=�������?��f�#f��ny<���\X�u�����5  �J��,��0   4�+���2˯�=l��ޱHy&7�,�ƹav�%�Ľ��W���O��yߟ2�U{t|<7�Ea���wM�Sl��>1�0k��\i�[V[c�x�����[�s��g?,��0!f�RՓ�ӹ���Y�����A��~�_�<m�Ra�k��ԥ'�+,���jk�|Sc�]�/�g�)2��?ǦT�֔t  ���YK��4   M�Ra�2��3���ŗ\6�)���[a��;����mw�;���N<�x�~�m�#�
��Xl�e➺��g>���USe�y�[�x�g��ZW��X<����Г��ڴښ��lSj��7��i�AG�s��+�>�m�\zk�7�(,W����U�<3�x.,O�Z4�l[���|   Цt�Z�/�a   hZ�
�O���x�a��V;��5-e��x*,?����9�'n�����6������Ӟ��pPX�s/�^�O��A1�z�{��)�^��/U7�wj]p��|mh��]}�A:W�N:�ʘe,���/��V=�:�����/flR�t�e�{SXΞx�Ә   ڔ�YK���Y   @�^�`���0۔�n)ð9攋����rLn����[���Xv�U��i��&��7ᤳ��t5�]q{<?���:\X>踘��[q���.x����荟��Xa���>���˭��O���`<�Ԛ0a����ٴ���4f�Ͻ��x��s��1�h=��O�q���O�-��K��?fl�#���w��D�8h��ԑ���ub�A�La��Oc>   hS�g-��0   4����JU~K9��m��ה�ar�W�a����뷑�v��7������m�9c�~�~����0f�����Ko�{��ϸ<fV[�{�g?\t�=��6ǜs��O����1��I���l�25���Î;'��;{#s��{�h����'����ro
���|�  @��=k����  �i]-,WV_{Øeؤ�5%�r㡰\}�~�[��L�Д��;2��UV_7��nSX�ꉪ3�0c�ˠ���&1�:��3�>��zBv:gv�爘�μ�xΦl��61��Xd��ڴ�.��l������7�/�l�ضK� �k��s����۞��Ma9SX  ��H����b  ��u��\�{��b�a��R�Ž5!�r�^X~�՟��6{�[����f1C�n����	��{D�@w),��Z�m��w?�V�<�ο⎸�~9䘳�y�p�/��P�^{��_����-��3�C�3*��MS�a��nx(k��a�A��]a���m
˽u�o]),W�))   �)ݳ��i   �����+�>�6�m�S�_��sOn��l�m�W�<����i+��F�ӄ�/�%f����㠣N�{��.�.fF���m1�\��}�˝���݄7~�o1C�l��y�m�o���U�;��MM�Q��˹���5[l�K�ضKox0�kZ����۟��Ma9�Ja���~�  @��=k����  �i]/,W���N��͉g]��O鼓�i����08霫➚�^�9����'�<M�a�����u�A�(,O��.�)�<6fVU7�_Yl�x�&����1K���������7WT���9�V�}+��:'���x�Ah���h����V;�<mSX�k��^,�l��2   ]��YK��4   M��r�΋��M�_���9�a-,��ݿ�����&�<ˬ1G.��ᘩ)��}x�A�(,O��=f��7�:fV�_qG�g?��]��LY���,���k��R=���!��m՟��1�7*$�7�챸솇b��m�ݮ1Oۮ���oкRX^Ua9z���c>   hS�g-��0   4mX
˕���!f6io���7�a-,{�C�~������mx��_�LM���7c�Eay�]tͽ1w�������c�=�&��7Ŝs���O7�=)��I�^������|��Ǟ��O{||<w�����/����1a����6����;��Ma�����9ۦ�  @�{�R}���   ж�FYXN�������|�d�eV�{�t���8�Lq��n�����&�{ȉ1K[�~R\U�N9��:^XN���ĳ�����*+���?���j���{�9�7��i�=�J��O�~k�x�~Xy�:��V]�t�����9�'g��A�e��bƶU����5m��<a��;�����S�AYh�%bζU��   ڔ�YK��4   M��r%�6[�G���J��0�'��A�K�N;�ژ�-�l����IW��H�Bw(,�M���M��)���XV~h��q��������7���~���3��?5&���x�&,���ū�>�hۊ��3N�ɒ�X��Q1c��w��<aZa�7��La  �.H����b  ��ca��'ߎ���\�75�y&7l���#�i���t�Ӗ}>!�j��+M�Y���ѻ��{�9�'f��C+�����6��\�����ߴw>�S��oǝ~i<�Ը�㹚���o}�ǘ�MW��X�7�.}�@��4el۠
�]y���w>�ZW
˫��^�7(/��   _J����b  ��ca�R�}S�a���?���R�����r+��ѴG��$�iˠ��>�ܫc�Aa���}��Ų+�3vEe�.h��ʕ�ox(fh�<��3��
+��=5θ��x�&���Fŋ��}�Ӧ�6�2����l�AaYa���La   �"ݳ��i   �6���t&��A����Ƽ�&�mJ��On����<��+W���b��TOxN����2+�<t��rv�3�ǟqY��bK�l]Q�n�{R��x��ӕ+��PI[����7���T�T�t��-��
�}O�3��{��پԥoPXVX�Ea9����e   � ݳ��i   �6���	���.Xy��c�a�ϒh:�䆩�<��+Ͽ�"1O��rW�ֆ��x&fb����o�����]�:�Xh�n<U�l��.#�o�~ƃ�_�i�wf�aƘ�-�o�u��o�||<����c�y�0�,�W��X�Ֆͷ�9f�,��Rqf��{�Lay��c�A����mSX  ��=k����  �i�),O?�q�+f�}��{���~qoc��=�a),���+O\'fjӤ�~��a�v���i������P<<�G#O?�+�o�yH���s�=o<WUO�����8�r̙q�MXa��c��l��1W���D�V�^N>窘��?�j�T�{���� (,+,����),  �_�{�R}1  @�FSX����ӯ�"�6�\zk��X��NnX
�U"�o��[�3����3fk�so�*�b��\X��ħ����o(ξ���/�n�x|�)�nq�Q��rb���G;�~`��f�+����Ϋ���s�������>,\x�x������ڲ�އ�\Mx��Sb����h�^��랲ܥ�),+,����-���1g��  �t�Z�/�a   hڨ
�CPr����b�aS}�~��h�cN������ϼ��e�����6�'�u�1����26��G�x~����������ע)Ǟzq�і�8%�jB���aJ,���x�AXc퍊G��$�l�m��T��AaYa����'��A�7(
�   �鞵T_L�   д��g�i�8�5])yL����JǛ�0��?�䘽-qj�նZ4�kCU�H�,���Zg�-��\}O|�ƻM��V�&M�>��r���S�s5a��w����l�1(3�<K����6i�ͷ�y�kAaYa���l�%��9ۦ�  @�{�R}1  @�FSX�J%i����h˸�a���+Ž�F:�䆡�<����m9��+b��-���1_[��屘��QX�Yg�����c���y?�_ӂ�^�y�6M��gb���z�51W�{ޘaJT��t�A�s��cަ\s�S1Gz� (,+,����-��1g��  �t�Z�/�a   hڨ
˳�g�j����&Ua#�훤cMn��g�s]QdR�6��1[�V[c���-;�~P���(,wC�D����۟�)�OӒSν:^�&���k1K[��>�j�=O�s�����+��V��9�þ�����+��m
�
˽t������|���   ��YK��4   MMay�Yg��]��ǿ��6S�U��8��zay������T=}1ek۠�^�S.Gayp6�|��«��˴��.�z5�I�,mi��5�_p]�1VG�x~<~�T����κ�ڹ��;�׶MaYa���l����9ۦ�  @�{�R}1  @�FSX���>�v�w?�e�<��q_'cr]/,�=�|1w��z�͘�m�o�혯M>�a��`(,���k���S��ٟF~���֤���e�Ӗ��~����c�J�]��{/����~��������A|� (,+,����),  �_�{�R}1  @�FSX�m�9�l�����t�M��u�1&������3��|�*��|m:�+b6Ca�}ˮ�J|/����ѧkִ��A�Ӗ[�1�j��ˮs��U�<��E�����ɸ��NaYa���l����9ۦ�  @�{�R}�_�   в�>�����s�g��V���4L�_ib�[��'7a��q��8Ἐ�m�VO��ڶ��G�|m�p��b6�O���ψ�Ǵn���׫i)K��x��׼���b�������e�vR��x֙��͏�|M�La��gc�A{��n��Xg��oP�\f���mO��y�   mJ����b  �����<ǜs��a�܊��}�=<&�����\X^m��c涽܇bX?�s�	1_�f�q��/�3�}].,o���=O��ӽO�S���{�>��ew=�f|O�e�(�u�T�>�v�֤��~.f��W[;��VXy��{��3),+,����),  �_�{�R}1  @�FSX�s�y��0�q���ކ�e7<�6�47���_|��c�Ax��?Ōm;��Sc���p��c>����rUBN���r+���U��>G�Ǵ���!^��-��1O����1[��;�Ҙe��?��x�ap�q��=�7
�
˽(,gK-�b��6�e   � ݳ���/  ���<�������ar���Ľ���۸�/���U��47h7��|�۶f�1��Î;;fl�AG��Ѿ�;\X���b�3���N��o�˴�ƻ��Q�]|���M�����G�,cu�mO���5�ݤxx�Gqo�EW
�W��H�״];TXN��Ɏ��\w�oP�RX~���c>   hS�g-��0   4m4���/���/�.�o�L7aB�ۗ���i~PN:����m���S�A8��cƶm���1��ray�CO��{9���㱺���{��{����4m��y�4���[}+f�������t�MW�q��q�²�r/
���ˮs�Ma  �.H����b  �����<ϼ���a��=�{&���n�[%�~r]-,���1o�Zd�o�9�;1c۪�֥|��Ӆ��N����F�m��eU�3�eZ�����kӴ.|�d��V^=f�ՇP�9��6�ڳx��{f
�
˽(,gK/�r��6�e   � ݳ��i   �6���|��aU�G���a��a'ǽ��N���媌��K�ܮ�+��1#��ra��~&}����m1�L3�cv��l�2-Yt��i��3�����>�ٚ4����,Sⶇ^��F�/�Hq�-��}+�e��^�SX�$��e   ��t�Z�/�a   h�h
��ο`�fU	;�u�\s��}��M��*���ZU��r��m[y��c�A�Ra���?�iW��>e��������e'�qY�˴���>^�6T�+�_`�Kٚ��O�5�Sb�-w��VS���.RXVX�Ea9[fy�e   �R�g-��0   4m4����_(��g^�E��y����}��L������B�
]*,���_�ra��#N��Gk������^�6��w?�f�4��~ߏ)q��/�s�U&�S<��{q��DaYa���La   �"ݳ��i   �6��r���iv�]r�}q��d��jO�5��ba��;��Ya�Ͷ��K��*K�H�:]X>�Ԙy,[b�x�Zh���>ƻ+nz8^�u�C/��cJm���<�l�fyb{��PXVX�+����4��eWX%�l��2   ]��YK��4   MMay����������<L�\w������'����^���3B�
�;�_�H��\X>��b��ރ/�cwٞ�2��|�U�ZЬ�o|(�S�ƻ��3z�Yq��@aYa���l�W�9ۦ�  @�{�R}1  @�FSX^p����x��Vߊ�&_�s�?�����<┘sv�u��q�TX���1#��ra��N����΋�ﲪ���2^�w�I�:Ь3/�1�S�+��&l��.�[?����.SXVX�Ea9SX  ��H����b  �����<-|��K/�>LF��Ĵ�A�z��c�A�}�#b�A�Ra�*�����˅�C�>#f�k��I<G�����^ƣ�w�#^�U�NH���x�����:+��F���Ľwմ^X�e�Cc��),�ֽ��j1g��  �t�Z�/�a   h��
�K����E�~K��U&�s��3B�
ˋ.�T�H�:]X>�̘yJ<�ʧ�]�ɖ;ƽ�G���5�Y�zb|?��5�=�7^�1���uw<��E
�
˽t�������|���   ��YK��4   MMay�Ŗ����M�L��o���*æ����a'Ō�Х�����3Ү.�=���yJ�s���<]VeN{o����Ӭw�?�����G�s�'�^vk�{�(,+,����-��Ę�m
�   tA�g-��0   4m4��i��G�|a��I�� Ϳ��1� t�i1� t��<�L3Ō��˅�Î;;f���w<W�=����^ƓA��^o㭊K��o(\z��}S=	9���֧(�Zf�x�ǓcO�8�K��{QX�VXy���m
�   tA�g-�?���   �{�����Q�j�%����Վ����x��<H�׵������e�ٶcN�(f�w��1'�y����M~�91��x�g(��x��Zn�U�^Ɠ9�'�-�n�S���y���-�x����~�����),?�5�;���b�A{��Oc޶����1ߠt���E�   mJ����b  �����������l��׍�b<H����)� u�1� t����G��9iO��ǟ3O���x&���=�����b�Yf��nˊ��s1���b����}<9�3���@aYa���l�Uֈ9ۦ�  @�{�R}1  @�FSX^|�e��xץ'��S��|�˯�z��PX�z���1'��ra��΋����O����z%�e<�n��T�ڔ�����w���r�~�ONaYa���La   �"ݳ��i   �6���K/gǻ����ǰK{��~���qP�T`��������Ĝ��˅�#O<?fV]+���f�c��/�;�e�������m[�F�TOZ^x�%�ON=�ARXVX�Ea9����r�Ma  �.H�����G�  �m������2+��i�^��0K���?�c�8(UI8��Ï;'f�g^�<�=?�pa��?;=�a<o���~q/���K�k�����c>���W^,���x�Ǔ�|����A٧#��n~4�kZW
�7��\�7h]*,�|�ҕ5=��1   �)ݳ��i   �6���Rˬg�{xl�.�*�q�R�A�Ra��cΌ�w�6�=].,}�1s?�z����]vɵ�Ž��>�}�g�|���~��b�����xr�o���4_X�󐘧m
˽),g
�   tA�g-��0   4m4�套])�NK6�d�xm�Q�� �<�,1� t��|�Q�ǌ�����)�=].,s�wb�~�l�����lқ��{F�|���Ƕ���k1��i���0^����q߃����܋�rV�N9ۦ�  @�{�R}�_�   вW>����2˯g�5�,�d�>�&�m��{ޘsF�2���3ʻ?��������.,�򝘹�^����>G��Uk��I�˰�n���>�t���l4��Ύ��xqБ��}�m��;TX����^*,�|���+*,�|�ҙ�����   �M鞵T_L�   д���]a�8;���#t�>�&�m�Xhјs�?���q�=䄘q&L�3�.��?��Ƈb�.;���^��,���ئ�ξ2f�Y\yG|?Ƌ�^�q�w���{QX�V��N��6�e   � ݳ��i   �6���r+*,�����h��}��K.s�Ig]3�^3��s�3�.���c��b�.��7�^������צ�;)f�yw=�F�~o��6;���&�e��^��UW_7�l��2   ]��YK��4   MMay��V��Ӫa�Z���AZc��b�A8��kb�A�}�#b�AXf�UbFڥ��זZvŘ��^tɸ�a��"����i�M���hG�w�Ͷ�)�7������{n�²�r/
˙�2   �E�g-��0   4mt��qvZV=�/]�a��3H��W�9g|���qv����q��Aa�����c�.�����2L�]aո�6�8�L����T���3�V]c��׶(,+,����UnSζ),  �鞵T_L�   д��WXy�8;�����zu]�� p��1� �s�-1� t����=�i��r�1�^�t٠�����z��}����_��h�M�L*Yl����~%�
�
˽(,g
�   �鞵T_L�   �4��)W�g�y�xͺ,�e�N=��sν�֘q6���1� q¹1#�RX���h˘��^z���^������ԶÎ;;�}�~��N}�fj���q�m�����R��),�ֵ��jk�s�Ma  �.H����b  ��),O����Y��}ҕ7=s��\3�F�o3B����2���7~YL�0!��Ͷ�)�e��qOm[}�c>���ޫa3��3o�ǸǦ),+,����),  �_�{�R}1  @����Ig_�[W�=ҽO�s�	g^3BUI��^�i���׻��;b�.;��k�^���cό��7~�/1#���3�u�<85�A�e��^���kns�Ma  �.H����b  ��),��n��]�����/��3OG=�b�A�J�����3�.��޾��1[�=���q/]v���ƽµ�?32x�r||φ��[��մ�\7�e��^�VX^}-�e   �R�g-��0   4Ma��Zo�x��&e�UW_7fm��G����+M��6��|�Oa���>��b�E�����s���e7�;)�e�;�Ę�n�����f�#�w]7�|�=5�#N�yڦ���܋�r��  @�{�R}1  @���k��װKR�A�s��cֶ�{�	1� ,���1c���h����),�.�iG��#N�{�Io�:�cV\e����x����n�e|���W?�{jҡǞ��MaYa������0�l��2   ]��YK��4   MSX��y/^�.I���+n�Y۶�~G�|�0ϼ�ǌm���>�}
ˣs�1�(��ŵ�?��U��>g�� <��ocF���Ύ�_�UO�N{iґ'���MaYa���La  �.H����b  ��),��y���cW�̃��cֶ�ہ1� �<�,1c۪2y�G��Go�5׏9�j�9�.^����^�h�5֋�����2�=���B��r+�����ǝvq��6�e��^�VX^c��bζ),  �鞵T_��|1   ���Q��,_�#N�ײR�.�w�c�6m��=c��}���|����ǌ��S�Ay�Ob�.�~�}�^���pG�� ,��b1#���/���s���{�5o�}�C�N>�ʘ�mUa9�kZW
�7��\�7hOu����JW
������   ڔ�YK��4   MSXn��[��砥�]�÷��y۴��;�lm{�G���ڶ�2+�|���؜q�u1k��}�Mq/]s�%7���r�5�Ĝ����:q����۴Xo����7�:��W��p1����+V]}ݘ�Ig\Ѝ�[�*,��A1O��{SX��  �t�Z�/�a   h�k+,7i�eW��t�R�.���;b�6m��61[�&�����m��{D��`(,��V���v�c��(�K��q�>(U7���z��kw���׶�����F]|���I]�0�²�r/
˙�2   ]��YK��4   MSXnV��'L��려�]���~�i��6�����~�m�E6����S1�����]5,�ۮ]�{�x#��/�{���u[o�-�k�t�Y��ru�����6�ު\y�#1_��{SXΪ����mSX  ��=k����  �i
�ͫ�R��려�]��Jc涬2q���mw=�z�׶�������PX�2W��X��e�{V�K�l���1����!1'���?�����~_ߦ;~�Xf��c�Az�?ƼM����b��]v��1_��{SX���x똳m
�   tA�g-��0   4Ma�ǝzq�����u��ǝ3�e�E����v����|m��S7�k
�Sn�������{&Žtű�^s���C�:��=������6��)U/]�#N87��m�~���V���O�%fm�u�?s�����5Ma�7��l�-w�9ۦ�  @�{�R}1  @��۳�n��۶��+���m�n�銏~�?1[�����3.�.fcp���r+��w�"�-Y����^���>����.�Y��e�?��t�E7��L��ߧɯחXh���A��0�r�m��g���t�C/�,m����c��),����m�ݮ1g��  �t�Z�/�a   h��r�V[s�x�۔ruɆ�ns�e�[����t��7�lmz���l���Թ勺c�.�y��^�b�%���e�g)���â�s>�,��վ��h��=my��N���=_?(g\p]�٦y�?fk�ϼ���Kn�����ܛ�r����ĜmSX  ��=k����  �i
��z�'���k�y�unK��%_sO�ݖ��a�զ�O� fk�zms1X
�S�3.���lPO:�}9>f��6�.f[l�K�W�A����JǮT�M��+n|(fm�K-s5��W?�Y�6�oYPX�Ma9�ʿ7���   Ƈt�Z�/~��r    Z6��r�e�����:�%e�A��/�ၘ�M{xL�֖3�s}��`���C�.�p�m��j	�S�~�2h7�7)f��Ϲ2��s/�%�g�%������_���7�� ]{��1kV]}ݘ�I]�9_}�#�k��{��ƻ�����]��m�����s��~�  @��=k����  �i�*,�Y��uR��9��Sb�6�t�1S���q����-�P���u��|�)��]�»�-f�i����Z���^�`�����Iż]���?*f�u���~~��*��sT6�b�83h�*-k�b���4��1O�?�옭i
˽),g�zB�ٶ�{?�  �6�{�R}1  @Ӟ{��=�W-��q��3���OY���w~��a�����ڴ����lm8�ԋb&�˅�#N87f��'�˺Z
?��3c�A[{��b�.����Y�<q��kmg��s��󥫿�x�����+�m�i�_�4m�E��y�T@S�&=��Ͽ���6��޺VX>��n����7c>   hS�g-��0   4�|\�G����+�2�6�l�x͛�rt�^�7m�mw�yڴ�2+�lM������S���u��\��R�.�e�C�^�춇^�{��d��v�1�|'f�w�?�\uˣqfJUBK��RU�Ns]��ST�|����i�L\'�iӮ{�5i�u6�YAa��������m���Wb>   hS�g-��0   4��X�GM>���y��bK,�ySR�.T)�e�g�)fkځG���].,t�1s׵��wj-���#O�M{�ͷ�9�킮?���.��+o�}��˭�j<��λ��8��-�L�܄�?����i�n�c�Ӧ�?�U}�E�1(
˽u��|��WĜm���I1   �)ݳ��i   �v�=߯ݣ&����<S�gߋ׼))CW��AqM�w�c��Lz�W1W�f�a��w3�/������ɧ)s�UO,N����="�e���퉘��\z���w~s����3����h�b[=�9�v�)�^3��
+��߆���)S�VZu͘�	���1� ),�ֵ��Wu�ߡ��x*�뇻{=�  �W�{�R}1  @Ӯ����=j��뿈��ǥ���{����~�дw������I�LM;��b���~߻.���cc�ap�	��=u�����2H���S�.Xm����~��1��\~�1�v��87���6��}�ǟ��_�{1�s���t��g����ϋ��4�<��l���s�3�4�o��廷>���rvݝOǜm���Gb��u磯�����   &�����[}1  @Ӫ��N��_U=8��?m��ҹ���cΌ�h�ӯ}���+o���T=M�_�O�Cw<��'��낦ʕmYsݍ㾺�*lNz�q/�rɵ�Ƭ]��&�t��"����1㗦�~��o�:�N��v�+�3���c�f[���ӽO��݆6?���?�瘯_^�����˭�=h]�PH�+�>���}O�s������V���/ϱ�*k��[��  ���U'S_L�   д=�;�v��T_ל�鯱��T:o����/�D�KS.�ၘ�yj�Ԥ�.�^�B�t���l��n1���S+�b�Ͷ�{�V��vŲ+�:�����-Փ[S�ɝz��q�v���x�d�eVy�q:� �x��1o�,�ؒ�m���wb���t落�_6�b�x�.�>Й2�Ó>�y۶�J��|�2�_Ŝm;��Sb��Q���z�馛.�   *_������0   4m��׭ݣ&U�%���?E�	�]w�Uwƽ4��O�9��Ɠ#'����t�����.Xo�-c�arΥ7ǽu�ǟ�2(m���s�=�Hi8�o�Ig_3M��?Kc�PL��~�ïĬ��ׁ�����O�-�jۑ'��������_pm�=h<�n�۶ŗZ6��~��1g۶�q��oJ=<��x�cO�(�   *�^�T_�Q�b   h�ǿ����O��_�܊��c���b�Yg��C?�s��7�*�	묿Y�І%�^>fj�CϽs�=W��p|�`��֊��MO���n2�eP�\g���k9����){�T��U��8�/U�-���K��'kP�}�1g����k�m�����i�%��٦�aǞ��%'�qi�>h�>�F�۶\$���Y<fmS�!e�R��iu���/�3�   *_������0   4��g߫ݟ�����1����g�{��|���|\L�0}�S��6�1CӞz��1OSN=��1�tf�Z�h��ֶw>��b���{��X�x񽿉��;y5����X���/�}�K�A�іO�=�;��t���cι����x�Ax��s�C���tζ��1_�.��oJCY�R=�>����|�~&�|��ŶߎY�v�So�|c��>����{  ���'K��4   M:�[k���\���84�*���aj�s�s/�%�	Wﱘ�I��wu�҄ͷ�9f���\��y�Yb�aT��K{�Ͷ�V�ˠ�sб1gW���c��(�eJUO��x������7�:��.���x�o��*k/��x̶}��Č�p����s���k���V��)ߔ��re�}�{����h�;���1㠜}�1g۶��>1�Xr���3�4S��O�%�   ���=e����  �I�l�y��������Cs�����^L�t�a��=����w�;��Il�M��os�5O�������������������e�v�� ���-Xp���˪rTVy�'���M^���Sν�Xq�����L3�\<����v�mO���܊����e<n��}��1��Z~������>f��<&f��~��b�������i�i/�V�|Iy��s�2�?���s��0�}z<f���Ί3   0�tOY�/�a   h�Ó>�ݛ~�馛n����Gs�Z,�&�ä*�,��qo�4��3�|����	����	7��L�@�M\k��~v�m�s��V]3��.�ၸ�A��c�a���׫'Wgz��_�����U��*?��;�7>X�}бc.)O��˿�an襘a��^n��O<�mj�{���|�R�S�A8��Sb�or�#����1��^#����}��3h�|\�;�ϒ�q����m�������{��������_�O�  �ɥ��R}1  @S�;�ڽ�hy¹�x4�͟��b����ǔH�6�?�r�[��t����M��@)C��~����t_�U��=�s/�%�V����.�e�Y��_��';�[�'E�gS��ڱz��c��Xb�劧^�i<~9���kj���F�\�t���Ye��w-�}⍘���Fm��ؾ5�z�x?��;5[b鸯A�R���0G�8H���T�:�����=\��U�>��7~Hʇ  �t_Y�/�a   hB�t�9皧vo:s�=o����Ks�y���~L�t�a��������k�+�����zb笳�3�S�d�t~��*����K���S�K_�?Z���n�ѯ�;�g��x����G����k֤���)�X-��R#J�h�*׉��F�$�t�A��s�;HUq���\_�t��G�4^}H�*SVky�}��	�����^,��Ax��ߓAj�FF��ωm�Ͷ�y�*�_|��1g�����G�����Qq   �toY�/�a   hB�����t�N>��x\�u�w��c�ұ��1�\��O~r<w?���M��!�������FU8{������f[}+�˶�q���A��_�|�)圖-��J#�t͚��G�+�RU�?��)߽�јcjT/N��6���O���.�ؒ��]��^{?��!���k�sZ?�"�o묿�����3r]}��
+O����E]r䃍i�   ����R}1  @�]qヵ{�)q�%7��Ӭ~|z:�0�}���>���&�s�C�T�t�~�����ǿ����X��:ǝvq��0�Jbm<��8줸�A���Gb�i�|,T<2�x��V=56e�;�ok��~?av��,���?�suE���ۄ�}���u�%���������۟�9�����Zl˕7=�   _'�_��i   �������{��j��S��k��Y[n�K|?F+s�m��vq��R}����&�{j5����[{�M�7~�O���������L3�<��iì_�i۩�~7�g��-�nι�)���x�ڒrM���M?�v׽�����t<W�\���b����[���.��
�u�����U�A����c�A���{b�A��W�5�a��y��x�eq   �K��,��0   �K����L��I�F���t>������c4��ƃmv�=�_��d�'S�sO����5��_6�l۾g�]�sf|o�l�V-^~���~�����u]*-�{�-1�b�Yg����J��/��uH��7����������^��E;�_�ð�qƙ��~���e�[9�~P�x��Z�Ax��b�A[i�5c�.x�城|@$�{xL�   |�t�Y�/�a   ��+�ӽh�����?�}<7�x��O�駟!��$o�h�鋓[m���V�j�|T=�;������������{;�_i�≗~�6̺��ϱ8�ԋ�~��\3�w�η@�O ����9��⭟�K<�XUŽt��q�a'�su�+�}�Բ+ƽ�9朻����>�[q�83(�7x��m��W�]l��1_TO�O���۟���n����   ���5K��4   S��{'����[`����9hƔ~Ex:�xr��'�}�˂/�
8ZU��zZs:v�T��tn��-�Mj��m�{���ˮ�?�qXU�Ӵ�aP}=~W�z�yWǌ��K-W<������f�9�췹枷8��Ӌ�y7��&���t�����cO�-��5������H�7�iX̿���}O��WYq�5�� m���1k&���b�Z�\]r�7��]p�9W��]������   L���f������?  `����-N���b��&��=�P��N=��۟�[�G�p�e�}�%g�9��+���i�=*^|�o�����mq��&L���g���������{�|����c�{;̶�n���۟,>��Ӱ�as�)�}�mwګ��*.�������<��3�'���a1�_�k0(�ο`�ڤ��՟���y�x�*�\�G���8��u6�<gj��������&�}��[*������U&�gm���)n��٘�)'�uy1�<��<]��^����M{���^��5��?���=   �h�{�R}1  �7y���^w��Se��`�b�g��sB�����Y�'�z�Ә���������c�G�?��H9&]�~�x��s.��xx�G#��w����P�����_��8�;�[m�[���ߴx⥟�kA7=��/��nx������WY=����,��6RƩ>PPf���-�n�E��0Ysݍ��O���O�{l�so|1�#����}Z��男M�σŗ\�Xc��F��]��M?=��΍�c�\t�]q]��>�Ž|��5׏�]Q����+�y��1�Ժ��G�t�ȷy�����]sw��;��{����h���zb�   c��;K��4   ���b�CNyj�6��s�?�����Ų+�Z�>ǜ��ˮZh�ŋ6�f�?b�s�q��'_0�t�ꉳ�sƈ�F�*��랤���:�tƓj�i�F����#O)"WOL�~~WO�]����V]�Xd�%�Yg�=��Ӛ�	�K-���ϰͶ�i�wĎ��W������>�a��؝�� }��/�馛.�iͷ�B#O6�x��G��M{nRUfL��YUO{킕W[+fϪ�k1����b��g�{횪������ﳃV�._o�-G�k��n�O/Փ�?����n�m1�������V\��d���v;�8�س���|&^������#�[�:H'�}e�   S"�{��i   ��N�������G�|����b�yF���i~���犉km��0�r�]G�����T%��~1e��Q�΃v���Ǽ�A�oӪ�גּ�:1�0�>���k��=v�X�>�3��aQ����7�uvw<�j��TЩ2��_g��6������~zi�[@���n�mz��,���Ș�mU���|s  ��J�����'�  �n��
�|��0:U�/]ׯJ�ӊ�ox�Xq���u&����p�3q��h?<�����K�:w��v�?fvi�m9����I�����x*�k����q_�w��=�;b������{ō[�{1a8�EG�t~m����/���O̮�A%��&c)�w�V����K�w�t��zu��P��^b�ecΦ-�Т�EW�s  ��J����b   �e��Kn��vrinZs��K.�|�>]��*k���iOt�����?�L��|�ł/s���6���.<��b��f���d�g)N>������y��^*��W�Ə�)�U���ߊ��{g�����.�h��F����W�sF�m�
�=�~�8l�u<n�m��.q?��2q�x��K��ʇ����s�*�Zf����_r���.*���c   �tOZ�/�a   �����#�#����������R��V]��#_㝮S��˾�m��@�(,�W�˕[�s���Ax�*N���b��W�9i�u6*ξ�Ƒru��eUi0�*���������.V_{�x��m�շ:��I���[a�xm�V�����1�Xl�yw��LIa��O�nZ�
˓���m�߿ͷީ���'�y  ��ҽi����   �D�'��O�^��Ͼ�Xa��BUd������3LK������m&L�?3�P=M��N+���1�8���k{�m�9�~�����x�G�t~��ZԎݤ馛��j�]��{=�UE�����o�_�f]}HaX��0|~��/�K����s����V\-�{�ˌ3�4��=������g_�E<   4%ݯ��i   `J����%�Z���?�� �ً���������0��'��u���U7?\����b6�i������'U���"�-��Ӫ��S�w�	�uw<���Ͼ��ϪĜ^;%���?W����ﴥ�]�v���eI��k�-����b����~�����P^�݊K���x����I�w�{�aq��O\q[q��qr��A�G�t~q�w�/��������/��a<   �)�o:��b   �R�4��AL����7��lq�9W��}��S��g�xmGc�E������P�������&���z����sʅ�N�0��y�?���&��6{���[���n���{������_����b��W�G[�{�������uf�y�b��k;�~`q�Y���>VR�꽭��AG�2���:���Ya��\==�*�{�0   @�����T_L�    S�z��:lV����#��e�_%���{���+���g��*���j�8��顑�$�t�s#�Uq��:��\=�:�������{S<��#?�o�呟��W��p�3#?�{�㑯������'g������k�T�����ţ?�x佹���F
��=����U���mͨ�.R���{ꭑ"�y��Z�s�M#����.�>�f��+?k��   @3R7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr������           c��ɥ�b          �%u�K��4          �K�&��i          ���M.��7�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/��|1          �X�nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_��ߖ           c��ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr������          �I�&��i          ���M.�Z�          `,R7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&��?��r           `R7�T_L�           ��nr����          zI��R}�g�          �"u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��O�          �E�&��i          ���M.��0          @/��\�/~���           ��nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&��i          ���M.�?+_          0��\�/~���o          �1I��R}1          ���ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����?�           c��ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}���          �"u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}����           ��M.�Q�          `,R7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_��|1          �X�nr����          zI��R}1          ���ɥ�����           �I�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}���           c��ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��/��           ��M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��_�/          ��M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}1          ���ɥ��/_           �A�&��i          ���M.��0          @/��\�/��|1          �X�nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&�ꋿ���          `LR7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}1          ���ɥ���          �E�&��i          ���M.��0          @/��\�/����          �1H��R}1          ���ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&�ꋿ)_          0��\�/�a          �^R7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����.           � u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr������           c��ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}����          ��M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&�ꋿ���           ��nr���7�          �"u�K��4          �K�&��i          ���M.��0          @/��\�/�a          �^R7�T_L�           ��nr����          zI��R}1          ���ɥ�b          �%u�K��4          �K�&��i          ���M.��0          @/��\���s/P�}uA���R@t��.��
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
 ��E/k�� �)L(X"��f,Z�%Q�۬���K��,_�ֽ�`���� ��F�3�|��{�n�ֵ�yu=��V�  ��&'��( ���������k�X�����=���c��=��ÿ��{�oV�c6��9�:�_�_Ql6����|���~�;��eü��D�k�gI�glX�i�0�s��0�Y4���Y*C����^s�kI��6�s}l=�l�N9O��di/�D4�ܾ���Ǟ��߷~N�u?�{-jE}/����°��=�{��(����Q�=����>���ݱa��]ҸJ������}����8��}͹��.����]�KY��c,�+���x��-�ˢ�%Q�KT3/��,����}��C?��1��ט�kɞ�0��h�?��Q�c�۷�g�q�>K�=O{���K�'k�<6�۸�  Ģ��d� p��/	k��D�/�\5�|�D�%��h�g����~���l�ǚ��Q�����_��{��%������:�1^�ZQ���z�zƺ�<�����c�c���(��(O$����[�ۨ/[[ZÞ~�í���B�������J������=��Y�(OI��(�y�����y�S�?[�{�^�}�����Z�K�XK9�xFK�N�Z�Ϸ��)��#�mkߦ���ھ�����}n��Q�V���ܿ�qn5�����g�9.I�z�����ОsM�֚>F�zD�K���<��^%��z���Ǟ���{��x��i�[�:D�ܟٚ�϶�O��ձ����~  8��69�7F�  \���`����k^�F����Kr��<=�x��K�A�K�Z[�H������e?/y����>�⋈��Ls��s[����<�l�<�>Q�%Q�H�#ʽ�-�g��>�y֢=��L�i��m1�k����=�jL�վ��(Oɖ����V=�}���̖{|����:և�gk�<���Y��f�Oq�[��T4ǵ�8�|6J��s��Y��U�4�5��Y�֬�t�{�8��$k�똵{sl�s#�o����%[��gH��ni�s���[�wɰ����Vk����-��m�c�8�� ��E��ɼ1
 ��E/k��8���h�W��d�/&�^8O�Q���^Ζ�_(��ޫ����/���������u:�9ǥ��)�F���Ҳ���Ai~�Y?����Z���-}~�~.���s���s=�k�[�9&�#�c��1k�y��/{Zc��(�d�o{�#�>m%���y����|�J�~��7[��%=�^��d���r_����~X�|������oK�c��[^�i�����ﱴ�%����<K�q�����Y�Ec������������2������Y���Q?N��b,ǮY�;����c֌瘚~k2�L�&�{K��!�����G���j���9�+K���G��d�yX���=j>?��-���ve�c��=M�.��Ҿ+Y�'����ާf��+[�OM�k�  �+�MN�Q0  �/z!X���a�|��F��$�QR�b�dil��K/p׾�-�m)w����{��R����˅���s�O�gKs6ֻ>K����=f�5��3Ks��Fq�z��2�Q��(ϒ(O���l�/ZE}���x��Md8zb���9�>��>:��Z��^�}�z׿�gK��d��Ի��S�]�=�\ӯ���bƶ�����=V2�sLO�����_d:���{ߣǱ>��昼��cZ�ۚ{�����R�"�㣼K���7��i�z����&��玭m�4﹕ư�Gk�5�k�4ǹQ�Rߖ�V�+Ys�l�Xz�Uo��7����w�����X�<�ֺ�ߥ�Zk�7��i�  ����d� p����5Z_�����y��)9�R��ҋإ�5/r���ܵ�Zӿ�ڼ����b���5_.,��9�����Һ/Ŏ�쟖��|���%Q�z���l)~�mK�9�?�˖b�Jy���,��Ԋ�e��=����G}.�j,-�?�������`����Iw	�Κ��gz��~Y�:GyJ������^�Z�y�4�j�]��[�Y�s��ji�,�嘞��fm�"�y��[��k{l�-y[�nm�(o��<�k\����Y����ϙ(wd|�E�c���3WQ�sYZ�(&R�����wi���{�W��$�#�Sci��r/�gQ\�5}��ߠu����+ʻd�Y��ϳ�qm�gI�}�.�g  �%�MN�Q0  �/zX��h��K�b�f\5cȟ�b=�[�[sL��(ǒ��n��z��%��),��R���=�(iݟ������ݓ�f,��o7ō�<Sc5��R�(��9�r,���86K�Ĺ�>�Y�~���?K����E��3�lܷ5���K���=�zN�G�:���T��d͸j�0����ƶz���sl1����}�vﯙ��ةi"���>���s���<�Y�����=�<-J���oK���D9K���S�ֵ�ٳ�X���c(�Gɸo=j�[��z�g,�y.�������\�gkԬ[7ō���Q��ֹ�<�5��Ҟ��_V�ǥ9�b����Y��S=�m�X�sR��5k��������9��ܧ����� ���&'��( �������3�_��e|����^Q��-_8��l�r���ho}Q��?�3֚s�Ey�i��s��޶�sM��-ֺE����KzϒAͳV�G�楷�[�]�1�8&�_��W＞R�>���wKZσ����c\Qα��P3/[���Z���G��(WI���w\5�]���X;�{�9�x���(6R;�5k�;W=q��O�(�E�y�K�\�<Q��TgL�}��c���(G�(_�%�=c׼=�u����3%Rz������>Q|�(�9l�G����g����v���\�Y�)�=��\c�y���X�s���.嘊rDzβ�y�o��r��ާ%w�Z������:v�ֽ�m��  \��69�7��gR   �N�"�Ʒ?��a�������^��jE9Kz�W3�{/���H?�Ej_,׬a����u>s?�<c[�Q�z.y��V�.Y;�Ӻ�5kt��3к��(O�-�;��{�3����:��gU�/��,�֫u�KϾnY���}���(�����yD�SQ�9���k����ӳ�Q���qՎ%���<�k��z�s�YT�g[�4����\��;�i�-�Z�"��U37Q\�=����Gy�����y8&�W�3?��y�Rע�9���ݓ{�Ϙ(_�u?���(�9�1G�k�GѼD�+���(6���1Q����f}Z��OEq%[��f��j�ekެu�Gj��X�l�GL˜E�5z�g��F��{  ._T����`  �_�"�F��c��rs���Qޒ�/&j_G��D�S�s��O��e����\�sS�~߶u�ۖ_�E�%Q�cj���O={�u��\S-{c��[;Q�Tm�-��Dy�L�]��rN����멞g*�����wݶ~��(~��Y8�S���g���Y���g\Q�����'�ֳ�k΢���g!�_�{v�=Z��e.����Z��#��h=ײ(�1-c����>���(_I��x*׾k΅��Fq%�k�2g=���Ԉr�CԷ��y��{�9�Ӳ�k��(��(��u���y�='�ة��j�ٚ7kݯ��3�e�-�?���}��Ϡ�Q\I�9 �e�j��y�K҇ �}��5��ӟ9˕۾��es�_>F9{D�K�<�Ծt�b��*g��(���̢<�Ԏ=��+�q�8w�XǶZ��=r�v�kףe��(�1����~����=�m=�ٹ�0��j�זqG�K�<%ӹ��}�R����F^�(W��g*�rEz΃(ϒ-���9��:�k�����s&GyJZǵ����so=�c�\5Z�k�<���=k�3WcQά&owL��e>Z�y�e���#�sL�Y��֜={p�+Y{����ע�\�]�ֱdQ�%Q�cz�w��\�����^5{t0ͽ&�����#g_R�lj��(���sn��v��k�?ŗ���cZǓ�쿖q��,�sL��g~�ڱ��ͱ(  �+�MN�Q0  �/z	xj[�l��Q��r���s��ڼKs��?�3��eM[_T����i��V�&��lݗQ�cZ�O�9&�/�c��Dy"��"�1U����QI����2�Q��(�1���2�K��<f[�I�{i�e�{��)��G�ިͻ�Zla�<,Y���Q���qE9"�k嘺�3�g={����z���Q����,=�ܭ��(iݏ=��z��e�2OQ|IM�OyG�J.�̙��[rik�{.���ֱdQ�����LD��DyN�g�[�h��-q�����W�c�噉�KZr����g�6o�������Zρ�uoEj�g�vL���ZҺO�(�1-���OV�6-�3��  p����d� p������_�n�9ݯ���p�K��ε/�^4�������nY۞�9��Ex˚���%g��:�e�[�9�����v�Ƣ<�D�%���k��(G$���yh9���c��%Q��tZ��K�2��ֽwL�~�<S=c����<d�}o������`�3���DyJZ�ug����g�Z�j�����iG�\�|S9�xy�z�z|��gB�<���k���U_R�����:�cQ��5��[�ߒK[����~�cɢ<%���r�8�X�г�Y�7�_{����6� ��b#Ql��GKά�gЖ�QKk��E񑖜�ֽi�˃�g�F�>ɢ<%Q����iG����  p����d� p����k����B7�/;�x�\#�SI�K�(>R�"{��>���Dy�z^T��uϋ�A�<sr��yO-}�c�[ω���9&�/��=�,�9���<c��JK�,ʱ$�35�g���$-{n��\S��c-�T�آ<%Q�HK����G���z��}n��?�(OI˸�\�(G�gN�г��gQ��$�Q�:��5虫��?��>����g@�3:hٗQ|IM�OyG�JN�[D�-�����yhK�)�r�D9j�b,{�=��<暳i�|�9�rFF�%Q�H�>�bKZ��(>���e[�oiw�l��2�A�~-�y�j׳��i�'Y�=Z�l���Z��Ey  �^Qmr2o�� �~�K�[� �K�ߒڗ�-/�[_8���Ǣ�cj���e�����K�˟�m뼗�/Ǵ�Kk�s���(�1Q|I�>D9���K���>nٳ�c_:[Z�eQ�%Q������%iY�AϳYZ�H�9���Z��;k�ş���X����{��<%-�⏉�K���9�C�z��E-�� �S�3�(�1=��S�uԇ��3�g@�~iٗ���S�<��|%k��=]�Z��<��%����Eyj�=�����cKk>ܣ���/��K���}Ŗ�>�-{���ZrG�<���Ҿl����A�~-�y�j�3�-i�'Y��z�=�'��s ��D��ɼ1
 ��E/kl�x/QKj_���Pm}�?�,�}ˋ�u�yQݲ?Z�^�:��\�:�u�{���`:ǹ�w�{�ŗԞ'Y�|G9J��R�[�e�<��"�i}��<��<Y�>���$�����|N����e��[��кWZ���z��2Y�1���DyJj����(GI�#�v�ҳ�={�u�fQ����ғ?�r�)�CI�^왇=��A뾬ݓ��Ӻ�[�},�W��^�p�k����:�,�s�)�j��)�[�<��s*Һ�%Q�H�^�bKj��=�ִ��������5�io�2��9�sV;�(��g�Z�8k]�=�'k�� �zE��ɼ1
 ��E/kl�x/QKj_�F�Ǵ�pnyI_3�5/�kװ�Eu���r�h}	��s�u�ZK��s�qj�/:����=�:�Q��(G������l����kI��t�]¾_�u�e[��4�Ǵ��{�-�?���YK�/e���kψ��� �SR;����(�rD.��zֳg�Gy�DyJ�ޛ=�ϱ�Q?JZ�x�3��9�z�e{�k�27K����(_ɚ{�ʵ����C�X�(�1Q|ɚ��{,{�9���9��Z��(�rDj�B[�G�=��3��1k��-���{9���O�>��>׮�g  �/�MN�Q0  �/z	Xc��{��[��{iy	�?;}��j�/��%rֺ?z�s̚��;��/-�����r�%�{�ŗԞ'Y�|�)�_}�����%�גi�����.a�k�컭��������{�-�?��9�S������{��<%���ٟ{��wj=��sEyJz���{�'�9~VE�(i��)΀��g_f9nz�<'�}��k�foD�Jα{\�Z��<�E9JZ��(�d�\��%�s=�\�f{����MEq%�sŞC�y�>̚=���}�r����s޴�̓�����^�3ػ�  \��69�7F�  \��%`�K~Q����}��ڽ���n��X[�������]Z�H��뼇�/:�<�D�%-����{����<�<s[���d����3S�0�^�c��Y�9c[��=��������C��[���Q���q��מ���Z�z�>���{��ӽ�s���Vԏ�ֹ>�<��sd+kϒ�=>�򕬹�58�Z��<��Y������K��\��HϹ�G?zE����Jj�X���˱繵�k��-���n�Z֜�=�ɩ�ǖ��5?�S�  �+�MN�Q0  �/z	X�_F�-�-_L���{�G�KZ�|��\���}����#�sL_Ҳ�z����=�d�產>-b[������Vc�9c[Ι=�v�>DZ����1���DyJj�u�g����g��:�9z�i��劜j����u�O1={�u_������(_ɚ{]�KX�����1�>=}_3W=c��S�ǵZϮ�����MEq%5{�g�������s,�s*���۳�5g�g�����k~��'  .WT����`  �_���%�(��[R��D�K������":���Eu���yi��e-n�:׈�X���^z�L��(��eOD�%�>w[��,�sJQ��丞}t��Ϻn��z~���3{��5�m��i}n�����g�)�W[r�3��zֳuNz΄�}�3�,�9��_+�GI�\�bz�L�lu����������=>�򕬹�%�����yhKV�,��l�%�sN=�]�vͲ(��%��>���ֽ�f���y��<�r���=�T��-�OϹV���^�3�����  ��&'��( ����q�/
���ܖ/&Z�^@��k��z���yq��v=n�:���Xr	���/Q�c����9��KN}�>٩�8�iI�\�p�u��y����c��?"����9�{��<%��bK�}��g=[��L�٧=cɢ\�S���~���u�<��{�iU:c��*}.ҳ�Q��5��$��{?�c�j����l�\��%�sn=k^�vݢؒ���ZQ�Jj�X��]�o{�<_�\�s���y�9C��%=��9�j�?�k~=�?�> `{Qmr2o�� �~�K����0�o�_LdQ�K�׮�rx����E������z�][���Jj�to={%�sL_�2'Q|I�@[�y���!���.��a�����e��9���Q�ۤ��]�3c��3��Ԏ+�-���5�Y�ֳ���٧=cɢ\�S���~���u�<�>G={f�.�!���~�-��1={|�+Ys�Kp�k����:���Y�y���U�X�<��g�kԬ]WҲ�����f����(Ϟz�X����\��y�j�t��g<=�Z�3����� �u�j��yc ���^ָ��QKn�K�-��ֵ�E����K�^�߆u^#ߒ(�)�|��9&�/�9OQ|I�@[�yN�ǩ�O�z�;��r�zƹř����-{jϱ����&��o���{��<%��bKz�~Mzֳ�,�G�>��>���-D�(i��yh}���9�d����+�G�?�g��|%k�uNװ{?�c�j����l�\��%�sIzƴ$��XSҲ�����f����(Ϟz�xε���?w�:C��%=��9Z����  p���d� p����5.�Ea�ߒ�/&z^�^����X�9��g�/���9�K�o�9���z�G��(���<���%�z�<̢\�������٩�f��u���ߖ3fϱ��}�y�T�>���?�(OI�N=_���T�=z�}���b�����u�{��8��ȥ�4WK�S��>�+Ys�s�����yhKV�,��=�r��K�����W%K���{Q�Jj�۳w��Ѣw_D�N�\?w{ֲ������>{���T� ��j��yc ���^ָ��QK��b��e�j^�ֳg�[��ǔ���y=�\��쩧�Q�c�����ػ�k�<�����?K�9z�%��� �3��Vk��������5j}VZ��ȩ� �SR;�K?�O�g=k��A�=z����9�\m!�GI�\��C���ϑH͘���l���(_ɚ{�ڵ����C�X��g���Y��F�X�<��w.���1��y��SR�<���m�wO�>�[;����^ky����<�5?�=�  p���d� p����5.�Ea�ߒ���Ȣ\�R�����}��ѻS��^�:o��K�,�u*=}��ŗԜ'���������}����$���,�q�m�T_�Fz�T�~�sl��C��h}NZ�����g�)�ץ���ֳ����X���g���=-k{��Z+�GI�\���Y�޹�=��-={|�+Ys�S�Ƶ��yhKV�,��=�r��K����1RZÞ{���[��SR�<�z�������?w�z�f�N5�(ϒ(�1{��`ϟ1  \��69�7�4} ��'z	X#�(��]���%��_晊b��k�Z^.G�c=/��w^��-�����K.�yh���J��|�г�<�D�%��I�3�-����F����,����{٩קGQ�V=��r��=�(~�m��i=�>{�A���v\�p�R�z�<O�ا=k۲7O5WkE�(i��S�=k�羬��)�� �W��^�r�k����:���Y����k=c��\��=1U�#=��=�Vԟ�ھF�Kz�n��5�\{��W[�g��Rs����O����Y�k~{�� ��E��ɼ1
 ��E/k\�¨�%�/n{^�����Rԗ�%�p��UO�Ƣ��kZ���k�[�Y�(�1Q|I�\��w�k/=󛵜[���$ʓ��:���K?zƵ�z������{lw��C�u>Z�����g�)�WO��(�m�3-���\�٧=gO��S��ZQ?JZ�gZπ��l��{D�S�<��|%k�u
׼{?�c�Z��(~I��3�(�9}��혞�J{�g�eQ�S��RR�<��o˳���=�����r��g�~�]�(��g<Q���=���d={�� ��D��ɼ1
 ��E/k\�¨�%�_L��H϶���'�u)_ˋߚ���D�%-c�׮�X�K�,�5��u�K4��_�l�gDy���KZ�q�cIϳrL�UZ���!���-�ǒ(Ϡge�z~k����[����l�{l��Ö���V��j}>֞]{�A���e\Q����6�s���zֳg���'�S�z6����j��%������sη�u�j�'.�?�u�Ǣ|%kk_����ֱd-�Bϳ��<z��9���-s���li���%k��T�U�������g�f�kxLΓ�Z����,ʷ���n1�=�H��kSmޱ(OI�=���l8�Zl��  \��69�7F�  \��%`�K~Q��������Ym�%���K��>F9�z�ݲ?��[�=/��<���f�\�=��� ʷ���*�/i�����L��z��Y��cG�X���<�C^���%[췬g>[���c;�~΋s�|�j=������g�)iW��x�������g�<D�J�%�gO�<�\�����L虇ֹn]ˬe�[΀ھ��+��>�+Ys��]�Z��<��%kyN�<z��9���={e�9�ˬeo3O��S��ר��y�z�02�Ai�z���>g����~�l��j��Mk�9��D9J�����,�b� �Qmr2o|�R   �N��ƽ�A�K����A�H_c�|/s�����7�1��u�C\�Z���1�1��8�:�g�[�;߿e[�|�9&�/i��o���_�ek�z��K��3ǃ��<&�Q����$�3v�5�K�x�X�l��j���)�娱v����U�ٰ��{�A���e\�:O�s�v/n�g={��z�(GI��(9�\����u_��C���ϑ�����u��9O�|%kϮ�܆���yhKV�,�K�<Kz��9���:����/�qZ����x������5Z����]����5}mY�����jϛ�=�E��{�V혲(����Yk�[�g{�O���b� �Qmr2o�� �~�K����0�oI�=/��s6��R���jD9�z^T��uM���w�z_�:��\_n��!O��Fy���KZΓA��FϽ-9���jM�a}k���޿F�gj�G�s?ˑ��m�q��ڵ�N1��{z�q|�K�S�����`��3�򔴎+�Qc��m�ck=�f�G��i�O˳г���^Q?JZ�gZπ=����{���֜Q��5���mX����ֱd��N�:d=g[�X�<�0����ҲWj�7���f���߻FK;������59Z�}j�Ze9>�9sNuL��O�x-�j��־G9��=?�i�ź p9���d� p����5.�Ea�ߒ���=/m�Z�m|��~�����)^$�c[_֏��S3�K^�S8����>����(�1Q|I����a��~�ku�>N�`���8�V�'�ֈr�S˗���}�׸f/M׵F�y��ئz�0��~=�ҞZυ���{��<%��j����9߯u��g=����~���r�/9�\���QҺ�{�����Hϳ�뙓�sd�+Y�W������k����3G���`�1d=�Ȣ\�0=�Z�xP{�E�S�s������Z�?�o�����ݱ5�JK?O�V�y�>3��������L�Z����r��2g�}�=?={siM �.Qmr2o�� �~�K����0�oI����c�s7}�}f��rv�?���%r��E~�K������������ �Gy����6n���(�1Q|Ik�Q�Z-s=�kyΦ�=z��ڸlW+�靃�=���q��ʰ^5����SӽS�ez�3V{���}朦�[��Y�{��<%=��Ԫ��t�z�bO=������ֹ�rDz��s�#�GI��3�s�?�)�����7V�����|5z���9�~�k����3G=��3wQ�H��|��O��%5sQ��SQ|���L�ܲ�ql���gk����y�>��g�Y�j�9��ڸ�>���Џ�������{D�����==���L�u��= �u�j��yc ���^�h}Q{JQ�DyJ���N_����^>��w���8��>����X���e��K�-��*�s�u>�i�z��9�~v�4�F�'��ŗ�~)�E�|���:z6Ɵ��U�>�s������l��ql�(�1�<�h�ޢ>��^��3�5����.�Y��~G�\��؎=�gǟ����9�{��<K�<%=c���9_�?Y�6��]���i��}Z��Gj�r��?�\����$�s�)�v=�j��zfr\��'߯7�`ܯ=�3���0�5�F����{5�ױ�^�*��3=��X�X��s���j�#�;�g�M�qL�d�g���jL�kDyJ���L�[��h���Uc���f�{��t�Z��8�Ʊ��?�g�}�ꙗ����rd��c{��{~j�0�f<  \��69�7F�  \����V=/jO�wLQ��5sת�e�^}�͛�~�1��L�9j-��n��;���|J�ֺ>Y�gI�'�3�5븴��3��)��3�Q�%Q��(G��9�C�:��#�e��v�n��=�C���=���K�<�z�O����k�kI�~�uI�Ϡg��������sV��k����-���9���\����|ιv^r|��=߻������>k[�k�s����Q�Z��ҽz�M��圎�sk���g�={�W����Q��-�Z�5{�ku��4�O^�q��>�}��Ռ���Y����`�������s  ���&'��( ����rP�2��z�����_N��s�c�a^�kU�����ڗ�K����|����N1���?���[�^�=O�\�z��i�3�c��c�3��~��l��oe�"����;Q�H��u��z�\�K���;�(W�(_�(�1��F�}g�^z�z����5��}֮A�\�9�z�����gH�/����9:fأ[�my���^[=ӭn�Z��E�K�\�(vI��O�����=͓��w�r\�_{�����-��M-�AS�wn[����s��m%�st��X��~nqv�լc�|D����}���D�K��� ��D��ɼ1
 �zm�R}���[��Y�w,�/qk��i��L�6��X��/�g�5�%K�[cm�Na�/��j�����K͜�9#�|����\�;�c9�%>{�k��la�5�Ϛy��ǚ�D�jm5��ֽyJ���Ӛ��r&G�5�1�����}Ԩ�ݦ��5c��[�7�]����p�=�;���k�05��yk�(�j�k��k_�5��ɿ�Y���Zs��0֭�����j�z�_K9֌��k��-�U��\�z���i*��%�s����5��wO�޳��f���/{�Ϛ��� p����d���W�*  W�_=�Y�����E� �̹�{�JM���</y�cɟ���/���-��ƿf�5c�ݓ�s�Yx���.�m�\�]�:���Ϣ�V�?���c��4���<���Z;�H�\/٣�9gt�c���jQ�^[�l����ڱԬ�V�T���ܾ՚D�kl��j����8�k[�3����Y��V��(�G�{x��{�f˹���"��E��C�����Y�7r�-��-����[��<l1��f�׎�ؾZ�����c�\-zk\���k��۫��T��1[�?�������%���ط�Vs<�ӷ-��<�|�����G���}���k�U?��m��r�M�����O���3tl~jǱ�X~  �KT����`  .�/{-�(�Jt�-��?�,���k��-�G�ܗ(����Kߞ�T;�[���ڹ���u>�s��g��}ﳳ��\�l5ޥg{�s�q��ploF�]+�;��)�>���i,�[k<g�}/�]����O�ܞR�ǭM����7��-�죭Ɯ�י�cϵ�yf�|v��kM�(�V����k]�{w�u�#��-�#��Z�WQ̒5gE��r.��3��"�Y+�ȿ��6�ߢwo�s�3���>�~];��m=��~.��-ʳ�ܗ�~��~�(W��sf�s�Ւ��+��V������u�>�Dy��k��bl�q p>Qmr2o�� �l{��]r��ѽ����5_P���>�����_2�쩖/.�}�����g˘z���Ԗ�|.���k1��?�_k�^�w^����u�>�����k��G�:���lY�(~����=��k�$�Mo���s�����=�z���<G���5��w����:�[칖�-��<wk����~dQ��D��׺��p�Y�d|��߷0�?�v�Z�U�ϟ��(G��{f[=#[����b�߻ekՎ-��vm��W[��,�k|6\�ܧq_���c�ck=��r�z� ʷ����ȱQ�%{���9֪FK�j��;�k�km�z�g$ʿָ��Cɖ� �����d�  ����؋��>|&ʱ�c����%ñ��9s\�/����=v������u>�a�K�?��V�>*���s���L��ð'�x����k�9���k��ϱ=4|Ƴ����eٰ.�sQX2쳭�V)�^�u����mZ����<��導�gs��s�jg�o�0?Y�ܱ9Z������+>smkRkim�s�����R�~F�{9֧�/�e���a�Y�  �IT����`          ���69�7�D�0          @��69�7F�           %Qmr2o��          J���d�          �D��ɼ1
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
          (�j��yc          P�&'�?���W����           �\�<�KNr���>j���~C�           �k��u�I�U~����{^��׆I           "�yZ���Z�^;j��/���פ           ��yZ���Z�^5j��y/|��kސ          *��i]r�k�x�����ȋ�$           �\�<�KNr��?>j����y~�           �k��u�I�U~�GG�<�?&          ���i]r�k���Q�=O���&          ���i]r�k���Q�=O��g�I           "�yZ���Z��1j��[����I           "�yZ���Z��>j������$           �\�<�KNr��O5����[�$           �\�<�KNr���|�p��|�S�$           �\�<�KNr��_3j��K��W�I           "�yZ���Z��d�p����/�           D>��>��&��\���g����O��a          �H�A��%'��<��Q�=�͟���$           �\�<�KNr��1j��?���$           ���;�I�/�*?𾣆{��{<��)          �F�A��%'(y��G�����=a          �H�A��%'�V�����y�~�0	          @$� O뒓\�|��O��?���~)L          0�k����I�Q~�����~�߽<L          0�k����I�Q~��|ﳟ&          ��g=����\����Y����i�g�          `,�O둓\����K>�u���7����           P�k����I�Q~�����O��'��           �r��9�5�o�>'9��'~ʧ��           �r��9�5�o��tr������          ����i=r�k��|�wr��w}�G��           �r��9�5�o�~[r���z����           �r��9�5�כ����K^&          �r��9ɵɳ������zv�           �5��:�$�&ϮoO>��_��0)          @�k��u�I�M�]_�|�o|�¤           Y�9��!'�6yv}fr�����O�y]J          �5��:�$�&Ϯ?�|�C>��ä           Y�9��!'�6yv=:9�໼���          d��xZ�����������o���I          �\s<�CNrmrx�.9������          �ݖk����I�I>z}Or�u���ar          �n˵����$�$��29�����ar          �n˵����$�$�>99x�G~t�          ��r���8�5�G��M~��&          �\k<�?NrM��뭒Y����           wS�1�j��\�\�^�}�w=3�	          p7��i�q�k���~��{Jx          �n�5�Ӻ�$�"/^_�~ҧ|Zx          �n�5�Ӻ�$�"/^*9|�����ޘ          $��xZw��Z���]����=��           wS�1��'������ �߾��          ��%�O덓\�\}�`r�����          wK�-��'������ �_��Ͻ���          ����g~�A��}������ �c��Û          wK�-��'����/��o�oq�_����          �ݐk�sm��8�5�M׫��$��]��          �rM��8ɵ��׿H=��          ���E_zPc|_�=n�>+9H��O�7          �\S<�3Nr�q��~�A���;gxS          �n�5��:�$�w]��$��� �1          p��Z�i}q�k����M�������          ����_�j���5����L>��?)�9          p�}��}�Am�}������ ỽۣ          �ۻ��#j���5���oMfI_���;           �N��8�-Nr����Aү���v          ��r�8ɵƫ��J��}b�	          �v����AM�}��x����A�w|�w
;          �N��xZW��Z��׃���$�3�v          �]r��8�5ƹ�x�뻒�|��$�          p����i=q�k�7�>;9�����y}�9          p����i=q�k�7�~2��O���          �C��j��\c������&���o�y��RG          �[)�O눓\[���5���>�>9�          p;��iq�k�7��Trp���{~o�)          �v�5��:�$�o~=4���Y?��a�          ��k���$��r����fO|ғ��          �-�
O뇓\S��������h�9          ��Z�i�p�k�w���n��y�A          �:��v8�5Ż^/Jn��_��R��#          pK<�˾�f��\K���������$�$          p�r��n8ɵĻ_�If7ޏ�D�Q          �<��?9��/���zArp�'}闇�          �K���'���d�����c? �,          p]��xP+|���]0�u⇟���          �u�5�Q�p�k�Oz=/9�����a�         ����'�����~$9����AG��}�7�4          p����A��}�v��ף�Yg�����SG         ���C������˵�g���t���%a�         �˖k����I�>��W��=��;          \��>��j���5�g��5�u�;����           ������Y]�}�f���3��N}�'~j8          �2��i]p�k��~}Rrб?�!7?��7�          .K��}�C��&��\+|��-�_N:����?���_N           .Z����'�F8�
_���䠃��~�          �,���8��/�_���ɬ��x��          .C���j��\#|Q�s��N�����          �ːk~�u�I����ӓ�����>��u����          ��ʵ�{��� ߗk�/�z��?$�}��|}88          �r���7�5��6�"��K:��?���          ��A�Cj��5�{}P2�����#�          ���5�Q�o�k�/�z^r��O���          8�\�;��Mr-��_��t�Az��+^��p�          �i���=��5���Z��������_��/	          �V����&�8�_���� �����ip          �y��;��Z��r��\�/��~�?          ��W��4��/� _��ϒ�A<�1�          8�Ǽ�{�������]8��[��i��         �}�Zި�7���*�g$��������          �ʵ����$��^���&�A=��~ �           `��7��Mr��U_?����Ǉ�           �#��N�z�\�{�ק'����E/'          �V�ݍjz�\�{+��$��+��9�d           �ʵ��z�$��ޚ�	�� ��-���/{u8!          �6r�n�ݝ��&����\�#���� ��_��pR          �m��io�k{s�ﭺ�N2�_��pb          �ur�nTÛ���[w�}򟒃���O�n��+iB          �M�Z�i�n�kzsmﭼ������4�           �O�эjw�\�{k�ߚ�Zr0�O��O'	          �kt�u�I���5�����d6����o~)M          �N�͍jv�\�{믷Iޘ���I�d          mrm�^7�5����N\���&�����	          ��ܨV7�5�w�zp��`�ܟ|8i          @�?��>��F��\��kx���?%������7��          e�7��Mr����g����������          e�wZ����\�{'��Jf��M��m�          �\���&�f�N_?�L�|�G��          �r�.7ɵ�w���d69_���N8�          ��\{��&�Vו��N&�m��a7/�χ
          ��\s�ko���I��uݿ�?�M�_���'          �M��6��Mr��kt}m2���}�s         ��.��F5�I��uM�wH~=9��?���'          �\k;��MrMn��u��%�I�_����	�O          �}��6��MrM��p�(9��G>���I         ��*��N�n�\��Z��L2��/zғÉ         ��&��F5�I��uU\ߑL�o�-�������	         ��"�����i�m�kp]���$�I���'nޘ&          �\S��&���p��d6�O��'          n�\K��&����x=(yEr0���������υ           �U��͵����$����[W����٤~��?�捿�&          �\C��&��ֵ�zj2��o��o          n�\;��&�$q���>yCr0�����r��7�F�           p[��\;;��Mr�m��ump}j2����K�.
          ��f6��Mr��k��_%���?���Å         �k�ke��$�ֺ6��-����~��~�pq          ���Z�i�l�kjsm�k��s�٤�����         ���5�Q�l�kj];^�Nf������B         ��ɵ�Q�l�ki];_�Mf������n^���.          \�\�kc���$�ҺNp}Q2[����O          �E���je�\C�:���d�_��q�p          p�r-lT#�|O�:����?$���o�67��K�         �K�k`s-�>6�5��v�u��S�٢|��c�"         ���5�Qml�kf]g��I2[�/x��	          �&׾F5�I��u��z����l���]�s󦴀          p�r�kT���\+뺀�q�l�����o����         ���Z�w��`V{_��u]��7��B}���pq         ��>�S?mV�zߗ&�������~�W�          �k\����ى�B�?���d�pO���	          N-׶F5�I���5����b2[��?�7�x�χ          ��kZ�w�ջޗka]Wp}u2[��}�          N�#>�fu���X�]�Jf�?~�_          ��kY��$׾���z����lA��7���7�ZZt          8�\�ն&��5׾�����lQ���|��p�          `k�v5װF���G'�+��z2[��x�cn~�~-�          ��\��kW����%�[p}S2[�?��>&�          ��\�ղ&���uK�%�Of��O|R�1          `�\�հ&��5׸�n����Lf��_���         �^�F5�]MrMk�mu����I��7��ۿ3�(          �*צF5���V�-��V2[��{��z��>?�0          P+פ��Ԩf5ɵ��;p=5�m�G>��o~�U�	7          ,ɵ��&5�UMr��]�Hf�C�����Y          �U�E�jT�\��c�oO^��6ğ��Ǉ          ������Y]�}?���U���#��d�1>�s>/�H          0�kO���$ת�U���Xm��/��_n(          �Ө���%.���D��~�?7          |�?�ְ��\��r����$�(7��i�n0          �\c՞ޗkS]�����l�<���y�wO��          �{rm�C�Y��}�&��:z=-�m��=��n���?n8          ��{�ޫ-�jN�\��r��$�Jf����n��C�7          �_�%}�w�]�:��rj�Eu���K���6��}��o��c/7           �W�!}�#1�/�/מ�T����]ɏ%���G>��%?���F         ��ɵ���4�-Mr�i�=u���wN~2�m�G?�=o~�g~>ܐ          ��f�1��^�z��r��#���z��U�l�=��p�׽1ܘ          \�\+�kF�Z�$ט�ZS�k������F�����7��7�
         ���5��V4�!Mrmi�1u�6�>0��d���F|�/�����         �ȵ��b�\S�kK]�ͯO��d����p�S�zM�a         ��&4׆F5�I�%�5�.�n�G'��y�{���K~�Í         ��˵��~�{����}T�r�~��$ڀ7�|�n~��^r�+i�         p=rh��jD��5�.�ɮ\��_z�6����<���F    ��ڻ����3��v�fӌ0��3Z�E��:�RV���6����|�ʲ�Xm7�J$�e����f�l1�J(qj�%A��lXb��q�Ծ���so�I.q�;w���5�������     �?�~���Ku�]���|��a���~��>���         @���;�鴟U��f�l�m���tԨ5�Ov[��         `���3�Ϫوf+j�ҷu����>������?p          V�l<���Аmh6�f�fc����KtU�C         ��e�Y5�]�	�6Ԭ�m��`�~���3�.�          ��l:�ֳK���~������c�;���         ��e�Y5�]��Ԭ�o�pG�~�m�}�/           +N6�U��%��l@�̆�kC��n���.m�s���          @��f3�ͪ���g��fr���F��~q���a          ��e���f�rv���l�oB�~�m�5�j�^c�@          X~�hf�Y5�]��44;(T?�)�[�P          x��ͬ��n��4t�=��};�S�         @�e�Y��]���-��m	�h_��>�����         ���3[̪��2;d�i6���05T�m�Ŗ�{Y>$          ���e6�U��%��l8͆�.	ՃhÆk?����A         ��l.�����.�l�ٝ���q�q'��          X�q��l0��V�l�o��V�I��/�9s_,         �P�me6�U{�%��l4ͬk[�B�`�n�n������xd          CX6��VV�e�l2��4�%6,\���q���         `(Ȗ�j,��3�L3[ƾ��qȡ���-x�|�          �Q���PVme7�	f�Í�	�cjc�|��ϙW         `0�f2�ɪ��������6w��au��݉��         ����n������>�����u���C��()         �@�md6�U;�M6�f�K��ck9�]p�E�         H���6�j&�dS�m�����,L���8������_+/         @�d��U#��!�J3[���ر�f��n��|�          �Q���@Vmd7�����vO��1v���         Пd�X���d3����F��C�0;�����ƛ���         `e��1[Ǫ��&[�l&�l%���f�iǑ�pT����c         �K�4f�X5��d������}2��۱���˯��|�          }![�l�ֱ��!�H3�;)Tw�}�?�=;���          ��.f�X��K�������B��;�^�#m��?,         @o�f1�Ūi�&��l �l ���z�����|��Z         �ޔ�b��U��M���<�� �f��P=��UW]��2���P          ,�l�Q���n�q�������c�=zL���+ʣ         ��"f�X���dӘm����.��_d�ؿjw�=�<"          �l�A���%\�i4�A���_��,r��#����G          ek��a�".!��l�l�lX8;Ta�5��N���         ��l�5��%d����}:�����oѮ����          CK6��V���Q�V��l���Cu,�c���]��(�         0�eC�-a�.!��l�̖�	aa���"{�_��˙�7�          �Z6��VM��A���l�[7��C��<�=0kvy�         ��-�l�����a6�ff=���g�:*�9��o��}�<V         ���M`��U3X��0�C3�����P�Ō?����Sϔ�         �߲��j�fchf�k;:<������>�=�ϗ�         �_�����j�fShf�B�G��P�EV_}�v̱ǷG�x�<n         �ʕ�_�~��U-��̆0[B3���	a~���b��o�v�]���         �[��e�W5�l��v�̬�72�^ՑZ�n��mW���q��          ��l����B���f+hf�ҷV�^�}���b���/�y�_XC         �we���^����&0�@3�~�u����P�Ŭ�����'����A         zG�y��e�W5|�l����̬��O��P��6�y�Qm�}3ˣ	         �L�x��e�W5{�l�������l@���`���Rv����_/�����[�!         ��]�w��Um޻ȶ/?3�A��SCu���ƨQm����{zy\        `���.[�l��]d��7��l�n�pY��`i��o�;���ʂ��        �P�-]6u��U��2�(|6���m�o��0.e�����C�ަ�6�<�         0Xe;�]�tUc�.���V/�=3�!�Q�p��eiĈ�䉓ڬ)3         t3g��r��U-�2d��m^6zff�m;�Û�:��m��L�4������˃         E�p��m�Ͷe3���e��-�����>��꨾����v��Oo�>�ty�        ��y䱧:�[6pU�����.�;33[�m���C��v�q�v�g�'�z�<�         ��d������nU���˶.;33륭
����.�Ο�B��������         V�lزe˦�j�z �l鲩33��M÷P�eZ��A�.���6�Wʏ         |P٨e���Z�kU����e3�휙����¯Bu���N�۹M<er�����         �T�h٤e�V5k=t�6.933�G��pb�����[��m�j�\zy�3���c         �/[�lβ=��j�z(۷�l���l l�01��a��>9�p��m��f=��         ��lɲ)˶,��={�q��-�733��"n	��ﱏ~t�6n�=ۤɧ��?��^]�V�A        `��F,[�lƲˆ�j�ާlٲi˶����F���T�>=��j��v�l;�o����q{��ʏ         �_6`قe�mX6bU;�>e���Z�kٰ���ۧ�1��v�>�˨5��|����vι��n�����W���1        `�˦+ۮl�����+ۯ�	[٢e��mZ6jfff���0.�f��C��6�p���q{�ogB���+��.?�         ��l����v+�l������,۳lвE333���
_
��m��P}h�۰a��&�l�v�e�v�!_o'O��.��v����3s�+?�         �#[�l����+[�l����F�j�>�lɲ)˶,�l����zm�	G�υ�c�kF��F��v�}l�;0�{�u���ƛni�g�ߞ������ʏ0        �@�mT6R3�{��Le;�U�T�Te[��U�^��lŲ�v,233�>ݟ������j�>X+�ȑ#�'>�~�r�����K���~�>��6�Nlg��/����kڵ���~�ӛ�-��ަ�uo���_��g>�~�����϶9s���8����+���*�#         ��l��9��(�l��I�6)�l��Y�v)�l��iʶ)�l��y�9�l����Z�>�X�`لe�����Y���a���!���P}�         Xy����+[�l���233���m88�n/��#        @��V+��l���ʖ+�.33�!�Qa��g8*��0+�.TO         ޑ�U6W�^e��-V6Y�fe�efff���a��p\87\����C�        Ȳ��F*[�l���ʆ*[�l����������`�W@�m�ذoȿ:%L	W�kÍaj�3L3�#��07������P�G         @w�es��Q6H�"e��mR6J�*e���R6L�2eӔmS6NG�l��}�6d�_F6[e�� ��48�f�    IEND�B`����� JFIF JJ  �� C 


�� C		�� ,," ��           	
�� �   } !1AQa"q2���#B��R��$3br�	
%&'()*456789:CDEFGHIJSTUVWXYZcdefghijstuvwxyz���������������������������������������������������������������������������        	
�� �  w !1AQaq"2�B����	#3R�br�
$4�%�&'()*56789:CDEFGHIJSTUVWXYZcdefghijstuvwxyz��������������������������������������������������������������������������   ? �S��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ������M�����+[XT��L�u%� }hz+�o�_�S/����s��L5X�>��t�0�fȄs��}��o���ZOx���>�L�#����56}B�
��Eob}h�
����)��G��-�����?�iq�Fn}�
�!���K_���%��S����Z���,v����8�!��� �է���_����V���Zr}�����k�Z�>�ʵ8�϶�u� ���Y%�Gox�E�d���&"}��`=�/�pod�#�»x?�s�j�&~�G���5���'g����Oĳ/%�;�Q��a��=��� ��o�݇F�6�e*��l#i��߭v�-��&��,\�����
��N|@f_�M�������$�֗&�n>2~�?����?#ɱ�I���ѥU�UUB�� JZ�Y�#'��D~j�����g���O�-�_x��� ����U���5�����#u�f�4����Mh��[Eud}j�����N��������;��P�o�;���!� ��_��_��Ug��'�U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~a� ð~)� ��� �u� ��ð~)� ��� �u� ���yE��;0��C����O���� �˯�F����O���� �˯�F���(�ϡه֪�������]�5������]�5~�QG�}�>�P��� �`�S� �� �� �2�� ��� �`�S� �� �� �2�� ����?��va�����;���� ��_��G�;���� ��_��_��Q��C��T?0� ��?� �?�� ��� �j?��?� �?�� ��� �j�<�����}j���� ���� A� � �e�� #Q� ���� A� � �e�� #W��g����U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~a� ð~)� ��� �u� ��ð~)� ��� �u� ���yE��;0��C����O���� �˯�F����O���� �˯�F���(�ϡه֪�������]�5������]�5~�QG�}�>�P��� �`�S� �� �� �2�� ��� �`�S� �� �� �2�� ����?��va�����;���� ��_��G�;���� ��_��_��Q��C��T?0� ��?� �?�� ��� �j?��?� �?�� ��� �j�<�����}j���� ���� A� � �e�� #Q� ���� A� � �e�� #W��g����U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~a� ð~)� ��� �u� ��ð~)� ��� �u� ���yE��;0��C����O���� �˯�F����O���� �˯�F���(�ϡه֪�������]�5������]�5~�QG�}�>�P��� �`�S� �� �� �2�� ��� �`�S� �� �� �2�� ����?��va�����;���� ��_��G�;���� ��_��_��Q��C��T?0� ��?� �?�� ��� �j?��?� �?�� ��� �j�<�����}j���� ���� A� � �e�� #Q� ���� A� � �e�� #W��g����U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~a� ð~)� ��� �u� ��ð~)� ��� �u� ���yE��;0��C����O���� �˯�F����O���� �˯�F���(�ϡه֪�������]�5������]�5~�QG�}�>�P��� �`�S� �� �� �2�� ��� �`�S� �� �� �2�� ����?��va�����;���� ��_��G�;���� ��_��_��Q��C��T?0� ��?� �?�� ��� �j?��?� �?�� ��� �j�<�����}j���� ���� A� � �e�� #Q� ���� A� � �e�� #W��g����U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~a� ð~)� ��� �u� ��ð~)� ��� �u� ���yE��;0��C����O���� �˯�F����O���� �˯�F���(�ϡه֪�������]�5������]�5~�QG�}�>�P��� �`�S� �� �� �2�� ��� �`�S� �� �� �2�� ����?��va�����;���� ��_��G�;���� ��_��_��Q��C��T?0� ��?� �?�� ��� �j?��?� �?�� ��� �j�<�����}j���� ���� A� � �e�� #Q� ���� A� � �e�� #W��g����U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~a� ð~)� ��� �u� ��ð~)� ��� �u� ���yE��;0��C����O���� �˯�F����O���� �˯�F���(�ϡه֪�������]�5������]�5~�QG�}�>�P��� �`�S� �� �� �2�� ��� �`�S� �� �� �2�� ����?��va�����;���� ��_��G�;���� ��_��_��Q��C��T?0� ��?� �?�� ��� �j?��?� �?�� ��� �j�<�����}j���� ���� A� � �e�� #Q� ���� A� � �e�� #W��g����U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~a� ð~)� ��� �u� ��ð~)� ��� �u� ���yE��;0��C����O���� �˯�F����O���� �˯�F���(�ϡه֪�������]�5������]�5~�QG�}�>�P��� �`�S� �� �� �2�� ��� �`�S� �� �� �2�� ����?��va�����;���� ��_��G�;���� ��_��_��Q��C��T?0� ��?� �?�� ��� �j?��?� �?�� ��� �j�<�����}j���� ���� A� � �e�� #Q� ���� A� � �e�� #W��g����U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~a� ð~)� ��� �u� ��ð~)� ��� �u� ���yE��;0��C����O���� �˯�F����O���� �˯�F���(�ϡه֪�������]�5������]�5~�QG�}�>�P��� �`�S� �� �� �2�� ��� �`�S� �� �� �2�� ����?��va�����;���� ��_��G�;���� ��_��_��Q��C��T?0� ��?� �?�� ��� �j?��?� �?�� ��� �j�<�����}j���� ���� A� � �e�� #Q� ���� A� � �e�� #W��g����U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~a� ð~)� ��� �u� ��ð~)� ��� �u� ���yE��;0��C����O���� �˯�F����O���� �˯�F���(�ϡه֪�������]�5������]�5~�QG�}�>�P��� �`�S� �� �� �2�� ��� �`�S� �� �� �2�� ����?��va�����;���� ��_��G�;���� ��_��_��Q��C��T?0� ��?� �?�� ��� �j?��?� �?�� ��� �j�<�����}j���� ���� A� � �e�� #Q� ���� A� � �e�� #W��g����U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~a� ð~)� ��� �u� ��ð~)� ��� �u� ���yE��;0��C����O���� �˯�F����O���� �˯�F���(�ϡه֪�������]�5������]�5~�QG�}�>�P��� �`�S� �� �� �2�� ��� �`�S� �� �� �2�� ����?��va�����;���� ��_��G�;���� ��_��_��Q��C��T?0� ��?� �?�� ��� �j?��?� �?�� ��� �j�<�����}j���� ���� A� � �e�� #Q� ���� A� � �e�� #W��g����U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~a� ð~)� ��� �u� ��ð~)� ��� �u� ���yE��;0��C����O���� �˯�F����O���� �˯�F���(�ϡه֪�������]�5������]�5~�QG�}�>�P��� �`�S� �� �� �2�� ��� �`�S� �� �� �2�� ����?��va�����;���� ��_��G�;���� ��_��_��Q��C��T?0� ��?� �?�� ��� �j?��?� �?�� ��� �j�<�����}j���� ���� A� � �e�� #Q� ���� A� � �e�� #W��g����U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~a� ð~)� ��� �u� ��ð~)� ��� �u� ���yE��;0��C����O���� �˯�F����O���� �˯�F���(�ϡه֪�������]�5������]�5~�QG�}�>�P��� �`�S� �� �� �2�� ��� �`�S� �� �� �2�� ����?��va�����;���� ��_��G�;���� ��_��_��Q��C��T?0� ��?� �?�� ��� �j?��?� �?�� ��� �j�<�����}j���� ���� A� � �e�� #Q� ���� A� � �e�� #W��g����U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~a� ð~)� ��� �u� ��ð~)� ��� �u� ���yE��;0��C����O���� �˯�F����O���� �˯�F���(�ϡه֪�������]�5������]�5~�QG�}�>�P��� �`�S� �� �� �2�� ��� �`�S� �� �� �2�� ����?��va�����;���� ��_��G�;���� ��_��_��Q��C��T?0� ��?� �?�� ��� �j?��?� �?�� ��� �j�<�����}j���� ���� A� � �e�� #Q� ���� A� � �e�� #W��g����U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~a� ð~)� ��� �u� ��ð~)� ��� �u� ���yE��;0��C����O���� �˯�F����O���� �˯�F���(�ϡه֪�������]�5������]�5~�QG�}�>�P��� �`�S� �� �� �2�� ��� �`�S� �� �� �2�� ����?��va�����;���� ��_��G�;���� ��_��_��Q��C��T?0� ��?� �?�� ��� �j?��?� �?�� ��� �j�<�����}j���� ���� A� � �e�� #Q� ���� A� � �e�� #W��g����U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~a� ð~)� ��� �u� ��ð~)� ��� �u� ���yE��;0��C����O���� �˯�F����O���� �˯�F���(�ϡه֪�������]�5������]�5~�QG�}�>�P��� �`�S� �� �� �2�� ��� �`�S� �� �� �2�� ����?��va�����;���� ��_��G�;���� ��_��_��Q��C��T?0� ��?� �?�� ��� �j?��?� �?�� ��� �j�<�����}j���� ���� A� � �e�� #Q� ���� A� � �e�� #W��g����U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~a� ð~)� ��� �u� ��ð~)� ��� �u� ���yE��;0��C����O���� �˯�F����O���� �˯�F���(�ϡه֪�������]�5������]�5~�QG�}�>�P��� �`�S� �� �� �2�� ��� �`�S� �� �� �2�� ����?��va�����;���� ��_��G�;���� ��_��_��Q��C��T?0� ��?� �?�� ��� �j?��?� �?�� ��� �j�<�����}j���� ���� A� � �e�� #Q� ���� A� � �e�� #W��g����U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~a� ð~)� ��� �u� ��ð~)� ��� �u� ���yE��;0��C����O���� �˯�F����O���� �˯�F���(�ϡه֪�������]�5������]�5~�QG�}�>�P��� �`�S� �� �� �2�� ��� �`�S� �� �� �2�� ����?��va�����;���� ��_��G�;���� ��_��_��Q��C��T?0� ��?� �?�� ��� �j?��?� �?�� ��� �j�<�����}j���� ���� A� � �e�� #Q� ���� A� � �e�� #W��g����U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~a� ð~)� ��� �u� ��ð~)� ��� �u� ���yE��;0��C����O���� �˯�F����O���� �˯�F���(�ϡه֪�������]�5������]�5~�QG�}�>�P��� �`�S� �� �� �2�� ��� �`�S� �� �� �2�� ����?��va�����;���� ��_��G�;���� ��_��_��Q��C��T?0� ��?� �?�� ��� �j?��?� �?�� ��� �j�<�����}j���� ���� A� � �e�� #Q� ���� A� � �e�� #W��g����U�?�v�?��?� .����v�?��?� .���O(��>�fZ�~b/�/�ͯ���M������'�� �}|t�� ��#��X�_��+��?��y�֪����?�_�� ��]Qv����[��fD����r�7�O��q 7�WP�ə���~��P��=+�s���_ۇ���s��WZ������H��h��?��u~� ��|X��)�?�5u���.,g?Rd�A� �W�՗��WD�DN���j��]��N��05����g�,c� �_�����寊���Bv�K������'A�j�'�����y���Gm��C�|o��O�O������7�c�g�{�ρ,,�n��4�{}��ʿ��x������Լ�<3�_C���G�@�ylՍr�.��ͣ�����p�|I���OK�T��,_��X\$�7ѐ�kF�o�`�_����(�.4RM>����9��@� ��~տ���V^4�}j�6
!��x��6;�r�3�ҸgF�?�6:#R2�Y��E~q�� ��x�FO��Լr�5������2��(�
� Z�{�o�χ_�ρ�g���m���ܩ�1��''�	Ebhw�QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE Q^e�����q�:��^<�E��][���o.�� ��r�21�A�E zmy����~�ڞ<�^���KEě�'�x�̒?�M~Z�~� ��|D���7��
�W��s7�w�]{F�2E��;w+ͼ��N�ɭ?��*���I{��g����j���&�bG�kjtjVv���+ɞ� ���TŮ4߄~���x���ܖѷ�?�S�|ըx'����n�� �ږ��,��D����N�'��(����2k�?�����o�����T���mKd��a�� ���5�R�z�q�<_H#�� ���L�Uǌ����Z��}���H&F��_�}G�_���ej �t�cK[u�����k���Zt)R�"qʤ��0��+��(�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� *+�X/�䷹�;�$^)T2�� �jZ(�>%~�� 	>$y�?�Wú��?l�XZ���<�>��z�gǟ�M��"�u��&�Y�ݼ�#y��l:yn�z�ۓ�_�4WL%�����5�Cf~x��
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
(��"�����[k�c���JI�OPA���7�w���^=�F��G����}�4/a+z���  �t���V5)B����#9A�,�������� �wx��9n�&�i�k�*�X�����������,� �B>~�K�gv|-�6\��uiI!��yxYǰ�㒀W=�/�^.��t�oN��tۅ�-��BH�{���W����}>i|K���H剼� �	�!�#�m�'9p�s��x8��P���ѧ�R�z�tW��� T��}Q<���Q���7���]_O �������bN�~~��?��៊���/�u�?hW�0��I�	T��áV ������袊 (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (��_�� <����-���[�"��k�0J�
uwlp�$ H ���{M�y�k��^�gMsyy*�1����@ z���������5+� |V�I�V�E�o���'m#����|y��N�-������� �����4���0��A�|C�q}7Fa�g�^���o��f���� �݊_�ߌ&�mƱ2qG̐)��؟���8�|,�Ms
��My��6� �=e����7�<�.#�� �Y��ݺ��ROT���k�?N��,`�������oo�8�tUP0��W�Q�
��T�*��
(����(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(�&���2�/���LZ������]j�B��q�c� -���0y��t�C��k�"�ZE�\�n�a�`<�N���L��ܘ��q��#9�C�� xM�V�y�k0jZe�f)�nP<r)�A�  ��#
��$tү*zn���C����Zx}N�t/�Ǻ� �7S4x�d��<��@?�.F~������4�W�C����~�/��,�o��L�;ч�$`g�����W�� �o��N���~��=����VZ�M]�#�S����?w_3V��˖h��8�W��EV%�Q@Q@Q@Q@Q@Q@Q@Q@Q@W���ߵ��?d��X�
�~"���� Yo$�c�.F�� e�[������)�Ox�����L�m��B`>��Pdn���$)�q�������$�޷vl|?d�<��?H���0�~yH��� �x+�� ࡟��o�Nx�%�e֠�V#+if� `}��Ř�ߤ��.�����>����٦Ƞ�u�f=Y���^��{��3��g��s����|��>�ŗ�qu&{�1̒7s�:� WsE�������[m݅QT ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ���k�>5��>[%�����XX��R�vI;��ݰ�{��+
�aZ<�F��*n�>S���w:-Տ¿����"qg���P$I�[ޖ� ��S�<?0�XGYYX2��e9z��/���|d���_������&�b�"���,�߀�`�1��D�O��շ�/��\E��&����'�dh��w& xy��Ϲ�<���{4�Ƣ�?`�C4w��:��U��\��EPEPEPEPEPEPEPE��l��Ꮐ?�x��ZiZ|y�m̧�C���x�N $ r��G�A�_�O�Ǌ|E ���t:V��6�q��_D�Q�J����?�� �
��S�Ǎ�M
9�������YZ���,�f��.>!� �K?h��o[�M;�����k}"�q�Y�ny�Y�1 ��D�m�]?Þ��N�lc��:�Vc՘��ǒI5�`�������:������_J�_��t=�7J��Cok�QG�$�I�I�5�E��M�QLAEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEP_/~�_���O�n�|ue�.��"s��>��O+��Jq�	���2>:� �n���_uKO�����c��Z>�����J߲N[�xV?������:�����c��(i�>9�}���-#�yg���J:�� =�;�y!k�� �]�����;/����ۯ���Z���H��V=dP>B~�~��KBXy�j�ER7G�Q\��EPEPEPEPEPEP{]��/�����:v�a��]�8H��;� 	��_�{㷊�ߴv����<e+A����8�����:3pyUڃ�'w��U�lk�x�|�ܷ���%��w��jG�Tln���}'�D��mg��|��x�TU�S�L7�:��7�W�y�zc�\;�N����q�=�����w�,�1��q_����6�b�\��p:  
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
(��?1i��?���1��K�:������WO���?����?/�$*>��� �?�O�
�9�I�;X�Y�����QĐ�?�2�u���� h��V�������hӮ\��@|���?�W���H��"���]o���e{K=fi$ӡ�l}�S�,<�1� ��P9c_+���j]l�gS�6{��
�(�8�
(��
(��
�t� ��|~����B�+�v���@�3wo�ڵY�T���h�G�SHJ�J���M��o�~<�6�5�M.q��p��=�FE� �W�g��=K�'�_|U�3�����x�3s��J����زڍ7Z���EI(Eɟs|�Oi�O�N��;m�5�;�.���W�n$E
;W�QE}�b��V��m݅QT ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( �����Ꭳ�����k��&���u�uo����!�A�\ a��nOZ�"�<��>����|) Qsu�e+�;���l��ؑ޹1T}�'����sL�{�a��a�G|�,�C.�m����hvO<�8;s�J���u�	� s��?�>%x��ƽ#Zì���t�yZ�+�x���x�'������<��(��*��4�+���V[x�Ye��Qe��@4�c� �������M�=?��d�����l�g*��J���g�*��#��X|�RH�ݩd|��P�� tm�肾�N�?��{j^.�G�E:���C.O�alU-b#����c_���e����7=���(�t�(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(�̿�'�ڟ���T�|2����mYr^B�gC��f��+�� �ώ�ϊ�i<�/\��P�9�"�}gv �����>/ď�� U�� �������|�N����7����#O��/���Û�������7&���%��H��OK�Vv��{Xy��y��QEp!_+�L�1���Vm��u_m�����|�;��Ϯ+��}� ��|N��������t����G�뫧؈GvĤ{M�@�1�b|7׼gq.5˿��1��A#�#8?��W��|%�,?��g�� �a��� �s#� ���~5���Ч�F�R\�r
(����(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(���a�����5�	���7�`��E~s~ž$��N� ��'�/fh4�J�o��6�;+Z9ϫ�s�ٍ~�W����sxS◄�y����R��x�)sl�Q�ٶ:� \���i�RS�v�ei��?sh��oĈ~/��_�!+�sI���W�r�`ȟU}�� ��h���muC�T� �IuM}�ڴ��彍�����m���a ��=~�~ҟ��T~�� �Z�yW^�u5�g�2�����?��� �Z�7����W$|[[C�B�u21�@>�T��uaa�+F&5��M��:�(��< ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��?ࠞ� ���s�/b�}�s�;1�>�$f?�
�N��e�<e�s@����K��26�!� Ы*����;�	r�H���?��?f}C�w��$�9ŭ��?������ҿ��>4����w��~M�ִ��$������I.�_���'��`�t����G�]�x�[����ě�X�7A� �
��'o�W��m��M�k��7�H����1G����� �񧙪|+�o�*�V�<���E�6M�׺|���!����-��j%\c��5i?����屽W.�����H�(���C�
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
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
(��#�(�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� �� ��HĿ�u?�%�������� ��/�ŝO� Id��*�3/�'���Y����K��� b����ZW�~_� �/�k� �,��U�~�Wn_��͊� QE��EPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEP��п�@~%� س�� �,��_��/����,��K%~ׁ�|q=L'�ϰ?����_�� �n?��ҿP+��	y� %�_� �b�� J�+��p���د���Ez79�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\�N����)�QpE:�.h�QE�m�(���E�S���6�u\=����K� �gS� �Y+���p�ho� ??�X�� ��J�=�2��z�O��`� �/?�k� �,\�U�~�������K��� b�����W��]���د�	E.�F�^��%�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM %�MM y��� $�g������_�����C� 
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
(���ߴG��� ���5?�$����'��� ���� bƧ� ����xY�����,��	m� '��.?��ҿQ�������N��,\�]�~�Wf�'>'��vѶ�Ez79F�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\<���1� �o�����K_�����Dɿ�M� �cS� �Ik�Ҽ,���a~}�� �� ������Wi_��o� ��������Wi_�;k����A(��F��Q(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh�?h��7� ���,j�I-~W�o�����X�� ��Z�2�0��zX_��c�,� ��A� b�����W�U~Z� �,� ��A� b�����W�Uv��s��QEz(QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE翴W�������O�$��k�;��� �}��� b��� ����xY�����,��	g� '���?��ҿR������N��+��]�~��n�'>'��m�m-��m�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����P�~�K� ��;��}S� I%��
�������~'د��$��_^a�D���>�� �Y��x���{��+��Ժ��� �X� ��x���{��+���ݢ�0?�0�m����sXm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�Xm������N�(�(�X��+�M��w������K_�����E� ��;��}S� I%����3�'���Y�G���N��+��]�~������N��+��]�~��^�&��Q]�0QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE矴W�������O�$��k�K��� �}��� b��� ����mx������,�'�	_� '���?��ҿS����	_� '	���?��ҿS��G�����9���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� ����_�!��h�ߎ�/T�`���_ƒDH)<�� ���?� E+�� �3�� ����� ��>2�C�<a� 	������,�C�X�"�*q��z�g��o�/Ý:�T��շ�t�g��&��(��۾��}���\D.�.���֧��5'�� �)^� ��_�U�ԟ� �x_� q�U�w$o�����k+GPE6���'�����~��R|� ������ �Q� I���W�� �g� _�4Q��?�C��{����ԟ� �x_� q�U�x+��>$]\�x[�z?�.m�I4:u�s<jN�8���Z�|Q�~��K�_�����1�>\���0��`�j��J�hK«h�ޚ+��#�A����O�^�&!��w6���i8|O�20{���{^�d����vbQKE1	E-��/xs�ޒ�����Ny�u2���B�<d�c�c\W�5?���)^� ��_�O���韴g��𾫩��6��G~�j��� !�1��_%x��	;lmh�[��=KM�r8+� |��*�~�n��m{���_��� ��xg� Q���� ��xg� Q�~S|p��>$��[�sI]CAF�kKs5�����h��@��&�R�%��iF�J��J韸���� ��xg� Q���� ��xg� Q�~�_���'��1k� ����u����X�Rq�I�A��EJ�T������g��5?���)^� ��_�]w�>&xO�U������mYRy4녙cf��O�k���_�J�jqE/�<a�h;�h, {�ؒc\�	Z�'�a��t�ٓC�t�?[��T�9���k"�B��=ɮ�U+I��9�q^�=��Z+��J)h���� J)h���� J)h���� J)h���� J)h���� J��jvz&�w����Y��0T�4R��O@ $�j�X�5�7�n��˩����,.4�f��%��,�� 1#4��7�57�� �)^� ��_�K� G����_�� �i� _4]�(|"��f�ֵv2���~[k̾"�+�]���u�� Xx�������L��[s�?��*����t�Rh�������� ���*��j/��R�/� �H���J�'��_k׺.�a>���Hb���B��~�� �8������ �����~�� �Q| � ����C� �Q� E����_�� �i� _�R�П��V=���j/��R�/� �H��?ᨾ�K��!� ��⦶���b��I;H��Q��?�B��{��?��_?��_� �� �T��P| c��/�ߎ�� ٫�!|7���t���[����xwU�R�i��(%���S��� �>�������I:|M��� {[���F���}P�g���d��5�g�O_��Q��/��������H�1�R��=>�p�g�&�k��V*���A�"����|,T��2������<C�V��;ė���޺+���� �P��>h�[���� G֭L�dM�ō}e���
{�Ok���H����!����c��b x�� u����e)�{z�ʄ��}�EU�u{{M��t��}F����Zʲ�*��VRA��i�%�PQKE %�PQKE %�PQKE %�PQKE %�PQKE %�PQKE %�PQKE %�PQKE %�PQKE %�PQKE %�PQKE %�PW⯏_|�O�x�������g{z�J��2���A�k����h��'��������K�t=R�(�0��'�Lq�jq�[���{�U��We�E�y��� O�{��W��e��� O�{��W��e���o�� ���qx�5U�X��Ă(u��{�$$�~h� �ʞ�N@���卫�(ٝ���J韸���� ��xg� Q���� ��xg� Q�~�Q��?�E}V=��O�j���R�3� �(�ƽ�>'�<g���:�m�iWA���	"�k;Xpp�G�_��~�����i� � �W_�Y=u�R�7�*�T�t�|����NQ(����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((���<���� �|��� b��� ����]~��џ�o� �W�?��Z�*�0��z_��d� �+� ��<C� b�����W�}~X� �+��<C� b�����W�v`��a�� �)�Wq�6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E |)�~�qx�M�����;n�n�Οl���Ϊ?媎[xs�����_��~(��^Ҽ�Ox�J�l�����a�-�Dim��������g��q�c�G��a�9{��:(����(�� �� �?���� g�Q^����ڃ$:ƞ��%A��L�=A+�r?f<;�-7Ś����Ũ�wЭŵ�)"0�#�:��~+����� �[�C���WA|#���]��]6��rzD�ec��\ק��{7��9+��\�s�6�u�h�)�Ph�Q@u>�V���������6�ky�<r!*�x �ƿ�m/�v�~7_iZZy~����$C�)�'��F��۞��G_ �Y<2�h|@������I �u����������J�7Tt���������.5[�k+H�{��VbN��@U�+��x��|�c��kM�?��k.�z�滺*7�=p
;(�%c�?~ӟm%]蚴wX�����ۚ��AY͚���)�W�p���E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E |� N�[��x7ß-���U��]�ȸi����{�d �C_��� �[��S�R�͵ݔ��D� �������u=L;��(��+�o�%���� �"/�+����	C� %/���� G
��� &5����)�WҞA����υ<iG����0���S� �Jk�?�� �M� ��<��_�q�}X�Rm=̶Żo��Odd�_XQYΝ:��q���?�>~�~3�����%�˧\�-^Ж���=U���9�*���t?���|+�U��u����	P�u8 �zdW��O�z��_�߃�������[�Ԅ�&��sٔ���#�xX�?�w���h����s���+��=��u���a�9��}"s�h37�f�u#}�a��� �rz8PÊ�v�'���k�o\�ϗsi)kI����g � �85�K^��.�к����2�Z��B�e����?,�����d��\��k��t�,���5u��qES�u�?��-���]E�_B��P����X{j�}�卢�E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�P�~�����'��K!�K�X7(m�ƦH�{�E��
��uȄ�.�o"�jk��0J�g��z4QEy'hW���� &����u� ����%~�~���17����� J�K�G�rb~{��+�<ѴS���N��(kzՏ��{�[S���N���nne8H�E,�}������Z���`xc�R]�[iS��O����HHmV�!Pq�	f� iڽ��
S�M	�7���U+6�q�OX�r=8w� ��
���O�����1i��Sxت��8�B�>X�>����+�<�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S���N��E:� m�(�S��:��?���'� د��$��U_�� �g����� ��O�$��
+��>(���g�_�J��8O� دq� �v���_�� �J��8O� دq� �v�����V�F��J)v�6��9���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���o���N������ ��~�m5���}��x���?�"��;�5�u�7�x+�]Gc��L�sH��p�Do".��(�2@�\�~��ſt?��5�w�5��yϊ.f��$@�Y\;M���s�\��H�8#��?u߄�6�<+�;O�j�|�\�2RE�!�܌`}n��N��c>��ԛ�Ts4QE`hQE ~�� �=k��m����żEe=����|�'���p����q�����լ�=6�[F�d�������C+� ����j?�;���Co�M)V=V�0�`�1��?p>�dt�O���s�g-�:�.W���(��h�k�8Ģ�i�i���ખ���n~������� ��~�_�T���ց��q� ��uψ��i|h���	� n.?k���V7�q_����� �� ���O��� �C=~�m5́��� #lOƄ��i�i�D��]�����]�����]�����]�����]�����]�����]�����]�����]�����]�����]�����]�����]����>s� ��(o�'�����OA_��!� ��7� z�� K��������^�����QE�a_s� �(�������_W�� �J�)�8� �D�8WV��1��6~��K�Ѵ�ў@�R�4m4 ���� X�gl,|!��xv]y͢�H��R�,9�f���Ѝ��p� ��h1k?���f�7K��gy	���П��G��jRFԛ����EW��QE��� ����*�]�x�V{�N%�,s�I�0Q��O��v��+�G�	���t�f�KV�λ��Y�猢}���H�M~��5�8I��W��W�,ؔR�4m5�s�E.�F�@	E.�F�@	E.�F�@	E.�F�@	E.�F�@	E.�F�@	E.�F�@	E.�F�@	E.�F�@	E.�F�@	E.�F�@	E.�F�@	E.�F�@	E.�F�@	E.�F�@	E.�F�@	E.�F�@	E.�F�@	E.�F�@	E.�F�@	E.�F�@	E.�F�@5_���� �M=�����?�� �&����1� g���^�EW�w~����i�� �����_��A���bo�� ��s� �sW���#�91?=�]�������]�����~���j/����3C?�uDh��6�y}��������7�~�?�g�� gI�j��z��h��%	.�� E�,݇��f���+_���MK�>%�k�R��1� ��ƿ8��$�+��,w���4��3���B�V�.o�g����V�i�b�#�%���I$��]����.?�bӿ��:᫹�� %�����Zw���^$~$zOc���]������Ģ�i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i���i�i�;��?���'� د��$��Q_��`?��� �W�?��Z�'��D�0�>�� �U� ��x���k��+���*��� �U� ��x���k��+���*����Q]�0QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE W�'��� 'o����� �޿f����� �������Co^v;�k���|o��C�	W� &���n?��һ��O�S����/�!��z<Lt�B}�>X��ǌ�����|��J��7�� ��q� ����]mJ
�mb*Iƫh�zo�.���+�i��-�hg��B�E"�O*��<�*
�/� ������]�S�}��Nm#ݮi��}£�>P�h�|���g���x��:2�g�Nj��
(��4
�~�V�~
�@Ҽ[�����>h\�*�#�� U���0@#���n.�M_F~���š|v�w���@������-'\o�� �RG�a�
��[�G���ٿ�$wr���%�2ìi񜒟�4c���$xn\��~�x^��U��kE�z��}
�[]Br���*��5�Xz�|��U��"�Q]F_�U?���=� cM���]��u��S�M����4�� �%�s�?�#Z_>8� �{��x#��� �!��d���	�� 's��o� ��d���7��F؟�(���N@��( ��( ��a���oO�	���tȏ�<]*�V����oI���('�o���[�����jY�����c��3B&�?�J��?�c\uqT�nΈQ������?�x%�<C�Cu)��0�ߓ�5�]����|G������?��k�IiI�ؖfc�I�I���� �_�x̞𦷯Ƨi}/N���S\^�~�N��En��}+���A�H�m�'�;��V}J(I>�9ף�zņ�j�Zu���}٭eY�I��O����$> �5M
W���Vr[��QQxgŚ߃5H�-W��o�9[�>�����Mc���<2���"����/� 4��'���6�x�HV�
��яP�?ї'xW�g���� tW�����R��Ӯ�wm��dy���ʞ��},E:���)���z(����(��(��(��(��(��(��?�?�h�7� z�� K������� ��� ɣ�����.��k��z��~�EW�u�}�� �� ����G� ��|/_t�'� �x���� ��]X_�DƷ���kEWў@QE W�~������A/���R?��W��{m|�� ��C���W��2�YE�M�n���G��;S��.�Ə�*(��\���(�v��e�ڻ��C�����6҆�	��z��� �n�RO~�:M��|:%�����`_�zq�W��{����y���Q^��QE QT��f�ú]֧�^���u�g���GJ:�18�@*;��l���T���F
�=I=|�@�O���.�_�z|w�ψ�$>Q�xa��g�O�G5��O��8��|�^.�F��;��&"� r%�'�EpT�BGS�yKW��C�~��	�7!�P��ᨥcMN)}U�X���7�[ɼ����+g�V�� W���;��d>𦩮���Z۟%���A�'5s�����]��|W��WF��#7�ú Ǡ2.T۞k��ṃ����1���xwƶ�������ێL�m�w
>���=�V����}���iװ������H��2�G�_^|� ��x��w6�o�Cx�@�V�l.����K�~O��kO;MX�a������Es�%xs�ׄ��K�]N=SI�G�F��F�yG�O<��WO^�i���U�QE(��q4v�I4�,0ƥ�IUU$�z (���x�J�͙��5;=*�z�_\$1��+��s�
U5��ׇ>���9�Oqp�p~�c?�����c�_��q�jϩ��[��oܜ�j3����
���ah�u�)k-�}K���5��Y�$�u�� �k՜~q�4oڻ��2Ei�#âG�V��m���ͼ���0�q2Em,�B� %������K�z�����^أr���>�
���� �����@Zn�e�Y�w��[��I�Om*�}I�W�?��%x�ᦦ����A�h7JۋY\2+�:�k�fWٿ
� �Z}�����"�5-*E��� �d$���?�k��6�z�%��~����o�4K�cX��M��i�.�,q �I5��� ⧄�*��^xG^���md�-��
x�����?nMS���-�7���o�q�&�	'M} �S& T�g�1< >�� �L^y���i��W�r�����%\qJ�U�L��Ù�}�EWq�W�_���?h/��J��R���di*n� �,c��K.�T�>vy}+���c��u-�ߎ�+[��l��en��ЃA�*��ӗ,�:#FSWL����?Ṿ:�B�� �[o�5G�77�O�(W��m� ƫ?�S��������W���77�O�(W��m� ƨ� ����� E
�� m���^�هէ�������o�F��o���u���կ��g��eU29�
 ��zEz|�3��;����Þ;�.~���J��˸�(ϝ�_A�־,���>>i�R�]��R������V��}1d㩊�9r�3xє��?hh�ŏ�no���P���� �Q� �������������2��>��=��� ���������?Ṿ:�B�� �[o�5Gש�a�i�?i�:���M�x��O�u�j鯵mCI�{��US$�ybϰ�E�B/�&r�g`��)���� �.������k���OU� �]��q����^F?��� C��(���N��/�G�M7�� �p�� ҹ��v�h� a� ������?��j�p?�~�&'�G�QE�hW�~��� �f�	��E�P����9 ��oo�FW#�d�^g�W~�^�'��(����!Ll�ӛ�V�q� <������^5�~!x��^�&�q�j�m�k���@EP8
  p �?�T��nuҢ嬶5�,|[�7Ư]���W�{7ʑ�V+x�%b��zu<�I$�E����E+h����� �\>�ŧ�Lu�U�X���a��ӛ]B��;�i�c�20p@<�Qf�=��F��Y� ����� E�� -��������%�� �o���3�������W���7G�_�(7��	m� ƫ���4�ֽ�+�m���<K�>������r"�I��xׅ p���ڶ���ir�ԣ*j���(����(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(�;�����(ح��$��M_�_�o������O�$��	k���Q=/�ϳ?��?�p�!� �Z�� J�+�R�*� ��?�p�!� �Z�� J�+�R��_�1�?|(���9�QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.QEp��(�Q@\(���EP
(���QE(�.���~��x� ��?�!�����a��?��<������ם���:��>�� �T� ɼ���ƛ��$���+�?�%O��ψ�i�� �KJ�2�p� �Wﰯʯ���C?	�ɼ����� L��O�K�s�G݅���ʬvp
��Y�!������ F�,��4��Z���a��6*
u����b�Q�w?�
+�k�ُQ��~!=�^u��"�h������C!y��	�!��2TxU|ܢ�'nzђ��
(�����?��?���+Mr��ǋ����S��ˆ�4ۗ9�I��HǞ�����]�:��%(�(��3� ��"� �z����m��ŷa�O�@F�{3|���Lm��ƣ�VE��Y�۵���F�T�y3���
���
�� &����m� ���̯�� ���o>� ���� I.�<G�U'��o�'���ׁ�ݿ� �����Ɵ�'���ׁ~�� �Aq_��́��� #\O�QEz'-(�.���u���=������ �u`�5�g��n�o�:3�AʌH�_������ �GЮ|���+I��p$��nv��D��"� ]�Ggv,�rY�I>����?���N����3HŘ�f9,NI4�S���&H�F�G`��2X�� ��>���� c};��r�� ���j�Ŵ�;Lc���O��c�j~]����
C~����ivp�Y��ii�	j8
�8 z
�>x��
�+�8#TV���z<�s+�ZB��j��hRT��SȩQ�W1�Y��Z֋�.�XҮ�l����F����A���e~���7�HKm<�q�=]�J���iH�d,��VS�5�E_7�@����m��A���m�ےy����Wҳ�RU �TU�2�F~8V׃|m�|=����V�F��Nb��}�=A��{��#�b�_?���O�O������ h}4Mi��<}gg�_�;�����N������>��f��_i� �/$�մ�����3������	�_�߳OǍ?���W����4���jV(��5ʁ�G}�˟�a�s^��,�G�Z�'��=R�(���EP
(���QE(�.QEs�?�(W��7�?ް� ��z�k��O�(W��7�?ް� ��z�k��U�zXo���Q^q��O�w�Jg�?����}�� �� ����E� ��ua���g�Q_Fy7
(���W��o�.��{ipUw�o�V�����c����_g�r��῀^��7�� U;K�y׳���1�{��2OJ�I���#Y����gź��v���e`�v%	�*�UÞk��UQ��[��'��9j(��3�
(�j��f�_�����Im�I�}KlQg�Z��>Qۖ親1sj1�Rj*��g�	w�v_
�9�|y���W~#�`��ai ��!o��Oz�v��z��t['K���M��-��aH�E
�� U��jT�8(�ǜ���QE�
(�f
���I4�x�C�i�=Gľ#��M��#�&���
��1$ ��H� ?j���ߴv�-�<�?�m��e��`ɏ�,�}����,zO۳������|7�^���$�m����NָoU2��9?�k���X�7��J�.U�-½��G������-�h���xv�>۫�Fp�  ��wbz����J�z� �r�-��?�Ŏ�-�Ǫ��f�g����	H�m��u5������V�$n���7�+�����i��^�e����G��̞��O5r���T�����+�I��,�x�R0U��A�OE}����o���f���>',�<L�׃�飒-��[�=v�� e�R	�����
#��m�2k��2���ê۰�G/>�\�q�ҿ+�qT�:�[3գ>xjz����F���s�z��+]is�]GG�Ȇ�?���;_��h>|Iо.xK�W�.�ޗ�D�th�vu9z�Nk�"��� �t���/����K���D���m|Fؤ��#>�C�5����%�"�>e̷?Y����O6�_��R_ښH^o�^�h�M~�� ��j�A� U_�
���w�-��:�'����:>�5��c"��f��5�3�j0���공֧�\�ws3uy�1��y�ʮ�]N�<9�3�fV� �� �_<m�x_G�I�j�Ikﺥ�7���`k���,�v���{�z��$PX��KV�c����HQ��H���Z񠔤���i6�կ��7��g����:l3�X[�r� n�_1��D�DS�=NI�K�ho-��$�	׎E�==EKE}Db������ۻ>f���|5��euu�i��/Ĭ����#�߷� �0�+���OJ�����;� �a'��Uc�JA���Z��<���#�x W�y���>����xk\�Rl,5Pe����"�N̿�2=㯅�Ex�Ίu�]���E~�� �#���E�� �[� ��O��}�}_ᯌ����e������0r��S�X���c�	/y�������ZK��d�����}���emi��Қ(�������N���������"�O�a�jF�@�y�;��2�9�^1�������b����Y��\h�7ѭ�ֳ�84˄�d��)EW��QE���#� ɳ|6� �,�z�y'�� &������� �z�}M?�>��'�0����_|��c�^)���_E�nVX�a��).3�� ��ʸ_��lz���6R������lۿ?�b���.Y&���(��X�B�(��o�/�M������
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
(���QE(�.QEs�� h��7���+j��I-~����� &��C��mS� I%��:�� OC�>�� �T� ��x���k��+���J��� �T� ��x���k��+���J��� ��Q]�0QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE W���� 'q�� ��� �޿g����� ���� ��_�Co^v;�k����l�K�	U� &���n?��Ҿ̯������ ������ Im+���� 
&U~6QEtG�o�:�/����|E�K��¨�m'\�3ٔ�`��� �&�f�A���u/	x�.�����l�3�X�e8*@��	��?f7���~�ш��]��ˣ�07����y� '�Ng[�C�ǚ;���NGg���E_����W\���{9t�R�f����v�R)�R>�B��
(��.h����5�-WK���Q��.-��b��AVS؂�ٟ�� �����~��;O�j���+�f������?��/ oŪ�|^��|B�<[�����F�5��]���!�wV�Veua�2�f5i�E�~�W��Uo�7���o� ��u���������g�|9q��]��`c��Y�<�d�O�A2�	���
�� &����m� ��CR�&�
zTI�� �=� ���	��� �����Ə�'������� �Aq_����o���1
(���*���Z�]�}:������'ݎ4R��� O�V���(���>��w�m��]� ��SJ\�BA�c�*��t��5N.O�Q�3H����>1^�v���x���\K��[�� Qj����O�3���Q_/)96��%ed�|Ї��8�?ҙwGu�X� � c�M���k��p����~�_�#;o�����T�旘��[?l(�����
��Ɔ�(�W�4g]������D��c���(�@?��+C�V�� ��q�"�X����+>�H�¾�� �w|p���_^\l�� ���L��*\��y=K/�'����[;��.��������H�C�F �؂+Jst䤺(�&���\W�_��/���X*�j�l7�6�%_�Ï»Z��ԕ����(� ��( ��( ��( ��(�?�(W��7�?ް� ��z�k��O�(W��/�?ް� ��z�k��U�z8��Q^q�����Z��_�5}B���$�cTk��0���d��TV����؉ǚ.'�'�|3��E���k���Y���_ i1���φ������c�z������9������"���	�]X��E�.1�Ng�'��������^ҭg���/u�� ������l�������U���%���<�o�o�?��8x�����j7C+#��	�ȣ*��8�$�\EWnN��J� ������ b��SC�_3��}"���\�C�_D�d��U*�ۿW��� a?�� Z��� �&�c���׶�oT�5��v�s]t��*+�gZ0ө���6��>8��wi�j�x[���}V�,Ip��o�����G��W��*������cºz�i�� 1$�yޒF�'8�� 0 ����8xQZnpԪ�nQEt��Q@|�� 
��'�_���z|�w�,�e�)�G��H>�Bdt2�ھ����)ħ��{��7��xf�=>%S����Y�����+�S��vݛя4ϕ袊���%�����x�4�#DY���@���>𦋡ۅi�P�F`m�5A��k��}��>4�Ȯ�s�>���c�k���|�L�Ľ�QE�G��I��߈Zy7��}�۾��W��@�4��w���b8��x���_��x�����6�)��%��,N�KI����GCL���������&��?<#����o��]��d�7��#c؊�j���	S�f��x��ɽ�mQnc\�\'�������g�N2<���G˟�RMr]#�[���V��GP��}����� �U�_���R�!�/�~��P���[;� }�Y���e�Ȫ�1�����(���Γ���/��,��]+�?���T�V�N�͈���H���;<7|7-��;�����  ��ϱ�k�_h�7�U�G3>��=:~����c���c���-�T��*�:��e��� �� ���3k����Αr>K�v�uV�a�U�#�������G�� �S�k�x^�↓j>ߤ�ռ��KflG)�䣝��8��?�J;�/�O�-3����.?ݹ����g��x�Ú���ۋ�7R����z7R�?#_��?���g?�:�����W�ƕq���i;��Z�Fd(�b9�����Qj�jEi��E��Y��Q^���|R�g��zd^藮�����L�]� `����������+��h(�� ����� ��>� ��� A�Z�%��� ��>� ��� A�Z����G���0�W���� ��Ƽ� �]Ur�������ן�!������+��h(�� ����� �ھ�
�� @�����W�O�Z� �^�_UO�^��/��QVIWU� �]��q����_ж��  ����� �	�禼�����v�QE䝡_�� ���������� ҩ��b�j?a_�4� ����?��M^������|(��(����
���
�� &���m� ��5�%|{� J�:g��v��"���/��4QE|��w � �|9� ��M� Ҩ넮����O�?�2i��U\>$)l��Z(����B�(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(�;�����~(ح��$��G_���w������O�$������Q;��3���%G��7��V�� һJ�U�ʯ�%G��7��V�� һJ�U����1QE�s�Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@~0~���w>?� ~�� Hm��~�?o��;�� �e� �6��c���N�?�϶?��������k��&������%����{���� I��Z��� 
&U~6QEtQ@~�߲	������ ��/m _�T�A� �� ���F�HA_��W�o�;�A������4[�k�/�1� ���nO�?��=�b�� ��?3��O�σ袊�Nࢊ(�c�ڊ� �m���\�.���ëخIQ�n#��3���r�v���� ?�,�A�1�GT�n��ӯ|GiqmuIc{+�WS�AZ������N��>������Y�Y�y�a�$J?��>�q�y�q�
�4�M�c*i�M�� �� ������ ��fk��	�� '}�/�� �}�~�W����� #�񠢊+�9B�2��^,7�<ᥐ����d��˰g�� ��~������V}K�������}����/��1ڗ�Ӈ^��Q^�}� ��O�o�#"�� �)���|�_R�6!���C��O�q� ~�����G�Χ����(�����(��� �2��!x�GE�.�� ��\�t�����·_�9����{��[QHg�G��O�~����Y�Mf��4'��꓌{�O���U��� ��M�m�bR�>p��u�R��������E�Uei���+��(�� (�� (�� (�� ���
� &����?��޿k�[�
� &����?��޿k��� z��aEW�u}�� �������qs��1`J�����������$�� ���� P��utΉ�?I&�ޓq�K��K�Q�{t*�F �A���j~̳��� ݴ�d�,�\�!9��1��dc=T������y� �o�:7ǯ����u�ؗ͵�Q���P|�W�IwVa޽��UX�npR��/#�v��!�Z�[�M[��!�6z��1�h���*�{�)pA�v�}����QHf��h�� ��╮��[�^��lP�ͷ-��t�#9e��K���v��Xx�E�������������7P���E=u���g����ԣ�K�{�i����ĭ�3�,�ğᐝ����9.1��+r�g-��^�2�G�MQ^��Q@Q@�T�'�F	)fc����� ��2x���$����ڕ��������� ��O�:�h?�o��g��\�=R��*��'�wa��(���Ӹ��=;�� ����]+M�l���(�"�'FwN ' ���� �������Y{� �k���[ۨm�]�L�.@��9��l� �"��� D�R� ��� ���V�4��*S��rg��7��O��� �e�� ����'����� ����1�/�'z����� �Q� E����ޥ� a� ����� ��3cK���mS���s�����g�E_���IR � S_��?��_���_���9G�1�/�'z����� �W=iU�nh���i�t�fxu�?��?���_���9G�1�/�'z����� �W?�����5�s���Z����� ��B�F��Z)�Q��V�k��5� `�����O���ω�}�����6��J�`���by(;W�E{XD�;4y���tp� ��|]� PM����n�}՜�1�*¿�-6�GԮ�/`{k�Y^	�a��IVR=A~���W���c��~"o���f�5Ikq�	�G���N=���c���5м<�_+�|%EW�zEP�|���?��%]g�:�Z;`\YM���_��p��a�b+�k�c��<�DZŦ�/�|d��tk�& r��q��~��02j[;����n�g���Y"�(�9�9�Et��J��1�J5=O�v��'�:� ��'�&��O���V����!#�n��es�?F��}�^�:����y����
(�� +���� $������ �/]=s�������%�e�����+�h(�� ��������� ��� A�Z�%������� ��� @�Z����/CƗ�¹o���K�c� `k�����-�[�Io��y� ��[1-��j(��X���(�ܿ�?�M������X�W�(��� ��~�
�� E��z��<i|L(���J�����?�� �&��z��5_���� �M<�������u
(����+���� �O�{� ^�� �L���_����������� Je�G�G�rb>{�Q^��|{� J� �p��b�� ��|� H� �o�� �b�� �7ψ��i|h����+�X+��� %������o��G\%w � ��9� ��M� Ҩ��!Kf~��E�'�QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE y��� &��C��mS� I%���������~(ح��$��C^>;�߆ٟf� �)� ��C� b�����W����� �S��x���k��+���z��� �񉶍��Wu�a6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6������ 'u�� ���޿h��{��� ���� ��_�Co^v;�k����l�w�	f��������������� �K_�6�C��+��o_`�N��&U~6&�6��]2mih��&ڂ�O��,�,�-⻳�����t���VS�<j�\���W�O����h5&9'�&�3:rKI9ck!<��'�(=J�>l���_��^
��+�[1{���z<m�dF�S���ҿhO��������´��2�ڈM��ےv�=f\���`�Cپh��J�Nug��tQEp�!EP�?�O�;� � �C� M�5�7��?��� �w�� �����k�r���߯�v#�Bm�m-���m~2� �Bch� k���0�+������7_�� �S-�?���B�SҬ����´9� �8�+���Mz�X��J�(��B��� �g���V_�
�� �V�S��?��%G�?o���/6�p��5����Χ���M�m�����&�6��E��~�"�������R�� �V�z��_?ڼY���̾�� 9�M|�����ER�� ����✙;4�#ܛ�#_�k���$��m���cf�S����TL��<� ־���F'�[��M�m����`&�6��E�M�m�����h�KE6Ѷ��.�?�P�� ����?��޿+�c�
� &����?��޿���z��aEW�u}�� �\�B��� �\�8������_�J����5Ն�4Lk|�3�F�Z+�yG�� �P��g���O��9h��E����,��ut嗹��W�U���� �;�Z?	�h|w��2���cq)���%0:$���2�6���Q� �����*}�|uEW�v�Ioq-��S�+�<L9cb��C9��-��ӑ~��V�T�Wƚ%���7+�%ҏG�ݎ� ���k�_��-o�?č+��o�ъOj�V;��!f+@��?��4O��|W��w�jP�co�Cѣq�ՁR;k������{�̭O��lt�h�KEv��mih��y�� �� ��?g�� �_S� �Y+���}~2i�X�C� 2n�+� � z�����Ux��'~fQEyg`�xeI#b��2��GC_�Ɠ}����E����>����_�-~��Ȟ<� f� k�e�i�cp{���������z�ZRG%h���h�KE{8�F�Z(�	����Qpmih��&ک�h�:��w��V��i�q43�\ x�Fee<E\�����O�!������
�mW���K�a��J����>��z7��Ǖ���#��0*�� ��_��� �NO5��������qᬄ��=I�'��� p���^�ȯ��T�㺝~�?4(��慨�gX�ҵ{�7S���qiu�X�uVS�5F���(�� +���	��^]����]�;�?U���P����4kw'����{��yQ��
��u��L�t����4넺��)"0`~�;��j�2�Fu ����Gmk��S��O���7��-��cז|�+��}��)�Z�����'m�\��� �g���^藮����?�L�]� `������`�? h���S�
(��?rd�� �����v��z���_�'�M����������o�^��/����_��� ��?��� �]]r�?�����ן�!��b[���QE|��Q@���p� �i�i� `+_�+�6ה~ɿ�m??�k� �
����o�G�/������V�$����J�?�� �&��Z���o���� �M<��������u
(����+����� ?��o�)������� �O�y� ^��L�����LGw�F�Z+۹牶�?� ��/�cm���6��&z��@� ��ɶ���k� ���|C�ԍi|h����+��X+��� %������o��G\%w� ��|8� ��M� Ҩ��!Kf~�m�m-�W<Q6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.m�m-\�F�Z(�	����Qpmih��&�6��E�M�m�����h�KE6Ѷ��.����� ��G��mS� I%��
����� ��~(� ح��$��?^6;�߆ٟg�)y���E� b������������ �R� ������k��+��՚��� ���F�u�s�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�_���P?k�� �e� �6��I_��P��������� Hm���z�X��p� �-Ff�� ����M�}���>� �Y�?�n���x��?��޾î��_���F�u�d7mi�PvѶ�E 7my�9�:�� ���{���g�ۓq�j�2ֳ��c���`�U#ר��T�+n.��}<s��gΉ�_�����gL��q��B���G�������������/�-|}�C����I�c���C�<6W���[K� �7�hأ�"�d`pA��|�j.����i�U�QE�j}� �� �����?��s_��k��	�� '�����g���� �� ��b>47mi�W�r��_���V���P^1�2���.�3� � *����>�Ə�����&mZ�/�-9Te�� �E�u����Þ�F���4��
(��p����� �l����zgu���o���]�~��u��/���3�O=ţ_6�T��e��j:T��>~�m�m:��S���<�o���F���Ե�|DԆ��� ���v�].�r}6���)=����f���W��1s�'4�(��=������u��4�<1��e�uk�������I=�4���/�'o�_�� �߇�|�����U���F�C����k3¾��w�4�N�����H����8�"��
կ��y"���.f�ݴm�QZ7mi�PvѶ�E 7mi�PvѶ�E |�� _��?�a� �����_�_�P��4?� �a� �����^&7����p� 
(���+�/�$�?�{� `�?�i��k�/�$���/� �.�k��X����i����+�(n���� ��o����/�Qy�j�'���P�{2�b��CE-�c?�>�־ �L��%���B|�+��[�f'˕~�GfV���������m� hφr%�QG�5Ɠr�8���� v@3т������}���P�mu�R�2�x�NXAA���Q�2�g�J����Q\��_[� �?jc�c��"�)�~eP����7�6���)��V����J$�*J����Ѷ�1� �u��� ��"�|Iu��:���ϖ����#�O��<�ƾЯ��QT��<�E�ٍ�F�u�S��s�ʡ�J:�B������3x/ƚ� ��A[�*�{u���� A��6���_�?����EM?ĐG����w��a�_1�� lW���������h���(���H�������|O���P.m��ٰyx�,s(� u�m� mҿ7k���kP��S@񖞆gӧ������3�F8'���[ѩ��gR<�h���F���/�����<O���#T�[�yq���;09R;Gjݯ���G�7mi�Sݴm�W���&��ᾫ��q�Cj�-�U��w�
{���1�&�U�һ�;��m��� �����_�^0�e���ͪ�j�5�h����Y	��2>�����>�0�՝:��y�R���9�x���ς��j��K�i6�^\-���A��	 �{����^���=OS𔗊�u������E2�����a���ToJ��IB.?g�F���52�����W����T?
~<C�jW�D�LK���p�q�6��E���5�JΝX���!O��G��6Ө���� h���_���[U��'ı&�]~�1�A��?�<��+�_��f��κ���5���N�ZŮZ��F��Tl�9��*��o��_x~�C��m�i7i�kK�ã_b:�9�A�:�h��hΊu�4�=�W�?�W�]�~���k�����]�\7�.��x�c��'~��^�p�9r��c%%tQEAG���/�R럲�����t�Z��?e%'��9����ǿ�K-.KٿR�ue[�\΄�Q�y�7�_bW�P�}&�����h����L�_�&�U=b��KI��� ��.�R?�l�3G��E+�F̬6��{J�S�
(��?r� d�� �f�i� `;��omx��w8����09H�?�W�W���������cv�+�a���������!묮C��V	<m4�0��zZI*���OJ�l�[���QE|��Q@�߲j� �4|4� ����+�6ה~ɿ�l� � �k� �
���z<i|Ln�6Ө�	)��� ����� �	�玿��[�Aw���� �_�y8��;��B�(�(�
����W�1?���zM� �2��~��¿�i�?��o�)����Ɉ�Q��h�N������A� ����c]�����*z��I� ��C�~�l�� U�Z?�"� Z¿�kO�G�mQ_6z�]��/�.�?�d���:����?�_��_��� �TuQ����?{6Ѷ�E}I��F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u�?�z� �<|Q� �[T� �Ik�v�x� i�7��?�+j��I-~W����w�g�� �J_�8حq� �vu��_���J_�8حq� �vu��]X?��>0��+��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
�@��5�u� ړ�=��#��?�B���#����t���5�kG��l�$n�3͌G� ��G�ҿ��j�Ě棫����j2]\I�n�سē^^:J�'n:���I�r���0�����^�S'�Ko��L�W�U�� ���i��d��s���ͮ&p���'��G��}��ں0�R��c*Ѵ�QE�s�Q@Q@Q@~y�E�c��� �h<�_3�:t�:� �� Y����Ci�ē��J�$n
�0�`x �ⱫMU�+4�����Q_[~�_�#|�G�%���A5i��24˃Ʉ�Fܔ=�^���5�Ӄ�'�d���������w�� �����k�z�����w�� �����k�z�|�߯���Q^��QE ~9~�_��~2�j�}��x�俲dR����Y���Q��̕���@|�?h����u�X�y��r�w*�Wےz�F����_�K�MK��&�kZ�M��RE�#0�C�W�����̶g�F�2���+���D����kR����e4���7� 㹮�J�G�ӳ���X��(����6Z�w�%�k��e�m�[^2�t��@v#�~$�;�^�_Q)%$x�8�0�#��<D</�3�I�-���=���N<���A^�_� �S�3[i>ѾY\+�Z��j�d�k>Xa�ܘ#���QYV�%6ʧi�~c�E�g��� ���^�P����[���O��ׇ����?�����ݯ�?eٟY��> G�[��^�+6�����<(���ڿxgL�o�t�F���J��Kkkh�F �>��9&�%){G�9+Բ�F�Q^��Q@Q@Q@Q@Q@4� �E�D�O�������\-� ���m~��S�mi����7������"������=�C�_�5�c���K����+��
���	3r��S�����h� �:� �
�^��� �g���¿�ͽ���״��62�
d�L�>�� z�z�ûU��]`�ת(����B�(��8?��ɳ����,ɆfX�Akp�¥�f�_�kwc_��SV�l��.�M�m�������-�]�,n�YXw?ƭ5V<��p����Q^��^��ן����I|-���o�-F�� ~2B�PU���󒋃qg����
(��f��k��9�v��m��j�]��[�@a؎�FA$w�����㶍�B�2ӼQ��Et@�Q���� o��s�=ԃ� ~׶~�� ����6�H�SM�^�
[�z�����D�+뒽5ن��g�0�O�i��yEfxg��_��?��%�:��
�mue$B8#��AEi׼y�_(�F���O�-�6�ͮ�M��5E��j��v
��3�_W�d�&���U�`U��A�"��H���FN-4;4W�� �'�7����-&|�J�i�v�I��jޛz����j���nqt��#֌���QEY�ÿ�i���|/⦚�����h�0]�َ�W����Q��9���� i�+�l�}��S��#��Z�$�U=ՇZ�y����ҟ�ݳ�G_���F�6�r�kI����������C�Yj�Z�y�[��4W懇୞#��U�~�z��2i�����.�� :����+�u+7�Þ
��9���\�zW�@��r=����+ns{	�?E�#|H�������:�ZV�h>id9gc�G.�(����1�Hx�����g���h���MA��l��p�9ۜp8Q�'��(�f�ƍhj�2�ֵp�G! ��%} �z���	��$˧��ſZl�H��#��̪F�Б��ķ�Mqʤ�R�p��*+�[�T��?����i^t�����u{��|˧r�� A�=�{W���q6���C�����?��?h������?�~�#E"�1GS�e8 �⿡X� i�:��3��i"���G��瞼�v�gnf~��Ŀ����	���1��-uXY��`a.G���0a��&�~�`�'��Yx���g�����Ah�^�p=� �~�����i���l��x�����f�Y�z8�ц�*��[�֥��=��(����t�xOV��n�Z^�l��0*��è=��~��ó�G�Zօrw\�w�YJ}Z9��꺥����j7�	icg�\\Jp�ƊY��@?�~�H�@��O��_U�� ۷�l�'N�{�yX�{��ۆ��;E�_�/�?����M�ז�g��<��:���Il�������z�!9(��d�����/�;���+������1�%���V2�;�O��
� �k�袾�1�J+��wp��*�~|X���O◌4Yci��ݮ1��3��8�R��� ���#��?��$2��+�oc�/ȷ(N���Y�u��k�*E�n,�!.h�QEfY�1� ��D>%���'H�K��O�P~�,�����3��������Z�����n���|%��k��,�Ƞ��"x݃���9��_�)?��Ӗ��z��FM�Ҧ�G㍟��{�1pJN�mJR�m#�J��� ��~���k���3Px~�Lk7�핽�ND
ìhFI��e����?च������ �n|3�܆��U�p/�c<P��T��%����?�>#�����v��-6�VԯeȆ��g�8��aT^��(f^lF#�~*\���.�(�,�
(��?r�dY�����ՇA�[����z�|������ ���t���Yʽգ�� �;OЊ�����"��zI�QZ�S�G��1�H��_�~�|b�~�M�=nG����܆?�XX��'��z�1�X���
(����+����M�%�=e�-��S���N�`� ���$�\���k$�<�=�唋�KJg � �g�z'���r�>�Z(����
���
9����'��e��uepq����� ڕ��r>��T�i�_���2ډ��r�$�Up�� ��G�=ˋ�g����
�~�6���fլ�]6��n`q�]N=Gp{�d���WK��T�C���=Jv���ip�N0�21?��j�i��袊��T� ���;W���|I����n�J���� U}�	YH0 y���_��ᧀ4i�����'!e"mQ��&�,�@��d��^��KMǚ��{X�m��?ٟB��t�-u?js����n1��7�d}��ܱ�a��� 7�_�~�y���Z-��^2���giܼ2��)��e$)��� ~������$��U���a��]��Ŋ�N�&︪:���_�V���[EooA(#�(�*���  �
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
(��<���� �x��� b��� �����~����o�V�?��Z���D��lϳ� ��|��"� �Z�� J���km~S�(� ��<E� b��������u`� �c����F�Z+��mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� �/~�_����C�����^ɂ\Y�3H�P��{��� �O�xj�q��ռ_z�D����?�<�e��[�����>0|Bּc��~)���YVYⱻ�HT�TC[��X�%� ��I� C�?�:�� �k�Ov���e���?���g�?�ST��a�K����-4+6&�y�Ǚ�Ǡ' d����k�G���/�u�� "�� ��I� C�?�:�� �k�XZ�w��J�N*�����?g�)f���K����>�ʱA�[�����,����ƾ�� �Q�$� ��Ɵ�i� ȴè����O��� �Z�xz���)U�5fzO�o���&�Y?�5M*S���;I�u���O�ƽK���o|]mDx;Ķ^!:w�����>W��fr]�� |����G���/�u�� "״~�_���f��mGZ����?��'�M�O��g�c>kg9�:w�)��\�X咧ou�϶���WQ��h�KE &�6��@	����P?��#����:����(�#Q���m��~�G �G"�� jٯX��~"I�\�/t����/�ʶ8&@a��~�W�|z��� �����o���(2�\ BJ���%$w�rb(���toJ��~G�'���O�?�� �}�~���/��c��������� ���ϨC.��5��ep�"ˌmrˌ�� +�:���ƛMu+Ӓ��h�KEw��m�m- �k�� jO�Oß���D7{4�Y!n������ �8�����>�ED���e�I����W�/����+h�/�^�ۛ{ȳ%��� z)1��0�y�B�)���G��g��w�_B����#�;��_�V� �W���Ky��v��ӱ'�W
om>�����.�J�j`䵆�t1	�G�_���%�ٷ����V7Ac�t��|������2v��d���G��)��iqO����/v�2��O�p�ׁ\0�B@>¾?�g���$��_�~%�}�g}�u�P� ?Z����K㽫� /O�s���SQNX�:(�T�*��}��k�
��isE�;+�k��k���3���#c��Fq��Z�����]���K��%��Q�5	<ɦ~ �
�:*�   �dӿ��� �)G��V<�q��F܃.!^���%Očbdo��B����%���� U?��)�z������>(��?foأƿ�=����j�3kW��/qn�i� k��<��� �9�|.��� V����LO� -��u�_��:9z��R�$�$X�E
��@� ;
ޖ�S�3�#�O�o|�e��<+����Z���O!t���;c��  :ݴ�W��J���bm�m-�&�6��@	����Pm�m- �j����xwG��u+��Ӭm亹��;b����`��¯V?��/k�k���h�u����`�,r�ѱBA�c��3��w��<�����b�����0� �0�)��3^[�+�
��:}��RK�:��!�l֖۱�w�>�j��:��'�^4� ��O�E��G���/�u�� "�x��H�J����O���O������7)%� ����J�i$��RI$��I$�k��կ�u�O��i� �֟��G�:��'�^4� ��O�E���'vt��Ҳ?)h�կ�u�O��i� �֟��G�:��'�^4� ��O�E���Q�b�-Y�5+�R��,.e����g��*�H�2��ЂڿT��	?�b��Z�-��?���1x�� �?���ꇷ���� ��iض�g�&�58�Fu�"%�)�>��n�f�O@���_���U�]<o#%.�n�e�;��5���?���1x�� �?���u�O��i� �֟��]���U���/bϨ�|V����&�R�~�m�X���y��I�6Ӹ��vז~�?����f��Z�/�[�K�϶����I }��Sh1�A��z�vǛ�sns���M�m�����?i��O��Q�� ��<�3Peɴ�Pv���ʰ�{���xKU��5O떏c�i�mso �]N8��#�#�_Е|� �� �#�y�E����Z�ھ���n-��E��Ҁ�L�H�q^@# pbp��ގ�M��=��Z+�k�G���/�u�� "�� ��I� C�?�:�� �k��uN��@����V���	?�b��Z�-��?���1x�� �?�����@���Y���S�5��_�ǈ�!q/�.�<�?y�~v�H*}�A�� �%����ǆo�z����=~xá�*���	?�b��Z�-��?���1x�� �?����MYY�	ʔ�ϳ�Ѷ����3�����@����G�&��R�/�d�?��]���CA��Q~ž.��5)o�_�.G?g֭�$�;%�m���n�9Q�EL�����H'�&�E(�ȡ���A��\ըF��sju?C�آ�]�3� �6~|J�{� ���%b��j	lُsnH�{W��1� �[�X��$�ý��Kk�+�m)���(?G?Z�'���;�Z�|uE{�������I���@;�i.�MM�� �?�=j,�|$��5+8���f��+eS�_�i���T�a��d�(�Yd`�Y��  �k��Jj�C/��W�hV�����.1�ܰDS�cߥ}��7�5�e�d���vԵ�P��:����8	�Ě规�-�FR���/�k� ���,<o�R��4�a5��nP���ʽʟ������Ł�~��
��(DQ�U �
}�R�Q�N	��݉����V�bm��r�լ�'�o�������+��������u{�ß`��7۾���	���gPۃ���Eyغr�⢎�Q�g���k�-u�S��5[V�圦9���#�����υ&��xvO��x}b�oa�U=@�n�Ҿl��� ��Q㶖�AK�j.s�5��b}ᐜdd�,�]�V"/F����T�A��[�� �;�5�um1����"<�?���� S�G�3*��-3�sen~�h�+�o� �*�$�;�o^��Ml3�H�e9� �0d����
���c��2��_Zȧ�+)��G���ѯ�9h�Tλ�����E��I��&������>w�sz� *dgb��M|�_K��	���\�,��Ѣ?��P�����woҾ��C� �����{�(��ۂ�F��$�5�a�}���.�zҼ��=:jɟ�	���]�Bx��6,aV�T�H����H��裖�:��7��_�g߇�~��{�����E[����� ;f������;m��cH��4�|��mj��bz��I'����(a�^��U�Ӡ�h�KEu�	����P�~џ �o�'���F�@��N�nkK�W�u9*������/�_
�K�o�W��V��v�lx�1΄����� ~� W�|h�࿏�:7����Tɶ���Z1�(��Gl��8��u[U��J����~Q_q|R� �U��C������JZ߷��Ǣ����������~;Y�d�y|�{�um"� �K)�y�R.�'z��<��{�� �N�^"�ׇ�|?yڦ� �+������������P���M�˄!��,U�l��vϙ �l��j��beZ�|I�6��~1���A��07����u�!�?礞�=�Ts_���/�+�\�ǅ�-g�[jqu&�r�W,�ݏ�  W��:��}�I��-��2�1���B8�Q�TQx����,�ί�ݼ��jv�YL�RD(�Igq�k֥��(���W7�=W���:��G��2� ��O�F��C���O�e�� #W��:�_� ~S�_�����2x�� -?���u�<e� �����G��� |��~���p�o�<Caq��:�qq�V{�
Ί�V er@ ���O�a��N�o�xS��k��.�5����[�P =	�����2x�� -?��߂�?�~|B�񎁬���R��XRJ���"D(�	���vэzv���j��+���������O¿�:�ޕ��iڵ�+˧�2�8���D=A�Y��~5~��;|D�<e��>&��/R(�M��H@�5�pn�3�<�eNt�w��yo�!~؟�P(�3xn�������3t���N�主U`�4E'dd�NI,00� �Q_�_��?���1x�� �?���u�O��i� �֟��^]L=z��Gljӂ�?)h�կ�u�O��i� �֟��G�:��'�^4� ��O�E���T��@����?c�����b�5k�>Mg�Z�+�Y��f�E,�d�'N7 ��W����?���1x�� �?���u�O��i� �֟��W5hK�;�*���P�o���O���;o%��ìkg�ZM�� .v�''�^� ��;�?���
t]Z�P��1i�'I�d��*YX0Ì�Z��Zn��x�)�}�6Ѷ����������[㍩񧃣�?Y��{3�]V%*�<	Tp��F��G�6���ծ��V��M�md1OiwG,N:�+����v���7�����տ�U���PT��jL��H�x�������ª��:3��n]%��OE~����	#g%���緃������� #u� �+��� �Ix�IH�����xhl%���Y�y�U���_��s�z�������?h/C��cOv�0�^�20��N�#��E3v����	_���7��x�]�<b��mB�g?�*��A_^�O����+}úM���[�Gie����Ory=뢞Mަ�S�/�q��� �?��?�|3���/���JD-���F��s� NI����W���#���؛h�KEP��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���<���_���?�+j��I-~��/�!� &��K��]S� I%������N�6��C�	G� '�/�n?��οV���	C� '�/�n?��οW+�� tc_�E>���2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��y��!� &��K��]S� I%����������)ث��$��1^F;�߇ٟh� �(��<E� b��������~Q� �(?��|E� b�����������¿�6�v�6�i�a�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�ma�S�Ѷ���)�h�@Xm�m�,6�v�6�E;mh�������N�F��h�m�mc�� i�7�_�*��I-~W�O�$�����U���Z��'�D���ϴ���� �q>"� �V�� J���z�(�����^"� �V�� J���mt�?�c_��]�m��(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh�?i/�7o�_�*��I-~
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
(��<���� �s��� b��� �r��e~��ҟ�n?�T��Z��+�#����_�	/� '�?�.��οY��g�	/� '�?�.��οY���3�QE�bQE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE y��+� &��S��MW� H����������*ة�� ���^^/�Ge����[��7��*\� �e�~��ɟ�$���o�� �T�� ��:�h��/��k|cvѶ�Eu���F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@�F�u ݴm�Q@m�K/�c��O�5_�#��+����� �r��� b��� �r��=yx���6g��I_�8�� ةs� ��u��_�� �I^h�� ةs� ��u�Ӷ�p��1��E.�6�Y��R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK����ߴ�����_�5_�#�������_�����*j��G-~W�����Cf}�� �� ������?�Yg_�5�-� �� ������?�Yg_�5Ӆ��o�(���1
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
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
(��B�EQ�X(��5QF�`��(�,QE��6��?���+ا����� W��Lɷ�V� �OV� �9k���7�:��Ϸ?����r$� �N�� K,���m~H� �#��<I� b����Y��ta��̪�CvѶ�Et��n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE�n�6Ө��7mi�Qp����(��F�u\�m:�.vѶ�E�h�N���ݴm�QE��O�e��+ا��������Mɷ�W� �OV� �9k���;�:��Ϸ����r$� �N�� K,��¿#� �����$� �N�� K,���mo����BQK����b%�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m y��5� &��_��=[� H���_�i��+� ا�����W���tvg���H��9اs� ��u��_���H��9اs� ��u��[��:��
(��L�QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`��(Q@X(���EP
*���ei�W�f���#i>�i�d���^������%�h;F��iŸ_��$Y�n�TJ�c��J����tW�� �ZxC�������� ��� i��k���;Q��+�O��W�'�e��li�A���=���.�"F5����o�h����?T�B[1:r[��(���b��&�T��6�=Aj[P���.�>�qy�F�}��7�&v��'����!F�"���O��.��G�����%e>�r*/�kO�;[� �� ����S�i����*���G��vwЫ,7P�� +(`\�[�QHv
(���QE��(�,QE`��(Q@X(���EP
(���QE��(�,QE`����!|T��m�6����v�o��۫�8��̠j�IE^CQrvGiEx��5��?��� ߈��u_�7h_5y��.�P�xa3�]����`mv�ڦ5!'h��NQWh�*(��"�EP
(���QE��(�,QE`��d�y0�&Ɠb�ڃ,p:�Fڇ-�Ex���I�.��tO�\/ފ{H���sU� �<!� @�o��C� �kmO����c�(���t[R�dK{��5� �YA �H��_��2�P��)�EP
(���QE��J��W����{%����K��sZF����TJ�a��Q�)l{��֞� �v�� ~!� ����%���=c�Y�,V���F�����޴Fq���Pq�֢�*ɰQE��(�,QE`��(Q@X(�z��e���_jQY��7<�6��Ҽ{T��4�����&���m�yBO��⢳�H�ٽK�9KT{m���5����T��ՃF�eO��T=�F�y�$W����
I��|N��l�aA��Hȏf~���s�IV�NI��M;X�h��,� k�>]Ic���Ž�l��h����'U��4�k������J�MT'��d�.̷Ey_��h��������$�F�dkDX��q�H7)����kO�;[� �� �ڏmO�^�}�l�������:,��������EW� 9³|ýt������(�,QE`��(Q@X(�ږ�k���y}q��+�I�`��ܚ��g����Z���}��)d?��X� �"��H�ٲ�NR��E|�'�d�v�֚��n��#m����O�4k�Ԋ�� �BѾ"i&�G��!�4��aoF��������E�6
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
(���QE��(�,QE`��(Q@X(���EP<���� �n��� b��� �r���_�7�5� &��_��=[� H����V�ꣳ>�� �E� ��x���+��,���=��� �� ���'��W?�Yg_�u��2��	����WI��h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�����o�x/[���e���;\g�?P@?�o�o�� �[տ��o� �U��/FiN���?;�Կf��+�g�q�� E5ymz���� %wL� �3� 覯7�TwV�>��g����r[��l��̇�\r�5��O��j,�eX1JF7��*�\}��������8�!i�.ɧ!^$p���l\Tyd���'�Y�|(���|5բY.tY�&ś+����5��{o�X��Z�&���e�E��� �F�9k�����[��:LK2��������nQq}�QjK����:��M5����`F�_?W�_��~�g����Ҿu�:��/VwC������� �/@� �}���Z��X�� �/@� �}���Zگ~_<h�(M�m����6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6�}F�-CO���C�4M���EY���������9z1���?8dP�0"���I� ��S� �s�į��t���^��$� ���� �9��bW���*��L�q_�_����6��^�扶���Pm�m- �h�KE &�6��@	����Pm�m-��c����Uo�"����c�ξj��l��k���� ��|�^%]*�y���i� ?�Sￆ� ��� ��� �k]6���8���m� �Z�M^��&y�P�h�KEABm�m- �h�KY�"�&��=&}OU�KK8FY��=����j*�i6�^9񅏀�3y�_7��\Gy�C�P{��d��k�_�5?x�]Ze���%A��#
 �z�k�� ��/~(k�B�I�$ZZg���z��:{���e��%i� _���F��RO�����R���]o�3�v׀~׶Q�{��O|s��H��W��ߵ���:���b���_�?���>V��>����� �o����z���w��׷��՞i|�R��mih��M�m���mih��F�Z(6Ѷ�� M�m���>2��>#\���WZdS0�4�Z�S���}Nr��5��&��/�+�t���Q��r�?
�q��Z�1�7�F���I� ��S� �s�į#	'*�Ow��q+��K��3��������=�-J�4,G�5���#�?K�����oSլ�&�J�t7�ۓ�w�>����]_úW�4�4�=I#9E��IB�7��e��鶩kim��,0�z 8��o���4VG�*�H�UK3$��������/��H�ԕ����зXù`������-|����hZm��w}�+8�L��5�]4i{$�՘թ�-n��~Қt7�	uIdP�k$2��r�U��c_W��G������J����� ŗ����5�Ϯ?d���� @J����� ���{���� �G^�^�~�/��_��_�&�6��LBm�m- �h�KE &�6��@~���K���
	�t�-��-N���냕�>���f����3Ge+!#��	���|h�Ok��$�ϓ� m�+�U� ��� ��/�J�p�r���w��KZn+���i�{��*\�4-2[��HG(\|6G�[��^��I�J�5�tĴ1L���)U��O� ���մ3��Z��i�D�r��*��Si�m��j��6��[/݆�1��J����~� ��R�<m��'�F�Z+��mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m��:�a�.}GS�������p��'��M����۲*x��>��汨I���7m���(=����7�� xV� T�eW��BD�Dq�R�0q���m���~��v��kmՏ٭���O1� �>������� $����2��-��I>��Q�J~��ўѶ���WY�&�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�������;��*<�q>9*��c�������|+��� �Uŋ��Ӈ�'�|�_�?��P>� �u���Z����/�� �s���� �-h�|2�~���Q�����6��]�(�h�KE &�6��@	����Pm�m- �h�KE &��{����k�@��Iv����Ԛ���&��?�k�H?�JW���?3������/� ������ ����m|��JG��#��+���� �� $c_��%������WI��h�KE &�6��@	����Pm�m- �k�����᷄��ܭ�͘��<�H�^���_�:�~�jZ��4B��3uGs���/�O�Z�ďK�ߟ-�ol�)}�{��渱�W$w���4��4�>��W�^����3}����˱Wq�Fp�Ӱ��mp� ���k����v�滎S�� ڻ�z���������.-��r�뜡��r=r}x��+���n�\^\8DE���=��2���P�{{�#���a��B�iѼ>��^�c���l�"�]�f�>���7��u�w�kj[��Z���	�x�Tg�� j}��W{wwgVVmih�16Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih�4����1����Z����_��A����o���o�#���z��[������	� '%�O����ο\��7�	� '%�O����ο]kl?�gS�E>�鹐�)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�����o���� @5�Y�'� �oV� �I��YU���.�ƽO���/٣�J�� \g� �M^[]�����x��o���Z�� ��那)C��]z��2�娙�UW�G���/���<]�*�[I֖h�qȧ!��ۋ~U�(�㧌,^�M.�G��l�l.!B�z匬�}�OC���Q}"N��L����H>� S� }WEnz�%�.Zi��<sEѯ<C�Z����ˈ�GR��O�}������	�z4m�Y���� y���X�X���g����l��2�.8U��5�WUj�m՘T��%~������ ��� �9� �Wε�g�� �� �s� 4���ɟ�/Vz1�#�~�x�D����� �-j׈5�h�Z��:���&�$o�ܓ�����D����� �-k���ƒ�kV>�B-�P\�(<4��A�/?�*��5}��ݳ͡id�9��y�O�mj=/A�m"��Q֯�i�8�s�z׭w����?��fׯ$�.��c����cv3�����^[�8��?�0��,�UC�ԃ�W�5�J��>i�ٵJ��������_|:�?�����4Ȥ0�kr���c�#t�}:���>�i�iv����qis�)�S�}�ᯌ��k�K������9?��x��|O%� �u]g��	�Xrz$���R�Taj�k�{�"����~5��?�j���`��T^^W=G�� _:�_|a����hV:��(g�hO�F~� ͎������o�`�44����oˌ���?��A�/��j��My#��eh,�ݞ ��L{�r��[�[^�4�4��nϦ�I�*Xj:|�Z��T�Wە�Q�e��Ѐ��}rq�_:i^(��uIO���y
<qJvnS�~���_�V�#��Ү�J�;69�����<C� a�����)�M7�,:��|�G��1/��6k[�H5�5e�a&Cǘ��<ؑ�^�_��W�H��s!�| �Ҿ�9��Z�!>jjOs�q囊>f�� ��Q��n<5�˷�Kc����N��@�G9��9g�_�������5��u~�l+l�\�gX��c���M&��5��=�$2n�q��_�� �>�������m�rv� ys�?��>uƣ�]��V2��>���#���F];�ח�����p�G�!�yg�H�;+�?f���Wo���n��������a|S�m�6�0�:U�F3$G;N2:�)���\����G�>�4�J�[y��9�ZD*%\(��q]���Ւ{5�ގI�J��ў��N �_3|`��n���G���#&9u@4���?��{c��/�[ƒ�W��ZZ�c��d�0e8+3!����*��O��f�mo�y��y��+
�%9�(�oJ
1��>����/�/E���[����|�#YC8C�Y��g=p1����^.���m>�m��=�$H��KAq������ ��}�kn�v����$�; 0|��_[���C���YO�!���i׊�iRК2u.�j{��?���6����gC��ڱ�C&:{��� �u�����y4_�˧Ŷ�D˞7�.����*��a>x)ҏ$�FQO��䑳R�B���x�`���J_\�O��J���JI�(I�C��S��_��%𿁓N���w���r��ɏ�U~�k�
�o��:������B�Z��ԕJ��.�e8(��H�¿���6�Z����H��yЫI��_3~K0Hc8�������V�:�Յ՜���adb*xa�Њ��cEUUF �|_�B�P�׏u-O��6������UI@��)=�l�{b�g%*z/�ã'Q8�S���`_��\����ݚ�2��J��"���Gc�Z��?շ���� ���E�-��I���n���&ve8_P6�~������k�Rr�yn�9$�f���ߛ�t���^��O���]gV�.�{���ooYd.��O�z�ٿ�I���u���,5� ��%���6�u
7�P���#��a�sJ.��ǥ�i&����$?<]�M�kV7��<�k�$�1݀?���^m�� j��K�'P���NwB��ãc�~�k��u����[Ȅ7B
�_��6�����T�Y4ߴ�����
�N� ���Ss:,Z���n�>�8�p;�A�@��c�]s�V�;�)1�G��_OM2[�$�6��R�ǰ$׫�
o�<�G�n(�>.|]��]�!d[�Z��Y����~�r>}�x����d����q�
��+	+K�X��c���'8�����ό5-bvb�JD*O܈p�>��k��MH|��`y��Y=��)���'��s�kc���0�w}N/�����1��ºƢ����,y8�c����_�q�_A��J�y/�i��cvG�=�~��kվ0[���� ��#�27�p�E|�Z��9�X�:2�A ���k/��({jW{����|M���|3����}̄����g��;��Rj���x7FՏ߻�I��a�Pk�_��8��5	ą�m�ڮx����'���&����vcB��{#���E�O��k�ˬͧ��Y��%!� ��r@��՟��`�c��kZ>�wqJ!�\��W ���;q�u?��{|MxFeA�S���?*���^���ýOM�_2�v�g�v�H��kR^Ǚ|[� _#X�~ח����������/P�J�mcEf$7^H��ѱ�c����־����,m�-d[\F��"�ea�#_�K��M�0X\�G��_)��ޛq���?�+��B�/�����VU�v�]��q�V�i�P|���FJǎ������w�� AZ������������� AZ�����ş���_Ï�՟~|0� �q���� 赬���V��O�n��ˈ *��ݱ���H9��AZ�_��$�2ܓ� lּ��o�K� �������&L�CRF�K8e��Sߠ�k�Ĵ��{m�>��xS\���]v-'I��sr��ј"(��wp���_R|=�e�xN�;�[�:ֿw�k{��6����H�r>����᎛��A[+@&��uxˆ�������]�*t�5��.�G7��E>���C+�O�����s��Ě��>�o Hl���l��zd��������� �G� _��\���9�8v��|a^����xUo�7��7�t��|�;QH�k�k�o�g�IZ����-sa�'~ߪ7�'�ܳ���>6���?�|G��m4!!���ϱ�d�H�;��k��t?��?�-��+�?k��4?��?�-�����^_���V����{,x������	a5ֱ�.��rǅ���nx��*�c�w�톋��K���%Π�{�O,�Ϡ����0�rRW�4�%;]�&x�C�F��l�R^.� ��Sy�ǨpH?�⺯��5������s_�.�&������ d�#�:ֽC�·V��7>ѽg�7w�U8�+�
�JTj4��M�?H!�'�$����Yzy��"�bTd����1�k� �~��,�|'� �W�Nx�_�t�YwZ���pD@fO�!~�k֭QQM�<�QuG�|]���./'�|#7٭c%$������_���ث~��_�~�]�n������E�3"�U��1$���^���f5b���'�8��	����-)o
E�DP��� `
�h��jt՗�j0�����^/��[�&RK�kv����A�O�H��Fk��t� ��T��Ԁ�w�r��*}Gp{�ʾ{����?i��O�	Ψ���'�I���Qs�mRR���@]O�~4����/;|��������]�9W��P�k��+��Þ*���$E��Ѹ��v���{=�⦬ˌ�v?8�ym5�'2h�dw��	������/�5�_�� �ױۙ$��2��7(ۺ0I�:q\ω��c���_�׮~�2m���&>����?�+i(���=D�"��������'�Ǉ�?km�]�)h��Y�4*H9<W��p�|Q�l/4� x��O�v�mơ/����aҾ��
�����r5�� l޻q:S�[���A�;=�Ϝ� �fx����o� 3�U�g����SᏇ�.&����t�N����������_���$�� �o���p��W��G�V��� �#�H���� G%|I_n��_�Gu� �?�::���k� _�Cz�������O� �ow� a?���k�� d��&����� @J���<i'��\�YW���8�N
�	v�G�+ӕOgIK�~H��\�⻿�󿌟���w�'��T1��e��8����T�|������|_��nיkh�`������ӿ5��Q�eDY����cG��K�l�����@�T����}�s���Z^ͨ�C�z�ş�ߍ�N��&��4+&L�FN���x���_
�&�|N���`_���Dwv�����=T�?_J��
�WW���|��~������5��O�N���E���k"�'��w ?KQƣ���O�u��=���:�}�s�ey�o�^>�W����ω� ���n����2�-�5#��⽆��qSVe�N;��Ν������p�I҆-��͓���&�?�^���N�^���!�7h����'��=�3���<x����1�пey6�R��e0� �MyXX���}�CѮ������� ���:�^*�� ��Uh�r�yjf@1ϥvZ��e��T�f�vɾI��� =MiW�?�����S��/�`����;d ?A�� �F�Oe7�4��'���F��<G��a�.f����8m_l�=]�<��Z�	��ڶ�����Auk�Ļ-�q	# 1`K�����<'f5/���n�C�p+�B��MN.s՚֛��a��n��x��/�&�-��淇l�G1/��c��8�y��>�b�����yb�Y�Pn-3�GM�{��П"��-5� ]����2��#� B5���[q�i��@�Q�?~3é����j�O�OKغ�ԡ�����j��2H��)�D�]r1�{_)�n��W��Ok?���O��ϨJ$F\d��w�� �.a�h�|r(ua��k�o��<�1� ].?�V��v<�ne�w��� 3��5��1�� ���~!�xo�7DӴ�E߉.t�$�yN� ��N���s��>� +W��G�m��.�(��\�ډo*i��;@\��>��y�E;y���Rj��׈�O^�ާ�\_ܓ�2y��¾���|]�x�������]6�8�)�;��e�R{���;V/��~Ok�j���������rF�ݐ��_��~�Ӿ�6�ᇅ������d3\��g � �=�\4'.c:��<G�Ƈ���dԴ������f��5���ln `�{b���g�?�k�� �c7� _F~�_�#�?��-�����g(��� $uCނ��� ��.���W���5�x���33ŗlr���f�q�yO�O��C��o�{��d��S��c��G�x��f���7�]?Q���H4��,�vP��POsU?j{t��k;(/�,�Н��]U�^y;���Qsr���w�� �^ �{�C%���n�3��׾���s��_B�����l��[���jP����qnd��� ��c'���_Z�ʾ�m<���T��H�v�)�U�ϯ�(��Q8_E���(����^��kz��b�xgl,�N�F��
�'��|#��{��������X�q1H=:_a�}�������dH��B$\�n��z5-m�B�G���A���kP'���"H��.@p�����K*A�#�q�,�� ԓ^+�$��j�'Pl{��k��޷�|)�������#���g�9��Z�9i�ۥ� τy���o��o�_�桨jRi�t���� iH��/8ʆE�#=���f���O��Zյ\�pbRG@1#��|�_m������|+�$ԚGti"�I3��S��9�
�j���VtV�;rh��~ x?U�?����e�g���f��	  ����>~�W�Z����Oq���a��'��s�G?�@�ry�lo���G�'���� F=x��3[�r�!���!��(�p�����lmR*p�[����x�P������h�.6��6+8E!J [�$��Ҿ�Vܠ��ל~�_�Gu� �?�::��+���\<���S�*�τ� �{���y�c��pѵ��� ˅R[ljGB9��y}u�&� �7�� ����W6
ri��Q�^N)[���� 	>���z��׊<K��e%��8~�<�_r��dP �ך��}���I#����QO������~ƶq%��2oX\�.=���e��<Ï���Ɵu�x����r�`~��,����ʒ}k��4-��R�'��q1�b���� ,��mw����C��+v���]��67����kͣS�O��ӱ�R��'�o��? 药hڍ�ݜ,�q�S*n �2@�^��x�g��	{K��8$������P�������� �C�����J�v����CCHZ�?|���e��-cObַI���z>��?
��o
���Y[�Z��2�$[�6c�6�:�9�G W�~�wo�2�I>U��>�P� 2k�kՕ�~�� 3�M����O�^,�|I�i��Z�[)�"� �<n�9���,kz�t��o����T2��k�dFP�����Vw���*�&� ���B��xr/x�F��v���Iz��@��5�P���g����QQ�և�����1��^�L����{�Gm ��Qԓ���W t��^k}y�^	��&����$��:��9�Ell`�lവ�`��qă
� �G�R��?���B��3�GA�W^>�j��3�.x�Տ�?gO�W�4�}[���-c�a�o�4`�Cz���>�>�_~ͷ&�������I�>� �c�+�j��Mʚos��yf�7�wċ�~mB�|���kj�Wǯ`:��"�A�?��|KՒ���\�#�O�c �¨Py9=['޷?h�I�o�W��!k=/�%�����ݑ� G��WT�������|{�3�
���"�M����ʨ�r�{O�� f7���֟�KM�*�o>݂ۆ��*⽾�>���>7��S�zN��O���E�xn�#�0�d��
O�q��
�{���m�u�T�6���,2�>���R}�r?Zbb�˵�2�.i���~�|<� �ß��� �k_���	���D�:�� E�^���'�G���S�ۜ�+�>:�\?�H`����a q��GY����ҽJ�<���v��]�R��"�M�l���rbj8�%՝"�'~���z�Ǎ_P���{%���1���6p��*�`�����M�Gğ�$���-A��o�O��vP9 �O�J�3ҼO�ĭSឹ��<���\�H~I��؎Ƿ�"���w�[A��b$�n<��\�a1X�N?��}:TS�*B��EԔ�+�X�'�4x�S�W��v�W�����O�H��"` 0 �+��� �;��5�+)F����J%���в��$C�^��]r�Q�g�+8b=v�q��?ޜ�(�.����f�����܏7������cM��;�~T��R�OF|u'��'����W�_��{y�-z��"Հ���3���>���8���oX���y�^�e�������?ʾ���t԰�K�Ƞo��i��ʏ�EsQ�y9T�t7�j1J���*|4��\Ykֵӧ(���H��@b�۱�G�_�C�����["�cB�̀*�(��fHȭO�j�g�K�3b�_��?�5���ˏk�~�j�g��f\w��}���ƫ����nҚ�S�F�&��?�k�H?�JWڶ7i�Y[�E�SF�/Ќ��_�ѿ�X�ߤ�%)�7�~d�]ܟ��w�o���� �F�W�U�/���)� �����k|7�o�FU� ��_��|K����3�M���$%-���?����
��O�G�>?x�-!�Y,l����)���?9��y#�`�p�5���M���sil��"	��ƽK�=�Q��.�@2*���w1�C�zs��O{��o8�z|_��W���:��|5.��_Qk��gY\*>�# R1�<g�g�_�*�=F�C�Uɺ���pj2��� H�O��#�GO��Kn���`q���T?B�W�e��W�۰�mO�߹�EX8��O�jo���c�q�W=}��2k'௉��g�=�g�r�y��-+��@񯛿ioK�oM��!���p�^�7�?/��q}��w��z0����:_�6���m�P�I�����ٹD�02��0_ g�3�
����<Z���i�������j[~`�.��`�����c�8�nXe�!=�:��C���~ xm�a�cGB[�vH�t�9\� W?��ϼ���?oij��#�߇�|K�K�{6���7r]�����Ǩ5�6���W�bմG�Cwnd��FB���|r0�z�R����j���X�[߫�6��O���u��~xgP��=2�S���c$����1!H�q�v�i���&�"�TZks�O���,��uk-I��0D��9p[j�"�^;~>���W�}� "��� _�� E�|�^}D�RQ]?��|Ћ}��O�)���xOM�Ѽa���O�-�����$cj�Q�zW�xᏋ-|/�i~(�e��R��^�ڌ�4*�@y� u��7�I/����� gj��T��5ԓ��o�K�|5�t�?�&���>Z���e,��7��W1��I�t�����u�V�AS��A�W�~���&�@�:dd���J�+�B���?�!^e9?ik��N�%��n���1\��<q�|=������b��B�~iD_�� ��+�o�c������bH~���!<����]��:q�wg|��df����~'kqX�y-�R�0i�2��  �r����zշ쓧df�_���+�>0�J�?�F�=�¼o�������a��4���3�
����N3��=[.�G��d|S�|R��\�7���e1�[+�2��t�T��G�OÏ��gĭu<��B���^�C�c���� �z�m�,^8]�}�2s�n������i'�>#i�d�e����Rx�����MF������|˭Mr�H�u>���� �n�����e ]����9V��Bk���������E>3�.�EYb��̥��FF����W�����S� %K��a�mO�����Ň��� [n��7�K�Y��U��d2o�&+�p�~n��`�o�|^����eu�lJ=�.!�s���z(9�^�;�ڛ��{`��,{�2@�~��O�^������#���.q˷v>��Tb�A�=VJ������������y��}F
�$$H�9=s_@~ο�+�j/��r_%�"��f�"���y �q�A�'�/
C�sS(�w��;�,���|6�6_�76�u2}���Ya�(�Q��cJ�R�������}���QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E̿i��6ߊ� �)j��G->��~ӟ�m�?�Rտ�Z�}k��[3��$G����?�R�� ��:�u�ȯ�$?����?�R�� ��:�u��� �O�(���2
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
��?��z��zM� �Ӭ�� ȷ�פ����� ^�ҟƽO���?ٟ�J��� \g� �M^Y^��3�^�?�� �)���� ]���g�TQEzǜQE |��a� ��?�� �+�:�3��� ���9� �W�u�O㗫=h��?B�� "^�� `����Ǜ�����6|���cw�������� ȗ��>�� E�|��I诤�Xԥe";�㹌��S� �)��oğ�Ʉٯ$s�	t[����M���<=up�*_Z�޸B��YO8�Z��'�?�k�����_0x;�RxO�Z^��l�R����~##��-'T��4�]B�U���1,R/B�f��!:w{�V��=6>~���nu+�n�<o-�̧t�O��w>���&���� ��N�y�����P��_d�v��;�=�w�T�=+�Zn��=WE�y&��˪f5f�w��9�q�i���R��}��Ƌa�>[J��9F�2���j�[D�W�N�5i�-�5 7��%����;��p�Y���q�h��i�[[�!.2!A�����?LW�1����E�������#k�سf Oֹ�S���׽������}��h���|h���j2�Ί}_R{*�����Y˦��T���&��H���VAc�y$g�z�����O��4�6!�����o�v�$��>>� ���?��� эX�#�(����З2����N��a� ��g� ^�� ��m|e�0�Z�� �i� �
�6�h� ?�S������f�G���Q�@��wLe��1�o!��G�	�=�>���π�{o��֢�mq���e4g����~�����Z��߳����	�ߜ����S������W�Ys��3�u���������z��i&0H��&)�ї����O��� ��MJ�fC��ڱ�C&:{��� \��CR�n�{�쯭䵻��I��R;^��#j��}V�X��k#�m�"�?�7�Z�j6��x�$���_�Ӷ��~v,I��?�|�l�n""F��H�W���W�߶��X�{VU&8����8�e� �Z�e�I{���1�
ݏ�c���#W_�Z�VM� ���~�:��d�MW��$cj=��TzӜ
��~6�����@ז����g�u�� x`�5�W�*4ߚ<�՚�s���.7��U�kK⏴�)��O�~��u�ӌ�c^�Xz�����N�;����� ��B�$�z/��ܫ��ch�L����
(����O������%�?��ĸ���^!��?V��= �%� �X�_B���n��*�״���9ܿ���W���]£}S� �z�JT���(fK�c�3�C)�dS��� g�E�/��P���5E��g�#}
�Ϩ5��;t�����?շ��4��:��gb����ޖO�m�5�O�^�-�����I���u�߲?��ڟ���� F%x���� x� :����G�O���� ���/�����b�z����R� �}��ro�k�O������K�A�_�ɿ���9?�7֞'���X����� � �կג� ��_B�H�k?��$�3�4����>[W�_���Z��y/���������7������%�}Y����x[.��g>���u������7����_�u/Gos��-���*�zt�ּRX�	^9��e=AE{o��ht]�r��UD,��\�~ ���0����S��4cu��;����7���]|S�.m�R�C4r�:���q�+�� �7� ��� )������ږ�m���׳-���%���&�eF�����L�<3�>�5���o}>�f[�/����|�lc8�ھf.Ř䓒k�3C�-<m�5d�l�IQ2l}�#�_ k�L����)Y�'x\�IɊ��W�t�������=ԑg�.�O�a�G��� �C󯥫�o�^6O�B���.�l�\�`��G�B�¾�F"���!��G�wQ�54�i�|�J���Ţ�dS$�'Y$nS��pG�A�Lϟ?l/�xk��&� �V�`��� l/�xk��&� �V�`���g��#ԧ�8� ]Y�0�٪j�Ѽ7���~�og���F��cP���O~+��O�/��[=f��{|ѓ�ȇ�!�"���~�~8�����^Kna�4 �pc�8�N}A5����^xoX��5Z�Yr#z�����ڷ=*�M��ʏ%JVKC��~(�����:M��)�������Z������c�_ej3��� ��`z���������%R*H�\+
(����)��?��� ���0W�W�~ӟ�H�����͈��~f�?����}}��,� �+O���� e��k�o�g�IZ����-a������lO½Fz�x'�}� "��� _�� E�{�x/�}� "~�� _�� E�k���~h��E����}��7�I/����� gj�R���� $��_��� ��Fi|�R���?l��=� _r���}K�`ȷ������ @������e�tGM?����Ͻ�� �0���xE��¿l�oxzߝ�k#�L����+�~� �0���xE���lB�ՕI�L��q�'�ɿ*�1�_����o#�+���l�n�b�d�����#�־�_��/`����}�� ��&���;�h<u�=6�$wkot��eP'��~5�1�Z{�J�Zkc̵��_T��.���S�5ڒ^Z<������U����i�����f�_3������8��kݫS񖛤��J�ew}OQ�Cnڪ	,� �^ںcN��Z����{�QEldQE ~vx��FMW���� ��z��� %"�������ɼM� #&�� _r� �f�g�K� ��}� `�?�du�`�8�?��q_�_���xW�u� "6�� a� �޽ּ/��� �I� ��� �o]���?��|?���g���_���$�� �o����J���o��_׷���]���J�n�?�W��������"��h��#�� ��|E\���/롽����}u�&� �7�� ����W3�b]7��{q���H}3��?��]?�� $����2��_�{�Iq�SU%-n^>�b�柭u���w�9�� �ߩ�ʱV���}ae�7ŷ�p\G�g\�ȼMЌ��x����߀^6��_��}�"�\�yF��@�A���3RO���yJ-4p��칫���mW��bȌ�Z��3��8��O����&�������������ۿk�>i�}q^�Xz������F�u#���7�o+��.ò����]*�8�I-o����{�QEndQE ~|x��G�� �B�� F5z�� %V?���W����G�� �B�� F5z�� %V?���W������ #��|2��O��῏�s�s�%�vJ�>�5�䯎?jɤ�N��f Ԡ�dn۔la��A�k|W�g>�#��c��7�T��ѷ?���������e����o���+����� C� ��/������ ���9� �W�u�'�6u�C��������+�l�%�/ ���L�(�X���󤜦��w��
���>��?|5,�.�m�'�٭x��� ^� ���+�|?�� b�Zv�9����]�����?�������%z߅����^���5~�|7�[���@t�:{�'����~��?���_��:ߧ�sZ�/�/��i������(��S¿k��������}}c�]ȍ�� �@路��ƫ�Y]����������� �� ѯX_����� ���kw���� �� ѯX_����� ���k��� ��3��ƿ�����g�� $������ Ы�:�3�_� �Mi� _S� �U�y�~������f��o�$~$� �u� ы_��?� ������� F-|-\�������(�� C����D�g��� ��׼u��+�� d?�5��� � �k^��O�x#�7���L��b�]�¨�$�^�����_�8l�6�w���~�3���泤X�ƭ��y�(��$�c�J����:]F��9�>�3b&�-� \C wl�F+��{�Y���5��p̀�"�	@��?��k�ٗ�Vz� ך5�S�1�)b�?��rSn�j>�|���M'-Y�?�hz��x�V��ծ�D탔�a�A������k��� ��?�� :�����G�'���� F=x���~A� ]��$?�o?���/����R�+οh��#�� ��z4_��+�h��#�� ��z�����y��#�*���M� �oy� a?���k��7�I����d� �#�0��?Ttbv����j��+�8(�oſ<5�k�L7�F6��h��z����+�u��~������Pc��l����1$��9�W�A|o��5k��-k�[1��x�p�\��:c�9⹏���v�8��mAk-:/��2%m�*�铓�����奥����橭�'�.�"|C����Z�]Ñ"��
��H.�:����0;ׄW�_���A����J�v�k%�#z-�}� ����M� �oy� a?��k������FO�:���#���~G��~��>��� %[������S�� �X���}�L��� %[������S��Y���}�^N���W�zu��/F}�_�X� ����M� ������� %3�� ��� C5�+���1��/��n~�� �X|?��� D�}�_��� �X|?��� D�}�]t?��g=o�E����(�k�jӾwIw3���ɯ@��� �i�������L�_�� <A`��.�B��۔�DV��It��Z�ڲ��qˣ(�H�7	�A؏�L���m����>�_�J�����lo��O���Ү�W��N\?�>� ����A>� ȃ���[� 赯Ϻ��w� "�?�o� �֞���O���t4QEv�^E�A|#������mj�HX�>3�L�����޽v�[�Aa���s�#{�#uɓ�8b{q��=��Ƥc4�.�BN�t>�����H."x'��<r)VVA��t�J�G�����KK�[tsB�YO�"�����KC��O"� g�
�MB��u�1��ƾC��ýg�ήlu{}��a����_U?��+̔%I����g�C鏁?��>Ɖ�m�]�7G2�%ҁ��gH�ҷ� h����C��dGǡ�3�¾@��>���݊��8�|�~ ����_4i5"��KI�c2����� �k�Ru0�O�9TU:�Km���O�?|E���w�MWC��$b��HR"��G=z���})�%x��Q�ā&g�]�c���{���8uI�]�k����K�o��f�����R���-�2H����i��X���� �7� �3� �W�u��i�Ю�}N_&��r���U�<]r�I+�s�����wJm@�tן�-go��n��T.�d�8�_~ѿ�X�ߤ�!+�K��,m�WX�dU�J�ddv<���� h��,Z�����2�ѿ���۷�w�o���� �F�W�ڄ����U��0�5�O���)� ������M���񫢛�����l�]�_�?7���+�Y��O�z��_��6mb+��kq�Y����;dN��"�eе�GN�J�kq$,�,Ez���h|�
��Ee�G�I�1 �m��*��r�I=��G2�kt{$� |]"2?ōiц
��A��"���c��� �� ۫�:l�$<�2���p��J�{�q*��3��S���]�c�� i��<�o#����noO^�������Ě�Ą�%��ؓ�K�_|�GŶ6����+���B�I]���Oc_|P�$���{O�J��wL�e�S�\���m�� C������H��u$����h�ub�=ʺ�~Y���Z���w���G�4�erR�a��h�����&�����S����U��tG"�R2��椼���z����袙�6���m���i�>��1<��� �CC� ��� �ھU��� k���?�� ?�-��J�� _�Dz�� ��>��� $��_��� ��wU��� �K���� �ڻ���璶>I��� �i��O�-x֕� !K?����e��?�i��O�-x֕� !K?���������������Ѻ���u����Z�S�I/ff>��~�����?'��%����Q�Z�/t����tb�(�_������o~�1�>/iY��a� ~��ԯ�f��?���9�}|�� 
�Z��տ�ώ� j��*_���j�}6f�ԭ%C�IQ����C������t�'˵�~�.�p� �<O�M�K�~; 9c�(&��z����S�m*w}�C�
6��Ԍ�ǿ�O��/�q���}�_��_�T��xC�ں�}Fsaw���-�?�=�� �����j���?�=�� �����j��?�����G������� ���� E���� #��� _��Ů���� ��� ^��\�?�s������G���_��S�_/����+�<���( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��(�?i��6ߋ�)j��G->���ӟ�m�?�Rտ�Z�}+��:ilϸ��� �r~%� �J�� K,��޿"?��� �r~%� �J�� K,��޶��u> ��+��(�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� +����M'�~����n$��8mˏ2Ge!B�S�J�m*�\�.����;��)c�Wr��2#�5��8�;� �{� 'n��cQNQq�SZn1jL���߁�&����"� Q�m�� y�?Dތ��l��}C� ���������� ��8�;� �{� 'n��sS�8IKCyՄ���z����n��\Eu���u?B*z�t��	�m�i�/ٯ�شR��w�H#�9�z�^ڜ�}�*���Z�e֟}�guE,{��R0FA~�{h
���k���%���-���+�l����`��0���^�������3�ÿ���v�� �Q� ���������� ��<<�n��U��F�¿��"��%��S4VP�,h�n���#�r�7©~!xv+�6/3Zӷ4q��F~�}x��Gz�|3�g���"�4}�w�T��S>��}+�������HKٻ��ښ-�x�F�Tb��0T���ƽ��M{�-gǨ�L۾�pH�OR�>��_V�����y+\j�Z���e��?R8o�5�G�$�5eou�\�c\E���y�k�4j�~�:�V�׼�*����O�-º�����ċ1��g�_�B�ԃ��׸|,�ï�bd�N�;�n�&YO^O8�ޮ�?�߇�n�h�dV��M�����'�ǵt����'v�IK�D���?k�:�/�άl�D� ù]�/��������)�s�8tu�A���3�^�|m����٥��[!��2��q^am�'�*�3ͫ\G���K�}����7���q}nu���JG)��ƞ>��D�|9o&��MK�&�C#&�/�,F{���<����V��IM��S��|��o�s�ƾ�������c��C��JB�,}X�c�I5�|@�9ᯉ$����^��[�G�.�C�APqZU�)���ΝX���G�_|sm���6z�強6���Ȑc~H��NqԊ�O�~+�ƺj�isi0N���p��$_����ӟ�q^������hn�i#!�uU�H�UA���eP� �oF2�m&eRJR�O�b��u�㇉o"C=�����ڎ7�d'#������}Y�/�^8���G�a�R�h�,~̝Eb��[�m�k:���u��e���dg����{���gV��|�M��,r��+� ��*Q�J
غ��I9lx��E�����(N�473Aj�]�	L����� �ڽ;�Y�us����Z�-Ƥ��n0�s����1��]W�g?xV�.��mR�3�R�H�� )�A�N����Kٷ9n� QT�ΔV��~ �2�����nN�=3�ϗ �[�?�k�� x_R�~�q���mw	�V�ó)��:��g��/Y����C~���WO�\��5�s�h�U*܋��|3�?�z�Ýc�CH�T��5��tS/���W���ajSX���3n��0%k�x��0��������4���X�Rs��qQ��FO�]o���C��	uc����9[���dS�§� iB�[r�Y)�o��f���uu�O;M�=I6���k�q�$��`=I�[���IE$�9[m݅QLG;��v�<){����e�R�&)*��BE|!��	�~��4�V٭�?�^̧��Z��Ox+C񵐵���u�����/�k��w�ΚU�+��G�^8���I���L�mx�n�U��/q��W�/�}����[��c�m���f�������<<��u�[�9�qQ��FO�]_�~x?�7���j���n�[�u>��O� �ӧU.[�9�o��g/�"�V��;�w�����*���_g��p3�9=���I���7���V�E6�\.�A���
ԮW��1��{�{L�t��c���,{T���a�}kyE�r@�I9sL�
F�#7��{췯��&K��V�����`�_r�\��^�� ���������� ��8�;� �{� 'n��sң:rRV:*U�D�:�x����/�����`}�G �m8U�$���NI5����s��$��m�#�i�[U�L���:/��=:���� �q�w� B��N�� ��uhΤ���Ռ#fxg������1Z���ˉ�p��:����g�3�G"���������;s� �+����MK�Ӭ"�,�cE��Tt$���馜`�.�=F�'%��G���Is�~�Ěu�I�_?�?��[�O9�V<��H�ω�#�"�2:���pA5�E41������U�uX���W����f�>�p�Gis����|�T��X2�� W��.���tզy�� k-sE���4�uƍv����+{����|
��u����v3�4oZȲ^-�'���B��1� Ԏ��i?���4ہ,��:��|���� !��ֽ_L���[��-a���a!�"� ��Ӝ�����f�W�տ�����Y[Eok �8�`*���_9~�_�.���5���P�1��ߎ�Z�J����Ug����z�߇?����]6-6��=oN�m�f��,k���<{q_Ax��{�W�.�]9��9y���}Ws�5�c�'�*��d�mV���P)��E?�rX?u�2�NK�G�<k��x�zRh�x��k)��N���@A���������,�k2]!\�3}�d<����0*߇|/��OO[-O�O���q�����ܒkFhR�"�n�E*˜d�]��*ov� ��唹�[#�o���v�h:]��7V�4��1d(P��N��}�� ���������� ��8�;� �{� 'n��qK9Iɵ��ЌT{����^�k����2�kiB���6�Q���9���cs����;}Bݖ�Ub�Bp1��L�W��g�σ��E�h�?����/ڦ|<3���U��)�/���t��\�����c��a��Uhʴ|�zRT�~����W����U���ֽx��˶������2�v�J�?�~� н� ��?�r��g��/���� ��S�I���J���I���n�o����V��V��QG���� d�x��v�\GƏ���᮳�Y)��Yb�uvF�{������AT��.p���~l:4n��U��}+�>
��_�/w��Zu�����jW��e_$`du�ҽ��� �#�羺����C�.,$��ՁI��MQ�����%ٷ�զ��O�C��]���rS�V���JsZ�G���?Xj$����$W��s�[���Fy�g9�{ך~�^"���4]"����᧖��B�n�I��_F�q�:W�I�:�=�Fw��Y�b~�q�� �Jڵ9N*����MF\���~����u�[�~�mos�Vq�S�o��=G\A�i���?� �^� �۟�9Z��3�����F��K�Օe�Tπ���9�*4�N��:�#Q+;�]����t-"��)���Y��6b�P��Оx�����	?g_��H���1,O�n:��iM� �q�w� B��N�� ���������5����c�?�lu��ZV�0�qkn �p^6RG�:��?�m�G�E��/_hӑʻ����r�U����Y��-����N��K�w��Tό��t'�v��%�"���q|��z�^$�ޣ�-b�K�m����#����������}X�i( �[LE2���{v �+�x@�ͪ��i��*��l����?^g?��ɦ.���*O���#*?8��k�P�x���BJ�Gq�`jw"O�Ǩ0ڲ��H���aA�7W����i5ߌ<Y#O�mQ@� ����c�I��t Z��_�?�t����:�>���y�/��>�����4���9&�b���(��©kՆ�c-�w��jY��>�O�]�+���-^MSX���}"���i�2`p�Oj�s[�ܨ���>�.��V����M;ȹ�,H�u鿳7��|;�)M��v�]�Il���T9*�xۏƾ�� �q�w� B��N�� ��?�~� н� ��?�r��Н&���R�*&��*��2��yt5��s�k	��F��)�V姚�3F�J��t�n���~��x/L���k���1i��}Į��1#�C\��3�ÿ���v�� �V��*��^FT�o�������,5�E�715�43[����pA57�3�ÿ���v�� �U��?��R�Q�����Y�����k�p\����4�M���V�j-7G?�Nx��O�k}��u��$��x�+���^���}Ϫ|���]_�h^u�ԭ4�}�u��rN�Oj�� ���������� 㕄�Nrrv��5�����E�� �'�h�uz�ޙ�pё@*;򧧷�{�<)k��-�]G�HLl9V��k�O����n���)�?m�� ��
���yy��?C���\��?<�a��O����V�n�\F~V��"�u=�� <������Z]KG�K��E ������;����<ib-5�6B�|�C!�Ve��0���|q3<wZ����W��4d����*A�,��В���7퉨5��������1�?�m�ẻ���-Rx�������Kyo١>��$���'<�]�~x;�71�Zi���3���5���������*�|�w��)+r�YQ[UԵ[=�K���l�igp��ժ�|U���>6�?���'�Z����h�Q�#�����5�ң�x�w�Z�Z��5{�`���h���9#�5�~ξ �������m�K:Vڪ�>\��$���_��?� �^� �۟�9G�3�ÿ���v�� �W*��V��Z5O��Ȳ��0ta��rp?����	��;SU�&k9�n�(O��ڷ��� @���:��"�`�$���c�������T�,�X��\����� º����KK���W��_{xC�ZO��;�cw	��2� ���_Џ�G�����j�ޗ܊0����=���y�߲O�e����0���,�2ZjRN;�NP����O�+�~7��ͧL�60�T�VF�2��'��s�6���P�źջAo�ao*�b8���;z�{s�����)㸇L7�Q������P�.}���>���V*�y�,v
������o��6��)�'H�1�v���М�f��;�Y��]�̆7\��#d{W��8�;� �{� 'n��:��Qr��JQ��g�5���]{O���?a��ᡲ�)Q�� ea�Ec�8�;� �{� 'n��^���#�z�����f��m�K���i�3�r^�4�N��¬�R�n��(������VhzF��Rߋ�<���ѨB2öKq�C_/��3~�� �&�Y<?�I�7�nI9?�ҙ� ���������� ��KRRrv��hF*+�[�q����I���6��E=��?x��pA�/�Y��[�?])���и����0/�rK�{��g�|�O��{mSL�~�}nŢ��S��A�<�TZ��o뚝֡}����R4�I���s�p�
�ԇ.�=9F��}w�*�������E�CP����ܰ��C�u�zWG� ���������� �cM��='P���м��iVh��e�mu ��&#�Es�&��B�ԍH�t3� h�i�O�Z���5��D�}�{��I� �������x�X��5M�W�e�T鸀 �\�g�8�;� �{� 'n��eR��K�ظU�#ʏ=��|Aa����E�\����H�0v��s�Q]��E�]_�,�푝m���`���d��jE�:�=�T�<?�у)�m��J�Y���!�5�)���e#z��%R�#�O��MF�2��گf�U�C�������]
�l�e)"�9ـ��d�8<��`տej��sԴ�c�oi:�c���k��O|?��1�5���1���O��V4�ԋ�Ƶ*�Kk�,|u��`��V��Վ)��Ek �mNDg���Sׯ���!�Dq�R�+�7�~�|m����٭��C#��9�^ao�'�*���է���K�}���C��N�صZ<���A��
|N�h������n/nYB,�b4;���p:s[��7q��h3u<���C�J���A��3�C��v�X�D0��0>����MhWl�����k���42[��J���Y`�G�{'�� �_t;�/S����IM�k�!ʀU���9��{� �~�O�׍{{i-������� �*O����~̞	���]<7z��ܩ�J�������9�N�ZrѝS�NkT_�;�x��K�x�Y�e��SN�3��%��9��^3���zM5UcP��� 
uw%d�8ۻ�R�ih�} ���F�u���S��F���En��y�z�5����6����CT��.��=����	�r�����W���R���� �l�q���l���_�������=�ί"��ԥ����kΧ��M�=��քӺ��x������i��f�<���%�Nv+m��SҾr��X�(V$EH�v�P�� zW�x��_�V��Iv�}�4�sEa2�y=HWF��`UT��˚$Ӭ��<���H������[�_U���-#�*�*��Ws����[��-�rGj���iO_��z�����?wI��x˵�.[̘�L�Q���x+F�֛��g��H�,~kǇ �r�F?�u�8�۳��_d|;�SW�׾"����E��k�1ȧ!�p{f��X��~$xz��E��+��$c���w�3����������;s� �(� �q�w� B��N�� ���Bt�&�:�ZN=�@�T��O7�]��Aw�������1_ x� V�]�ƽ�Z��k�٥������+�Y>x>_E���3��q����L�K�7n߻�<g�� ���������� �u�ά���  �U#M;�,�֭<=�K@���`�YYYj��d��,+��f��%�'Yca�t �=A�8� �q�w� B��N�� ���+�'�:P�t[O�Yi<�1���S�$����FP�,��J3�2<C���Sq�*x�J������cb���|t>�=�|�c{6�}own�;�$Ycq�YNA��~�W�x��t�W��鬦��d;�M6Ac봂�������53xV\��24?ړ�>��S�{=Qcm�[��p9��m�=2Exg����6}+Ě��i�E��f��&�d?��z��|7�3x+÷�t�]j�!܋�ʮ���UP~�"�?�	�?�M�k���O��Em Em�r��t"��:���S�)�M�������=�-A��[]r;e���Ɯ�p��q�]���Z;�Ö�8�t�8?���v^���jWN�m,���,c� ��9S�ތ���MRԿEWI�W��]�և���R�P��Ed����_O�)⟆z�5�+U� k�ӕ�;Y6�v9u#��:g��a)ۗK3ZrQ�7c��[��C�����0���P����w��d��}:���?�n����c�g����9-��b;=�]'�?g/���-�ɘ�ɧ�c����4��O�:}��3�z�)ϓur��b)�j*F�H�4���7̙⟳�û�x��T�N� �IX|� �z���a�+��U�YXV ���N�e���X�ְ�ZD0�T~�V�y���.f|C��=�ßM4��7r��FU3ɉ���?y]h��^�N��p8�)�8eaЊ��t�M^�[;�h��ed�t�=��)�?e��XR� K�����H��¸Q~�;x�{�����gb�j���r���Cp`��˵�~�J��D~#��<Ik�x���/�oc�򎃞_ۧa����ً���-�ެ�r��N�TU�r+�m�⳷�"Ha�B�q�UU ���BWR��c�R��`�x�ayf�b�Y݂��ROJ�O�v�i�o�����qj�$i,g*�S �S�}���h�:�c��,��i�e��x��Q����G�������;s� �+:ԧQ�ZȺU#M>������/�4+���om{�7���ʿAlo�+d�������VX\:��#����g��/���� ����|=��n�����]3~�I7m�>�u=+Z1�8�:���2<#���Ks���~�n�A"��Gɍ����`Lg�|�_�졔�2�?~Ξ	�E��>M.��^M:O(�G�+�xg{�c�խ#�<�Qk��a��+(��hWdRI)�`��/�q�z�_��.�����2�t�7�D�B#�L�(=��wz��;��iu[�a��P������t���a���ah���\}O��<�ѧ9�[�Jp��֣�=�Ú=��c�ik�������׍~�_g�u�~"Ѡ3j��츷�|��9���} �t���MUVft�黣�]��F��^���������i�8�]$���h��c���E}9�/����ėW�g����ʑ���c�A5����~�ew��.�丌)��c��8ѫUi�{��'���Y�
�=�^M?P��<�1�t������ |��6�rNo5+�7׌I3L�NO$�?^���_h^	�6�&��m��.� �9�7�kr�a][�g,�͢VG���^��7Æ��6�m:�n�5-VV����5�~�זx��k�_�/����ZT��gJ��}v����W)M�=N�u�c�.���7���<�}�H��Ԭգ�0����7��H�88�w�]��>�⋖mN�W�x_��F�؀� �l���W�x�;�kI��9��P�Hc>�������:g�4y��^��v�/�ɒ#� ���4�-Y�'�Gȟ�׉,�G�'6Q��gi�I]���u��μ��o��C.3庾>�5���3�ÿ���v�� �Q� ���������� �Ȱ�������լw:�c�].�� N����d'GC��zW���e�������kvW�^�u*?�$}H�+���
�P��B���u$~S���dʒ0�GP+��ʐUci����ϯ����'T���X܃4da�����5�|>,Ѯ41���jt������{n���/����>��!�a��۱�S�+���@�d��7}�Z��h��+5�#V�y��p�.g���/'��!j���w��`���J(�{��{��� �����n�"֠0��1췶�|�񞥽�;���o���d��t�k��wt�Y�'�?�]�h�Z�X����-�e�K1
�d�p|U�Ik�> ��s%��wP�oXX2�$:�v?
����U��iz��ڬ. Ž�pʐG t5��8�;� �{� 'n��*��V�l�J�iݽ��V�r�z��Ì�S$�����-7Z��4��;�n,�M�x��3��\�3�ÿ���v�� �V���w�N��h��G���E�?i���NW�/���EU(N�\tIFrL�C����G�KU���K�dX�Fr�U l�9��3�7��~#ү&�+{���U�?���?�~� н� ��?�r��g��/���� �!B�$�u{��6�hJ<��=OԭukHb���e%��)QVk��o���+��?�	t�ҏ:I7�>�u=+��@�
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��<���� �l��� b��� �r���_�o�9� &��c��-[� H��ʸ1;�������	� ''�_����ʿ^6��G��NSĿ�)\� �e�~�V��Χ�&�6��]&Bm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����Pm�m- �h�KE &�6��@	����P��Ӌ� ��c��-[� H��ƿ�� �s�M�����Z����_υpbwGM-����NSĿ�)\� �e�~�W�?���SĿ�)\� �e�~�m��|u>!�S�Ѷ���N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m��N�F� m�m�2��?���,� إ�����W�!�N��c_ş��o�#������n��[3�?�$�����R�� ��*�z��_�$�����R�� ��*�z��|u> ��+s0��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��( ��(�� i��6��?�)j��G-=������m�Rտ�Z�{���:)lϸ� ��_�r�%� �F�� K,���m~B� ��|K� b����YW��kC�"��&�6��]Z�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�����P�e�O/�c_ş�uo�#������?i��6��?�(���G-=uň��gܟ�H�9_� أs� ��U��_���H�9o� أs� ��U��Z������Et�m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\��{�M�����:����_�]Bߴ���Oů�uo�#������n��[3�_�#� �����Q�� ��*�}��/�#� �����Q�� ��*�}�(�T���(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(��(�3��� ��~-أ�������-�O� ɴ�Z� �GV� �9k��:��z{s����K� b����YW��~@����K� b����YW��mG�"��QElfQE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE y��A� &��k��_� H�瞿��ڃ�M�����:����_�=q��������� '/�_�n��ʿ`+�� ��� '/�_�n��ʿ`kj?>!(���3�Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((����Z((���<���� �h��� b��� �R���_���A� &��o��_� H��璸����}�� |� ���/��7?�[e_�5��� {����7��7?�[e_�[kZ?	>!(��F���J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� J)vѶ��]�m���mh(��F� �ڃ�M�����:���K_�%C��� �4|[� �GW� �)k��:��z{t�� ��<M� b�����W�~>� �� ��<M� b�����W�kG�"��QEoc0��(�QE�(��,EQ`
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
(�� QEX�(��QE̿j�6��(�� �E-<U���?�l� � �Q�� �Z�x����ޞ���G��9�ءu� ��U�	��?��?�s&� �B�� Kl�����=��F�Z+b�F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih�0����1��������K_�C� ���?ſ�5�"���
孺6���O���Nc���(]�m�~�m�ǯ�#���ω��P�� ��*���)|$�q6Ѷ����M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� M�m���mih��F�Z(6Ѷ�� �ډ��.ء�� ����W�E�Qɳ|\� �CW� �)k�ݮJۣj{u�� ��|M� b�����W�5~=���M� b�����W�=kK�&{���Elf6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E 6�u �)�Ph�Q@��E y�E� &��s��_� H���v��Oڏ�M��������K_��rV�S������ '5�o�.���ʿa������ '5�o�.���ʿa�Z_	3�(����(�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� (�� �ڏ�M��������K_��D����/���5�"����䭺6���_��}����4-�������g�W�F�_�?�K�����$��� mF}.`?�-~�U�~�3�n�F�N���v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.v�6�u\�4m4�(��h�i�Qp�ѴӨ��7i�i�QE�n�F�N���ݦ���E�MM:�.�~�y?�Ŷ=?�՗�Z�{��]F�:���k���&�?��a�g�w�� }B����_��lxf�m���2��r��S��?�_�D�����3y�i��J�� �_�~��5���1� ?�RW�ț�Z� �� Q�m������-]/���7i�i�Q[�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�ݦ���E�����c�ǿ�3�~���� ��G��_�?�� ���7��~G��k�[n�z[�7��_c� �Q����ƞ,�-��z}����r���UG�ݠ���}s$J�<Oyr�˜��_�#\�w5��ǿ�2� ¯� ��i���ʹ�惪��m�_��'�_�5�9� ������o[�!�h
���������/ʿS<�H|g�����5k{����,j��UR�)�^�u�ga�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ea�S��,6�u��N����)�PE:��h�Q@Xm�(��Ec�O�,w�����Z2���<D�0��Eo6G��"¾�� �gxk��؛�L�f����C��6�gS� |�|?� ��h��O�/#a�,�5)W=|�4'��I���<��|"����Ò&���XZ�1��-�?�n?�q���������?���c��N�.��*�0'�UR+ۿ���0� ���?�kvd�-��i�3���W�A� �F��Em� �H>����~"[G�w�Z&�%>��,���,�񯘿��~?���C�_�.��]]/�V<�71c�=���Ӧ�!Kc��}�s�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���E�S��2�}\QO����)�QpE>�.(��E�e�(���M�E�݂"��c����~7��k� ��
C�x)� K�Þ�l��H>�ܾ~��o~%~����F�J��m*��+[�S����壷�^x1���jᓻ�ж3�I��x�ú��~�e��k-��z9�ɍ~.�0�믃?�_��k/��j��$�7�>�i!q���?گ�J�I��4�� f�(����bx����<U�8�o��>��ˑ�{�"��T~��P����Y�wk*�mq��r��2{EM]�EPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEP^9�b|F_�?��_��M�zL����q?�"#�"��ױ��_�S�w�O�_x�m�k:���ҩ�ɷM����0#�*�;+�nQ� �!|8d������d��C����\� ���6�k� �t�-� �O� x�X|��Z��{����S#���_JW!�W�������/
|D��uφ����˭���YQ�S_����|�[m�࿌�rޙ5�N�#�n�d� ����4�������-��O\M�j���e�ـ=�ahNO|��{E~H�"�(�x�g��kA��Z����c���з�_S����+�º��d���(�)h�$M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M����� M��?j�����
/c�-:f�K��m�1F��ᅋ^�1�X�� �տ�⥷���1����t}:I��C�����FO�R�/�
�9���_
��k�Ɵ�u{���K{�y2^�Ҹ=�cW���ƣ�i��ݍ�m����K��kQ ¢(T{ OE��QE��� �~��&��Z_���b��nb�=�ߕ�_������q_���V��	��#�&�N��Y�}k/����������
��?�Ưٚ�Wӭ��!�c��je�������fS
���_�H��C�� �?�ڕ�}[�r��A�-%��H�˔��F;V�޶"K����h�]BQK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R���i��w����Z����,q��f'� O�@�_�Y��b��%���^���DSȉ	���?����M}���� ��.xCC���>�������0��
��#��~^|ѯ�� �P��A�A$����P�	Xt�fU��=�r"F���5��\�ww7J�(����Q@�wP�Ѭ�ȥ92���=A�a�Cw� �� ��I"�4~fX�~Ѣ]���.Gl�o�_�u����i� ���-|u��y�&�H{���/>���Bq����3�$�婁�����coyi4w6����I�*����jj�c�	?�F����}�\��!�hX����m9�����c�Q�}ӴWR����h�mm�!�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�QpE;h�h��6�v�F�E�m�����)�EE�S��6�.h�mm\�N�(�(����Q�Qp_� �W� hd�g�F? �;<A�2������A����v��&>�}��j�6�u�_O���O<��Ԏ5��� 
�X��u/�)g��^U����d�-���UiK�%��+9�Ȩ�����$���~��m�MR�����[�2(����� ��ב"zW������g��1�Z�ƱEJ#EUP:  ���(�� (�� *+�Xo�����;�y��$2�du#H<G5-�[�c�:��C��ִ8f��SJo�#�E�3m�Ա�z|�3Фnz��������x_I��v��F�k��g!�u����;+�� ������]��#��@/�hnG�#���}%P �c'�_� �&� jyt�����-�x$W����tv�n	i�z�G��?�+X=lȒ�~��K�Ѵ��bQK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK�Ѵ�QK����և�/ᖽ�oM��U����Hx������4�W����� ��v�	t��߉#���X�3�\����o�	o�-���S��ݧ��$w�"��l�M���!��=\�|)��%����]j� |u��֙x���$�#�[� dPG������Y;�d���(�QE QE QE W���3�dԾ�W�><xK�Xj��=�*Ciڲ���p�)�뢶~��\��_���ok~�-���lַP�U�2��e8ea� Ҁ<?�;����ړ���Hv����=oOR3�Q��y����:��s��}��_�K_���PK�G$G)��մ�s�t<��<vtu��$���[�O��7����:���@�6��픖6?��:k�2�1�lj�Eb�QE��(�,QE`��(Q@X(���EP
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
��xV�W��#�e��1&�b�&GݒNU{���!s�� ����������u�r�%����X��#����x��
�r�\WS��_�w�� f��'�4�.{�W����oo\6^y�T��;W��EdhQE QE QE QE QE |��q~�zO�i�m9R+O�+%Ɓ��% n�C� <���![��?�?�N� ��V��>#�|���i%Ơ���@�mP-�G$�E#u=�n�ƿf����
u�����>���*���
�ڽ�
��x�q��m�)�gqq�_���L�ۨ������-O�*;T�R��k�qk#��)?}F>���G6�JwW3�]�m�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�QK��� �R�m %�h�@	E.�6�W�~�ߴ��+�'��>����7�h�I|5��g��>��8��s��׆��>�|e�����4���O!��ƿ��x�8 �����G� ��~��^�.��k4�^�D���G٦�ݛ,p��ĥa�s����e�{���٫�^���j~��77�\.�{�[e����p*�o�Tq�1�"�EUT` : +��s��A�O�m��];B�-����NHQ�f?��Ifc�,I�]%s�Q@Q@Q@Q@Q@Q@Q@���S��A�/~7�)����/۵�'N$�����/#��@�������O/���G���7����=�8�����£�t(uye�p_�Ddu�0U�A��� ���O��$�_a���`�څ����D�6����y �?�}ʌ�D���*+��7�
�~�]��|_-���.�,m�!��G2�:	02яv^2�������R�4m4�J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J)v�6� J��}�����Ox�R�I���L�7S�����@�_x��� �#�x��:�6���e���l*�����$� &��>||��� (�ե|>�}�]A�h�-��l�Wh8k��
 >�ڻ��h��F��>(�I��� C����	�N�kZJ�O�����@����Ǯ22I-��7�ξ��~�x7°e#����"�=�� <���z*�J���d� �ɿ�����}���X�0��L�;�L�	'$�}ư��(��@QE QE QE QE QE QE QE ɡ���TYb�J�8ʰ#�S������	Ũ|/Ԯ~.�����-e�~��饖}-���M��|�F⣘����� �G�>4A���M�g��D6Z��G�� =��z9��'m~�W�_��� �_�j?>�%���]j>������K��ѿ��R�Bj����c�� �N��Ѷ�s����>�� d��Mڟ6�&���Ǜ�c����&��B��Y����X��ä��*�Â �+t�3�-QTEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEP\g�ό��{� ��V�I�-�c�'���(���c���G�O~����|.�����ֳt�t�́��a��x27�N�Wh>��� U��5��� ٞӤ�%���� c�.5��fa��1�Q�ĥa�q��)|Z� ��|p����<-k)��K.E��9�ח���؀�K~�~�� �O��d�.������v����<`O} ��jI�8��$�� ���o�_����
�6��C�K�Jp��`0e����*��z�a��QE�(�� (�� (�� (�� (�� (�� (�� (�� (�� (�� ����� �q�O�����'���¿�o���}H��]*��`J��0 �>~�� ���>��k�7D�5/Z>�4+��h��Oe)ʴm�۝��
���ï-��?f� ��>o��G[�L���o��,�� 2���S�l��>|n�Wǿ��_k��Μ�YU٭�#>\џ�7��Q��]�~(�Y��~:� �6�|�8�6�w��K�ں�h^,�C�[��9RpU�p>���৞���?�����"�<���[��ec�1?�����֊W!��Z)h�$J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h���� J)h�������|�+�sx��������D~sf[�>\Q��G�P}N4�W�߶G������ ����8
�Kx�O�$�w�� [ ��������E�C>'~�~!?~h����jt�NR�����$d$C*��w9��_� ��>����0�k�O��ko���G�3q������� 3r�Z]Ϛ� e��O�G������{Wԭ|%}(�������n�b8q�|m�S�b>|9���	X�c�-���Y.�l��j��f=Y��I�]q�1�"�EUT` : )ՙAEPEPEPEPEPEPEPEPEPEPEPEPEPWV��Z�msw�!�HeP��F
�<A�~t~��H����&�8����2I�َ�6��8���v>���aM~�Q@��	?n���� �"��-P�tK</�F�J]�ő�Z�� �	t�
W�~�~Ͽ���_�cG^ף}E=·{�o�ދ'p� i/�z��>���Y<=�Y��N;�hLKn�c|2<m�Ґ{t�˟��	!�υ�����f�u���9��Mk��j��d��e*��w���RhV?Y��?��UO��5������6O�\O<c֭1�Dp�)8F=K��'�!�S�0��,Vox��Q�	�]*c�_C��`s� WК�I2OU��EP�6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\m�(���Eq�S��.6�u��N����)�PE:��h�Q@\mꫪj�:�q�jW��}��y�]]J�E���� =�r�VԵ;=O���.౲�C$�72,qƣ�31 �k��#�
����w:g�"?|B�A=�t�[զ#2��#N1�WŶz?�O� 0��YZ�|&� [.�?@���932�O6A�J�$�>������7�������[\����e�`lu�pӐ�N�3+毄߱����� �x��N��xj���$D� ���}	S�<��e��%��?���[�:G��1��֥ V���͹$G}ǀ@^���d�eE�;��� �f� fx#DK{��%�u�o�ps��p8�!�W��E 
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��
(��<���������E��<-i�Kl��Pom���tÁ�v��~l�t� �;���w����(�^Kg���	��S����\.��n��+�ך(�;��Q���\�������lv���a��RT�ۗ��g<Wܿ� ��>.$6ڎ���X|�~%�9��?�d'Ҿ���ß
�P��Ѽ]��3ĺ[�6���N�� yw��`�ƾ��� i�q�&�� �ι}�=A�aas��z�a*}w�aj���}�e}o�ZEuiqմ�9�p��zÂ>�5~,j�_���M�{�Z� H�����k�y1�d�a��~h��]��o�,g�O�e|K�E��d�o.i�i��G]�U���
�|�V?]诐>� �T�|D�a�u��߿G�-
ǻ�&��`{�_¾����/�A��g��I������Pr�j�#~�Z)��R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E- �R�@	E-6GX�gv�2�� S@Ex�ď�+��x�� ���F��y�ecq�ې}Poq��_$|P� ��x3H����W�S�o5���� x"�����4��~�מ�W��>���O�|e�xt�2%��ᮥQ���G��T��1{�Y~׿�5��w��k.F(��6ɭ`�<bK�%���Uһ��F���6�]_ⷍ-<<'2��Ŏ�!�RDj����Q��,v�� ೚e�����BMJnUu�� ����=:9Z�='���g� � ԼE>�o�g�̊�^'O�aRG����8�r#����~�����	|�=·�5}n�x��a��.�3�/Z�"���|3�=�$��)6���_�^ �WP���cl�'������_o�Y[�p�ZA��#�P"F�`*�� ;
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
�$|�{����_�%ڲ�1�4�I��}Ԋ�h���)�~����𞡮x�V+y��h3�ٓ�W��0� �7|b�5�����m�k��L�ͥ��2��+�R� �I]S����s�����߂e�O��LO�~�1�Z���g��>��xT�/���Y�au��ʏ��_�U��
�-��&���k6���P�I���pE;�X�������w�>j�y�6�uX�8�,ש�w�
��ZU�o�$�=~ߥ��� ~d���Ŀ���o�w�_D_��셉��)^U�/�$���3O�{i�����2S�adli?�SٻW��DV���J���C�ֺ�/ۋ��x�,�aA� �ע#�>|��c�	�w��y���A�]�Kyr��2�����_�A���/�=J{{K�^D���h��#��X������q�>/x$�Mz�?���� Y�O����� ���x��6:>� &[��� =O�Q\SpH���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh���������� � �Oe� �h� ���'� E��?�S����oh������_����+����2?��%�)+>�;�M�w�\�y� �z�/����â���7$fQ���^�ῃ�6��ι�\� �944s0��a��A�g�w�>)��c�}{��W��Ug--\���� V+�6z-�'��F��"�#� � �y�$������3�%պ� �P5}u�� ��X��7�1�r��!P}�w\�,r�$� ��|(�ܺ'�<[���k��co��]�5�#���]�VO|.���f�5I.��)q���ڞ� �g~����|7��&^�jZ���謁�?%�\������.�xs�ǅtYפ�z=�r�.q�M�,����o�����&� �K��:`� �IR��p������_���7��J��c��1�yH=��i��+��q��Rһ�e����l��>��,�Ǚc������ �U�����'O���L�> ��o���w������$�j}�})E  ����m"������%�0 DA�pOE QE QE QE QE QE QE QE QE QE QE QE��awk '{i=1;while(i <= NF){col[i]=col[i] $i " ";i=i+1}} END {i=1;while(i<=NF){print col[i];i=i+1}}' inter_result/snpmatrix.txt |sed 's/IID //g' > inter_result/snpmatrix2a.txt
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
