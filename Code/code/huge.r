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




load("C:\\Users\\MSI\\OneDrive\\æ¡Œé¢\\èŠ±æœŸå¯å¡‘æ€§\\data230605.RData")


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
write.table(ph,file = "C:\\Users\\MSI\\OneDrive\\æ¡Œé¢\\GWA_shiny_app\\data\\input\\SMR_type1\\phe_FT10.txt",col.names = T,row.names = F,quote = F)


ph=data.frame(family=1:nrow(ft1016),id=ft1016[,1],AT5G55080=p10[,"AT5G55080"])
write.table(ph,file = "C:\\Users\\MSI\\OneDrive\\æ¡Œé¢\\GWA_shiny_app\\data\\input\\SMR_type1\\phe_AT5G55080.txt",col.names = T,row.names = F,quote = F)

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
write.table(ph,file = "C:\\Users\\MSI\\OneDrive\\æ¡Œé¢\\GWA_shiny_app\\data\\input\\SMR_type1\\phe_AT2G33600.txt",col.names = T,row.names = F,quote = F)



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
  relmat[upper.tri(relmat)] <- t(relmat)[upper.tri(relmat)]#ä¸Šä¸‰è§’å˜æˆä¸‹ä¸‰è§’
  svd <- svd(relmat)
  Z <- svd$u %*% diag(sqrt(svd$d)) #å·¦è¾¹çš„çŸ©é˜µä¹˜å¥‡å¼‚å€¼
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



‰PNG

   IHDR  ,  ,   ^uè   sRGB ®Îé   gAMA  ±üa   	pHYs  2À  2À(dZÛ  ÿ¥IDATx^ìİ	”­eu áá2Š ˆ" "“È £8
ˆ‚™DTAÃ$ˆŞ âePQæI¦ ¢ˆ‚ bœ5Q4&&Ñt¢±Mì¬îô”ît:æï½«KRç|›Ÿ[poÕ©:Ï^ëY«×îúNóŸ¡¬•÷Œ1ófV„Â^á°p|X.·„;Â=áğ‡áÛááOÃÏÂ/Âß‡ÿş5t         ![£l²=Ê)[¤l’²MÊF)[¥l–²]Ê†)[¦lš²mÊÆé¸ÍS¶OÛ‡l¡²‰2ÆcÌÌÓÃ¶á€pR¸$Ü
~ş%Tÿ#         `.Ë6*©l¥²™Êv*ªl©²©Ê¶*+cŒ1ÆôLş+ ­Â¾!ÿ•Ğásááğ_CõC        €—­U6WÙ^eƒ•-V6YÙfùKÍÆcÆfV;„·„üø`øu¨~x        °ôd«•ÍV¶[ÙpeË•M—1Æ3gçùa¿pZ¸-ü$T?        ˜=Ùveã•­W6_Ù~cŒ1#7†ÃÃáá¿‡êÛ2µÒÊ«tk®µN·ÁÆ[t[mÿòn—İöëö|ı‘İAo>±{ëq‹»ãN»´;yñ5İû>tCwÚÇ>ÕyÑgºs®¸«ûÈÕ÷u\ÿåî’›¾Ş]yÛw»kîøAwıİÖİxßO»[¿üóîÓı         W¶FÙe{”R¶HÙ$e›”R¶JÙ,e»”S¶LÙ4eÛtà›Oœh²yÚr»—M4PÙBeUµR3 °lÁ²	Ë6,1cŒ1fFçÅáøÿªæoCõk©Yq¥•»uÖß¸ÛfÇ]»=ö;¢;üèEİ	g\Ñ~ş-İy×~±»òÓßënºÿgåÿ         ˜Ë²ÊFê¼kîŸh¦²Ê†*[ªlª²­ÊÆªj¯–²lÅ²Ëv,2cŒ1f©Íêá5á¬ğ¥ğ¿BõÃèq[~áÂîÙën0ñW‘wßçĞî£Nø—Cg_vgwÕg¾_ş        àßek•ÍU¶WÙ`e‹•MV¶YÙhUíÖ”-Y6eÙ–ec–­™1Æ³D³bØ'\ÕšÇm­µ×ë¶ßyî€Ã›ø—>~â+åO         –lµ²ÙÊv+®l¹ªÆë	Êæ,Û³lĞ²E3Æc™†ÂÂoBõƒdZVYõ©İ¶|ñÄvà¨Ïí_zGwÃ=?)        0ó²éÊ¶+¯l½²ùÊö«jÂ‡lÑ²IË6-5cŒ1c6Oû…ËÃOCõÃb‰-·Üò?¨ö=ø˜î”³¯í.ûÔ7Ën         Œ¾lÀ²Ë&,Û°lÄªvlš²UËf-ÛµlØŒ1ÆÌÃÙ*¼'|1T?–Øj«¯Ùí°Ëİaïx·ø’Ïv·~ùçå-         æ¾lÄ²Ëf,Û±lÈª¶lš²eË¦-Û6cŒ1sx^>~ªü%òœçnÔíº÷Áİ;O>¯»ğ_)         0>²%Ë¦,Û²lÌªöl²qËÖ-›7cŒ1s`¶g„‡CõÁş˜şÌµ»]_óÆî„Ó/ï®¹ãáò‡         üV¶fÙœe{–ZÕ¦-¡lßNÙÂcŒ¡yAXş8Tàió­wê}ûûº_uoùÃ         –T¶hÙ¤e›V5kKèB¶qÙÈcŒ™…Ù8¼7|'TÔ½Ö|Öºİî¯=¤;éÌw7Üó“ò         <QÙ¨e«–ÍZ¶kUÓ¶²•Ëf.Û9cŒ1Ëp–G†¯„ê¹×ÛìÒ~ô¢î¼k¿XşP         €e-¶lÙ²i«Z·%]¶tÙÔcŒYJ³C¸$ü·P}ø>ªÍ¶Ú±;ò]gvWŞöİòƒ         fK¶mG¾ûƒ­[ÕÀ=†lê²­ËÆÎcÌã˜UÂÑá›¡ú }T›n¾}wÄ1§w—}ê›å<         ŒšËoùÖDû–\ÕÆ=†lí²¹ËöÎcÌcÌKÃUáŸCõ¡ZÚè[w‡½ãıİ%7}½ü         €¹"[¸lâ6Şl›²™ë‘í]6xÙâcŒ™2«…ãÂ‡ê´´âJ+w‡ujwá•Ø         0×]ôÉ¯N´rÙÌU-]lò²ÍËFÏcÆv6
ç‡ÿªËÆÂVìvßçĞnñ%Ÿ-?˜        `¾Êv.ºléªÆîQd£—­^6{Æ36³s¸)TŒ¥M6ß®{û‰ênºÿgå‡0         Œ‹lé²©Ë¶®jîzÜ²á3Æ˜y;û‡Bõ!Øxò*OéöÚÿÈîÜ+ï)?p        `Üec—­]6wU‹÷(²åË¦ÏcæÍ~ª½Ææ[ïÔ{êùİ­_şyùá
         Êæ.Û»lğª6ïQdÛ—Ÿ1ÆÌÉY3œş.Tr–_¸°Ûûoë>zíå)         °d²ÅË&/Û¼ªÙ+dë—Í_¶Æ3ò³V8/ü&TjşÌµ»C:µ»şî?+?4        €Ç'Û¼lô²Õ«¾B¶Ù fhŒ1#7«‡sÂ¿„êClÀoÙ½óäóÊH         `éÊf/Û½ªé+d˜M`¶Æ3ë³jø@øŸ¡úĞ°Í»v§{}ùa         ,[ÙğeËW5~…l³ÌVĞcf|V‹Â?†êCjÀ+ö<¨;÷Š»Ë?         `feÓ—m_Õü²Ìf0ÛAcŒYæó;áäğëP}(=b¹å—ïö=ø˜îò[¿]~Ø         ³+¿lı²ù«ZÀ!ÙfC˜-¡1Æ,“9>ü2TB^÷»Gw×ÜñƒòÃ         -ÙüeûW5…l	³)4Æ˜¥6‡„¿Õ‡Î€×ğ–îÊÛ¾[~˜         £-Àl«F°ma6†Æó¸ç%á¾P}Èxõ¾oê.½ù›å‡         0·d˜m`Õ²5ÌæĞc–x®Õ‡Ê€]÷>¸»ğ†‡Ê+         `nËF0[Áª!,d{˜¢1ÆôÎ)áŸBõAòˆ—ïñ†îc×=X~8         óK6ƒÙVMál³E4Æ˜f?Õ‡Ç#vxéİ‡¯º·ü0         æ·l³%¬Ã!Ù$f›hŒ1¶w…êÃâëoøÂîä³®)?|         €ñ’Ma¶…Us8$Ål1c8+„Cõáğˆ'¯ò”îˆcÏ(?l         €ñ–a¶†Uƒ8$›Ål1c2„¿
ÕÂ#öÚÿÈîºÏÿIù         ²5Ìæ°j‡d»˜£1fÏZá¡úxÄ6;îÖ}èã÷–*          •l·}Éîe›8äúM£1fÍÛÂÕÂ:ëoÜø+Ë         €%‘-b6‰U«8E6Ù6cæÁlîÕ›}Â“´\wØ;Ş_~h          <Ù&f£Xµ‹Sdã˜­£1fÎ{Â¿…ê>a‡]öì.¾ñkå         À‘b¶ŠUÃ8E¶Ù<cæĞì¾ª7õ„§<uõî§|´üp          Xš²YÌv±j§Èö1HcÌˆÏBõ&~Ä+÷úİîÚÏı°ü@          X²]|Å•mãl!1#8›†BõÆğ¬ç¬ß¼øšòC          `&¼çƒWO4Uë8Åƒ!ÛHcÌˆÌÛÂ?‡ê;áµÕİò¥¿)ßø          3)›Æl«æqŠl#³‘4ÆÌâ¬®Õ›tÂ†›nÕqÁmå›         `6eã˜­cÕ@N‘­d6“Æ˜=ÃÏBõÆœpĞ›O,ßÜ          £$›Çª…œ"›Él'134g…êÍ8á9ëmØ~ş-å         `eû˜dÕFNqv0Æ,ÃÙ<<ª7à„İö>¸»ùÿP¾‘         FY6ÙBVä_ÙTc–ò¼3ü&To¼n¥•Wé=õ‚òÍ         0—d™mdÕLNÊ¦2ÛJcÌRšËBõf›ğ¢^Ñ]zó7Ê7,         À\”md6’U;9E6–Æ˜'0Ï…ê6áà·R¾I         æƒl%«†rŠl-³¹4ÆLs^ş!To¬n½çmÒyÑíå         `>Éf2ÛÉª©œ”Íå>Á³„ó¾P½™&¼úu‡w·~ùçå         `>Êv2Êª­œâÔ`Œé™Â'Cõšğ¶ãÏ.ß„          ã [Êª±œ"[Ìl21C³uø~¨Ş8İZk¯7ñçÌ«7         À8É¦2ÛÊª¹œ”Mf¶™Æ˜É9,üŸP½aºívzUwİ?*ßp          ã(ÛÊl,«örR¶™Ùh3ösn¨Ş$^è»Ê7          ¿ìö;äØ²Áœ"[McÆv>ª7F·üÂ…İ»]\¾±          øwÙ\f{Y5™“²Ù4f¬æáP½!ºçm´y÷‘«ï+ßP          ´²½Ì³j3'e»™§1ó~¶?Õ¡{Ù«öï>õÀ_—o$          ]6˜ÙbVæ¤?Ùr3ogğ_Bõè~Û{Ë7          K.›ÌªÕœ”-çÁ˜y7G†êE?á˜÷~¬|Ã          0}ÙfVÍæÙv3ofQ¨^èİ*OY­[tŞMå         €Ç/Íl5«†sR6ÆÌù¹,T/ğnİõ7î>vİƒå         €'.[Íl6«–sR¶ÆÌÉY1Üªv·åv/ë®ûüŸ”o          –l6³İ¬šÎIÙ|fûiÌœ™§…‡Bõ‚î^¾ÇÊ7          ËN6œUÛ9)ÛÏl@ùyVøN¨^Èİ~‡[¾          Xö²å¬ÏIÙ€fjÌÈÎóÂCõîŞò{‹Ë>          3'›Îªõœ”-h6¡ÆŒÜlş2T/Üî¤3?^¾à         ˜y'}ğª²ùœ”Mh¶¡ÆŒÌl~šìÂ…+t‹>rcùB         `ödã™­gÕ€†lC·	ÆÌúìş!4/ÔUV}jwæE·—/p          f_¶Ù|V-hÈF4[Qcfmvÿ#4/Ğ§­±fwÎw•/l          FÇ¹WÜİ­şôg6=è¤lE³5fÆgŸğ¯¡ya>óÙëu½öò         ÀèÉö3ĞªÙŒ¾63cs`¨^Œİ:ëoÜ]|ã×Ê2          £+ĞlA«FtR6¤Æ,óÉ:¾zvl¼Ewå§¿W¾€         }Ù‚fZµ¢“ö	Æ,³Ù=äŸôn^|›n¾}wİ?*_¸          ÌÙ„fZ5£![ÒlJYê³sø¡yám¹íK»ïûËò         ÀÜ“mh6¢U;²)Í¶Ô˜¥6Û†Í._ˆ·~éoÊ*          sW6¢=Ñr¶¥Ù˜ó„g³ğ‹Ğ¼ĞòO}ûËÊ          óW¶¢ÙŒV-iÈÆ4[Sc÷</üeh^`l¼Ewİ?*_˜          ÌÙŒf;Z5¥![Ó‚1Óg‡†æ…µÎúwW~ú{å         €ù'ÛÑlH«¶4dsší©1K<Oß	Íê™Ï^¯»øÆ¯•/D          æ¯lH³%­Óíi6¨Æ<æ¬
Íiõ§?³ûèµ”/@          æ¿lI³)­ZÓj¶¨ÆôÎ¡y­²êS»s¯¸»|á         0>Î¹â®‰¶´jNC¶¨Æ<ê\šÎÂ…+tg^t{ù‚         `üd[ºp…›îtR6©Æ4³(T/˜nÑGn,_h          Œ¯lL«ötR¶©Æ<2G†ê…ÒôÁ«Ê          dkZ5¨“²Q5fÁ¡ztoù½Åå          ~+›ÓªE”­ªãÙ2ü—Ğ¼8ö;äØò          Ã²=­šÔ­j6«fçáÇ¡ya¼|7”/$          x4¯ØóÀ¦K”Íj¶«fÌæĞ¼ ¶Üîeå          K¶¨U£²]5c4ŸÍaİõ7î®»óGå‹          K¶¨Ù¤V­jÈ†ÕŒÁœšÀ*OY­ûØu–/          XRÙ¤f›Z5«![V3ç°P=ñİ¢ón*_0          0]Ù¦VÍê¤lZÍ<œ­Ãÿ	Í“~Ì{?V¾P          àñÊFµjWC6­Ù¶šy4+„ï‡æ	?ø­§”/          x¢²U­Ömk6®fÌ'CóD¿ìUû—/          XZ²Y­ZÖ«™sjhàçm´y÷©şº|Q          ÀÒ’Íj¶«UÓŞÌ}BóÄ.¿pa÷‘«ï+_          °´e»škÕ¶†×3çùáBó¤¾{ÑÅå          –•lX«¶5dóší«™cóPhĞı9¶|          À²öúCßÕô­“²}5sh.Í¹İN¯*Ÿx          ˜)Ûíüê¦s”¬™óÎĞ<k­½^wİ?*Ÿt          ˜)Ù´®µös›ŞuR¶°f„góğ›Ğ<yg^t{ù„         ÀLË¶µj^C¶°ÙÄš¯†æ‰{Ûñg—O4          Ì–l\«ö5<ÌÎÙ¡yÂ^ıºÃË'          fÛ«÷}SÓ¿N:+˜š=CóD­÷¼Mº[¿üóòÉ         €Ù–­ëzlÚt°“²‘5#0«†Ÿ…æI:ó¢ÛË'          FE6¯U²‘ÍVÖÌò\š'èà·R>¡          0j²}­šØ­¬™Åy[h˜íğŠò‰         €Qµõ‹_Ùt±“²™5³0›†OÈJ+¯Ò]zó7Ê'          FU6°ÙÂ÷±!›ÙlgÍÏƒ¡yB=õ‚ò	         €Q—-lÕÈ†‚™Áù@hˆİö>¸|â          `®È&¶jeC6´ffÇĞ<ÏYoÃîæşCù¤         À\‘Ml¶±U3²¥5Ëx¾š‹úù·”O          Ì5ÙÆVÍlÈ–Ö,ÃyOh.üAo>±|¢          `®ÊF¶jgC6µfÌ&áßÂÀßpÓ­Ê'          æºle‡ûÙMm¶µf)Ï¡¹àg\p[ùä          À\—­lÕĞ†lkÍRœ·…æB¿öÀ£Ê'          æ‹lf«–6dck–Â¬şs¸ÀÏzÎúİ-_ú›òI         €ù"›Ùlg‡{Úm¶¶æ	Îõ¡¹À'/¾¦|B          `¾Év¶jjÃ'‚ys@h.ì+÷úİò‰          €ù*Úª­ÙÜšÇ1+„¿
ô)«­Ñ]û¹–O          ÌWÙĞfK;Ü×†ln³½5ÓœCsAßyÊGË'           æ»li«Æ6d{k¦1Û…æBî°Ëå…         €q‘MmÕÚ†lpÍÎ]aà>éIËußøµò¢         À¸È¦6ÛÚáŞ6dƒk–`Í<ìï//8          Œ›lk«æ6d‹kc~.Üºëo\^h          WÙØw·![\Ó3§„æÂø+Ë‹          ã*Ûª½Ùäšbş)\°m_²{y         `Üek;Üß†lr³Í5CsEh.Ø‡>~oyq         `Üek[5¸!Û\3e^šµ×şG–          øÿ²¹­ZÜ®™œûÂÀzò*Oé®ûüŸ”          øÿ²¹Íöv¸ÇÙèš˜CBs8öŒò‚          ƒ²½­šÜ­îØÏ_„³ş†/,/$          PËw¸ËÙêõšsòY×”          ¨eƒ[µ¹!›İ±œß	¿d‡—îY^@           _¶¸Ã}nÈf7Ûİ±›“CsA>|Õ½åÅ          úe‹[5ºá”0V³bøu¸/{õå…          –L6¹ÃnøO!Ş±™E¡¹»îÁò¢          K&›ÜªÕÙğÅ¬ş1\€]÷>¸¼`          Àôd›;Üë†lx³å÷óĞ\€ox¨¼X          Àôd›[5»áÌ0¯gõğ?ÃÀõ¾o*/          ğø¼úu‡4»“ş)dÓ;oçœĞ<ğKoşfy‘          €Ç'İªİÙôÎËY+üKxÀ¯9à-å          ˜lu‡ûİğC¶½ónÎÍ¾ò¶ï–          xb²Õ­ŞğÑ0¯fÍğ›0ğ@_÷»G—          X:²ÙîxÃ¿…l|çÍœärË/ß]sÇÊ‹          ,Ùìf»;Üó†l|çÍü]x€û|LyA          €¥+Ûİá7ü}˜sthàå·~»¼          ÀÒ•ínÕô†l}çüü0<°WìyPy!          €e#Şá®7dë;§çõ¡y`ç^qwy          €e#Şªíû‡9;„´Í»–           X¶²åî{C6¿srvÍ:õÜëË          ,[ÙòVoÈöwÎÍaàl°ñ–å          fÆó7Ùr ñtS˜S³QhÈ;O>¯|Ğ          ÀÌÈ¦·j}C6ÀsfÎàéÏ\»|À          ÀÌÊ¶w¸÷Ù Ï‰Y-üï0ğ 9êÔòÁ          3+ÛŞáŞ7dœ-ğÈÏqaàÎ/¿pawıİV>X          `feÛ»pá
Íï¤lG~ş8Üñ½ßğ¶ò          ³#ßáî7d<ÒóÒĞÜñ^û@ù          €Ù‘oÕş†l‚Gv®
wxó­w*           0»6ßfçöwR6Á#9«„wøØSÏ/          0»²õîC6ÁÙÜîì“WyJwë—^>8          `veë»ÊªOh€'e<róÍ0pG÷ÚÿÈò          £!›ßá8d<R³Chîè¹WŞS>(          `4dó[µÀ!á‘™KÂÀÜdóíÊ          Œ–M7ß~ ”ğHÌòá¿…;øö?T>          `´dû;Ü‡l„³õ92Ü¹…+¬ØİtÿÏÊ          Œ–lWXq¥&xR¶Â³>_	wl÷}-          0š²î‚C¶Â³:‡æ-¾ä³åƒ           FÓâKïhºàIÙÏÚ¼7Ü¡WZ¹|           Àh[q¥'´Á“²µùN¸C‡ujyç         €Ñ–-ğp²•yAhîĞEŸüjyç         €Ñ–-pÕ‡l‡g|…;²ñfÛ”w          ˜²	î„C¶Ã3>îÈaïxy§         €¹áğ£4Â“ş8Ìèlš;rÉM_/ï4          07d\µÂ!â›ÓÃÀØtóíË;          Ì-›n±Ã@+<)â›‡ÃÀ8â˜ÓË;          Ì-Ù÷Â!â™…æ\~Ë·Ê;          Ì-ÙWÍpÈ–x™ÏÃÀ7Şl«Ë;
          ÌMÙwÃ![âe>?ßøÈw°¼“          ÀÜ”ğp7²%^¦³Uh¾ñ•·}·¼“          ÀÜ”pÕ‡lŠ—Ù¼'|Ã-¶Ù¥¼ƒ          ÀÜ–­ğp?²)^fóÅ0ğ?zQyç          €¹-[áá~8dS¼LæÉ¡ù†ç]ûÅòÎ          s[¶ÂUC²-^ê³_øFk>kİò          óC6ÃÃqÈ¶x©Ïåaàí¾Ï¡å          æ‡İ_{È@C<)Ûâ¥>?ßè¤3?^Ş)          `~Èfx¸#Ù/Õyah¾Ñ÷ü¤¼S          ÀüÍpÕ‡lŒ—Úœ¾ÁæÛì\Ş!          `~Ù|ëZâIÙ/µùBø‡¾ı}å          æ—l‡‡{âñR™ÃoÂÀ7øğU_(ï          0¿|øª{ZâIÙgkü„gŸ0pãOæÚå          æ§lˆ‡»â­ñ‹ÃÀïº÷Áå           æ§lˆ‡»â­ñ‡ÃÀŸpúåå           æ§lˆ‡»â­ñšÕCsÃ×Üñpy'          €ù)âª-Ù?îyM¸ÁuÖß¸¼          Àü–-ñp_²9~ÜsV¸Á]÷>¸üæ          Àü¶ÛŞ´Å“²9~Üó¥0pƒï<ù¼ò›          óÛ;Oùè@[<)›ãÇ=ÿ+Üà…ŸøJùÍ         €ù-[âá¾8dsü¸æÅaàÆV[}Íò          ã!›âáÎ8d{<í9>ÜĞ»ìY~S          `<dS<Ü‡l§=·…:ìï/¿)          0?zÑ@c<)ÛãiÏß†Z|ÉgËo
          Œ‡lŠ‡;ãíñ´fÃ0p#Ë-·|wë—^~S          `<dSœmñpo²A^â9<ÜÀ¶|qù         €ñ’mñpo²A^â¹"ÜÀ¾S~3          `¼ìwÈ1­ñ¤l—x¾nà”³¯-¿          0^²-îC6ÈK<ÿ=ÜÀeŸúfùÍ          €ñ’mñpo²A^¢y~8¼ÊªO-¿          0²1îC¶È9û…ƒ/ØòÅå7          ÆS6ÆÃİqÈù1ç´0ppı(¿	          0²1îC¶È9·…ƒGxnùM          €ñ”ñpw²E~ÌùI8¸øÒ;Êo          Œ§lŒ‡»ã-rï¬šƒ7Üó“ò›           ã)ãª=Ù$?êì¬µözå7           Æ[¶ÆÃıqÈ&ùQç-aàÀö;ïQŞ8          0Ş²5îC6É:„~\yã          ÀxËÖx¸?Ù$?ê<œpÆå          ã-[ãáş8d“ü¨óë0pàÂO|¥¼q          `¼ek<Ü‡l’ËY-|ñò–7          ²9îC¶ÉÍl¾ğÙënPŞ(          @Êæx¸CÙ&7³oøÂ­¶yy£           )›ãá9d›ÜÌqaàwßçĞòF          R6ÇÃrÈ6¹™ÂÀrÔ©å          ¤l‡;ämr3Ÿ_xÜi—–7
          ²9îÃ¡™‡ÃÀ}Ùå          ¤l‡;ämr3ÿ5|áUŸù~y£           )›ãá9d›<0O_´âJ+—7          0U¶ÇÃ=rÈFù‘Ù6|Á:ëo\Ş          ÀTÙ÷È!åGæ€0ğÛì¸kyc           Se{<Ü#‡l”™“ÂÀì±ßå          L•íñp²Q~d.	_pøÑ‹Ê          ˜*Ûãá9d£üÈÜ¾à„3®(o          `ªl‡{äò#ó•0ğ§ŸKyc           Se{<Ü#‡‡Â#ó'aàÎ»æşòÆ           ¦:ïÚ/´È“²Q~d~¾àÊO¯¼1          €©²=î‘C6ÊÌ¿„/¸éşŸ•7          0U¶ÇÃ=rÈFybVÿŸ+­¼JyC           •l‡»ä­ò‚¦,&¬¹Ö:å           T²Aî’C¶Ê¶Ÿ²˜°ÁÆ[”7          PÉy¸K;„{MYLØr»—•7          PÙjû—4É“²U^pØ”Å„]vÛ¯¼          €J6ÈÃ]rÈVyÁqSö|ı‘å           T²Aî’ÃñaÁâ)‹	¾ùÄòF           *½ùÄ&yR¶Ê.Ÿ²˜ğÖã—7          PÉy¸KÙ*/¸eÊbÂq§]ZŞ          @%äá.9d«¼à)‹	'/¾¦¼          €J6ÈÃ]rÈVyÁ=SŞ÷¡Ê          ¨dƒ<Ü%‡l•<0e1á´}ª¼          €J6ÈÃ]rÈVyÁNYL8ó¢Ï”7          PÉy¸KÙ*/øö”Å„s®¸«¼          €J6ÈÃ]rÈVyÁ¦,&|äêûÊ          ¨dƒ<Ü%‡l•üé”Å„®ÿry#           •l‡»ä­ò‚ŸMYL¸ä¦¯—7          PÉy¸KÙ*/øÅ”Å„+oûny#           •l‡»ä­ò‚¿Ÿ²˜pÍ?(o           ’òp—²U^ğS®¿ûÏÊ          ¨dƒ<Ü%‡l•üÓ”Å„ïûiy#           •l‡»ä­ò‚²˜pë—^Ş          @%äá.9d«Ü,Ë           èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]Ş_          0U›Úeu           OÕ&‡vY          èSµÉ¡]ŞşÕ8           0U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Ÿùê          ˜–ªMí²:          Ğ§j“C»¬          ô©ÚäĞ.«Ã           }ª69´Ëê0          @ŸªMí²:          Ğ§j“C»¬          ô©ÚäĞ.?ó‡q           `ª69´Ëê0          @ŸªMí²:          Ğ§j“C»¬          ô©ÚäĞ.«Ã           }ª69´Ëê0          @ŸªMí²:          Ğ§j“C»¬          ô©ÚäĞ.«Ã           }ª69´Ëê0          @ŸªMíò³ñÅ           ÓQµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]~ökq           `ª69´Ëê0          @ŸªMí²:          Ğ§j“C»¼#¾          `:ª69´Ëê0          @ŸªMí²:          Ğ§j“C»¬          ô©ÚäĞ.«Ã           }ª69´Ëê0          @ŸªMí²:          Ğ§j“C»¬          ô©ÚäĞ.«Ã           }ª69´Ëê0          @ŸªMí²:          Ğ§j“C»¬          ô©ÚäĞ.«Ã           }ª69´Ëê0          @ŸªMí²:          Ğ§j“C»¬          ô©ÚäĞ.«Ã           }ª69´Ëê0          @ŸªMí²:          Ğ§j“C»ü\|1          ÀtTmrh—Ÿûú¯           ¦¥j“C»¬          ô©ÚäĞ.«Ã           }ª69´Ëê0          @ŸªMí²:          Ğ§j“C»¬          ô©ÚäĞ.«Ã           }ª69´Ëê0          @ŸªMí²:          Ğ§j“C»¬          ô©ÚäĞ.«Ã           }ª69´Ëê0          @ŸªMí²:          Ğ§j“C»¬          ô©ÚäĞ.«Ã           }ª69´Ëê0          @ŸªMíòÎøb          €é¨ÚäĞ.«Ã           }ª69´Ëê0          @ŸªMí²:          Ğ§j“C»¬          ô©ÚäĞ.«Ã           }ª69´Ëê0          @ŸªMí²:          Ğ§j“C»¬          ô©ÚäĞ.ïüF           ˜†ªMí²:          Ğ§j“C»¬          ô©ÚäĞ.«Ã           }ª69´Ëê0          @ŸªMí²:          Ğ§j“C»¬          ô©ÚäĞ.«Ã           }ª69´ËÏÇ          LGÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Ÿÿf           ˜†ªMí²:          Ğ§j“C»¬          ô©ÚäĞ.ïŠ/          ˜ªMí²:          Ğ§j“C»¬          ô©ÚäĞ.«Ã           }ª69´Ëê0          @ŸªMíò®oş          À´Tmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]Ş_          0U›Úeu           OÕ&‡vY          èSµÉ¡]Şı­8           0U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vyO|1          ÀtTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—÷|;           LCÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vyo|1          ÀtTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vyï·ÿ          `Zª69´Ëê0          @ŸªMí²:          Ğ§j“C»¬          ô©ÚäĞ.«Ã           }ª69´Ëê0          @ŸªMí²:          Ğ§j“C»¼÷;q           `ª69´Ëê0          @ŸªMí²:          Ğ§j“C»üB|1          ÀtTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—÷Å          LGÕ&‡vY          èSµÉ¡]V‡          úTmrh—÷}7           LCÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›ÚåıñÅ           ÓQµÉ¡]V‡          úTmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa          €>U›Úåıßû5          À´Tmrh—Õa          €>U›Úeu           OÕ&‡vY          èSµÉ¡]V‡          úTmrh—Õa   H»úóÎ¿æ®	\{÷„¯»gÒ½.úƒ/ÌiÕcŸo>vÕóÒå7}©|¼0ú?³?¯/¾ş¾	—|âşî’ÆW~†\ù©‡º«oûzwİg¿Õ]ÿ¹ïvŸ¼ëûİM÷>Üİrÿ»O?øçİg¿òÓîÎ?üëîîoüm÷…ïüª¼ö 0Tÿ{z.©   0ûª69´Ëê0   |ğü›ß!ç«×¼şğòÌ'ÏYwƒò±Ïu{íwhùxa=uµÕË÷	K×“´\·ÚÓÖèÖ^÷yİÆ›½¨ÛvÇWt/Õ¾İŞÑ|äqİÛ?³{ÏwgÿÉ‰`üÚÏ|³»ë¿(Ÿ3 ˜Inºeù³m.xÖÚë•	   ˜}Õïò¡]V‡  `œ‚åtôIg•×a¾Xg½ç—{®‡Ø–”`y´­ñŒµºM7ß¦{Ùî¯ëŞpØ1İ»Şû¡î¬‹nî®ùô×'şzsõœÀÒ´Ñ¶*FÍ‚e   ]Õïò¡]V‡  `Ü‚åô¡Ë>]^‹ù`çnX>æ¹.ÿ¢iõxaå_ı­Ş'Ìù›wÛëİ±§œÛ]zÃËç ˆü/T?ƒæÁ2   Œ®êwùĞ.«Ã   °ø‚ñ–ó¯“Ş|ïËë1×­»şFåcëËğïV[ıéåû„¹iå'¯Úí°Ë«º·¾ë÷»¯»§|Î`:6yáÖåÏœ¹àÙÏynù˜   €ÙWı.Úåã‹  `Ø¾ı»S_ÑuÜİ!o=¡ÛçGvÛï´ÛÄ_ê]n¹å›ß/ç‹-¶~Iy=æºĞO>ãâîÇŸÙ|äqİëßxT·×~‡v¯ÜãõÕ³º£hõßí°óî÷ÿØ“Ïé>rùgÊÇãè¬oê{UwÒitÇ¼çìî'5ñşÆ7ÿ^÷š×Ş½t·}ºm¿K÷´5Q¾¿mk®µv÷»G¼»»òæ/—Ï? <–ÿÉ‰ÿpàáÇv{¼în³-·/æŒ¢–«Ç   Ì¾êwùĞ.«Ã   ĞçŞoı²;çâ[ºƒŞô®î9ënĞü®9×½ö€#ÊÇ=Ÿİ|ïÃİïŸ{u·ÿÁoï6İ|›òºÌ´µ½n·ëûwGŸğÁîü«?ßİı_”÷˜¾ë?ûíî}g]92ûËÌsÏ·Ú¡{÷{?Ô}öKQ>¿ °¤nğÏ»Eç\ÕíôŠ½ÊŸ9£B°   £«ú]>´Ëê0   LÇi¾¶ÛjÛšß9ç²üËcÕcù<ó/0W×fYÚpÓ-»7½ã”î’OÜ_Ş/`Ù8â§N¼ÿª÷%£-ÿRşU·|µ|^`:.¹ş¾‰4Xı¼™mk¯³~yŸ  €ÙWı.Úeu   ãŞw^·ÚÓÖh~÷œ«>|ùíåã'ŸüüwvL·Â
+–×hiÙ|ë»Ó>tMy€™sä1ï/ß£Œ¾ü¯ÜõõŸ—Ï+ LÇş‡¼£üY3›Ë   0ºªßåC»¬  Àãuıg¿İm¶ÅvÍïŸsÑj«?½»å?*ç¸ùø§êÖYïùåuz¢}Û‰å÷fÇïŸ{u·üÂ…åûuÔ<mgt½`«	›lö¢î[lÛ½p«º-¶~ÉÄ_şßzû—vÛîøŠnûvë^¼Ë«º—¼l‰ÿä}ş¿·yñË»-·Ù©ÛlËí»6Úlâ¶ªï1—¬ñŒµºãßÿÑòy€éØå•{—?kf‹`   FWõ»|h—Õa   x¢vÍÍï sQÆlÕãG7ßûp·éæÛ”×éñ:é´ÊïÌ®Ó>|mù5û½ñ¨òş?^÷}çWİÍ÷ü »ô†û»“N¿pâ/Ln³ÃËºÕ×X³üş£jó½Ø% €'ä²¾XşŒ™-‚å_ww~õ¯&şáUşC¬üGY[m·swİíß(¿   fRõ»|h—Õa   Xvİó€æ÷Ğ¹hŸ7¼¹||ãè®¯ıÍÄÿ¼ºNÓuÌIg•ßG¼ã½å{w”,í`¹Oş¥ùwrn÷òİ÷íV{Úåı5o}×ï— –Ä®{î_ş|™k¯û¼ò>“×pDs]Ë   Œ‚áßW'µËê0   ,-/İmŸæwÑ¹èØ“Ï)ß8º÷[¿ìÖ\kíò:-©ŒªÛFËZÏ^·|Š™–‡ş‘ë&âåê~’Wí}ĞÄ_® ô9ù—”?[fÃ¸Ëö¯Ë   Œ‚êwÖĞ.«Ã   °4åª¶út®ùÈåŸ)ß8:é´Êk´$VYõ©İ_ùiy»Àh9ä­'”ïãQ1›Áòoİpç÷º7¾ù÷º'=i¹ò>‚^°ÕÄ_ˆ®î? <š[ïÿqùse6Œ{°¼Í‹_^^Á2   £ ú5´Ëê0   ,M7ßûp·Ö³×i~'k¶Æ3ºÛüóò1£­¶İ©¼Nåo>®¼=`ô\rı}åûxTŒB°ü[·|áGİş‡¼£¼Ÿ£`ÅWêŞöÇËû f³-·/®Ì´ç¬»AyÿÆÁûÎº²¼&I°  À(¨~gí²:   KÛùWİÙüN:e¤[=¾qtöÅŸ*¯Ñc¹ù”·Œ¦UŸ²Zù^£,ÿÖ9ßÒ­óÜËû;
_ô±ò~@å¯Ş¯üy2ÓÆ9Xîûß‚e   FAõ;kh—Õa   XŞsÆEÍï¥sÑ>o8²||ãæ¾ïüª[~ù…å5z4/ŞåUåm£k“Í^T¾ŸGÁ(ËéŞoş²Ûk¿CËû<
{Uy¿`ØAozWù³d¦k°ü–cß_^ß,  0
ªßYC»¬  À²òÆ7×ün:½û½*ß¸ÙaçİËëóhNXt~y;Àèzå¯/ßÏ£`Tƒåß:ôm'–÷{ä_‚®î3 Lõ®SÎ-Ì´uÖ{~yÿæ³ïú~·Â
+–×ã·Ë   Œ‚êwÖĞ.«Ã   °,íòÊ½›ßOç¢ó®¼£||ãäè—×æÑÜşàŸ—·Œ®QnG=XNï~ï‡Ëû>ÛVZùÉİÅğ…ò>Ào->ÿÆòçÈLÇ`ù¥»íS^‹©Ë   Œ‚êwÖĞ.¿øGq    fĞç¿ö×İzÏÛ¨ùu®YãéÏœx,ÕcWßö‡åµ©<óYÏ)omÇrNùÁrqŸGÍëzKyÿgÛ3ùìîÚÛ¿^Şg Hç_}gù3d¦MËÅı›¯NøıóËë0ìºÏ|£<   3©ú5´Ëê0   ,k—ßø@ó;ê\ô¢íw)ß8Y~áÂòÚÛaçİËóÀh;uñåå{zÌ•`9mòÂ­ËÇ0Û6ŞìEåı€tÑÜ[şü˜ië<wÃòşÍG7Şıınå'¯Z^‡a‚e   FAõ;kh—Õa   XÖ>ıàOšßQçª}zkùÇÅs7Ø¤¼.Ã|Ó±åy`´-¾p4şSğ•¹,_~ÓƒİïüÎï”c¶òÖÊû —}ò‹åÏ™6NÁòKwÛ§¼Á2   £ ú5´Ëâ‹  `¦İ>‚åô{§~¸|œã`çW¾¦¼&ÃŞıŞÿÇŞ_¿[RŞiãöKãîî.»;‚»»»KpnÁİ=‚kˆ@±±'“™gf’Œ~_ıê­Úùò¤C],öî^U«Öîó‡ó˜9nÖ§êºkÑ{w…kÕº>ÎİvÛ#¯Ç?Ó]°×ÇÅÌ]uÆE7Å}tÁm¿30m»ÿÙãï¶-¹ô
1ßxsö%·Æıª°œ   mJ÷¬¥úb  €¦·Ârå–^Š{ï8üÔx=¾ê¢kïó@·=4é“øgº†­°\Ùyƒã^m•Õ×y˜¶uåïÓBaùé×YÌ6Ûqÿ_Ga  €.H÷¬¥úb  €¦ÇÂò|,\¼õÙ¿Çıgg_ú½x=¾êº;Ÿ‹ó@·=öâÏâŸé.ÆÂò#Ïÿ8î¥?ş¼˜€i×£/ş4şÎhÛ’Ë¬ó'›o»KÜ{/
Ë   tAºg-ÕÓ0   4m<–+ël¸EÜïxvã=ÏÇkñUw=ñvœºí©W?¦»`Ë•-·ß=î§î|ìÍ˜€i“Âr;;ı²¸ïo¢°  @¤{ÖR}1  @ÓÆka¹²çÇÆ=WU¹-]‡¯zü¥Oã<ĞmÏ¼ñ«øgº†µ°|Ó½/ÄıtÁ–Ûí30mRXnŞ­¾÷<
Ë   tAºg-ÕÓ0   4m<–+g|çÆ¸ïñèÁI?Š×à«¾ÿÖoã<ĞmßûoâŸé.ÖÂreâšÄ=uAUNK™˜öt¥°¼Ô²ã³°üæÏşT,³ÜÊqÏ£¡°  @¤{ÖR}1  @ÓÆ{a¹rÛ#¯Ç½7O½öEÜÿW½ôá?Æy Ûºüóz˜Ë\uWÜSì{èI13 Ó…åfí¼ÇÁq¿£¥°  @¤{ÖR}1  @Ó¦…Âòü.÷>ŞŒö½¬–æn«>lşLwÁ0–+/¶dÜ× Í6ûœÅ?ıcÌÀ´Ea¹9g]|KÜëX(,  ĞéµT_LÃ   Ğ´i¡°\YwÃ-ãşÇ“şôqï“›~úâ,Ğ}/ôûøçº†½°\åOûê‚SÏ¿6f`ÚÒ•ÂòÒË­ó«ûhä)íu,–  è‚tÏZª/¦a   hÚ´RX®ìuĞpêFã›şc{õ´Î4tŸÂrs.üî=q_]°ì
«ÆÌ L[–û¯ú»Õr+MŒû+…e   º İ³–ê‹i   š6-–+g^ts¼ãÅÌ³Ì÷ı¥¹æ7Îİ§°ÜœÇ_ú4î«+î|ìÍ˜€i‡Ârÿm¸évqSBa  €.H÷¬¥úb  €¦¶°üøËŸÅõatÏ“ïÄk1(,Ãø¥°Ü¬ùæ_(î­N8óŠ˜€iGW
ËË,·rÌ7lvŞãà¸¿)¥°  @¤{ÖR}1  @ÓF[X®^Ûå¯Ì‹Z¬vÆ‹o*,Ï9×<qè>…åfmºÕNqo]PeK™˜v(,÷Ï¡Ç÷65–  è‚tÏZª/¦a   hÚX
Ë•=86şóa³ŞÆ[ÿÕu/fy–¸ß/),ÃğRXnÖÑ§\÷Ö€ÂrœqÑMq_SKa  €.H÷¬¥úb  €¦µ°\Yeâ:ñ5Ãfïƒÿ«}3Î4SÜë—–ax),7ë¦{_ˆ{ëŠûy?æ`Ú °<õ¾{ë“qOı °  @¤{ÖR}1  @Ó¦¤°\™nÂ„øºasÖÅ·Ôö6Ìf˜aÆ¸Ï/Í¡°CKa¹Y¯}ò/qo]qêù×ÆÜ L:SX^~•˜¯ëª&Í0cï{¥©¡°  @¤{ÖR}1  @Ó¦´°üÒ¿‹¯FNúQmÃjÂ„éã¿¤°ÃKa¹y3Ï<KÜ_l½ã13 Ó…å)wÛÃ¯³Î6{ÜO¿(,  ĞéµT_LÃ   Ğ´)-,W®½ãÙøÚa³àÂ‹Åı£é¦›.îñK
Ë0¼–›7ßüÅıuÁr+MŒ™˜6(,O™»Ÿx»˜kîyã^úIa  €.H÷¬¥úâ{å‹   m/Œ²°œf+‡^|ı°Y“mâş†MÚÛäæœk8tß+/,§ÌÃfÉ¥Wˆûë‚ùX8f`ÚğXG
ËË®°jÌ×E<ûa1ÿ‚‹Ä}ôÛ£/ü$f   €6¥{ÖR}1  @Ó¦¶°\Ù`ÓmãÌ°Ùçàâş†IÚ×äæ˜sî8tŸÂróV™¸NÜ_L˜0}ÌÀ´Aaylªñ"‹-÷Ğ…e   º İ³–ê‹i   šÖÂr¥zroš6ç^v[Üß°H{šœÂ2/…åæ­·ÑVq]QıÎN¹ÿ–GïÉW~Şú·&(,  ĞéµT_|ï‹r    ZöÂû£,,‡Ù¯JsÃèñ—?ûi?“),‡9 û^ù¸ã…åyØl¹ınq]ñğ>‰¹ÿºRX^nÅÕb¾®xà¹‹Å–X&foÒ£/ş$æ  €6¥{ÖR}1  @ÓúYX¾û‰·ãì°Yh‘Åãş†AÚÏä–ax),7oç=ûëŠ[x)æ`üSXşf·?òZ1Ï|ÄÜMSX   Ò=k©¾˜†   iı,,WN=ïš8?l6Ül»¸¿®K{™ÜìsÌç€îSXnŞ¾‡÷×—^ÿ@ÌÀø§°ÜÛõw=WÌ<Ë¬1s–  è‚tÏZª/¦a   hZ¿Ë•ívŞ'cØìwèIq]–ö19…e^
ËÍ;ì¸sâşºâ‚+ïŒ¹ÿºRX^~¥‰1ß ]vÃC1k›–  è‚tÏZª/¦a   hZ…åÊâK.3l.ºæŞ¸¿®J{˜œÂ2/…åæísÈ	q]qÅMÄÜ Œ
ËÙyWÜs¶Ma  €.H÷¬¥úb  €¦5UX®¤ã£gßøUÜ_¥ü“›mö9ãcóÔkŸ7Şóƒâ’ë(Î¾ôÖâ¤³¯*?á¼‘§rï²÷aÅ>_yÒÅi\W\|íıÅõw¿¸çÉ·‹§^ıEñîçÿ	ßDa¹yßŞó¸¿®¨~–¤Ü|½—>ü]qção—ßøpqÎ¥ß+N8ëŠâĞcÏ.ö<ğØ‘Ÿ×ÕÏíÃ?wd½úy^½ö‡?ù×x,Fç­Oÿ4òûîöG^+.½şâÔó¯-9õ¢â £Îø?×ıÀ#O+=í’âìKnyoîöƒV~?|îÕÅq§_:’§RınN¯ƒ.RX®;å¼«cÆAPX   Ò=k©¾˜†   iM–Ÿxå³x¬a³ÈbKÆıuQÊ?9…å±©
¢U‘jÿÃO)6ßæÛÅ2Ë¯RÌ8ÓLñÚÖ„	Ó+­ºV±Ó§wMñ½‡^)^ıäŸãùar
ËÍÛzÇ=âşºâöG^¹ù³ª${ÖÅ7;îz@1q­Š¹ç?^ÇÑXréŠ­vØ£8ã;7*¡}ƒÇ_úYqÁ•w»í{ÄH‘1]ÏÑ¨~?®¸ÊšÅ·v;pä?Õû™Î7¥ª¢òWÏY•Ók¡‹–ÿZõˆ”oPü®    Ò=k©¾˜†   iM–+]so<Ş°Ùh³íãşº&eŸœÂò7«Y=-yu6×°)U9®*}ŞùØ1(,7¯úYŸö×<÷aÌ=-«:½ÿa'Ë­¸Z¼fı²èâKr«'0§ÓšêçÑ™İT¬½şfñzõËJ«­]vÜ9Sı»ñŞ§ßÇ?êäïÄ×Cu¦°¼òê1_[ªû×7Û.f$…e   º İ³–ê‹ï}ñ?   ĞºŞÿ_µ{Ô$ÍV×~5¥ªBTÚ_—¤Ü“ûsa9ÏNËyşÇ#_U¿àÂ‹ÇëÖ¶êÉuÇviñôë¿Œy™6½òñ?Å_ºàÏ…åœ{˜´ıA…±úáOş-æÖ¼ğŞ?Ç~Y±ä2+ÆëÔ´‰kmXœ}é­Å»ŸWÿ g¯ªß—»ïwT¼.M«ş|^tÍ}1W/Ï¿û÷#
JÇ<àˆSãL¿ÜtßC©_?kÒ±‡Á»¿ø¯¸ŸA{ìÅŸÅÛöçÂrÎØ´;{³Xbéåc®A«>¨{ÛÃ¯µ"]   ¨¤{ÖR}1  @ÓÚ(,WV™¸n<î°¹üÆGâşº"eÜl³Íç¦EUÑìü+ï,ÖİhËx­º¢z‚ÜWİ÷À´Ea¹y3Ï<KÜ_T¨H™§%·?òÃb‡]÷×g\x±âÔó¯‹YÇ›G_üi±íNûÄëĞ¶ª¬xò¹Woÿü?bÖ¯Zã­ãq*ÕS³ÓL¿,°Ğ¢ñ¼]wŞwÄıŒÕœsÍßu—ŞğPÜÏ Më…å¿{O1İ„	1Ó´&]   ¨¤ûÈR}1  @ÓÚ*,Wf™u¶xìaS=Ù1í¯RŞÉ),ÿÙi\?R4K×¨«–Zv¥âÜËn‹ûaÚ °Ü¼fœ)î¯ÖÙ`‹˜yZpÿ³v¦,›T?Ÿ«À¤ìãÁ	g]ÑÉ’àÜóÌWuÒ…ÅkŸüKÌıÒ‡ÿXl½ãqöKÛí¼oœí—®|sÃXõ«°<×ÜóÆãwÂro+¬¼FÌ×¤Ã?/f™V¥k   •tYª/¦a   hZ›…å—?ú}<ö°Ydñ¥âşº åÜ´^X®Ì¶Ìò+Çk3,–]aÕq]Œãë),7ëÑ~÷Ö»îÓì“`»èñ—?+¾½ç¡ñztÑæÛ|»xîÍßÄ½£êú¯¹Ş¦q¯]2ël³Ú9õââ £Ïyêõ>‡œXÌ0ÃŒñõ“ÛdËoÅ½÷Ë"‹-ÏÛuıú{ÆÜóÎßu—İğpÜÏ M‹…åêƒ[í°GÌ1-K×
   *é>²T_LÃ   Ğ´6Ë•ëî|.Øl´ùqƒ–²NnZ-,ßùØ›ÅÄµ6ˆ×dXU_Çı«ï‰ûe|RXnÖ7=÷ÖÓÚN=ïÚ¡üêÿê‰²]s_ÜÓ0¹ægF`œö8¬µŞfqÿı²èâKÇóv]¿~ŞÌ7ÿBñø]§°Ü[[…åîùÁÈES†i]º^   PI÷‘¥úb  €¦µ]X®~ÂøøJß85îoRÎÉUO!LsãÙ)ç]¯E“¦Ÿ~†¸Ş„õ6Új¤À’öÎø¢°Ü¬£NşNÜ[W<õêç1÷xóôë¿,6Ûzçx¦ÔÌ3ÏRì°ëşÅ'?òAïy¾8ó¢›F~o±İ®”*÷<ğØ¸¿apÒ9ß{V\eÍxúåª[Ÿ(=öœ‘Ï–_ibÌĞE\yWÜÏX]yËãÅ!ÇœUl¾í.#ß‘ÎÕE—ßøHÜÏ u¥°ÜôŸ›Ê±§]ÏÍŸ¥k   •tYª/¾_¾   Úöâ(Ëivjl¸Ùvñ<ÃæšÛŠû””qrUa9ÍG¯ò/Åv;ï¯C?­2qİâğãÏ+n}à¥âñVüğÇÿ:rş—?ü]ñèó?..½şÁâà£Î(–\z…8ß3Í4sqÎ¥·Ö®ãË«.,ï}àq1ó0Ùn§}âŞº`ñ%—™Ç›s/»mä› Ò5˜UA²úĞÊ?ù·x¾/½ıéŸŠÃ=§˜nºéâq¦ÔVÛï^¼÷ùÇsvÕß¹1î¥ª÷¶*£ïºÏ#×{¿CO.¶ùÖ^Åël\,¶Ä2q¦imÿÙzâ¥O‹ãO¿¬Xmõb®¨
Ë)ÿÔzô…ŸÇœzq±òÄuây»¢*,§üƒVı=3åm[UXNùú¡º7İlëâyù‹tí    ’î#KõÅ4   MTa¹²ÀB‹Æs›×~ôÏqƒòMnZ),ßúÀËÅâK-¯A¿l¾Í·GŠ#éü_çª[/–Zv¥x¼~Øq×ı‹×>ùßñÜ?…åfUÅÅ´·.ØyƒcæñdßCOŠ{ŸRÛî´wñÖÏşÏõu{ãWÅÆ›ï7¥ª2î¤·ÿ&¯kÎ»üö¸‡©U•”oº÷…xÎÉ=óú/‹³.¾¹Ør»]Gş¾’ÕoóÎ·`ÌÒ†‹¯¹¯³O]¾ğª»cæ~ªÎ±Ìò+ÇóÚ7=3Úx/,_{Ç3Å<ó-ÏÉ_K×   *é>²T_LÃ   Ğ´A–+é\ÃfÑÅ—{„”or³Ì:[œOª§âM7aBÜ?Ì9×<Åw¾{O<÷h¼û‹ÿ,öØï¨xì~Xdñ¥F
éÜ7…åæ<øÜGq_]qõ÷Œ¹ÇƒêCU¡5í{JuÒ…ñ\£ñŞ/ş«XwÃ-âq§TõšêÉºé|]Q}Ğ'eŸË.¿JqÛC¯ÆóÆ©ç];òw¬tì~™yæYâ¹ÛT} *e¤©ù{ÎX¼ù³?ëlĞß?oı °ÜÛJ«®óMêçv:Yº†   PI÷‘¥úb  €¦º°|ÏoÇó›M¶üVÜ_ÛR¶É÷ÂòYİ÷İ/UAãÑçÏ=V‡{v<G¿sÊEñ¼/…åæÆåq_]°Ğ"KÄÌãÁıÏ|P,»ÂªqßSjïƒç‹şø_‹U&®?¥º\Z®Dİï§ÿWO¸~ã§ˆç«êwûRË¬ÏÓéœm[k½Mc¶Ai«°\©JË/ºDÌ1(WŞüXÌ:hã¹°\ı="‹,]C   ¨¤ûÈR}1  @Ó]X®œvşuñœÃæ°cÏ‰ûkSÊ5¹ñ\X®
ºiÏı²àÂ‹”ºÒ¹§ÔIg]ÏÕ/‡wn</ÃIa¹9]|Âç—:òô˜yØİ|ß#¿“Ò§ÔZëmÏ5%šø@UWKË;í~PÌ;¥ª'f§óL­<-oj¥sµ­z¢oÊ6(]İ^a¹Òôß‡ÆJa¹7…åÁK×   *é>²T_LÃ   Ğ´.–+Ûí¼o<ï°¹éçãşÚ’2Mn¼–«²xÚo¿Ì9×<Å¤·~Ï=µ6İj§xÎ~©
^é¼…åftùºVúşÇ1÷0«ÊÊ3Ï<KÜï”š0aúâÉW~Ï7¥š(QV¥åçŞøu<ß œuñÍ1ç”Zaå5Š·?íï‡{&wÍíO‹,¶d<÷”Jç„V^=æ„¶Ëïüü?ŠÙf›#f„«ny<æ´ñ\XŞuïÃã¹ÈÒ5  €Jº,ÕÓ0   4­+…åÊ2Ë¯Ï=lúıŞ±Hy&7ó,³Æ¹avö%·Ä½öÓWŞÏİO¿úyßŸ2úU{t|<7ÃEa¹§œwMÜSl·Ó>1ó0k¢¬\iê[V[c½x¾©±şÆ[Çsµíg?,¦›0!fœRÕ“‚Ó¹ú©úY¸ı·û÷A·÷~ñ_ñ<mëRaó¢kî›Ô¥'İ+,÷¶òjkÇ|ScÇ]÷/æg¾)2¡Ï?Ç¦TÊÖ”t   ’îYKõÅ4   MëRa¹2ıô3Äó“Å—\6î­)ÏäÆ[aùº;ûì§mwÚ;»ŸN<ëŠxî~Úmß#â¹
ËÍXl‰eâºàşg>ˆ™‡USeåyæ[ xógŒçœZWŞüX<çÔÚïĞ“ãùÚ´ÚšëÇlSjíõ7çiÊAGsŒÕ+ı>¿mç\zkÌ7ƒ(,Wß‘²ÂU·<3Úx.,OZ4æl[õş¤|   Ğ¦tÏZª/¦a   hZ×
ËO½ú‹xşa³éV;Åı5-e™Üx*,?ôı‹9çš'î³Ÿn¼çñüı6ÿ‚‹Äó÷ÓÏÍpPXî¿s/ı^ÜOì´ûA1ó°zğ¹{šü)ç^ÏÙ/U7wj]påñ|mh¢ˆ]}ãA:W“N:ëÊ˜e,ıá/ã±ÛV=ñ:å„‹¯¹/flRõtî”e–{SXÎxéÓ˜   Ú”îYKõÅ÷Y   @Ë^ü`”…å0Û”ê©n)Ã°9æ”‹âşš”rLn¤°æ†Í[Ÿş©Xv…Uãûiù•&Æó7á¤³§¾t5ç]q{<?İ÷ê:\X>è¸˜¹ë–[qµ¸Ÿ.xô…ÇÌÃèŸı¡XaåÕã>§ÖÒË­ÏÙO—Şğ`<÷Ôš0aúâÑÏÙ´µÖÛ4fšÏ½õ›x®¦sêÅ1Ïh=üƒOâqÛöÆOÿ-æ„K®»?flÒ#Ï³Âw¿÷DÌ8h¿Ô‘ÂòÄub¾AéLaùåOc>   hSºg-ÕÓ0   4­‹…åJU~K9†Ím¿÷×””arÕWñ§¹aóí½ûë·‘§v†ó7áİÏÿ³˜mö9c~š~úŠû0f Û–ûëœKo{é‚ãÏ¸<fV[í°{Üg?\tõ=ñœı6ÇœsÇóO­·Ø1¯I×İùlÌ25ı÷‹Ã;'æ;{#sÚø{Àh¢°üÚ'ÿ³‚Âro
ËÙ¯|ó  @›Ò=k©¾˜†   i]-,WV_{Ã˜eØ¤½5%rã¡°\}…~Ú[ªòLÊĞ”İ÷;2æè·UV_7ŸnSXîŸê‰ª3Ì0cÜË ­¹î&1ó°:ø˜3ã>û¡zBv:gvİçˆ˜¡Î¼è¦xÎ¦l°É61ÇÔXdñ¥â¹Ú´Ã.ûÅlßäú»¿7‹/¹lÌØ¶K¯ ækÚìsÌó´íêÛŒùMa9SX  €¿H÷¬¥úb  €¦u¹°\™{ùba²øRËÅ½5!rÃ^X~òÕŸ³Î6{Ü[¿­µşf1C“näõ˜¥	»í{DÌ@w),÷ÏZëm÷Ñw?ñVÌ<ŒÎ¿â¸Ç~9ä˜³ây›pë/ÅıPı^{îÍ_ÅóöÛ-÷¿3ôCõ3*³MSòa·Ënx(kºòa½A–—]a•˜§m
Ë½uíƒo]),W÷))   ´)İ³–ê‹i   šÖõÂò+ÿ>æ6›m½SÜ_¿¥sOnØËlºmÜW<ñü˜¡i+¬²FÌÓ„³/¹%f ›–ûã £N{è‚Ó.¸.fF“Şşm1ç\óÄ}öË¿Ïİ„7~úo1C¿ló­½âyûm»oïÏßU©;³MMúQÌÖË¹—ß5[l»KÌØ¶Kox0ækZ›ÏëåšÛŸŠùMa9ëJaù©×~ó  @›Ò=k©¾˜†   i]/,Wª¯êN™†Í‰g]÷×Oé¼“›i¦™ãÜ08éœ«âšò½‡^‰9š¶ÏÁ'Ä<M˜a†‹çŞüuÌA÷(,O½³.¾)æï‚½<6fVU7í³_YlÉxŞ&­°òê1K¿ÜóäÛñ¼ıÔä7WT¿§Ó9ÛVı}+åû:'ŸûİxœAhúÏÍhª°¼éV;Å<mSXî­k…å^,æl›Â2   ]îYKõÅ4   M†ÂråˆÎ‹¹†MÓ_õŸÎ9¹a-,ÿàİ¿ùÊü´§&Ì<Ë¬1G.¿ñá˜©)»î}xÌA÷(,Oó®¸=fï‚õ7Ù:fVç_qGÜg?¢à]ı¼LYú¥é§,ßşÈkñ¼ıR=åôŸı!»mÕŸ©”1Ô7*$Û7øì±¸ì†‡b¾¦m¹İ®1OÛ®¹ıé˜oĞºRX^Ua9zúµÏc>   hSºg-ÕÓ0   4mX
Ë•·Ø!f6ioı’Î7¹a-,{¯Cã~š²ÑæÛÇmxîÍ_ÅLMºóñ7cºEayÊ]tÍ½1w¬¶æúÅËıcÌ=Œ&½ó7ÅœsÍ÷ÚO7Ü=)¿Iç^ö½˜¥Ÿš|Êò¡ÇÏÙO{||<wÛîöƒ˜/Ùÿğ“ã1a§İŠÛ6¨ÂòÖ;îó´Ma¹·î–9Û¦°  @¤{ÖR}ñƒòÅ   Ğ¶—FYXN³ƒ°ø’ËÆ|ÃdÉeVˆ{ë‡t¾ÉÍ8ÓLq®Ën¹ÿ…¸—&í{È‰1K[Ú~R\UĞN9è–×:^XN™»àÄ³®ˆ™» *+¿òÑ?ÆÜÃj§İŒ{í§9ç7»i·=ôJÌÓOÛ~k¯xî~Xyâ:ñœıV]§tş¶öÛ9ª'g§ùAØeïÃbÆ¶Ußöò5m‡<aúº;‰ùí‰–S¾AYh‘%bÎ¶U…å”   Ú”îYKõÅ4   M¶Âr%å6[ï°GÜÛÔJçšÜ0–'®µAÜK“N;ÿÚ˜¥-›lù­˜«IWÜôHÌBw(,Mõ´òM·Ú)æí‚ñXV~hÒÇq¯ı¶Õö»Åó7íé×~óôÛ3¯Ï?5&½õ›x®&,±ôòÅ«ÿ>æhÛŠ«¬3N®É’øXí±ÿQ1cÛõw‚®<aZa¹7…åLa  €.H÷¬¥úb  €¦caùŞ'ß‡Íé\÷75Òy&7l…åó¯¸#î£i×ŞştÌÓ–}>!æjÒò+MŒYè…åÑ»èê{Š9çš'fí‚êC+Õû™²³6®\©œ›Îß´w>ûSÌÓoÇ~i<ÿÔ¸ì†ã¹š²á¦Ûo}úÇ˜¥MWİòXÌ7¹.}Ë@õ³4elÛ 
Ë]yÂôõw>óZW
Ë«­±^Ì7(/ª°   _J÷¬¥úb  €¦ca¹R•}SÎaóĞ÷?Šû›Ré“¶Âòr+¬÷Ñ´GŸÿ$æiË şı>åÜ«cºAaù›}ÉÍÅ²+¬3vEeÔ.hëéÊ•Ëox(fhÃ<óÎ3õÓ
+¯Ï=5Î¸ğ†x®&­±öFÅ‹ïı}ÌÓ¦õ6Ú2æûÒëlçAaYa¹…åLa   ş"İ³–ê‹i   š6šÂòt&ÄÙAÛáÛûÆ¼Ã&ímJ¥ãOn†‡§°<¨§+WŞüÙb¦¶TOxN¹š¶Ô2+Æ<tƒÂrvÿ3ïÇŸqY±ÈbKÆl]Qªn¼{RÜÃxĞÖÓ•+ƒüPI[¤¹ëñ7ãù§TõTêt¦-¹Ì
Å}O¿3µå¦{ŸÙ¾Ô¥oPXVXîEa9ëÊï…e   º İ³–ê‹i   š6šÂò„	ÓÇÙ.XyµµcæaÒÏ’h:şä†©°<¨§+Ï¿à"1O›ªrWÊÖ†ëîx&fbğ–ÿìoÿ¶¸ì†‹]ö:´Xh‘n<Uñ›l±í.#¿oÓ~ÆƒÇ_üiÜwf˜aÆ˜¡-ëo²uÌÕoû||<ÿ”Úó€câyÚ0ó,³WİòXÌÕ–Í·Ù9f«,ºøRqf––{éLayÍõc¾A©ş§œmSX   Ò=k©¾˜†   i£),O?ıq¶+f›}˜{˜ì°Ë~qoc•=¹a),òéÊ+O\'fjÓ¤·~³µa«v™¼i©°üâûÿP<<éG#O?ñ¬+ŠoïyH±úÚsÎ=o<WUO¡¬öö8rÌ™qÿMXaåÕc†¶l³ã1W¿õ»DÛVî^N>çª˜­·?üjÌT™{ùâÌ (,+,÷¢°œ),  À_¤{ÖR}1  @ÓFSXôÓ¿ÉÓ¯ı"æ6ç\zkÜßX¤ãNnX
ËU"åoÃÆ[ì3µé½Ïÿ3fkËsoş*æb°º\X®ÊÄ§œûİâŒo(Î¾äæšÓ/¸n¤x|Ì)‡nqĞQ§ûrb±ÇşG;í~`±áfÛ+®²ÆÈÎ«§ú§s“êÉÏıø™>,\x±xš°ı·÷Ú²ëŞ‡Ç\Mxô…ÇSbİ·ˆçhÛ^óµáë²Ü¥¿),+,÷¢°œ-ºøÒ1gÛ–  è‚tÏZª/¦a   hÚ¨
ËCPr½äÚûböaS}­~Úßh¥cN®ëåóÊıÏ¼³·e§İŠ¹Ú6È'‡uÒ1ƒÕåÂ26ïüGx~ñşçÿßÃñèò×¢)ÇzqÌÑ–8%æjBõ€”aJ,¿ÒÄxAXcíŠGŸÿ$ælÒm½óTÒëAaYa¹—®–'®µAÌ7(
Ë   ğéµT_LÃ   Ğ´Ñ–gšiæ8Û5])yL­´·ÑJÇ›Ü0–÷?ìä˜½-qjÌÕ¶Z4ækCUüH™,…åîZgƒ-Šï\}O|ßÆ»M¶üV¼&M©> ”r´åèS¾s5a«íw‹¦ÄëlÏ1(3Ï<KñïŞ³6i£Í·yÒkAaYa¹…ål±%–‰9Û¦°  @¤{ÖR}1  @ÓFSX®J%i¶‹ÖÛhË¸‡a²ô²+Å½F:Şä†¡°<ßüÅìm9ñÌ+b®¶-½ÜÊ1_[®ºå±˜‹ÁQXî–Yg›½ØëÀc‹ûŸy?¾_Ó‚§^ıy¼6Mºîgb–¶œzŞ51Wæš{Ş˜aJTåçtAÛsÿ£cŞ¦\sûS1Gzí (,+,÷¢°œ-¾ä²1gÛ–  è‚tÏZª/¦a   hÚ¨
Ë³Ìg»j¡…û&Ua#íí›¤cMnúégˆs]QdRî6Å1[ÛV[c½˜¯-;í~PÌÅà(,wCõDáêçÄÛŸş)¾OÓ’SÎ½:^£&İşÈk1K[ª÷>åjÊ=O¼sŒÕûß+­ºVßö9ÕÃ¾šáıÏÿ+¾¶m
Ë
Ë½t¥°¼úÚÆ|ƒ¢°   ‘îYKõÅ4   MMay–Yg‹³]õêÇ¿û6SòUéé8“ëzay»÷‰¹ÛT=}1ekÛ Ÿ^ÿS.Gayp6Ú|ûâÂ«îŠïË´¬º.éz5é¡IÇ,miûƒ5§_p]Ì1VGx~<~—Tø”½ßÎºè¦Ú¹«¿;¦×¶MaYa¹…ål‰¥–‹9Û¦°  @¤{ÖR}1  @ÓFSX®¾ö>ÍvÙw?îeØ<óúq_'cr]/,Ï=Ï|1w›îzüÍ˜­m›oóí˜¯M>÷aÌÆ`(,·ïâkï‹ïÿS¼óÙŸF~§¤ëÖ¤çŞøeÌÓ–¶ÿ~±ã®ûÇc•Jº]´õ{/ø»¸‡~ºåşş›îùA|Í (,+,÷¢°œ),  À_¤{ÖR}1  @ÓFSXmö9ãl×ÃÓ¿ÉtÓM÷öuÒ1&×åÂò¼3·í‘|óµ­*©¥|m:ñ¬+b6Ca¹}Ë®°J|/øŸ‘§Ñ§kÖ´êÏAÊÓ–[ï1æjÊÒË®sŒÕU·<ßEÕş¯¹íÉ¸ñNaYa¹…ål‰¥—9Û¦°  @¤{ÖR}ñƒ_•   Ğ²—>üæÂòìsÌg‡ÁVÛï÷4L–_ibÜ[’æ'7aÂôq®8á¼˜¹mÏVOïùÚ¶ÇşGÅ|mÚp³íb6ãµO–áà£ÏˆïÇ´n¯×«i)K›îx´ı×¼òÑïb–±¨‘İeûvRÜËxÖ™ÂòÍÄ|MëLaù®gc¾A{âån–×Xg£˜oP–\f…˜³mO¿şyÌ   mJ÷¬¥úb  €¦¦°<ÇœsÇÙa±ÜŠ«Æ}“=<&îí«Òìäº\X^mÍõcæ¶½Ü‡bX?ìsÈ	1_›f˜q¦â½/ş3æ£}].,o³ãÅ=O¾ÕÓ½O¿SÜÿÌ{Å>Ñew=ñf|O¦eƒ(ˆuáT÷>õvÌÖ¤ëï~.f«•W[;¿ËVXyõâ{½÷3),+,÷¢°œ),  À_¤{ÖR}1  @ÓFSXs®yâì0™q¦™âŞ†Ée7<÷6¹47¹®–_|ÿïcŞAxû³?ÅŒm;ğÈSc¾¶İp÷÷c>Ú×åÂrUBN™¿Îr+®ÓU³Í>GÜÇ´êÅş!^§¦-´È1O›üş‡1[“;ıÒ˜e¬ö?ìäxüapØqçÄ=7
Ë
Ë½(,gK-³bÌÙ6…e   º İ³–ê‹–/  €¶½<šÂòÜóÆÙarßÓïÄ½›ŞûÛ¸¿/¥™ÉU…å47h7ß÷|ÌÛ¶f˜1æ„Ã;;flÛAGóÑ¾×;\XŞ÷bæ¯3ˆ§ÔN­í¿½oÜË´èÆ»¿¯QÓ]|©˜§M °¼ÇşGÅ,cuÍmOÆã‹5×İ¤xxÒGqoãEW
ËWŞüHÌ×´];TXNùíÉ–×\wã˜oPºRX~æõÏc>   hSºg-ÕÓ0   4m4…å¹ç™/Î›Ó/¸.îo˜L7aBÜÛ—ÒÌä¾i~PN:ûŠ˜·mÕÓÄS¾A8ê¤ócÆ¶m¾ÍÎ1íëray¿COŒ™{9áÌËã±ºì²ëˆ{™Ö{ÚÅñú4mşyÚ4ˆÂò¦[}+f«şøŸãñ‡ÉtÓMWœqáõqãÂ²Âr/
ËÙÒË®s¶Ma  €.H÷¬¥úb  €¦¦°<Ï¼óÇÙaôí=‰{&«¬¾nÜ[%½~r]-,ï¸ëş1oÛZdñ˜o9å;1cÛª§Ö¥|´¯Ó…åÃNŠ™¿ÉF›m×eUé3íeZ²Õö»ÅkÓ´.|¨d…åV^=f™Õ‡PÒ9†Í6ßÚ³xéƒˆ{f
Ë
Ë½(,gK/·rÌÙ6…e   º İ³–ê‹i   š6ªÂò|ÄÙaUıGö´Ïa²ÿa'Ç½¥×N®«…åªŒ•ò¶­KåÜ®–+ïñŸ1#íêraùë~&}“çßùm1ÓL3ÇcvÕël÷2-Ytñ¥âµiÚÌ3Ïó´éÁç>ˆÙš4×ÜóÆ,Sâ¶‡^çFó/¸HqÕ-Å}+…e…å^ºSXŞ$æ…e   ø‹tÏZª/¦a   hÚh
ËóÎ¿`œfU	;íu˜\sÛµ}¥×M®ú*õ¯ÎZU†­r¥¼m[yµµcÆAèRaùÑç?‰iW—Ë>e…åÊÅ×ÜÙe'œqYÜË´à¥÷ÿ>^“6T¿+æ_`áKÙšöæOÿ5¾Sb“-wŒçVSú„÷.RXVXîEa9[fy…e   øRºg-ÕÓ0   4m4…åùæ_(Î³g^ûEÜë°yõãü«}¥×L®‹…åª›²B—
]*,ª¸Ä_ëraùÀ#N‰™Gk—½Çí²Ç^˜6‹üw?ñf¼4ëñ—~ß)qËı/Äs³U&®S<ğÌ{q¿ÃDaYa¹…åLa   ş"İ³–ê‹i   š6šÂrõ•ÛivØ]rİ}q¿ÃdÂ„éÿjOé5“ëbaùº;Ya£Í¶‹¡K…å*KÊH»:]X>òÔ˜y,[b™xì®Zh‘Åã>Æ»+nz8^šuÛC/Ç÷cJm¸é¶ñ<Ãl†fyb{Úï°PXVXî¥+…åµÖÛ4æ”eWX%æl›Â2   ]îYKõÅ4   MMay…³ãÁş‡÷<LÖ\wãÿ³ŸôÏ'×ÅÂò^³Â–Ûí3B—
Ë;î²_ÌH»º\X>è¨Óbæ±øŞƒ/ÅcwÙ÷2|ÎUñZĞ¬Ëo|(¾SêÆ»¿Ï3zìYqÏÃ@aYa¹…ål¹W9Û¦°  @¤{ÖR}1  @ÓFSX^páÅâìx±éVßŠû&_–sÒ?ûª¯îĞ<â”˜svÜuÿ˜qºTX¸Ö1#íêraùà£N™ÇêˆÎ‹Çï²ªø™ö2^íwØIñ:Ğ¬3/º1¾S£+åØ&l±İ.Å[?û·¸ï.SXVXîEa9SX  €¿H÷¬¥úb  €¦¦°<-|íıK/÷>LFû´Ä´ÿAÚz‡İcÎAØ}ß#bÆAèRa¹*‚¤Œ´«Ë…åC>#fk®»I<G—½ıÙã^Æ£­wÜ#^šUıNHïÇÔxïóÿÿù:+®²FñèóŸÄ½wÕ´^XŞe¯Cc¶),÷Ö½Âòj1gÛ–  è‚tÏZª/¦a   hÚè
ËKÄÙñ¤ËEÀ~Kû¤U&®sÂŞ3B—
Ë‹.¾TÌH»:]X>æÌ˜yJ<ùÊ§ñ]¶É–;Æ½ŒGÕ×Ó5 Yûzb|?¦Ö5·=Ï7^Ì1çÜÅuw<÷ŞE
Ë
Ë½t¥°¼öú›Å|ƒ¢°   ‘îYKõÅ4   MMay‘Å–Œ³ãÍM÷LŠûoÒŞ©*Ã¦œƒ°ßa'ÅŒƒĞ¥ÂòÜóÌ3Ò®.–=ö¬˜yJsÉÍñ<]VeN{oª¿¤ıÓ¬wİ?¾ı°ÇşGÅs'ç^vkÜ{×(,+,÷¢°œ-¿ÒÄ˜³m
Ë   tAºg-ÕÓ0   4m4…åiééªGŸ|a¼ãIÚ÷ Í¿ÀÂ1ç tÔi1ã t©°<ãL3ÅŒ´«Ë…åÃ;;fÛí´w<W—=ıÚÏã^Æ“AÿÌ^oã­ŠK®»o(\zıı}S=	9½ığÖ§(–ZfÅx½Ç“cO»8î¿K––{QXÎVXyõ˜³m
Ë   tAºg-Õ?üÕÿ   ­{ùÃßÕîQ¿jÑ%–³ãÕ»ì¯Ãx‘ö<HÕ×µ§œƒğçÒeÎÙ¶cN¹(f”wşï1'íyı“ïM~Ü91óÔxëg(æ›¡x¾®ZnÅUã^Æ“9çš'î½-›nµSÌÅÔyü¥Ÿ-²x¼æãÉ~‡÷ßİ),?ó5­;…åçb¾A{òåOcŞ¶­³Áæ1ß t§°üEÌ   mJ÷¬¥úb  €¦¦°¼ØËÄÙñl•Õ××b<Hû¤êé½)ç uÒ1ã t­°üÊG¿‹9iO§ËÇŸ3O­ëîx&¯Ë=ö¬¸—ñbæYfûnËŠ«®s1õîæ½bŞùŒ×}<9ä˜3ãş»@aYa¹…ålÅUÖˆ9Û¦°  @¤{ÖR}1  @ÓFSX^|Éeãìx×¥'ÿöSÚë |ğË¯ızªPXşzÏşğ—1'íéraùˆÎ‹™ûá€ÃO‰çì²Ûz%îe<˜nÂ„¸ç¶T…Ú”‹ş¸ïéw§‰Òr—~ßONaYa¹…åLa   ş"İ³–ê‹i   š6šÂòK/gÇ»ª´“®Ç°K{”·~ö‡˜qPºT`ª²¤ŒƒòèóŸÄœ´§Ë…å#O<?fî—•V]+·«fŸc®âƒ/ş;îe˜½û‹ÿˆûm[ÊFÿTOZ^xÑ%ãµON=ïš¸ÿARXVXîEa9«¼Ÿr¶Ma  €.H÷¬¥úâGå‹   m¯Œ¢°¼ä2+ÄÙiÁ^¯É0Kû”·?ıcÌ8(UI8å„Ã;'f”g^û<æ¤=?ìpa¹é?;=÷a<o—í°Ë~q/ÃìÿKÜkÛŞıìßc>úçéW^,¿ÒÄxıÇ“ï|÷®¸ÿAÙ§#…å«n~4ækZW
Ë7Üõ\Ì7h]*,§|ƒÒ•5=ûú1   ´)İ³–ê‹i   š6šÂòRË¬g§{xl¼.Ã*íqRÆAéRaùĞcÎŒå…wÿ6æ¤=].,}ò…1s?zîÕñÜ]vÉµ÷Å½«×>ş}ÜgÛ|€¢¯~ô»bİ·ˆïÁxrïoÅıÂ4_XŞó˜§m
Ë½),g
Ë   tAºg-ÕÓ0   4m4…å¥—])ÎNK6Ød›xm†QÚß Í<ó,1ç t©°|ğQ§ÇŒƒòúş)æ¤=].,sòwbæ~Ûl«âù»lÒ›¿{Fï|ö§¸Ç¶İñğk1ÍØi·ãû0^¬¼ÚÚqßƒ °¬°Ü‹ÂrVıN9Û¦°  @¤{ÖR}ñ£_—   Ğ²W>úæÂò2Ë¯g§5‹,¶d¼>Ã&ímæš{Ş˜sF2Â‡Ÿ3Ê»?ÿ÷˜“öüğÇ.,Ÿò˜¹ß^şğ³Í>GÌĞUk­·IÜË°šnºéâ>ÛtÙõÄl4ç°ãÎïÅxqĞ‘§Å}·mŸƒ;TXùš¶Ë^*,‡|ƒöä+*,‡|ƒÒ™Âò¿ˆù    MéµT_LÃ   Ğ´Ñ–—]aÕ8;­©ş#tº>Ã&ímXhÑ˜s?ı’˜qö=ä„˜q&L˜3Ò.…å?»üÆ‡b†.;ñ¬Ëã^†Ñ,³Î÷Ø¦“Î¾2f£Y\yG|?Æ‹Ç^üqÜw›––{QXÎV™¸NÌÙ6…e   º İ³–ê‹i   š6šÂòr+*,éÒëî×h˜¤}ÒâK.sÂIg]3Â^3ÂìsÌ3Ò.…å¿Øcÿ£b.»ç‰7ã^†ÍÜóÌ÷×¦ı;)f£yw=şF§~o÷Ó6;î÷Ü&…e…å^–³UW_7æl›Â2   ]îYKõÅ4   MMayù•V‹³ÓªaÿZô´§AZcbÎA8õ¼kbÆAØ}ß#bÆAXfùUbFÚ¥°ü×–ZvÅ˜¥«^tÉ¸a³Ğ"‹ÇıµiƒM¶‰ÙhGõwÇÍ¶Ş)¾7Ãîú»{n‹Â²Âr/
Ë™Â2   üEºg-ÕÓ0   4mt…å‰qvZV=…/]«aö3HÛî´WÌ9g|çú˜qvÚíÀ˜q”»Aaù¯İõØc–.«œö2L–]aÕ¸·6Í8ãLÅ¿üï˜öTîÓû3ÌV]c½¸×¶(,+,÷¢°œUnSÎ¶),  ĞéµT_LÃ   Ğ´Ñ–WXyõ8;­«ŠÜézu]ÚË pøÉ1ç œsé-1ã t©Èıí=i—Ârİ1§^ótÙ Š€ı²Öz›Æ}µíÖû_ˆùh×M÷L*YlÉø«Û~%îµ
Ë
Ë½(,g
Ë   ğéµT_LÃ   Ğ4…å)W•gy–xÍº,íeN=ïê˜sÎ½ìÖ˜q6ßæÛ1ã qÂ¹1#íRXÎÖİhË˜©Ë^zÿïã^†Á–Ûï÷Ô¶Ã;;æ£}¯~ü»N}Èfjí²÷¡qŸm˜ÖËÕ¤R¶),÷ÖµÂòjk®s¶Ma  €.H÷¬¥úb  €¦),O›î¯Y—¥}Ò•7=sÂé\3ÂF›o3B—ŠÜÓ2…åì¹7~YL˜0!æêªÍ¶Ş)îeì¾ïqOm[}íc>§úšŞ«a3ÓÌ3oöÇ¸Ç¦),+,÷¢°œ),  À_¤{ÖR}1  @Ó–§ŞIg_¯[W¥=Ò½O¾sÂ	g^3BUIáÖ^Œi—Âò×»àÊ;b®.;õ¼kâ^ºîĞcÏŒû„7~ò/1#ƒóÀ3ïu¦<85õA…e…å^–³‰kns¶Ma  €.H÷¬¥úb  €¦),÷Çnû¯]¥üƒôÁ/ÿ«3OG=ê¤bÆAèJù¢òÒû3Ò.…åŞ¾µÛ1[—=ğìûq/]vÚù×Æ½Âµ·?32xûr||Ï†ÅÆ[ì÷Õ´®\7…e…å^ºVX^}-…e   øRºg-ÕÓ0   4Ma¹ÖZoÓxıº&e´UW_7fmÛÁGŸóÂò+MŒÛ6ßÅ|´Oa¹·÷>ÿb¡EùºªúsöÒe7İ;)îeö;ôÄ˜‘n¸äºûŠÙfŸ#¾w]7ï|Ä=5í€#N‰yÚ¦°¬°Ü‹Âr¦°  @¤{ÖR}1  @Ó–ûkş×°KRîAÛsÿ£cÖ¶í{È	1ß ,¶ÄÒ1cÛÖÛh«˜ö),³.•iGëÀ#N{éªIoı:îcV\e˜‘îxúµŸën´e|ÿºî©W?‹{jÒ¡Ç³´MaYa¹—Î–×Ş0æl›Â2   ]îYKõÅ4   MSXî¯y/^Ã.I¹íü+nYÛ¶û~GÆ|ƒ0Ï¼óÇŒm«¾>å£}
Ë£sÈ1İ(üÅµ·?÷ÒU³Ï>gÜÇ <ÿÎocFºå°ãÎï_—UOˆN{iÒ‘'³´MaYa¹…åLa  €.H÷¬¥úb  €¦),÷ßy—ß¯cW¤Ìƒöğ¤cÖ¶í´Û1ß Ì<ó,1cÛª2yÊGû–Goµ5×9»j9ç.^ûÑïã^ºhÕ5Ö‹û„ªà™2Ò=·ŞÿB±Ìr+Ç÷±‹ñ­ÇvqÌÒ6…e…å^ºVX^cbÎ¶),  ĞéµT_ü¸|1   ´íÕQ–Ó,_ï #N×²RŞ.˜wşcŞ6mû­=c¶¶}ğÅÆ|ƒğìëŸÇŒ´ï–SæAyì…ObÎ.Û~ç}â^º¨úpGÚÃ ,¸ğb1#İôá/ÿ«Øsÿ£ã{Ù5o¾}ÜC“N>ûÊ˜¥mUa9åkZW
Ë7Şõ\Ì7hOu¨°œòJW
ËÏığ‹˜   Ú”îYKõÅ4   MSXnÎæ[ï¯ç ¥¬]°Ã·÷yÛ´ÙÖ;Ålm{ıG¿ùÚ¶Ô2+Ä|†ÂòØœqáu1k—}ñMq/]sÎ%7Çüƒrñ5÷Äœüµª„:qÍŠµ×Û´Xo£­Šõ7Ş:¾®WŞôp1ßÅ÷³+V]}İ˜½Ig\ĞŸ[ƒ*,ï¼ÇA1OÛ–{SXÎ–  è‚tÏZª/¦a   hÚk+,7iéeWŠ×tRÎ.¸ğÊ;bŞ6m°É61[Û&½õ«˜¯m»ï{DÌÇ`(,İVÛïóvÙcÏÿ(î¥K™ôqÌ>(U7åä¯İzÿóµkwÙõÄ×¶áåş¡ÓF]|©˜»I]ù0€Â²Âr/
Ë™Â2   ]îYKõÅ4   MSXnVõ®'L˜¯ë ¤œ]ğƒ·~ó¶i­õ6‰ÙÚöè~óµmE6ê–Çî‡ŸüS1×ÜóÆÌ]5,åÛ®]×{x#æä/î{êíÚu[o£-ãkÛtâY—×ruÁ¬³Íó6©Şª\yÓ#1_Ó–{SXÎª¿¿§œmSX   Ò=k©¾˜†   i
ËÍ«R—®ë ¤Œ]±ÂJcæ¶¬2q˜«mw=úzÌ×¶—Şû»˜ÁPX2WßúXÌÜe‡{VÜK—l¼ùö1û ì²×!1'ñäË?‹×îÙ×~_ß¦;~µXf¹•c¾Azç³?Æ¼M¹ôºûb¶]vİı1_Ó–{SXÎÖßxë˜³m
Ë   tAºg-ÕÓ0   4Ma¹Çzq¼¶ƒòuÅáÇ3·eáE–ˆ¹ÚvıÏÄ|mêÂS7ùk
ËSnŸƒ‹¹»ì¦{&Å½tÅ±§^sÒËüCÌ:Œö=äøâÅ÷ş6ş³)U/]·#N87¾¾mï~ö§ÎV¿ôæOş%fmÊu·?s´í¢ïŞó5Ma¹7…ålÓ-wŒ9Û¦°  @¤{ÖR}1  @Ó–Û³ÓnÆëÛ¶”­+¾îém™nºéŠ~õ?1[›ºğ„Å3.¼.fcp–§Îr+®³wÕ"‹-Y¼ıéâ^ºà¡ç>ˆ¹é¸Ó.Y‡Íe×?ğötÖE7Ä×L‰êß§É¯×—XhÑøúA©0r¶múégˆùštûC/Ç,m»àŠÛc¾¦),÷¦°œm¹İ®1gÛ–  è‚tÏZª/¦a   hšÂr»V[sıxÛ”ruÉ†›ns·eÒ[¿Š¹ÚtöÅ7ÅlmzáßÆlÂòÔ¹ï©·cö.Ûy÷ã^ºbñ%—¹eÆg)ô¥¬Ã¢ús>Ë,³şÕ¾ÖÛh«¾=myòãNîêï=_?(g\p]ÌÙ¦yæ?fkÒÏ¼³´íœKnùš¦°Ü›Âr¶ıÎûÄœmSX   Ò=k©¾˜†   i
Ëízã'ÿ»˜kîyãunKÊÕ%_sOÌİ–»ûaÌÕ¦£Oº fkËzms1X
ËSï„3.ù»lPO:}9>f¤6ß.f[l»KÜWõA–ôú±JÇ®TçM¯¤+n|(fmËK-s5ééW?‹YÚ6¨oYPXîMa9ëÊ¿7Ãş   Æ‡tÏZª/~ü›r    Z6êÂr˜eÊÜúÀñ:·%eêšA–º/»á˜©M{xLÌÖ–3¿s}ÌÅ`½ñ“–Cæ.Úp³mãºjÂ„	ÅS¯~÷2h7ß7)f´“Ï¹2æíºs/½%îg±%–¯Ÿéø_ºûñ7âÌ ]{û“1kV]}İ˜©I]ù9_}¸#åkÚÎ{ó´íÆ»Ÿ‹ù­ú]ò¶m¤°òÊûs¶í±~ó  @›Ò=k©¾˜†   i¯*,ÄYİ¯uR®9ğˆSbö6œtö1S›¶Şq÷˜­ó-°PÌÄàu¹°|ô)ÆÌ]ôÂ»¿-fœi¦¸®Z“­ã^º`‰™íáIÅ¼]õèó?*f™u¶¸—~~¦*À§sT6Şbû83hƒ*-k·b¦Í4óÌ1O›?şì˜­i
Ë½),gûzBÌÙ¶Ÿ{?æ  €6¥{ÖR}1  @Ó{ã‹Ú=êW-¼è’q–©3¨¯ÒOYºæùw~³·aŸƒ‹™Ú´öú›Ælm8öÔ‹b&¯Ë…å#N87fîª‹¯½'î£ËºZ
?ôØ3cŞA[{ıÍbŞ.úà—ÿY¬<qí¸‰kmg¦ÔsÎÏó¥«¿÷xœ´‹®¾+æmÒiç_³4máE—ˆyÚT@S¶&=óÚÏ¿¶´ß6…åŞºVX>ô˜nüº÷É7c>   hSºg-ÕÓ0   4í‘|\»GıªÙç˜+Î2õ6ÜlÛxÍ›”rtÑ^ó7m‹mw‰yÚ´Ô2+ÄlM«êøÃÿSÌÄàu¹°\•…Ræ.Ûe¯Câ^ºì¶‡^Š{¤ªd˜²vÁ1§|'fîšwİ?æ¯\uË£qfJUBKçùRUœNs]ĞöSTï|ôµ˜£i«L\'æiÓ®{³5iu6ŠYAa¹·®–«ô¤œm»ıáWb>   hSºg-ÕÓ0   4í®ÇX»GM>üÕÇy¦ŞbK,¯ySR†.T)£e­gœ)fkÚGœóĞ].,tä©1s×µıówj-³ÜÊ#OãM{¤Í·Ù9æí‚®?ñò¸Ó.¹+o±}œ™Ë­¸j<×äÎ»ìÖ8Û‹-¹LÌÜ„·?ı·˜¡i›n¹cÌÓ¦¶?¼U}»EÊ1(
Ë½u­°|ÒÙWÄœm»ù¾I1   ´)İ³–ê‹i   švÓ=ß¯İ£&¯ıèã<Sïgß‹×¼))CWí¼çAqMšwşc–¶LzëW1WÓf˜aÆâ…w3Ñ/¾ÿ·ñ½ë‚êÉ§)s×UO,Nûé²İ÷="îe®¹í‰˜µ–\zùâùw~sÚÅ×Ü3©‰²õhb[=…9ÍvÁ)ç^3÷Û
+¯Ïß†êéÆ)S›VZuÍ˜­	ßùî1Ã ),÷ÖµÂò…Wuãß¡ëîx*æë‡»{=®  ÀW¥{ÖR}1  @Ó®¼éáÚ=jòìë¿ˆóôÇ¥×İ¯{Òù»ê±~÷Ğ´wş§˜§Õ×I§LM;úäbºãÙ~ß».Øû ccæapÄ	çÆ=uÙÅ×Ş÷2HÕÓéSÖ.XmõŠ·~ö¯1÷ \~ãƒ1ë—vİç°87µª§6§ó}ÕÇŸçí½_ü{1ûsÅÌıtøñgÇó·áÈÏ‹™Ú4÷<óÅlıöàsï3Î4˜o•èå»·>óšÂrvİOÇœm»òæGb¾©uç£¯¯‰ÿ   &÷ÕûÕÿ[}1  @Óª¯ıN÷©_U=8ÍÓ?m•æÒ¹»ìĞcÎŒûhÒÓ¯}³´á‚+o™šT=Mó£_ÿOÌCw<öâ'ñıë‚¦Ê•mYsİã¾ºª*lNzë—q/ƒrÉµ÷Æ¬]±Á&ÛtæÛ"®şŞã1ã—¦Ÿ~†âoÿ:ÎN­ívÚ+3¹å¾ÄcÚf[ïóöÓ½O½Ïİ†6?ÄÖË?ùç˜¯_^ÿä÷ÅÒË­Ï=h]üPH¥+¨>’òÊ}O¿s¶íŒ¯‹ù¦VõÄó/Ï±â*kß³[¿ÿ  è–ÉïU'S_LÃ   Ğ´=ö;²všT_×œæé¯±‰¦T:o—½ÿÅ/ºDÜKS.»á˜¥yjÌÔ¤ó.ÿ^ÌB·t¥“l½Ãn1ó°èÊS+ÇbÃÍ¶{¤V³vÅ²+¬:òäş”½-Õ“[S¶ÉzŞÕq¶vß÷ˆxÎd©eVy¢q:Î xÖå1o¿,²Ø’ñ¼m¹ÿ™wb®¶İtï¤˜¯_6Ùb‡xŞ.¨>Ğ™2ÚÃ“>ŠyÛ¶üJ«Å|ƒ2é­_Åœm;ğˆSb¾©Q• ¿zé¦›.¾   *_½ü¿ÕÓ0   4mÕÕ×­İ£&UÁ%ÍÓ“?E«	éœ]wáUwÆ½4å€ÃO9ÚĞÆ“#'·úÚÄtÏí¿ßÃ.Xo£-cæarÎ¥7Ç½uÙÇŸ÷2(mÿ¬sÍ=ïHi8åoÚIg_3M®é?KcıPLŸ~ûÃ¯Ä¬ı²×ÇÄó¶åíOÿ-æjÛ‘'óõÃŞÏÙ§_pmÌ=h<ónÌÛ¶Å—Z6æ”~ùŸ1gÛ¶Şq÷˜oJ=<éÃxcO½(¾   *é^²T_üQùb   hÓÇ¿şŸ‘¯O÷©_µÜŠ«ÆcĞßã‹b–Yg‹ïC?¤sƒõ7Ş*î§	ë¬¿YÌĞ†%—^>fjÊCÏ½sĞ=WŞôp|»`åÕÖŠ™‡MO¹ï·ën2îePÖ\g£˜³k9úô˜¿){îTÌñUÿƒ8ß/UÑ-·—K®½'kP}ı1g¿Üõèkñ¼mª²µi±%–Ù¦ÖaÇÏ×%'œqiÌ>h÷>ñFÌÛ¶\$æ¤…Y<fmS•!e›Ré÷iuïşá/ÿ3¾   *_½—ü¿ÕÓ0   4éÁgß«İŸöòúş1‡ş»ñ®gã{Ğé|ÃàÑ|\L˜0}ÜS¿Í6û1CÓzù§1OSN=÷»1İtføZğ®hªØÖ¶w>ı·bŞùŒ{ìªùX¸xñ½¿‰û„;y5æì¢ÕÖX¯øŞ/Ä}ôKõAƒÑ–O=å;ñıtÆü™cÎ¹‹¼õ«x¼AxÿÿsöCõ„ëtÎ¶í´Û1_Û.¼ò˜oJCY¹R=½>å´®ü|­~&¤|ƒ´Å¶ßYÛvßSoÅ|cµÛ>‡Åã·ñ{  €á–î'KõÅ4   M:ï²[k÷§½\çÓñ84£*“¦÷aj¥s‹s/½%î©	Wï±˜¡I§wuÌÒ„Í·Ù9f »º\ôšyæYbæaTı®K{ì²Í¶úVÜË ìsĞ±1gWí¾ïÅcÏÿ(îeJUOéİxóíâù’õ7Ù:§ß.ºú®xşo²â*k/¼ó›xÌ¶}ğÅÄŒıpõ­Æs¶íâkîùÚV½ï)ß”–²re÷}{´ïŞòhÌ;ïşï1ã œ}ñ1gÛ¶ÿö>1ßXrôñØ3Î4SñÖOÿ%Î   À—Ò=e©¾˜†   Iël°yíş´—ª€”CsöØïÈø^Ltaòí=Šûê·íwŞ;¿Il²MÌÒosÎ5O§–Éèì´ûñıìŠñôïÔ‡Ÿ÷Øe‡vÜË ¼ı³-Xp‘˜³ËªrTVyó'ÿ÷õM^ùğŠSÎ½ªXqåÕãñ¿ÎL3Ï\<ıê§ñ˜ıvímOÄ£±ÜŠ«“Şüe<n›}íç1ßÔZ~¥Õâùá•ş>f„½<&f­÷~ş§b»öŠÇîªêi½i/ƒVı|Iy¡ús˜2Ê?ş§˜s¦æ0‡}z<fåğãÎŠ3   0¹tOYª/¦a   hÊÃ“>¬İ›~“é¦›n¤Œ“GsÆZ,ÿ&éÃ¤*¾,µìŠqoı4ıô3ï|ö‡˜¡	Ïığó˜£	7ÜõLÌ@·M\kƒø~vÅm¾s«•V]3î³Ë.¿á¸—A¸ò¦‡cÆa±æº×«'Wgzşí_¯ü»‘ßU™´*?ğÌ;Å7>Xì}Ğ±c.)OîüË¿¯anè¥˜a´–^n¥O<ó²˜mj{éÍñ|ƒRèSÎA8ğˆSbÆorç#¯ÁÓ1“ê^#­·­ú}›ö3hû|\Ì;ÕÏ’”qª¦¬m«şùı¿‹{éõ¡ÌêŞçã_ÿOœ  €É¥ûÊR}1  @Sö;ô„Ú½éhyÂ¹ñx4çÍŸşïbÁ…ïÇ”Hç6·?ürÜ[¿tÖåñüM¨Ê@)C¿~Áµñüt_õUàé=íŠs/½%æVÿƒ¸Ï.›e–Y‹Ç_˜ò';ö[¯'EògSûôÚ±zè¹÷c±Xb©åŠ§^şi<~9úŒ˜kj¬¹ÎFñ\ƒtÍ÷Ye«íw-î}â˜õ«ªFm¸éØ¾5¢zÒx?ÿ¾;5[bé¸¯AëR‰½ú0GÊ8H×İñTÌ:«­±ŞÈ=\ÊùU÷>ùæ7~HÊ‡  ­t_Yª/¦a   hBõt¦9çš§vo:sÍ=oñÁÿKsîyüø~L‰tüaô«îŒûë§Ùç˜kä+¦Óùû©zbç¬³Í3ôSõd¾t~º¯*³¤÷´KºúúS£K_¿?Z«®¾nñÑ¯ş;îgÖßx«˜“ÿG±ñÛÇkÖ¤êéÈ)ËX-ºøR#JçhÚ*×‰™¦Fõ$àt®A«®sÊ;HUqù¬ï\_ÜtÏ÷G4^}H¢*SVkyê}À„	Ó¯ÿè‹^,şóAxÿÿß“AjãFFëˆãÏ‰m£Í¶y¡*à_|Íİ1g¥úó³û¾GÄÙÉí¹ÿQq   ’toYª/¦a   hBõÕåéŞt´N>ûŠx\šuá•wÄ÷c¬Ò±‡Õ1§\÷ØO~r<w?µñÒM¶Ø!›áĞÖ¸§FU8{÷çŒù‡Ùf[}+î·Ë¶Ùq¸—Aøş_Œ|Ø)åœ–-½ÜJ#ÍtÍšôêGÿ+æ™RU©?§)ß½åÑ˜cjT/Nçê‚6şÓO½ò³‘ı.²Ø’ñŸÂ]½^{?é•ÿ!æ”‰k®sZ?"ßoë¬¿ÙÈıÄŞ3r]}ÓÑ
+OŒ¯ıªE]räƒi¯   ¤ûËR}1  @¿]qãƒµ{Ò)qÎ%7ÅãÓ¬~|z:î0Û}ßÃã>ûé–û&Ås÷CõTÇtÎ~ª¸úúÇ¿‹ç§ûºX¾ù:ÇvqÜÃ0«Jbm<½ß8ì¤¸ŸA¸êæGbÆiÕ|,T<2é£x­šV=56eš;ì²ok…º~?avŞù,ŞüÉ?ÇsuE—ªÛ„»}íÿìuñ%—‰¯„¶ËøßäºÛŸŒ9©úı˜²ZlË•7=÷   _'İ_–ê‹i   úéÅ÷ş¦˜{ùj÷¤SêÒkïç¡Y[n»K|?F+sØm¼ùvq¯ıR}¥óïü&{j5ñÕö“[{ıM‹7~üOñÜ‡Õ×û×ëÊL3Ï<òçiÃ¬_öiÛ©ç~7îgúõ-ÃnÎ¹æ)îúíxÚ’rM­ê‰ÑM?‘v×½ç×ßùt<W—\ı½Çböñàê[ı«½.µÌ
ñuƒ°õ»ıU¶AÛ÷ãcÎAºäÚ{bÖAûèWÿ5òa½”y˜œxæeq   ĞKºÇ,ÕÓ0   ôKõô²ªL’îI§FõÔÄt>šµìò«Ä÷c4ÒñÆƒmvÜ=î·_Ößdë‘'S¦sO©­¶ß5«_6ÚlÛ¾g¦]‡sf|o»lÙV-^~ÿïâ~†Ùû÷Ûu]*-Ÿ{é-1ã´b–Ygû«§ÉJÊÖ/»ìuHñı7¾ˆçÕÓÛÓù¦Æ^ÏÕE;î²_ÜÃ°šqÆ™Šï~®¶Ïe–[9¾~Pîxø•ZÆAxæµÏb¾A[iÕ5cŞ.xêåŸ|@$å{xLÜ   |“tŸYª/¦a   è‡ê+Ó½h¿ì²çÁÅ?ù}<7ÍxüÅOŠé§Ÿ!¾ß$o¼hâé‹“[mõúVÄjº|T=‰;—áğñ¯ÿ§Ø÷àãâ{;–_iµâ‰—~÷6ÌºôôÏ±8öÔ‹â~á¬ï\3wóÎ·@ãO ­”¯ß9úŒâ­ŸşK<ÿXUÅ½t©qÀa'ÅsuÕ+ü}±Ô²+Æ½›9æœ»¸ıá—ã>—[qÕ83(Õ7x¤œmúèWÿ]l³ã1_TOÏO¹»àºÛŸŒ™»n³­¾÷   £‘î5KõÅ4   Sãæ{'µö¸ç[`¡â²ëî‹9hÆ”~Ex:ÖxrÀá'Ç}÷Ë‚/Ÿ
8ZU‰³zZs:v¿TÅítn†Ã-÷Mjä‰øm›{ùŠË®¿?îqXU…Ó´×aP}=~WzÚyWÇŒãÕK-W<öüâµ„™f9æì·¹æ·8äèÓ‹y7æø&×ßùt±Îú›ÅcO-·Û5¯ë™ôÑHÙ7íiXÌ¿ÀÂÅ}O½÷WYq•5âÜ m¿óŞ1k&½ùËbâZÄ\]rî¥7Çü]pÊ9WÆÌ]µúÚôíÃ   L›Òıf©¾ø£ßş?  `ª½õ³-N¿àšbù•&Öî=ÛP•ìN=ï»ÅÛŸş[ÌGpæeñ}è%g¼9éì+âŞûiç=*^|ÿoãù“Şımqä‰ç&LˆÇë‡gš©¸àŠÛâùé¾{Ÿ|³Øû cã{;Ì¶Ün—âÚÛŸ,>®şÓ°ïasÌ)Æ}‹mwÚ«¸½*.‡½µ©úàÇ<óÎ3'«¯½a1é­_Åk0(óÎ¿`ÌÚ¤êï¥ÕŸ›îy®x³*ã…\•GŸÿ¸8ë¢Šu6Ü<gjí±ÿ‘ñ¼Ãâæû&Å}ƒê[*¯¼öõ¥U&®gm•Õ×)n¸ëÙ˜¹)'uy1×<óÅ<]´Ë^‡ŒüùM{´ó¯ø^ÌÜ5Õ?üÕÅ=   Àh¥{ÎR}1  À7yòåŸ—^wßÈSe×Ù`³bÆgªİsBõô¾ªĞY•'ŸzõÓ˜şøöÇ÷àë¤cŒG÷?óÎH9&]ƒ~Úxóí‹s.¹¹xxÒG#¥äwñ§âÏşP¼úñïŠÇ_ü¤8ó;×[m¿[œí§µ×ß´xâ¥ŸÄkA7=ûú/ŠËnx Øÿ°“ŠWY=¾¯ãÉ,³Î6RÆ©>PPf‡ùƒ-ën¸EÜã0Ysİ‹ãO¿äÏO¿{lÃso|1’#åª¿¥}ZõÄç”·MÕÏƒÅ—\¦XcFş]ıÿM?=øˆÎ×cØ\tõ]q]¶Û>‡Å½|ÕÄ5×ó]QıïŒ¯+yıç1ÿÔºûñGtşÈ·y¤óƒê‰Ô]swñü;¿‰{”ïŞúhÌÛûzbÌ   c•î;KõÅ4   ½ğ£bßCNyjÕ6ßÚsä?”¯¹ÎÆÅ²+¬ZÌ>ÇœµûË®Zh‘Å‹6İfä?bïsğqÅÑ'_0òtàê‰³‡sÆˆ´F§*ü¤ë¤ùñ¬ú:útÆ“jiïFõóìÀ#O)"WOL®~~WOµ]“­‹•V]³Xd±%‹Yg›=¾—Óšê	»K-»âÈÏ°Í¶ŞiäwÄ»îWì´ûÅîû>òaœªØ®ó }ÿÍ/Šé¦›.îiÍ·ÀB#O6İx‹íGşM{nRUfL¹†YUO{í‚•W[+fÏª¥k1¬®¹í‰búégˆ{íšªà›öŒåï³ƒVı._o£-Gîkª‚nÚO/Õ“Ä?şœ‘Ÿnºm1×ÜóÆó³åV\­ØdËŠv; 8ìØ³Šëï|&^‹¶ÜşĞË#ï[Ê:H'Ÿ}eÌ   S"İ{–ê‹i   ª¯N÷‘ãÑûŸÿG¼|³Şû›bîyF÷Õúi~¼»ñçŠ‰km¯Ç0Ûrû]GìœöÌàT%Üô~1eª‚QºÎƒvÁ•·Ç¼ãAÚoÓª¿ï¬²ú:1Ï0©>œğôkŸÅ=vÅX>Ã3ÆõaQıî¿÷É7ãuvw<òj§ßTĞ©2¦ì_g­õ6‰ÇêºêÛÒ~ziã[@º¦únºmzçç,öØÿÈ˜¯mU™û‘|s  À”J÷ ¥úâ'å‹  à«nœ†
Ë|şñ0:U©/]×¯J³ÓŠËox Xq•Õãu&Õ×çßpç3qŞh?<ÀèÜşĞKñ:wÁ·vÛ?fvi¯m9óÂë†òI£Õï–ëîx*î©k¶ØöÛq_šw¾Š=ö;bäéæÕë«ÿ{Å[ï°{1aÂ„8ÓEGŸt~mïãÍã/ü¨“OÌ®¾A%åı&c)ÓwÉVÛï÷ÓKõw¹t¬ñ¬zuºƒPı¼^b©ecÎ¦-¸Ğ¢ÅEWßs  ÀÔJ÷¢¥úb   …eÆâœKn×vrinZsÁ·K.½|¼>]¶â*kŒ¼ÇiOtÇÜóÌß?¦L—Ë|ñÅ‚/s³´×6½ù“.<âäbÖÙfùºdæ™g)N>ûŠ¸®Úy÷ã^*Ûí´WñÆÿ)ÎUŞùôßŠ‹¾{g±ÉÛÇù.Øh³íFªœòW‡sF¼m«
Ç=÷~Ì8l²u<n×m¹İ.q?½¬2qíx¬ñ¬K…åÊ‡¿üÏâ”s®*–Zf…˜·ß_r™â¸Ó.*Şıùc   è‡tOZª/¦a   ¨¢šî#Ç#…åşØÿ°“âõıRš™V]çÓ#_ã®S—ì¸Ë¾Åm¾÷@÷(,÷W—Ë•[î›s³´ÏAxï*N¿àšb™åW9iu6*Î¾øÆ‘ruÊŞeUi0í©*›¦×ßû›â´ó¯.V_{Ãx¼¶m¶Õ·:ÿó¢IÕß–[aÕxmšV¯ùŞã1×Xl¼yw‹ğ½LIa¹‹OÆnZ×
Ë“«şımêß¿Í·Ş©¸ö¶'ây   ßÒ½i©¾˜†   ¦D¯'¦×Oë^ùàï‹“Ï¾²Xaå‰ñšBUd©Âöòû3LK®¾õ±‘ßm&Lˆ?3ÛP=Mûà£N+ıÁÇ1ã°8ü¸³k{›mö9Š~òûøúÑxò¥ŸGŸt~±úZÔİ¤é¦›®Øjû]‹»{=æšUEúêÒõê·oí¶_Ÿf]}HaX¥ı0|~ğÖ/‹K¯»¯Øsÿ£ŠåV\-ş{ßËŒ3Î4ò–ı=±¸ìúû‹g_ÿE<   4%İ¯–ê‹i   `J¼ûÙŠ%–Z¶ö¿?øß ¾Ù‹ïş¶¸ü†¦¸¨0¥ª'ˆîuÀÑÅU7?\¼öñïb6€iİû¿ø÷‘'Uî±ßÅ"‹-öÓª«¯SìwÈ	Åuw<ó£“Ï¾¢¶ÏªÄœ^;%Şøñ?WßúèÈï´¥—]±v®©õeIùÒkï-Şùôßbşüô÷~‹ÄôÓÏP^ûİŠK®¹§xë§ÿÏãIõwò{ÿaqİíO\q[qÂ—qr±ÏAÇGt~qæw®/®¸ñÁâÖûŸ/™ôa<   ´)ıo:¥úb   ˜RÕ4÷¿AL½êéË7ŞõlqÊ9W»î}ØÈSÓæg¾xmGcáE—ùêôªèP•¾÷ÀÅïü&€Şzîı‘™sÊ…ÅN»0ò„ßyæ?şüı&³Î6{±şÆ[”Ğn¹ïã¶{î¥·ÔöşîÏÿ_Û“Şüb¤ìW½G[ï°{±äÒË×Îÿufy–b•‰k;ï~`qÒY—ü>VR›ê½­®ÛAG2òíé:¥–Ya¤ \==»*¡{’0   @·¥ÿ§T_LÃ    S£zâ×:lV¬¹ÎÆ#…Ÿe—_%¾±{ùı¿+ñ“âÁgßù*úêÉjÕ8¯¼é¡‘§$ßtÏs#ëUq¼ú:èªø\=ù:€şªşúâ{S<ıê§#?§oèå‘ŸÑÕWıßpç3#?Ÿ{şã‘¯ÿ¯üñ¯ÿ'g¼ª®ÏäÒkšT•Ÿû×Å£?øxä½¹úÖÇF
±·=øâÈûUıŞômÍ¨ş.R•Èï{ê­‘"ùy—İZœsÉM#Åÿêï.÷>ùfñÔ+?k´Ä   @3R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾øãòÅ           c‘ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.Õü7å           À¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/ş¤|1          ÀX¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_üÉß–           cºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾ø“¿ı          ŒIê&—ê‹i           —ÔM.ÕZ¾          `,R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹?ı»r           `R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}ñgå‹          Æ"u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅOË          ŒEê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/~ú÷å           À¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.Õ?+_          0©›\ª/~ö÷ÿo          €1IİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾øÙ?”           cºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}ñçå‹          Æ"u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}ñçÿ«           ƒÔM.ÕQ¾          `,R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_ü¼|1          ÀX¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úâç¿ûÿ           ŒIê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}ñ‹òÅ           c‘ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ/ş±           ƒÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ_–/          ‹ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úâ/_           ŒAê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/şª|1          ÀX¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹¿úıÿ          `LR7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úâ¯Ë          ŒEê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/şúŸÊ          €1HİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹¿)_          0©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾ø›.           Æ u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾øÛòÅ           c‘ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}ñ·ÿüÿ          “ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹¿ıßå           À¤nr©¾ø7å‹          Æ"u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\ª/¦a          €^R7¹T_LÃ           ½¤nr©¾˜†          zIİäR}1          ô’ºÉ¥úb          è%u“KõÅ4          ĞKê&—ê‹i           —ÔM.ÕÓ0          @/©›\úÿ³s/Pû}uAàÿ³R@t‰­.ÚÒ
tR°¼¤¨™:Ú„¦ÕL«©ïNšS©£c“ÌJI-5'ÒœÌi°QF'TÀ¥$&ˆ„€Š "‚×.3Í;{ÿüæ<ç|Ÿ}öŞçœçò¾Ÿ³Ögéÿï÷ìÛÙïûç«óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞø¢ôa          €Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆıt
           hÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7şxú0          @‹¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
    (ùöï|æÍ×Ã7İóeOşÛ7Ÿ÷ÙŸ=“Û³ü™(ÇTş\Ÿï}ë•×tØCÇöO6ìŸÚ=  Ñï˜¹Íß  °¨69™7FÁ  \¶üruü²õ”†¢‘qáÀ/z£{·Ê}Œr—e«3-ÊÍõ±._4ß­ò:G¹¹|Ci´®5ı~–óFŸ«ıİ.Šmå÷ÈıÔ¬õ1Ãß Q^ nŸègA¯(?—+ZÃV~Ÿ»›j~×ôû$  l/ªMNæQ0  —m«b®­mY@åoå‹	¸[iQn®ıpù¢ùnåâë´åï ãÿƒ·–Öšßï¢¸V~ÜŞšBå)ëp7D?zEù¹\Ñ¶òûÂİ³õß  @½¨69™7FÁ  \—ü26ä­ÑØsÈ}Y[Œ4Œ+Ê_Ã‹g¸.kÏ±('×Ë~¸lkF¯ıÓZúâ??«ù3ÃçkŸİ!&ú·’šıã÷ÈË²´'ÆkZ»nÖàîÈ?×³Úß1"Q^.›ßçhíƒ’áï  `½¨69™7FÁ  \¯5/ò³èEmÎ¹æK¡qB¯(ï_LÀuê=o¢\\?ûá²åß¢ù_²ÅïœFÍG¿?FŸ›ÊÏwÏï®­¿ãE9–ø=r;KçxtÔıŠL îÚŸSQ.®G´¦Kü>w·ôü]‘E¹  €vQmr2oŒ‚ ¸n½/h³¥/ı{s¯ı’àÅ,Àeè-€Œrqı|éxùzŠFòºF¹¸,5çñ±ß·¢Ïn)ºç1~<Ÿšó!úû£vÍ¬ÀİãïÅ»Éïs,©ù½3å  ÚEµÉÉ¼1
 àºõ~y“E‘—Àk¾(è“/&àzEÏô’(·C´ŞK¢<ü¦áËşÚŸùKzŠrL”‹ËQû»W´–k~­5½g‰ß#Ï£öï…(¶å\‰â}ägÎùÈ%˜ş,¨åázø}îvÚòoÓÚß=§¢\  @»¨69™7şø+ €[(z[ãŞâ _äËüåa’åªå+Ys/à¼¢gzI”‡Û!Zï%Q~åà‹ş–Ÿù%İËA..G´n‘h-{
JZMï¹$ÊQâ÷ÈuZö@}î˜­Î2 løûßùÈ%˜ş,¨åáºDëZâ¼ºl[ÿmÚówiå  ÚEµÉÉ¼1
 àúE/`k´¼ n)DëıÂ ÊUâ‹	¸^CAD‹(·C´ŞK¢<>[[ùõü> `ù²µ|Ù­eo±@‹é=—D9Jü¹N4§Çìlkü³ŞùÈ%ğ÷âİ­k‰óê²mı·iÏß ö  l'ªMNæQ0  ×/z	[£õqÏ—DYO¡R”§ÄKg¸^=gK”‡Û!Zï%Q»nú®‚e‰Öì˜h-[ŠzöOÏŞò”ø=²_ëšF9¢ÏE¬œ†çKãïÅ»)Z×çÕåÚëoÓqÎş. €íDµÉÉ¼1
 àúE/ak´¾ n)H™Šò•D9J|1×ËĞŒEë½$Ês—E„[})ÜZœ˜ùbørµ®ç±µŒ>;5ü®ÖrÏŞßï¢\%~ì×ú3<ÊQû7ÆVçpÜô™v>r	ZÖdQ®K´®%Î«Ëıî¿Õït§ø»  ˆEµÉÉ¼1
 àúE/bkô¼ òÔh-VŠr”xñ×ËĞŒEë½$Ês—Est/…
–/Wm¡èàØZ.í‹éïi5ûhÍ¾‰ò•ø=²_4Ÿ%Qléw­Î0à¸èg‚ó‘KàïÅ»)Z×çÕeŠÖjËßëjş®°7  `{Qmr2oŒ‚ ¸~ÑËØ=/ˆ{¾(ÊZ_G9J¼|†ëÕs®Dy¸¢õ^å¹«=O[})\ó…ğ”‚åËÕzş–Ö2ïi¾üß¥½—óMcrÛÚı:ÎW#÷!ÊCYÏyåìµ€eÇgç#—`ú³¡F”‡ë­k‰óêò{v÷øİ.ÿÎ8¾G¾wæ÷H  ØGT›œÌ_œ> Àí3~!Û"¿´ò•L_ ·ˆòÅ—ä—ĞQàòû«$ÊÃí­÷’(Ï]Tz–z~æG4•äß¢\œ_ëù{-kõ½Äï‘}zÎƒ(p^¥gÙùÈ%ğ÷âİ­k‰óê²”Û­ş6  Î'ªMNæQ0  ×/zù[£çñš‚å–"—(¾Äp½|ÍX´ŞK¢<wÍRáàV_
/İ'¢`ùrµ¿
–ë9¢<Ày•~8¹ş^¼›¢u-q^]Sım
  œOT›œÌ£`  ®_ôò·FÏbËÀÖ|ÍX´ŞK¢<wIMÑàV_
×ÜkJÁòåj=,3ÖsDy€óYú9à|äø{ñnŠÖµÄyuNù·)  p>Qmr2oŒ‚ ¸~ÑËß=/ˆ{
-_Dñ%¾˜€ëåhÆ¢õ^å¹+j.oõ¥pÏï
–/Wëù«`™±ó ÊœGÍÏ ç#—Àß‹wS´®%Î«ó;õß¦  ÀùDµÉÉ¼1
 àúE/kô¼ î)L´|yÅ—øb®—/ ‹Ö{I”ç.hù™¼Õ—Â=¿(X¾\­ç¯‚eÆzÎƒ(pzµç¿ó‘KàïÅ»)Z×çÕyãoS  à|¢ÚädŞ pı¢—¿5z^÷&ZŠ\¢ø_LÀõò4cÑz/‰òÜvùgj4Çlõ¥pÏï
–/Wëù«`™±ó ÊœN~n[Î~ç#—Àß‹wS´®%Î«ó9×ß¦  ÀùDµÉÉ¼1
 àúE/kô¼ î)L(X"¾€f,Zï%QÛ¬ç™ÙêKáß,_®Ö½¤`™±ó ÊœFÏ3ë|äø{ñnŠÖµÄyu=ÏçV›  çÕ&'óÆ( €ë½ü­Ñó‚¸õÿkÆXËı¢ø’=¿˜ÈcÎò=¢ñÃ¿Ÿë…{¾oVêc6üÛ9ú:î_©_Ql6ÄíÙï|ÿÒÿ~ê¹;•ÒeÃ¼äÏDñk»gI”glX¯iÜ0†s®ã0×Y4ö¡Y*Cÿ­Ï^sİkI”ç6Ês}l=–lµN9O”¿di/ÏD4¶Ü¾¿·¡Ç‡áß·~N¡u?{-jE}/ÉóåÙÂ°¿‡=İ{ø÷(ş’õœQ=ŒçüØ>æşŸİ±aœã±]Ò¸Jë°Õœâ½ò}Çı‹ú8íç}Í¹£û.ÉıŠò]£KY‹½c,+ŠÛÛxîõ-ŠË¢˜%QKT3/ÇÖ,·Ÿ‹ş}ù¹C?õ1‹â×˜ŞkÉó0ƒh†?÷™Qêc–Û·êgÎqì>KÎ=O{¨ûKØ'kã<6¾Û¸¶  Ä¢ÚädŞ pı¦/	kõ¼DÌ/£\5¢|ÇDñ%ùÅh”g±æ~œâål¾ÇšµäQş­ôÎá¸_ãá{Ìí%¯óŞÖì£é:­1^ãZQÖñœzó½zÆºÕ<×È÷êéc–c·˜Ï(÷’(O$ŠíåŞ[Û¨/[[ZÃ~ÛÃ­ûí”ÏB–ïõ£äÔçJ­ŞçºÅô=óéY÷(OIŸ(ÏyôÌû©÷ySœ?[{Í^Ï}Øêùò·ZêKÍXK9¶xFKÏNëZäÏ·Îÿ)îÑ#Ïmkß¦¶ìëÚ¾ÔÈ÷ˆî}n—´QîVÇÎÊÜ¿–qn5¦­ó±g§9.IëzòÜë–ÿçĞsMï±Öš>FùzDùKö˜‡<è^%¹§z¾²ÜÇµÊÆ{ªÅxÿíiè[Ï:D¶ÜŸÙš¹Ï¶ìO”¿Õ±ş´§~  8¨69™7FÁ  \¿èå`—‡k^ÄFù‰âKr¿¢<=¶x¾õKğAëKâZ[¿HÎù¢û¬µe?/y÷¶å>Úâ‹ˆ¾Ls¬Ùs[î«Èó½Å<—lñ<Ö>QÎ%QHÛ#Ê½§-×gÉÒ>ëyÖ¢=ÑûLä¸i®­m1ßkŸƒ­­=ƒjLï¹Õ¾í™Ë(OÉ–ûêÎüV=Ï}«µÏÌ–{|‹ùò¶:Ö‡Ögk«<‘èÙY›·fîOq[ÌéT4Çµò8·|6JÖôs—¶Y”³U×4ïš5ŞëYÈÖ¬Át®{Æ8¿$kÖë˜µ{sl‹s#Ço±·¢Ü%[ÎÃgHô¼niËs®¥¯[ŞwÉ°¶ºçVk²õÏ×-úåmõcÍ8÷ü ÀùEµÉÉ¼1
 àúE/kô¼8ŒòÔh½W”£d‹/&–^8OïQóù­^Î–î•_(çŸŞ«õ…ş–/î£üƒ¡¿ãÏ×öu:Æ9Ç¥®ó)ä¾FãŒÇÒ²‡ÖÌAi~Y?µÕşŸZšÃé¼-}~~.íïás¹¯µs=kå[å9&#Ïcë¾É1köy¯¨/{Zcş÷(®d¼o{â#ã>m%÷­öy¨ıüÒ|Jë~ïİ7[š§%=ç^”§dº¶½r_£üƒé~Xú|ÏØ÷ûõoK½c­é[^ßişš¸µóŸï±´Æ%Óı’õ<K¥q¬íãôÙYó¬Ecœâ­Íá°÷¦÷Êÿİ2üÙÖş®Y×¹Q?Ní×b,Ç®Y›;ÎÕÒ÷cÖŒç˜š~k2ÄLç&ÿ{K¾©!ö’”ÆÍGíçj¥½9İ+KŸ§G”³d‹yXšóé=j>?·-Ôö±veã¸c¢¸=Mç.ÿ÷Ò¾+Y»'³šùœŞ§fò¿¯İ+[ÎOMŸk¬  —+ªMNæQ0  ×/z!X£õ¥aş|”§F”¯$ÊQRób½dil¥üK/p×¾œ-õm)wş÷–Ì{ÏãR—úºöË…¥şsO¡gKs6Ö»>KıŠ¬‰=fë5ìİ3KsóFq­zïÓ2çQü’(Ï’(Oš±l½/ZE}ÚÓÒx—öMd8zbÙê9ô>ÙÒ>:÷ÊZÛ^Ñ}Çz×¿çgK”§d‹ıÔ»Öì½Sé]»=ë\Ó¯¥ù‹bÆ¶˜ÿŞù›î™¥=V2ÎsLOşñü¬é_d:şüß{ß£Ç±>çæ˜¼ï£ØcZúÛš{­šñîíR×"’ã£¼K†³²7ş˜iÿzÕö«´&ãÃç­mÉ4ï¹•Æ°´Gkæ5Šk±4Ç¹QÜRß–ÆVå+Ys¯lÍXzç¯Uo—ú7Å¢Ïï©wßÓóûæXÍ<–Öº¦ß¥øZkç§7ş˜iÿ  ¸¢ÚädŞ pı¢—5Z_‚¶¼àëyÙå)9öR¾ÆÒ‹Ø¥Ü5/r£¸¥ÜµóZÓ¿±Ú¼‘Ò©ıb Ôß5_.,ÍÃ9×ù–ßÒº/ÅõìŸ–üƒ|Ÿ¸%Qÿz¬™ïl)~ímKû9ß?ŠË–bÇJy‰ò,‰òÔŠòe¹ï=ûùòúG}.Ùj,-ë?Èıíéó’µÏÁ`ÍóÕÌIw	–Îš©gz¬õ~YÏ:GyJö×Òó·¿Õ^ßZÏy°4µjî]³®[åYÒsçji,÷å˜õæfmÿ"ãyïé[µk{lÜ-y[ænm³(oÉ÷<…k\‹–ûòY²Çó°ÕÏ™(wd|¾EÆcÌóÔ3WQŞsYZ³(&Rš‡èóµ–æwi½–â{÷W”«$÷#ÊSci–r/ÅgQ\‹5}¬éß u×şÓ+Ê»dÍY·´Ï³šqm•gIÍ}¦.ıg  —%ªMNæQ0  ×/zX£åhïËKÿb¢f\5cÈŸ‰b=ı[ê[sLë‹ê(Ç’¥şn±ßz÷Ó%¯ó),­ÿR¿—Æ=å(iİŸ‘ÜÇñËÿİ“·f,©¹o7Åµ<Sc5ÏÃRî(æ˜Ö9r,‰òÔ86KÏÄ¹µ>“Yï~™ªÙ?Kòü÷EÎÙ3¦lÜ·5ã©ÙÃKı¿Ô=ÕzN®GÏ:×ÌÿT”§dÍ¸jæ0Š›ŠâÆ¶z†·Ôsl1ÚûÖîš}¹vï¯™«Ø©i"½÷™>ù¿ÇsŸóö<ûYúµõ=†<-JÏş¸oK¢ñ•´äD9KÖîıS¸ÖµèÙ³ÑXÇıÈc(ÍGÉ¸o=jï[»§zæg,Êy.¥¹©Á±\ÑgkÔ¬[7ÅõœµQ’Ö¹Ô<û5ÏûÒíí_VÓÇ¥9b©ï çYíÙS=÷m×XÍsR»¾5k™­£ù‰Æ9³Ü§š¹ˆŒû ÀíÕ&'óÆ( €ë½¬Ñòò3Š_²æe|”¯¤÷^Q®±-_8·¾lrŒµ¼ho}Qİó?Ê3ÖšséEy‹i©s®óŞ¶èsM±-ÖºE©ÿÑçKzÏ’AÍ³V»G–æ¥·¯[ä]Ê1å8&Š_å©­Wï¼RÍ>›ªİwKZÏƒ©Òüc\QÎ±ÚıP3/[­Á–ZŸåµÏGÏç˜(WI”§¤w\5ã©]÷¥µX;÷{è9¶x¢¼‘(6R;5kĞ;W=q‘¨O‘(¶EéyíKô\ä<QşìTgLÍ}¢¸c¢ø’(G­(_É%=c×¼=ûuêØó°ô3%Rz¶–´Œ¥å>Q|­(ß9l½G¥Ñg—Ôô­v½–ö\ÏYå)é=¯¢\cµy­ÍXïs¶Åü.å˜ŠrDzÎ²Şyë¹o‰r•ÔŞ§%wíZ¬™§ù™:vÿÖ½”m±æ  \–¨69™7¾øgR   ·Nô"°Æ·?ı™a¾±ü™‘÷^–ùjE9KzîW3®{/œƒØH?ÅEj_,×¬aÖú¢ºu>s?¢<c[¬QËz.y÷V³.Y;ÅÓºÖ5ktÌÒ3Ğº÷³(O-ç;«é{ë3±ÕóĞ:¯µgUÅ/‰ò,‰Ö«uïKÏ¾nYƒ’Ú}©™ß(®¤õÛêyDñSQÜ9µ¿kŸ‘ıÓ³ÆQ’qÕ%ŠÔ<×köûzÖsíYT»g[×4Êéíï\ÿ;iØ-ùZæ"Š¯U37Q\‹=îÑóüGy¦öœ÷šy8&ÊWÒ3?§õyêR×¢å9”îİ“{ÍÏ˜(_¤u?­™£(ß9ä1GıkİGÑ¼DŸ+©Û(6²Çï1Q’óªf}ZúÅOEq%[õ±fÆj÷ekŞ¬uÏGj÷ğXël¹GLËœEñ5zæg¬´F§˜{  ._T›œÌ£`  ®_ô"°FéåcÖûrs‹—QŞ’Ö/&j_G±ÇDñSµsÅÅOõ¬e”ç˜Úù\ÚsSÓ~ß¶uŞÛ–_¾E±%Qcjú©ÙO={¿uŸ¢\S-{cëı[;QìTmß-ãâ—Dy–L÷]ëùrN­óŸõîë©g*«ßÖó wİ¶~¾²(~ªåY8…SÍ÷ gÿôÌY”§¤g\Q©–¾ï±'÷Ö³kÎ¢–ûµî›–g!Š_Ò{v¢=Z›³e.¢øµëZ»Ï#µ÷h=×²(Ï1-c¨ûÖ>·îï±(_IÏùx*×¾kÎ…šûFq%½kİ2g=÷ˆòÔˆrCÔ·©y™Î{í9ÇÓ²¿kŸÇ(ö˜(¾¤uÏÕç–y­='¢Ø©ÚñjûÙš7kİ¯‘3´eî³–ó-Š?¦¥ï­}ôÌÏ æQ\IÏ9 Àe‹j““yãKÒ‡ ¸}¢5şÕÓŸ9Ë•Û¾¡ães–_>F9{DùKò½£<ÇÔ¾tbÙ*gÃ(î˜š±·æÌ¢<ÇÔ=Ëû+ÊqÌ8wÍXÇ¶Z“±=rî¡vÍk×£e³(Ç1­¹³Ú~÷ìıÖ=šm=ßÙ¹Ö0Šj×–qGñK¢<%Ó¹Èÿ}îRåù÷¿F^³(W«g*‹rEzÎƒ(Ï’-Ÿ‡Á9÷Ö:ßkŸ•ıÓs&GyJZÇµõùœí‘so=ë™c¢\5Zökë<µœ«=kĞ3WcQÎ¬&owL¿¤e>Zæy¬eßôÜ#ÊsLë¹YÓ÷Öœ={på+Y{îïéÚ×¢÷\¨]“Ö±dQ%Qczæ«w¢\§ÖÒ÷¼^5{t0Í½&ö˜–õÚ#g_RûljŸ‘(ö˜­sn™¯vµk•?Å—´ì×cZÇ“µì¿–qµî½,ÊsL¿¤g~²Ú±ÔîÍ±(  ×+ªMNæQ0  ×/z	xj[¼l‹îQÒòr¸ö¥sëçÚ¼KsÕó²?Ê3Öó¢ºeM[_T·Ìí¸ï·i÷V»&µılİ—QcZ÷Oå9&Š/Écò”Dy"­û"Ê1U³‡ÇÏQIËóÅÓ2§Qü’(Ï1ÓıÖ2æKÑú<f[Iµ{i¬eı{Æå)©½GëŞ¨Í»ÕZlaú<,Yû¼ì½Q’ÖqE9"­kå˜º¤3«g={÷ëÙĞzŸÖüQ’¹,=¥Ü­óå(iİ=óĞz´®eÖ2OQ|IMÿOyGùJ.éÌ™Šú[rikÑ{.Ôî×Ö±dQ’Öç­õLD¹–DyN­g[æh¼Æ-qãû•´®W”cªå™‰âKZr×îİÖg¼6oÍÜÖîŸ½æ´æÌÌZÏ¬uoEjçg¬vLÙøùZÒºO²(Ï1-ıôÌOV»6-ó3ˆò  p½¢ÚädŞ pı¢—€§_˜nñ’9İ¯¤ååpíKÕÖÎµ/é—^4×æ‹òŒõ¼¨nYÛÕ9¦õExËšÔö©%g¶Õ:ï©e½[Ö9Š´½v­Æ¢<ÇDñ%­ıßk¾³(G$ŠÛãyh9«¢øc¢ø%QÈtZŸÿKÑ2÷ƒÖ½wLË~´<S=c‹ò”ìñ<dµ}o™½ÕÎÅ`í3³÷şDyJZÆugş©ô¬gëœZ÷j”£¤õìiGÏ\¢|S9ÿxy¾zæz|ß­gBÏ<´­k™µÌU_RÓÿÖıİ:ïcQ¾’5÷Ú[Ôß’K[‹ç¡å~­cÉ¢<%­÷ˆrÔ8ÅXöĞ³ÆYí¹7Î_{µô©6ç Ê‰b#QlÉÏGKÎ¬ögĞ–çQKkû—Eñ‘–œƒÖ½iÙËƒg«Fë>É¢<%Q’ùiGíş‹ò  p½¢ÚädŞ pı¢—€kåƒüB7Ë/;·x©\#êSIËKÕ(>Rû"{Ğò’>Š´äDyÆz^T·¬uÏ‹êAÍ<sr›ÖyO-}Üc[Ï‰ıå9&Š/¹¤=åˆ,õ9Š‰´<cµçJKÎ,Ê±$Ê35İg­ıº$-{nĞú\S»îc-ÏTÏØ¢<%QHK¿³–¾GñçĞzş®}nöŞ?ƒ(OIË¸ö\ç(G¤gNöĞ³½gQ”«$ÊQÒ:–Ö5è™«ìÔ?«¢>”´öïg@Ë3:hÙ—Q|IMÿOyGùJN½[Dı-¹´µØûyhKå)‰r”D9jœb,{è=û³<æš³iø|í9ÖrFFñ%QHí>bKZÇ(>ÒòÌe[ÎoiwílÉÙ2æAí~-éy¾j×³õÌiİ'Yë=ZçlÏùÉZûŸEy  ¸^Qmr2oŒ‚ ¸~ÑKÀ[¼ ŞKÔß’Ú—é-/Ò[_8·äâÇ¢˜cjúÙó¢ºeô¼¨ËñK÷ËŸ¹më¼—¨/Ç´Kkİs®ôìŸ(Ï1Q|Ií>D9‰âK¢‘Ò>nÙ³­c_:[ZóeQ%Q±éëé×%iYÓAÏ³YZóHË9»÷ØZò·ô;kÉÅŸÃôÙX²öÙÙ{ÿ¢<%-ãŠâ‰âK¢‘9ÙCÏzöœE-ÏÖ ÊSÒ3–(Ï1=ù³S¯uÔ‡’Ö3ág@Ï~iÙ—­ù£S§<‹£|%kÏı=]ûZìı<´%‹òÓó¬Eyjì=–½ôıcKk>Ü£ö›æ/‰âK¢‘Ú}Å–Ô>-{·å™ËZrGñƒ–<­çĞÒ¾lÍ×Ò×Aí~-éy¾j×3Š-iİ'Yë¹Öz=ç'»Ös €íDµÉÉ¼1
 àúE/klñ‚x/QKj_¨·¼Pm}œ?å‰,Í}Ë‹åšuìyQİ²?ZÆ^Ò:çÇ\Ë:ï¡u­{ú˜ç`:Ç¹­w¼{ÑÅ—Ô'Yë|G9J¢‘RŸ[æ·eìƒ<Óç"çi}¶ã<µ¢<Y´>½ıº$Óù®Ñû|Nµîù¬eÎ÷[ËóĞºWZú¾Õz¬Õ2YÏ1¶÷şDyJjÇÕÚÿ(GI”#²v¶Ò³={¿uŸfQ’½ÇÒ“?‹rí)êCIë^ì™‡=ÏâAë¾¬İ“µóÓºÇ[ç},ÊW²æ^§pÍk±÷óĞ:–,ÊsÌ)çjï±ì)ê[«<şÖs*Òºç¢%QHí^ˆbKjó¶ì§=•Ö´¥µãË÷ö5çioÖ2æÁ9ösV;¾(¶¤gŞZÖ8k]ç=ç'kíå àzEµÉÉ¼1
 àúE/klñ‚x/QKj_ÜF±Ç´¾pnyI_3÷5/—k×°çEuëşˆrôh}	‰òsîuŞZKÿ²sôqjï/:¢ø’–=Ø:ßQ’(G¤ÔçèóÇlñü­õkI”§tî]Â¾_£uße[¹4¯Ç´œ³{-Š?¦¥ßYKß/e¶¿kÏˆ½÷Ï ÊSR;®Öıå(‰rD.á¬ÎzÖ³gïGy–DyJöŞ›=ùÏ±ÎQ?JZûxŠ3 õ9Ízöe{”kĞ27K¹¦Öì(_Éš{Êµ®ÅŞÏCëX²(Ï1Q|Éš¹Ú{,{ê9“é9«ÆZûå(‰rDj÷B[²GŞ=•Ö3úü1k­-ôìñµ{9ÛëíOë>Éö>×®ıg  —/ªMNæQ0  ×/z	Xc‹Ä{‰ú[²×{iy	œ?;}œÛjÖ/¦ç%rÖº?zïsÌšıå;‡/-ÖÊ÷ŒúrÌ%œ{ÑÅ—Ô'Yë|ï)ê_}ö˜–±ï%ê×’i¼¯£Ï.aœkôì»­õ¥¹´œ…{-Š?‡–9ÙSëù»öÙÙ{ÿ¢<%µãêÙŸ{‰úwj=ëÙsEyJzöéŞ{³'ÿ9~VEı(iíã)Î€ç´g_f9nz¿<'µ}îék¶foDùJÎ±{\ãZìı<ä¾E9JZ…(¾dÍ\õŒ%Ês=ë\²f{÷ù¢şMEq%µsÅCéy>Ìš=±…}ÕrÖ³×ÚsŞ´œÍƒû´ÌÛ^ó3Ø»ÿ  \¾¨69™7FÁ  \¿è%`K~Qõ·¤æ}ÏËÚ½ô¼Ìn±ÅX[÷ÇóÛóÅË]ZçHë—­ë¼‡/:¢<ÇDñ%-û®§ï{‰ú×ú<ô<s[‹úµdŸŸ»è3S—0Ö^µcÛêYï9c[ÎÂ=ÇÖÓ÷½œãçC¤õ[ûÜì½Q’ÚqµÎ×¢şZÏzÖ>¯ƒ{ôìÓ½ïsª½¿VÔ’Ö¹>Å<ìùsd+kÏ’=>ˆò•¬¹×58çZìı<ôŒ­öYèéûš¹êK”ç\Öî³HÏ¹µG?zEı›ŠâJjöXÏŞİË±ç¹µk­-œëçnÏZÖœ¡=ÏÉ©Ç–ûì5?ƒSÍ  —+ªMNæQ0  ×/z	Xã’_Fı-¹-_L¬µå{öGÏKZ¿|¹ë\õ£ä}œêù¢#ÊsL_Ò²çzú¾—¨=Ïd”ç”¢>-b[×ã’–ô¬ëVcí9c[Î™=Çv×>DZŸ™–ó1²÷şDyJjÇuégş©õ¬gíó:è9zöiïùåŠœjï¯õ£¤u®O1={¦u_öÚêéÙãƒ(_Éš{]²KX‹½Ÿ‡1Ö>=}_3W=c‰òœSÔÇµZÏ®­öı¢şMEq%5{¬gïîåØóÜós,Ês*çú¹Û³–5gègçØŞ÷Ùk~§š'  .WT›œÌ£`  ®_ô°Æ%¿(Œú[RóÅDëKôšœ—¢ô":»çEuïşèyi½¤e-nó:×ˆÆXÒò…Ä^zöL”ç˜(¾¤eODñ%§>w[Ÿ‡,ÊsJQŸ–ä¸}t­ÏÏºnµ÷z~´œ3{­5÷mûùi}nÖÎÉŞûgå)©W[rê3ÿÔzÖ³uNzÎ„}Ú3–,Ê9ÕŞ_+êGIë\ŸbzöLë¾luì¬ÍíùŞÇşı˜=>ˆò•¬¹×%º¤µØûyhKVû,œêlôŒ%ÊsN=ë]£vÍ²(¾¤%÷¢>•Ôì±Ö½»fßöêy¾¢<§r®Ÿ»=ÏTÍÅ-éOÏ¹VÓÿÁ^ó3èéÿ¹Ï  ¶Õ&'óÆ( €ë½¬qÉ/
£ş–Ü–/&Z•^@çşkÜó¢zÍşèyq½¤v=nã:·ˆÆXr	ãïÙ/Qc¢ø’–9‰âKN}î¶>Ù©û8õiIÏ\ãpÎuíùyÒòÅçc»ë?"­ÏÎÚ9Ù{ÿ¢<%µãŠbKÎ}î­g=[ç¤çLèÙ§=cÉ¢\‘SÌÕ¢~”´ÎuÏ<´{şiU:cÇã*}.Ò³ÇQ¾’5÷º$—¸{?­cÉjŸ…çlÍ\õŒ%Êsn=k^£vİ¢Ø’½ÎÆZQŸJjöXëŞ]³o{õ<_ç\«sõ·çyª9C£¸%=ãé9×jú?Øk~=ı?ç> `{Qmr2oŒ‚ ¸~ÑKÀ—ü¢0êoÉ_LdQK×®ôrxú’¹çEõÚıÑóòzÉ][çÑøJjæto={%ÊsL_Ò2'Q|IË@[èyÎı³!êÓŞ.ùçaäœëÚóó¤eßï9¶ÜQÛ¤õü]û3cïı3ˆò”Ô+Š-ééû5éYÏÖ³¨ç¹íÙ§=cÉ¢\‘SÌÕ¢~”´ÎuÏ<´>G={fë¹.õ!ÏÙô~¹-úì1={|å+Ys¯KpÉk±÷óĞ:–¬öYèyÎÖÌUÏX¢<— gİkÔ¬]WÒ²ßöõ©¤fõìİ(ÏzúXûìîá\ıíy–jöt·¤g<=çZË3¹×üöü Àuˆj““yc Àõ‹^Ö¸ä…QKnËKòš-½ÖµçEõû£ç¾K–^ ß†u^#ß’(Ï)õ|Ñå9&Š/©9OQ|IË@[èyNİÇ©¨OµzÆ;ˆúr©zÆ¹Å™õœë-{jÏ±õäòÜ&­çoËùÙ{ÿ¢<%µãŠbKzú~MzÖ³õ,ê¹GÏ>İû>§˜«-Dı(iëyh}öü9²déŞÇæ+·GŸ?¦g¢|%kîuN×°{?­cÉjŸ…çlÍ\õŒ%ÊsIzÆ´$ºÏXSÒ²ßöõ©¤fõìİ(ÏzúxÎµêéï?w÷:C£¸%=ãé9ZÖùšÆ  p¢ÚädŞ pı¢—€5.ùEaÔß’š/&z^Ö^ÒÕôÿXÏ9öœgë/ã¢ûÎ9ÖKoÉ9¿ÜÊzöG”ç˜(¾¤æ<´ö½%÷z¾<Ì¢\§õ§Æğ÷œÙ©×fuİêœë™ß–3fÏ±õô}«y»T§>ÃöŞ?ƒ(OIí¸N=_—îÏTÏ=zæ}ïûœb®¶õ£¤u®{æ¡õ8×ÏÈ¥ó¡4WK±S­ó>å+Ys¯s¹–µØûyhKVû,ôô=‹rÕèK”çÒôœW%Kûã”û{QŸJjúÛ³w·øÑ¢w_D¹Ná\?w{Ö²æİóìÛû>{ÍÏàTó ÀåŠj““yc Àõ‹^Ö¸ä…QKöúb¢åeíj^—Ö³gì[ï—ÙÇ”Öåš×y=ó\óüì©§ÏQc¢ø’–ùØ»ïkõ<ÙÖÏ‹¨?K¦9zÖ%»–³ ÷3êÉVkº÷»çØöîû5j}VZÎÇÈ©Ö ÊSR;®K?óO­g=kŸ×AÏ=zöéŞ÷9Å\m!êGIë\÷ÌCë°çÏ‘HÍ˜–æ©õléÙãƒ(_Éš{Úµ­ÅŞÏCëX²Úg¡§ïY”«FÏX¢<—¨w.‰î1¸¶yŒúSRó<îıÜm¡wOÔ>¿[;õÏİÁ^kyªñìı<î5?ƒ=Æ  p¢ÚädŞ pı¢—€5.ùEaÔß’½¾˜È¢\§RÛç¥ñ÷Œ}ıÑ»S¥ñ^ã:o©çKƒ,Êu*=}òÅ—Ôœ'ƒ¾ŸòÏŞçá”}œŠú³$ÊÓû,ìqömíT_¢FzöTË~Úsl½ÏC”ë¶h}NZÎÇÈŞûgå)©×¥Ÿù§Ö³µÏëX”§¤gŸöœ=-k{ª¹Z+êGIë\ŸâèYËŞ¹®=–ò·-={|å+Ys¯SºÆµØûyhKVû,ôô=‹rÕèK”ç’õŒ1RZÃ{´Á[ŠúSRó<zïöèíã¹Öê”?wÇzæ©fN5(Ï’(Ï1{ÍÏ`ÏŸ1  \‡¨69™7¾4} €Û'z	X#¿(Œò]‚¨¿%ûÉ_æ™Šb—œkZ^.Gñc=/ª÷w^¯è-¢¼ƒèóK.ùyhÑóåJ–ã¢|§Ğ³¢<ÇDñ%µçIÖ3ß-ù·õ¡F”ë¢¾,‰òôœ{Ù©×§Gï—¨Q®V=óÚr¾ì=¶(~Émùùi=×>{ïŸA”§¤v\×pæŸRÏzö<O§Ø§=kÛ²7O5WkEı(iëSœ=k¹ç¾¬éÿ)öø ÊW²æ^§r­k±÷óĞ:–¬åYˆâ—ôk=c‰ò\º=1UÚ#=çãš=¾VÔŸ’Ú¾F±Kz÷n¯¨5¢\{ëÙW[ÌgÏóRs†öäíO”§¤õYÜk~{ÿŒ àòEµÉÉ¼1
 àúE/k\ò‹Â¨¿%µ/n{^¨¶¼°İRÔ—È%¼pò·ÌUOŸÆ¢œƒkZç­õÎkë—[êY¯(Ï1Q|IË\ôÎw”k/=ó›µœ[Šú²$Ê“åç:úü’K?zÆµÕzöìù–ùÜ{lwùçC¤u>ZÎÇÈŞûgå)©WOÿ³(×mĞ3-Ïë õ\èÙ§=gO”ç˜SÍÕZQ?JZçºgZÏ€µlë–{DñS§<‹£|%kîu
×¼{?­cÉZ…(~Iï¹Ö3–(Ï9}şí˜ñJ{¤gÏeQ®SˆúRRû<öÌoË³·…Ş=ĞûŒ­Ñr¶ègÏ~®]Ç(¶¤g<Q’Ö=¸çüd={ôû €ıDµÉÉ¼1
 àúE/k\ò‹Â¨¿%µ_Lô¼HÏ¶š«œ'÷u)_Ë‹ßš¾åÏD±%-cò×®ÃXÏKî,Ê5¸–uŞK4¶­_€l¥gDy‰âKZ÷q”cIÏ³rLÎUZ»Şç!‹òí-êÇ’(Ï geçz~kô¬é–ç[”¿¤ålÙ{l½ÏÃ–ó—÷äVùÖj}>Ö]{ïŸA”§¤e\Qü’µó6–sõÌÉzÖ³gï÷Ü'ÊSÒz6´®é©æj­¨%§˜‡ÖıŞsÎ·ÌuËjç'.Š?¦uŞÇ¢|%kîµ·k_‹½Ÿ‡Ö±d-ÏBÏ³Öú<zÆå9‡¡ï-s›õìli£˜%köùTÎU»¢¾”Ôö³gïf­kxLÎ“ûZÊ×ÛÇ,Ê·§¾n1—=ÏHíŞkSmŞ±(OIë=öœŸl8ÛZl±î  \¨69™7FÁ  \¿è%`K~Qõ·¤ö‹‰¶Ymş%ÃİÒKàÖ>F9¦zÆİ²?Æù[â=/º£<ƒñf§\ç=õÌç Ê··­×*Š/iİ½óİó¬L÷zôïYïócGıXåœû<ØC^—¨Ï%[ì·¬g>[öÑŞc;÷~Î‹sı|˜j=¿ÖÎÃŞûgå)iWëœ¶xÇóıû©õ¬gï<D¹J¢%­gOë¾<å\­õ£¤õLè™‡Ö¹n]Ë¬e®[Î€Ú¾·+­ó>å+Ys¯½]ûZìı<´%kyNñ<zÆå9‡¡ï={e9î™Ë¬eo3OôïSãû×¨ãyÍzÖ02¬Ai­zû˜õ>g½òı¢~”l½ŸjÕÎMkî9ò”D9JöœŸ¬ç,Ùbİ ¸Qmr2o|é«R   ·Nô°Æ½…A¾Kõ·äŞAH_cí|/s—úÚúÒ7Ê1Õó¢ºu¼C\ËZŒµö1Ê1ÅÔ8Õ:ï©g½[ô;ß¿e[÷|å9&Š/iƒoøÇí_Ğekçz¼ÎKóİ3Çƒ–µ<&ÏQÔ‰ú°$Ê3v®5ÚKÏx¶XÇl¼ïjµ¬ÿ)Æå¨±v‡çğ’öUëÙ°¶ï{ïŸA”§¤e\ç:OÆs·v/n¥g={ûŞz¯(GIë³å(9å\­õ£¤u_÷ÌCë°çÏ‘Öş×ö½uÿ­9O¢|%kÏ®½Ü†µØûyhKVû,¢K¢<KzÆå9‡ñ™Ô:¿ƒ–ñ/í‘qZ¬ÙëÙx¿×ÎÃøş5ZúÅ×è]ÃÁ°–5}mY÷©µıÌjÏ›=µEÿÆ{ªVí˜²(ş˜–¼Ykß[óg{ÏOÏşÜbİ ¸Qmr2oŒ‚ ¸~ÑKÀ—ü¢0êoIË=/Ó½s6¾çRñıjD9¦z^T·uMì ö…wÍz_ú:ïí\_nåØ!Oôï‘¾Fy‰âKZÎ“A”§FÏ½-9ÆëÒjM³a}k¿üšŞ¿F”gjèG«s?Ë‘óm«qôì¥ÚµÏN1¶{zçq|ÏKÚS­ÏÅÚó`ïı3ˆò”´+ÊQcÍüm‘ck=ë¹fÿGùi½OË³Ğ³§«^Q?JZç¢gZÏ€=´æ®í{ËşËÖœQ¾’5÷ÚÓmX‹½Ÿ‡Ö±d­çNë:d=g[ÏX¢<ç0£ŞıÒ²Wjæ7Š«±f¿÷äß»FKî½;èÙÃÙø59ZÖ}jÍZe9>ç©9sNuLõÌOÍx-ãjïÖ¾G9–ì=?Ãi±Åº p9¢ÚädŞ pı¢—€5.ùEaÔß’–—Ã=/mÇZçm|¿¥~öô­æ…ò)^$c[_ÖóS3æK^çS8õø³ñ>Ûû‹(Ï1Q|IÏúõŒa°ö~µkuê>N÷`ô™È8¦V”'ÅÖˆrSË—§ƒ¥}’×¸f/M×µFËy°ÇØ¦zÆ0¶æ~=ÏÒZÏ…µıß{ÿ¢<%­ãj·±9ß¯uÿí©g=×ô¿å~­÷‰rÅ/9õ\õŠúQÒºŸ{æ¡õØóçHÏ³åë™“sdå+YºWîÿšşôºk±÷óĞ3GµÏÂ`ï1d=ãÈ¢\ç0=“ZçxP{¶E±S½sšõìùñıZÆ?¾o–¾õìİ±5ÏJK?O½VÓy‰>3µÇÏİÜ÷¥ÏLûZ£õü‰rÅÓ2g­}ì=?={siM ¸.Qmr2oŒ‚ ¸~ÑKÀ—ü¢0êoIë‹øêcµs7}}f¬ç…rv¬?¹½ç%r¶öE~ïKõšşÖÎÿ¥®ó©ô®ı ÇGy§¦û¬6nĞÓÏ(Ï1Q|IkÿQ®Z-s=kyÎ¦±=zŸ½Ú¸lW+ÊéƒŞ=±—q”öÊ°^5ãÜúŞSÓ½S£ezî3V{Ïéù}æœ¦ı[²öYØ{ÿ¢<%=ãŠòÔª½ßt¾zæbO=ëÙó¼ÕîÙÖ¹ŠrDz×àsÕ#êGIë³Ó3­s?å)©ëÚı7Vš£¾¢|5zÆåŸ9õ~½k±÷óĞ3G=ëØ3wQHÏ¢|çÍOô¹%5sQÚãSQ|­ÚûLûÜ²³ql–ñgkû¬öy™>‹ÑgYójû9Ú¸>–öÂĞ¥õÜú¾‘–{DñÇÔî½Ö==¶÷üL÷uÚ= Àuˆj““yc Àõ‹^Öh}Q{JQ—DyJ¢­òN_¸æÿ^>×ÌwÏåÁ8ÿ±>´¨éïX”£çeôÒKï-úÕ*ßsËu>¥iÿzÿ±9˜~vÉ4¾F”'’ûÅ—ô~)ÍE«|ïéæ:z6ÆŸ«±U‡>sçÿÎù§ılÏql­(Ï1Ñ<ÖhÇŞ¢>–ë^·á3Ó5ôì¡é.éYŸš~G¢\­òØ=ÑgÇŸ»­ó½ö9Ø{ÿ¢<K¢<%=c™Êó9_Ş?Y´6ãÏ]‚ÜÏi—¬}ZîÅGj×rÍş?Ç\õˆú±$ÊsÌ)æ¡v=ÇjïÑzfr\û'ß¯7×`Ü¯=ó3îûØ0†5ÏF¯»ºµ{5ë×±µ^å*©Ù3=çÅXïX¶­sÍø§jæ#Š;¦gÿMåqL÷dîgí¿ñçjLãkDyJ¢­òL÷[şïh§óUc«µÖfœ{èçt½Z÷è8¶Æ±ü¹?Ãg¦}ê™—ùörd©¿c{äœÚ{~jÇ0¶f<  \¨69™7FÁ  \·ñËÜV=/jO¡wLQ®’5s×ªöeÿ^}êÍ›ã†~å1”öLŸ9j-½ônÍ×;öµë|J§Öº>Y”gI”'Ò3ş5ë¸´·Ô3×Ù)ûØ3—Q%Q’(GŞ9ßCÏ:×#eœ£v­nÓŸ=çC¯Úù=µ¨¯K¢<µzöOëÜõ®k”kIÏ~íuIçÏ g®·øû£ö¾µsV³kŸá½ÏÎ-œâÙ9ÅĞó\ÖÎõÏ|Î¹v^r|şïš=ß»ÎãÜÓşÖ>k[ºkÑs¯åŠôÌQïZöì«Ò½zòM×åœ­skÿ–æ¤g¼={¤WëŞêİQ®’-öZ­5{ò”kuª½4¾O^‡qš>ì}†ÕŒ¯¦ÏYíëíë`ïù©™“©Ös  €ËÕ&'óÆ( €ëÖórPû2õÔzÇÔóâó_N´Îs”ca^¢kUšãÒËêÚ—ŞKëÑòò|ì×ù”N1ş¬´?éí[í^è=O¢\µz¾¸iÕ3×c§ècÏ3±÷~ôî‹líÜoeÍ"µãêİ;Q®H»¤uıÇz÷\‹KıùĞ;ö(W­(_(×1½ÏFï³}gş^zæz«ç¡æŞ5÷ª}Ö®AÏ\­9ÛzœâÙÙûgHÅ/©ëŞ9:fØ£[ämy¶¢ø^[=Ó­nÃZäÏEñK¢\‘(vIí³éO¾ßøÉÿ=Í“ÿ½w®r\_{†¯‘ïõ-úÅM-ÍAS£wn[ôÌï³Øs¯Òm%Ïstï§X«Ş~nqvÕ¬cï|D¹–Ôî‘­ú}ÎıåŠD±KÖüŒ àòDµÉÉ¼1
 àzmñR}‹— [‹úY£w,½/qkôôi«şLï½6ïÒX–ò/Åg¥5ñ%Kı[cmßNaë/¦j¾˜‰ô®KÍœ¯9#×|‘²ÅÙ\Ò;×c9Ç%>{î‡±µk´Åla«5¬Ïšy«¹ÇšñDùjm5‘Ö½yJ½ãîÓšıÓr&Gñ5Î1®µÏç©í}Ô¨ùİ¦´®5cÈñ[ô7Ê]ã”ë¿÷™pŠ=Ó;†¬ækÆ05¾ßykç(Ûjµk¿‡k_‹5÷©É¿æYˆòÕZsßÈ0Ö­òÖÌİÖjÖz©_K9ÖŒ«¦kôö-ÊU£÷\Úzïõöi*Ïå%÷s«¾Õì™5û¶wOÖŞ³”¿f¶Ø/{ÏÏšµò p¢ÚädŞøÒWı*  Wî_=ıY÷¾œßêEğ çÌ¹£{JMÑÁ’</y­cÉŸò­ñ›/™ãû-‰òµˆÆ¿fŒ5c©İ“ùs¿YxôÿÇ.õmÍ\]Ú:ŸÚãÏ¢ıV’?ŸÕî™c¢½4äßê<‰ò×Z;¾Hë\/Ù£9gt¯c†ıåjQ»^[ì½lÍŞØÒÚ±Ô¬×VÏTÎíáÜ¾ÕšDùkl±§jæöÔò8·k[í3ï·ÅşYºçV÷æ(ºGÉ{xª§{æfË¹îÓ"çˆòE÷ªCÇôÈ÷Y»7r-æê˜-ú¸´[ì™Ò<l1†¬f®×åØ¾ZÓÿı‘c¢\-zî»¥k\‹üïkûÛ«¹ıTÏÂ1[Ü?ßœ§’ïõ%ÍûÒØ·ÓVs<ÖÓ·-<–|ïÖûçÏGùÖÈ}‰îµÆkµU?×öm©×rÆMó×ö»·Oƒ¿Å3tl~jÇ±äX~  ®KT›œÌ£`  .Û/{-½(ŞJtï-´ö?¿,İâÅk–×-ºG­Ü—(ï’ÜÿÒKßıT;[ÍİÔÚ¹œº¤u>‡s¿gÖö}ï³³¤ô\•l5Ş¥g{söq¯ıploFŸ]+;º×)õ>ÓÇæi,Š[k<g½}/é]“¼·êOÍÜRÔÇ­Mï¹×óç7ú÷-ôì£­Æœï½×™ßcÏµ¬yf¶|v³œkM²(ïV¢ûõØk]Ç{wËuß#ú÷-ù#½ãZÚWQÌ’5gEr.¹¤3êšÖ"ŠY+È¿Çó6Îß¢wo¦sÙ3¶“åµ>ç~];Çä±m=®­~.ôô-Ê³…Ü—è~Çä~ç˜(W«¥sfs®Õ’Şù+Í×Vãô¬ÓÖÏuë>Dy×÷k«çbl‹q p>Qmr2oŒ‚ ¸l{¾”]rª‡Ñ½·ĞÛÿ5_Pô¼è>¦¥ù¾µ_2´ì©–/.Æ}ÇåûõÌgË˜z´ÌïÔ–ë|.çËşk1ì•ü?£_kí^ìw^£µ÷®u>îµíÑè³kåñG÷:µ–¹lY³(~­ñœåÿ=úÌk×$ÏMo¿í½s‹úºµé=÷z¾³ñ<Gÿ¾…5û¨wìù§:ó[ì¹–ƒ-<wkúšç‹~dQş­D÷ë±×º÷pïYºd|èß·0ä?¦vşZöUËÏŸü¹(G«–{f[=#[º–µˆb×ß»ekÕ-’ç°vmÇîW[ş·,ßk|6\ŠÜ§q_Çí­ó“cßck=ıÊrßz× Ê·…ñœ·ÈãÈ±QÎ%{¯ÏØ9ÖªFK¿júÒ;ÎkÖkm¿z÷g$Ê¿Ö¸½ÏCÉ–ã àô¢ÚädŞ  ±á¤ìØ‹ÙÜ>|&Ê±…c÷îİû%Ã±¼½9s\/ÅæÏ÷=vïŞû÷º”u>§aìKã?åºÜVÃ>*Í÷ğ™sÍ÷ğL”öÃ°'¢xÎëØÚkæ9®ççk÷Ï±=4|Æ³¹­¥ŸeÙ°.ùsQX2ì³­÷V)ï^ûu¯±œÊmZ‹ÛèØ<ÖÎå°Ãgs¾ès—jgôoÙ0?YşÜ±9Zš§­ûõ+>smkRkim²s­ÏØĞÏR‡~Fñ{9Ö§¡/·eßÔî“aÜY”  ®IT›œÌ£`          €’¨69™7şDú0          @‹¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7şÄ«S      ›ØBó  <yIDAT     @ƒ¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ2}           ET›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ—ıì¯          4‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞø²Ÿı5          €&Qmr2oŒ‚          J¢ÚädŞøòôa          €Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæ/ÿ¹           Ğ ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7¾"}           ET›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2o|ÅÏ§           €Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæ?•>          Ğ"ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæ?ıš_          hÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆW¦          ´ˆj““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2o|åkS           @ƒ¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞø3éÃ           -¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼ñg^—           DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2o|Uú0          @‹¨69™7¾êu¿          Ğ$ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yã«á7           šDµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞø³éÃ           -¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oüÙ×§           €Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oü¹ôa          €Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yãÏ½şß          4‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼ñç~1           4ˆj““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7ş|ú0          @‹¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'óÆ(           $ªMNæQ0          @IT›œÌ£`          €’¨69™7FÁ           %Qmr2oŒ‚          J¢ÚädŞ          ”DµÉÉ¼1
          (‰j““yc          PÕ&'ü?“†›W¿î×Â           ‘\ƒ<­KNr­ò¿>j¸ç¯~C˜           ’k§uÉI®U~à£†{^úŠ×†I           "¹yZ—œäZå^;j¸ç…/ù™›×¤           €¹yZ—œäZå^5j¸çy/|ÙÍkŞ          *ääi]r’k•xù¨áçüÈ‹Â$           ‘\ƒ<­KNr­ò?>j¸çÙÏy~˜           ’k§uÉI®U~àGG÷<ã™?&          ˆääi]r’k•ø¡QÃ=Oÿ×ß&          ˆääi]r’k•øşQÃ=Oûg„I           "¹yZ—œäZå1j¸ç[ÿ÷§‡I           "¹yZ—œäZå>j¸ç›ÿùÓÂ$           ‘\ƒ<­KNr­òO5Üóÿô[Ã$           ‘\ƒ<­KNr­òÿ|ÔpÏ×|İSÃ$           ‘\ƒ<­KNr­ò_3j¸çKÿöW„I           "¹yZ—œäZå¾dÔpÏçşÕ/“           D>çó>ÿ &ù¾\«üÀgîùäOı‹a          €H®AÖ%'Ÿ™<ğøQÃ=ÿÍŸşïÂ$           ‘\ƒ<­KNr­ò1j¸çƒ?ôÃÃ$           ‘ş;¨I¾/×*?ğ¾£†{şà{<ææµ)           F®AÖ%'(yàG÷¼Óïş=a          €H®AÖ%'¹Vù‡îyë·~›0	          @$× Oë’“\«|ïúOÉÁ?¾òç~)L          0–k§õÈI®Q~óõóÉÁ~ôß½<L          0–k§õÈI®Q~óõï’ƒ|ï³Ÿ&          ûg=÷ ù¾\£üæëYÉÁ¾íiÿg˜          `,×Oë‘“\£üæëK>ğußøÍ7¯ı¥ÿ           P”k§õÈI®Q~óõ”äàOøâ'‡É           Æríñ´9É5Êo¾>'9øÀ'~Ê§…É           Æríñ´9É5Êo¾ştrğû£“          ŒåÚãi=r’k”ß|½wrğw}·G…É           Æríñ´9É5Êo¾~[rğ·z«‡†É           Æríñ´9É5Ê×›’ƒıØK^&          ÈrÍñ´9ÉµÉ³ëÉÁ¿óÿzv˜           Ë5ÇÓ:ä$×&Ï®oO>øµ_ÿÔ0)          @–k§uÈI®M]_™|ğo|á“Â¤           Y®9Ö!'¹6yv}frğÁû„O¹y]J          É5ÇÓ:ä$×&Ï®?™|ğC>ôÃÃ¤           Y®9Ö!'¹6yv=:9øà»¼Ëï“          d¹æxZ‡œäÚäÙõ°äàƒoù …I          ²\s<­CNrmrx½.9øğüĞÂÄ          Àİ–k§õÇI®I>z}OrğußøÍar          ànËµÆÓúã$×$½¾29øìÏıëar          ànËµÆÓúã$×$½>99xÜG~t˜          ¸Ûr­ñ´ş8É5ÉG¯÷M~ïÃ&          î¶\k<­?NrMòÑë­’YĞËæÂ           wS®1j“\“\¼^’}Çw=3¼	          p7åãiİq’k‘¯‘~ùß{Jx          ànÊ5ÆÓºã$×"/^_˜~Ò§|Zx          ànÊ5ÆÓºã$×"/^*9|ìûàÍëŞ˜          $¹ÆxZwœäZäÅë]’ƒÀ‡=ìíÂ›           wS®1Ö'¹¹êú•ä øß¾à¥á          €»%×Oë“\ƒ\}ı`rà©ÿìÛÂ›          wK®-Ö'¹¹úúÚä Á_ş¬Ï½ù…”          ¸ÛşÒg~ÎA­ñ}¹¹úú¸ä ÁcßÿÃ›          wK®-Ö'¹¹úú/“ƒoñoqóš_üõğ†          ÀİkŠsmñ´Ş8É5ÈM×«“ƒ$ßñ]ßŞ          ¸rMñ´Î8ÉµÇÍ×¿H=ñ‹ŸŞ          ¸ğE_zPc|_®=n¾>+9HôüO„7          î†\S<­3Nríqóõ~ÉA¢ßñ;gxS          ànÈ5ÅÓ:ã$×w]¿‘$ûç¾ ¼1          p»åZâi}q’k»¯ïMşı§ü£ğæ          Àíö•_õj‹ïË5Çİ×ßL>şã?)¼9          p»}ìÇ}âAmñ}¹æ¸ûúÈä á»½Û£Â›          ·Û»¾ë#j‹ïË5Çİ×oMfI_ü“¯;           ÜN¹†8ª-NrÍñªëÉAÒ¯ûÇßv          ¸rñ´®8ÉµÆ«¯¯JìÇ}bØ	          àvú˜ÇÂAMñ}¹ÖxõõÑÉAâw|Çw
;          ÜN¹†xZWœäZãÕ×ƒ“ÿœ$Æ3Ÿv          ¸]ríğ´8É5Æ¹Öx“ë»’ƒ|Á¿$ì          p»äÚái=q’kŒ7»>;9¸Áş ¾y}º9          p»åÚái=q’kŒ7»~2»ÉO½êÂ          ·C®j‰“\c¼éõ²äà&ßøÔo¹yı›RG          €[)×Oëˆ“\[¼ùõ5ÉÁ>î>9ì          p;äšáiq’k‹7¿şTrp£ßı{~oØ)          àvÈ5ÃÓ:â$×o~=4™İìY?ğÃaÇ          €ë–k…£â$×ïrıëäàfO|Ò“ÃÎ          ×-×
Oë‡“\S¼Ûõ¹ÉÁÿÈhØ9          àºåZáiıp’kŠw»Ìnú‚yØA          à:åá¨v8É5Å»^/Jnú¥_ö÷R§ş#          pK<ùË¾â fø¾\K¼ûõ¤äàÆğ$ì$          prğ´n8ÉµÄ»_If7ŞıDØQ          àº<ÿ…?9«¾/×ŸäzArpó'}é—‡          ®K®Ö'¹†ød×“ƒ¼ïc? ì,          p]Şï±xP+|ß’“]0™uâ‡Ÿÿâ°Ã          ÀuÈ5ÁQ­p’kˆOz=/9èÄ¾øÉa§         €ëğÄ'ı­ƒáû~$9ùõùÉAGŞû}Ş7ì4          pŞç½ßAğ}¹vøä×£’YgşÍ¿ğæSG         €ëòCÿößÍêƒïËµÃg¹›tæğ%aç         €Ë–k§õÁI®>ÛõW“ƒ=ô¡;          \¶‡>ô­jƒïË5Ãg»Ş5™uê;¿ûûÂ           —ééßıÌY]ğ}¹fø¬×3“ƒN}ü'~j8          à2åài]p’k…Ï~}RrĞ±?ä!7?ûÚ7…          .K®ı}ÈCŞê &ø¾\+|öë-’_N:÷÷şş?¸ùÅ_N           .Z®ıÖ'¹F8×
_Äõ”ä ƒïû~ï          ¸,ï÷Ø8¨¾/×_Ìõ¾É¬“Ïxæ†          .C®ùj“\#|Q×s’ƒNş÷Ÿöá           €Ëk~§uÀI®¾¸ëÓ“ƒ¾íÛ>ìæu¿ôáÀ          €óÊµ¾{ØÛÔ ß—kƒ/îzëä?$}Ê×|}88          à¼r­ï´ş7É5Á¹6ø"¯¯K:ü‡?èƒÃÁ          çõAäCjïË5Á{}P2ëô³Ÿó#á          €óÈ5¾Qío’k‚/úz^rĞéOÿŒ¿          8\ã;­ûMr-ğÅ_Ÿ™tüAzĞÍ+^õºp           ÀiåÚŞ=øÁ5¿÷åZà‹¿–üûä ó_ğÄ/	          œV®íÖû&¹8×_ÅõÉÁ Şñßéæip          Ày½Ó;ıîƒZßûrğÕ\¿/™â~Ê?          œÆWıƒ4«ó½/× _ÕõÏ’ƒA<ú1ï          8Ç¼ç{ÔøŞ÷ÍÉÕ]8™æ[¾õiáÀ         €}åZŞ¨Æ7ùÀä*¯g$ƒù£ÿõãÂÁ          ûÊµ¼ÓúŞ$×ü^íõß&³A=ãû~ œ           `¹†7ªíMrÍïU_?–êÏìÇ‡“           ì#×ğNëz“\ë{õ×§'³Á½àE/'          ØV®İjz“\ë{+®×$ƒû+Ÿõ9ád           ÛÊµ»ÓzŞ$×øŞšë	ÉÁ ßò-ßòæ¥/{u8!          À6rÍn®İÖó&¹Æ÷Ö\¿#ùÏÉÁ ÿÒ_ùìpR          €mäšİio’k{sï­ºşN2ì_üŠpb          €ur­nTÃ›äÚŞ[w½}òŸ’ƒÁş…OÿnŞğ+iB          €MåZİiın’kzsmï­¼œÌı¼¾4œ            O®Ñjw“\Ó{k¯ßšüZr0èOşÔO'	          è“kt§u»I®åÍ5½·úú¢d6øçşÈo~)M          °N®Íjv“\Ë{ë¯·IŞ˜şñÿIád          mrmî´^7É5¼¹–÷N\ŸŸÌ&áşÍóÂ	          êäšÜ¨V7É5¼wæzpòºä`şÜŸ|8i          @?ûç>ö F÷¾\»›kxïÔõ?%³Éø¾ïÿ7áÄ          e¹7ªÑMríî»ş‹äg“ƒÉø¨ş“áä          e¹wZŸ›äšİ\»{'¯ÏJf“òMßòmá          ±\ƒÕæ&¹f÷N_?‘LÊ|G‡“          Ärî´.7Éµºwşú˜d69_ú·şN8‘          À¡\{Õä&¹V×•®ïN&çmßöa7/åÏ‡
          ü¦\s›ko§õ¸I®Ñuİ¿Ş?™MÒ_ø´Ï'          øM¹æ6ªÅMr®kt}m2›¨ï}ösÂ‰         €».×ÚF5¸I®ÍuM®wH~=9˜¬?ö¸'          îº\k;­¿MrMn®Íu×ç%³Iû_ú¿¦	ıO          À}¹Æ6ª½MrM®«p½(9˜´G>êİÃI         €»*×ØNën“\‹ëZ¸şL2›¼/zÒ“Ã‰         €»&×ÖF5·I®ÅuU\ß‘LŞoù-¿åæßşè‹Ã	         €»"×ÔæÚÚi½m’kp]•×û$³Iüãõ'nŞ˜&          îª\SÕÚ&¹×Õpııd6‘Oùê¯'          n»\KÕØ&¹öÖÕx=(yEr0™¿í·ıö›—ıÔÏ…           ·U®¡Íµ´ÓúÚ$×ÜæÚ[WÇõ§“Ù¤~ìã?áæ¿š&          îˆ\CÕÖ&¹æÖµâzj2›Øoú–o          n›\;ÕÔ&ÿ$q­¼Ş>yCr0¹ïüÎïróú7şF¸           p[äšÙ\;;­§Mrm®µump}j2›äÏøKŸ.
          Ü¹f6ª¥Mr­kÃë_%³‰ş?ş¯Ã…         €k—ke£Ú$×Öº6¾Ş-ù“ƒÉ~¯÷~Ÿpq          àÚåZÙiıl’kjsm­k‡ës“Ù¤ÿµ¿ñ…á         ÀµÊ5²Qíl’kj];^ÏNfÿíßñİáB         ÀµÉµ±QÍl’ki];_Mf“ÿûŞõİn^û‹¿.          \‹\›kc£šÙ$×ÒºNp}Q2[€ÿ„O          ®E®‰je“\Cë:áõŒd¶_óÿq¸p          pér-lT#›|Oâ:ñõÈä?$‹ñÖoı67ÏáKÃ         €K•k`s-ì´>6É5³¹vÖu†ëS“Ù¢|Øıcá"         À¥Ê5°Qml’kf]g¼şI2[˜/xÂ‡	          —&×¾F5±I®•uùz›äåÉlş]ßsó¦´€          p©rÍkT›äÙ\+ëº€ëqÉl‘õî¿ÿæoú÷áÂ         À¹åZ×wÿı`V{_®‘u]Ğõ7“ÙB}ò§ü…pq         àÜ>åS?mVÿzß—&®¼•Ììï~ÅW…          ç’k\£Ú×äÙ‰ëB¯?üçd¶pOÿ®ï	          N-×¶F5¯I®…Í5±®¾şb2[¼‡?ü7¯xåÏ‡          §’kZñˆwÕ»Ş—ka]Wp}u2[ÀÇ}Ä          Nå#>ò£fu®÷åX×]ÏJfù?~Î_          ö–kY£×$×¾º®ìz—äõÉlAÿÑ7ü“›7ıZZt          8‘\ÃÕ¶&¹æ5×¾º®ğúÉlQô İ|ÿşp¸          `k¹v5×°Fµ­ÉG'®+¾şz2[Ø÷xôcn~á—~-Ü          °•\³škW£šÖä¯%®[p}S2[à?ûç>&Ü          °•\³Õ²&¹ÆÕuK®%ÏOfı…O|R¸1          `­\«Õ°&¹¶5×¸ºnÑõÉLfşÕ_ûõá         €^¹F5ª]MrMk®muİÂëñI´è7ÿòÛ¿3Ü(          Ğ*×¦F5«÷åšV×-¾şV2[ø·{»ßzóœç>?Ü0          P+×¤æÚÔ¨f5Éµ¬®;p=5™m€G>êİo~êU¯	7          ,Éµ¨¹&5ªUMr«ë]ÏHfáCÿ«¿ùå´Y           U®EjT“\»êºc×oO^”Ì6ÄŸÿ˜Ç‡          ù˜ı¸Y]ê}?äÚU×¼Ş#ù¥d¶1>ûs>/ÜH          0•kO£šÔ$×ªæšU×¾şXm›/ÿ»_n(          äšÓ¨õ¾Ç%.×Ÿ”DäæŸ~ó?7          |Ó?ûÖ°õ¾\£êr½ùúü$Ú(7ÿòißn0          î®\cÕŞ—kS]®ÙõÕÉlÃ<øÁ¾yúwO¸Ñ          ¸{rméCòYİé}¹&Õå:z=-™mœ‡=ìín¾÷™?n8          îï{ÖŞ«-jN“\‹êr¯‡$ÏJfèíßşnóCÏ7          ·_®%}‡wø]³:Óûrj®Eu¹¯·K›Ì6Òï}øÃoÿc/7           ·W®!}ø#1«/½/×æT—«úú]É%³õÈG>êæ%?ñÓáF         àöÉµ£¹†4ª-MrÍi®=u¹š¯wN~2™m¬G?æ=o~êg~>Ü          Ü¹fô1ïù^³zÒûr­é#—«ûz÷äUÉlƒ=ö±póš×½1Ü˜          \¿\+škF£ZÒ$×˜æZS—kõõŞÉë“ÙFûı°›7¼é7Â
         ÀõÊ5¢¹V4ª!Mrmi®1u¹6»>0ùÕd¶áòF|Í/¼éæ—ıÿ         àÈµ¡…bå\SškK]®Í¯OşŸd¶ñûşpóS¯zM¸a         ¸¹&4×†F5£I®%Í5¥.×n×G'Ñæ»yÌ{¾×ÍK~ò§Ã         ÀåËµ ~Ì{†µ¢÷}Târí~ı™$Ú€7|ä£n~ôÇ^ró+iÃ         p=rh®jDïË5¤.×É®\Ÿÿ_zÏ6ãÃñˆ›ç<÷ùáF    øÿÚ»÷§¿ç3ãvfÓŒ0º¥3Z§E“î:¬RV©¶é6Øíì¶ã|ªÊ²ãXm7»J$¥e°Œ°ºfê°l1ªJ(qj¢%Aœ•lXb«ÒqèÔ¾÷º¾soÜI.q‹;wîÃã5óø¾ïëã‡ç     ô?Ù~®¿ÁKu¡]²Ìú|»„a©æºë~¬İ>íŞò         @ÿñó;îé´ŸU²ÍfÔl¥m‡ğ›°ÔtÔ¨5ÛOv[ùÃ         `åËÖ3›ÏªÙˆf+j¶Ò·u˜–ú¡>¼ıøºŸ”?p          Vl<³õ¬Ğmh6¢fıfcÂ“¡úÁ¶KtUùC          ïeÛY5Ÿ]²	Í6Ô¬ßmÃğ`¨~¸íô3Î.ğ          ôl:«Ö³K¶ ³~»…¡ú·c;¡üá         °âeËY5]²ÍÔ¬ßoÍpG¨~Èmï}ö/           +N6œUÛÙ%ÛÏl@ÍÌ†‡kCõƒnŸÿÂ.mÎsóÊÇ          @ïÉf3ÛÍªéì’Íg¶Ÿfrç…ê‡İFÓ~qßÌòa          ğÁe«™ÍfÕrvÉÖÓlÀoB¨~àmÍ5×j×^cù@          X~Ùhf«Y5œ]²ñ44;(T?ô)ü[ùP          xÿ²Í¬šÍn²í4tÛ=¼ª};ñ¤SÊ         @Ïe“Yµš]²åÜ-˜Úm	Õh_ıÚ>íåùËÇ         À»Ë3[ÌªÑì2;dËi6è÷‘05T¡m±Å–í®{Y>$          ––íe6˜U›Ù%ÛÍl8Í†Ô.	ÕƒhÃ†k?¼èâòA         ğl.³½¬šÌ.ÙlšÙª‡Ñqìq'´ñ          XÚqÇ«l0»ÉVÓlÈoßğV¨Iû¥/·9s_,         ÀP”me6–U{Ù%ÛÌl4Í¬k[…Bõ`ÚnØnºù¶¶àõxd          CX6•ÙVVÍe—l2³Í4³%6,\ª‡ÓqæÙç–         `(È–²j,»É3›L3[Æ¾ªÔqÈ¡‡µß-x³|„          ƒQ¶“ÙPVme7ß	fÖÃ¿	ÕcjcÆ|ªóÏ™W         `0Éf2ÛÉª©ì’ÍåÁÌŞç6w„êauœøİ‰åÃ         ²•¬Ên²µÌæÒÌ>ÀÎÕëØu×İÚC³Ÿ()         À@”md6’U;ÙM6–fÖKŞÕck9²]páEåƒ         H²‰Ì6²j&»dS™m¥™õòş,LÕÃë8àÀƒÛËó_+/         @–d¶U#ÙÍ!›J3[›ªØ±éf›µn¼¹|È          ıQ¶Ù@Vmd7“‚™õÑvO‡ê1vüã„Ë         ĞŸdóXµİd3™í¤™õñF†‹Cõ0;¶ŞúÓíÆ›¦–         `eÊÆ1[Çªì&[Él&Íl%îĞğf¨iÇ‘ÿpT›ÿêåc         èKÙ4fÛX5İd™¤™õ“}2ÜªÛ±ÑÆ·Ë¯¼º|ø          }![Æl«Ö±›©!ÛH3ë‡;)Tw‘}÷? =;÷…ò          ¬Ù.fÃXµKÈÒÌúù¶ÓBõˆ;Ö^û#mÊù?,         @oÊf1ÛÅªiì&ÛÇl Íl íøğ¿¡zÔãÆíÙ|è±öZ         €Ş”b¶ŠUÃØM¶Ù<šÙ İfáúP=ğUW]µ2éÔòP          ,l³Q¬ÚÅn²qÌÖÑÌÁ¿Õcï=zL»ô²+Ê£         ĞÙ"f“XµŠİdÓ˜m£™²­.ÕÃ_dìØ¿jwŞ=½<"          •l³A¬ÚÄ%\²i4³A¼¯„_‡ê,røø#ÚÜçç•G          ek˜ÍaÕ".!ÛÅlÍlˆlX8;Ta‘5ÖÕNûşå         †¶l³5¬Ä%d³˜í¢™Á}:Üªã°Èæ›oÑ®¸êšòØ          CK6…ÙVÍá²QÌVÑÌl•¿³Cu,Ùc½Ú]÷Ì(         0¸eC˜-aÕ.!›ÄlÍÌ–Ú	aa¨Ç"{ï³_›ñË™íµ7â          ƒZ6ƒÙVMá²AÌÑÌl™[7œªC²˜<¨=0kvyœ         €-Ál«†°ía6ˆff=Şöág¡:*‹9ô°o´‡}²<V         ÀÀ’M`¶U3XÈÖ0›C3³åŞŞá‰P™ÅŒ?âÈöøSÏ”Ç         èß²Ì°jÙfchfÖk;:<ª£³˜£>¶=óŸÏ—Ç         è_²ùËö¯jÙfShf¶BöGá›á¥P¡EV_}õvÌ±Ç·GŸxº<n         ÀÊ•_¶~ÙüU-à²Ì†0[B3³¾…	a~¨ÒböİoÿvÇ]÷–Ç         è[ÙôeÛW5…l³ÌvĞÌ¬Ï72œ^Õ‘ZÌn»mWÿÇõqìş          ô±lø²å«¿B¶Ùf+hf¶Ò·Vø^ø}¨Öb¶Úê/Úyç_XC          we³—í^Õô²Ì&0Û@3³~·uÂéáíP±Å¬·ŞÇÛÉ'µæ½ÜÆA         zG¶yÙèe«W5|…lÿ²ÌĞÌ¬ßïOÂ‰áÅPµÅ6¬yäQmÆ}3Ë£	         ôL¶xÙäe›W5{…lı²ùËöÏÌl@îğğ`¨ÜRvúÜÎí_/¼¨½úÚ[å!         —Í]¶wÙàUmŞ»È¶/?3³A³¿SCuô–²Æ¨QmüßÙî¼{zy\        `¨ËÆ.[»lîªï]dË÷7ÁÌlĞn‡pY¨`i»í¶oçœ;¥½üÊ‚òà        ÀP‘-]6uÙÖUÍİ2ü(|6˜™™mÎo„ê0.eøğáíàC¾Ş¦Ş6­<Â         0Xe;—]¶tUc÷.²ÑËV/›=3³!»Qá¨p¨eiÄˆíä‰“Ú¬)3         t3gÍî´rÙÌU-İ2d“—m^6zffÖm;†Ã›¡: ¥m·ıL›4ùÔöĞìÇËƒ         E¶pÙÄm³Í¶e3·ÙŞeƒ—-™™½Ç>÷†ê¨¾«í·ÿËvê÷Oo>ştyÈ         ¿yä±§:í[6pU÷²µËæ.Û;33[mÎÿªCû®vÜq§vúgµ'z¶<ğ         °²dÛöƒÓÏì´nU÷²©Ë¶.;33ë¥­
·‡êø.ÓÎŸÿB›ü½ÓÚô”‡         V´lØ²eË¦­jİz ºlé²©33³¸MÃ·ÂŒPäeZıÚAÚ.»üª6ï¥WÊ         |PÙ¨e«–ÍZ¶kUÓÖÙÊe3—íœ™™­„Â¯Bu¨ßÓNŸÛ¹M<er»ûŞûÊ         ôT¶hÙ¤e›V5k=tÈ6.933ëGûópb˜ªşÖ[ïãmÿj—\zy›3÷…òc         ÿ/[³lÎ²=Ë­jÓz(Û·ï†láÌÌl lË01ÌÕaï‘Í>9ºpàÁmÊ¶f=Ü¾         †¬lÉ²)Ë¶,³ª={²qËÖ-›733ÀÛ"n	ÕÁï±~t6n=Û¤É§¶©?ŸÖ^]øVùA        `àËF,[±lÆ²Ë†¬jËŞ§lÙ²iË¶ÍÌÌáF„½Â”ğT¨>=¶Új«µvøl;æ¸o¶«®şq{äñ§Ê         ı_6`Ù‚e–mX6bU;ö>e«–ÍZ¶kÙ°™™ÙÛ§Â1á¦ğv¨>ïË¨5×ì|¨ûÆøvÎ¹çµÛn¿³ÍûïWÚëñ1        `åË¦+Û®l¼²õÊæ+Û¯ª	[Ù¢e“–mZ6jfff‹ö¡0.œf…êC²Ü6Üp£öåq{´ogB»ô²+Ú³.?„         ôlµ²ÙÊv+®l¹ªÆëÊæ,Û³lĞ²E333ëÑÖ
_
“ÂmáõP}h–Û°aÃÚ&›lÚvÙe×vğ!_o'OœÔ.¾ä²vû´»Û3s+?         ¼#[«l®²½Ê+[¬l²²ÍÊF«j·> lÉ²)Ë¶,³lÍÌÌÌzmŸ	G‡Ï…êcÔkFŒÑFÓvÛ}lç;0ù{§uşÒçÆ›niÓgÜßúõœöÛùÊ0        À@–mT6R3î{ ÓLe;•U¶TÙTe[•UÕ^õ²lÅ²Ëv,233³>İŸ†ıÂùáğj¨>X+ÔÈ‘#Û'>±~Ûr«­Ú¿¸Kû»¯~­>şˆ6áŸNlgõ/¿ºòªkÚµ×İĞ~òÓ›Û-·ŞŞ¦İuoûÅô_µûg>Ô~äñöäÓÏ¶9sÿ«½8ïåöò+¯¶ß*ÿ#          »l²9Êö(¤l‘²IÊ6)¥l•²YÊv)¦l™²iÊ¶)§l²yÊ9¨l¡²‰ªZ©>X¶`Ù„e–˜™™Y¿ÛÆa¯ğÏ!ÿªæ±P}Ø         Xy²íÊÆ+[¯l¾²ı233°ûã°m88œn/…ê#        @ïÉV+›­l·²áÊ–+›.33³!±Qa‹°g8*äñº0+ü.TO         Ş‘­U6WÙ^eƒ•-V6YÙfe£efffËØÚaëğ•p\87\î‡ÂïCõ        È²ÊF*[©l¦²Ê†*[ªlª²­ÊÆÊÌÌÌú`ùW@…mÂØ°oÈ¿:%L	W†kÃaj¸3L3Ã#áé07ÌóÃÂğ‡PıG         @wÙes”íQ6HÙ"e“”mR6JÙ*e³”íR6LÙ2eÓ”mS6NG‡l²}Ú6då_F6[e•ÿ êô48ßfË    IEND®B`‚ÿØÿà JFIF JJ  ÿÛ C 


ÿÛ C		ÿÀ ,," ÿÄ           	
ÿÄ µ   } !1AQa"q2‘¡#B±ÁRÑğ$3br‚	
%&'()*456789:CDEFGHIJSTUVWXYZcdefghijstuvwxyzƒ„…†‡ˆ‰Š’“”•–—˜™š¢£¤¥¦§¨©ª²³´µ¶·¸¹ºÂÃÄÅÆÇÈÉÊÒÓÔÕÖ×ØÙÚáâãäåæçèéêñòóôõö÷øùúÿÄ        	
ÿÄ µ  w !1AQaq"2B‘¡±Á	#3RğbrÑ
$4á%ñ&'()*56789:CDEFGHIJSTUVWXYZcdefghijstuvwxyz‚ƒ„…†‡ˆ‰Š’“”•–—˜™š¢£¤¥¦§¨©ª²³´µ¶·¸¹ºÂÃÄÅÆÇÈÉÊÒÓÔÕÖ×ØÙÚâãäåæçèéêòóôõö÷øùúÿÚ   ? ıS¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š‚öúßM´–êîâ+[XT¼“Láu% }hz+åoŒ_ğS/€¿üûsâ¿øL5X²>ÁátŸ0ìfÈ„sşŞ}«âo‰ßğZOx–éì>øLĞ#•¼¸®56}Bí½
¢ìEob}hö
¼×âí)ğ«áG˜¾-ø…áİâ?½iq¨Fn}ñ
“!ü¿õK_Ûö©%µıSÄñé´Z…ÀÒ,v‡ìë°8ô!Šè¼ÿ µÕ§òäñ_ììÇV·ÑíZr}„’ÇıòkªZÕ>˜Êµ8îÏ¶üuÿ ƒøáY%‹GoxÂEádÒôß&"}Úá£`=Âšğ/Ápodó#ğ—Â»x?¹s­j­&~±Gãşû5«áø'gÁïªûOÄ³/%õ;öQŸ÷aòÆ=õÿ üøoàİ‡Fğ6ƒe*ôl#iïãß­vÇ-ªş&‘Ï,\ÈøÚóş
¯ûN|@f_éMã…şÂğü·$ßÖ—&²n>2~İ?ã¨øÚĞ?#É±‹I†Øâ¯Ñ¥UUUBªŒ JZéYµ#'Œ—D~jÍğ›öÑñgü„üOã-»_xÓåÿ ¾ÁÇåUöı¢5¯Şßø#uûf¿4ø¿Mh­–[Eud}j¡ùŠßğNŒ·ßññâõóµ;–şPšoü;âŸıü!ÿ —_ü_§”UgĞó'ëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~aÿ Ã°~)ÿ ĞÁÿ øuÿ ÈÔÃ°~)ÿ ĞÁÿ øuÿ ÈÕúyEÙô;0úÕCóşƒñOşƒşÿ ÀË¯şF£şƒñOşƒşÿ ÀË¯şF¯ÓÊ(şÏ¡Ù‡Öª˜ğìŠôğş]ò5ğìŠôğş]ò5~QGö}Ì>µPüÃÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘¨ÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘«ôòŠ?³èvaõª‡æü;âŸıüÿ —_üGü;âŸıüÿ —_ü_§”QıŸC³­T?0ÿ áØ?ÿ è?àÿ üºÿ äj?áØ?ÿ è?àÿ üºÿ äjı<¢ìú˜}j¡ù‡ÿ Áø§ÿ Aÿ ÿ àe×ÿ #Qÿ Áø§ÿ Aÿ ÿ àe×ÿ #WéågĞìÃëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~aÿ Ã°~)ÿ ĞÁÿ øuÿ ÈÔÃ°~)ÿ ĞÁÿ øuÿ ÈÕúyEÙô;0úÕCóşƒñOşƒşÿ ÀË¯şF£şƒñOşƒşÿ ÀË¯şF¯ÓÊ(şÏ¡Ù‡Öª˜ğìŠôğş]ò5ğìŠôğş]ò5~QGö}Ì>µPüÃÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘¨ÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘«ôòŠ?³èvaõª‡æü;âŸıüÿ —_üGü;âŸıüÿ —_ü_§”QıŸC³­T?0ÿ áØ?ÿ è?àÿ üºÿ äj?áØ?ÿ è?àÿ üºÿ äjı<¢ìú˜}j¡ù‡ÿ Áø§ÿ Aÿ ÿ àe×ÿ #Qÿ Áø§ÿ Aÿ ÿ àe×ÿ #WéågĞìÃëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~aÿ Ã°~)ÿ ĞÁÿ øuÿ ÈÔÃ°~)ÿ ĞÁÿ øuÿ ÈÕúyEÙô;0úÕCóşƒñOşƒşÿ ÀË¯şF£şƒñOşƒşÿ ÀË¯şF¯ÓÊ(şÏ¡Ù‡Öª˜ğìŠôğş]ò5ğìŠôğş]ò5~QGö}Ì>µPüÃÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘¨ÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘«ôòŠ?³èvaõª‡æü;âŸıüÿ —_üGü;âŸıüÿ —_ü_§”QıŸC³­T?0ÿ áØ?ÿ è?àÿ üºÿ äj?áØ?ÿ è?àÿ üºÿ äjı<¢ìú˜}j¡ù‡ÿ Áø§ÿ Aÿ ÿ àe×ÿ #Qÿ Áø§ÿ Aÿ ÿ àe×ÿ #WéågĞìÃëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~aÿ Ã°~)ÿ ĞÁÿ øuÿ ÈÔÃ°~)ÿ ĞÁÿ øuÿ ÈÕúyEÙô;0úÕCóşƒñOşƒşÿ ÀË¯şF£şƒñOşƒşÿ ÀË¯şF¯ÓÊ(şÏ¡Ù‡Öª˜ğìŠôğş]ò5ğìŠôğş]ò5~QGö}Ì>µPüÃÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘¨ÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘«ôòŠ?³èvaõª‡æü;âŸıüÿ —_üGü;âŸıüÿ —_ü_§”QıŸC³­T?0ÿ áØ?ÿ è?àÿ üºÿ äj?áØ?ÿ è?àÿ üºÿ äjı<¢ìú˜}j¡ù‡ÿ Áø§ÿ Aÿ ÿ àe×ÿ #Qÿ Áø§ÿ Aÿ ÿ àe×ÿ #WéågĞìÃëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~aÿ Ã°~)ÿ ĞÁÿ øuÿ ÈÔÃ°~)ÿ ĞÁÿ øuÿ ÈÕúyEÙô;0úÕCóşƒñOşƒşÿ ÀË¯şF£şƒñOşƒşÿ ÀË¯şF¯ÓÊ(şÏ¡Ù‡Öª˜ğìŠôğş]ò5ğìŠôğş]ò5~QGö}Ì>µPüÃÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘¨ÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘«ôòŠ?³èvaõª‡æü;âŸıüÿ —_üGü;âŸıüÿ —_ü_§”QıŸC³­T?0ÿ áØ?ÿ è?àÿ üºÿ äj?áØ?ÿ è?àÿ üºÿ äjı<¢ìú˜}j¡ù‡ÿ Áø§ÿ Aÿ ÿ àe×ÿ #Qÿ Áø§ÿ Aÿ ÿ àe×ÿ #WéågĞìÃëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~aÿ Ã°~)ÿ ĞÁÿ øuÿ ÈÔÃ°~)ÿ ĞÁÿ øuÿ ÈÕúyEÙô;0úÕCóşƒñOşƒşÿ ÀË¯şF£şƒñOşƒşÿ ÀË¯şF¯ÓÊ(şÏ¡Ù‡Öª˜ğìŠôğş]ò5ğìŠôğş]ò5~QGö}Ì>µPüÃÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘¨ÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘«ôòŠ?³èvaõª‡æü;âŸıüÿ —_üGü;âŸıüÿ —_ü_§”QıŸC³­T?0ÿ áØ?ÿ è?àÿ üºÿ äj?áØ?ÿ è?àÿ üºÿ äjı<¢ìú˜}j¡ù‡ÿ Áø§ÿ Aÿ ÿ àe×ÿ #Qÿ Áø§ÿ Aÿ ÿ àe×ÿ #WéågĞìÃëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~aÿ Ã°~)ÿ ĞÁÿ øuÿ ÈÔÃ°~)ÿ ĞÁÿ øuÿ ÈÕúyEÙô;0úÕCóşƒñOşƒşÿ ÀË¯şF£şƒñOşƒşÿ ÀË¯şF¯ÓÊ(şÏ¡Ù‡Öª˜ğìŠôğş]ò5ğìŠôğş]ò5~QGö}Ì>µPüÃÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘¨ÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘«ôòŠ?³èvaõª‡æü;âŸıüÿ —_üGü;âŸıüÿ —_ü_§”QıŸC³­T?0ÿ áØ?ÿ è?àÿ üºÿ äj?áØ?ÿ è?àÿ üºÿ äjı<¢ìú˜}j¡ù‡ÿ Áø§ÿ Aÿ ÿ àe×ÿ #Qÿ Áø§ÿ Aÿ ÿ àe×ÿ #WéågĞìÃëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~aÿ Ã°~)ÿ ĞÁÿ øuÿ ÈÔÃ°~)ÿ ĞÁÿ øuÿ ÈÕúyEÙô;0úÕCóşƒñOşƒşÿ ÀË¯şF£şƒñOşƒşÿ ÀË¯şF¯ÓÊ(şÏ¡Ù‡Öª˜ğìŠôğş]ò5ğìŠôğş]ò5~QGö}Ì>µPüÃÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘¨ÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘«ôòŠ?³èvaõª‡æü;âŸıüÿ —_üGü;âŸıüÿ —_ü_§”QıŸC³­T?0ÿ áØ?ÿ è?àÿ üºÿ äj?áØ?ÿ è?àÿ üºÿ äjı<¢ìú˜}j¡ù‡ÿ Áø§ÿ Aÿ ÿ àe×ÿ #Qÿ Áø§ÿ Aÿ ÿ àe×ÿ #WéågĞìÃëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~aÿ Ã°~)ÿ ĞÁÿ øuÿ ÈÔÃ°~)ÿ ĞÁÿ øuÿ ÈÕúyEÙô;0úÕCóşƒñOşƒşÿ ÀË¯şF£şƒñOşƒşÿ ÀË¯şF¯ÓÊ(şÏ¡Ù‡Öª˜ğìŠôğş]ò5ğìŠôğş]ò5~QGö}Ì>µPüÃÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘¨ÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘«ôòŠ?³èvaõª‡æü;âŸıüÿ —_üGü;âŸıüÿ —_ü_§”QıŸC³­T?0ÿ áØ?ÿ è?àÿ üºÿ äj?áØ?ÿ è?àÿ üºÿ äjı<¢ìú˜}j¡ù‡ÿ Áø§ÿ Aÿ ÿ àe×ÿ #Qÿ Áø§ÿ Aÿ ÿ àe×ÿ #WéågĞìÃëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~aÿ Ã°~)ÿ ĞÁÿ øuÿ ÈÔÃ°~)ÿ ĞÁÿ øuÿ ÈÕúyEÙô;0úÕCóşƒñOşƒşÿ ÀË¯şF£şƒñOşƒşÿ ÀË¯şF¯ÓÊ(şÏ¡Ù‡Öª˜ğìŠôğş]ò5ğìŠôğş]ò5~QGö}Ì>µPüÃÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘¨ÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘«ôòŠ?³èvaõª‡æü;âŸıüÿ —_üGü;âŸıüÿ —_ü_§”QıŸC³­T?0ÿ áØ?ÿ è?àÿ üºÿ äj?áØ?ÿ è?àÿ üºÿ äjı<¢ìú˜}j¡ù‡ÿ Áø§ÿ Aÿ ÿ àe×ÿ #Qÿ Áø§ÿ Aÿ ÿ àe×ÿ #WéågĞìÃëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~aÿ Ã°~)ÿ ĞÁÿ øuÿ ÈÔÃ°~)ÿ ĞÁÿ øuÿ ÈÕúyEÙô;0úÕCóşƒñOşƒşÿ ÀË¯şF£şƒñOşƒşÿ ÀË¯şF¯ÓÊ(şÏ¡Ù‡Öª˜ğìŠôğş]ò5ğìŠôğş]ò5~QGö}Ì>µPüÃÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘¨ÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘«ôòŠ?³èvaõª‡æü;âŸıüÿ —_üGü;âŸıüÿ —_ü_§”QıŸC³­T?0ÿ áØ?ÿ è?àÿ üºÿ äj?áØ?ÿ è?àÿ üºÿ äjı<¢ìú˜}j¡ù‡ÿ Áø§ÿ Aÿ ÿ àe×ÿ #Qÿ Áø§ÿ Aÿ ÿ àe×ÿ #WéågĞìÃëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~aÿ Ã°~)ÿ ĞÁÿ øuÿ ÈÔÃ°~)ÿ ĞÁÿ øuÿ ÈÕúyEÙô;0úÕCóşƒñOşƒşÿ ÀË¯şF£şƒñOşƒşÿ ÀË¯şF¯ÓÊ(şÏ¡Ù‡Öª˜ğìŠôğş]ò5ğìŠôğş]ò5~QGö}Ì>µPüÃÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘¨ÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘«ôòŠ?³èvaõª‡æü;âŸıüÿ —_üGü;âŸıüÿ —_ü_§”QıŸC³­T?0ÿ áØ?ÿ è?àÿ üºÿ äj?áØ?ÿ è?àÿ üºÿ äjı<¢ìú˜}j¡ù‡ÿ Áø§ÿ Aÿ ÿ àe×ÿ #Qÿ Áø§ÿ Aÿ ÿ àe×ÿ #WéågĞìÃëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~aÿ Ã°~)ÿ ĞÁÿ øuÿ ÈÔÃ°~)ÿ ĞÁÿ øuÿ ÈÕúyEÙô;0úÕCóşƒñOşƒşÿ ÀË¯şF£şƒñOşƒşÿ ÀË¯şF¯ÓÊ(şÏ¡Ù‡Öª˜ğìŠôğş]ò5ğìŠôğş]ò5~QGö}Ì>µPüÃÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘¨ÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘«ôòŠ?³èvaõª‡æü;âŸıüÿ —_üGü;âŸıüÿ —_ü_§”QıŸC³­T?0ÿ áØ?ÿ è?àÿ üºÿ äj?áØ?ÿ è?àÿ üºÿ äjı<¢ìú˜}j¡ù‡ÿ Áø§ÿ Aÿ ÿ àe×ÿ #Qÿ Áø§ÿ Aÿ ÿ àe×ÿ #WéågĞìÃëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~aÿ Ã°~)ÿ ĞÁÿ øuÿ ÈÔÃ°~)ÿ ĞÁÿ øuÿ ÈÕúyEÙô;0úÕCóşƒñOşƒşÿ ÀË¯şF£şƒñOşƒşÿ ÀË¯şF¯ÓÊ(şÏ¡Ù‡Öª˜ğìŠôğş]ò5ğìŠôğş]ò5~QGö}Ì>µPüÃÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘¨ÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘«ôòŠ?³èvaõª‡æü;âŸıüÿ —_üGü;âŸıüÿ —_ü_§”QıŸC³­T?0ÿ áØ?ÿ è?àÿ üºÿ äj?áØ?ÿ è?àÿ üºÿ äjı<¢ìú˜}j¡ù‡ÿ Áø§ÿ Aÿ ÿ àe×ÿ #Qÿ Áø§ÿ Aÿ ÿ àe×ÿ #WéågĞìÃëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~aÿ Ã°~)ÿ ĞÁÿ øuÿ ÈÔÃ°~)ÿ ĞÁÿ øuÿ ÈÕúyEÙô;0úÕCóşƒñOşƒşÿ ÀË¯şF£şƒñOşƒşÿ ÀË¯şF¯ÓÊ(şÏ¡Ù‡Öª˜ğìŠôğş]ò5ğìŠôğş]ò5~QGö}Ì>µPüÃÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘¨ÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘«ôòŠ?³èvaõª‡æü;âŸıüÿ —_üGü;âŸıüÿ —_ü_§”QıŸC³­T?0ÿ áØ?ÿ è?àÿ üºÿ äj?áØ?ÿ è?àÿ üºÿ äjı<¢ìú˜}j¡ù‡ÿ Áø§ÿ Aÿ ÿ àe×ÿ #Qÿ Áø§ÿ Aÿ ÿ àe×ÿ #WéågĞìÃëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~aÿ Ã°~)ÿ ĞÁÿ øuÿ ÈÔÃ°~)ÿ ĞÁÿ øuÿ ÈÕúyEÙô;0úÕCóşƒñOşƒşÿ ÀË¯şF£şƒñOşƒşÿ ÀË¯şF¯ÓÊ(şÏ¡Ù‡Öª˜ğìŠôğş]ò5ğìŠôğş]ò5~QGö}Ì>µPüÃÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘¨ÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘«ôòŠ?³èvaõª‡æü;âŸıüÿ —_üGü;âŸıüÿ —_ü_§”QıŸC³­T?0ÿ áØ?ÿ è?àÿ üºÿ äj?áØ?ÿ è?àÿ üºÿ äjı<¢ìú˜}j¡ù‡ÿ Áø§ÿ Aÿ ÿ àe×ÿ #Qÿ Áø§ÿ Aÿ ÿ àe×ÿ #WéågĞìÃëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~aÿ Ã°~)ÿ ĞÁÿ øuÿ ÈÔÃ°~)ÿ ĞÁÿ øuÿ ÈÕúyEÙô;0úÕCóşƒñOşƒşÿ ÀË¯şF£şƒñOşƒşÿ ÀË¯şF¯ÓÊ(şÏ¡Ù‡Öª˜ğìŠôğş]ò5ğìŠôğş]ò5~QGö}Ì>µPüÃÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘¨ÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘«ôòŠ?³èvaõª‡æü;âŸıüÿ —_üGü;âŸıüÿ —_ü_§”QıŸC³­T?0ÿ áØ?ÿ è?àÿ üºÿ äj?áØ?ÿ è?àÿ üºÿ äjı<¢ìú˜}j¡ù‡ÿ Áø§ÿ Aÿ ÿ àe×ÿ #Qÿ Áø§ÿ Aÿ ÿ àe×ÿ #WéågĞìÃëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~aÿ Ã°~)ÿ ĞÁÿ øuÿ ÈÔÃ°~)ÿ ĞÁÿ øuÿ ÈÕúyEÙô;0úÕCóşƒñOşƒşÿ ÀË¯şF£şƒñOşƒşÿ ÀË¯şF¯ÓÊ(şÏ¡Ù‡Öª˜ğìŠôğş]ò5ğìŠôğş]ò5~QGö}Ì>µPüÃÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘¨ÿ ‡`üSÿ  ÿ ƒÿ ğ2ëÿ ‘«ôòŠ?³èvaõª‡æü;âŸıüÿ —_üGü;âŸıüÿ —_ü_§”QıŸC³­T?0ÿ áØ?ÿ è?àÿ üºÿ äj?áØ?ÿ è?àÿ üºÿ äjı<¢ìú˜}j¡ù‡ÿ Áø§ÿ Aÿ ÿ àe×ÿ #Qÿ Áø§ÿ Aÿ ÿ àe×ÿ #WéågĞìÃëUÌ?øvÅ?úø?ÿ .¿ùøvÅ?úø?ÿ .¿ù¿O(£û>‡fZ¨~b/ü/âÍ¯Í¿áMßôÎúèı'õÿ ‚}|t±ÿ é#óÇX¸_ı¦+ôâŠ?³èy‡Öª™ğşË?µ_„ÿ äâ]QvôşÌñ[ÃùfD«°Ãûrø7çOøòq 7ˆWP÷É™óôÅ~‘ÑPòÚ=+ësìÎë_Û‡öÓøs…ÖWZ»·¶³áHÊûh°©?÷Õu~ÿ ‚Ó|XĞî)ğ?†5uŒá–Ù.,g?Rd‘Aÿ €WÜÕ—®øWDñDN³£Øjğã]õªN¸ú05„²µögø,cê ğ_üÛáş£å¯Š¾ø‹Bvá›K¹‚ıúåÌ'AŸjú'áßüƒöyøåGmñËC»|oâßOØO¬’äæ¾ñ7ìcğgÅ{ÍÏ,,än¦4–{}ÂÄÊ¿˜Åx·Œ¿à—¾Ô¼É<3â­_C‘¹ßGä@úylÕrË.­¬Í£Š¦÷ĞıpĞ|I¤ø«OKıT²Ö,_îİX\$ñ7ÑkF¿o¿`ï_µÕüâ(î®£å.4RM>ó÷Šø9®‡@ÿ ‚„~Õ¿³ÄV^4†}jÊ6
!ñ†”Çxï¶å6;Ÿrì3ëÒ¸gF¥?6:#R2øYûƒE~qüÿ ‚ÒxÄFOˆÔ¼rÜ5şšßÚŸï2€²(ö
ÿ Zû{áoÇÏ‡_ì…Ïüg£ømÜĞÙÜ©1şÜ''ü	EbhwÔQE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE Q^eñÃö’øqû:èÚ^<ñEŒ][Øîóo.ºÿ ªrì21»AêE zmy÷Æß~è¿Ú<ñ^ŸáèKEÄ›®'ÆxŠÌ’?…M~Zü~ÿ ‚¸|Dø±ª7…ş
èW³ºs7¦w¬]{F€2E‘;w+Í¼ûüNøÉ­?‰¾*øŠçI{Æóg’şàßj—ïÄ&ÚbG÷kjtjVv‚¹©+Éÿ ñÓşTÅ®4ß„~’‰®x›œöÜ–Ñ·â?¦Sµ|Õ¨x'ö¤ıµn’ÿ ÅÚ–¨Ú,æDÚôÆÃN'ƒª(ÏûÉè2kî?„¿²ÏÃoƒº‡¡›TŸímKdú‡a„ÿ €Õë5ìRËzÕqÁ<_H#âÿ †ßğLèUÇŒõëïÜZÒÌ}ÛèH&Fú†_¥}Gà_„¾øej ğ·†tİcK[u¿ûÒ™¿k­¢½Zt)Rø"qÊ¤çñ0¢Š+ È(¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š *+«X/­ä·¹†;ˆ$^)T2°ô ğjZ(Á>%~Ãÿ 	>$y³?‡WÃº„™?lĞXZœæ<Ï>©ŸzùgÇŸğMßø"ùu‡&‹Y’İ¼Ø#y…ôl:ynÂzüÛ“é_¤4WL%»ÆÏÈŞ5êCf~xøş
3ûI~Ëú¤:Ä;+ØGÇØü[¥Ó(ÆLWcæ÷ŸÍ×Şß³ÿ ü3àÏÆ¹-ôí^úO‡~!”…~ u[iöè|‡¨?–OaSx‹Ã:G‹´¹tÍsL³Õôù~ı­ô4mÿ `Fkä¯Œ_ğMŸø«Ï¾ğ-üÔ[-ö)‹\Y9ô;ãïĞ°–¼Š¹lã­7sº¨½$¬~¸Á<wPÇ42,ĞÈ¡ÒHØ2²‘AA%~xâí;ÿ ÿ ¾‰`¼¼ÿ „R9xµ¹'PĞçÉèrXãî˜œãÒ¿@fŸø+7Ã‹ÿ dÑüt£á·‰¤Ây—²ïÓ'nŸ-Æ•N%
MìkÉ”ei+3µII]tQQZİC}mÅ¼±Ü[Ì‚Hæ‰ƒ#©à‚9È©jFQE QE QE QE QE QE QE QE QE QE QE QE QE QE QE ƒãxwá¯†o<Câ­jË@Ñ,×t×·óã_A“Õ@£$ &¾fı¯ÿ à£~ ı—ÒçB±)ã/€WûÎ`#²ld©vu ËŸEu~bÜAñóş
3ã5Ö¼C¨I‡ ”ˆî%VƒI°-âï|pq¹Øu«Œ%QòÅ]“)(«³éÚƒş	­I?†>é²[‰[Èÿ „£R·İ<„ğ>ËlGÊsÑ¤œıÀy¯	øcûüGøñ®Iã‹:æ¡¤%ûùó> æ}Vë>¡Éò¿ày#û˜¯¯şşÉ^ømúuˆÕüDkš‚†Ÿ'¨Œt‰İçI¯i¯v†\—½[_#Î©ŠoHğŸà/¾	éßfğ…”Ì»f¿y—Sÿ ¿+|ÄgøFv½Š+ÙŒTU¢¬ÜØQEB
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€"»´‚şÖ[k¨c¹·•JIÈOPAàŠù7ãwü¯Á^=ûF£à¹G‚õ¦Ë}4/a+zúÇÿ  àt×ÖôV5)B²´ÕÍ#9AŞ,üÊøñßö‡ÿ ‚wxŠú9nü&òi‘k*äXÛÊîŸü¤œº•ú‹û,ÿ ÁB>~ÔK™gv|-ã6\¿‡uiI!ïöyxYÇ°Ãã’€W=â/é^.Ñît­oN¶ÕtÛ…Û-­äBHÜ{ƒüûWÂ´ü¦}>i|KğâHå‰¼ÿ ì	ç!ã#m¦'9p®sèİx8Œ¾P÷©ê¿Ñ§ŠRÒzµtWã¯ì»ÿ Tñ·Á}Q<ñÊÃQ×ô»7­©Í]_O àùÊØóÔïbNù~~³ü?ø‹áŸŠŞ²ñ/„u»?hWƒ0ŞÙI½	TªÃ¡V ƒÁ¼ƒ¸èè¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢¸_Œÿ <ğÀ·¾-ñ¶¯•¥[Œ"Ÿšk™0JÅ
uwlp¹$ H êõİ{Mğ¾y«köÚ^—gMsyy*Å1’ÌÌ@ zšüı±¿àªÚï5+¿ |V–IöVñE¼oöëæ'm#ÆèÔö|y‡°Nş-ñÛöø¿ÿ ñğğ†4é´ßÃ0–ßA·|Cq}7FaÔgå^Š¥¹o¬¿fŸÙÂÿ ³İŠ_ºßŒ&mÆ±2qGÌ)û‹ØŸ¼İÎ8¸|,ñMs
µ£Myû6ÿ Á=eº¸‹Å7Ï<¬.#ğÿ šY™‰İºéÇROTıãÕkï?N´Ò,`²±µ†ÊÎÅoo8ÔtUP0°«WÓQ¡
´äT©*ò
(¢º‚Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( &øñû2ø/ö€ÒLZåŸØõ˜”‹]jÍBÜÂqÀcÿ -ı†ã®0y¯…t½Cã—ükâ"êZEÛ\ønîaæ`<ºNª£¢L™ıÜ˜éÈqÎÖ#9ıC¬ÿ xMñV‹y¤k0jZeäf)ínP<r)ìAÿ  Œ×#
ú­$tÒ¯*znÇöCı¸¼ûZx}N•t/ÛÇºÿ Ã7S4xÆd…°<è³ü@?ˆ.F~¯ÂßÚö4ñWÀCñáæ¤Ú~/ÚÕ,åo·éL¹;Ñ‡Í$`gŸ¼ŞÈËWÙÿ °oüëNøØö~ø£=¦‡ã¶ÄVZ¨M]º#¤SŸîı×?w_3V”èË–hõá8ÔW‰úEV%…Q@Q@Q@Q@Q@Q@Q@Q@Q@W‚ş×ßµ÷…?d‡çXÖ
ê~"¾ Å Yo$–cü.FçÇ eˆ[ö ı©üû)øOx¦ãíÓLÑmÜB`>êƒ÷Pdnğ ÷$)üq¸“âßü³ã$úŞ·vl|?dû<Àû?H€œù0¦~yHÆ‰ Âx+Àÿ à¡Ÿ¯¼oãNx´%—eÖ ªV#+if‡ `}¹ÜÅ˜üß¤şğ.‡ğÛÂö>ğæŸ™¤Ù¦È ˆuõf=Y‰ä±ä“^¦ë{óø3µgîÇsàïÁŸ|ğŒ>ğÅ—‘qu&{¹1Ì’7sè:À WsEô±ŠŠåŠĞò[mİ…QT ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¯ˆ¿kØ>5÷>[%¦·Ì÷šXXîÏRğvI;•èİ°ß{íÚ+
ÔaZ<³F©*nñ>Sı€à¦w:-ÕÂ¿·ÒÆ"qg¦ø£P$I´[Ş–ç ü¢SÈ<?0ıXGYYX2°Èe9z×å/íûÙ|d¶¹ñ_„¡†ÃÆñ&ébá"ÔÔû,¸ß€Ü`¯1ûÁDµO‚ºÕ·Â/ŒÓ\EáØ&û¯¨'Ñdh‚ãw& xyıÏ¹ò¸Œ<ğò´¶î{4êÆ¢º?`è¦C4w¤±:ËŠUâŸ\¦ÁEPEPEPEPEPEPEPEÂüløÏá€?µx¶ôZiZ|y©mÌ§îCŸ¼ìxêN $ r¿µGíAá_ÙOáÇŠ|E ¹¾›t:V6¡qŒ„_DŸQêJƒøáà?øÿ ş
ñÃSñÇï¦M
9€½¼ŒŠÊYZ©Èß,Ùfù‰.>!ÿ ÁK?h‹½o[M;Ãö„òÉk}"ÇqÙYá¥nyÆY·1 ÒøDømá]?Ã±NÒlcòâ…:ŸVcÕ˜œ’Ç’I5ê`ğÙóÏáüÎ:õ½šåå¯ø_Jğ_‡ìt=Æ7J±ˆCokáQGó$òIä’IÉ5©EôËMä…QLAEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEP_/~Ù_²¯ÇOøn­|ueÂ.¥"sığ>ëŸ÷O+õ•Jq«	­Œœ2>:ÿ ‚nşß÷_uKO‚ß®å´Òc›ìZ>­¨’¯¦Jß²N[‘xV?êÏåÆÏ×:ü’ı·¿cøş(i×>9ğ}š§Œ-#İygÇöœJ:ÿ =”;°y!k¹ÿ ‚]şŞÒø¢;/ƒõøÛ¯‘áİZñğ×H£ÎV=dP>B~ğ~ğ¾KBXyò½jER7Gé¥Q\ÆÁEPEPEPEPEPEP{]Óü/¢ßëµä:v—aÜİ]Ü8Há‰³;Ğ 	¯Ã_Ú{ã·Šà¢ß´vá‹ˆ<e+A£ÚÌ8‡úëû€:3pyUÚƒæ'w´ÁU¿lkßx¬|ğÜ·÷ˆÏ%¯¯w–jGŞTlnäÀş}'öDıšmg¿ƒ|ÜxÃTU“SºL7–:­º7÷W¹y²zc¸\;ÄNİæªªqó=àÏÁİàw,¼1áøq_¼¸º6îbé\úœp:  
îh¢¾²1QJ1Øñmİ…QT ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¯ÏßÛËöWoÜMñkÀñIhÉ0¸Ö-lò­™È¼o+óc~:?÷~Twñ^[ËoqOªRH¤PÊêF
x ÕÏ^ŒkÃ’F´ê:ræFOüöÚ·ı§>/†üIv«ñ+Ãöê/ƒ¥ !Ví¯A ƒÀpÙø'ñëá·Š?a¿Ú/Ä‡—Oa£ÍvÓé’®Y-ß–ÎQŸš6RÀ÷‘ÕI¯Ù¿Ù¯ö€ğÿ í1ğGñÇ‡ÛË[•òo¬Xåì®Ô6ú?Ä¬­Ş¾B¥9R“„·G¹)®dzQY”QE QE QE QE òGüoö¿_Ùá	±Ğ®Txÿ Ä©%®•±†ë(ÀÄ—d³ÕØu
ÕôÏ<q£|5ğnµâ¯Ş%†‰¤ZÉywpçî¢v=I M~Aqâ/ø(Ïíu¨x‡ZY­¼8$š Ø\lDVÊÃí’2:³ÈøàÕÂ.¤”c»&RQWg¥ÿ Á=fÙn§ÿ …¹â¸yågşÅ†èfbH{¶ÏRNBŸ÷›û¦¾øªúnŸm¤iö¶6P%µ¬KF0±Æ *¨€ Â¬WØP£PG‡R£©.fQEtQ@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@ŸÅ_†z?Åÿ êŞ×"ße}Õ•T‚AÊJ™èÊpGäx&¾ı“ş8x‹ş	ïûNjñƒÊ|¨Î–ZÄkŸ(ÆOî5ë´œrPºıà1úO_-şŞŸ³¯ü-Ï‡Ÿğ“èÖŞgŠü;Hªƒæº´å¤‹Ùyu÷ÜŞ¯/‡ö°çèìÃÕä—+ÙŸ©¶·PßZÃsm4wó ’9¢`ÈêFC)AÎEK_œğHÏÚãşÏ¿ÁÏ^n×ü;›¡Í+ ntñ€a­	<ùæGd&¿G+æ\(¢Š (¢Š (¢Š (¢¼Ëö’øá¦~Î¿|QãÍKd‡M¶?cµsµ]?Ë]s†r¹Ç!CÔùÙÿ „ı¨ZÕ´ß~§¼WŞ û9Üd˜àÛZ`rJäHGrÑwS^“û%|‡à/Â[:xTx‹Q{«Ê9>qD÷c\/¦wâ¯ÿ aŸ†Z—ÇŸšçÅŸÈú¢i×­z×?jÔä%ÃÛ<ïÇbcí_¥5ôu/m.»f*¥ß"
(¢½³Ï
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€?1i¯ë?±ïí1 üKğ:ı†Âæ÷ûWOÚ—À?éÍşÃ?/÷$*>é¯Úÿ ‚?´O
ü9ã¯I;XµY¼’ÁŞQÄ¾?‰2Ÿuô¯ÿ h¯ƒV¾êş¸—¬¿hÓ®\¨º@|¶ú•?ì»WÎßğH¯Ú"óá§Å]oàŠe{K=fi$Ó¡¸l}—Sˆ,<œ1ÿ ÀâP9c_+¡ìj]lÏgS6{£ö
Š(¯8ê
(¢€
(¢€
ütÿ ‚¸|~½ø±ñ›Bø+áv’öÏ@3woÏÚµYÀTĞùháGûSHJıJøıñƒMøğoÅ~<Õ6´5“M.qçÎp°Ä=ŞFEÿ Wãgìà=Kã'Ç_|Uñ3µûéóËx÷3söJà³÷ÚÛØ²Ú7Z¢‚êEI(EÉŸs|øOiğOáN…á;m5¬;ï.¯¹šWún$E
;W QE}œb¢”VÈğmİ…QT ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¯ÍßÛûá£ğã…ñkÂÏ&ı¡uíuoÁµÔá!ÒAè\ aşÔnOZı"¯<ı >Áñ£áˆ|) Qsu™e+Ë;”ù¢lö€Ø‘Ş¹1T}µ'½èÔösLú{öaøçaûG|ğ·,‚C.£m¶úÕ×hvO<à8;sÕJõêuù	ÿ søñ?‚>%x›àÆ½#ZÃ¬´—ÚtœyZ„+¶x€ÏŞx“'ş½ıëõî¾<÷Š( Š*ëØ4Û+‹»©V[xÚYeáQe˜Ÿ@4ùcÿ ªøèÅüğM¸=?·µdºıèí£lÛg*é™ÇJôïÙgá*üø#áíX|­RH¾İ©d|ÆêPÁÿ tmè‚¾ğN¡?í«û{j^.¿G—E:„šÃC.O—alU-b#ğ§Şc_§•ïe´·ªı7= ‚Š(¯tó‚Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Ì¿Û'ÃÚŸìáûTè|2¯Ûî£×mYr^BêgCêífüæ+öÿ áÏôÏŠğï‹´i<Í/\°‡P·9É"Ú}gv Šüèı½>/ÄÙÿ U½‚ úŸ†Ïö¼|´N¹ôòË7Õºßø#OÇã/‚ºçÃ›ùËê»óìÃ7&Êà³»%çĞH‚¾OKÙVvÙê{XyóÁy¡´QEp!_+ÁL¾1ÂŸı‘üVmçòu_mğıûğ|â;ñÍÏ®+êšü}ÿ ‚Ò|N›Äß¼ğŞÁÚtÑìüğGÎë«§ØˆGvÄ¤{Mï@Á1şb|7×¼gq.5Ë¿²Û1òïA#ë#8?õÌWÚÉ|%ğ,?¾øgÂĞ —a»²ÿ s#ÿ Àœ³~5Ö×ÙĞ§ì©FƒR\ór
(¢º‚Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( ®­a¾µšÚâ5š	‘£’7¤`ƒìE~s~Å¾$—öNÿ ‚…'…/fh4JúoÎÒ6Ã;+Z9Ï«‹s“Ù~Wçü»ÁsxSâ—„üy¦—µ—R·òx¸)slÊQóÙ¶:ÿ \ëÉÌióRSìváeiò÷?sh®àoÄˆ~/üğ_!+sI·½‘W¤r´`ÈŸU}Ëÿ ®æ¾hõ‚¿muCûTÿ ÁIuM}Ú´„×å½ÉÊıÈm·Èôa Ÿ=~Ê~ÒŸ¿áT~Ïÿ ¼Z²yW^‰u5³gé2°ŒûÈÈ?üÿ ‚Zø7ÏñüW$|[[C¦Bäu21’@>T÷Õuaaí+F&5¥ËM³ô:Š(¯±< ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¯›?à ÿ „ÓösÕ/b}æs©;1Ê>$f?î
úN±üeá¸<eás@¹ÇÙõKì¤Ü26É!ÿ Ğ«*°ö”;—	rÉHäà?Œ?f}CÂwï»ğ–­$¡9Å­Çï£?÷ğÜ¢ŠûÒ¿à>4›ÁŸ´wŒü~M¹Ö´‡Ì$õº´”¸öI.á_²Õñ'ĞÁ`ütŞı‘ÛGŠ]²x“[´ÓİäÄ›îXı7Aÿ 
ñø'o„WÃ³m…ùM³k—÷7ìHç¼…ü1GûŞõ•ÿ Áñ§™ª|+ğ”o*İVâ<õŞÑEü6Mù×º|ğÏü!¿ü£Ù-®j%\c÷¦5i?ñâÕëå±½W.ÈáÅËÜHï(¢ŠúCÊ
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€?6¼3x~ÿ ÁR´{˜›ìÖ·*„3tQ¢¡dÿ €tßM¾Õû§_„_ğR-2ãÂ?´†<Ub|™®tØfI1ÿ /ó??‚˜«÷#Â¾ ·ño…ô}rÓ›]NÎØyÏÉ"_Ñ…|f"<•¥3Ş¥.h&~0ÿ ÁUïÄÛ«HğÚ±²éº^‘°†Y^\}Ò~‚*¬jªªT`( zWç/ÆKƒñ#ş
½¨–ıà·ñ]´x=1ekü—¯Ñºö2Èû²‘ÃŒz¤QE{gQE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE |/ÿ NĞÄŞğ°›kÛ«BØÿ ˆşA?­~şÂş%ÿ „³ö@øK|_Ì1èöE³mÁ·?ú*¾ÿ ‚’h£Tıœ¾Õ·'MÖ-nsé¸IşÕ¯£?à’~":çì_áûBÛ¿²5=BÄ{fs>?ò=|¶aWo½c
ïLüñøM7ü%ŸğS¯ê~/øIüAv?İÿ J	ùn_Ê¿J«ó+ös­~Ù#¿—çÅ©Üîÿ i¦POş>kôÖ½Lµ~åúœ˜¯âQ^©ÄQE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE x7íÓ§ÿ i~Ë8P2Ñ%¬ãÛeÔ,@j—üCâÅ§ƒfÏi÷…KÂUs"n=ÚZqùƒù×_ûWZ­çìßñ6£Í'â£pıE~{~Ì/<)à;ûKyü¤}NIH÷1D?öZù¼Í~ö/ÈõpŸGMÿ áo·~Ó ¸ëÿ KÉsõ¹€ìÕúu_˜ğLù/šÿ ı‹7úUi_§•èåÿ Àùœ¸¯âQ^™ÈQE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE yÇí%û>|H_O_7åŸé_>	ñ'ö>•,;öî˜¿ş:£úWì7í	ÿ $âWı‹:Ÿş’É_ˆ›ˆèkçs?>‡©„øYõ÷üşKæ¿ÿ bÍÇş•ZWéå~aÿ Á0ä¾kÿ ö,ÜéU¥~W~_üêsb¿ˆQEzg QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QEÀ~ĞŸò@ş%Ø³©ÿ é,•ø‰_·´'ü?‰_ö,êúK%~"WÏf}S	ğ³ëïø&ü—ÍşÅ›ı*´¯ÓÊüÃÿ ‚`ÿ É|×ÿ ìY¸ÿ Ò«Jı<®ì¿øÔæÅ(¢ŠôÎ@¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š(€ı¡?äüJÿ ±gSÿ ÒY+ñ¿nÿ hOù ¿ìYÔÿ ô–JüD¯Ìş8ú¦ág×ßğLù/šÿ ı‹7úUi_§•ù‰ÿ ¿ÿ ’ù¯ÿ Ø³qÿ ¥V•ú}]¹{ıÇÌæÅeú+Ò¹È2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢àyïí	ÿ $âWı‹:Ÿş’É_ˆ•û{ûBÉø—ÿ bÎ§ÿ ¤²Wâ|şgñÄõ0Ÿ>¿ÿ ‚_ÿ É}×ÿ ìY¸ÿ Ò«Jı>¯Ìø%ÿ ü—İşÅ›ı*´¯ÓêíËÿ ó9±_Ä
(¢½#(¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š óÿ ÚşHÄ¿ûu?ı%’¿«ö÷ö„ÿ ’ñ/şÅOÿ Id¯Ä*ğ3/'©„øYöü÷şKö¿ÿ bÍÇş•ZWê~_ÿ Á/ä¿kÿ ö,ÜéU¥~ Wn_ü™ÍŠş QEéEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPŸşĞ¿ò@~%ÿ Ø³©ÿ é,•ø_¸´/ü‰ö,êúK%~×™|q=L'ÂÏ°?à—¿ò_µÿ ûn?ôªÒ¿P+óş	yÿ %û_ÿ ±bãÿ J­+õ»pÁùœØ¯â¢Ez79ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\ÑN¢‹€Ú)ÔQpE:Š.h§QEÀmê(¸¢E´S¨¢à6Šu\=ı¡ä€üKÿ ±gSÿ ÒY+ğş¿p¿hoù ??ìXÔÿ ô–Jü=¯2øâz˜O…Ÿ`ÿ Á/?ä¿kÿ ö,\éU¥~ ×å÷üïşKö¿ÿ bÅÇş•ÚWêÓ]¸àœØ¯â	E.ÓFÓ^Ê%»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM %»MM y÷íÿ $âgı‹Ÿş’É_‡µû‡ûCÿ 
âgı‹Ÿş’I_‡•àæ_OO	ğ³ì/ø%×ü—íşÅ‹ı+´¯Ô:ü¼ÿ ‚]É~×ÿ ìX¸ÿ Ò»JıC®ìğNlWñŠ(¯@å
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€<ûö†ÿ ’ñ3şÅOÿ I%¯ÃºıÄı¡¿ä€üLÿ ±cSÿ ÒIkğî¼Ëã‰éá>}…ÿ ºÿ ’ÿ âû.?ô®Ò¿Pëò÷ş	sÿ 'âû.?ô®Ò¿Qk»üŸüA”Sè¯@å¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¸Ê)ôPE>Šã(§Ñ@\eú(Œ¢ŸEq”Sè .2Š}ÆQO¢€¹çŸ´7ü‰Ÿö,júI-~Wî?íÿ &ÿ ñ3şÅOÿ I%¯ÃŠğs/'¥„øYöüçşNÄö,\é]¥~¢×å×üãşNÄö,\é]¥~¢í5Û€ş	†'øEMMz'(QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4 QFÓFÓ@m4m4çß´?ü›ÿ ÄÏû5?ı$–¿+÷öˆşÿ âoı‹Ÿş’K_‡àæ?OKğ³ìOø%¿üœˆ?ìX¸ÿ Ò»JıF¯ËŸø%¿üœˆ?ìX¸ÿ Ò»JıF®ÜğN|OñŠ(¯BÇ(QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹çß´Gü›ÿ Äßû5?ı$–¿«÷'öˆÿ “ø›ÿ bÆ§ÿ ¤’×áµxYÇÒÂü,ûş	mÿ 'âû.?ô®Ò¿Qö×åÏü×şNÄö,\é]¥~¤Wfø'>'øƒvÑ¶Ez79Fí£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\<ı¢ş1ÿ âoı‹Ÿş’K_†µû—ûDÉ¿üMÿ ±cSÿ ÒIkğÒ¼,Çâ‰éa~}‹ÿ µÿ “€ñı‹úWi_©ùoÿ µçöñı‹úWi_©;k³üŸüA(¥ÛFÚôQ(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mhÏ?hù7ÿ ‰¿ö,júI-~Wîoí¿ñßìXÔÿ ô’Zü2¯0ø¢zX_…ŸcÁ,ÿ äà¼Aÿ b½Çş•ÚWêU~Zÿ Á,ÿ äà¼Aÿ b½Çş•ÚWêUvà‚sâˆQEz(QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QEç¿´Wü›ïÄïûõOı$–¿k÷;öŠÿ “}øÿ b¾©ÿ ¤’×áxY‡ÅÒÂü,ûş	gÿ 'âûî?ô®Ò¿Rö×å¯üËşNÄö+Üé]¥~¥×nø'>'ø‚m£m-è¢m£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´P~ÑKÿ ûñ;şÅ}Sÿ I%¯Ã
ıĞı¢¿äß~'Ø¯ªé$µø_^añDô°¿>Çÿ ‚YÉÁxƒşÅ{ı+´¯Ôºü´ÿ ‚Xÿ ÉÁxƒşÅ{ı+´¯Ôİ¢»0?Á0Ämí¢¢»îsXmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xmí¢¢‹…†ÑNÚ(Ú(¸Xó¿Ú+şM÷âwıŠú§ş’K_…õû¥ûEÿ ûñ;şÅ}Sÿ I%¯ÂÚñ3Š'£…øYöGüÃşNÄö+Üé]¥~¦×å—üÃşNÄö+Üé]¥~¦×^ø&âQ]ç0QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QEçŸ´Wü›ïÄïûõOı$–¿k÷KöŠÿ “}øÿ b¾©ÿ ¤’×ámx™‡ÅĞÂü,û'ş	_ÿ 'âûî?ô®Ò¿S«òÇş	_ÿ '	âûî?ô®Ò¿Së³üGñ¢–Šï9„¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š óïşĞ_ü!­İhúß´/Tµ`³Ùİ_Æ’DH)<Îÿ †¤ø?ÿ E+Âÿ ø3‹ÿ Š¯øßÿ êĞ>2üC×<aÿ 	¥£ßêÒ,²CöXç‰"§*qòç¯zùgãüoâ/Ã:ãTğõÕ·tègÆ&†ô(ïäÛ¾ˆÌ}«‚¥\D.Ô.˜Â”¾Ö§è‡ü5'Áÿ ú)^ÿ Áœ_üUğÔŸÿ è¥x_ÿ qñUøw$oŒ¬§k+GPE6¸¿´'ü¨ßê±î~ãÃR|ÿ ¢•áüÅÿ ÅQÿ IğşŠW…ÿ ğgÿ _‡4Qı¡?åCú¬{Ÿ¸ßğÔŸÿ è¥x_ÿ qñUÒx+âç‚>$]\Ûx[Åz?ˆ.mI4:uäs<jN§8Ï¯ÁZë¾|Q×~øóKñ_‡®…Œ™1±>\ñŸ¿ƒº0àş`€j£˜Jë™hKÂ«hÏŞš+ƒø#ñ›Aøñğ÷Oñ^&!›÷w6ÀÉi8|Oî20{‚à×{^Ìd¤®švbQKE1	E-ã/xsáŞ’šŸ‰õ«Ny„u2Å‚Bî<d…cc\Wü5?Áïú)^ÿ Á”_ãOı¢¾éŸ´g€àğ¾«©İé6ğßG~·j¬ûÑ !†1‰å_%xş	;lmhå[‘Ò=KMÛr8+ÿ |šå©*Ñ~än¡m{ÎÇÖ_ğÔÿ ¿è¥xgÿ QğÔÿ ¿è¥xgÿ Q~S|pı>$ü[İsI]CAFûkKs5°ÉÀßÀhû@ÉÀ&¼R¸%«iFÇJÃÂJéŸ¸ŸğÔÿ ¿è¥xgÿ QğÔÿ ¿è¥xgÿ Q~×_ğËá'‹¾1kÿ ØŞĞîu«Õ¥ò€XáRqºI…AìEJÇT“²ˆş­«gì¯ü5?Áïú)^ÿ Á”_ã]w‚>&xOâU½ÔşñâmYRy4ë…™cf€ÅO€kóÇÁ_ğJŸjqE/Š<a¥h;¹h, {ÙØ’c\ı	Zû'öaı˜tßÙ“CÖtı?[º×Tš9æ–æ‹k"•B“Ç=É®úU+Iûñ²9§q^ë»=ªŠZ+°ÀJ)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J­ªjvz&™w¨êÙØYÂ÷0TŠ4RÌìO@ $ŸjµX¾5ğ¼7ğn½áË©¤·¶Ö,.4ùf‡Ñ%£,¹ã 1#4í 7ü57Áÿ ú)^ÿ Áœ_ãKÿ GğƒşŠ_…ÿ ğiÿ _4]Á(|"êßfñÖµv2ÚÂà~[kÌ¾"Á+ü]¢ÙÍuàÿ Xx™“•±½€ÙLÃÑ[s¡?ïÁ*˜˜ëÈt¨Rhû—ş‹áı¿ÿ àÒş*øj/„ôRü/ÿ ƒHøªüJñ'†µ_k×º.·a>—ªÙHb¸´¹BÄ~ ô ‚8¬Úäúüÿ ”ßêÑî~ãÿ ÃQ| ÿ ¢—áüCÿ ÅQÿ EğƒşŠ_…ÿ ğiÿ _‡RşĞŸò¡ıV=ÏÜøj/„ôRü/ÿ ƒHøª?á¨¾ÑKğ¿ş!ÿ â«ğâ¦¶³¸¼b¶ğI;HùQı¡?åBú¬{Ÿ¸?ğÔ_?è¥ø_ÿ ÿ ñT«ûP| cñ/Âß«ÿ Ù«ñ!|7«É÷t»Öú[¹ş”“xwU¶RÓi—‘(%íÜÔSúüÿ ”>­çíü´§ÂI:|Mğˆÿ {[¶ÍëFÇã—Ã}PgñÂ×dôò5«gşO_ƒôQı¡/åª®çô¦øHÖ1öRÊû=>Íp’gò&´kùéV*Á”AÈ"º½âç|,TèŞ2×ô­½§<CòV¢Ì;Ä—…ìÏŞº+ñëÀÿ ğP¯>hÒ[ø–Õÿ GÖ­LıdM’Å}eğ‡ş
{à¿Ok§øßH¹ğ…ì„!¿…¾Óc™b xÁÿ u€îÕÕe)é{z˜Ê„ãæ}§EUÒu{{M·ÔtËÛ}FÂáİZÊ²Å*ŸâVRAâ­×iÎ%´PQKE %´PQKE %´PQKE %´PQKE %´PQKE %´PQKE %´PQKE %´PQKE %´PQKE %´PQKE %´PQKE %´PWâ¯_|®O£x‡ÆÚªÀËg{z‘J”2’¤çAük½¯“¿hø'ş‘ñïâ©ãøKït=Rú(£0ı'…Lq¬jq¹[ {ÖUãÓWeÁE¿y»ÿ Oğ{şŠW†ğeøÑÿ Oğ{şŠW†ğeø×åoíû øßöqx¯5UƒXğİÄ‚(u«ù{È$$ˆ~hÛ õÊÌN@ğÚòå«Ë(Ù‘ÃÂJéŸ¸ŸğÔÿ ¿è¥xgÿ QğÔÿ ¿è¥xgÿ Q~ÑQı¡?åE}V=ÏÜOøjƒßôR¼3ÿ ƒ(¿Æ½Ã>'Ò<g¡Úë:£m«iWAŒ–’	"“k;XppÊGá_Ïİ~ÍşÁ¿òi¾ ÿ ®W_úY=uá±R­7Œ*ÑTãtÏ|¢–ŠôNQ(¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€<ïöŒÿ “|øŸÿ b¾©ÿ ¤’×á]~êşÑŸòoÿ ìWÕ?ô’Zü*¯0ø¢z_…Ÿdÿ Á+ÿ äá<Cÿ b½Çş•ÚWê}~Xÿ Á+äá<Cÿ b½Çş•ÚWê…v`¿‚aˆş Ú)ÔWqÌ6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E |)û~Çqx«M½ø™à;n¿nºÎŸl¿ñûÎª?åª[xs÷‡Íù¡_ĞÕ~(şÚ^Ò¼ûOxïJÑl£Óôèîa™-¡Dimâ•öÀ»±ÀàgŠñq´cŞG©èaê9{Œñ:(¢¼£´(¢Š öÿ Ù?ö–Ôÿ gˆQ^––ëÂÚƒ$:Æ¼ïŒ%AÓÌL’=A+Ær?f<;â-7Åš†³£ŞÅ¨éwĞ­ÅµÔ)"0È#ü:Šş~+ì¿Øö¶ÿ …[¯CğûÅWA|#ªÜ¡]ÌØ]6åÏrzDç¯ec»€\×§ƒÄ{7ìå³9+Òæ\Ësõ6ŠuîhÚ)ÔPh§Q@u>×V±¸²¾¶ŠòÎâ6ŠkyĞ<r!*Êx Æ¿¿m/vß~7_iZZy~Ôá¦›$C³)‡'®ÇFœíÛµûG_ ÁY<2²h|@‘€ğÜİØI °uĞ§–øúšàÆÁJ—7TtáäÔíÜüëÓì.5[ûk+Hš{«™VbN®ì@Uä+÷öxø¤|øc¥økM‚?¶ùk.§z£æ»º*7¹=p
;(¿%cı?~ÓŸm%]èš´wX÷„‡ë¯ÛšçÀAYÍšâ¤ôˆÚ)ÔW®p¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E |ÿ Nø[§Üx7ÃŸ-­£Uµ½]ÜÈ¸i ‘ãÜ{ìd ×C_›ûÿ ±[¯ÙSÄR‘Íµİ”£ñ¸Dÿ Ùëñú¾ªİu=L;¼Š(®¤+îoø%ü”¯ÿ Ø"/ı+ášûŸş	Cÿ %/Çö‹ÿ G
êÂÿ &5¿†ÏÓ)ÔWÒAÍø›áÏ…<iGâèúÚ0Á…ŒSÿ èJkç?Šÿ ğMÿ …¾<µš_Ãqà}X‚Rm=Ì¶Å»oÏOddé_XQYÎ:Ÿ¹qœ£³?¾>~Í~3ıµä³ñ%šË§\³-^Ğ–¶º°=UñÕê9¯*¯ßŠt?‹¾Õ|+âU¹ÓuŠ¾ş	Pöu8 ûzdWâÆO…z¯Á_‰ßƒµ‚²İéÒí[ˆÔ„&£‘sÙ”ƒÇ#µxXœ?±wÌôhÕöŠÏs‹¢Š+„é=“öuı©¼aû9ëÂ}"s©h37úf…u#}aıäÿ rz8PÃŠıvø'ñ«Ã¼kâo\…Ï—si)kI€¢„g ô ‚85øK^Áû.şĞº§ìëñ2ÏZäŸBºe·Õôõ?,ğ÷€şúd²Ÿ\Œk¿‰tŸ,¾šµ5u¹ûqESĞuÍ?Äú-¯¥]E¦_B—×Pœ¤±°X{jõ}Ñå¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPš~Ò¶ø‰ğ'ÇÌK!›KX7(m³Æ¦H˜{‡E¯Ã
şƒuÈ„Ú.¡o"şjkùò¯0Jñg¡…z4QEy'hWìßìÿ &›àúãuÿ ¥“×ã%~Ï~Áãş17ÀõÆëÿ Jç¯KüGèrb~{Õê+İ<Ñ´S¨ ÑN¢€(kzÕ†ô{í[SºËN±înne8H£E,Ì}€×ä¥ñ÷Zøçû`xcÄR]İ[iSø—O·±°óHHmVå!Pq’	fÿ iÚ½ëş
SûM	ş7†îşU+6¿qğOXír=8wÿ €ï
ø¿àOü—‡Ÿö1ißúSxØªüÕ8ìB>X¹>§îıê+Ù<ñ´S¨ ÑN¢€E:Š mê(´S¨ ÑN¢€E:Š mê(´S¨ ÑN¢€E:Š mê(´S¨ ÑN¢€E:Š mê(´S¨ ÑN¢€E:Š mê(´S¨ ÑN¢€E:Š mê(´S¨ ÑN¢€E:Š mê(´S¨ ÑN¢€E:Š mê(´S¨ ÑN¢€E:Š mê(´S¨ ÑN¢€E:Š mê(´S¨ ÑN¢€E:Š mê(´S¨ ÑN¢€E:Š mê(´S¨ ÑN¢€E:Š mê(´S¨ ÑN¢€E:Š mê(´S¨ :ı£?äŞş'ÿ Ø¯ªé$µøU_ºÿ ´gü›ßÄÿ ûõOı$–¿
+ÅÌ>(†ágÙ_ğJßù8Oÿ Ø¯qÿ ¥v•ú¡_•ÿ ğJßù8Oÿ Ø¯qÿ ¥v•ú£´×VøFŒJ)vš6šî9„¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¯ÆoÛëşNÛÇßïÙÿ é½~Ím5øËû}ÉÛxûıû?ı"·¯;ü5êuá¾7èx+Ø]GcëÛL¶sHğÇpÈDo".Šİ(’2@ä\õ~ƒşÅ¿t?Úö5ñw†5ŸôyÏŠ.fÓõ$@ÒY\;M²‘¸s†\ÊHÈ8#áï‰?uß„ş6Õ<+â;O²jÚ|¾\Š2REê²!ÀÜŒ`}n•åN“„c>ŒîŒÔ›Ts4QE`hQE ~›ÿ Á=k…ñ¦“mğÇÅ÷Å¼Ee=ôïÍíº®|–'¬‘€p‰ª’ßq×óá£ë¾Õ¬õ=6ê[FÎd¸·ºŠÉŠC+© ııj?Ú;ÀíÒCoãM)V=VÉ0¾`è·1î?p>ëdtÚO·„Äs¯g-Ï:½.WÍ (¥ÚhÚkÒ8Ä¢—i£i ¯Œàª–Áşøn~ñø–üÖçÿ ‰ö~Ó_ÁT²¿³Öïâ‹qÿ ’—uÏˆşi|høËş	ÿ n.?kÑóşV7ıq_²µøáÿ ÷ÿ “·ğOû—ÿ úC=~Èm5Íşõÿ #lOÆ„¢—i£i¯DäŠ]¦¦€Š]¦¦€Š]¦¦€Š]¦¦€Š]¦¦€Š]¦¦€Š]¦¦€Š]¦¦€Š]¦¦€Š]¦¦€Š]¦¦€Š]¦¦€Š]¦¦€>sÿ ‚„(oÙ'ÆÇÑìşOA_û!ÿ ñˆş7ÿ zÃÿ K ¯ÆúğñßÅ^‡¥†ø¨QEça_sÿ Á(ä¥øãşÁèá_WÜÿ ğJù)8ÿ °Dú8WVøÑ1­ü6~™ÑK´Ñ´×Ñ@”Rí4m4 •ùñÿ Xøgl,|!ãûxv]yÍ¢ŞH£‡R­,9÷fˆô¯Ğ¦¾pÿ ‚…h1k?²‹f‘7K§Ëgy	şë˜ĞŸûâGüëŸjRFÔ›ŒÑøçEWÍ¸QEúÿ Àø¼ş*ø]«xúV{ÏN%µ,s›I‹0Qşì‚OÁÔv¯µ+òGş	§âÉtÚfÓKVıÎ»¦İYºçŒ¢} ¨òHüM~·í5ô8IóÒWé¡åW,Ø”Rí4m5Øs‰E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@	E.ÓFÓ@5_ùŞ×ÿ ĞM=õıêª²ï?ë‹ÿ è&¿êò1ÿ gçúø^¡EW’w~ÏşÂòi¿ÿ ë…Ïş•Í_ŒûAû¯üboÃÿ úásÿ ¥sW£ş#ô91?=êŠ]¦¦½ÃÍŠ]¦¦€¾~ı±¿j/ÙÓÀ¶3C?uDhô«6Ãy}šá×û‰Ø¼Ø7Ô~Ò?´g‡ÿ gI¬j¬·zµÀhôÍ%	.åÿ EÈ,İ‡’üfø™ñ+_ø¹ãMKÅ>%½kİRù÷1è‘ û±Æ¿ÂŠ8óÉ$×+ì×,wüª4¹ß3ØçõBçVÔ.o¯g’êòæVšiåbÏ#±%™‰êI$“ï]ÀŸù.?ìbÓ¿ô¦:á«¹øÿ %ÃáçıŒZwş”Ç^$~$zOc÷†Š]¦¦¾¨ğÄ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ¢—i£i ;ı£?äŞş'ÿ Ø¯ªé$µøQ_»´`?ğÏÿ ìWÕ?ô’Zü'¯ñDô0¿>Ëÿ ‚Uÿ ÉÂx‡şÅkı+´¯Õ*ü­ÿ ‚Uÿ ÉÂx‡şÅkı+´¯Õ*ëÁÃñ…Q]Ç0QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE Wã'íõÿ 'oãï÷ìÿ ôŠŞ¿fëñ“öûÿ “¸ñ÷ûöúCo^v;økÔëÃ|oĞûCş	Wÿ &÷âûn?ô’Ò»¿ÛOöS¶ı¢¼/ôˆ!‡Çz<Ltë‚B}®>XÚÈÇŒ’…‰ä|ğŸğJ¿ù7¿ÿ ØÓqÿ ¤–•ö]mJ
¥mb*IÆ«hşzoì.´«ë›+Ûi¬ï-¤hg·¸B’E"’O*À‚<‚*
ı/ÿ ‚†şÈğ”Ø]üSğ}ŸüNm#İ®iğ§ü}Â£ş>Pùh€|Ãø”g‚‡æ…x•©:2åg£Nj¢º
(¢°4
ì~üV×~
ü@Ò¼[áéü»ë>h\Ÿ*æ#ÃÃ U‡£‚0@#¢šn.èM_F~ñüøÅ¡|vøw¦ø»@“ıäÕØ-'\o…ÿ ÚRGÔaÃ
î«ñ[öGı¦µÙ¿â$wr™®¼%©2Ã¬iñœ’ŸÃ4c§˜™$xn\Œä~Íx^ÓüU¡ØkEÜz†—}
Ü[]Br’ÆÃ*Ãğ5ôXzê´|ÑåU¦é¿"ıQ]F_ÁU?äŞü=ÿ cM¿ş’]×Ùuñ§üSşMïÃßö4Ûÿ é%İsâ?…#Z_>8ÿ ‚{ÉÜx#ıÛÿ ı!¿d«ñ·ş	ïÿ 'sà÷oÿ ô†â¿d«›ü7ëşFØŸŒ(¢ŠôN@¢Š( ¢Š( ¢¾aı¤¿oOü	»ºĞtÈ‹<]*öV²…·µoI¥çğ('×oùéñ[öÖø·ñjY’ûÄóèšc“3B&Ò?ºJî?ßc\uqTé»nÎˆQœõØıƒñ?Åx%™<Câ½Cu)¨ê0Àß“°5Á]şØ¬›|GĞØç¹˜È?ñĞkñIiIØ–fc’IêI®‹Âÿ ü_ãxÌğ¦·¯Æ§i}/Nšäé”S\^›~ìN«EnÏÚ}+ö¥øA­H‰mñ'Ã;ÜáV}J(I>Ÿ9×£ézÅ†¹j·Zuí½ı³}Ù­eYıIøâOø‹Á“$> Ğ5M
Wû±êVr[±úQQxgÅšßƒ5Hõ-W¾Ñoã9[>ááÀ”ƒMc¤Ÿ½<2èÏè"Šü¬ø/ÿ 4ñçƒ'†ÏÇ6Ñx×HV¸
¶÷ÑPê?Ñ—'xWègÁÚÀÿ tW¿ğ®·RÂ¹Ó®—wmdyéşĞÊÄ×},E:¿Ôå)ÓÜôz(¢º‚Š( Š( Š( Š( Š( Š( ?à¡?òhş7ÿ zÃÿ K ¯ÆÚı’ÿ ‚„ÿ É£øßıëı.‚¿kÃÇz–à~¡EWœu…}Ñÿ Ÿÿ ’™ãûGÿ £…|/_tÁ'ÿ ä¦xãşÁÿ èá]X_ãDÆ·ğÙúkEWÑ@QE W‰~ÚËŸ²ÏÄA/İşÏR?ŞÆWõÅ{m|ßÿ ××Cı”üWí’ê2ÚYEîMÄnÃşøG¬«;S—£.ŸÆÇ*(¢¾\ö‚Š( vı†e–Ú»áëC÷ÍÔÊİ6Ò†ı	¯Úzü„ÿ ‚nøRO~Ô:Mğ|:%…İü‡Æ`_üzqùWëİ{¸û·êy¸ŸQ^ÈQE QTõfÃÃº]Ö§ª^ÛéÚuªgººGJ:³18ë@*;›¨líäâT‚Æç’F
ª=I=|û@ÁO¡Óî.ô_…z|w¯Ïˆµ$>QÇxaà°ôgÇOºG5ğ·ÄOŒ¾8ø±|×^.ñF¥®;ˆ®&"ÿ r%Â'üEpTÆBGSªyKW¡ûC®~Òß	ü7!Pø‹á¨¥cMN)}U‘X¶¿¶7Á[É¼¤ø¢+g–Vï¦ WãßÃï‚¾;ø¬d>ğ¦©®ÅÃÜZÛŸ%¡áAö'5sâìûñá]Û|WàİWF°Ş#7’Ãº Ç 2.TÛkŸë•mÌ£¡¯Õá·1û…áŸxwÆ¶ÆãÃÚö™®ÛLºmäw
>¥ÖÕ=ºV± ßÅ}¦ŞÜi×°ÑÜÚÊÑH‡Ô2Gá_^|ÿ ‚“xÛÀw6ÚoCxÓ@ÈVºl.£ú‡àKô~O÷ÅkO;MX‰a¤µ‹¹ú­Esş%xsâ×„¬üKá]N=SIºGÃFãïFêyGåO<ÄWO^’i«£“U QE(¢™q4v°I4Ò,0Æ¥ŞIUU$’z (õŸ®x‹JğÍ™¼Ö5;=*ÑzÜ_\$1ø+àÚsş
U5½å×‡>öÆÆ9¼OqpÄp~Íc?òÑÁ²ôcğ_‹¼qâjÏ©ø“[¿×oÜœÜj3ö‰Àö
óêãahêuÃ)k-Ú}Köµø5¤¹Yş$øuÈÿ ŸkÕœ~qî©4oÚ»àî½2Eiñ#Ã¢GáVæõmóíûÍ¼×áä0Éq2Em,®Bª %˜€ÔÕíKÃz¾¾¡¥^Ø£ræİã>ä
åúôÿ ”ÛêÑî@Zn©e­YÇw§Ş[ßÚIÊOm*É}I­Wà?ş%x¯á¦¦º‡…¼A¨h7JÛ‹Y\2+û:çkfWÙ¿
ÿ àªZî§}Çş‹Ä"á5-*Eµ•ÿ ë¤d$ú®Ğ?ºk¦6ÒzË%ğê~øƒÄo…4KİcX¾‡MÒì¢i®.î,q êI5‰ğÿ â§„¾*Øİ^xG^³×ímdÏ-›¹
xëŠü±ı­?nMSö‹Ó-ü7£é’øoÂqÈ&	'M} ÁS& Tògœ1< >ÿ ‚L^yø‹iŸõW¶rãıèåû%\qJ¥UìL¨¸Ã™î}ñEWqÌWÆ_¶÷ˆ?h/øãJ¹øRšãødi*nÿ ²,cºäK.âT£>vy}+â»ßÛcãæ›u-­ßµ+[˜›lÍen‡ĞƒA®*˜¨Ó—,“:#FSWLı¡¢¿?á¹¾:ÑB½ÿ À[oş5Gü77ÇOú(W¿ømÿ Æ«?¯SìËú´ûŸ´ôWâÇü77ÇOú(W¿ømÿ Æ¨ÿ †æøéÿ E
÷ÿ m¿øÕ^§Ù‡Õ§Üı§¢¼ßöoñF©ão€ş×u»¶¿Õ¯ô¸g¹¹eU29±
 ü…zEz|É3•«;ñí½âÚÃ;Ó.~¦¸şJ¿ì‹î€¹Ë¸•(Ï_AÖ¾,¼ı¶>>i×RÛ]øëRµ¹‰¶ÉÖVèè}1dã©Š9rÉ3xÑ”ÕÓ?hh¯ÅønoŸôP¯ğÛÿ Qÿ ÍñÓşŠïşÛñªÏëÔû2ş­>çí=ø±ÿ ÍñÓşŠïşÛñª?á¹¾:ÑB½ÿ À[oş5G×©öaõi÷?iè¯:ıœüMªxÓàOuİjé¯µmCI‚{›–US$…ybÏ°¯E¯B/™&rµg`¢Š)ˆ«ªÿ È.óş¸¿ş‚kùì¯èOUÿ ]çıqı×óÙ^F?ìüÿ C¿Ô(¢ŠòNà¯Ú/ØGşM7áÿ ıp¹ÿ Ò¹«ñv¿hÿ aÿ ğûş½î?ô®jôp?Ä~‡&'àG½QEîhW~Ó´ÿ †fß	›ÍE×PñÒìİ9 –áºooîFW#Ødñ^gûW~Ş^ø'ï†ü(Ğø‡Ç!LlîµÓ›ÖVæqÿ <Çü‡ò·Æ^5×~!xó^ñ&©q¬j÷mºk«–Ë@EP8
  p ¯?ŠTıØnuÒ¢å¬¶5ş,|[ñ7Æ¯]ø›ÅWí{7Ê‘®V+xÁ%b‰…zu<’I$×Eá¶äîÏE+h‚»Ÿ€ÿ ò\>ØÅ§éLuÃUÍX¼ğî³aªéÓ›]BÆâ;«iÔc•20p@<úQf˜=èFŠüYÿ †èøëÿ Eóÿ -¿øÕğİè Şà%·ÿ ¯oëÔû3Îú´ûŸ´ÔWâÏü7GÇ_ú(7Ÿø	mÿ Æ«ô÷ö4ñÖ½ñ+ömğ‡‰<K¨>«­Şı³ír"£I²òx×… p¨£Ú¶£‰irÅÔ£*jìöš(¢ºÌŠ( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( ;ı£äŞş(Ø­ªé$µøM_»_´oü›ßÅûµOı$–¿	kÆÇüQ=/ÂÏ³?à•?òpŞ!ÿ ±Zãÿ Jí+õR¿*ÿ à•?òpŞ!ÿ ±Zãÿ Jí+õRºğ_Â1Ä?|(¢Šî9®QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .øÃû~ÉÜxÿ ıû?ı!·¯Ùêüaı¿?äî<şıŸşÛ×şõ:ğß>Òÿ ‚Tÿ É¼ø‡şÆ›ı$´¯³+ã?ø%Oü›Ïˆìi¸ÿ ÒKJû2ºpÿ Â‰Wï°¯Ê¯ÛûöC?	õÉ¼áÔçÿ L´·O“K¹sÙGİ…ÏİÇÊ¬vp
ú«YŞ!ğş›âÍÿ FÖ,âÔ4»èZŞæÖa”–6*
u¨ªÑåb§QÓw?Ÿ
+İkÙQı›~!=¤^uï„õ"Óhú„ƒ’ƒïC!y‰	ş!µ°2TxU|Ü¢á'nzÑ’’º
(¢¤ ¯µ?àŸ?µ×ü+MrßáÇ‹¯¼¿êSìëË†ù4Û—9ÚIû±HÇÊí»€ÎÕñ]¥:’§%(‘(©«3ú ¢¾"ÿ ‚zş×âmğÏÅ·a¼O§@F•{3|ú…²Lm²Æ£¯VEÉåYÛµô´êF¬T¢y3‹ƒ³
øÏş
­ÿ &óáïûmÿ ô’î¾Ì¯Œÿ àªßòo>ÿ ±¦ßÿ I.ë<Gğ¤U'ï£ãoø'·ü×¿İ¿ÿ ÒŠı•¯ÆŸø'¿ü×~—ÿ úAq_²ÕÍşõÿ #\OÆQEz'-ÂŠ( .ùóûuşÜ×=æ£ğÛáÕÿ “u`Ö5Ûgù¢noŒ:3AÊŒHö_ÛËöœ€ÿ “GĞ®|¯ø+IüÖp$¸önv§ûDŸà"¿ ]ÚGgv,ìrYI>µåâñ?»†ıNÚù½é3HÅ˜–f9,NI4”S¡…î&H¢F’G`ªŠ2X“€ ¯ï>Öı€ÿ c};âÜrøÿ Æö¿jğÅ´Æ;Lc„¾™OÎòc¬j~]¿ÄÙÏ
C~ ØØÛivpÚYÛÅiiâ‚	j8
ª8 z
å>x…¿
ü+á8#TVŸ¼›z<Ûs+ıZBíõjì«éhRT —SÈ©QÎW1üYáÇZÖ‹â.×XÒ®—l¶·‘‡F÷ç¡ˆäAøëûe~Ìïû7üHKm<Íqá=]ãJ¸˜å—iH÷d,¼÷VS×5ûE_7Á@¾Ãñöm×î„Aµ•Öm¤Û’yú™ÏÕWÒ³ÅRU ßTU2·F~8V×ƒ|m¯|=ñ®»á½VëFÕíNbºµ}¬=AìÊ{©È#‚bÑ_?¶¨õOØOØçöÊÒÿ h}4Mi¡Ó<}ggµ_–;ä™¡¿ŞN£¨Èéôå>¾ñf«à_i Ñ/$°Õ´éÖâÚâ3‚®§õ¡‚	ƒ_¶ß³OÇ?öˆøW§ø¦Ö4´¿ÛjV(Ùû5Ê¹G}¤ËŸáas^îí,·G›ZŸ'¼¶=RŠ(®ó–áEP
(¢€¸QEÂŠ( .QEsç?ø(Wüš7?Ş°ÿ Òëzük¯ÙOø(Wüš7?Ş°ÿ Òëzük¯üUèzXoú…Q^qÔ÷OüwşJg?ìşğµ}Óÿ ÿ ’™ãûEÿ £…uaßÃgé½Q_Fy7
(¢€¸WçüoâŒ.¾ø{ipUwÖoãVûœàÜæcƒşÉï_gürøßá¿€^»ñ7ˆî U;KØy×³à•Š1ê{Š2OJüIø¡ñ#Yø¹ãígÅºôŞv¥©Îe`¿v%	ú*¨UÃkÎÆUQ‡³[³¯'Ìö9j(¢¼3Ñ
(¯jı•fÍ_öø‰Imá»Iµ}KlQgıZò>QÛ–è¦ª1sj1ÜRj*ìûgş	wğv_
ü9Ö|y¨ÛùW~#•`±Şai °ô!o¨Ozûv¨èz‡†t['K´ÇM±-­­aH£E
ª  UêújTÕ8(£ÇœùäØQE©
(¤f
¥˜àI4ÌxïCøiá=GÄ¾#¿MÑì#ó&¸“ò
£«1$ £’H¿ ?jÛÄß´vµ-š<º?‚mæİe££`É»,ä}÷ïº¹Àç,zOÛ³öª¸øåã©|7¡^çÀš$åmü£òßNÖ¸oU2§û9?ÄkåšğñX‡7ÉJ.UÍ-Â½¯öGıæı¢¾-Úh³ù‘xvÅ>Û«ÜFpË  ÔöwbzÍÎÚñJızÿ ‚rü-‹À?³Å³-ºÇªø¢fÔg˜Â	HûmÇıu5†—µ¨“ÙVŸ$n¤ü7á­+Áú–‹¢iöú^•eŠŞÒÕGÀÌ¤’O5rúÆÛT³ÎòŞ+»IĞÇ, xäR0U”ğAOE}‘äÜüoı¸¿f”ıŸ>',Ú<L¾×ƒÜé£’-˜æ[–=v’ÿ e‡R	¯›ëöş
#à¼mû2k—‚2ûÃóÃªÛ°€G/>\qê£Ò¿+çqTÕ:–[3Õ£>xjz¿ìëûFø›ösñ¤z¾+]is•]GG‘È†î?ı•Ç;_àh>|IĞ¾.xKñW‡.ÅŞ—¨D¿thÜvu9zNkğ"¾ºÿ ‚tşÑü/ø¡‚µK€¾ñDËù‡m|FØ¤ğ#>¹Cü5®»„¹%³"µ>eÌ·?Y¨¢Š÷O6á_ğR_ÚšH^o„^¹hÎM~òÁ €ÉjìAÿ U_ï
û£âw-şü:ñ'Š®‚¼:>Ÿ5ç–Çc"©õfÂş5ø3âj0ñ§®ê³µÖ§©\Éws3uy‹1üÍyØÊ®ä]N¼<9Ÿ3èfVÿ Ãÿ ê_<m¢x_GŒI©j×Ikïº¥7¢¨Ë`k½·ö,ñv›àÚ{Àz®­$PXıªKVšc…§‚HQ‰ìH¼·Zñ ”¤“Øï“i6Õ¯€³7‚¿g¿ÛÙè:l3ëX[Írâ n®_1İÉDÏDS€=NIõK›ho-ä‚â$	×E¬==EKE}Db¢¬–‡ŒäÛ»>føïû|5ø½euu¤iĞø/Ä¬¤Å¤Ä#ß·› Â0õ+µ¹ÎOJü·øÕğ;Å üa'‡¼UcäJA’ÚòZŞî<ıøŸ#Ôx Wïy·Çï>ı ¾Şxk\…Rl,5Pe²¸ÇË"ŸNÌ¿Ä2=ã¯…ExèÎŠuÜ]ÇáE~‡ÿ Á#æçâ¤E¹ÿ ‰[ÿ À°Oò¯ƒü}à}_á¯ŒµëÖßeÕô¹ÚŞâ0r¹†SİXÀ÷ûcş	/yåøÃâ¦ÖØZK÷dìõåá}Úêçemi¶Òš(¢¾ˆò®ò¯üàN‡ñà¹âÈôè"ñO‡a±jF²@¤y±;¼»2Ã9Á^1“Ÿª«Šøİbš§ÁYÉÊ\hñ7Ñ­äÖ³©84Ë„¹d™ø)EWËÈQEû‡û#ÿ É³|6ÿ °,úzíy'ì“ÿ &Ïğ×şÀ–ÿ úzİ}M?‚>‡‹'ï0¯“à¢ô_|ÕücŸ^)ğê¥Ú_EÍnVX¤a÷”).3œã šúÊ¸_ÚlzÇÁˆ6R’ãÃ÷ñŸlÛ¿?…b§˜á.Y&Áš(¢¾XöBŠ( ÜoÙ/şM£á¯ı€í¿ô
õªòoÙ;şM§á¯ı€­ô^³_UOà^‡‹'ï0¢Š*É¹SVÿ ]çıqı×óÕ_Ğ¶«ÿ  »Ïúâÿ ú	¯ç¦¼ŒÙùş‡v¨QEäÁ_´¿°—üšoÃïú÷¸ÿ Ò©«ñj¿i¿a?ù4ï‡ßõíqÿ ¥SW¥ş#ô91?
=ê¼ëö‡ĞüQâ_‚¾-Ò¼n‰ï,ü›²İ-´›Ë.q!eòîç"½Šö¤¹“G¥gsñµ¿à¿‰>
Œ“É'X²ÿ ãÕÅüVı•>'|ğÜ:÷Œ¼:ºN•5ÊÙ¤ëm>eef¶9ôFçâ¿q«ã¯ø*oü›‘ÿ c%·ş“Ü×•WNœ“zÄJRI£ò’Š(¯ ï
¿áı÷Åö›£i°ı£QÔnc´¶„°]òÈárH,@É8ªİ|çãŸÃ¡ÿ Sÿ ¥QÕE]¤'±êğï?ô$Çÿ ƒ‹ş=Gü;Ïãßı	1ÿ àâÇÿ WìíıFŸwı|;ë3ì~6Ã¼ş=ÿ Ğ“ş,øõ~—~Ç¿5ÿ …?³Ÿ„¼-â{!§k–kûE°™%	¾òiæF*r§ƒßÖ½–ŠÚ–eÍg:Ò¨¬ÂŠ(®£…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEsÎÿ hïù7ŸŠö+jŸúI-~×îçíÿ &óñCşÅmSÿ I%¯Â:ñ±ÿ OC³>Ìÿ ‚Tÿ ÉÃx‡şÅkı+´¯ÕJü«ÿ ‚Tÿ ÉÃx‡şÅkı+´¯ÕJëÁÿ Çñ…Q]Ç0QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE Wãíùÿ 'qãÿ ÷ìÿ ô†Ş¿g«ñ‡öüÿ “¸ñÿ ûö_úCo^v;økÔêÃülûKş	Uÿ &óâûn?ô’Ò¾Ì¯¿à•«ÿ ï®ûøçÿ Im+ìšéÃÿ 
&U~6QEtGÆo„:Ç/‡ºŸ„|EûKµİÂ¨óm'\ìš3Ù”ş`²œ† ş&üføAâ¿u/	x.òÔï†áò® lìš3İXîe8*@ıé¯	ı®?f7ö’ø~öÑˆìü]¦£Ë£ê07‘“‡şy¾ 'øNg[‡CÚÇš;£¢NGg±ø«E_ñ¨øW\¾Ñõ{9tıRÆf·¹µv¼R)ÁR>µB¼Ó
(¢€.hºÕ÷‡5‹-WK»–ÃQ²™.-® b¯ŠAVSØ‚¯ÙŸØÿ ö °ı¤~¤·;Oéj°êö+Àfè·ùæøÎ?…·/ oÅªí¾|^×ş|BÓ<[áÙü»ËFÛ5»å]ÀØó!wVêVeuaëº2òf5iûEæ~õWÆğUoù7ŸØÓoÿ ¤—uô·Á¿‹šÇ‡Úg‹|9qæÙ]®Ù`cûÛY†<Èd™OàA2¬	ù§ş
­ÿ &óáïûmÿ ô’î½œCR£&
zTIŸÿ Á=ÿ äîü	ô¿ÿ ÒŠı–¯Æø'Çüç€ş—ÿ úAq_²õ†øo×ü1
(¢½”*®©ªZèš]æ£}:ÚØÙÂ÷Ï'İ4RÌÇØ OáV«å¿ø(ßÄÙ>şÎwÚm¬¾]ÿ ‰®SJ\˜BA’cô*›ıt¨©5N.O¡Q3HüÌı¢>1^üvø¹¯x¶äºÛ\KäØ[¹ÿ Qj„ˆ“ëO«3õæÔQ_/)96Ùì%edŞ|Ğ‡‰¾8ü?Ò™wGu¯XÇ ÿ cÏMßøîkƒ¯pıˆôáª~Õ_á#;oû÷’ìµTÕæ—˜¥¤[?l(¢Šú“Æ
ÄñÆ†¾(ğWˆ4g]é¨é÷Œ¸ÎD‘²cõ­º(ß@?ú+CÄVÙş Ôíq"êX±şë‘ı+>¾HöÂ¾¨ÿ ‚w|p›áÇ_^\lğÿ ŠÊØLß*\óöy=K/é'°¯•ê[;¹¬.¡º¶• ¸…ÖHåC†F ƒØ‚+Jstä¤º(ó&™ıÑ\WÁ_¯Å/„¾ñX*Òjºl7í6Ğ%_ÁÃÂ»Zú„Ô•Ñã½‚Š(¦ ¢Š( ¢Š( ¢Š( ¢Š(ç?ø(Wüš7?Ş°ÿ Òëzük¯ÙOø(Wüš/?Ş°ÿ Òëzük¯üUèz8…Q^qÔöïü¿ZÓô_‰5}BúÚÅ$ÒcTk™–0ÇÎÄd×ÄTV´çìæ§Ø‰Çš.'ô'|3åüE¤ õkè‡şÍYš‡Æ_ i1—¾ñÏ†ìÓû×½¼cõzü¢½¯¿å9¾¬»Ÿ¶Ş"ı³ş	ø]XİüEÒ.1ÛNg½'şü«×Ïß¿àª^Ò­g¶ğ‡/uÛÿ »î¬µªŸïlÈãØìú×æU”±µ%¶…¬<çoñoã?‹¾8x¡õßêÒj7C+#ä‚Ù	ÎÈ£*ô÷8É$ó\EWnNìèJÚ ¢¾ı•ÿ bÏÚSC—_3Òô}"ÚàÛ\ÛCÜ_DàdŒìU*ÛÛ¿Wßßÿ a?„ÿ Z¤Ğÿ á&ÖcÁşÑ×¶ÜoT‹5Áèvîs]t°µ*+ìŒgZ0Ó©ùïû6şÃ>8øíwi©jÓx[Á­‡}Vò,Ip™éoåÉşùùG©èWşü*ğ×ÁŸÙøcÂºzØi–ÿ 1$î–yŞ’Fş'8äû 0 ®¢½Š8xQZnpÔªênQEt˜…Q@|¿ÿ 
øİ'Â_“éz|ÆwÅ,úe»)ÃGÜÜH>ˆBdt2ƒÚ¾ ¯Èø)Ä§ñ¿í{£Å7™§xfÖ=>%Sòù¬¢Y›ë¹Âúæ+“SÙÒvİ›Ñ4Ï•è¢ŠùÓÔ%µµ–öêxÉ4Î#DY‰À@¾ğü>ğ¦‹¡Û…i–PÙF`m5AÁkğ«à}ˆÕ>4øÈ®ásâ>¾»®cık÷®½|ÒLáÄ½QEëGñóIçÀßˆZy7¿}˜Û¾ÓøWàå@ş4´ûwƒõÛb8šÂxÏãë_ÏÅxøıâÎü6Ì)ğÍ%¼É,NÑKIÊÀäGCL¢¼£°ıÚı¾&Œ?<#âÂÁ®o¬”]ãµÄdÇ7áæ#cØŠôjøƒş	SãfÕşx«ÃÉ½ômQnc\ò±\'é¾ü×ÛõôôgÏN2<Š‘å“GËŸğRMr]#ö[ÕàŠVûGP³´}§—Ìóÿ ÈUù_®¿ğRı!õ/Ù~öáP°Óõ[;– }ĞY¢Éüeñ¯Èªò1¿ÅùØ€(¢ŠóÎ“ô§ö/ı¿,õ«]+À?®ÖÏTV×NñÍˆ®€ÀHîıÙ;<7|7-÷½;õõßì«ÿ  ñÁÏ±økÆ_hñ7ƒU‚G3>ëÍ=:~ì“ûÄÜcÀû¤cÕÃã-îTûÎ*”:Àıe¢¹ÿ øÿ Ãß¼3kâêÖúÎ‘r>K‹vÈuV•aU€#¸®‚½„ÓÕGÁÿ ğSï€k®x^Ãâ†“j>ß¤í³Õ¼µæKflG)Àä£¤ú8ìµå?ğJ;Ï/ãO‹-3ş·Ãí.?İ¹„ìõúgâØxÃÃš…ªÛ‹­7R¶’Òæz7R¬?#_şÆ?±—g?:ş¿­İéWšÆ•q¦ÚËi;´ÎZâFd(åb9äàšóêQj¼jEiÔéEìÜYö½Q^ÊÌ|Røgâàzd^è—®¹ŠòLü]ÿ `‹¿ıô¥³Üü¢Š+äÏh(¢Š ıÆı’ÿ äÚ>ÿ ØÛÿ A¯Z¯%ı’ÿ äÚ>ÿ ØÛÿ A¯Z¯©§ğGĞñ¥ñ0®Wâ°İğ·Æ ôşÆ¼ÿ Ñ]Urß¿ä—øÃşÀ×Ÿú!êå³Üü¢Š+äÏh(¢Š ıÈı“ÿ äÚ¾Ø
×ÿ @ëå²ü›WÃOûZÿ è±^¯_UOà^‡/‰…QVIWUÿ ]çıqı×óÓ_Ğ¶«ÿ  »Ïúâÿ ú	¯ç¦¼ŒÙùş‡v¨QEä¡_´ÿ °§üšÃïúö¸ÿ Ò©«ñb¿j?a_ù4ÿ ‡¿õí?ş•M^–øĞåÄ|(÷Š(¢½³Î
øïş
™ÿ &ç¥ØÉmÿ ¤÷5ö%|{ÿ Jñ:gıŒvßú"â¹ñÂ‘­/”4QE|Ñëw ÿ äº|9ÿ ±“Mÿ Ò¨ë„®ïàü—O‡?ö2i¿úU\>$)lÏŞZ(¢¾¨ñBŠ( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( ;ı£¿äŞ~(Ø­ªé$µøG_»Ÿ´wü›ÏÅûµOı$–¿ëÆÇüQ;ğÛ3ìßø%Güœ7ˆìV¸ÿ Ò»JıU¯Ê¯ø%Güœ7ˆìV¸ÿ Ò»JıU®œğŒ1QEÜs…Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@Q@~0~ßßòw>?ÿ ~Ëÿ Hmëö~¿?oïù;Ÿÿ ¿eÿ ¤6õçc¿†½N¬?ÆÏ¶?à•Ã³®µïâkŸı&µ¯±ëã¯ø%ˆÇìç«{ø’çÿ Ií«ìZéÃÿ 
&U~6QEtQ@~ßß²	ø·¢Éãÿ Ùñ–—/m _ŸT¶A ûÓ û½ÙFŞHA_•ıWæoü;öA¾¸ø£àÛ4[¹k–/“1ÿ „¤nOÌ?…ˆ=åbğÿ òò?3¶…O°Ïƒè¢ŠòNà¢Š(ßcßÚŠÿ ömø‚²\™.üª²Ã«Ø®IQÑn#óÑ3Óø—rğv²ı›ÿ ?Ö,¼Aû1øGTÓn¢¾Ó¯|GiqmuIc{+¶WSÜAZüµ®şããN»ğ>…÷¯ö½ÏYY°yî´aé$J?¸í>üq†yŞqÕ
Î4åMìc*iÉM•ÿ ÷ÿ “½ğı¿ÿ é¾æ¿fkñ—ş	ñÿ '}à/ûÿ Ó}Í~ÍW£şõÿ #“ñ ¢Š+Ğ9B¿2¿à¬^,7Ÿ<á¥”°Òå¿dÓË°gßÿ ¯½~š×ãßüƒV}Kö¬×íØüº}••²ı/ó”×1Ú—©Ó‡^ùóQ^é}ÿ ñµûOíoà·#"ï¤ÿ É)×ùµ|ã_RÁ6!şÔúCş¯O½qÿ ~ˆşµµâGÔÎ§ÀÏ×ê(¢¾”òŠ( çÿ â2ˆş!x¡GEÕ.€ÿ ¿Í\ít¿ä¢ø«şÂ·_ú9«œ¯•{³Ú[QHgëGüÅO¯~Íòé’ÊY´MfæÖ4'îÆê“Œ{–OÈ××Uù÷ÿ ‘ÔMâm‰bR>p½u¸Ròü«ô¾‹ù©EUei°¢Š+¤Ä(¢Š (¢Š (¢Š (¢Š ùÏş
ÿ &‹ã÷¬?ôºŞ¿kö[ş
ÿ &‹ã÷¬?ôºŞ¿kÄÆÿ zàaEWu}Ëÿ ¥±·½ø‘ãqsªé1`J€ığõ¯†«îßø$Øÿ ‹ã³ÿ P¨ôutá¿‹ß?I&ğŞ“q‘K¥ÙKŠQã{t*ÊF Aøéûj~Ì³şÎÿ İ´èdë,×\ä!9ËÛ1şòdc=T©äîÇìåyÿ Çoƒ:7Ç¯†š§„u¥Ø—Íµ»Q—´¹P|¹WèIwVaŞ½ŠôUXÙnpR©É/#ğvŠè¾!øZø[ãM[ÂŞ!µ6z¾›1†hú©ã*ê{«)pA®v¾}¦™ê…QHfı”h›ÿ ÙÇâ•®´¦[^íµÖlPçÍ·-÷Ôtó#9eéüK×ív‡®Xx›E±Õô«¸ï´ÛèæÚê”–7PÊÀúE=u÷Ïügö¥şÃÔ£øKâ{Ìi÷²ğıÄ­Ä3±,öÄŸáÉÓçÜ9.1èá+r¿g-™É^Ÿ2æGéMQ^Ñç…Q@Q@šT‚'’F	)fcĞÔ×óÿ ñÅ2xãÇŞ$ñ¥‹êÚ•Åñİ×÷’³ãÿ ¯İOŒ:³h?üo©§g¡ß\®=RİØ*ü¯'şwa–ì(¢ŠòÓ¸øâ=;Áÿ ¼®êó‹]+M×l®îç(Î"‰'FwÂ‚N ' ÇúÑÿ ÷ğş‡ø¿ğY{ÿ Ækñ’ÖÖ[Û¨má]óLë.@Ë€9÷¯lÿ †"øåÿ DïRÿ ¿°ÿ ñÊì¡V¥4ÔÌ*S„šrgéü7ßÀOúâÿ Áeïÿ £şïà'ıñà²÷ÿ Œ×æ‡ü1Ç/ú'z—ıı‡ÿ Qÿ EñËş‰Ş¥ÿ aÿ ã•Óõšÿ Éø3cKùÒmSöôøs¦İÃ¢g’E_ìËŞIR ÿ S_î?ğÄ_¿èê_÷öş9Gü1Ç/ú'z—ıı‡ÿ W=iU­nhíäÍiªtïfxuî?ğÄ?¿èê_÷öş9Gü1Ç/ú'z—ııƒÿ W?²©ü¯î5çsÛà”ZáµøÅâí ¾óBûFÜğZ)ãQøâVıkõ¿5ÿ `ßÙÏâ—ÂOÖÚÏ‰ü}£èòé×6²İJñ²©`¬ íby(;WéE{XDÕ;4yõìçtpÿ ¾‹|]á PMªéÒÃnÒ}ÕœÑ1ö*Â¿µ-6ëGÔ®¬/`{kÛY^	áa£‘IVR=A~ı×æWü£öc—Ã~"oŠ³fÒ5IkqÆ	÷G…Ÿ’N=öÅcŒ¤å5Ğ¼<ù_+ê|%EWŒzEP |øíã?€Ş%]gÂ:«Z;`\YM—µº_îËpŞÄa‡b+õköcıµ<ûDZÅ¦»/‡|d‰™tk©& rÖîqæû~ğî02j[;Éôû¸n­g’ÚêY"š(ñ¸9¬9EtÑÄJ‹ò1©J5=OèvŠø'ö:ÿ ‚…'Š&°ğOÅ¨íµVÛˆä!#¹n‹ÇesĞ?FïƒË}í^í:‘ª¹¢y²ƒƒ³
(¢µ +˜ø£ÿ $ÏÅßö¼ÿ Ñ/]=sä™ø»şÁŸú%êe³Üü¢Š+åh(¢Š ıÈı’¿äÙşÿ ØÛÿ A¯Z¯%ı’¿äÙşÿ ØÛÿ @¯Z¯§§ğ/CÆ—ÄÂ¹oŠßòK¼cÿ `kÏıõÔ×-ñ[şIoŒìyÿ ¢ª[1-ÏÀj(¢¾Xö‚Š( Ü¿Ù?şM§á§ı€­ôX¯W¯(ı“ÿ äÚ~Ø
×ÿ EŠõzúŠ<i|L(¢Š²Jº¯ü‚ï?ë‹ÿ è&¿zş†5_ùŞ×ÿ ĞM<õäã¾ÏÌîÃu
(¢¼£´+ö«öÿ “Oø{ÿ ^³ÿ éLÕø«_µ°¿üšÃßúõ›ÿ Je¯GüGèrb>{½Q^Ñç…|{ÿ Jÿ “pÓìbµÿ Ñö|ÿ Hÿ “oÓÿ ìbµÿ Ñ7Ïˆşi|hü¢Š+çX+»øÿ %×áÏıŒšoş•G\%w ÿ äºü9ÿ ±“Mÿ Ò¨ê£ñ!Kf~óÑEõ'ŠQE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE yßíÿ &óñCşÅmSÿ I%¯ÂıŞı£¿äŞ~(Ø­ªé$µøC^>;â‰ß†ÙŸfÿ Á)ÿ äá¼Cÿ bµÇş•ÚWê¾Úü©ÿ ‚SÉÃx‡şÅkı+´¯ÕzëÁÿ Çñ‰¶´´WuÎa6Ñ¶–Š.m£m-\ÛFÚZ(¸	¶´´Qpmih¢à&Ú6ÒÑEÀM´m¥¢‹€›hÛKE6Ñ¶–Š.m£m-\ÛFÚZ(¸	¶´´Qpmih¢à&Ú6ÒÑEÀM´m¥¢‹€›hÛKE6Ñ¶–Š.m£m-\ÛFÚZ(¸	¶´´Qpmih¢à&Ú6ÒÑEÀM´m¥¢‹€›hÛKE6Ñ¶–Š.m£m-\ÛFÚZ(¸	¶´´Qpmih¢à&Ú6ÒÑEÀM´m¥¢‹€›hÛKE6Ñ¶–Š.m£m-\ÛFÚZ(¸	¶´´Qpmih¢à&Ú6ÒÑEÀM´m¥¢‹€›hÛKE6Ñ¶–Š.m£m-\ÛFÚZ(¸	¶´´Qpmih¢à&Ú6ÒÑEÀM´m¥¢‹€›hÛKE6Ñ¶–Š.m£m-\ÛFÚZ(¸	¶´´Qpmih¢à&Ú6ÒÑEÀM´m¥¢‹€›hÛKE6×â÷íıÿ 'uãÿ ÷ì¿ô†Ş¿h«ñ{öşÿ “ºñÿ ûö_úCo^v;økÔêÃülûwş	f¿ñ™õñ×şˆ·¯°ö×Çÿ ğK_ù6íCşÆ+¯ıo_`×Nşê&U~6&Ú6ÒÑ]2mih¢à&Ú‚ûO¶Õ,®,ï-â»³¸¡šŞt¨À†VSÁ<jÅ\ÆßÛWöOºıœüh5&9'ğ&±3:rKI9ck!<ä•'ï(=J¶>l¯ßï‰_ô‹^
Ôü+â[1{¤êì‘z<mÕdFşS‚¨ïÒ¿hO€úïìóñûÃÂ´ÖÙ2éÚˆM±ŞÛ’vÈ=f\ü¬ä`ŸCÙ¾hìÏJNug¹ætQEp!EPÑ?ğOù;ï ÿ ÜCÿ M÷5û7¶¿?àÿ òwŞÿ ¸‡ş›îkör½¼ğß¯ùv#ãBm£m-èÜåm~2ÿ ÁBchÿ kÏ–ä0°+ôû¸ş•û7_ÿ ğS-é?µåÑB£SÒ¬îÁÇŞÂ´9ÿ È8ü+ÏÆëMzXŒùJŠ(¯ôB¾¬ÿ ‚gŒşÔV_ö
¼ÿ ĞV¾S¯¬?à™%Gí?o»©Ò/6ıpŸÓ5½âÇÔÎ§ÀÏ×M´m¥¢¾–ç&Ú6ÒÑEÀş~ş"¿™ñÄíıíRèÿ äV®zµ¼_?Ú¼Y­Í×Ì¾ÿ 9ÖM|£Üö–ÁERúÿ ‰šóâœ™;4Å#Ü›¬#_£kàßø$¾†müãıcfïS·´¾TLøü<ÿ Ö¾ó¯¡ÂéF'—[øŒM´m¥¢ºî`&Ú6ÒÑEÀM´m¥¢‹€›hÛKE6Ñ¶–Š.Î?ğPÁÿ ‹ã÷¬?ôºŞ¿+öcş
ÿ &‰ã÷¬?ôºŞ¿ëÃÇzàaEWu}áÿ ™\üBñéÿ ¨\ú8×Áõ÷ü_şJìş5Õ†ş4Lk|ı3ÛFÚZ+è®yGÇÿ ğPÙgş÷‚O¼9hÆÓEåïì×,Ñàutå—¹—’W“Uı×äÿ ü;öZ?	üh|wáË2¾×çcq)„Óï%0:$œ²ö2ğ6çÉÆQÿ —‘ùÔ*}†|uEW’v…Ioq-ÄSÁ+Á<L9cb¬ŒC9õû-ûşÓ‘~ĞßVßT™WÆš%¾§À7+Œ%ÒGÁİ ÏÑÛkğ_àÆ-oà?Ä+ÅúoÑŠOjìV;¸‰!f+@¯Ü?†¿4O‹Ò|Wáë¡w¥jP‰coâCÑ£qÙÕR;kßÂ×ö‘å{£Ì­O‘İlt»hÛKEvÜçmih¢àyÏí ­ÿ óñ?gŞÿ „_Sÿ ÒY+ğ‚¿}~2iÇXøCã› 2n´+è ÿ zİ×ú×àUxøïŠ'~fQEyg`è¥xeI#b®„2°ìGC_ĞÆ“}­¥ÙßEşªê™>Œ ç_Ï-~âşÈ<ÿ fÿ kşeÂiÉcp{ù¶äÀÄû“ïøz˜ZRG%h™ëûhÛKE{8ÛFÚZ(¸	¶´´Qpmih¢à&Ú©«hö:ö—w¦êV°ßi÷q43Û\ xåFee<E\¢“¿¶Oì!¨ü’óÅş
ŠmWÀÌÆK‹a—ŸJÏ÷»¼>Õz7÷Ç•ıÉÍ#ªº0*ÊÃ ƒÔ_şÕÿ ğNO5÷ŠşÃªœËqá¬„·¸=I·'ˆØÿ pü‡±^‡È¯„ûTşãº~“?4(«Úæ…¨øgX»Òµ{7S´ÅqiuXœuVSÈ5F¼³´(¢Š +ô¿ş	Ûû^]øËÊø]ã;ï?U·€P‰’æ4kw'«¢‚Ê{ªyQŸÍ
Óğ¿‰uøL×t‹—³Õ4ë„º¶)"0`~™;Öôjº2æFu ¦¬ÏèGmk–øSãûOŠŸü7âÛ-¾¯c×–|·+óÆ}Õ÷)÷Zê«é®®'mÛ\ÇÅÿ ‹gâïû^è—®¢¹Š?òLü]ÿ `‹Ïıô¤ô`·? h¢ŠùSÚ
(¢€?rd‘ÿ Ïğ×şÀvßúzŞÚò_Ù'şMŸá§ı€í¿ôõºúŠoÜ^‡/‰‰¶¹_Šëÿ ·Æ?ö¼ÿ Ñ]]r¿?ä–øËşÀ×Ÿú!ê¤ôb[Ÿ€”QE|©íQ@¹Ÿ²pÿ Œiøiÿ `+_ı+Ö6×”~É¿òm??ìkÿ  
õŠúŠoÜG/‰‰¶´´V—$©ª¯üJï?ë‹ÿ è&¿Zş†õoùŞ×ÿ ĞM<•äãşÏÌîÃu
(¢¼“´+ö·öñ‰ÿ ?ëÒoı)–¿«ö·öÿ “Oøyÿ ^“éLµèàˆıLGÂwÛFÚZ+Û¹ç‰¶¾?ÿ ‚¤/ücmıŒ6¿ú&zû¾@ÿ ‚¤É¶ÙØÁkÿ ¢§®|CıÔi|hü™¢Š+æÏX+»øÿ %×áÏıŒšoş•G\%wŸ ¿ä»|8ÿ ±“Mÿ Ò¨ê£ñ!Kf~ôm£m-õW<Q6Ñ¶–Š.m£m-\ÛFÚZ(¸	¶´´Qpmih¢à&Ú6ÒÑEÀM´m¥¢‹€›hÛKE6Ñ¶–Š.m£m-\ÛFÚZ(¸	¶´´Qpmih¢à&Ú6ÒÑEÀM´m¥¢‹€›hÛKE6Ñ¶–Š.m£m-\ÛFÚZ(¸	¶´´Qpmih¢à&Ú6ÒÑEÀM´m¥¢‹€›hÛKE6Ñ¶–Š.m£m-\ÛFÚZ(¸	¶´´Qpmih¢à&Ú6ÒÑEÀM´m¥¢‹€›hÛKE6Ñ¶–Š.m£m-\ÛFÚZ(¸	¶´´Qpmih¢à&Ú6ÒÑEÀM´m¥¢‹€›hÛKE6Ñ¶–Š.m£m-\ÛFÚZ(¸	¶´´Qpmih¢à&Ú6ÒÑEÀM´m¥¢‹€›hÛKE6Ñ¶–Š.m£m-\ÛFÚZ(¸	¶´´Qpmih¢à&Ú6ÒÑEÀM´m¥¢‹€›hÛKE6Ñ¶–Š.œşÑËÿ óñGşÅmSÿ I%¯Â
ıàı£ÿ äŞ~(ÿ Ø­ªé$µø?^6;â‰ß†ÙŸgÁ)yı¡üEÿ bµÇş•Ù×êÆÚü¨ÿ ‚Rÿ ÉÃø‹şÅkı+³¯ÕšêÁÿ ÇñÛFÚuÜsÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@Û_‹¿ğP?k¯ÿ ¿eÿ ¤6õûI_‹ğPöºøƒŸËÿ HmëÎÇzXŸpÿ Á-FfÛÿ û®¿ôM½}ƒ¶¾>ÿ ‚Y°?³n¡ƒœxŠë?÷æŞ¾Ã®šÂ‰•_ÛFÚuĞd7miÔPvÑ¶E 7myí9û:èÿ ´‡Ã{ôÇg«Û“q¥j…2Ö³ã¡îcºËÜ`õU#×¨©”T—+n.èş}<sàgá¿‹µ_ø†ÉôıgLœÁqó†‚B¤Á‡G°«õóöêı‘áøóáñ/‡-|}¤Cû€íI·cıñÉCë•<6WòâŞ[K‰ 7†hØ£Ç"•d`pA¡µ|õj.Œ­ĞõiÔUÆQEÎj}ÿ ÷ÿ “¿ğıÄ?ôßs_³»kñş	ëÿ 'àûˆé¾æ¿g«ÛÀÿ úÿ ‘çb>47miÔW rÛ_ŸğV·üP^1‰2Ÿ¿Ò.Â3ÿ ÿ *ı¯ı±>·ÆÙ÷Äú´&mZÚ/í-9TeÄ °E÷ußü¹ëÃ›F´åË4ÏÄ
(¢¾põ‚¾¢ÿ ‚lÏäşÕzgu…êäoı–¾]¯~ıƒu¯ì/ÚÇÀ3ÂO=Å£_6ÚTşúeü«j:T©>~Ôm£m:ŠúSÇ¶™<‹o’¹ÂF¥ĞÔµÏ|DÔ†ğÿ Ä×ìv‹].êr}6ÄÍı)=Ïçîâf¸¸’Wûò1sõ'4Ê(¯•= ¢ŠŞğƒuˆ4Ñ<1¤ÆeÔukÈìáÈìãìI=€4÷Ğ×/ø'o‚_Áÿ ²ß‡å™|¹õ«‹U×öFãCø×ÒÛk3Â¾²ğw†4NÊÓô»H¬­ÓÒ8Ğ"È
Õ¯¦„y"£Øñ¥.fØİ´m§QZ7miÔPvÑ¶E 7miÔPvÑ¶E |áÿ _øÄ?½aÿ ¥ÖõøÍ_³_ğPïù4?ÿ ½aÿ ¥ÖõøË^&7ø‹Ğôpÿ 
(¢¼ó¨+ï/ø$¿?¼{ÿ `¸?ôi¯ƒkï/ø$·ü”/ÿ Ø.ık§üX˜Öøúi¶´ê+è(nÚÀñÿ €ôo‰¾Õü/âQy¤j'‹¡ÇPÊ{2°b ×CE-ôc?ş>üÖ¾ üLÔü%¬©B|ë+À¸[»f'Ë•~ GfV«Îëö“öĞı™mÿ hÏ†r%„QGã5Æ“rÜ8ùíØÿ v@3Ñ‚™Ïãåœú}äö·PÉmu´RÃ2•xİNXAA¯ŸÄQö2òg©J§´™Q\¦Á_[ÿ Á?jcğcÇÂ"¼)à½~ePò–Âñ°«7²6¿¦º)ÏÉV”æéÉJ$Ê*JÌş‰¶Ñ¶¾1ÿ ‚uşÔÿ ğ³ü"¿|Iu»Å:¸ûÄÏ–¿³¼™#áOªí<Æ¾Ğ¯£§QTŠ’<™EÁÙÛFÚu¡SÛÇsÊ¡â‘J:·BÁüûøãÃ3x/Æšÿ ‡îA[*ş{uİŒ‡ÿ A¯è6¿¿à£_À?´¦­¨EM?ÄGªÂØùw‘²aŸ_1ˆÿ lW›à¥ØëÃËŞhù‚Š(¯ô¿H¿à”¿¢¹Ñ|OğŞòP.m¥şÙ°yxØ,s(ÿ u„mÿ mÒ¿7k¼øñkPøñS@ñ–†gÓ§ıı¶ìˆ–3õF8'¡Áí[Ñ©ìæ¤gR<ñhıêÛFÚÄğ/´‰Ò<O İÍ#T·[›yqƒ´õ;09R;Gjİ¯£½õG’7miÔSİ´m§WñÇã&‰ğá¾«âıq·Cj»-íU°÷w—
{±À1è&ÔUØÒ»²;½´m¯†ÿ àŸ´§Å_^0ñeŸ‰şÍªøjØ5é¿hü·²šY	Ú2>òæ­’¡>ö0ÜÕ:Š¤y‘R‹ƒ³9ïxÊÇáÏ‚õ¯j‘ÜK§i6²^\-ª—ËA–Ú	 œ{Šã¾şÑ^ı£´=OSğ”—ŠšuÀ¶¸·¿„E2’¡•¶†a´òÏToJ·ûIB.?gŠFïø¦52¸µÔWæü¯ãT?
~<C¥jWßDñLK¦ÌÎp‰q»6îàE“ÛÍ5JÎXÅìÍ!OšõGìÚ6Ó¨®³Ãÿ hïÙÁ_´†–[Uƒû'Ä±&Û]~Î1ç AÀ•?Ù<á+“_“¿fŸşÎºø²ñ5™§Nì¶ZÅ®ZÖèFşÇTlî9¯İ*Çñoƒô_x~ïCñ—m¬i7i²kK¸Ã£_b:‚9A®:ØhÕÕhÎŠuœ4è=ôWÔ?¶Wì]©~Îú£kúªøò]±\7Í.ìx†cÜ‹'~‡¾^¯p•9rÈôc%%tQEAGëÇü/ÄRëŸ²ı­¤¹tZîÆ?e%'Çç9¯«¶×Ç¿ğK-.KÙ¿R¸ue[ß\Î„÷Q¼yŠ7å_bWÒPş}&§ÆÆí®âhğ‰¢ë¿L¹_Î&ÑU=bËûKI½´ÿ Ÿˆ.ÚR?­lö3GóÇE+«FÌ¬6²œ{JùSÚ
(¢€?rÿ d…ÿ Œføiÿ `;ı½omx÷ìw8¸ı—ş09H‰?ï’WúW±×ÓÓø¡ãËâcv×+ñaâÕøËşÀ·Ÿú!ë¬®CãÄV	<m4ò¤0¦‰zZI*¨ò’OJ©lÉ[Ÿ€ÔQE|±íQ@¹ß²jÿ Æ4|4ÿ °¯ş€+Ö6×”~É¿òlÿ ÿ ìkÿ  
õšúz<i|LnÚ6Ó¨­	)êËÿ »Ïúâÿ ú	¯ç¿¡İ[şAwŸõÅÿ ô_Ïy8ï³ó;°İBŠ(¯(í
ı®ı…Wş1?áçızMÿ ¥2×â~×şÂ¿òiß?ëÒoı)–½ñ¡ÉˆøQîûhÛN¢½£Ï¶¾Aÿ ‚£¯üc]§ıŒ¿ú*zû¾Iÿ ‚Cç~ÌlØÿ U­Z?é"ÿ ZÂ¿ğ¤kOãGämQ_6zÁ]çÀ/ù.ß?ìdÓôª:àëºøâ?_Üô_éÄÿ àTuQø¥³?{6Ñ¶E}IâÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚuç?´zÿ Æ<|Qÿ ±[Tÿ ÒIkğv¿xÿ iù7Š?ö+jŸúI-~Wø¢wá¶gÙÿ ğJ_ù8Ø­qÿ ¥vuú³_”ßğJ_ù8Ø­qÿ ¥vuú³]X?áâ>0¢Š+¸æ
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
ü@ı´5¨uÿ Ú“â=Ô½#Õ©?íB‹Á£#ğ¯Øï‹ßt¯ƒß5ßkG®›lÒ$nØ3ÍŒGÿ ´íµG×Ò¿µíjïÄšæ£«ßÉæßj2]\IŒn‘Ø³Ä“^^:JÊ'n:¹§ğIİrƒ¾0ÑÃæâÓ^ûS'¢Ko©üLùWÜUùÿ ßøÙiğ³ãdÚ¯s‹â¨Í®&p‰ÒĞ'ŒÒGõ‘}ëõÚº0²R¤—c*Ñ´ØQEØs…Q@Q@Q@~yÁE¿cö»ÿ üh<Ä_3Ä:tà:Ş Àÿ YôßıòCi’Ä“ÄñJ‹$n
²0È`x â±«MU+4„ÜÑüíQ_[~Ş_²#|ñGü%³øA5iÕÆ24ËƒÉ„úFÜ”=¹^À·É5óÓƒ§'¬d¤®¢¿à¿òwşÿ ¸‡ş›îköz¿à¿òwşÿ ¸‡ş›îköz½|ğß¯ùQ^ÊQE ~9~ß_³äŸ~2İjš}¾ßx™ä¿²dR‰ÌğÀY·ı×QØ×Ì•ûÁû@|Ñ?h†š…u„XäyÖÛröw*ÉWÛ’î¥‡züFø•ğß_øKãMKÂŞ&²kZÂM½REê²#0äCøWƒŠ£ìåÌ¶g§F§2³Üæ+ºøâDğ‡ÆÏkR¿—»e4Íé7ÿ ã¹®€JG¸Ó³¹»ÕXşŠ(¯ı“ş6Züwø% kÂåeÖmá[^2Ùt»@v#°~$Î;ƒ^Ã_Q)%$xí8»0¯#ı®<D</û3üI¾-°¶‹=¢·£N<‘øæA^¹_ÿ ÁS¾3[i>Ñ¾Y\+êZ´ë¨j«dÅk>XaşÜ˜#ş¸ŸQYV—%6Ê§i¤~cÑEóg®úÿ ºıŸ^ëP¾ø³¬[‘¸“OĞÕ×‡‚³Î?İŒ„´Ÿİ¯˜?eÙŸYı¥> G¦[ù–^±+6¯ªâ³÷<(íË¯Ú¿xgLğo‡tíF³ÃJÓàKkkh†F ÷>§©9&½%){G²9+Ô²åF¥Q^Ñç…Q@Q@Q@Q@Q@4ÿ ÁE®D²O‹œ§±Œà\-ÿ ²×ãm~«ÁSümi£üÑü7ç§ö³«Ç"ÛîŒ£³¾=´Cş_•5ácêüK½À¢Š+„é
û£ş	3r«ñSÆÖäüÒh© Ë:ÿ ¡
ø^¾®ÿ ‚gøÒÛÂ¿´Í½•ÔÂ×´»62Ç
dÊL£>§É z’zèÃ»U‰•]`Ï×ª(¢¾ŒòBŠ( ¿8?à¥²É³¹“âç†,É†fXüAkpÂ¥Şfá_ßkwc_£õSVÒlõí.óMÔm£¼°¼…íî-æ]É,n¥YXw?Æ­5V<¬Òp•ÑüñQ^åû^şÍ×Ÿ³Åôè’I|-©¹Ñoî-Fè˜ÿ ~2BŸPU¿‹áµó’‹ƒqg¬š’º
(¢¤fçük«ü9ñv“âmèÙjú]ÂÜ[Ì@aØêFA$w¯ÛïÙÛã¶ûBü2Ó¼Q¥²Et@‡Q°–³¹ oŒûs•=Ôƒ× ~×¶~Éÿ ´¦§û6üH‹SM÷^Ô
[ëzŸõ°ç‰ÏDÉ+ë’½5Ù†­ì¥g³0­Oi¹ûyEfxgÄÚ_Œ¼?§ëš%ô:–“
Ïmue$B8#ÓÜAEi×¼y_(ÁF¾¿ÅO‚-â6ÜÍ®øMú5EËÉjÀ„v
²Û3ë_WÓd&£‘Uã`U•†A¨"¢¤H¸¾¥FN-4;4WÒÿ ·'ì»7ìûñµ-&|¯JóiìŠvÚIÏjŞ›z¯ªú•jù¢¾nqtäã#ÖŒ”•ĞQEYşÃ¿¶iıŸõ|/â¦šçÀº„Şh’0]ôÙ‘W©ŒãæQÏ‡9úËáÿ i+Ñlõ}şßSÒï#ÛŞZÈ$U=Õ‡ZşyëÓşşÒŸ¾İ³øG_–ÚÊFß6™r¢kI©¸ı¥ÚŞõßCé®YjZ”yõ[Ÿ»4Wæ‡‡à­#µµU×~éz•È2iú„–Šà.’ÿ :£âïø+Œu+7‹Ã
Ñô9˜çß\ÉzWİ@Œır=«Ğúİ+ns{	ö?Eş#|HğïÂŞø—Å:œZV‘h>id9gcÑG.ç(ä×ã¯í1ûHx—ö²ø‘g½¬ğhÑÏöMAŒîl»ÀpÓ9Ûœp8QÓ'Îş(üfñ§Æhj2ñÖµp¹G! ö%} Ïzû£ş	Ãû$Ë§µ¯Å¿Zl’HÏü#ö§ÌªFÛĞ‘Ä·÷MqÊ¤±RöpÑ¨*+š[ŸTşÊ?àı¾i^tµÉÇÛu{ˆğ|Ë§r†î¨ Aê=Í{W­¨¥q6äîÎCã—ö—Â?Ùã?hĞï¢ÇûÖî?­~Ç#E"º1GS•e8 â¿¡Xÿ iè:•3ö‹i"ÇûÈGõ¯ç¼¬vñgnf~ÏşÄ¿´µ¿í	ğ¶õ…1ĞÑ-uXY²ó`a.G´€ú0aÓú&¿~ü`ñ'ÀßYx§Â÷g¾ƒä–ŠæAh¤^êp=Á Œ~ÇşÎµ„i‹½álµëxÔêïûûfèY¿z8ãÑ†Ä*‹–[˜Ö¥Êî¶=Š(®ó˜Æñ—„tÏxOVğæ³n·Z^©lö·0*ÃèÃ¨=ˆµ~ø»Ã³øGÅZÖ…rw\éw³YJ}Z9êµıêº¥¦‡¥Şj7÷	icgÜ\\Jp±ÆŠY˜Ÿ@?…~üHñ@ñ·ÄOøˆ_UºÔ Û·ılÍ'Nß{¥yXë{¯©Û†¾§;Eô_ì/û?ÍñËãM„×–Şg…ü<ñê:œŒ¹IlÅ÷‘—‘ıÕzó!9(®§d¤¢®ÏÔ/Ù;áÛü+ıüáùã1Ş%€ººV2Í;Oº™
ÿ ÀkÖè¢¾š1åJ+¡ã·wp¢Š*„~|XğûøOâ—Œ4YciúÅİ®1¹3¨ü8®R¾²ÿ ‚“ü#ŸÀ?§ñ$2èş+…oc/È·(N™õÈYıu¯“kæ*EÂn,ö!.h¦QEfYû1ÿ óñD>%ı”ü'HãK’ëO¸P~ã,îêı³’3ø×Òø¿ûş×Z‡ìÇâèn¬¤Ö|%ª•kë˜,±È …š"xİƒ‚§†9¾ü_ø)?ÀÖÓ–äëzšÌFM™Ò¦óGãŸøõ{´1pJNÍmJRæm#êJüÃÿ ‚ˆ~×Éãkéşø3Px~ÒLk7öí•½™ND
Ã¬hFIèÌeÉÎı§?à¤š·ÄÍ÷Âÿ ¬n|3¡Ü†ŠçU¸p/®c<P¤ˆT¸%ˆî¼ƒòÇÃ?„>#ø¬úóèv›ì´-6ãVÔ¯eÈ†Ú£gÁ8åßaT^¤ú(f^lF#Ú~î™­*\¾ôÎ.Š(¯,í
(¢€?r¿dYÅÇìÉğÕ‡A¢[§ıò»¥zõ|ùûëÑëÿ ²™t–‘ÜYÊ½Õ£¸ à;OĞŠú¾š›¼"üzI…QZSÖG£ß1è°HñÓ_Ï~ü|bñ~øMã=nG®Ÿ£İÜ†?ŞXXõ'ñ¯Àzò1ÏX£»³
(¢¼³´+ö³ö˜Mû%ü=eè-®ñSü«ñN¿`ÿ àšş$‹\ı•ôk$”<š=ıå”‹İKJg ÿ Àgz'ûÇèrâ>êZ(¢½³Î
ù—ş
9¥¶£û'øe¾Åuepqéö˜ãÿ Ú•ôÕr>ÚüTøiâ_Şë2Ú‰ü²r¿$ŸUp­ÿ ¬êGš=Ë‹å’gà­âÏ
ê~ñ6©áıfÕ¬õ]6áín`q÷]N=Gp{‚d×ÌÀWKğÇT‹Cø“á=JvÙ­ipìN0©21?®jŠiÙÜè¢ŠøöTÿ ‚†ø;WğŸ |IÕáñ—n°Jå­ïÑ U}Ê	YH0 yœ§ã_ü‹á§€4iãğ…ßü'!e"mQÒÒ&Ç,¬@şêd™^µôKMÇšç“ì§{Xè¿mÚæ?ÙŸBÑítˆ-u?js¤±Øİn1¥¢7ïd}¤»Ü±çaìßş 7Å_†~ñyÑî´¬Z-ÚØ^2´‘«giÜ¼2°Ã)à•e$)È’ÿ ~ø«öâı£$¿ñUä×ÖaÖû]½ÎÅŠÙN»&ï¸ª:ÍÎÓ_²V¶°Ù[EooA(#(Ô*¢€  
’ªÜşÏBêEA(õ%¢Š+°ç
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€<ëöÿ “xø£ÿ b¶©ÿ ¤’×àå~ñşÒòoìVÕ?ô’Zü¯ñDïÃlÏ³ÿ à”|şÑ"ÿ ±Zãÿ Jìëõkm~SÁ(ÿ äâ<Eÿ b­Çş•Ù×êÕu`ÿ „cˆøÄÛFÚZ+¸æmih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š ò/~Ö_şø‚÷Cñ¬´İ^É‚\Y´3Hñ±PÀˆ{ñÏÿ ÁO¾xjÚqáôÕ¼_z¼D¶ö¦Ö?í<Øeáé[¿¿à¿>0|BÖ¼c¬ë~)¶ÔõYVYâ±»¶HT„TC[³…X×%ÿ £øIÿ C?ğ:Óÿ ‘kŠOv¢‘ÑeÖçÂ?´—ígã?ÚST‹ûaÓKğõ«‰-4+6&ßyÇ™Ç ' dçÄëõkşGğ“ş†/àu§ÿ "Ñÿ £øIÿ C?ğ:Óÿ ‘k‚XZÒw–çJ­N*Èü¥¯¶?gø)f¿ğßK³ğ÷Ä>è¶Ê±A©[¸„¶ì,ÀÆâ­êÆ¾‚ÿ ‡Qü$ÿ ¡‹ÆŸøiÿ È´Ã¨şĞÅãOü´ÿ äZªxzôİâ)U¥5fzO†o¿&µY?á5M*S÷­õ;I¡uú…OàÆ½KáÏÆo|]mDx;Ä¶^!:w—ö¿±±>W™»fr]ÿ |šù“şGğ“ş†/àu§ÿ "×´~Î_²‡„fâámGZ¿şÛû?Ú¶'†MO›³g—c>kg9è:wô)º÷\éXå’§ouÏ¶´´WQˆ›hÛKE &Ú6ÒÑ@	¶´´P?‹¼#¤øïÃ:Ÿ‡µÛ(õ#Q­îm¥¤~„G €G"¿ÿ jÙ¯Xıš~"I¤\™/tİÓé£/ƒÊ¶8&@aîá…~âW|zøáÿ Úáİï…µøöoıí•ò(2Ù\ BJŸş%$wÍrb(ªÑÓtoJ§³~Gå'ü×şOÀ?÷ÿ Ó}Í~Ïí¯Ë/Ùöcø•ğöÄğÔŞ ğ¥¥Ï¨C.µ»5“ƒep‰"ËŒmrËŒàå€ +õ:³ÁÅÆ›Mu+Ó’°›hÛKEwœÂm£m- ›kÃÿ jOÙOÃŸ´¿…D7{4ÏY!n´‰—¿• ş8‰íÔ‘ß>ãED¢¦¹e±I¸»£ğWã/Àüñ+hŞ/Ò^ÍÛ›{È³%­Òÿ z)1†ú0îyıB)ğ‹ããGñ•g¬éwÚ_B²ÆŞø#‚;Èí_üVÿ ‚WøÄÓKyàv÷ÁÓ±'ìW
om>‹¹„‹õ.ßJòj`äµ†§t1	üGÁ_³¯í%âŸÙ·Åï«è—V7AcÔt›‚|›¸ÁÈÉê®2v¸ädõƒúGà¿ø)ÁiqO¬ê—…/v2Ò÷O–p¾×\0ÏB@>Â¾?ñgüãƒ$ŸÙ_Ø~%„}Ãg}ä»u™Pÿ ?ZóûØKã½«ÿ /Oıs»¶ıSQNXŠ:(èT•*š¶}±ñkş
àisEà;+Ïk¸†k˜ÒÊ3êûñ#cû¡Fq÷‡ZüÎøñ]ø£ãKÄş%¿“QÖ5	<É¦~ ã
ª:*¨   ¯dÓ¿àŸÿ µ)G€ä·V<Éq¨ÚFÜƒ.!^¥à¿ø%OÄbdoøƒBğå±Æï%¤¼œÀ U?÷İ)ûzû ²§³>(¯¢?foØ£Æ¿´=ôïğj¸3kW‘ç/qn‡iÿ kîç<¾şÿ Á9ş|.šŞÿ V·›ÇÌLO¬ -‘‡u·_”:9zú–RŞ$Š$XãE
¨ƒ@à ;
Ş–ïSî3#¤OáoÂŸ|ğe—†<+§®Ÿ¦ZŒòO!t²·ñ;c“ô  :İ´´Wª’JÈâ½õbm£m-Ä&Ú6ÒÑ@	¶´´Pm£m- ›j–µ«ÙxwG¿Õu+•´Ó¬mäº¹¸“;b‰³±Ç` ŸÂ¯V?Œ¼/kãk¾’h¬u‹ôùä·`²,rÆÑ±BA€cŒ‚3ØÒw¶ƒ<Šóöàøb…¤ø‹§0ÿ ¦0Ï)ü–3^[ñ+ş
ğ·Ã:}ÂøRKÆ:–Ò!ÙlÖ–Û±Áw”>ˆj·ü:á'ı^4ÿ ÀëOşE£şGğ“ş†/àu§ÿ "×x—²HèJŠêÏÎOß¼OñûÇøŸÅ7)%Ë †ŞÖİJÁi$ˆãRI$œ’I$’k€¯Õ¯øuÂOú¼iÿ ÖŸü‹Gü:á'ı^4ÿ ÀëOşE®…­'vtªôÒ²?)h¯Õ¯øuÂOú¼iÿ ÖŸü‹Gü:á'ı^4ÿ ÀëOşE©úQıbå-YÓ5+½RµÔ,.e³¾µ•g‚â*ñH¤2²‘Ğ‚Ú¿TáÔ	?èbñ§şZò-ğê?„Ÿô1xÓÿ ­?ùŸÔê‡·Ãüÿ ‚¦iØ¶ºgÅ&ñ58”Fu­"%’)ğ>ü‘n×fàO@½Ğú_íëğU…]<o#%.¬n¢eö;¢ò5æğê?„Ÿô1xÓÿ ­?ùøuÂOú¼iÿ ÖŸü‹]±úÌU™Í/bÏ¨¾|VğÅİ&çRğ~»m¯XÛÍöy¦¶I´6Ó¸œëv×–~Ï?³†¿f¯êZ†/µ[ëKûÏ¶Êú´±I }Š˜Sh1…AïÍz¥vÇ›•snsÊ×ĞM´m¥¢¬“Ë?i€ºOíğ¿Qğ½ÿ —ş<ı3PeÉ´ºPv¿®ÓÊ°î¬{àÄxKUğŠ5Oë–c«i·mso å]N8õ¨#‚#ƒ_Ğ•|ÿ ñÿ ö#øyûEøª×ÄZüÚ¾“«Ãn-¤ŸEš¾Ò€åL¢HŸq^@# pbpş×ŞçM¼š=ÅZ+õkşGğ“ş†/àu§ÿ "Ñÿ £øIÿ C?ğ:Óÿ ‘kƒêuN¯¬@ü¥¢¿V¿áÔ	?èbñ§şZò-ğê?„Ÿô1xÓÿ ­?ù©Õ¬@ø»öYı´¼Sû5Èú_ÙÇˆ¼!q/™.“<¥?yà~vİH*}æ¿Aüÿ ø%âëä»ñÇ†oÖz½” ©=~xÃ¡ğ*âáÔ	?èbñ§şZò-ğê?„Ÿô1xÓÿ ­?ùºéÃMYY£	Ê”İÏ³öÑ¶–ŠôÎ3šø‰ğç@ø­àıGÃ&°RÒ/“d‘?§ø]ªºCAùûQ~Å¾.ı5)oâ_ø.G?gÖ­ã$Â;%Âõmş×İnÇ9QûEL¸·ŠêŞH'&†E(ñÈ¡•”ŒAêí\Õ¨F²×sju?CùØ¢¿]¾3ÿ Á6~|J–{ÿ ùŞÖ%bÌÚj	lÙsnHô{WÉŞ1ÿ ‚[üXĞî$şÃ½ĞüKkÕ+£m)ú¤€(?G?Zò'…«—;£Zê|uE{ı÷ìñæÁ˜Iğúæ@;Ái.ï™MM§ÿ Á?ş=j,|$Çî5+8ÂûfÏä+eSù_ÜiÏçÏTøa’âdŠ(ÚYd`©Y˜œ  êkî‡ßğJj“C/Œ¼W¥hV„åàÓï.1ıÜ°DSîcß¥}Ÿğ7ö5øeğd¾ĞôvÔµåP¿Û:»‰îİ8	ÕÄšè§„©-ôFR¯í©ò/ìkÿ óºÕî,<oñRÁí4äa5†nP¬³‘Ê½ÊŸºŸôÌòßÅÃ~“Ç
Å¢(DQ…U €
}ìR¥Q´N	ÍÍİ‰¶´´VÆbm¯çröÕ¬¯'·o½úƒŠş‰+óáüÏÅú×íu{ñÃŸ`ğ—7Û¾Ûı¼	“ØägPÛƒò«‚EyØºr¨â¢ºQ»gÂ©áŸk×-uSºÒ5[Vßåœ¦9ûØô#¡ıõ“áÏ…&ğÊxvOèòx}béoaµU=@‹nßÒ¾lø¥ÿ ÑøQã¶–çAKßj.s5üÛb}áœddÏ,ã¬]ÍV"/FşÁTüA Ù[éÿ ü;‰5Úum1ÖŞå½Ş"<¶?î”Õíÿ SøGö3*èŞ-3ã‹sen~¿hÆ+çoÿ Á*ş$è’;øo^ĞüMl3´HÏe9ÿ €0d÷ò¼¦ûö
øñ§Ìc‡÷2ó€Ğ_ZÈ§ß+)ãëG´ÄÓÑ¯À9hËTÎ»ö¦ı¾üEñûI›Ã&ŞğŒŒÄ>w™sz° *dgbõîM|¥_Køş	Ññ×\”,ş¶Ñ¢?òÛPÔíÀü£woÒ¾øCÿ £Óôû‹{ï‰(şÕÛ‚úF†­$ú5ÃaÙ}•ûÖ.•zÒ¼‘§=:jÉŸü	ıŸü]ûBxÂÃ6,aVóT™Hµ²»Høëè£–ì:‘û7ğà_‡gß‡Ö~ğü{öşöòşE[ÙÈ¥Ë„ ;fºŸøÃßü;m øcHµÑ4‹|˜ímjäõbz³ìI'¹­Úõ(aÕ^¬â©UÔÓ ›hÛKEu˜	¶´´P—~ÑŸ ôoÚ'á¥ï…õF×@ı£NÔnkK€WÇu9*ËÜĞàÅ/Š_
üKğoÆW¾ñVœúv§lxÏ1Î„ü²Äİî ~ÿ WŸ|høà¿:7‹ô¥»TÉ¶¾‡İZ1ş(¤ÆGlƒ•8¸±u[U¹ÑJ¯³Ñì~Q_q|Rÿ ‚UøßCº–ëº‰ôşJZß·ØîÇ¢ó˜ÛıíËô¯ºı…~;YÜdøy|Î{Åum"ÿ ßK)­y£R.Î'z©³<Šú{Âÿ ğN^"•×‡ì|?yÚ¦¥ õ+‘¿ñÚú—àçüÇÂ¸‡Pøƒ®MâË„!¿³,U­lÁôvÏ™ úl÷ªj³ébeZê|Iû6şÊ~1ı¤¼Aäé07Ãöìî½uò!Õ?ç¤Š=²Ts_©·ü/ğ+ö\ñÇ…ü-gä[jqu&÷r›W,­İä  W³è:›á}ÓIÑì-ô½2Ò1½¤B8¢QÙTQx£Ãöş,ğÎ¯¡İ¼‘Újv“YLğRD(ÅIgqkÖ¥‡(¾¬âW7ä=WêÇü:‡áGı2ÿ ÀËOşF£şCğ£ş†Oàe§ÿ #W›õ:§_Ö ~SÑ_«ğê…ô2xËÿ -?ùøuÂú<eÿ –ŸüGÔê‡Ö |§û~ÚÇöpûo‡<Caqªø:şqq›V{ˆ
ÎŠÄV er@ õôƒáOíağ³ãN©o¥xSÅŞkÆÒ.›5´ĞÎ‚[‡P =	áğê…ô2xËÿ -?ù»ß‚°?€~|B³ñ¬ø’óRµŠXRJâİá"D(Ù	àñÍvÑzv‹µj’§+µ¹ô¦Úò‰ŸµOÂ¿„:µŞ•â¯ÛiÚµª+Ë§¬2Í8¡—åD=Añ¯Y¯š~5~À¾;|DÔ<e¯ë>&´Õ/R(ä‡Mº·H@5pnŠ3–<×eNt¿w¹„yoï!~ØŸğP(ş3xnóÀş°¹Óü3tËöíN÷ä¸»U`Â4E'dd€NI,00£ üQ_«_ğê?„Ÿô1xÓÿ ­?ùøuÂOú¼iÿ ÖŸü‹^]L=zšGljÓ‚²?)h¯Õ¯øuÂOú¼iÿ ÖŸü‹Gü:á'ı^4ÿ ÀëOşE¬¾§T¿¬@ü¥¯£?c¿Úöóöb×5kİ>MgÂZ³+ŞYÂÁf†E,±dí'N7 ¼ŒWÙßğê?„Ÿô1xÓÿ ­?ùøuÂOú¼iÿ ÖŸü‹W5hKš;“*Ôä¬ÏPğoíİğOÆÂÖ;o%…íÃ¬kg©ZM›É .v•''³^ÿ ¶¾;Ò?à–¿
t]ZËPƒÄ1i­'IÑd½´*YX0ÃŒZû½Zn¥¿x)ò}6Ñ¶–ŠØÌøïöèıŠ[ã©ñ§ƒ£?YÁ²{3„]V%*î<	Tp¤ğF‘€Gå6µ¢êÕ®´ÍVÆãMÔmd1OiwG,N:«+ƒõ¯èv¼·ã7ìËğëãÕ¿üU¾ŠãPTÙ«jL‘ÀH¿xÊû—•ç×Âªš:3ªn]%±øOE~ø·ş	#g%Á“Ãç·ƒµ¾­§¬­ÿ #uÿ Ğ+˜³ÿ ‚IxšIH»ø…¤Ãxhl%ãèYyßU­ü§_¶‡sàzôŸŸ³÷Œ?h/C¢øcOv€0û^©20µ³Ní#ãÇE3vúğçş	_ğïÃ7İx«]Õ<bÑòmB‹g?í*“ÂA_^øOÁú€ô+}ÃºM¦‹¥[ŒGie©ÀêOry=ë¢MŞ¦†SÄ/²qŸ³ÿ À?şÏ?í|3 §/úÛíJD-ìøæFô•s… NIô­´´W¯¨«#¶İØ›hÛKEP„ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€<çö_øÇŠ?ö+jŸúI-~×ï/í!ÿ &ïñKşÅ]Sÿ I%¯ÁªññßNü6ÌûCş	Gÿ 'â/ûn?ô®Î¿V«ò›ş	Cÿ 'â/ûn?ô®Î¿W+§ÿ tc_ãE>Ší¹Î2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢ày¿í!ÿ &ïñKşÅ]Sÿ I%¯Áªıçı¤äİş)Ø«ªé$µø1^F;â‰ß‡ÙŸhÿ Á(äâ<Eÿ b­Çş•Ù×êå~Qÿ Á(?äâ|Eÿ b­Çş•Ù×êöÚéÂÂ¿Æ6ŠvÚ6×iÏa´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£ma´S¶Ñ¶€°Ú)ÛhÛ@Xmí´m ,6ŠvÚ6ĞE;mh¢¶´†ÑNÛFÚÃh§m£mcÍÿ iù7Š_ö*êŸúI-~WïOí$¿ñß¿ìUÕô’Zü¯'ñDïÃìÏ´¿à“ÿ òq>"ÿ ±Vãÿ Jìëõz¿(à“üşÑ^"ÿ ±Vçÿ Jìëõmtá?„c_ãŠ]´m®Ó(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mhÎ?i/ù7oŠ_ö*ê¿úI-~
×ï_í$¿ñ¿ÿ ìUÕô’Zü¯'ñDîÃìÏ´ÿ à“¿òq^"ÿ ±Vçÿ Jìëõ~¿(?à“¿òq^"ÿ ±Vçÿ Jìëõƒmtá?„a[ãŠ]´m®Ó(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mhÎ?i/ù7_Šö*ê¿úI-~	×ïoí$¿ñ¿ÿ ìUÕô’Zü¯#ñ#»³>Ôÿ ‚NÉÅx‹şÅ[Ÿı,³¯Öüÿ ‚NÉÅx‹şÅ[Ÿı,³¯ÖéÂÿ Æ·ÆQEuœáEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPœ~Òòn¿ÿ ìUÕôZü¯ŞïÚOşM×âŸıŠº¯ş‘Ë_‚5åc>$waögÚŸğIÏù8¯Ø«sÿ ¥–uúÃ_“ßğIÏù8¯Ø«sÿ ¥–uúÃ]X_á˜ÖøÂŠ(®³œ(¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š óÚOşM×âŸıŠº¯ş‘Ë_‚5ûßûIÿ ÉºüSÿ ±WUÿ Ò9kğB¼¬gÄì>ÌûWş	7ÿ 'â?û.ô²Î¿Xëòsş	7ÿ 'â?û.ô²Î¿Xë«ü3Ïß
(¢ºÌ.QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .QEp¢Š(…Q@\(¢ŠáEP
(¢€¸QEÂŠ( .yÇí)ÿ &éñOşÅMWÿ Hå¯Àêıñı¥?äİ>)ÿ Ø©ªÿ éµø^V3âGnf}­ÿ ™ÿ “ŒñıŠ—?úYg_¬•ù7ÿ ™ÿ “ŒñıŠ—?úYg_¬•Õ…ş•oŒ(¢Šë1
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€<ãö”ÿ “sø©ÿ b¦«ÿ ¤r×àe~ùşÒŸòn?ìTÕôZü¯+ñ#²†Ìû_ş	/ÿ 'â?û.ô²Î¿Y«ògş	/ÿ 'â?û.ô²Î¿Y«£ü3ßQEÖbQE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE y¿í+ÿ &çñSşÅMWÿ Hå¯Àªıõı¥äÜş*Ø©ªÿ éµø^^/âGe™ö¿ü[ŸÚ7Äö*\ÿ ée~³í¯ÉŸø$¯üœoˆÿ ìT¹ÿ ÒË:ıh®œ/ğÌk|cvÑ¶Eu˜ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@ÛFÚu İ´m§Q@mûK/üc—ÅOû5_ı#–¿+÷ãö–ÿ “rø«ÿ b¦«ÿ ¤r×à=yx¿‰”6gÛğI_ù8ßÿ Ø©sÿ ¥–uúÑ_’ÿ ğI^hïÿ Ø©sÿ ¥–uúÓ¶ºp¿Ã1­ñ‰E.Ú6×Yˆ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´æß´·ü›—Å_û5_ı#–¿ë÷çö–_øÇŠ¿ö*j¿úG-~W—‹ø‘ÙCf}³ÿ “ÿ “ñıŠ—?úYg_­5ù-ÿ “ÿ “ñıŠ—?úYg_­5Ó…şoŒ(¢Šë1
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€<ßö–ÿ “qø«ÿ b¦«ÿ ¤r×à-~ı~Òßòn?ìTÕôZü¯/ñ#²†Ìûgş	'Ïíâ?û.ô²Î¿Zö×ä¯üKşN;Äö)Üÿ ée~µ×Føf5¾16Ñ¶–Šê±ˆ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹ ›hÛKE6Ñ¶–Š,m£m-XÛFÚZ(°	¶´´Q`mih¢À&Ú6ÒÑE€M´m¥¢‹æ¿´ºÿ Æ8üUÿ ±SUÿ Ò9kğ¿?iù7Š¿ö)ê¿úG-~×›‹ø‘ÙCf}µÿ ’ÿ “ñıŠw?úYg_­uù+ÿ ‘ÿ “ñ'ıŠw?úYg_­µÑ†ş•oˆeú+ªæ(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.™şÒÿ òn?ìSÕôZü¯ßÿ ÚcşMÃâ·ıŠz·ş‘Ë_€æbş$vPÙŸmÿ Á$äãüIÿ bÏş–Y×ëm~IÁ$?äãüIÿ bÏş–Y×ëmta¿†gWâ
(¢ºŒBŠ( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( 6ı¦?äÜ>+Ø§«éµø_¿ÿ ´Çü›‡Åoûõoı#–¿ +ÍÅ|Hë£³>Üÿ ‚Gÿ ÉÈx“şÅ;Ÿı,³¯Öêü‘ÿ ‚Gÿ ÉÈx“şÅ;Ÿı,³¯Öêßü3*¿QEÕ©‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…‚Š(£P°QEj
(¢BÁEQ¨X(¢Š5QF¡`¢Š(Ô,QE…6ı¦?äÛş+Ø§«éµüÿ WôûLÉ·üVÿ ±OVÿ Ò9kùş¯7º:èìÏ·?à‘ßòr$ÿ ±Nçÿ K,ëõ¿m~Hÿ Á#ää<Iÿ bÏş–Y×ë…taßîÌªüCvÑ¶EtÜÄnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀnÚ6Ó¨¢à7miÔQp¶´ê(¸ÛFÚu\í£m:Š.vÑ¶E»hÛN¢‹€İ´m§QEÀóOÚeãş+Ø§«éµüş×ôûMÉ·üWÿ ±OVÿ Ò9kùü¯;º:èìÏ·à‘¿òr$ÿ ±Nçÿ K,ëõÂ¿#ÿ à‘œşÒ$ÿ ±Nçÿ K,ëõÇmo‡øªüBQK¶µÒb%»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m y§í5ÿ &ßñ_şÅ=[ÿ Hå¯çò¿ _Úiãş+ÿ Ø§«éµüıWŸŠİtvgÛßğH¿ù9Ø§sÿ ¥–uúå_‘¿ğH¿ù9Ø§sÿ ¥–uúå[áş:«Ş
(¢ºL¬QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
*¦­¨ei·WŸf¸¼ò#i>ÏiùdÀÎÔ^äöåú·í%£h;F§á¯iÅ¸_µØ$YúnTJ¤c¤™Jœ¥±ëtW‰ÿ ÃZxCşÚßıø‡ÿ Ñÿ iáúk÷âş;Qí©÷+ÙO±í”WŠ'íeàöli«Aô–º=ö„ğ.¹"F5°ÊÇo¢h‡âøÚ?TªB[1:r[£Ò(¨íîb¼&‚Tä’6¬=Aj[PşÉÓ.¯>ÍqyäFÒ}Õ7Ë&v¢÷'°­²»!Fú"İãµO…ì.ŞçG×í§ŒáâšÖ%e>„r*/økOĞ;[ÿ ¿ÿ ñÚÇÛSîiì§ØöÊ*®—¨G«évwĞ«,7P¤è² +(`\µ[™QHv
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Šâş!|T°ømä6£¥ê×vò©o´ØÛ«Å8ÚÌÌ j™IE^CQrvGiExŸü5§„?è­ÿ ßˆøíu_ş7h_5yôí.ÓP‚xa3³]Ç®ĞÀ`mvç‘Ú¦5!'h²¥NQWhô*(¢´"ÁEP
(¢€°QE‚Š( ,QE`¢ŠdÒy0É&Æ“b–Úƒ,p:ëFÚ‡-ôExİ÷íIá­.å­ïtOÚ\/ŞŠ{H‘ÇÔsUÿ á­<!ÿ @íoşüCÿ ÇkmO¹§²ŸcÛ(ª±ˆt[RÙdK{Èâ5” áYA €HÎ­_­ö2ßP¢Š)ÁEP
(¢€°QE‚ŠJñıWöğæ‡{%¥ş‰âKˆÎsZFõÁ”TJ¤a¤™Q§)l{âğÖÿ  v·ÿ ~!ÿ ãµêŞñ%·‹ü=c¬YÇ,V·‰æF“€‘ÈŞ´FqŸÂÂPqÜÖ¢Š*É°QE‚Š( ,QE`¢Š(Q@X(ªz¶¯e éó_jQYÙÂ7<Ó6ÕçÒ¼{Tı©4†ÔÇÃÚ&¡â†m¨yBO÷â¢³•HÅÙ½K9KT{móßü5 ²Ô×Tğ…ÕƒFÛeOµæT=ÁFyö$W£¯Æß
Ià™|N—¬l£aA·‰HÈf~ñúã‚sIVƒNIè†éM;Xïh¯Ÿ,ÿ kí>]IcºğåÅ½‰lˆî„’ë³h†ê÷'U´×4ÛkûÖæÒá‘JMT'ŠñdÊ.Ì·Ey_Š¿hÁº´ö¦‡â$FdkDXæÚq¹H7)ê¡¬økOĞ;[ÿ ¿ÿ ñÚmO¹^Æ}l¢¹‡¿´ï‰:,š¦™ÔñÌĞ»EWÜ 9Â³|Ã½tõ·™‚Š( ,QE`¢Š(Q@X(ªÚ–¥k£ØÍy}q¥¬+ºI¦`ª£Üšñİgö¤Ğã¿Z“â–}‰å)d?ìğXÿ ß"³•HÅÙ²ãNRÕ×E|÷'ídúv¢Öš¯ƒn´ù#m²Æ×½Oª4kÏÔŠõÿ üBÑ¾"i&ûGœ!Û4®ÙaoFÔ­©Ÿ©ÊÒÑE¡6
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¬¿k¿ğéßÿ gßêYQö]6:vÉ…ÈÎ3“í^]}ûSxkL¹{{ÍÄ—	÷¡Ö$uúƒ.Eg*‹³eªr–©ËExŸü5§„?è­ÿ ßˆøí{ªE®iZŒ
éÜ):,€
ÊóÏ­Teü,™EÇF\¢Š*…`¢Š(Q@X(¢ŠÁEP
(¢€°QX>2ñt~ÒF¡.™©j‘ù=23 Á;˜0£}Åyü5§„MÖÿ ïÄ?üv³uafËTå%tl¢¼‹Aı§<-â-jÇK¶°ÕÒâòe‚6–‚c€I}+×jã%%x’ââìÂŠ(¦+Q@X(¢ŠÁEP
(¢€°QE‚Š( ,W	ñCâö‘ğ¾Öµ£Şj ˜l¡ 1ø˜Ÿº¹ã<ûµÅx/ö¦Ò<I¬Ã§êš\š'Á#¸û@š=Ç sµv~k5V\—Ô¿g.^kh{…QZ`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP<×öšÿ “nø¯ÿ b­ÿ ¤r×óó_Ğ7í5ÿ &İñ_şÅ=[ÿ Hå¯çæ¼üVèê£³>Şÿ ‚Eÿ ÉÉx“şÅ+Ÿı,³¯×=µùÿ ‰ÿ “’ñ'ıŠW?úYg_®u¾à2«ñ	¶´´WI›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &ÚÁñö“o­x/[´¹‰e‰ìå;\g•?P@?…oÖo‰ÿ ä[Õ¿ëÒoı ÖU¿‡/FiNüñõ?;ëÔ¿fŸù+šgıqŸÿ E5ymz—ìÑÿ %wLÿ ®3ÿ è¦¯7üTwVş>¯ñgÃıÆÖr[êÚlÆáÌ‡Õ\rù5ñÄOÍàj,®eX1JF7ÆÃ*ß\}Á¯¿ëä¯ÚÒŞ8ş!iò.É§!^$påü«l\Tyd»ØË'¬YÈ|(ø¹ª|5Õ¢Y.tYı&Å›+ƒÕĞ×5öæŸ{oªXÛŞZÈ&¶¸eŠEèÊÃ şF¿9kíßÙîê[¯„:LK2¬±‚º²¸ ¶nQq}ñQjK©äµîŸ:ç‡ïM5¼‘»Ë`Fï£_?WÑ_¶ü~øgş¹ÏüÒ¾u¯:§Ç/VwCà¡úà…ÿ Š/@ÿ °}¿ş‹ZÚÛXŞÿ ‘/@ÿ °}¿ş‹ZÚ¯~_<hü(M´m¥¢ ¡6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Õ}FÊ-CO¹¶™CÃ4M«‚¤EY¦Éş­¾†³«ü9z1«¦š?8dP²0"½»öIÿ ‘ûSÿ °sèÄ¯›ıtŸïç^Ûû$ÿ Èı©ÿ Ø9¿ôbW›…ş*ùşLôq_õ_™õÚ6ÒÑ^±æ‰¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m-ó—íc·ğÕàUoš"Øä®ŸcŸÎ¾j¯§lùøkş»Íÿ  ­|Å^%]*Íyşˆõiÿ ?×Sï¿†ÿ ãÂÿ ö·ÿ Ñk]6Úæ¾É8ğ¿ıƒmÿ ôZ×M^ô¾&yøP›hÛKEABm£m- ›hÛKY¾"ñ&á=&}OUºKK8FYÜõ=”äö“j*ìi6ìŠ^9ñ…€ü3y¬_7îá\Gy–C÷P{“údö®kà_Œ5?xµ]Zešéï%Aµ #
 ¦zkåÿ ‹ß/~(k‚BÛI¶$ZZg şûz±ı:{Ÿ¡¿e¯ù%iÿ _³ìµËF£«RO¥´ûÑÑRŸ³‚]oú3×v×€~×¶Qè{—O|s´¦HüÔWĞàßµçüŠ:ııÔb¿„ş_š?ñÏò>V¯º>ø´şÿ ¯oı™«ázû§àwü’×·şÎÕi|¿Rñãmih®óŒM´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€>2ı¡>#\øÃÆWZdS0Ò4ÉZ¢Sò¼‹ÃÈ}Nr°÷5¯û&ÛÇ/Ä+÷tñéîQˆårè?
òq‹ëZƒ1É7àF½öIÿ ‘ûSÿ °sèÄ¯#	'*ŠOwÈôq+–›K¥¿3œı£íî¡øµ«=Ä-J±4,G‚5ïÈ#ğ®?Kğ·«øoSÕ¬í&›J±t7¹Û“œwÚ>½ëï]_ÃºWˆ4Õ4Ë=I#9E¼·IBı7Š³e§Úé¶©kim­²,0Æz 8ºÂo©—Ö4VGç*«HÁUK3$šû—à‡†ïü/ğÏH±Ô•¢»ÃÊĞ·XÃ¹`§Ğàôíšé-|áûï¶ÛhZm½îw}¢+8ÖLúî5±]4i{$õÕ˜Õ©í-n‡—~Òšt7Ÿ	uIdPÏk$2ÆÄrÌUãğc_WÛ´Gü‘ıéşJø’¼Úÿ Å—õĞí¡ü5êÏ®?dñŸ‡ö“ÿ @Jöµâÿ ²ü“{¿ûÉÿ  G^Ó^Ä~ú/Èó_Äı_æ&Ú6ÒÑLBm£m- ›hÛKE &Ú6ÒÑ@~Ò®¼Kâû
	™t-ü¿-N“¼Íëƒ•˜>´ßÙfæø Ğ3Ge+!#îœ¨Èü	üëÎ|hÆOk¬Ç$ßÏ“ÿ m½+öUÿ ’œÿ õá/óJòp’rª¤÷wüKZn+§ù•ië{˜¾*\É4-2[ÃäHG(\|6Gá[¿²^›¨IãJö5‘tÄ´1LøùË)UúğOÿ ®¾Õ´3Ä¤Z¦i©D‡r¥Ü*ƒêSiúm¦“j¶Ö6°Ù[/İ†Ş1 ÙJ‡³¨ç~ÿ ‰ËR¯<mÛğ'ÛFÚZ+¨çmih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥ª:æ¹aá½.}GSºÎÎÜòÈp°õ'°šM¨«±¤Û²*x»Å>ğõæ±¨I²Şİ7mîíü(=Éâ¸ï€ş7Ôş xVÿ TÕeW—ûBD‰DqíR¨0q“Éæ¾møÍñ~ëâv®©kmÕÙ­ÉåÏO1ÿ Ú>‡¤ûŸìŸÿ $ŞïşÂ2è×-®­I>–ıQ½J~ÎïÑÑ¶´´WYÎ&Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Úù§öÀ³†;¯İ*<‰q>9*¥cù×Òõóí‰÷|+õ¹ÿ ÚUÅ‹ø©Ó‡ş'Ş|Û_ ?ÇüP>ÿ °u¿ş‹Züş¯Ğ/‡ò øsşÁÖÿ ú-hÂ|2ù~£ÄüQùş†şÚ6ÒÑ]§(›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Úø{ãæŸ›ñkÄ@¡¤Iv¨ÀÜñ«ÔšûŠ¾&ı£?ä°k¿H?ôJW‹Ş?3·¼—‘Îü/ÿ ’áŸûÁÿ ¡Šûëm|ğ»şJG†ì#ş†+ïªèÂÿ æÿ $c_øŸ%ú‰¶´´WI€›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›k‰ø¹ñßá·„§¾Ü­¨Í˜¬¡<î“Hô^§ğë_Ç:Ò~è’jZ´â4B¼É3uGsúõñ/ÄOˆZÄK©ßŸ-Éol§)}”{úæ¸±¹W$wüª4¹Ÿ4¶>ÎøW­^ø“áî‰©ê3}¢úæòË±WqÜFp Ó°®¯mpÿ ¿ä“økş½¿öv®æ»Säÿ Ú»Ãz…·Œ­µ¦äÓ.-Ò”r¨ëœ¡ôÎr=r}xï‡ô+ïëºn\^\8DE©ô©=«ô2òÎßP·{{¨#¹·aâ™«B©iÑ¼>ÎÚ^“c¦´ƒlí’"ß] f¸>«ïó7¥îuûwËkj[ÓíZÒÆÚ	ÌxãTgşñ j}´´W{wwgVVmih¤16Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih 4ı¦—ş1·â¿ıŠZ·ş‘Ë_ÏÅA´×ü›oÅûµoı#–¿Ÿzóñ[£ªÌûş	ÿ '%âOû®ô²Î¿\ëò7ş	ÿ '%âOû®ô²Î¿]kl?ÀgSâE>Šé¹Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE•›âùõoúô›ÿ @5«Y'ÿ ‘oVÿ ¯I¿ôYU»—£.ŸÆ½OÎÚõ/Ù£şJî™ÿ \gÿ ÑM^[]§ÁûÏØxêÎoØÛêZÊÇ ŠŞé‚£)C¸’]z÷¯2‹å¨™èUWƒGİÕñ/íâØ<]ñ*ö[IÖ–h¶qÈ§!¶ä±Û‹~Uê(ã§Œ,^ÊM.ÓG¶”l•l.!BêzåŒ¬À}®OCı’üQ}"NşÃL„ıí¬ÓH>Š Sÿ }WEnzÍ%….Zi¶õ<sEÑ¯<C«Zé¶µÅåËˆãGR êOµ}õàŸÇàï	éz4m¼YÀ±³ÿ yº±üX“X¾èág±®µl—÷2Ü.8Uö‰5İWUj”mÕ˜T©í%~ˆù“öÂÿ ßÿ ×9ÿ šWÎµôgí‡ÿ ¾ÿ ®sÿ 4¯œëÉŸÇ/Vz1ø#è~„xşD½şÁöÿ ú-j×ˆ5ëh÷Z¦¥:ÛÙÛ&ù$oĞÜ“À¦«øşD½şÁöÿ ú-kçÚËÆ’İkV>‚B-­P\Ü(<4÷Aú/?ğ*õñ5}Úİ³Í¡idö9ßüyñOÄmj=/A’m"ÆâQÖ¯¶i‰8ÜsÎz×­wÒşÉğ?‡Úf×¯$ñ.Íşcò˜ÎÜcv3ÆíŞøí^[û8é±ê?ô0¶ë,àUCÔƒøWÛ5ÏJš>iêÙµJåÈøÏáïÇ_|:Ö?³õ‰§Õ4È¤0Ïkråå‡ËcÈ#tœ}:××Ú>­i¯ivºŒËqis’)¡Sı}«á¯Œ–ëkñKÄèƒöéşú9?©¯xı’|O%ÿ †u]gİö	–Xrz$™ÈúRàTajÊk–{†"š‹æ‰ë~5ñ–›à?Üjúœ…`åT^^W=G©ÿ _:é_|añÃÆÑhV:Àğ¦(gÛhO˜F~ÿ ÍÀ¨ö­¿ÛoŒ`û44»ˆû¾oËŒûíÎ?ùãAÔ/´½jÆëMy#¿Šeh,îİ õÏL{Örªå[–[^Ö44©İnÏ¦üIû*Xj:|²Zø‡T›WÛ•—Q‘eÛĞ€¡†}rqï_:i^(ñuIOÔî´ë›y
<qJvnS‚~ë„_ VÍ#ÛÄÒ®ÉJê;69ùùãïù<Cÿ aıÕˆû)®M7ü,:ö|ŞGÖ¾1/Äİ6k[äH5Ë5eŒa&CÇ˜£·<Ø‘ë^£_şËWÅHÖs!Ç| ßÒ¾Æ9ÁÇZô!>jjOsqå›Š>føÿ ñÛQµÖn<5áË·²Kc²îö‰NèÕ@îG9Èè9gÃ_ÙÎ×Ç·ñˆ5›óu~¾l+lë•\œgX¼càšãM&µ¨5Æï=®$2në»qÎõ_‚ÿ î>Ç«Æ÷ºï‘“™m²rvÿ ysü?ˆô>uÆ£½]ÙİV2‚µ>†§Ä#âÀÛF];Ä×—¹ıÊÌpÆGİ!·yg®HÎ;+Õ?f†¸øWo¹™™nç·ûÙşµÛÜa|SğmÌ6÷0ê:UüF3$G;N2:«)ÁÁä\Ÿì÷á½GÂ>Ô4JÚ[y­õ9•ZD*%\(¹ê§q]”¢éÕ’{5§ŞIÉJšïÑšÌN ä“_3|`ı¥nšòãGğŒÂ#&9u@4‡¡ç€?Úê{c©ô/Ú[Æ’øWáûZZÈc¼Õdû0e8+3!ü°¿ğ*øóO¶ûf¡moœy²¬yú+
Õ%9û(èoJ
1ö’>„ğì÷/Ä/E¯ø»[ÔÚêù|Ø#YC8C÷YÙÃg=p1Á×¨ë^.ıŸüm>•mªË=´$HÌKAqä„ü§±Ç ƒÍ}›kn–v°Á…Š$ª; 0|µû_[ªø»C˜ìYOá!Çş„i×Š£iRĞš2u.ªj{ßÃ?ˆ¶ü6š• ògCåÜÚ±ËC&:{ƒÔÿ u•ñ÷ì»ây4_ˆË§Å¶©DË7¨.§ëÁğ*û»a>x)Ò$œFQO¢ªä‘³RÌBªŒ’x¾`ø½ûJ_\ßO¤øJ²ÙÆJI©(ICåçî¯¿SÛè_´ç%ğ¿“N´Çw«¹ƒrœˆÉ®U~Œkä
ÕoµÍ:Ùşä×ÆßBÀZóêÔ•JÊ.Çe8(ÃÚHúÂ¿³„ş6ğ¼Z¿‰õíHë±yĞ«Iæ_3~K0Hc8¯µñˆşë×VÖ:¥Õ…Õœí‹‡adb*xa‘ĞŠûıcEUUF …|_ûBøPğ×u-OìÒ6•¨Ëöˆ®UI@ÍË)=ˆlğ{b§g%*z/êÃ£'Q8ÏSßşü`_‰Ú\Ö÷±¤İšƒ2ÇÂJ§"ÜğGcZôÙ?Õ·Ğ×Ëÿ ²„õEñ-î¼ğI–¶­n²º&ve8_P6ò~•õŸêÛèk®Rr£ynÓ9$”fÒØüß›ıtŸïç^…ğOâŸÃ]gVÕ.â{™ÄÃooYd.„öO·zóÙ¿×Işñşuë³‡,5ÿ ˆ%üãì6Íu
7İPøÉ#ßçaÓsJ.ßğÇ¥ˆi&åßõâ$?<]¥M¯kV7ÖÚ<Ák$õ1İ€?‰‡ã^máÿ jŞ½K½'P¸°NwBäìÃ£cÅ~‚k–ëu¢ê¸²[È„7B
‘_œì6±”êÇÙT÷Y4ß´‡¼¸¾
üNÿ …›áSs:,Z£ˆnÑ>é8Èp;½A¯@¯˜cë¶]sÄV¹;Ş)1î¬GşÍ_OM2[Ã$²6ØãRÌÇ°$×«Ş
oª<ùG–n(à>.|]°ø]¥!d[ÍZà³YîÇüı~½r>}ğ“x¿ö†ñd¶º–·q™
ù×+	+KœXÁÁcĞ“Á'8®â‹§ñÏŒ5-bvb³JD*OÜˆpŠ>ƒkè¿ÙMH|«ß`y·ŞY=ö¢)«šà¥'ˆ›sÙkc®¢ö0´w}N/â§Àù¾é1øÂºÆ¢±Àê·¥,y8ñœc¸®—à_íq¯_AáßJ¯y/Éi€¾cvGÇ=~‡kÕ¾0[­×Âÿ Æã#ì27â£pıE|òZÍÑ9XØ:2ğA Šı…k/‡°({jW{Ÿ£õÈ|Mø•§|3ğù¿¼}Ì„¥­¢¶gşŠ;ßRjøÄğ•x7FÕß»µIıüa¿Pkã_Ş8—Ç5	Ä…¬mÚÚ®x§‡ûÇ'ñ•¾&«§îÇvcBÓŞ{#£ĞüEãOÚÆk¥Ë¬Í§éìY£´%!† ºÎr@‰äÕŸŸ`øc¤ÙkZ>£wqJ!œ\²ïW ÊÊ;qÍu?±åŒ{|MxFeAôS¼Ÿä?*õŸ^¸ñÇÃ½OM³_2ùvÏgvıHÈükR^Ç™|[ÿ _#XÔ~×—¦ß×Ìù£á—í¯ø/P†JîmcEf$7^H×ûÑ±äcû§ƒíÖ¾ÆÓï Õ,mï-d[\F²Å"ôea#_K¡êMª0X\G—ö_)¼İŞ›qœ×Ş?ô+¿øBÓ/¿ãîÚÕVUÎv·]¿†qøVØiÊP|İëÅFJÇşØòğ×ıw›ÿ AZùŠ¾Ÿı°¿äá¯úø›ÿ AZù‚¼ú¿ÅŸ¯èÊ_ÃõÕŸ~|0ÿ ’qáûÛÿ èµ¬¿‰ï¼V¦óOñn­áËˆ *Òå’İ±“¹ÕH9÷ÏAZ¿_†Ş$à2Ü“ÿ lÖ¼¯âoÄKÿ ‰º£øÀ‡í&LCRFÄK8eıßSß ÎkÕÄ´ï«{mÙ>‡ŠxS\ø‡ã]v-'Iñ¹srç–şÑ˜"(êìwp£üó_R|=øeªxNæ;İ[Æ:Ö¿wå•k{‹—6 û’Hìr>•¥ğÃá›ğÇA[+@&¼—uxË†™¿¢Ãú“]*tÕ5«».¥G7¢²E>ŠÚæC+åOñÜÖsøÄš®›>o Hl¬¥“lªÈzd“šúº¼§öœÿ ’Gÿ _èÁ\˜”œ9»8vùì|a^Åğ·á¿ÄxUo¼7âŸì7ÎtßÚ|ã;QHçkÇkìoÙgşIZ×ìßû-saà§'~ßª7­'­Ü³ğ‡À>6ğ«?Š|Gıµm4!!íÓÏ±÷dœH ;ŠækÏùt?úş?ú-«ß+Á?kïù4?úş?ú-«£îĞå^_™åVş¿‘ò­{ï‡¾,x‚ïÀúƒ¼	a5Ö±®.ïòrÇ…İò£ænxçšğ*ûcöwğí†‹ğÇKºµ€%Î †{™O,í¸Ï õõ¬0ñrRW²4­%;]Ÿ&xóCñFƒ­lñR^.£ ó—Sy¥Ç¨pH?âº¯„ß5êÖÖ÷·s_è.Á&¶™‹˜”ÿ dô#®:Ö½CöÂ·VÑü7>Ñ½g™7wÁU8ı+æ
ÅJTj4ºòª°Mõ?H!•'‰$ƒÆêYzy•˜"–bTd“ÀÎü1¼kÿ ‡~‰,Ú|'Ô ùWûNxÒ_øtûYwZ´†ÊpD@fOÏ!~ŒkÖ­QQMö<êQuGŸ|]ı¥¯./'Ò|#7Ù­c%$ÔÔò‡ËÏİ_öºØ«~ı_Ç~‡]ñn·©›ëäó¡E3"Uœ1$õÀÇ^µóæf5bÆÔğ'8ü	€şµú-)o
E…DPª£  `
ä£hœêjtÕ—³j0Ğø¾ëÄ^/øã[&RK‹kv¶óÖ÷AØOÊHôäFkê‡tÿ ‰‹T²ıÔ€ùwÌrĞÈ*}Gp{Ê¾{ı®­Ò?i†“O¿	Î¨şÊş'“Iø„úQsömRR¹ãÌ@]Oä~4°õ›§/;|‚´Š¨®è§Ñ]÷9WüPøkñÅ+’ûÃ*şÉÓ$E¶şÑ¸‡Ìv¢‘Í{=œâ¦¬ËŒœv?8µym5¨'2h¥dwÌ	äõæ»„ñ/‹5ë«_êÿ Ø×±Û™$›í2Á¹7(Ûº0Iä:q\Ï‰¿äcÕëî_ı×®~É2mø‰¨&>öšçò’?ñ¯+i(§Ûô=Dœ"Úïúíğ§Áş'ğÇ‡õ?km«]Ï)h®îYš4*H9<W‹üpĞ|Qğµl/4ÿ xŠóO»v‹mÆ¡/™‘’ÏaÒ¾©¯
ı®”Â¤r5ƒÿ lŞ»q:Sæ[«‘ÍA·;=Ïœÿ áfxÃş†½oÿ 3ñUögÁËûSá‡î¯.&»º–ßt“NåİÎæä±ä×Áõ÷_Àßù$¾ÿ ¯oıªp­ÚWòıGˆVµŒÿ Ú#şHş¿ô‡ÿ G%|I_nşÑ_òGuÿ ¤?ú::øŠ¹kÿ _×CzÃù¿ĞúçöOÿ ’owÿ a?ôëÚkÅÿ dßù&÷Ÿö“ÿ @Jİı¢<i'ƒ¾\‹YWº‹‹8™N
‚	vğGâ+Ó•OgIKÉ~HàŒ\æâ»¿Ìó¿ŒŸ´”öw×'„¥T1ú¦eºˆ8ş÷åëTş|¸ø› §ˆ|_­ên×™khÖ`ÒììÎƒØÓ¿5óÔQ™eDY€üëôcG±KÒl¬âŠŞ‰@áTü«šŒ}ªs©©ÑZ^Í¨ÃCã­zëÅŸ³ßŸN°Ö&’Ğ4+&LïFNàƒxà×Ó_
ş&Ù|NğïÛ`_³ŞÂDwv¤äÆøê=Tö?_JñÛ
İWWğÜà|íÈ~”ı×û5øOüN³µÜE¶¦k"ö'ıw ?KQÆ£¤İ×OĞu œ=¢Üû:Š}ßsŒeyÅo‡^>ñW‰£¼ğÏ‰ÿ ±ôñn±›íˆ2à±-¶5#¡õâ½†ŠÎqSVeÆN;íÎ­ê·’ù÷pÜIÒ†-½ÃÍ“ÉÉÏ&º?…^ñŠ¼NÖ^Õ±õ!ä7h’•İ'œ=«3Çßò<x‡şÂú1«Ğ¿ey6üR½e0ÿ ĞMyXXª’Š}¿CÑ®Ü·õ©ôÁÿ ø¯Â:Œ^*Öÿ ¶æšUhírÜyjf@1Ï¥vZö½eáëTÔfövÉ¾Iò¤œ =MiWÌ?µ§¥›SÓü/…`…İÈï;d ?A“ÿ ¯FµOe7Ù4áí'©ÇüFı¢<GãÉaÓ.fĞôœâ8m_l®=]Ç<ú¯Zô	şËÚ¶…÷ˆõAuk˜Ä»-Ùq	# 1`Kßïõ¯Ÿ<'f5/èö„nŞC¹p+ôB±¡MN.sÕšÖ›„”a¡ñn±®xÓà/Œ&Ò-õ¹æ·‡l‘G1/Ñ‡c·¸8äy¯£>üb°ø£§ºybËY·Pn-3GMè{®ÓĞŸ"ı°-5ÿ ]ù¤¶’2ßî¸#ÿ B5ãßü[qàiºÄ@‚Qæ¨?~3Ã©úŒş•j¸O’OKØº”Ô¡Ï©÷íÔj¶–2H¼Ä)æDÛ]r1{_)ünÒüWğ®ûOk?ø‚óO¾°Ï¨J$F\d†àw¯¬ ™.ahÛ|r(uaÜkçoÛş<¼1ÿ ].?’V˜¯v<İne‡w•ºÿ 3Æô5ëø1›ÿ Š¯¡~!üxoø7DÓ´ÙEß‰.tè$–yNÿ ³îNöÏŞsÔõ>ÿ +W×³G€m¬ü.(¾ˆ\êÚ‰o*i†æŠ;@\ôÎ>˜•y©E;yšÕåƒRjçË×ˆµO^µŞ§¨\_Ü“Ÿ2yôôÂ¾«ı—|]¨x“Á—–š„ò]68Š)¥;›ËeÈR{àƒø;V/ÅÙ~Okòjº³ÓÖäîÎërF­İª·_îã¯~ÃÓ¾ü6‡á‡…†š³‹«¹d3\ÜÀg è ó=ë\4'.c:óŒ’å<GãÆ‡â†¿dÔ´ïø‚âÂöfŒÃ5ü¢ln `Æ{b¼‡şgŒ?èkÖÿ ğc7ÿ _F~×_ò#é?öú-ëäúâ¨Üg(§¢ÿ $uCŞ‚“İÿ ™ô.—ñ»ÄWŞĞü5áxîµ¯Ë»é33Å—lrÙÉÁfàqøyOÄOøËC¾ŠoÅ{çÜd¤÷S‰Ãc¨†G¦x¯¨fïéÚ7Ã]?Q¶€íH4··,ØvP öPOsU?j{t—ák;(/ì,‡ĞÃù]Ué¾^y;¿Àç£Qsr¥¡ówÃÿ ‹^ ø{¨C%ì³ØnŸ3–‰×¾û§Üsõé_Bøëö„‰lôû[¶³¯jP¬‘¢¡qnd£«ÿ ³Ğc'Ğü_ZşÊ¾Óm<úêÀT»šHväª)ÀUôÏ¯à( çQ8_E¯˜ë(Á©ØğŸ‰^ø…kz§Œb¿xgl,³N²F¤ó·
Ä'Ó°|#ãí{À÷ÑÜéŒÖáXƒq1H=:_aü}·ãáˆ„ŠdHëÇB$\øn¹§z5-mªBòGèÃïAãïØkP'•ö…"H³Ÿ.@pËùËĞK*AÉ#¬q ,ÌÇ Ô“^+û$Ş¾jŸ'Pl{ˆk¶øŞ·Ò|)ñéáÚàÛò#ëåîgş9º½Z•9ióÛ¥ÿ Ï„y§ÉçoÄño‰_´æ¡¨jRit²³åÿ iH É/8Ê†E÷#=ø®ÍfÔ×ìOø¿ZÕµ\™pbRG@1#ñ…|“_mşÎú–¥ª|+Ó$ÔšGti"†I3¹¢S…ç¾9Ø
ä¡jÜŞÓVtV½;rh–~ x?Uø?ãÓâÔe±g·¼¶f‰	  ½«Ö>~ĞW÷Zµ¯‡¼Oqöµ¹a®¡'úÅsÀG?Ä@İryÏloÚëşG­'şÁãÿ F=x®•3[ê–r©!£™×!¬(Ôp©ËÒöülmR*p¿[£óÇßxóP±¾ñ»ªé’hÖ.6’Ê6+8E!J [æ$úãÒ¾ VÜ úŒ×œ~Ñ_òGuÿ ¤?ú::ëÅ+Âı\<½ô—Sâ*õÏ„ÿ ü{âÏÍyáŸcéëpÑµ¿Ûî Ë…R[ljGB9ëÅy}uû&ÿ É7¼ÿ °ŒŸúW6
riöıQÑ^N)[¿èËÿ 	>øëÂzõÕ×Š<Kı³e%¹8~ß<û_rÛdP ó×šõš}éÇİI#»¶ÆQO¢ÀòŸŒß­~Æ¶q%ş»2oX\Ÿ.=ñÉÏe˜ã<ÃøãÆŸu¬x‡ÆöÚrÎ`~šâ,ÉåæÊ’}kÈş4-òüRñ'ö€q1»b›ûÇÿ ,ñí³mw²¥©Cãë‹+v‘´é­]î“’€®67±ÉÇâkÍ£SÛOßëÓ±ÛRÊé'ÆoÙî? è¯hÚÍİœ,«qáS*n  2@Æ^õÆxãg‰¼	{Kùµ8$°»ºîPûÄú“ö€ÿ ’CâúçşJøv²œ­CCHZ­?|ığŸ‰ì¼eáë-cObÖ·I¸ûÊz>à‚?
­ão
ÍâíY[ëZƒ2È$[­6cä6¶:¯9ÆG W~Êwo†2ÆI>Uüª>…Pÿ 2kÙkÕ•§~¶ÿ 3ÏMÅéĞøOÅ^,ñ¯…|I©iø»Z–[)Ú"ê Ø<nî9«ßş,kzt«oÄú´ºT2¸kÙdFP§‚¥¹çVwÆÏù*Ş&ÿ ¯¶şB²¾xr/xÛFÑçv»…Izíê@÷À5äP”Üág­ÑéÔQQ•Ö‡¬ø³Åş1Úİ^øL¾²ğÄ{‚Gm ¦QÔ“ÒöW täó^k}y¦^	íî&´º²$Ê:‘î9¿Ell`Ólà´µ‰`¶…qÄƒ
ª ğGÅRßâ?‰£B¢ê3áGAó“W^>ÎjÎ÷3£.x´Õ£?gOŒW4}[—ÏÕ-cóaºo½4`€Cz°Èç¸>Ù>ß_~Í·&ßâö ãÍI£>ÿ ºcı+íjô©MÊšosŠ¤yfâ7âwÄ‹†~mBí|û™—kj§WÇ¯`:“ş"¾Añ?ÄÏ|KÕ’›ùØ\Ê#‡Oµc ±Â¨Py9=['Ş·?hïIâo‰W¶Ë!k=/ı%Ï‡úÃõİ‘ÿ GöÓWTø· £®ä‰Ş|{¢3Ô
óÔŞ"ªMû¿¡ÛÊ¨Órê{O‡ÿ f7ğş—ÖŸâKMñ* o>İ‚Û†Çİ*â½¾÷>«Ê>7øãSñzNâO²ø‡E’xnÙ#˜0d‹õ
O§q×í
ù{ö¾Ñm­uT6ÜİÅ,2>ğŒ©R}şr?Zbb£Ëµö2¡.ië¹óİ~|<ÿ ‘ÃŸö·ÿ Ñk_Ÿµú	ğïşDØ:ßÿ E­^á—Ëõ'ñGçúÔSè®ÛœÃ+Ì>:ü\?ôH`°úİğa q•‰GYïÉÀÏÒ½J¾<ıªvø ]ŞRØÄ"ÏM¹lãñÍrbj8Á%Õ"¥'~ğ¯Àz—Ç_P½ñ½{%ß1¼Íò»6p©»*€`öôôíüMğGÄŸí$Õüâ-AÒÜo—O•ÁvP9 ²O÷Jı3Ò¼Oá·Ä­Sá¹öı<¬°È\ÚH~I—ÓØÇ·Ó"¾Çøwñ[Aø•b$Ón<«Õ\Ía1XıN?‰Ú}:TSŒ*BĞÒEÔ”á+ËX'û4x«SñWÄÍvûW¼’ööãOËHà»"` 0 è+éªñÿ ü;›Á5û+)F‰©éòJ%•¥Ğ²Œ•$Cí^£â]røQÕgæ+8b=v‚qø?Şœ¹(§.—¿ŞÌfœªµ¶·Ü7øÕñÊßá¼cMÓÒ;Í~TÜRİOF|u'²ş'ßÇşøWÄ_´­{yâ-zûû"Õ€£à3Ÿà>êñÉ8ôãšòoXºñ¯y©^Èeºº•¥‘©?Ê¾ÆıštÔ°øK¦È oº–iß»ÊÑEsQ½y9TÙt7«j1J³Åş*|4Ôş\YkÖµÓ§(¿™¶HäÆ@b Û±½Gà_ÇCãïø’ë["×cBÑÌ€*İ(ëÇfHÈ­OÚjİgøK¨3b_®ğ?‘5ò†õËkÚ~«jÅg³™f\wÁä}ãñ¥Æ«…ıßÈnÒš—SôF¾&ı£?ä°k¿H?ôJWÚ¶7i¨Y[İEÌSF²/ĞŒç_şÑ¿òXµß¤ú%)ã7Š~dá]ÜŸ—êwáoü”ÿ ØFıWßUğ/Âßù)ÿ °Œú¯¿k|7ğ¾oòFUÿ ‰ò_©Ç|Kø•¦ü3ĞMı÷ï®$%-­°ó?ô¹íù
ùŸOñG>?xÊ-!µY,l¥ËÉ©)¤¨?9èây#¥`üpñ¤5ø‡©Mæ²´sil¹à"	ïŸÆ½Kö=ÓQ¦ñ. @2*Ãú€w1şCò®zsúÅO{áìo8ûz|_×äWøû:¯|5.¿á_QkËógY\*>ó# R1×<gšgÁ_Ú*ö=FÛCñUÉº¶™„pj2¬ H‰O÷#¾GO£¼Kn·Õ`q”’ÖT?B„Wçe“¡WİÛ°ãmOŞß¹úEX8ñ¶›àO«joˆ£ùc‰q¾W=}Ïè2k'à¯‰äñgÃ=ögßr‘y˜õ-+“î@ñ¯›¿ioKâoM¥Ç!şÏÒpˆ^²7×?/üºq}š´wÕÎz0öŞè:_‰6øßãmËP“IµºªÚÙ¹D02Ìì0_ g3Ğ
Üø¹û<ZøÂ¯iš¥åäö¬¿j[~`Ä.äÀ`‘ÁÏ¯ìcŞ8ÕnXeà°!=·:äşCõ¯¤~ xm¼aà½cGB[«vHËtŞ9\ş W?³æ¡Ï¼¿ğ?oij¼»#ãß‡ß|KàKè³{6§¦–7r]¿ì’‡éÇ¨5ö6“âñW…bÕ´GCwnd¶ûFB‡ÁÀ|r0ÜzøRğş¥¤j¯¦ŞXÜ[ß«ì6ïŞO°ïøu¯µ~xgPğŸÃ=2ÇS »c$ÍıèÃ1!HìqÔvÍi†”§&ô"¼TZksçOñ“,÷‹uk-I®¥0D–’9p[j«"…^;~>µäµõWí}ÿ "†‡ÿ _Çÿ Eµ|«^}D£RQ]?É°|Ğ‹}Ì÷Oü)ø™®xOM¾Ñ¼aıŸ¦Oè-¿´î£ò×$cj©QÎzW¬xá‹-|/¬i~(ñeô·RÆğ^éÚŒÍ4*½@yã u¿ğ7şI/†¿ëÛÿ gjî«×T£Ï5Ô“±ñoÅKÏ|5ñtÚ?ü&šåä>ZÍÆşe,ê7õøW1¦üIñtš…ª·Šu¦V•AS¨ÌAíW¡~ÖÊâ&@Á:ddûşòJñ½+şB–õÙ?ô!^e9?ikõıNù%ìïnŸ¡ú1\÷<q§|=ğìú¶¤Çbü‘B¿~iD_Ëğ šé+ãoÚcÆÒø“âÚbH~Á¤!<õÎş]õë:q´wg|ò×dfø£ãŒ~'kqXÁy-ŒR¬0iö2Ôî  ÄrıºñèzÕ·ì“§dfã_¾şÛ+Ÿ>0¾J¾?ºFâ=÷Â¼oàˆ¿ø¹áÔa“4¿÷Â3Ô
û–³£N3§Í=[.­GòÇd|S¦|Rñ¿Â\é7„—Ée1†[+×2ÆÀtTÈÁGõOÃˆÚgÄ­u<˜¥BâÕÎ^ÇCêcßóåÿ ÚzÑm¾,^8]¾}´2s·nñÚÏıŸüi'„>#iêdÛe¨°³¸Rxù¿ÇàMF¬¹½œİú|Ë­MrûH«u>ÈñŒŞ Ñn´õ¿»Óe ]ØËåÍ9VíÓBkãï‰ú‡Œ¾øÂëE>3×.ãEYb˜ßÌ¥‘†FFş£‘øWÚÕñßíSÿ %KşÜaşmO¥¤·¿èÅ‡Öéÿ [n“ñ7ÅKªY™üU¬ùd2oÔ&+·pÎ~n˜¯`ñoÄ|^º¼±øeu‡lJ=ì.!’sşû´z(9Ç^¸;ØÚ›ëë{`ÛÒ,{½2@Í~ƒøOÃ^ğı–‘§Æ#¶¶Œ.qË·v>ääšTbêA©=VJœ•–§çö¡¡¥êÓÇyö‹}F
È$$H®9=s_@~Î¿õ+Íj/ë×r_%À"Îêfİ"¸ØÌy €qAã¿'í/
CñsS(¡wÅ¶;Ÿ,ş•È|6º6_¼76í»u2}¼ÀéYa¥(ÔQó·ãcJÑR¦ååÔûşŠ}êÜó†QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸEÌ¿i¯ù6ßŠÿ ö)jßúG->õı~ÓŸòm¿?ìRÕ¿ôZş}kÏÄîš[3îø$Güœ—‰?ìR¹ÿ ÒË:ıu¯È¯ø$?üœ—‰?ìR¹ÿ ÒË:ıu­°ÿ Oˆ(¢Šé2
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
Ìñ?ü‹z·ızMÿ  Ó¬ßÿ È·«×¤ßú¬ªÿ ^ŒÒŸÆ½OÎÊõ?ÙŸşJö™ÿ \gÿ ÑM^Y^§û3É^Ó?ëŒÿ ú)«ÎÃÿ ]ÚßÃgÚTQEzÇœQE |Íûaÿ Çï†?ëœÿ Í+ç:ú3öÄÿ ß×9ÿ šWÎuâOã—«=hüô?B¼ÿ "^ÿ `ûıµñŸÇ›§ºø·â6|ü³¬cw¢¢ü«ìßÿ È— Ø>ßÿ E­|‹ûIè¯¤üXÔ¥e";Øã¹Œà¨Sÿ )®ÜoÄŸ©É„Ù¯$sß	t[ÏøêÃM°Öî<=up²*_ZîŞ¸BØùYO8ÇZúş'Œ?è­kŸ”ßü‘_0x;ÄRxOÅZ^±ÜlîR£ø”˜~##ñ¯Ğ-'T¶×4Ë]BÊUÒæ1,R/B¤fª„!:w{¡V”£=6>~¾ı‘nu+¹n¯<o-İÌ§t“O§—w>¥ŒÙ&»ƒÿ ›áN«yı¹ı¨·PˆŒ_dòvÙ;Û=ûw¯T¬=+ÆZnµâ=WE³y&»ÓËªf5fÎw÷†9àq´iÂœ“æRœå}‹úÆ‹aâ>[JÒÛ9F”2Ÿ¯¿jó[DøWğNå5i¬-í5 7Á™%Ääú¢;¿ïp­YøóñqşhÛiû[[¾!.2!AÁçœÓ?LWÈ1µ÷‹üEºº’æşşá#k‰Ø³f OÖ¹çS÷œ´×½ÜÚıËÍè}ÅühÔæÓ|hº–‡j2ÎŠ}_R{*‚ŞıëçYË¦ø‹T´å¯&‚êHáóºVAc’y$g¯zûëÂ°ğO‡í46!¼‚Øù¤oâvõ$×Â>>ÿ ‘çÄ?ö¸ÿ ÑXâ#É(«İëúĞ—2“µ¶ıNëöaÿ ’µgÿ ^Óÿ èöm|eû0ÉZ³ÿ ¯iÿ ô
û6»hÿ ?×S’¯ñòíğfïGÕîüQ¤@ÓéwLe»1“o!ûÌG÷	ç=‰>Õá÷Ï€ü{oãÅÖ¢ëmq¦ŞÉe4g™¤…~ƒ†ÁãZóŠß³·çê	§ßœ»éçå‚SşÇ÷·İúWèYsÓÕ3®uåèùçÀ¿µ¯‡z²ßi&0Hó­ß&)—Ñ—úõ«íO†´ÿ ‰MJÌfCåÜÚ±ËC&:{ƒÔÿ \ŠøCRÓn´{éì¯­äµ»ŠI«†R;^Ùû#jÅã}VÉXık#¯mË"…?“7çZáj6ùÆxˆ$¹úš_¶Ó¶­á»~v,Ií’Ê?¥|ùlÅn""F„‡H½W£ÜWÒß¶ŠòXø{VU&8¤–ÙÎ8€eÿ ĞZ¾e®I{µû1Ö
İ­cøâù#W_‹ZáVMÿ É™¬~Ë:¯ˆdMWâæ§$cj=å£ÊTzÓœ
ô‚~6‡ÇôÙÖ@×–±­­Òg•uÉÿ x`ş5ŞW©*4ßš<èÕšßsÁ¼û.7ƒüU¦kKâ´›)„¾Oö~Íàu¼ÓŒc^óXz·Œ´İÄN‰;Èú–¦Ì †ÜB¨$»z/ÏøÜ«‚ŒchìL›“¼·
(¢¬“åOÚöéßÅÚ%¹?»Ä¸îÒô^!£Ş?V²º= %ÿ ¾Xé_BşØn­*“×´‘±À9Ü¿ŸÏùWÍõã]Â£}Sÿ ‚z‘JTÒîÒ(fKˆc–3¹C)õdSëÌÿ gßEã/‡ÖP¼¡µ5E¥ÂgæÂŒ#}
Ï¨5é•í;tØòµÙîÙ?Õ·ĞÒ4É‰:¬’gb“‚ØëŞ–Oõmô5O‚^…-ÏÍù¿×Işñşuíß²?üÚŸıƒ›ÿ F%xŒßë¤ÿ xÿ :öïÙşGíOşÁÍÿ £¼ü/ñ£ü™èbşz¯Ìú³Rÿ }×ıroäkóOõõ¯ÑÍKşA÷_õÉ¿‘¯Î9?Ö7Ö'ø‹ĞXúùû ÿ ÈÕ¯×’ÿ èÁ_BüHºk?‡ş$™3½4éÈÛ×>[WÏ_²üZ÷ıy/şŒô¯‰ô¶Öü7ªéË÷®í%€}YşµÔÓx[.Ïóg>Š½ßuúµí¿şë7ğíüº_u/Gos±í-…•Ë*òztíÖ¼RX	^9«£e=AE{oì¯ãht]è—ráÕUD,ÇÎ\à~ ‘õ¹0ü²ŸS¦¿4cuĞô;¿ÙçÅ7ö²Û]|SÖ.måR’C4r²: ƒq‚+Ÿÿ †7ÿ ©¿ÿ )Ÿıº¾“ªÚ–¥m£é÷×³-½¥º%•Î¨&»eFíŠ¬öLæ<3á™>ü5şÊ¿o}>Úf[/Êİ÷œ|¹lc8ëÚ¾f.Å˜ä“’kô3CÖ-<máˆ5d™lµIQ2l}§#‘_ kÚLÚµ¦Ü)Y­'x\İIÉŠ¿´WìtáíÈíÜúö=Ô‘gñ.OïaG¨”ÿ èCó¯¥«áo‚^6OüB°¾¸.ÆlÛ\·`üGèBŸÂ¾èF"ºÊÃ!—G­wQ—54ûiı|J±å¨üÅ¢ŠdS$é¾'Y$nS‘pGàA­LÏŸ?l/ùxkş¾&ÿ ĞV¾`¯§ÿ l/ùxkş¾&ÿ ĞV¾`¯¯ñgëú#Ô§ü8ÿ ]Yê¾0øÙªjÑ¼7¤¤Ú~•og¥ÕÏF¹‘cPÈè¾İO~+†ğOŒ/¼â[=fÁ¿{|Ñ“òÈ‡ï!ö"¾³Ó~é~8ø¡èæí^Knaš4 ÇpcÌ8ë’N}A5ñö»¢^xoX»Ó5ZËYr#zãÔ úÚ·=*ÎMëÜÊ%JVKCôÂ~(±ñ—‡ìõ:Mö×)¸èİÔûƒ‘ZõñÇìóñcş_ej3Ğõ‹ — `zÀö¯±ÇÌŠô£%R*Hâ”\+
(¢¨¯)ı§?ä‘ßÿ ×Äú0W«W•~ÓŸòHõúïşŒÍˆşù~fô?ˆ¾‘ñ}}û,ÿ É+Oúı›ÿ e¯kìoÙgşIZ×ìßû-a„ø¥éú£lOÂ½Fzıx'í}ÿ "†‡ÿ _Çÿ Eµ{İx/í}ÿ "~‡ÿ _çÿ EµkŠşù~hÏüEóü™ò¥}×ğ7şI/†¿ëÛÿ gjøR¾ëøÿ $—Ã_õíÿ ³µFi|¿R±£Î?lùü=ÿ _rèòİ}Kû`È·áïúû“ÿ @¯–«’¯ñeıtGM?áÇúêÏ½¾ÿ É0ğÇıxEü«Â¿l†oxzß‰k#L—ÇşÊ+İ~ÿ É0ğÇıxEü«ÇlB¾Õ•IL¶®qÀ'¿É¿*ô1»_ÌãÂïo#ç+ıºßlÍnŞbâdûÉÈù‡#‘Ö¾¯_€¾/`øµ­}¦ÿ äŠù&¾ïø;ãh<uà=6ñ$wkot™åeP'ëÃ~51œZ{£Jò”ZkcÌµÙ_TñÉ.«ñóS–5Ú’^Z<¬«œàœàUßşÌà¿iºÚøŸífÎ_3ÈûÍã»Í8ëèkİ«Sñ–›¤ø—JĞew}OQÑCnÚª	,ÿ İ^ÚºcN’’Zş¦¤äš{”QEldQE ~vx›şFMWş¾åÿ ĞÍzÏì—ÿ %"ûşÁ’èÈëÉ¼Mÿ #&«ÿ _rÿ èf½göKÿ ’‘}ÿ `É?ôduä`ş8ú?Èôq_õ_™õÕxWíuÿ "6“ÿ aÿ ¢Ş½Ö¼/öºÿ ‘Iÿ °€ÿ Ño]˜¯á?—æ|?ñÏògÉÕ÷_Àßù$¾ÿ ¯oı«áJû¯àoü’_×·şÎÕ]¥òıJÄnŠ?´Wü‘İéş¾"¯·h¯ù#ºÿ Òı|E\µ¿‹/ë¡½áüßè}uû&ÿ É7¼ÿ °ŒŸúW3ûb]7—á{q„ÜH}3û°?¯ç]?ì›ÿ $ŞóşÂ2èÖ_í{£Iqá­SU%-n^>b‚æŸ­uâƒûwô9¨ÿ Ûß©òÊ±V¤£µ}aeğ7Å·–p\Gñg\ËÈ¼MĞŒùx¯“«íß€^6‡Æ_ìÌ}§"Ú\¦yF¾…@üA¨ÃÆ3ROëş¥yJ-4pÚÇì¹«ø‰¢mWâî¦bÈŒŞZ¼»3×§8ÎåOğ·ìªŞñ&™«¯Š¼ö²¹ãÊşÏÛ¿k·>iÆ}q^ıXz÷Œ´ÏêÚF™u#µö©7“o+¹¸.Ã²çßë]*•8ÉI-oøœî¤äš{”QEndQE ~|xûşGŸÿ ØBãÿ F5zì³ÿ %V?úó›úWŸøûşGŸÿ ØBãÿ F5zì³ÿ %V?úó›úW‘ƒøãèÿ #ÔÅ|2õıO±ëá¿×sñsÄ%ÎvJˆ>‚5¾ä¯?jÉ¤üN–÷f Ô dnÛ”laõùAük|WÙg>í#‡øc“â7†TôşÑ·?”€×ßõùığÖe·ø…á©î®£oŸûø+ô·Ãæÿ Cÿ Äù/Ôù›öÃÿ ß×9ÿ šWÎuô'í…6uÏC“•¶•ÈíËı+Àlí%Ô/ µ“Lë(êXœùšó¤œ¦ÒîwÅÚ
ı¿>Î÷?|5,‡.úm¹'şÙ­x·í‰ÿ ^ÿ ®—É+ß|?¥ÿ bèZv9¶ñÁŸ]ªô¯ı±?ãËÃõÒãù%zß…úş§â^Ÿ¡ó5~|7³[‡ş@tø:{ 'ù×çí~ƒü?Ëà_¹:ß§ısZœ/Ã/—êi‰ø£óıê(¢»SÂ¿k¯ù´Ÿûıõò}}cû]È¤ÿ Ø@è·¯“«Æ«üY]êÓş®¬û›àü’ÿ ×ÿ Ñ¯X_µü’™ÿ ëîækwàü’ÿ ×ÿ Ñ¯X_µü’™ÿ ëîækÑÄÿ ü¿3‚‡Æ¿®‡ÆÕögì¿ÿ $Óş¾§ÿ Ğ«ã:û3ö_ÿ ’Miÿ _Sÿ èUyú~¨ßğÇ×ôfçÇoù$~$ÿ ¯uÿ Ñ‹_×İ?ÿ ä‘ø—ş½×ÿ F-|-\øâü—æÍ(æÿ CêßÙşDİgş¿ÿ öš×¼uàò+Áÿ d?ù5Ÿúÿ ÿ Úk^½ãOÚx#Ã7ÚÕöL©»bğ]‰Â¨÷$^›’Œ–Ö_’8lå6—wùœ‰~ü3ÒÌÚæ³¤XØÆ­½åyŞ(³×$úcŸJàõßÚ:]Fòü9Ñ>Õ3b&’-ª \C wlF+Á¼{ñYø‹¬5ö«pÌ€Ÿ"Õ	@¾Š?™êké¿Ù—ÀVz‚ ×š5“SÕ1›)b×?‡¥rSn«j>ì|‰¥M'-Yà?´hzõƒx«VşÕÕ®íDíƒ”€aå¯AŒü ık‚³ÿ È?ë¢ÿ :öïÚëşG­'şÁãÿ F=xŸü~Aÿ ]ù×$?‹o?Ôéû/Ğı‡ıRº+Î¿h¯ù#ºÿ Òız4_ê“è+Îh¯ù#ºÿ ÒızŸáËúêyø#â*úëöMÿ ’oyÿ a?ô¯‘kë¿Ù7şI½çı„dÿ Ğ#¬0Ÿ½?Ttbv¯èÏj¢Š+Ğ8ÂŠ( oÅ¿<5ã£kšL7ÒF6¤¹häÓzØöÎ+ÌuŠ~ÚÜé°Š÷PcûÈlœ²î“1$ãâ9àWûA|o½Õ5kÏè—-k¦[1†êx›pã†\ˆ:c¿9â¹ÙÏÀv8ñÑmAk-:/´¼2%mÀ*‘é““ô®Ôö“å¥¥úæ©­º'.¾"|Cø©ø—Z¹]Ã‘"É›
”ûH. :•ä±ì0;×„WÜ_´íøAâ—şJøv¹k%#z-Ê}ÿ ÈúëöMÿ ’oyÿ a?ô¯k¯ı“ä›ŞØFOı:öªö#ğÇÑ~G˜÷~¯ó>øÙÿ %[ÄßõößÈSşÿ ÉXğ×ı}ì­LøÙÿ %[ÄßõößÈSşÉYğ×ı}ì­^Nø´ıWæzu¿‡/F}Ù_üXÿ ’™âûMÿ ¡šûê¾ø±ÿ %3Äÿ ö›ÿ C5®+âÏô1Ãí/—ên~Îÿ òX|?õ›ÿ D½}½_şÎÿ òX|?õ›ÿ D½}½]t?„¾g=oâ¿EúŸ¾(¸k¿jÓ¾wIw3œõåÉ¯@ı™ÿ ä¯iŸõÆıÕÈüLÑ_Ãÿ <A`à.òB¹î¬Û”şDVïÀIt¿‹Z¹Ú²ÈğqË£(ıH¯7	ñAØ‚Lû¾mı±¾ï…>·_ûJ¾’¯›¿lo»áO­×şÒ®ÌWÀ½N\?ñ>ÿ Èù®¿A>ÿ ÈƒáÏû[ÿ èµ¯Ïºıøwÿ "‡?ìoÿ ¢Öá—ËõOÅŸèt4QEvÁ^EûA|#›â‘£¥¨mjÁHXó>3ÉLúƒÈú‘Ş½v¹[ÏAañËÂsÀ#{Ë#uÉ“ï8b{qè¤ç=±ŠÆ¤c4¡.»BN™t>ºµšÎâH."x'Š<r)VVA¡©tİJëG¾†òÆâKK¨[tsBÅYO¨"¾Ùø¡ğKCø—O"ÿ gë
¸MBäã uş1úÆ¾CñçÃ½gáÎ®lu{}¡²a¸˜¦_U?Ğò+Ì”%Ië÷„g‹Cé?¿á>Æ‰­m‹]7G2Œ%ÒÉÇgH£Ò·ÿ h‹¦µøC¯ÎdGÇ¡•3ŸÂ¾@ğ¡>“ãêİŠÍì8Ç|¸~ ‘ø×Ú_4i5ï…¾"´‰KIöc2¨ïåøÿ Çk²Ru0ÎOò9TU:ñKmƒëèO„?|EãÙêwÄMWC¶ß$bÂÜHR"ç™G=zµóİ})û%xÚ‡Qğ½Ä&gû]¨c÷¸Ô{ŒùÖ8uIÆ]k·¦Köoñµföš‡ÄíRúÕÈ-Ì2Hƒ‘•iÈëXŸğÆÿ õ7ÿ å3ÿ ·WÒu•âŸişĞ®µ}N_&ÒİrØå˜öUÉ<]r£I+´sª•‰’øwJm@Ót×Ÿí-go¹›nİûT.ìdã8é“_~Ñ¿òXµß¤ú!+íK´Ô,mî£WXçdU‘J°ddv<ô¯Šÿ hßù,ZïÒı•†2üÑ¿™¦ÚÛ·êwáoü”ÿ ØFıWŞÚ„ÆŞÆæUäÇ0ü5ğOÂßù)ÿ °Œú¯¾äŒM£«§ñ«¢› Òóü‘l«]ö_›?7¦Ë+»Y˜’O½zÇÀ_êş6mb+ÆşkqºYïÄÙÜ;dN˜ïµæ¾"ÒeĞµíGN™JËkq$,û,Ezìëãh|ñ
»Ee©GöI1 £mÃğ*äÃrÊI=Ÿô¼G2‹kt{$Ÿ |]"2?ÅiÑ†
²ÌA‡ı"¹Ïøcú›ÿ ò™ÿ Û«é:l²$<’2ÇÌìp¤ŸJô{´q*³Ù3øSğñ¾ø]´c¨ÿ i´<âo#ÊÆàÜnoO^õğ÷‰®ûÄš­Ä„™%»•Ø“K“_|øGÅ¶6Ò´ôÃ+Ùù²B²I]û—ÕOc_|PĞ$ğÏÄ{O‘J„»wLŒeîSù\˜­ám¬ÿ C£´»ßüÏHı’u$µñö£hçubÛ=Êºœ~Yü«ëZüıøwâÇğG4­erRŞaæªõhÏÊãşù&¾û±¾ƒS²‚îÖUÚtG"†R2®ºæ¤¼¿áÎzÑå¨üÉè¢™É6ñ«”m­´çiô>†·1<ö¾ÿ ‘CCÿ ¯ãÿ ¢Ú¾U¯ªÿ kïùô?úÿ ?ú-«åJñªÿ _×Dz´ÿ ‡ë«>ëøÿ $—Ã_õíÿ ³µwUÂüÿ ’Ká¯úöÿ ÙÚ»ªöç’¶>Iı­ÿ ä¢ißöOı-xÖ•ÿ !K?úìŸú¯eı®?ä¢ißöOı-xÖ•ÿ !K?úìŸú¯Ÿñíçù¤¿‡òıÑºüğñuÃİø³ZSºI/ff>åÍ~‡×Á¿ü?'†¾%ëö…QîZâ/tïÏ…tb¾(¿_ĞË´—§êo~Í1‰>/iYş§aÿ ~š¾Ô¯Šf¹„?´€9Ğ}|¦ÿ 
ûZº¨Õ¿ˆÏÿ jŸù*_öãójò}6f·Ô­%C‡IQ”û†½Cö ›Ìø±t¹'Ëµ…~Ÿ.­pÿ ô<OãM‰Kî“~; 9cø(&¼êzÖÓù¿S¾m*w}¿Cô
6İ“ÔŒ×Ç¿µOü•/ûq‡ùµ}‰_şÕ_òT‡ıxCüÚº±}Fsaw—§ê-Ğ?ä=¦ÿ ×Ìú¯ÑjüéĞ?ä=¦ÿ ×Ìú¯ÑjÓü?ŸùˆøşGÆ´çü•Ëÿ ú÷ƒÿ EŠàüÿ #–ƒÿ _ğèÅ®óöœÿ ’¹ÿ ^ğè±\‚?äsĞëşıµÇGøëü_©×Sø_/Ğı¢Š+Ö<À¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š(Í?iÏù6ß‹ö)jßúG->µışÓŸòm¿?ìRÕ¿ôZş}+ƒº:ilÏ¸àÿ òr~%ÿ ±Jçÿ K,«õŞ¿"?àÿ òr~%ÿ ±Jçÿ K,«õŞ¶Ãüu> ¢Š+¤È(¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š +ø£ãM'Â~ÕŞúún$µ’8mË2Ge!B¯SÉJéµm*×\Ó.´ûè¼û;¨Ú)cÜWr‘‚2#ğ5ÀÃ8ü;ÿ ¡{ÿ 'nøåcQNQqSZn1jLø†»ß¾&²ğŸÄí"ÿ Q™mìÿ y“?DŞŒ Ÿl‘Í}Cÿ ãğïş…ïü¹ÿ ã”Ã8ü;ÿ ¡{ÿ 'nøåsS¡8IKCyÕ„ââîzõ¶£n³Ú\Eu¬¸u?B*zátş	ğ¾¯mªiš/Ù¯­Ø´RıªwÚH#£9‚zŠî«¹^Úœ}Š*¦­¥ZëšeÖŸ}ŸguE,{ŠîR0FA~‡{h
×Ôùkö´ñ%†¯â-ÂÊê+©l¡Ïä°`ŒÌ0¤ø^õàõö÷ü3Ã¿ú¿òvçÿ Qÿ ãğïş…ïü¹ÿ ã•ç<<ÛnèíUà’FçÂ¿éş"ğ‡%ÔS4VPÅ,hà´n¨«Ç#½r´7Â©~!xv+İ6/3ZÓ·4q³F~ò}xÈüGzé|3ğgÁŞÖ"Õ4}ìwñT—íS>Èè}+µ®ê‘ö±´HKÙ»ÄüÚš-æx¥FTb¬0T Æ½áÇM{áœ-gÇ¨éLÛ¾ÇpHØORŒ>îí_VøÓàï„üy+\jšZ‹Öûe³¥?R8oø5ÄGû$ø5eou—\çc\EƒíÄyıkŠ4jÓ~ë:åVœ×¼*ûö˜ñO-ÂºéúßîÄ‹1¸g©_•BãÔƒù×¸|,øÃ¯ÇbdûN¡;ïn&YO^O8çŞ®ø?áß‡¼nÑhšdVŒãMËÊşÅÛ'ÙÇµt•Ù¸ë'vÎIK›D¬“?k:î/é—Î¬læ²Dÿ Ã¹]‹/×æñ¯¶¸–Îâ)ásÑ8tuê¬Aüëô3Å^Ò|m¤¾›¬Ù¥å«Á[!‘»2°äq^amû'ø*±3Í«\GœıK”}²¨õ®7‡š›q}nuªÑåJG)áŸÚÆ>·‹DĞ|9o&¸ÈMKÌ&‡C#&Ü/¯,F{•àş<ĞõøÃVÓõIMÅüS±–|æ“óoüsŸÆ¾õğï…ôŸéëc£ØC§ÛJB¸,}XõcîI5Î|@ø9á¯‰$ú­´‘^¢ì[ËGÙ.ßCAPqZU£)ÙŞìÎXÆêÖGÈ_|smğïÇ6zÍå¼·6±¤‘Èc~HÈ€NqÔŠûOÀ~+›ÆºjÒisi0NäÛÇpÀÉ$_ÂäòçÓŸ­q^ı™¼áËä»hnõi#!‘uUÑHïµUAúŠõePª €oF2„m&eRJR¼Obø‹uğ£ã‡‰o"C=Œº„ñİÚ7Æd'#ı¡œş½}Yá/è^8ÓÒïGÔaºR¹h·,~ÌEbÇğ[Âm«k:•î™«uªÊe•¯•dg¨´{õ÷®gVı–|©M¾õ,rÒç+ÿ ‘ë*Q©J
Øº’…I9lx¿íE¬èúÇÄ(N—473Aj±]Í	L˜‘Á åÚ½;öYøusáıïÄZ„-Æ¤«´n0Âs¸ö1ìï]W…g?xVé.–ÊmRâ3¹RHúí )üA¯NéÀª£KÙ·9nÿ QT©Î”VÈç~ ø2ÛÇŞ¿ÑnNÏ=3¸Ï— å[ğ?¦kàÿ x_Rğ~µq¥ê¶Ímw	ÁVèÃ³)î§±¯Ñ:Âñg´/Y­¶·¦C~‹÷ÁWO÷\Ãğ5¨s¾hîU*Ü‹•ì|3à?ˆzÏÃcûCHT°Û5¼£tS/£ş£W±ÍûajSXˆ­ü3n—ì0%k¦x÷×0 ş«¶¸ı’ü4ÅÒïX·RsåÇqQíóFOë]o‚şøCÀ·	uc¦ı¢ù9[»ÖódSê¿Â§Ü iB[r¹Y)Óoš×fÁêëuuãO;Mâ=I6Ç‹ƒk÷qü$ñÇ`=I¯[¢ŠìIE$¶9[mİ…QLG;ñÁvŸ<){¢İü‚eÌRã&)*ÃèBE|!âÏ	ê~Öî4­VÙ­î¡?ğ^Ì§ºŸZı¬Ox+CñµµÖôØu—îŸî¸Ã/àk’µwÍÎšU¹+ØøGÁ^8ÕşëI©è÷LÀmxÜnUşë/qúúWµ/í‰}öŸğŒ[›ÜcÍmåç×fÌş«µ¹ı“<<Ìéu¬[©9ÇqQíóFOë]_ƒ~x?Á7ÜÙéŸj¾Œån¯[Íu> •O¸ ÔÓ§U.[Ù9Óoš×g/ğ—Â"ñV¼¾;ñ«w°®™§•*¶Èİ_gğäp3Ï9=«ÔüIâè÷7úÔV°E6é\.âA¤ú
Ô®WÆü1ãë«{{LûtÖèc¾Ñ,{Tœ‘„a}kyEòr@ÁI9sLø
Fİ#7©Í{ì·¯Øè&KëˆíVîÍà‰å`ª_r°\ä^ûÿ ãğïş…ïü¹ÿ ã”Ã8ü;ÿ ¡{ÿ 'nøåsÒ£:rRV:*UDÓ:Ïx’ÃÃÔ/¯®¢¶`}†G »m8Uõ$öùêNI5÷¦¹ğsÁş$¶Ómõ#íi°[UûLËåÆ:/Êã=:œšÉÿ †qøwÿ B÷şNÜÿ ñÊuhÎ¤ù©ÕŒ#fxgì§âÇö÷·1Zı²ÏË‰¦p¡œ: ägò¯¯3G"¼ÛşÇáßıßù;sÿ Ç+¾Ñô›MKµÓ¬"ò,ícEâÛTt$“øšé¦œ`£.‡=F¥'%ÔùGöøIsá½~ãÄšu»I£_?™?–¼[ÌO9ôV<ƒêHôÏ‰Ç#Ã"É2:ÊÊpA5úE41ÜÂğÍËŠUãuX ƒÔW–ø‡öfğ>½pÓGis¤»·ö|ÛTŸ÷X2  W°Î.ğØëtÕ¦y„ÿ k-sEÓãµÖ4¸uÆv­À˜Á+{¹ÚÀŸ|
×ÒuïşÓÌv3À4oZÈ²^-¹'ÍÁÈBçï1ô Ô•İi?²¿‚4Û,ßÚ:š|«»ÿ !ªŸÖ½_LÒìô[¬ì-a³´ˆa!"ı ®ˆÓœµªïåşfœWğÕ¿®„¶¶±Y[Eok 8Ô`*€á_9~Ó_î.îÅÚ5¹›äP‚1–ßĞZúJŠº´ÕUg¹æé½Ízõß‡?´–½à]6-6êÚ=oN„m‰fÇ,kıĞø<{q_Ax³ö{ğW‹.æ]9´ë§9y´çò·}Wsï·5ƒcû'ø*Îád–mVõæîP)úìE?­rÂ•X?u2©NKŞGß<kñğxÂzRhxÄ÷k)•áN…Œ›@Aôôæ½ûÀ¾´ğ†,ôk2]!\É3}éd<³Ÿ©ı0*ß‡|/¤øOO[-OƒO¶”…q¸ú±êÇÜ’kFhRâ"nE*Ëœdƒ]‘*ovÿ ¯¸å”¹š[#æoÚëÄvÒh:]µÜ7Ví4³Ç1d(PØèNÕó}½ÿ ãğïş…ïü¹ÿ ã”Ã8ü;ÿ ¡{ÿ 'nøåqK9IÉµ©×ĞŒT{¾øŠÇ^øk ­µÔ2ÏkiBÆÈ6á‡QÓõ¯9ı«¼cs¡Ãâ˜Ş;}Bİ–ŞUb¸Bp1êÊLúW¦øgàÏƒ¼¬Eªhú?Øïâ©/Ú¦|<3‘ĞúUü)ğ·/¡¼×t¿·\Ã”ö‰cÂäœa¤×UhÊ´|ÎzRTå~‡À•õWìëñ²ßUÒâğÖ½x°êË¶ÒâáÀÇÙ2ˆvõJí?áœ~ÿ Ğ½ÿ “·?ürøg‡ô/äíÏÿ ¬©S©Iô±¥J¨I¥ªÚnŸo¤éöÖV‘ùV¶ñ¬QG¸ª£ dòxêÍv¡\GÆÜøÃá®³§Y)’ñ£YbŒuvF´{úšíè¬êAT‹‹.p’’è~l:4nÈêU”à«}+Û>
şĞ_¼/w£êZuÍØ´öÏjW’Àe_$`duëÒ½×Æß ¼#ã«ç¾ºµšÂúC™.,$´‡ÕI÷ÆMQğ¿ìÕà¯ß%Ù·ºÕ¦ŒîOíC¢Ÿ]ªªârS£VÙØéJsZ£Gàî¡â?Xj$×Ùí¢Õ$W±Ósò[À €FyËg9ï€{×š~×^"°¸Ó4]"¨¦¾á§–ÜBàn¦Iïé_Fíq:WœIû:ü=šFwğşY‰b~ÛqÔÿ ÛJÚµ9N*ÛüŒ©MF\òÜø~¾Üıuû[á~mos—Vq´SÀoŒ‡=G\AüißğÎ?ÿ è^ÿ ÉÛŸş9Zø3àïêñêšFöKèÕ•eûTÏ€ÀƒÃ9¥*4çN÷ê:µ#Q+;û]ø’Æít-"Şê)îà’Y§Š6bÈP¡±Ğxö¯›ëî	?g_‡³HÎşË1,OÛn:ŸûiMÿ †qøwÿ B÷şNÜÿ ñÊç–¤¤äÚÔÚ5¡¢¯¡cà?ˆlu¯†ZV×0Ëqkn p^6RGÌ:Œã?müGğE¿Ä/_hÓ‘Ê»à˜õrU¾™àûYşø-àßëêšNöKøwåûTÏŒ‚ät'µvõß%í"Ôúœq|’¼z^$ğŞ£á-bãKÕmÖòÃ#ìÀ÷±«àˆú×Ã}Xßi( ³[LE2†{v ‚+îx@ñÍªÁ®iß*ıÇl¬‰şë©?^g?ì—àÉ¦.—šÄ*Oú¸î#*?8Éık…P©x³³ÛBJÒGqû`jw"OÛÇ¨0Ú²µËH›¤aAü7W¤üğ­i5ßŒ<Y#OâmQ@Ù Á¶‡²cøIÀàt ZÜğ_Á?øt¹Óôß:ù>íåãy²/ºç…>ê®îºáŸ4İÙË9&¹b¬‚Š(­ŒÂ©kÕ†c-æ£w•´jY¤™Â>½Oµ]®+Åüã-^MSXÑşÙ}"ª´¿i™2`p®Oj‰s[İÜ¨òßŞ>Ö.’ûV¾¹îM;È¹ô,Hşué¿³7ˆ¬|;ñ)MıÄv±]ÚIl²ÌÁT9*ÀxÛÆ¾†ÿ †qøwÿ B÷şNÜÿ ñÊ?áœ~ÿ Ğ½ÿ “·?ür¸èĞ&­¡ÕR´*&™é*ÁÔ2Êyt5óŸísâk	ôF‚î)¯Vå§šÜ3F¡JÀtÉnş•í~ø x/LºÓôk±Ùİ1i£ó¤}Ä®ÓË1#C\¿ü3Ã¿ú¿òvçÿ Võ¡*‘å^FT¥o™Ÿ×Û³¯‰,5†Eœ715íš43[‡ÓÄ½pA57ü3Ã¿ú¿òvçÿ U½à?´R×Q°Ğü‹ËY°Éö¹Ûkp\ƒøŠš4çM»ìÇV¤j-7G?ûNxšÃOøk}¥ÉuÛï$ßxó+–Û×^¾õñµ}Ïª|ğ³©]_Şh^uİÔ­4Ò}²uÜìrN€Oj­ÿ ãğïş…ïü¹ÿ ã•„¨NrrvÔÖ5¡ò£Šı’¼E§ÿ Â'¨hïuz‚Ş™–pÑ‘@*;ò§§·­{<)kãê-ç]G´HLl9Và€k’OÙÓáän¬¾Ã)È?m¸ÿ ã•é
¨ â»yy¡É?C—›–\Ğõ?<¼aàıOÀúäúV«nĞ\F~VÇÉ"öu=Áÿ <ÔŞñæ±ğ÷Z]KG¸K²E İ«ı×Çê;û»ÅÑ<ib-5½6Bå|ÀC!õVeü¯0ºı“|q3<wZ½ª“‘W•ß4dşµÃì*AŞ,ëöĞ’´‘Ä7í‰¨5‡–·ØÀ˜İ1?îmÏáº»ïƒÖõ-RxÁ™õËÈü»KyoÙ¡>‹ü$ƒ°'<“]ƒ~x;Á71İZiÆòú3”º¾5”ú€ ûšô*ì„|Ów‘Í)+rÁYQ[UÔµ[=ÎK«û¨lí£igpŠÔÕªã|UğƒÂ>6Õ?´µ­'í·¾ZÅæı¦hşQœ#Üö¨—5½Ò£ËxøwÅZ„Z·Š5{è`¹»–hÉşë9#ô5Ü~Î¾ ²ğïÅµˆímæŠK:VÚªÌ>\“Ó$ø×Ò_ğÎ?ÿ è^ÿ ÉÛŸş9Gü3Ã¿ú¿òvçÿ W*¤ÓVĞë©Z5O©éÈ²¢º0taÊrp?¾§Äß	µ´;SUµ&k9ànÇ(O£ÔÚ·¼ğÿ @ğ½Ì:‡Ø"¸`ò¯$›ˆçc½«¢®ÉÁT,XÉÂ\ÈüèšÖÿ ÂºğîŞKKû•šW¬¤_{xCÇZOŒ¼;±cw	£2— ÀØå_Ğ­Gã†şñäjºŞ—ÜŠ0“Œ¤ª=©ÙÅyÄß²Oƒe²Şë0¯÷â,Î2ZÂœjRN;›NP¨ÓØñOÚ+Ç~7ø€Í§L·60‹T™VF³2ã'ïŒ×sû6ü¹“PƒÅºÕ»AoÍao*á¤b8”Ê;z{sëıŸüá)ã¸‡L7÷Qœ¬úƒù¤P¼.}öæ½>ÍóËV*•y—,v
ù‡ö¼×ìo®´6Şæ)î­¼é'HÜ1vÀ¡±ĞœúfêÚ;ËY­æ]ğÌ†7\‘•#d{WÃ8ü;ÿ ¡{ÿ 'nøå:ÔåQr­…JQƒægÄ5÷¿Â]{OÖşè?a»†á¡²†)QÕ eaÔEcÃ8ü;ÿ ¡{ÿ 'nøå^Ğşø#Ãzµ¶§¦èŸf¾¶mñKö¹ÛiÆ3†r^âŠ4åNéìÂ¬ãRÍnòŠ(®“çÚëÄVhzF‘ÔRß‹£<£‚Ñ¨B2Ã¶KqŸC_/×Ü3~Îÿ ®&’Y<?ºI³7ÛnI9?òÒ™ÿ ãğïş…ïü¹ÿ ã•çKRRrvÔîhF*+¡[öqñ†©ğ¿I²‚æ6¼²E=¾ñ½?xÄ½pAÖ/íY¯Ø[ü?])®£şĞ¸ºÒÜ0/µrKØ{û×g |ğO…õ{mSLÑ~Í}nÅ¢—íS¾ÒAÈ<ÔTZ·ÀoëšÖ¡}¡ù÷—R4³Iö¹×s’pø
é«Ô‡.‡=9F¹ğµ}wû*øŠÂçáùÒEÌCP¶¹‘šÜ°±°CÜuäzWGÿ ãğïş…ïü¹ÿ ã•cMøà='P¶¾´Ğ¼««iVh¤ûeÃmu ƒƒ&#½Es¥&ŞÏBêÔHÙt3ÿ hÏiúOÃZÆâæ5½¾DŠ}ã{üêI® šø¦¾ïñÁx«X¸Õ5MíW÷eûTé¸€ á\ÀgÃ8ü;ÿ ¡{ÿ 'nøåeR„êK™Ø¸U„#Ê=ı‘|Aa“¬éİEó\¬ñÂîH¥0vƒ×sQ]‡íE§]_ü,™í‘m®¢š`¿óÌdô…jEû:ü=†T‘<?µÑƒ)ûmÇÛJôY¡æ!š5–)«Æàe#zŠè”%R—#ßOÀÅMF§2ØüÚ¯føUûCëĞáğòèÃ]
Ål•e)"–9Ù€­¸dğ8<õ¯`Õ¿ej—ÏsÔ´åc“oi:ùcè¦k³ğOÂ|?ùô1èŒ5äÄÉ1õùO¢àV4¨Ô‹ŞÆµ*ÂKkŸ,|uÒü`÷ºV½âÕ)õ˜Ek ùmNDgßæÏS×¯§—Å!ŠDqÉRü+ô7Å~Ò|m¤¾›¬Ù­å«ÀC#Œ¬9é^aoû'ø*µ™çÕ§ŒıK„}²¨CÃÍNñØµZ<¾öçAğ¿ã
|NºhôıòÚÊÚn/nYB,¼b4;»œäp:s[Ÿ¼7qâï‡ºæ•h3u<‰¼ÊCüJãñ­ÍAÓü3¦C§év‘XÙD0Ä0>§ÔûMhWlãÏÔã‹ä’k¡ùµ42[ÌñJ¨ÅY`©G­{'Àÿ _t;ı/S°¹»†IMÄk·!Ê€UƒÇÊ9ïÅ{ÿ ~øOÇ×{{i-ûıû«»ÿ ¼*O¹÷¬Ÿ~Ì	ğıò]<7z³¡Ü©¨J®€ûªªƒô9ÉNZrÑS©NkT_ø;¬x‹ÆK©x§Yße§êSNÓ3•Š%Ïï9–Ï^3ŒôÅzM5UcPªªŒ 
uw%d’8Û»¸R‡ih£} üêñFŸu¤ø“S³½Fê™En¹Üyüzş5³ğ×â6¡ğËÄCT±.ĞÅ=¼„…•	ìr¾ÀøğWÃîëR·–Şü ¿l³qŒ£ l‚ñ¬_şÌşğíÒÜ=­Î¯"œ¨Ô¥ƒşªªkÎ§‡©M«=éÖ„ÓºÜó¯x·Ç¼ªÜÙi øfÖ<¡œÉ%éNv+màÇSÒ¾r¯ÒXá(V$EH•v„PÓ zW“xƒö_ğV½¨Iv‹}¥4‡sEa2¬y=HWFÇĞ`UTÃÉËš$Ó¬’³<§àÆH¼¡ÍáÈô[Í_U»¼ó-#µ*Ë*‚WsáŠú®[¸í-÷rGjŠ»¤iO_˜ãzäüğ‹Ã?wI¤ÙxËµ¯.[Ì˜LôQşè­x+FñÖ›†¹göëHå,~kÇ‡ €rŒF?u¥8ÂÛ³™ò¹_d|;ñSW¶×¾"ø†şÎEšÖk·1È§!Ôp{fğŸX¶Ğ~$xzşòE†Ö+¥ó$c€€åwè3šúËşÇáßıßù;sÿ Ç(ÿ †qøwÿ B÷şNÜÿ ñÊä§Btå&´:§ZN=Ï@—T²ƒO7Ò]À–Aw›–öúîÎ1_ xÿ Vƒ]ñÆ½¨Z¶ûk›Ù¥‰¿¼¥Îâ+íY>x>_EáæÒ3£Çqö¤¶ûLÜK‚7nß»¡<g“ÿ ãğïş…ïü¹ÿ ã•u¨Î¬“íÿ  ŠU#M;Ÿ,üÖ­<=ñK@½¾•`µYYYj¦ôdÃ,+îØfâ%’'Yca•t ‚=A¯8ÿ †qøwÿ B÷şNÜÿ ñÊì¼+á'Á:PÓt[O±Yi<¯1äùS–$öõ­éFP,ŒªJ3—2<Cö øSq«*x³J¦šü»ècbƒîÈ|t>Ø=|Ëc{6›}ownæ;ˆ$YcqÕYNAüÅ~W™x¯ötğWŠîé¬¦Òîd;M6Acë´‚£ğ¹§‡’—53xV\¼³24?Ú“Â>ŠëS{=QcmŠ[»–p9ØÀmÁ=2ExgÆïëŞ6}+Äš•³iúEç›—fçæ&Òd?ïûzŸ |7û3x+Ã·Étğ]jÒ!Ü‹¨Ê®€úíUP~‡"º?ˆ	ô?‰M¦k›¥OßåEm EmÛr‚„t"®¥:•¾äSœ)ËM‚ëì‚Ÿ¼=â-Ağá’[]r;e¶ïÆœ•pÆqé]•ğÀZ;‹Ã–ó8ïtï8?ƒ±¥v^¥èjWNÓm,û¶°,cÿ ®9S½ŞŒš³MRÔ¿EWI€WÌµ]õÖ‡ãÏê¶R´PÛŠEdÈúõé_O×)âŸ†z5í+UÖ kÆÓ•Ö;Y6˜v9u#æÆ:g Öa)Û—K3ZrQ¿7cŸø[ñËCø…§Á÷0éÚàP²ÙÌáw·¬dıà}:Ôó?µnµ£¯cÓgšÕŞæ9-áîb;=ó]'ˆ?g/ëìÒ-´É˜äÉ§Êc‚ ü²4ÏÙOÁ:}ÒË3êzŠ)Ï“ur¡×b)ıj*F¥Hò4‹ƒ…7Ì™âŸ³ŸÃ»ŸxÚ×T–N• IX|¯ åzœàŸaî+ìÖU‘YXV ªºN‘e éñXéÖ°ÙZD0ÂT~®Vñ‚„yŒ¤å.f|Cñ»á=ßÃŸM4»è7rµFU3É‰½íê?yî›©]h÷Ğ^ÙNö×p8’)£8eaĞŠıÔtÛM^Î[;ëh®íed†tŒ=Á¯)Ö?e¯ê—XRÿ K“ÀÛùH­Â¸Q~ã;xÉ{ç›èµö­gb±jš¾¥r«´Cp`İîËµ†~˜JÛğD~#ı¡<Ik®xŠ³ğ–™/›oc‘òƒ_Û§aŒšìôÙ‹ÁÊÏ-½Ş¬Êr«¨N÷TUèr+Õmíâ³·"HaB¤q¨UU  ®¨BWR¨ïcR­`šx­ayf‘b‰Yİ‚ªROJøOãv½iâoŠíıŒÉqjÒ$i,g*ûS ÷SÍ}©âïhŞ:Óc°×,şİi¢eÍxğàQèÇó®GşÇáßıßù;sÿ Ç+:Ô§Q«ZÈºU#M>ìø÷Àš¤/4+û–Ûom{²7¢‡ŸÊ¿Aloíµ+d¸´¸ŠêİÆVX\:°ö#ŠóÏøg‡ô/äíÏÿ ®›Á¿|=ğın†ƒ§ı€]3~şI7mÎ>ûu=+Z1•8ò³:²Œß2<#öøKsö÷ñ~“nÓA"¨GÉ€À—„`Lg½|å_¤ì¡”‚2ó?~Î	ñEÃÜ>M.áÎ^M:O(øGà+šxg{Ãc¢Õ­#Ã<ûQkŞÓaÓõ+(µëhWdRI)Š`£ /†qŸzß_ˆ.ı¢µáİ2Ñtˆ7óDÆB#îLç(=øÍwzì£à«;…’iu[ôa¸¹P§ë±şµêº‡tÏéÉa¤ØÃahœˆ¡\}O©÷<ÖÑ§9İ[™JpğÖ£´=ÓÃš=™c•ikÅûÔû§Ü×~Ò_gñuœ~"Ñ 3j–‘ì¸·Œ|ÓÄ9ì¼ñÜ} ¯t¢µ«MUVftæé»£ó]”«Fàƒ^ğËãö½ğŞÔiâ8õ]$­­Ãhó×c€úE}9ã/ŞñÅÄ—WÚgÙï¤å®ì›Ê‘©åcîA5ÉÛşÉ~†ew»Ö.˜ä¸Œ)öùcõ®8Ñ«î³¦UiÍ{Èã'ı¤¼Yñ
â=Â^M?Pºù<á1¸tØªª©íÿ |Ÿü6¶rNo5+‡7×ŒI3LİNO$ƒ?^õ¥á_h^	³6Ú&™„m÷Ù.ÿ ï9Ë7âkr»a][»g,¥Í¢VG’şÓ^»ñ7Ã†–Ê6šm:ánÚ5-VVÀöŸ 5ñ~“×–x§ökğ_Š/óì÷ZTÒÒgJ¨Œ}v²°€ÉW)MÊ=Nšu”cË.‡•ü7ı¥­<à}çH¸¼Ô¬Õ£¶0²ˆ¤’7œåHÎ88¯wø]¿…>Ûâ‹–mNúW¼x_…µFÆØ€ì ÇlšÉğWÀxú;ëkI¯ï£9ãPHc>ª ßÙø›Ã:gŒ4y´½^Ûív•/˜É’#• õµÓ4¯-YÏ'ì´GÈŸ´×‰,¼Gñ'6QİÃgi³I]Á™ˆuÆïÎ¼®Öo³İC.3åº¾>‡5öÏü3Ã¿ú¿òvçÿ Qÿ ãğïş…ïü¹ÿ ã•È°õ¹®·¹Õí Õ¬w:¿câ].Şÿ NºŠîŞd'GC„zW–şÑeñö“­¤ÃækvWÊ^·u*?Ú$}Hî+µğÂ
øP–÷BÒşÃu$~SÉö‰dÊ’0îGP+®®ÊUci“¦îÏ¯ø¼ã­'T¸‰ÔXÜƒ4da‚ı×õÁ5÷|>,Ñ®41¬¦§jt¢›ş×æî{nµƒã/ƒşñä­>©¥!¼aµÛ±ŠSõ+÷¿à@×d™7}¿Zû¿h‹ú+5”#Vœy¦³p©.g¡óÅ/'¾!jú­¨w‚â`–àJ(¼{œ{×Ğÿ ³Á¹ünş"Ö 0ê×1ì·¶|Öñ¥½½;©½ğoÁŸød¸Ót¤kÕéwtÆY¸'…?î]½hªZ½XªÕöš-ÖeK1
 d’p|UûIkÖ> ø¡s%…ÄwPÁoXX2–$:ãv?
ûÄÓüU£Üiz¥¿Ú¬. Å½“pÊG t5ÃÃ8ü;ÿ ¡{ÿ 'nøå*ÔçVÉl‡J¤iİ½ÙñV—r¶z•¥ÃŒ¬S$‡¯Ñ-7Z±Ö4Øõ;¸n,¤MâxÜÆ3Éí\ü3Ã¿ú¿òvçÿ V¶ğwÂN‡©höºG•¦êEÔ?i˜ù›NWæ/‘ƒèEU(Nœ\tIFrLùCö€×ì¼GñKU¹ÓçK«dXáFr¬U lã9ö®3Ã7Ñé~#Ò¯&â+{¨¥÷UÁ? ¯³?áœ~ÿ Ğ½ÿ “·?ürøg‡ô/äíÏÿ ¬!B¤$§u{ßõ6•hJ<ºö=OÔ­ukHî¬®b»¶e%…Ã)QVkœğoÃßø+¨ô?ì	tÊÒ:I7Î>ûu=+£¯@â
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€<Óöœÿ “lø±ÿ b–­ÿ ¤r×óå_Ğoí9ÿ &ÙñcşÅ-[ÿ Hå¯çÊ¸1;£¦–Ìû‡ş	ÿ ''â_û®ô²Ê¿^6×äGüşNSÄ¿ö)\ÿ ée•~¼VØ€Î§Ä&Ú6ÒÑ]&Bm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´Pm£m- ›hÛKE &Ú6ÒÑ@	¶´´P™şÓ‹ÿ ÙñcşÅ-[ÿ Hå¯çÆ¿ ÿ ÚsşM¯âÏıŠZ·ş‘Ë_Ï…pbwGM-™÷üşNSÄ¿ö)\ÿ ée•~¼Wä?üŸÚSÄ¿ö)\ÿ ée•~½m­¨|u>!´S¶Ñ¶º†ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m ÑNÛFÚ mí´m 2ı§?äÚş,ÿ Ø¥«éµüøWô!ûN¯üc_ÅŸûµoı#–¿úàÄnš[3î?ø$üœ§‰ìR¹ÿ ÒË*ız¯È_ø$üœ§‰ìR¹ÿ ÒË*ız­¨|u> ¢Š+s0¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š( ¢Š(Ìÿ iïù6¿‹?ö)jßúG-=ÕışÓßòmìRÕ¿ôZş{«‹º:)lÏ¸ÿ à_òr%ÿ ±Fçÿ K,«õóm~BÁ äå|Kÿ bÏş–YWëåkCà"§Ä&Ú6ÒÑ]Z‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´P‰¶´´PeûO/üc_ÅŸûuoı#–¿Êş…?iïù6¯‹?ö(êßúG-=uÅˆİ¶gÜŸğHù9_ÿ Ø£sÿ ¥–Uúù_ŸğHù9oÿ Ø£sÿ ¥–UúûZĞø©ñ¢EtÜmê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\óÚ{şM«âÏıŠ:·ş‘Ë_Ï]Bß´÷ü›OÅ¯ûuoı#–¿šâÄnŠ[3î_ø#ÿ üœ·‰ìQ¹ÿ ÒË*ı}¯È/ø#ÿ üœ·‰ìQ¹ÿ ÒË*ı}­(üTø‚Š(­ŒÂŠ( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( Š( 3ı§ÿ äÚ~-Ø£«éµüô×ô-ûOÿ É´üZÿ ±GVÿ Ò9kùé®:û£z{sÁ¿äå¼Kÿ bÏş–YWëı~@Á¿äå¼Kÿ bÏş–YWëımGà"§ÄQElfQE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE QE y—íAÿ &ÓñkşÅ_ÿ Hå¯ç¿¡ÚƒşM§â×ıŠ:¿ş‘Ë_Ï=q×İÓØûŸşûÿ '/â_ûnô¶Ê¿`+ñÿ şùÿ '/â_ûnô¶Ê¿`kj?>!(¥¢¶3ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€ŠZ((¥¢€<Çö ÿ “høµÿ b¯ÿ ¤R×óË_ĞßíAÿ &ÑñoşÅ_ÿ H¥¯ç’¸ëîéì}Ïÿ |ÿ “—ñ/ıŠ7?ú[e_°5øıÿ {çö—ñ7ıŠ7?ú[e_°[kZ?	>!(¥ÛFÚÜÌJ)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ J)vÑ¶€Š]´m ¢—mh(¥ÛFÚ óÚƒşM£âßıŠ:¿ş‘K_Ï%CŸµÿ Æ4|[ÿ ±GWÿ Ò)kùã®:û£z{tÁÿ äæ<Mÿ b×ş–ÙWì~>ÿ Áÿ äæ<Mÿ b×ş–ÙWìkGá"§ÄQEoc0¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQE ¢Š(°QE€(¢Š,EQ`
(¢‹ QEXŠ(¢ÀQEÌ¿jù6‹ö(êÿ úE-<UışÔ?òlÿ ÿ ìQÕÿ ôŠZşx«¾èŞÇİğGŸù9Ø¡uÿ ¥¶Uû	¶¿?à?òs&ÿ ±Bëÿ Kl«ö´£ğ“=ÄÛFÚZ+bÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih 0ı¨—ş1ŸâßıŠ¿ş‘K_ÏCÿ µü›?Å¿û5ı"–¿
å­º6§±÷OüãşNcÄßö(]ém•~Ãm¯Ç¯ø#¿üœÏ‰¿ìPºÿ ÒÛ*ı†­)|$Ïq6Ñ¶–ŠØÌM´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š M´m¥¢€mih ÛFÚZ(6Ñ¶–Š óÚ‰ãş.Ø¡«ÿ éµüïWôEûQÉ³|\ÿ ±CWÿ Ò)kùİ®JÛ£j{uÁÿ äæ|Mÿ b…×ş–ÙWì5~=Áäæ¼Mÿ b…×ş–ÙWì=kKá&{¢Elf6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E 6Šu Ú)ÔPh§Q@¢E yíEÿ &ÍñsşÅ_ÿ H¥¯çv¿¢OÚşM—âçıŠ¿ş‘K_ÎİrVİSØû¯şëÿ '5âoû.¿ô¶Ê¿aëñãşëÿ '5âoû.¿ô¶Ê¿aëZ_	3Ü(¢ŠØÌ(¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š óÚşM—âçıŠ¿ş‘K_ÎİDŸµü›/ÅÏû5ı"–¿ºä­º6§±÷_üü}—ö£ñ4-÷¿á»ñ–gúWìFÓ_?ğKø—şÙŞ$¶éÿ mF}.`?û-~ÄUÒ~é3ÜnÓFÓN¢¶¹vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.vš6šu\í4m4ê(¸ÚhÚiÔQp´Ñ´Ó¨¢à7i£i§QEÀnÓFÓN¢‹€İ¦¦E»MM:Š.–~Õy?³Å¶=?áÕ—ó³”Zş{´í]F‘:Ûúıkúı±&û?ì¯ñaºgÃwÉÿ }BÃú×â_Áïlxfæm»¶Ş2ãˆ­rÖİSØú?ö_øDà¤Úş‰÷3y¯iÛë›Jøÿ È_¥~Ãí5øùà1ÿ ?üRWıÈ›ÆZ° ñÿ QÜmüüáù×ì-]/„™î7i£i§Q[ İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E İ¦¦E¡ó÷íóªcşÇ¿î3·~š¶ùÿ ®³GşÏ_Á?¾ÿ ÂÊø7¬ê~G›åkó[nÆz[Û7şÏ_cÿ ÁQµ¦şÆ,·-´êz}°÷ÅÔrãòˆşUGşİ ­ì•}s$J´<Oyr¬Ëœ¼_Î#\µw5†ÇÇ¿¶2ÿ Â¯ÿ ‚«iúğıÍ´Úæƒªú›mÒ_ÌÇ'ç_°5ù9ÿ ¤ğÌÚíào[æ!¨h
Èü¶¶¸‘‰úíš/Ê¿S<âH|gàİÄøò5k{øöôÛ,jãƒURê)£^ŠuĞga´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢Ea´S¨ ,6Šu†ÑN¢€°Ú)ÔPE:ŠÃh§Q@Xmê(¢EcàOø,wˆ¾ÃğÂZ2¶×Ô<D³0şòEo6Gıõ"Â¾‰ÿ ‚gxkşØ›á¼L›f¼†êşCıï6îgSÿ |…|?ÿ £ñh¸ñOÃ/#a­,¯5)W=|ç4'éäIùšı<ı|"ŞøğãÃ’&ÉôÏXZÌ1Ş-º?‹n?qÔøã±ñüÛÁÚ?şø©cÜúN·.Ì*—0'é›UR+Û¿àŸŞ0ÿ „Ûö?økvdß-ƒi3ÊıšWAÿ €F§èEmÿ ÁH>ÂÇı~"[G™w¥Z&µ%>Ìë,‡şı,£ñ¯˜¿à~?şÖøCã_Ë.ù´]]/¢V<¬71cØ=»Ÿ«Ó¦í!KcôŠ}Õs”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸE”Sè¢à2Š}\QO¢‹€Ê)ôQpE>Š.(§ÑEÀeú(¸¢ŸM’Eİ‚"Œ–c€©¢à~7şÜkÿ Ãş
C£x)Ò KÍÃ£lŒ’H>€Ü¾~†¿o~%~ÂöíûFÁJïüm*›‹+[İSÄä±å£·Á^x1şè¯Ûjá“»¹Ğ¶3¼I Úx«Ãº¦‹~e¥k-Âz9£É~.Á0õë¯ƒ?¶_ˆşk/äÜj÷Ú$±7Ê>Ùi!qøâÆ?Ú¯ÛJüIı½4Ùÿ fø(µ‡¬bx¬¯®¬<UÆ8oòî“>®ğË‘×{Š"ìîT~ÅÑPØŞÁ©YÁwk*ÏmqË¨r®Œ2{EM]†EPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEP^9ûb|F_…?³Ä_ù¾MÌzL¶–­œq?î"#×"ŸÀ×±×ç_üSâwöOÃ_xŞm³k:„š•Ò©çÉ·Mª­ìÏ0#Ş*™;+nQÿ ‚!|8d±ø™ãé¢ùd’ÛC´—Ô§\ÿ Àí«õ6¾kÿ ‚tü-ÿ …Oû xÆX|›ıZÔë—{†½ÑóS#±—ş_JW!¸Wçü—àùñÂ/
|D´‡uÏ†ïÚÂñ”Ë­Î±öYQ×S_£•Àü|ø[mñ³à¿ŒürŞ™5¬Nı#›nèdÿ €È¿à4óÏü—âçü-ÏÙOÂ’\Mçj Şe²Ù€=òahNO|×Ó{E~HÁ"ş(İxãg‹¾kA¬ÛZ…¥†Şc†şĞ·™_S”Ÿúà+õÂºâîŒd¬ÄÚ(Ú)hª$M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¢–Š M¢¿?j‹«Ûş
/cà-:f—KµÔmü1FÙòá…‹^È1İXÜú ¯Õ¿Úâ¥·Á‚ş1ñÅÎÂt}:IàC…–àFOûR²/ü
¿9¿à_
®¼kñ—ÆŸµu{‘£ÛK{©y2^İÒ¸=ÙcWş»ŠÆ£èiçëİŒm½¬K­¼kQ Â¢(T{ OE QEø™ÿ ğ~£û&şİZ_ÄßÂb´Önbñ=®ß•å_“şÛÍíq_®ñVã¯	èŞ#Ò&ûN•«YÅ}k/÷¢‘©úàŠùãş
›û?·Æ¯ÙšûWÓ­ÌŞ!ğc¶³je¤·‹¨Çı³fS
õä_ğHŸCÆÿ õ?†Ú•É}[ÂrùÖAÛ-%„ÌH×Ë”°öF;V´Ş¶"K©÷İ»hÛ]BQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Ríªº®¥i¡éwšıÂZØÙÂ÷,q¢–f'Ğ Oá@›_ğYŸbĞü%ğ®Âã^¹ÖõDSÈ‰	İ³?šØõ‰M}ûüÿ …û.xCC¹·û>·¨Ãıµª«0¹¸
Ûİ#Æë~^|Ñ¯à ÿ ğP‰¼A©A$†øëPÊ	XtËfU··=¾r"FŞí×5û—\’ww7JÈ(¢Š‘…Q@ÏwPÉÑ¬ĞÈ¥92²‘‚=Aøa­Cwÿ áÿ ‚ƒI"Å4~fXÔ~Ñ¢]ÀîÆ.GlÉoé_ºuğ¿ü›öiÿ …¿ğ-|u£Úy¾&ğH{·ò×/>œØûBq×ËÂÊ3Ğ$˜å©õåõ¾§coyi4w6—¬ĞÍI€*À ‚jjøcş	?ûF‰ß¦ø}«\ù!ğhXí÷Ÿšm9îˆõòÎcöQ­}Ó´WR•ÕÌƒh§mmî!´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´QpE;h£h¢à6ŠvÑFÑEÀmí¢¢‹€Ú)ÛEE´S¶Š6Š.h§mm\ÑNÚ(Ú(¸¢´Q´Qp_ÿ ÁWÿ hdøgğF? é—;<Aã2ĞÊüĞéèA™§˜vÆê­&>í}»©jš6›u¨_O­•¬O<÷¶Ô5™˜ö 
üXĞìu/ø)gíî^U¸ù’õdÑ-›…ÏUiKí%Áì+9ÊÈ¨«³îïø$—ìô~şÏmãMRÓÉñä[Õ2(ß™Ëÿ ËË×‘"zWÜÕ••¾›g¥¤1ÛZÛÆ±EJ#EUP:  ©ëœØ(¢Š (¢Š *+«Xo­¥·¸†;‹yÇ$2¨du#H<G5-ø[ñcÃ:ßüCöàµÖ´8f—ÂSJo¬#ÎEæ—3mšÔ±ãz|Ê3Ğ¤nzŠı˜ğŠô¿x_Iñ‰v—ÚF©kå­Äg!ãu§Øàò;+Æÿ à¡²Èı¨¾]Ùé#øÏ@/©hnGÍ#…ıí¶}%P í½c'_ÿ Á&ÿ jytû«Ÿ-¸x$W–ãÃÍtv´n	iì°zî‘G¯˜?º+X=lÈ’ê~¡ÑK´Ñ´ÖæbQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´Ñ´ĞQK´×ñ£âÖ‡ğ/á–½ãoMåéºU¹—ËÜHxÏñ»£ë“À4ñWü«öŸÿ „Àvß	t½ºß‰#êïĞXâ3Ü\ßÀğâ½oş	oû-Ÿ€¿SÄúİ§“ã$w×"áílğM¼Œ©!ŒŒ=\÷|)ûü%ñíûû]jÿ |u»ğÖ™xº®¬$Ã#‹[ÿ dPGüóˆç–ı½®Y;³d¬‚Š(©QE QE QE WãÏü3ödÔ¾üWÓ><xK¦Xj—é=ì–*CiÚ²¶õ¸p«)¿ë¢¶~øû\çÄ_‡úÅok~ñ-’ê±lÖ·P‚Uº2Ÿáe8eaÈ Ò€<?ö;ı¦´ÏÚ“àæŸâHvúı¦Û=oOR3ĞQ–şyÉ÷Ôú:©¯s¯Ä}ûÆ_ğK_ÛâÃPKGÂ—$G)ÆÔÕ´§s²t<ØÈ<vtuÎÖ$şÑø[ÅO¼7¦ëúü:¦©@—6—–í”–6?Àò:k¢2º1”ljÑEb°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,QE`¢Š(Q@X(¢ŠÁEP
(¢€°QE‚Š( ,øçÿ ı¡õoÚËã¦‹ğoáÏ›ªèšn¢,¢Õ·&§©±(Ò‚81Æ(nŸë8 ×Ô¿ğS¯Û~
ø¾xVùWÆş#·e¹š1&™bÃ&Gİ’NU{¹¸!sÿ Ÿı‹ÏÃÇñ‡ÅöŠu»rº%¥ÂàÙX¸¾#´“Äxşû
Ærè\WSë¿Ù_öwÑÿ f‚ú'‚4Ï.{¸Wí¦¡oo\6^yÇTŠª;W®ÑEdhQE QE QE QE QE |åûq~ÈzOíiğm9R+Oé+%Æ©¿% n†Cÿ <¤Úô![¸??ğNÿ ÚÃVı›>#İ|ø“i%Æ ö°@ìmP-†G$àE#u=ˆnŒÆ¿fëóóş
uû·ÆÍŠ>³İã½*Û
üÚ½ª
ÖxÔqİÔmê¨)§gqqÑ_œŸğL¯Û¨øÒÎÏáÄ-Oş*;Tò´Ríşkè”qk#² )?}F>ğù¿G6×JwW3Š]´m¦QK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞQK¶´ ”Rí£m %»hÛ@	E.Ú6ĞWŠ~Öß´æ‰û+ü'»ñ> ±Şë7ÛhúI|5åÉg¸>ó°è8êÊsñâ×†¾ü>Õ|eâËõ°Ò4ø÷òO!ûÆ¿Äìxñ8 ‘øóáŸüGÿ ‚ª~Ôİ^Í.“ák4ò®^ßDÓ÷±GÙ¦“İ›,p«òÄ¥a¥s¬ı„¿e½{öÛøÙ«ü^ø¦³j~·¿77’\.¯{Á[eóÅİÀp*§oíTq¬1ª"„EUT` : +øsğóAøOàmÂ±];BÒ-ÖÚÖİNHQÉf?ÄÌIfcÉ,Ië]%sšQ@Q@Q@Q@Q@Q@Q@”ßğSØçAÔ/~7ü)³–ÜÇ/Ûµı'N$·ÆşÜ/#Ÿš@¼©ùÇ±ëğO/ÛâÚGƒÀ¾7»Šßâ=„8†åğ‹¬Â£™t(uyeãp_¿Ddu¬0U†A•ùÿ ı‚µO¾$_a¸²Ğ`¸Ú…†™•—Dœ6ï´Á·‘y «?ì}ÊŒ¹DÕÏÖ*+ãØ7ş
¥~Ò]·„|_-¾“ñ.Ö,má!ÕÕG2Â:	02Ñv^2ìı¦º¹˜”Rí4m4ÀJ)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š J)vš6š Jçü}ãíáƒõOxŸR‡IĞôØL×7S€¬Äà’@š_x÷Ãÿ ü#©xŸÅ:¤6‡§Äe¸»¸l*À¬Äà$ &¿¾>||øÿ (øÕ¥|>ø}¥]Aáh§-§élÅWh8kûÖ
 >áÚ»™¾h”¹F•È>(üIø•ÿ Cı¡ôï	øNÎkZJÇO±”‘…¾@’úí‡ÈÇ®22I-ûû7şÎ¾ı˜~Øx7Â°e#ıõö£"=ıÉ <ÒÜãz*€JÀı’dÿ şÉ¿¢ğşŠ«}­İ—X×0²ßLş;ä„Lğ	'$’}Æ°ÜĞ(¢Š@QE QE QE QE QE QE QE É¡âŠTYb‘Jº8Ê°#ÜSè ÈïÛëş	Å¨|/Ô®~.ü·¹·Ó-eû~¡¡é¥–}-ÔîûM¦Ş|°Fâ£˜ñ‘òğ¿ûÿ ÁG¬>4A§øâMÜg•D6Z¤„G³ =’öz9û¸'m~‰Wæ_íëÿ »_Üj?>Ù%¹ó]j>·é‡ÌÒÚÂKÜÇÑ¿‡ÃR“BjçèÍùcûÿ ÁN¥ğÑ¶øsñÆêá>Ìÿ d³ñMÚŸ6ß&ô›ƒÇ›Ôcçîãõ&ÎòBÖ«Y£¹¶X¦…Ã¤ˆÃ*ÊÃ‚ ‚+tÓ3±-QTEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEPEP\gÅÏŒøà{ÿ øËV‹IÒ-Şc™'–(“«»c…äàGûO~Öşı•|.º‡‰îÍÖ³tŒtİÍº¼aßÁx27¶NşWh>ø×ÿ UøÕ5íÔÿ ÙÓ¤Û%ÃşÌÑ cŸ.5ÈófaÛï1ä•QòÄ¥a¥qß¾)|Zÿ ‚¢|p³ğŸ„ôé¬<-k)’ÓK.E­„9Ú×—’àşØ€’K~²~Éÿ ²Oƒ¿dß.‹áø…ö·vªú¶½<`O} øäjIÛ8ÉÉ$ÿ Ù×öoğ_ìÃğş
ø6ÃÊC¶KíJpêş`0e™ÀäõÂŒ*ƒ€z•a¹ QE€(¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š (¢Š øßöÚÿ ‚qøOöœ³»ñ'†Ö×Â¿”oŠ¦Û}HÂ]*¼`JáÜ0 >~Öÿ àş>¹økñ7DÔ5/Z>Ù4+©›h¤œOe)Ê´mÉÛÎ
œšıÃ¯-ı ?f¿ şÓ>oøãG[ÅLµ¦£oˆï,œÿ 2à•÷S•lÀÓ>|nğWÇ¿Áâ_k–úÎœøYUÙ­¤#>\ÑŸš7‡¯Q‘ƒ]Õ~(üYıš~:ÿ Á6ü|Ş8ğ6­w¨øKÌÚºîŸh^,äC¨[ò™9RpUƒp>Ùı“à§øïö?øËìşñ¼Ÿ"Ç<¸Óï[ òecò1?òÍùèœÖŠW!£íZ)h«$J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–Š J)h ¢–¸‹Ÿ|ğ+ÂsxÆúı®…¦¦D~sf[‡>\QšGöP}N4ÛWÂß¶Gü÷ÂßÖÿ ÂŸÚ×Åş8
ÑKxOÓ$èw°ÿ [ şâœ÷˜´ü±ûEÁC>'~Ö~!?~hº¶¢jt‡NRú¶¦‡‚$d$C*§¦w9ôì_ÿ ŸÑ>µ‡‹ş0ÇkâOÆÂko¡éöG¨3q‰äŸêÁşÿ 3rìZ]Ïšÿ e¿ØOâGí¹âÕø§ñ{WÔ­|%}(™ïïúº‰n¤b8qÀ|m„SÕb>|9ğÏÂ	XøcÂ-®¡Y.ØlìÓjŒõf=Y‰ä³Iä“]q¬1ª"„EUT` : )Õ™AEPEPEPEPEPEPEPEPEPEPEPEPEPWV°ßZÍmswó!HePÈêF
°<AÆ~t~×ğHÏøóí&ø8öŞ×Û2IáÙİ6äã8„›v>œÇĞaM~Q@ˆŸ	?n¿°ÿ ‰"ğÅ-PÖtK</öF¼J]ÛÅ‘óZÜó½ è	tã
W­~¡~Ï¿µ—Ã_ÚcG^×£}E=Î‡{ˆoí¿Ş‹'pÿ i/½zÅï‚>øñáY<=ã¯YøƒN;ŒhLKnÄc|2<mşÒ{t¯ËŸÚş	!ãÏ…º³ø»àf¿u®ÛÚ9¸ƒMkŸ²jö„dşæe*²àw±ªRhV?Yè¯È?‚ğUO‰¿5£á¼ñ6Oä\O<cÖ­1ÚDp«)8F=Kšı'ø!ûSü0ı¢,Vox®ÏQ¼	¾]*cä_CëºÃ`só WĞšÑI2OU¢EP®6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê(¢Eq´S¨ .6ŠuÆÑN¢€¸Ú)ÔPE:Šãh§Q@\mê«ªjÖ:q¨jW–ú}…ºy“]]J±E¬ÌÄ =ÍrÅVÔµ;=O¸¾Ô.à±²·C$×72,qÆ£«31 êkâÚ#ş
ÑğÓáŒw:g€"?|B¹A=»tÈ[Õ¦#2õÎ#N1¼WÅ¶z?íOÿ 0ÖÖYZå|&ÿ [.í?@´ÁÚ932ŸO6AŸJ—$Š>«ı¨¿à®ğ7Úô„–°øÇ[\ÆÚíÈeÓ`luŒpÓİN„3+æ¯„ß±ÏÇÏÛÿ ÆxïâN¯£xjà‡æµ$Dÿ «±µ}	S–<º¿e¿ø%¿Ã?€¦Ó[ñ:GñÆ1áÅÖ¥ V¯ÁıÍ¹$G}Ç€@^•ö•dÛeEû;şÊÿ ¿fÿ fx#DK{¹%î³u‰o¯psûÙp8Ï!”W®ÑE 
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€<³ãŸìÃğÓöÒE—<-i«Kl·ÔPom»şîtÃvä©î~lütÿ ‚;øÓÀwÏâ‚¾(“^Kgó Òõ	–ÏS„Œ‘å\.ØänœŸ+ñ¯×š(ñ;ÀğQ¯Úö\ÖÓÂß´¯ÃlvµŸŠa’ÓRT¹Û—ûî²g<WÜ¿ÿ à§ß>.$6Ú¸ŞÖX|Ö~%Û9ï¶ä?ßd'Ò¾©ñßÃŸ
üPĞäÑ¼]áİ3Äº[ç6º¥ªN€ÿ ywµ½`Æ¾øáÿ iøqã&ÿ áÎ¹}à=A‰aas›ëzàa*}w°aj¹˜¬}ãe}o©ZEuiqÕ´Ê9¡pèêzÃ‚>•5~,j²_í…ûİMà{Zÿ H‰‹´Ş¼kÛy1Îd²a¹½~hˆó]‡Ãoø,gÄOÜe|KğEˆd·o.iíi—ªG]èU°ô
•|ÈV?]è¯>ÿ ÁT¾|DòaÔu«ïß¿Gˆ-
Ç»¾&ˆ¼`{±_Â¾ ğ‡Ä/üA²ûg…üI¤ø×ó´›è®PrŒj®#~ŠZ)€”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E- ”RÑ@	E-6GX£gvŠ2ÌÇ S@ExÏÄÛ+à§Âxåÿ „‹âF‡ÌyİecqöÛ}Poqø_$|Pÿ ‚Ñx3Hó­üàWÄS•o5‰’Æÿ x"ùÃØì4®†~×üWı¡>ü²ûO|e¥xt”2%­Äá®¥QÔÇæGí÷Tõù1{ûY~×¿¶5äÚw€ìµk.F(Ğø6É­`<bKÖ%ãÖUÒ»ï…ğF¿ˆ6¿]_â·-<<'2âÖÅ£!îRDjßí“éQÍØ,vŸÿ à³šeˆ¸ÓşøBMJnUu¯“ ú¥ºÌ=:9Zğ='àÏígÿ Õ Ô¼E>¥oág“ÌŠó^'OÒaRGÍº¨ó8ãr#™³Í~ üı€ş	|û=Î‡á5}n­xƒ·a‡ñ.á²3ï/Zú"¡¶Æ|3û=Á$¾ü)6š§_‰^ ‹WPÊÓcl°'Ìïşµ˜Ñ_oÙY[é¶pÚZA­¬#ŠP"F `*¨à ;
ŠC
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
(¢€
á¾$|ø{ñ‚ØÃã_è%Ú²ê1É4îIè}ÔŠîh ‚ş)Á~øÁ¥¸ğ¡®xí²V+yşÛh3ëÙ“ò”WËş0ÿ ‚7|bğ5çö€¼m¢kïÌL²Í¥Şû2ûú+öRŠ üI]SöùısøêúÎß‚e·OÀÓLOµ~Œ1íZšüãg‚ï>ÁãxT’/õ‘ÜYÏauøáÊû÷_´U•â
è-³û&¹£Øk6¿óÃPµIÓşùpE;±XüÒğßüÓÃwˆ>jšyş6ÓuXî¿8¿,×©øwş
ûğZUûoü$ú=~ß¥«ÿ ~d’½ëÄ¿°¿ÀoûwÂ_D_©Óì…‰üà)^Uâ/ø$Ÿìé®3Oë{iÚÌìış2Sæadli?ğSÙ»WÚâDVÎ†ïJ¾‡‹CÖºë/Û‹à¡x¾,øaAÿ ×¢#ù>|Ÿñcş	ğwÁºy¼ÓüAã]ÅKyrßÚ2ü•Ïë_üAı˜ü/á=J{{Kı^D¡šh‰ı#ù˜Xıƒö¾øq÷>/x$×MzÙ?›Šµÿ YğOş‹€ğ§²ÿ ãµüùx“Á6:>ÿ &[†Ûÿ =OòQ\SpH£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…éş³àŸı ÿ áOeÿ Çhÿ †¬ø'ÿ E‡À?øSÙñÚşoh£™…èş_ÚËàŒ+–øÁà2?Ùñ%›)+>Û;àM¾wü\ğyÿ ®z¼/ü˜×óÃ¢éÑê7$fQşÉá^Ÿá¿ƒÚ6±·Î¹¾\ÿ Ï944s0±ûa«ÁA¿gwÚ>)èòcş}{ıW¬ÁUg--\Ûøºÿ V+ü6z-Ø'éæF‚¾"ø#ÿ ÿ øyñ$Áı§¬øš3û%Õºÿ èP5}u ÿ Á¾XÇ—7Ş1Õr»Õ!P}¿w\Ì,rŞ$ÿ ‚Ğ|(°Üº'„<[«ºôkˆí­co¡ó]¿5¯#ñüÓ]¸VO|.ÓìøfÖ5I.÷)qãé¸×Úÿ ‚g~ÍŞØÑ|7·Ô&^²jZ…İÎïª¼¥?%¯\ğìëğ¯À.’xsáÇ…tY×¤öz=¼rş.qüMÌ,ÈïøoÛãËà­ú&ÿ ¡KÂÍ:`ÿ ÓIR£ßpúÒÃûşÙ_´„‹7õJÒÊcŞ1ñyH=­ãi·–+öÌqÀàRÒ»ùeğãş…lŞ>ø™,ËÇ™cáË÷Äóÿ ÑUõÏÂßø'OÀ„şL¶> ²×oãÁûwˆÉÔ°èÛ$Ìj}Õ})E  ±±¶Óm"µ³·ŠÒÖ%Û0 DAèpOE QE QE QE QE QE QE QE QE QE QE QEÿÙawk '{i=1;while(i <= NF){col[i]=col[i] $i " ";i=i+1}} END {i=1;while(i<=NF){print col[i];i=i+1}}' inter_result/snpmatrix.txt |sed 's/IID //g' > inter_result/snpmatrix2a.txt
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
)‹      r‰0æŠàb```f`afd`f2XCCÜt-X˜€FN$ø@û0ğ à'pA   shinyUI(fluidPage(
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
