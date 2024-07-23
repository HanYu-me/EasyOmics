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
  message("Estimating the kinship h2 of each trait")
  cmd_grm=paste(gcta,"--bfile",paste0(inter,"bfile"),"--make-grm  --thread-num 10 --out",paste0(inter,"grm"))
  system(cmd_grm)
  summ=data.frame(va=1,vp=1,h2=1)[-1,]
  for(i in 1:(ncol(phe_data)-2)){    
    cmd_h2=paste(gcta,"--reml","--grm",paste0(inter,"grm"),"--pheno",phe,"--mpheno",i,"--out",paste0(inter,"h2_",i))
    system(cmd_h2)
    h2=data.frame(fread(paste0(inter,"h2_",i,".hsq"),fill=T))
    h2=h2[c(1,3,4),2]
    summ[i,]=h2
    message(paste0("Trait ",i," estimated"))
  }
  write.table(summ,file="info.txt",append=T,row.names = F,col.names = F,quote=F) 
  cols=colorRampPalette(brewer.pal(12, 'Set3'))(length(3:ncol(phe_data)))
  write.table(paste0(out,name,"phenotype_summary.png"),file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  kk=length(3:ncol(phe_data))
  if(kk>20){kk=20}
  png(file=paste0(out,name,"phenotype_summary.png"),width=kk*0.6,height=13*0.6,units = "in",res=600)
  par(mfrow=c(4,1))
  par(mar=c(3,5,3,1))
  boxplot(phe_data[,3:ncol(phe_data)],col=cols,frame.plot=F,ylab="Phenotype value")
  barplot(summ$vp,names=names(phe_data)[3:ncol(phe_data)],col=cols,ylab="Phenotypic variance")
  barplot(summ$va,names=names(phe_data)[3:ncol(phe_data)],col=cols,ylab="Additive variance")
  barplot(summ$h2,names=names(phe_data)[3:ncol(phe_data)],col=cols,ylab="Kinship heritability",ylim=c(0,1))
  write.table("what",file="info.txt",append=T,row.names = F,col.names = F,quote=F) 
  dev.off()
  pdf(file=paste0(out,name,"phenotype_summary.pdf"),width=kk*0.6,height=13*0.6)
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
  message("The genetic information was saved in Phe_info.txt")
  
  message("The summary of traits was saved in Rawtra_phenotype_summary.pdf")
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
  png(file=paste0(out,name,"phenotype_cor.png"),width=11.5*0.6,height=10*0.6,units = "in",res=600)
  pheatmap(phe_cor,cluster_rows = F,cluster_cols = F,na_col = "white",
           border_color = NA)
  dev.off()
  pdf(file=paste0(out,name,"phenotype_cor.pdf"),width=11.5*0.6,height=10*0.6)
  pheatmap(phe_cor,cluster_rows = F,cluster_cols = F,na_col = "white",
           border_color = NA)
  dev.off()
  if(name=="Plasticity_"){
    write.table(phe_cor,file=paste0(out,"Phe_Pla_cor.txt"),quote = F,col.names = T,row.names = T)  
  }else{
    write.table(phe_cor,file=paste0(out,"Phe_cor.txt"),quote = F,col.names = T,row.names = T)  
  }
  message("The correlation between each traits was saved in Phe_cor.txt and Rawtra_phenotype_cor.pdf")
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
suppressMessages(library(lme4))
suppressMessages(library(tidyr))
suppressMessages(library(data.table))
suppressMessages(library(pheatmap))
suppressMessages(library(png))
suppressMessages(library(RColorBrewer))
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
}else if(type=="Multi_Traits"){
  genetic_data_fun(phe=phe,vcf=vcf,name="Rawtra_",out=out)
}else if( type=="Single_Trait"){
  phe_data=data.frame(fread(phe))
  inter="inter_result/"
  gcta="../../home/software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"
  ##estimate h2
  message("Performing PCA analysis with genetic data")
  cmd_grm=paste(gcta,"--bfile",paste0(inter,"bfile"),"--make-grm  --thread-num 10 --out",paste0(inter,"grm"))
  system(cmd_grm)
  cmd_pca=paste(gcta,"--grm", paste0(inter,"grm"),"--pca 2 --out",paste0(inter,"pca_single"))
  system(cmd_pca)
  message("Estimated the kinship h2 of the single trait")
  cmd_h2=paste(gcta,"--reml","--grm",paste0(inter,"grm"),"--pheno",phe,"--mpheno",1,"--out",paste0(out,"h2_single_",1))
  system(cmd_h2)
  
  ##load pca
  pca=data.frame(fread("inter_result/pca_single.eigenvec"))
  eig=data.frame(fread("inter_result/pca_single.eigenval"))

  pc1=pca[,3]
  pc2=pca[,4]
  explainability=c(eig[1:2,1])
  
  if(sum(phe_data[,3]%in%c(1,0,-9))!=nrow(phe_data)){
    message("This single trait is a qauantitative trait")
    pdf(file=paste0(out,"Single_Trait_result.pdf"),width=15*0.6,height=7*0.6)
    par(mfrow=c(1,2))
    plot(density(phe_data[,3],na.rm=T),frame.plot=F,main="Phenotype distribution",lwd=1.5,col="skyblue")
  }else{
    message("This single trait is a qualitative trait")
    pdf(file=paste0(out,"Single_Trait_result.pdf"),width=8*0.6,height=7*0.6)
  }
  plot(pc1,pc2,col="skyblue",frame=F,
    xlab=paste("PC1 ",explainability[1]),ylab=paste("PC2 ",explainability[2]),
    main="Genetic structure",pch=20)
  dev.off()

  if(sum(phe_data[,3]%in%c(1,0,-9))!=nrow(phe_data)){
    png(file=paste0(out,"Single_Trait_result.png"),width=15*0.6,height=7*0.6,units="in",res=600)
    par(mfrow=c(1,2))
    plot(density(phe_data[,3],na.rm=T),frame.plot=F,main="Phenotype distribution",lwd=1.5,col="skyblue")
  }else{
    png(file=paste0(out,"Single_Trait_result.png"),width=8*0.6,height=7*0.6,units="in",res=600)
  }
  plot(pc1,pc2,col="skyblue",frame=F,
    xlab=paste("PC1 ",explainability[1]),ylab=paste("PC2 ",explainability[2]),
    main="Genetic structure",pch=20)
  dev.off()

  h2=data.frame(fread(paste0(out,"h2_single_1.hsq"),fill=T))
  se=h2[4,3]
  h2=h2[4,2]
  message(paste0("This trait with a h2 ",h2," (se=",se,")"))
}




#phe="data/input/phenotype/phe.txt"

