# Function ----------------------------------------------------------------

omicQTL_fun=function(omic,type,gtf,vcf,norm,out,threshold){
  source("code/All_basic_function.R")
  osca="../../home/software/osca-0.46.1-linux-x86_64/osca-0.46.1"
  suppressMessages(library(data.table))
  suppressMessages(library(MatrixEQTL))
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
    message("Z-score phenotype data")
    data=cbind.data.frame(data[,1:2],
                          apply(data[,3:ncol(data)],2,zscore))
  }
  
  ##matrix eQTL analysis
  suppressMessages(library(MatrixEQTL))
  ##make exp file
  message("Converting phenotype format to MatrixEQTL format")
  data=data.frame(t(data))
  data=cbind.data.frame(row.names(data),data)
  data=data[-1,]
  names(data)=c("gene",as.character(as.numeric(data[1,-1])))
  data=data[-1,]
  ###separate the data into five file to parallel analysis
  mark=split(1:nrow(data),cut(seq_along(1:nrow(data)),breaks=2,labels = F))
  for(i in 1:2){
    sub_data=data[mark[[i]],]
    write.table(sub_data,file=paste0(inter,"omic_",i,".txt"),
                row.names = F,col.names = T,quote=F,sep="\t",eol = "\n")
  }
  
  ##make snp file
  suppressMessages(library(dplyr))
  #BiocManager::install("VariantAnnotation",type="binary")
  #BiocManager::install('snpStats')
  suppressMessages(library(VariantAnnotation))
  suppressMessages(library(snpStats))
  #convert.vcf(vcf.file=vcf,genotype_file_name = paste0(inter,"snp_for_eqtl.txt"))
  write.table("maked",file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  ##make cvrt(PC1-40) file 
  plink="softwares\\plink\\plink.exe"
  cmd_cvrt=paste("../../home/software/plink","--vcf",vcf,"--pca 20","--out",paste0(inter,"pca_40"))
  write.table(cmd_cvrt,file="info.txt",append=T,row.names = F,col.names = F,quote=F)
  message("Performing PCA analysis with genetic data")
  system(cmd_cvrt)
  pca=read.table(file=paste0(inter,"pca_40",".eigenvec"))
  row.names(pca)=pca[,2]
  pca=t(pca[,-1:-2])
  row.names(pca)=paste(rep("pc_",20),1:20)
  write.table(pca,file=paste0(inter,"pca_40_eqtl.txt"),quote = F,sep="\t")
  
  message("Loading data")
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

  message("Converting PCA format to MatrixQtl format")
  message("Selecting 20 PC as covariant")
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
  suppressMessages(library(parallel))
  message("Performing Omic QTL mapping")
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
  message("The Omic QTL result was saved in qtls.txt")
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
  
  if(grepl("gff",gtf)){
    if(type=="Transcriptom"){
    suppressMessages(library(stringr))
    system(paste("/usr/bin/python3 code/gff_format.py",gtf,paste0(inter,"qtlgff.gff3")))
    gtf=data.frame(fread(paste0(inter,"qtlgff.gff3"),fill=T))
    message("As you provided gff or gff like file, cis-trans QTL plot will be plotted")

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
    suppressMessages(library(stringr))
    system(paste("/usr/bin/python3 code/gff_format.py",gtf,paste0(inter,"qtlgff.gff3")))
    gtf=data.frame(fread(paste0(inter,"qtlgff.gff3"),fill=T))
    message("As you provided gff or gff like file, cis-trans QTL plot will be plotted")
    if(sum(gtf$V3=="gene")>1){
      gene=data.frame(str_split_fixed(str_split_fixed(gtf$V9[which(gtf$V3=="gene")],";",2)[,1],":",2)[,2])
      gene$chr=gtf$V1[which(gtf$V3=="gene")]
      gene$BP=(gtf$V4[which(gtf$V3=="gene")]+gtf$V5[which(gtf$V3=="gene")])/2
    }else{
      gene=gtf$V9
      gene$chr=gtf$V1
      gene$BP=(gtf$V4+gtf$V5)/2
    }
    
    
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
  }
  chr_mean=aggregate(BP~chr,data=gene,mean)
  chr_mean=chr_mean[chr_mean[,1]%in%unique(vcf_data[,1]),2]
  #chr_mean=chr_mean+max(chrlen$BP)
  
  
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
        snp_pos[snp_pos[,1]==i,2]=snp_pos[snp_pos[,1]==i,2]+sum(chrlen[1:(i-1),2])
      }      
    }
  } 
  snp_pos=snp_pos[,2]
  ##plot result
  message("Plotting cis-trans QTL plot")
  pdf(paste0(out,"cis-trans_plot.pdf"),width=10,height=10)
  plot(snp_pos,gene_pos,
       frame.plot=F,pch=20,xaxt="n",yaxt="n",ylab="Gene position",xlab="QTL position",col="skyblue")
  axis(1,chr_mean,unique(vcf_data[,1]),las=1,cex.axis=1)
  axis(2,chr_mean,unique(vcf_data[,1]),las=1,cex.axis=1)
  dev.off()
  
  png(paste0(out,"cis-trans_plot.png"),width=10,height=10,units="in",res=600)
  plot(snp_pos,gene_pos,
       frame.plot=F,pch=20,xaxt="n",yaxt="n",ylab="Gene position",xlab="QTL position",col="skyblue")
  axis(1,chr_mean,unique(vcf_data[,1]),las=1,cex.axis=1)
  axis(2,chr_mean,unique(vcf_data[,1]),las=1,cex.axis=1)
  dev.off()

  message("Plotting QTL occurrence numbe plot")
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
  }else{
    message("As you don't provide a gff or gff like file, only output the qtl summary infomation")
  }
}


# Code --------------------------------------------------------------------
source("code/All_basic_function.R")
suppressMessages(library(parallel))
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
message("Converting VCF format to MatrixEQtl format")
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
#tryCatch({omicQTL_fun(omic=omic,type=type,gtf=gtf,vcf=vcf,norm=norm,out=out,threshold=threshold)}, warning = function(w){message("")})

#data_convert_fun(vcf)
