suppressMessages(library(data.table))
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

options (warn = -1)
# Display the manhattan result --------------------------------------------
##This function was used to illustrate the single trait, cojo result, and multi-trait result
plot_fun=function(result,out,show_peakloci=T,color_manh=c("tomato","skyblue"),threshold=5e-8,corrected=T){
  if(!file.exists(result)){
    message("No significant MR result")
    message("Maybe the over setted threshold")
    message("Please Check the .log file")
  }
  data=data.frame(fread(result))
  message("Plotting result")
  ##check the file format to continue right downstream analysis
  if(ncol(data)==10){
    #For single trait and multi-trait
    names(data)=c("CHR","SNP","POS","A1","A2","AF1","BETA","SE","P","N")  
  }else if(ncol(data)!=10){
    message("loading cojo result")
    #For cojo
    data=data[,c(1,2,3,4,4,5,11,12,13,9)]
    #Note, the A2 col isn't the real A2, because the cojo result file don't contain this cols
    names(data)=c("CHR","SNP","POS","A1","A2","AF1","BETA","SE","P","N")
  }
  data=na.omit(data)  
  data[data$P<=1e-50,"P"]=1e-50
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
  message("Plotting manhattan")
  pdf(file=paste0(out,"manhattan.pdf"),width=12,height=6)
  plot(x=data[,3],y=-log10(as.numeric(data[,9])),
       cex=1,pch=20,col=cols,frame.plot=F,xaxt="n",yaxt="n",
       xlab="Chromosome",ylab=expression("-Log"["10"]*"(P)"))
  chr_mean=aggregate(POS~CHR,data=data,mean)[,2]
  ##show the threshold
  if(threshold=="Bonferroni"){
     #threshold can calculated
    threshold=0.05/nrow(data)

    message(paste0("Assigning threshold as ",-log10(threshold)," based on Bonferroni"))
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
  axis(2,seq(0,max(-log10(as.numeric(data[,9])),na.rm=T),2),seq(0,max(-log10(data[,9]),na.rm=T),2),las=1,cex.axis=1)
  dev.off()

  png(file=paste0(out,"manhattan.png"),width=12,height=6,units="in",res=600)
  plot(x=data[,3],y=-log10(as.numeric(data[,9])),
       cex=1,pch=20,col=cols,frame.plot=F,xaxt="n",yaxt="n",
       xlab="Chromosome",ylab=expression("-Log"["10"]*"(P)"))
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
  
  if(show_peakloci==T&length(snps)!=0){
    points(top_data$POS,-log10(top_data$P),cex=1,pch=17,col="red")
    text(top_data$POS-0.04*max(data$POS,na.rm=T),-log10(top_data$P),top_data$SNP,cex=0.8)
  }
  
  axis(1,chr_mean,c(unique(data$CHR)),las=1,cex.axis=1)
  axis(2,seq(0,max(-log10(as.numeric(data[,9])),na.rm=T),2),seq(0,max(-log10(data[,9]),na.rm=T),2),las=1,cex.axis=1)
  dev.off()
  
  message("Plotting QQ-plot")
  ###QQplot
  pdf(file=paste0(out,"qqplot.pdf"),width=8,height = 8)
  lambda=qqplot_fun(data$P,plot=T,frame.plot=F,pch=20,col="skyblue")
  dev.off()
  png(file=paste0(out,"qqplot.png"),width=8,height = 8,units="in",res=600)
  lambda=qqplot_fun(data$P,plot=T,frame.plot=F,pch=20,col="skyblue")
  dev.off()
  
  ### corrected manhattan
  if(lambda$estimate>1.1 & corrected==T){
    message("The inflation factor is bigger than 1.1")
    message("Performing p-value correction")
    data$P=pchisq(qchisq(data$P,df=1,lower.tail=F)/lambda$estimate,df=1,lower.tail = F)
    write.table(data,file=paste0(out,"mlm_corrected.mlma"),
                quote = F,col.names = T,row.names = F)
    ##The below code are same as above manhattan. 
    ##If changed the above code,can directely copy to here, except the filename.
    message("Plotting corrected manhattan ")
    pdf(file=paste0(out,"manhattan_corrected.pdf"),width=12,height=6)
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
    axis(2,seq(0,max(-log10(as.numeric(data[,9])),na.rm=T),2),seq(0,max(-log10(data[,9]),na.rm=T),2),las=1,cex.axis=1)
    dev.off()

    png(file=paste0(out,"manhattan_corrected.png"),width=12,height=6,units="in",res=600)
    plot(x=data[,3],y=-log10(data[,9]),
         cex=1,pch=20,col=cols,frame.plot=F,xaxt="n",yaxt="n",
         xlab="Chromosome",ylab=expression("-Log"["10"]*"(P)"))
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
    
    if(show_peakloci==T&length(snps)!=0){
      points(top_data$POS,-log10(top_data$P),cex=1,pch=17,col="red")
      text(top_data$POS-0.06*max(data$POS,na.rm=T),-log10(top_data$P),top_data$SNP,cex=0.8)
    }
    
    axis(1,chr_mean,c("1","2","3","4","5"),las=1,cex.axis=1)
    axis(2,seq(0,max(-log10(as.numeric(data[,9])),na.rm=T),2),seq(0,max(-log10(data[,9]),na.rm=T),2),las=1,cex.axis=1)
    dev.off()
  }
}



snp_finder=function(gctadata,threshold){
  message("Finding leading SNP of GWAS peak")
  snps=c()
  gctadata=na.omit(gctadata)
  window <- seq(0,max(gctadata$POS),5e04)
  gctadata=gctadata[gctadata$P<threshold,]
  if(nrow(gctadata)==0){
    message("No SNP passed the threshold")
    return(c())
  }else if(nrow(gctadata)==1){
    message("Finished leading SNP finding")
    return(gctadata$SNP)
  }else{
    message("Finished leading SNP finding")
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



