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

  png(file=paste0(out,"TWAs_manhattan.png"),width=10,height=6,units = "in",res=600)
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


