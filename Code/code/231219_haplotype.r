vcf="../vcf.vcf"
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
