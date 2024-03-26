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
