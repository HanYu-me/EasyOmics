##blup
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
pheno # this is the paireise difference