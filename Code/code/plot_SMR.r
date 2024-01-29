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
    par(mar=c(5,5,5,5), xpd=TRUE)
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
    par(mar=c(5,5,5,5))
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
   par(mar=c(5,5,5,5))
   pal=colorRampPalette(c(rgb(0, 0, 1),  rgb(0, 0, 0)))
   breaks <- seq(min(snpCorr^2), max(snpCorr^2),length.out=2000)
   #image.scale(snpCorr^2, col=pal(length(breaks)-1), breaks=breaks, horiz=FALSE, xlab="", yaxt="n")
   dvd = (max(snpCorr^2) - min(snpCorr^2))/5
   pos = seq(min(snpCorr^2), max(snpCorr^2), dvd)
   axis(4, at=pos, label=sprintf("%.2f", pos), las=2)
   mtext(expression(italic(r)^2), side=4, line=3.5, las=2 )
}

