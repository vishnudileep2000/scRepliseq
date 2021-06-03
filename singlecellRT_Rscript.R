## singlecellRT_Rscript.R

#1.	Load R packages.

setwd("/Path/To/Scripts_and_binnedData")  # set path to scripts and binned data

source("./scRT_functions.R")
library(copynumber)
library(zoo)
library(circlize)
library(gridExtra)
library(lattice)
library(ggplot2)
library(mixtools)
library(fBasics)

#2.	Identify cells with complex karyotype aberrations or complete chromosome loss by plotting whole-genome coverage at 1Mb bins.

RT=read.delim("scRT_96cells_mm10_1Mb.txt", header=T)
RT=RT[RT$CHR!="chrM",]

for(i in 3:ncol(RT)) RT[,i]=RT[,i]* (1000000/sum(RT[,i])) #Normalinze to 1 million reads#

#Calculate chromosome lengths#

chrs=levels(as.factor(as.character(RT$CHR)))
chrlen=NULL
for(chr in chrs){
  RTb=subset(RT,CHR==chr)
  temp=data.frame(CHR=chr,Len=nrow(RTb))	
  chrlen=rbind(chrlen,temp)
}
chrlen$c.Len=cumsum(chrlen$Len)

#Plot all datasets and visually identify cells with complex karyotype changes#

RTb=RT
par(mar=c(0,2,0,2),mfrow=c(10,1),ask=T, omi = c(.3,0,0.1,0))
for(i in 3:ncol(RTb)){
  plot(RTb[,i],type="p",pch=19,cex=0.2,ylim=c(-20,800))
  abline(v=chrlen$c.Len,col="red")
  legend("topleft",names(RTb)[i],text.col="red",bty="n")	
}

complexkaryotype=c("C5","C8","G1","G6","H2","H4")

#3.	Function to Segment genome based on copynumber

segment_vd=function(data,my_gamma=3,my_kmin=10,my_fast=F){
  pcf.seg=NULL
  for(chr in 1:20){
    cat("\nCHR=",chr)
    segmentation=NULL
    datab=subset(data,CHR==chr)
    datab[,2:ncol(datab)]=na.approx(datab[,2:ncol(datab)],na.rm=F)
    datab[,2:ncol(datab)]=na.locf(datab[,2:ncol(datab)],na.rm=F)
    datab[,2:ncol(datab)]=na.locf(datab[,2:ncol(datab)],na.rm=F,fromLast = T)
    datab$CHR=NULL
    
    segmentation=pcfPlain(datab,gamma=my_gamma,kmin=my_kmin,return.est = F,fast=my_fast)
    segmentation$chrom=chr
    pcf.seg=rbind(pcf.seg,segmentation)
  }
  pcf.seg=pcf.seg[order(pcf.seg$sampleID,pcf.seg$chrom,pcf.seg$start.pos),]
  return(pcf.seg)
}


#4.	Load the 50kb binned single-cell data frame

RT=read.delim("scRT_96cells_mm10_50kb.txt", header=T)

#5.	Filter cells with low sequencing reads (<250,000 reads/cell) and complex karyotype

read=stack(colSums(RT[,3:ncol(RT)]));names(read)=c("mapped_reads","cell")
remove=read[read$mapped_reads<=250000,] 
remove=as.character(remove$cell)

remove=c(remove,complexkaryotype) # Add cells with complex karyotypes#
RT=RT[,!names(RT) %in% remove] ## removes all cells identified above from fur-ther analysis#

#6.	Convert read counts of all cells to reads per million (RPM) to control for sequencing depth variability

for(i in 3:ncol(RT)) RT[,i]=RT[,i]* (1000000/sum(RT[,i]))

#7.	Remove Mitochondria DNA and Y chromosome reads

RT=RT[RT$CHR!="chrM" & RT$CHR!="chrY",] 

#a.	Re-format single-cell data for copynumber package
RT[,1]=as.character(RT[,1])
RT[RT$CHR=="chrX",1]="chr20"
RT$CHR=as.numeric(substring(RT$CHR,4,6))

#8.	  Read population bulk replication timing dataset
RTbulk=read.delim("cas.129.bulkRT.txt", header=T)
names(RTbulk)[3]="cast.129"
RTbulk=RTbulk[RTbulk$CHR!="chrM" & RTbulk$CHR!="chrY",]

#a.	Re-format bulk data for copynumber package
RTbulk[,1]=as.character(RTbulk[,1]) ##Re-format data for copynumber package###
RTbulk[RTbulk$CHR=="chrX",1]="chr20"
RTbulk$CHR=as.numeric(substring(RTbulk$CHR,4,6))

#9.	Aligning bulk and single cell replication timing data
RT=merge(RT,RTbulk,by=c("CHR","POS"))
RT=RT[order(RT$CHR,RT$POS),]
RTc=RT
RTc$cast.129=NULL

#10.	Identify bins with extreme values (mean RPM>99 percentile and mean RPM<1 per-centile) and mask them in all single cells.

#a.	Load sorted G1/G2 control cells

g_con=RTc[,c("CHR","POS","H5","H6","H8","H9")]
g_con$mean=rowMeans(g_con[,3:6])
g_con=g_con[,c("CHR","POS","mean")]

#b.	Identify outlier 50kb bins based on G1/G2 data

badwins=subset(g_con, mean<=quantile(g_con$mean,0.05) | mean>=quantile(g_con$mean,0.99))
badwins=badwins[,c("CHR","POS")]
badwins$ID=as.character(paste(badwins$CHR,"_",badwins$POS,sep=""))

#c.	Identify outlier segments based on G1/G2 data (ie. multiple consecutive 50kb bins that have unusually low values)

gseg=segment_vd(g_con[,c("CHR","POS","mean")])
gsegframe=interpolate.pcf(gseg,g_con[,1:2])
names(gsegframe)=c("CHR","POS","mean")
goodsegs=subset(gsegframe,mean>quantile(gseg$mean,0.05))
goodsegs$mean=NULL

#d.	Mask outlier bins regions with NA in the single cell data. 

RTc$ID2=as.character(paste(RTc$CHR,"_",RTc$POS,sep=""))
L=RTc$ID2 %in% badwins$ID
RTc[L,3:ncol(RTc)]=NA
RTc$ID2=NULL

#e.	Remove outlier segments from single cell dataset.

RTc=merge(RTc,goodsegs,by=c("CHR","POS"))
RTc=RTc[order(RTc$CHR,RTc$POS),]


#11.	Normalize read counts of single cells by dividing coverage data of each single cell by the coverage of G1 and G2 control cells.

base=(RTc$H5 + RTc$H6 + RTc$H8 + RTc$H9)/4
RTc$base=base
for(i in 3:ncol(RTc)) RTc[,i]=RTc[,i]/RTc$base
RTc$base=NULL
RTc[,3:ncol(RTc)]=do.call(data.frame,lapply(RTc[,3:ncol(RTc)], function(x) replace(x, !is.finite(x),NA))) #replace non-finite with NA


#12.	Center and scale data from all single cells to give them an equal interquartile range.

for (i in 3:ncol(RTc)) RTc[,i]= scale(RTc[,i],scale=FALSE)

rescaleDatasets = function(RT) {
  for (i in 3:ncol(RT)) { # 4:ncol(RT) if probeID present
    RT[,i] = RT[,i] * 1.59 / IQR(RT[,i],na.rm=T) # Scale to BG01ESC.R1 IQR
    cat("Scaling dataset", i-2, "/", ncol(RT)-2, "\n")
  }
  return(RT)
}


RTc=rescaleDatasets(RTc)

remove=c("H5","H6","H8","H9")
RTc=RTc[,!names(RTc) %in% remove]  ###Remove control cells.

#13.	Smooth data by applying a median filter with a span of 15 windows. Perform median filtering one chromosome at a time.

chrs=levels(as.factor(as.character(RTc$CHR)))
AllData=NULL

for(chr in chrs){
  cat("\n",chr,"\n")
  chrData=RTc[RTc$CHR==chr,1:2]
  full=RTc[RTc$CHR==chr,1:2]
  full$ID2=as.character(paste(full$CHR,"_",full$POS,sep=""))
  
  for(i in 3:ncol(RTc)) {
    cat(" ",names(RTc)[i]," ")
    med.temp=RTc[RTc$CHR==chr,c(1,2,i)]
    med.temp=na.omit(med.temp)
    up.lim=quantile(med.temp[,3],0.99,na.rm=T)
    down.lim=quantile(med.temp[,3],0.01,na.rm=T)
    med.temp[med.temp[,3]>=up.lim | med.temp[,3]<=down.lim,3]=NA
    med.temp=na.omit(med.temp)
    med.temp$smooth=runmed(med.temp[,3],15)
    med.temp$ID=as.character(paste(med.temp$CHR,"_",med.temp$POS,sep=""))
    L=!full$ID2 %in% med.temp$ID
    med.temp$ID=NULL
    omitted=full[L,c(1,2)]
    add=cbind(omitted[,1:2], matrix(data = NA, nrow = nrow(omitted), ncol = 2))
    names(add)=names(med.temp)
    med.temp=rbind(med.temp,add)
    med.temp=med.temp[order(med.temp$CHR,med.temp$POS),]
    chrData=cbind(chrData,med.temp[,4])
  }
  AllData=rbind(AllData,chrData)
}
names(AllData)=names(RTc)
RTs=AllData

#14.	Segment S-phase single cell data using copy number package
pcf.seg=segment_vd(RTs,my_fast = F,my_gamma = 3,my_kmin = 5)
segframe=interpolate.pcf(pcf.seg,RTs[,1:2])
names(segframe)[1:2]=c("CHR","POS")

## Check segmentation by plotting few cells using custom function in scRT_functions.R file

check.seg(RTs,segframe,samples=c("B9","B3","B2","A11"),binary=F,chrs=c(1,2,16),x1=2e6,x2=200e6)

# 15.	Binarize single cell data using

copyframe=segframe

#function to plot results of mixture model
plotmix=function(linep){     
  d=density(test,adjust=0.5)
  xa=min(d$x)-0.1; xb=max(d$x)+0.1
  ya=0; yb=max(d$y)+0.1
  
  plot(d,xlim=c(xa,xb),ylim=c(ya,yb))
  abline(v=linep)
  par(new=T)
  d=density(rnorm(1000000, mean=model$mu[1],sd=model$sigma[1]),adjust=0.5)
  plot(d$y*model$lambda[1]~d$x,type="l",col="red",xlim=c(xa,xb),ylim=c(ya,yb))
  d=density(rnorm(1000000, mean=model$mu[2],sd=model$sigma[2]),adjust=0.5)
  lines(d$y*model$lambda[2]~d$x,col="green",type="l")
  legend("topright",paste(sample,round(mean.diff,2),round(skew,2),sep=" "),bty="n")
}


#Test 100 different thresholds to find best threshold for binarization for each single-cell#

c.noskew=NULL;c.lskew=NULL;c.rskew=NULL;
samples=names(copyframe)[3:ncol(copyframe)]       #list of single-cells
point2Cs=NULL                                     #vector to store the best threshold for each cell
for(sample in samples){                           #Iterate over each sample
  seg=NULL;model=NULL
  cat("\n",sample,"\n")
  test=segframe[,which(names(segframe)==sample)]
  test=test[abs(test)<2]                      #mask bins with extreme segment-ed values for binarization threshold calculation
  model <- normalmixEM(x=test, k=2)               #Mixture model with 2 compo-nents (Replicated and un-replicated)
  
  skew=skewness(test)      #Mixture models fails in cells in very early or late is S-phase. Calculating skew in the data and difference of component means in the mixture model to identify these cells.
  mean.diff=abs(diff(model$mu))    
  cat("\n","skew - ",skew)                       
  cat("\n","Diff - ",mean.diff)
  
  #testing=rbind(testing,data.frame(sample,mean.diff,skew))
  
  #Setting range of the thresholds to iterate over#
  sweep.min=min(test)
  sweep.max=max(test)
  
  cpoints=seq(sweep.min,sweep.max,(sweep.max-sweep.min)/100)  #Thresholds to iterate over
  cpoints.cors=NULL
  seg=data.frame(test=test,testb=NA,testc=NA)
  
  if(mean.diff>0.7){       #If difference in component means are above 0.7 (may need to be adjusted empirically) then there is a clear replicated and un-replicated fraction in the data
    plotmix(model$mu)
    for(cpoint in cpoints ){     #Iterate over all Thresholds
      seg[,"testb"]=seg[,"test"]-cpoint        #subtract threshold from data##
      seg[seg$testb<0,"testc"]=min(model$mu)           #Binarize with compo-nent means instead of 1s and 0s
      seg[seg$testb>0,"testc"]=max(model$mu)
      cpoints.cors=c(cpoints.cors,sum(abs(seg$test-seg$testc)))  ##Calculate Manhattan distance between Binarized and Un-brinarized data
    }
  }
  
  
  
  if(mean.diff<0.7){     ##If difference in component means are below 0.7 (maybe need to be adjusted empirically) then the cell is most likely in very early or very late S-phase. Classify based on skew of data
    if(skew<(-0.2)){ # Cell is in late S-phase, binarize using 5 and 50 per-centile of the segmented data
      c.lskew=c(c.lskew,sample)
      cat("\n","skew to left-late S")
      plotmix(quantile(test,c(0.05,0.5)))
      for(cpoint in cpoints ){
        seg[,"testb"]=seg[,"test"]-cpoint
        seg[seg$testb<0,"testc"]=quantile(test,0.05)
        seg[seg$testb>0,"testc"]=quantile(test,0.5)
        cpoints.cors=c(cpoints.cors,sum(abs(seg$test-seg$testc)))
      }
    }
    
    if(skew>0.2){ # Cell is in early S-phase, binarize using 50 and 95 percen-tile of the segmented data
      
      c.rskew=c(c.rskew,sample)
      cat("\n","skew to right-early S")
      plotmix(quantile(test,c(0.5,0.95)))
      for(cpoint in cpoints ){
        seg[,"testb"]=seg[,"test"]-cpoint
        seg[seg$testb<0,"testc"]=quantile(test,0.5)
        seg[seg$testb>0,"testc"]=quantile(test,0.95)
        cpoints.cors=c(cpoints.cors,sum(abs(seg$test-seg$testc)))
      }
    }
    
    
    if(skew>=(-0.2) & skew<=0.2) {  #Cannot determine, binarize using 25 and 75 percentile of the segmented data. These cell may have noisy data, consider discarding later #
      c.noskew=c(c.noskew,sample)
      cat("\n","no skew")
      plotmix(quantile(test,c(0.25,0.75)))
      for(cpoint in cpoints ){
        seg[,"testb"]=seg[,"test"]-cpoint
        seg[seg$testb<0,"testc"]=quantile(test,0.25)
        seg[seg$testb>0,"testc"]=quantile(test,0.75)
        cpoints.cors=c(cpoints.cors,sum(abs(seg$test-seg$testc)))
      }
    }
    
  }
  
  
  point.2C=cpoints[which(cpoints.cors==min(cpoints.cors,na.rm=T))]   #Find the threshold at which the distance between binarized and un-binarized  is minimum
  point2Cs=c(point2Cs,point.2C[1])
}

#Binarize each single-cell using the specific ideal threshold calculated from above step and label with 1s and 0s#

for(i in 3:ncol(copyframe)){
  cat("\n",names(copyframe)[i]," ",point2Cs[i-2],"\n")
  copyframe[,i]=copyframe[,i]-point2Cs[i-2]
  copyframe[copyframe[,i]<0,i]=0
  copyframe[copyframe[,i]>0,i]=1
}

#Check binarization by plotting few cells using custom function in scRT_functions.R file

check.seg(RTs,copyframe,samples=c("B3","D2"),binary=T,chr=c(1,16),x1=2e6,x2=160e6)   

#16.	Rank cells in S-pahse

rank=data.frame(samples=names(copyframe[,3:ncol(copyframe)]))
rankfun=function(x) length(x[x==1 & is.finite(x)])/length(x[is.finite(x)])*100
rank$rank=apply(copyframe[,3:ncol(copyframe)],2,rankfun)
rank=rank[order(rank$rank),]

#17.	Final formatting of binarized single cell data

#Merge bulk replication timing data to single cell data
copyframe=merge(copyframe,RTbulk,by=c("CHR","POS"))
copyframe=copyframe[order(copyframe$CHR,copyframe$POS),]

#Add missing genomic bins to the binarized data#
copyframe$ID=as.character(paste(copyframe$CHR,"_",copyframe$POS,sep=""))

RT=read.delim("scRT_96cells_mm10_50kb.txt", header=T)
full=RT[,1:2]
names(full)=c("CHR","POS")
full=full[full$CHR!="chrM" & full$CHR!="chrY",]

full[,1]=as.character(full[,1])
full[full$CHR=="chrX",1]="chr20"
full$CHR=as.numeric(substring(full$CHR,4,6))
full$ID3=as.character(paste(full$CHR,"_",full$POS,sep=""))

L=!full$ID3 %in% copyframe$ID
copyframe$ID=NULL

omitted=full[L,c(1,2)]
add=cbind(omitted[,1:2], matrix(data = NA, nrow = nrow(omitted), ncol = ncol(copyframe)-2,))
names(add)[3:ncol(add)]=names(copyframe)[3:ncol(copyframe)]

copyframe=rbind(add,copyframe)
copyframe=copyframe[order(copyframe$CHR,copyframe$POS),]
#Visualizing the binarized data
#Binary heatmap ordered according to rank in S-phase#(Type "q" to quit)
for(chr in 1:20){
  
  RTb=subset(copyframe,CHR %in% c(chr))
  RTb.mat=RTb[,as.vector(rank$sample)]
  names(RTb.mat)=as.vector(rank$sample)
  p.binary=heatmap_vd2(RTb.mat,c(1:ncol(RTb.mat)),cols=c("snow2","red"),lim1=0.5, lim2=0.6,col_bias=1,binary=T,legend=T)
  print(p.binary)
  l=readline();if(l=="q") break   #Hit q and <ENTER> to exit loop#
}

#i.	Write binarized data
system("mkdir Output_mainscript")
write.table(copyframe,paste("./Output_mainscript/scRT_binary.txt",sep=""),row.names=F,sep="\t")

# ii.	Remove outlier binarized cells based on the hamming distance to other cells with similar S-phase ranking
wd=getwd()
setwd(paste0(wd,"/Output_mainscript"))

copyframe=read.delim("scRT_binary.txt", header=T)
bulkRT= copyframe$cast.129
copyframe$cast.129=NULL

#Calculate S-phase rank#
rank=data.frame(samples=names(copyframe[,3:ncol(copyframe)]))
rankfun=function(x) length(x[x==1 & is.finite(x)])/length(x[is.finite(x)])*100
rank$rank=apply(copyframe[,3:ncol(copyframe)],2,rankfun)
rank=rank[order(rank$rank),]

#Find outliers based one hamming distance between cells#
hammingdist=function(vec1,vec2) sum(abs(vec1-vec2),na.rm=T)


RTb.mat=copyframe[,as.vector(rank$sample)]

dist.mat=as.data.frame(matrix(NA,ncol(RTb.mat),ncol(RTb.mat)));row.names(dist.mat)=names(RTb.mat);names(dist.mat)=names(RTb.mat)
for(i in 1:ncol(RTb.mat)) {
  for(j in 1:ncol(RTb.mat))
  {
    dist.mat[i,j]=hammingdist(RTb.mat[,i],RTb.mat[,j])
  }
}

#plot hamming distance matrix
dist.mat.p=heatmap_vd2(dist.mat,c(1:ncol(dist.mat)),cols=c("red","orange","yellow"),col_bias=1,binary=F,legend=T)
print(dist.mat.p)

jpeg(file=paste("scRT_outlier_cells_remove.jpeg",sep=""),width=5000,height=5000,res=350)
print(dist.mat.p)
dev.off()

#Visually identify and  remove outliers###

remove=c("G12","D4","E11","G4","C10","H2","H4")
copyframe=copyframe[,-which(names(copyframe) %in% remove)]

#Write final binarized single cell data after outlier removal
copyframe$cast.129=bulkRT
write.table(copyframe,"scRT_binary_corrected.txt",row.names=F,sep="\t")
