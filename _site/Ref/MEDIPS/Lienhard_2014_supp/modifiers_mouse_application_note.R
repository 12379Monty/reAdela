#!/usr/bin/env Rscript

#load the required libraries
library(MEDIPS)
library(BSgenome.Mmusculus.UCSC.mm9)
packageVersion("MEDIPS")
#[1] '1.11.16'

#the input files (bam-alignments)
#this script assumes, the files are in a directory called "bam_dir"
bamfile_B6_Min_ad=c("bam_dir/ApcMin_Adenom_MeDIP10-2.bam",
  "bam_dir/ApcMin_Adenom_MeDIP2-1.bam",
  "bam_dir/ApcMin_Adenom_MeDIP4-1.bam",
  "bam_dir/ApcMin_Adenom_MeDIP4-2.bam",
  "bam_dir/ApcMin_Adenom_MeDIP9-3.bam"
)
bamfile_B6_Min_no=c("bam_dir/ApcMin_Normal_MeDIP4-3.bam",
  "bam_dir/ApcMin_Normal_MeDIP9-2.bam",
  "bam_dir/ApcMin_Normal_MeDIP9-4.bam")

bamfile_B6_wt=c("B6_MeDIP10-1.bam",
  "bam_dir/B6_MeDIP16-8.bam",
  "bam_dir/B6_MeDIP4-4.bam",
  "bam_dir/B6_MeDIP4-5.bam")

bamfile_B6_Min_no_input="bam_dir/ApcMin_Normal_Input4-3.bam"
bamfile_B6_Min_ad_input="bam_dir/ApcMin_Adenom_Input2-1.bam"
bamfile_B6_normal=c(bamfile_B6_wt,bamfile_B6_Min_no)
#parameters
BSgenome="BSgenome.Mmusculus.UCSC.mm9"
uniq=TRUE
extend=300
shift=0
ws=250
wd="/my/working/dir/"

############################################
## PREPROCESSING: reading the input files ##
############################################
date()
#all APC-Min sample files are processed
B6_Min_ad=lapply(X=bamfile_B6_Min_ad,FUN=MEDIPS.createSet,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,window_size = ws)
B6_Min_ad_input=MEDIPS.createSet(file = bamfile_B6_Min_ad_input,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,window_size = ws)
#the normal tissue sample files are pooled and processed
B6_normal=lapply(X= c(bamfile_B6_wt,bamfile_B6_Min_no),FUN=MEDIPS.createSet,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,window_size = ws)
B6_normal_input=MEDIPS.createSet(file=bamfile_B6_Min_no_input,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,window_size = ws)
#Calculation of the CpG density
CS = MEDIPS.couplingVector(pattern = "CG",refObj= B6_Min_ad[[1]])
date()
#time to read and preprocess the data: about 25 min on AMD Opteron 6380 using 1 CPU and ~6 GB RAM (peak)


######################
## Quality Control ##
######################
#calibration plot
for (i in 1:length(B6_Min_ad)) {
  png(paste(wd,"plots/SupplFig2_calibration_B6_Min_ad_",i,".png",sep=""))
  MEDIPS.plotCalibrationPlot(MSet=B6_Min_ad[[i]], ISet=B6_Min_ad_input, CSet=CS)
  dev.off()
}

for (i in 1:length(B6_normal)){
  png(paste(wd,"plots/SupplFig2_calibration_B6_normal_",i,".png",sep=""))
  MEDIPS.plotCalibrationPlot(MSet=B6_normal[[i]], ISet=B6_normal_input, CSet=CS)
  dev.off()
}

#small version
i=3
png(paste(wd,"plots/Fig1D_calibration_B6_Min_ad_",i,"_small.png",sep=""),400,400)
MEDIPS.plotCalibrationPlot(MSet=B6_Min_ad[[i]], ISet=B6_Min_ad_input, CSet=CS)
dev.off()


#saturation analysis
for (i in 1:length(B6_normal)) {
  png(paste(wd,"plots/SupplFig1_saturation_B6_normal_",i,".png",sep=""))
  sr=MEDIPS.saturation(file=bamfile_B6_normal[i], BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,window_size = ws, nit=10, nrit=1, empty_bins=TRUE, rank=FALSE)
  MEDIPS.plotSaturation(sr)
  dev.off()
}

for (i in 1:length(B6_Min_ad)) {
  png(paste(wd,"plots/SupplFig1_saturation_B6_Min_ad_",i,".png",sep=""))
  sr=MEDIPS.saturation(file=bamfile_B6_Min_ad[i], BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,window_size = ws, nit=10, nrit=1, empty_bins=TRUE, rank=FALSE)
  MEDIPS.plotSaturation(sr)
  dev.off()
}
i=3
png(paste(wd,"plots/Fig1C_saturation_B6_Min_ad_",i,"_small.png",sep=""),300,300)
sr=MEDIPS.saturation(file=bamfile_B6_Min_ad[i], BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,window_size = ws, nit=10, nrit=1, empty_bins=TRUE, rank=FALSE)
MEDIPS.plotSaturation(sr)
dev.off()


#######################################
## Differential Methylation analysis ##
#######################################
date()
#Create the result table and apply the statistical test for differential methylation
results_adenoma_vs_normal=MEDIPS.meth(MSet1 = B6_Min_ad, MSet2 = B6_normal,
  CSet = CS, p.adj = "bonferroni",  diff.method = "edgeR", 
  prob.method = "poisson", MeDIP = T,  CNV = F, type = "rpkm")
date()
#time for creating the result-table and the statistical test: about 15 min on AMD Opteron 6380 using 1 CPU and ~20 GB RAM
dm1=results_adenoma_vs_normal$edgeR.adj.p.value<0.1
dm2=results_adenoma_vs_normal$edgeR.p.value<0.01
write.table(results_adenoma_vs_normal[which(dm1),], paste(wd,"suppl_tab1_dmr.tsv",sep=""), sep="\t", quote=F, row.names=F, col.names=T)

#the result table contains, for each window 
#-the genomic locations, 
#-CpG density, 
#-read count per sample, 
#-normalized read count per sample 
#-mean values per group
#-logFC between groups
#-p-value and adjusted p-value

#an example is given in Suppl Tab1 (DMRs)

names(results_adenoma_vs_normal)
# [1] "chr"                                "start"                             
# [3] "stop"                               "CF"                                
# [5] "ApcMin_Adenom_MeDIP10.2.bam.counts" "ApcMin_Adenom_MeDIP2.1.bam.counts" 
# [7] "ApcMin_Adenom_MeDIP4.1.bam.counts"  "ApcMin_Adenom_MeDIP4.2.bam.counts" 
# [9] "ApcMin_Adenom_MeDIP9.3.bam.counts"  "B6_MeDIP10.1.bam.counts"           
#[11] "B6_MeDIP16.8.bam.counts"            "B6_MeDIP4.4.bam.counts"            
#[13] "B6_MeDIP4.5.bam.counts"             "ApcMin_Normal_MeDIP4.3.bam.counts" 
#[15] "ApcMin_Normal_MeDIP9.2.bam.counts"  "ApcMin_Normal_MeDIP9.4.bam.counts" 
#[17] "ApcMin_Adenom_MeDIP10.2.bam.rpkm"   "ApcMin_Adenom_MeDIP2.1.bam.rpkm"   
#[19] "ApcMin_Adenom_MeDIP4.1.bam.rpkm"    "ApcMin_Adenom_MeDIP4.2.bam.rpkm"   
#[21] "ApcMin_Adenom_MeDIP9.3.bam.rpkm"    "B6_MeDIP10.1.bam.rpkm"             
#[23] "B6_MeDIP16.8.bam.rpkm"              "B6_MeDIP4.4.bam.rpkm"              
#[25] "B6_MeDIP4.5.bam.rpkm"               "ApcMin_Normal_MeDIP4.3.bam.rpkm"   
#[27] "ApcMin_Normal_MeDIP9.2.bam.rpkm"    "ApcMin_Normal_MeDIP9.4.bam.rpkm"   
#[29] "ApcMin_Adenom_MeDIP10.2.bam.rms"    "ApcMin_Adenom_MeDIP2.1.bam.rms"    
#[31] "ApcMin_Adenom_MeDIP4.1.bam.rms"     "ApcMin_Adenom_MeDIP4.2.bam.rms"    
#[33] "ApcMin_Adenom_MeDIP9.3.bam.rms"     "B6_MeDIP10.1.bam.rms"              
#[35] "B6_MeDIP16.8.bam.rms"               "B6_MeDIP4.4.bam.rms"               
#[37] "B6_MeDIP4.5.bam.rms"                "ApcMin_Normal_MeDIP4.3.bam.rms"    
#[39] "ApcMin_Normal_MeDIP9.2.bam.rms"     "ApcMin_Normal_MeDIP9.4.bam.rms"    
#[41] "ApcMin_Adenom_MeDIP10.2.bam.prob"   "ApcMin_Adenom_MeDIP2.1.bam.prob"   
#[43] "ApcMin_Adenom_MeDIP4.1.bam.prob"    "ApcMin_Adenom_MeDIP4.2.bam.prob"   
#[45] "ApcMin_Adenom_MeDIP9.3.bam.prob"    "B6_MeDIP10.1.bam.prob"             
#[47] "B6_MeDIP16.8.bam.prob"              "B6_MeDIP4.4.bam.prob"              
#[49] "B6_MeDIP4.5.bam.prob"               "ApcMin_Normal_MeDIP4.3.bam.prob"   
#[51] "ApcMin_Normal_MeDIP9.2.bam.prob"    "ApcMin_Normal_MeDIP9.4.bam.prob"   
#[53] "MSets1.counts.mean"                 "MSets1.rpkm.mean"                  
#[55] "MSets1.rms.mean"                    "MSets1.prob.mean"                  
#[57] "MSets2.counts.mean"                 "MSets2.rpkm.mean"                  
#[59] "MSets2.rms.mean"                    "MSets2.prob.mean"                  
#[61] "edgeR.logFC"                        "edgeR.logCPM"                      
#[63] "edgeR.p.value"                      "edgeR.adj.p.value" 

#MA-plot
png(paste(wd,"plots/Fig1F_maplot_small.png",sep=""),400,400)
smoothScatter(results_adenoma_vs_normal$edgeR.logCPM,results_adenoma_vs_normal$edgeR.logFC,ylim=c(-4,4),xlim=c(-3,4), xlab="avg log methylation", ylab="methylation logFC", main="MA Plot")
abline(h=-4:4, col="grey",lty=2)
abline(v=c(-2,0,2,4), col="grey",lty=2)

points(results_adenoma_vs_normal$edgeR.logCPM[dm2],results_adenoma_vs_normal$edgeR.logFC[dm2], col="orange",pch=".")
points(results_adenoma_vs_normal$edgeR.logCPM[dm1],results_adenoma_vs_normal$edgeR.logFC[dm1], col="red",pch=4)

dev.off()

######################################
## Annotation and selection of ROIs ##
######################################

#here we use annotation tables as obtained from UCSC and other sources.
#for gtf files, only the first(chr), third(start) and 4th(end) column is of interest. 
#for an example of biomart database excess, please see the vignette. 
#read annotation tables
cc=c("character", "NULL", "NULL", "numeric", "numeric", "NULL", "NULL", "NULL", "NULL")
m_cgi=read.table("anno_dir/model-based-cpg-islands-mm9.txt", header=T)
m_exon=read.table("anno_dir/RefSeq_mm9_known_genes_chr_filter.gtf", sep="\t", colClasses=cc)
m_exon$nr=1:dim(m_exon)[1]
m_intron=read.table("anno_dir/RefSeq_mm9_known_genes_introns.txt", sep="\t", header=T)
m_prom=read.table("anno_dir/RefSeq_mm9_known_genes_promotor.gtf", sep="\t", colClasses=cc)
m_prom$nr=1:dim(m_prom)[1]

results_CGI=MEDIPS.selectROIs(results=results_adenoma_vs_normal, rois=m_cgi)
results_exon=MEDIPS.selectROIs(results=results_adenoma_vs_normal, rois=m_exon)
results_intron=MEDIPS.selectROIs(results=results_adenoma_vs_normal, rois=m_intron[,c(1:3,5)])
results_prom=MEDIPS.selectROIs(results=results_adenoma_vs_normal, rois=m_prom)
results_CGI_prom=MEDIPS.selectROIs(results=results_CGI, rois=m_prom)

head(results_CGI)


reg=c("all regions", "Introns", "Exons","Promoter" , "CpG Islands" , "CpG Island\nPromoter")

DMR_rel=matrix(NA, 6,2, dimnames=list(reg, c("hypomethylated", "hypermethylated")))
f=!is.na(results_adenoma_vs_normal$edgeR.p.value)
DMR_rel[1,1]=sum(results_adenoma_vs_normal$edgeR.p.value[f]<.01 & results_adenoma_vs_normal$edgeR.logFC[f]<0)/dim(results_adenoma_vs_normal)[1]
DMR_rel[1,2]=sum(results_adenoma_vs_normal$edgeR.p.value[f]<.01 & results_adenoma_vs_normal$edgeR.logFC[f]>0)/dim(results_adenoma_vs_normal)[1]

f=!is.na(results_intron$edgeR.p.value)
DMR_rel[2,1]=sum(results_intron$edgeR.p.value[f]<.01 &results_intron$edgeR.logFC[f]<0)/dim(results_intron)[1]
DMR_rel[2,2]=sum(results_intron$edgeR.p.value[f]<.01 &results_intron$edgeR.logFC[f]>0)/dim(results_intron)[1]

f=!is.na(results_exon$edgeR.p.value)
DMR_rel[3,1]=sum(results_exon$edgeR.p.value[f]<.01 &results_exon$edgeR.logFC[f]<0)/dim(results_exon)[1]
DMR_rel[3,2]=sum(results_exon$edgeR.p.value[f]<.01 &results_exon$edgeR.logFC[f]>0)/dim(results_exon)[1]

f=!is.na(results_prom$edgeR.p.value)
DMR_rel[4,1]=sum(results_prom$edgeR.p.value[f]<.01 &results_prom$edgeR.logFC[f]<0)/dim(results_prom)[1]
DMR_rel[4,2]=sum(results_prom$edgeR.p.value[f]<.01 &results_prom$edgeR.logFC[f]>0)/dim(results_prom)[1]

f=!is.na(results_CGI$edgeR.p.value)
DMR_rel[5,1]=sum(results_CGI$edgeR.p.value[f]<.01 &results_CGI$edgeR.logFC[f]<0)/dim(results_CGI)[1]
DMR_rel[5,2]=sum(results_CGI$edgeR.p.value[f]<.01 &results_CGI$edgeR.logFC[f]>0)/dim(results_CGI)[1]

f=!is.na(results_CGI_prom$edgeR.p.value)
DMR_rel[6,1]=sum(results_CGI_prom$edgeR.p.value[f]<.01 &results_CGI_prom$edgeR.logFC[f]<0)/dim(results_CGI_prom)[1]
DMR_rel[6,2]=sum(results_CGI_prom$edgeR.p.value[f]<.01 &results_CGI_prom$edgeR.logFC[f]>0)/dim(results_CGI_prom)[1]

png(paste(wd,"plots/Fig1G_mouse_dmr.png",sep=""), 600, 600)
barplot(t(DMR_rel)*100,  main="Fraction of differentially methylated regions", beside=T, ylab="Fraction of 250 pb windows [%]")
legend("topleft", c("hypomethylated", "hypermethylated"), fill=grey.colors(2))
dev.off()
#small version
png(paste(wd,"plots/Fig1G_mouse_dmr_small.png",sep=""), width=600, height=400)
par(cex=1.1)
barplot(t(DMR_rel)*100,  main="Fraction of differentially methylated regions", beside=T, ylab="Fraction of 250 pb windows [%]")
legend("topleft", c("hypomethylated", "hypermethylated"), fill=grey.colors(2))
dev.off()


##########################
## BiSulfite Validation ##
##########################

bsv=read.table("BisSeq_regions.tsv",sep="\t", header=T)
#the bisulfite validation table has the following format:
head(bsv)  
#associated_gene_name  dir   chr     start      stop length_bp BS_b3
#1                Cd133   ns  chr5  44493538  44493931       394    53
#2                Dusp6 hypo chr10  98728708  98729061       353    64
#3          Cela1, Ela1 hypo chr15 100516782 100517155       373    87
#4                 Eps8 hypo  chr6 137530962 137531369       407    66
#5               Fer1l4   ns  chr2 155863070 155863425       355    77
#6                  Fos hypo chr12  86817212  86817563       351    90
#  MeDIP_rms_b3 MeDIP2_rms_b3 BS_n5 MeDIP_rms_n5 MeDIP2_rms_n5 BS_ad5
#1         0.04          0.48    56         0.09     0.7015023     37
#2         0.33          1.59    59         0.14     0.5410639     19
#3         0.09          0.86    97         0.08     0.7175639     29
#4         0.12          0.75    59         0.08     0.4177759     24
#5         0.35          0.83    78         0.18     0.4406465     79
#6         0.38          1.39    95         0.25     0.9783356     54
#  MeDIP_rms_ad5 MeDIP2_rms_ad5
#1          0.06      0.5409102
#2          0.10      0.4085561
#3          0.04      0.4572876
#4          0.01      0.2246999
#5          0.34      0.7871628
#6          0.25      0.9602527

columns=c("CF","edgeR.logFC","edgeR.p.value","CNV.log2.ratio","ApcMin_Adenom_MeDIP9.3.bam.rms","ApcMin_Normal_MeDIP9.2.bam.rms", "B6_MeDIP10.1.bam.rms", "ApcMin_Adenom_MeDIP9.3.bam.rpkm","ApcMin_Normal_MeDIP9.2.bam.rpkm","B6_MeDIP10.1.bam.rpkm","ApcMin_Adenom_MeDIP9.3.bam.counts","ApcMin_Normal_MeDIP9.2.bam.counts", "B6_MeDIP10.1.bam.counts")
mr.roi=MEDIPS.selectROIs( results_adenoma_vs_normal , rois=bsv[,c(3:5,1)], columns=columns,summarize="avg")
#combine bisulfite and medip-seq tables
if(all(mr.roi$ROI==bsv[,1])){
mr.roi=cbind(mr.roi, bsv)
}

BS=rep(c("BS_ad5","BS_b3","BS_n5"),3)
sample=rep(c("APC Min Ad","B6","APC Min No"),3)
method=c(rep("rms",3),rep("rkpm",3),rep("count",3))

png(paste(wd,"plots/SupFig3_mouse_BS_validation.png",sep=""), height=700, width=1000)

par(mfrow=c(2,3))
color=rep("black", dim(mr.roi)[1])
color[mr.roi$CF>= quantile(mr.roi$CF,.75)]="red"
color[mr.roi$CF<= quantile(mr.roi$CF,.25)]="blue"
for(i in 1:6){
  val=mr.roi[,c(columns[i+4],BS[i]) ]
  val=val[!is.na(rowSums(val)),]
  co=cor(val[,1],val[,2])
  plot(val[,1],val[,2], sub=paste("correlation:", round(co,3)), main=sample[i], ylim=c(0,100),ylab=paste(sample[i], "BiSulfit [%]"), xlab=paste(sample[i], "MeDIP-Seq", method[i]), col=color)
  abline(lm(val[,2]~val[,1]))
}
legend("bottomright", legend=c("lower quantile CF", "25%-75% CF", "upper quantile CF"), fill=c("blue", "black", "red"))
dev.off()
#small version



png(paste(wd,"plots/Fig1E_mouse_BS_validation_small.png",sep=""), height=400, width=400)
i=1
par(cex=1.1)
  val=mr.roi[,c(columns[i+4],BS[i]) ]
  val=val[!is.na(rowSums(val)),]
  co=cor(val[,1],val[,2])
  plot(val[,1],val[,2], sub=paste("correlation:", round(co,3)), main=sample[i], ylim=c(0,100),ylab=paste(sample[i], "BiSulfit [%]"), xlab=paste(sample[i], "MeDIP-Seq", method[i]), col=color,pch=20)
  abline(lm(val[,2]~val[,1]))


legend("bottomright", legend=c("lower quantile CF", "25%-75% CF", "upper quantile CF"), fill=c("blue", "black", "red"))
dev.off()

