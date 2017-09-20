
# WORKING DIRECTORY
# ==============================================================
# Setting the enviroment
setwd("C:/Users/User/Desktop/OMICS DATA ANALISIS/09 TFM")


# LIBRARIES
# ==============================================================

# Required packages
library (FactoMineR)
library(factoextra)
library(factoextra)
library(TCGAbiolinks)
library(reshape2)
library(RTCGA.clinical)
library(survival)
library(survminer)
library(gplots)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggbio)
library(GenomicRanges)
library(plyr)
library(devtools)
library(genefilter)
library(pheatmap)
library(edgeR)
library(Biobase)
library(limma)
library(rafalib)
library(RColorBrewer)
library(rgl)
library(dendextend)

# DATA SOURCES
# ==============================================================
# Dataset of gene configurations 
setwd("~/omics/code to TFM/data/source data")
matrix<-read.table("matrix_res2_2", header=T, row.names=1, 
                   colClasses = "character", check.names=F)
configlegend<-read.table("leyend2import.csv",  header = F)

# Expression data:
setwd("~/omics/code to TFM/data/Expression data")
expression<-read.table("BLCA_rna.txt", header=TRUE, 
                       row.names=1, sep="\t", check.names = FALSE)
exp<-expression[c(2:20532),which(substr(colnames(expression), 14,16) == "01A")]
colnames(exp)<-substr(colnames(exp), 6,12)
colnames(exp)<-gsub("\\.", "-", colnames(exp))
exp.filt<-exp[, intersect(names(matrix), names(exp))]
exp.filt.names<-as.vector(sapply(rownames(exp.filt), function(x) 
  unlist(strsplit(x, split = "|", fixed = TRUE))[1]))


# Clinical data:
setwd("~/omics/code to TFM/data/Clinical data")
patients<-names(matrix)
clinical.TCGA.query<-GDCquery(project="TCGA-BLCA",
                              data.category = "Clinical",
                              barcode = patients)
GDCdownload(clinical.TCGA.query)
clinical_BLCA<-GDCprepare_clinic(clinical.TCGA.query, clinical.info ="patient")

load("clinical_BLCA.RData")
rownames(clinical_BLCA)<-patients
clinical.BLCA<-as.data.frame(clinical_BLCA)
clinical.BLCA$age<-round(-(clinical.BLCA$days_to_birth/365),0)
clinical.BLCA$agecat<-cut(clinical.BLCA$age, c(seq(0,100,10)))
clinical.BLCA.subinfo<-data.frame(clinical.BLCA$diagnosis_subtype, 
                                  clinical.BLCA$neoplasm_histologic_grade, 
                                  clinical.BLCA$stage_event_pathologic_stage, 
                                  clinical.BLCA$gender, clinical.BLCA$race_list, clinical.BLCA$age, clinical.BLCA$agecat)
rownames(clinical.BLCA.subinfo)<-rownames(clinical_BLCA)
clinical.BLCA.subinfo[clinical.BLCA.subinfo == ""]<-"NA"
clinical.BLCA.complete<-clinical.BLCA.subinfo[complete.cases(clinical.BLCA.subinfo),]
prop.table(table(clinical.BLCA.complete$clinical.BLCA.diagnosis_subtype))
prop.table(table(clinical.BLCA.complete$clinical.BLCA.stage_event_pathologic_stage))

#Mutation data # processed in bash sed 's/\b[12345678]\b/FALSE/gi' matrix_res2_2
#| sed 's/\b1[01234579]\b/FALSE/gi' | sed 's/\b2[0245]\b/FALSE/gi' | sed
#'s/\b9\b/TRUE/gi' | sed 's/\b1[68]\b/TRUE/gi' | sed 's/\b2[136789]\b/TRUE/gi' |
#sed 's/\b3[0123456789]\b/TRUE/gi' | sed 's/\b40\b/TRUE/gi' >
#matrix_res2_2_mutations
setwd("~/omics/code to TFM/data/source data")
mutations<-read.table("matrix_res2_2_mutations", header=T, row.names=1, 
                      colClasses = "character", check.names = F)


#CNV data #processed in bash sed
#'s/\b[12]\b\|\b13\b\|\b2[48]\b\|\b3[089]\b/loss_2/gi' matrix_res2_2 | sed
#'s/\b[6]\b\|\b1[146]\b\|\b2[267]\b\|\b35\b/loss_1/gi' | sed
#'s/\b[578]\b\|\b2[019]\b\|\b3[36]\b/NULL/gi' | sed
#'s/\b[349]\b\|\b1[258]\b\|\b23\b\|\b31\b/gain_1/gi' | sed
#'s/\b1[079]\b\|\b25\b\|\b3[247]\b\|\b40\b/gain_2/gi' > matrix_res2_2_CNV
CNV<-read.table("matrix_res2_2_CNV", header=T, colClasses = "character", 
                row.names=1, check.names = F)


#Methylation data #processed in bash sed
#'s/\b[1368]\b\|\b1[068]\b\|\b2[18]\b\|\b34\b/low/gi' matrix_res2_2 | sed
#'s/\b[2459]\b\|\b1[17]\b\|\b2[79]\b\|\b3[02]\b/medium/gi' | sed
#'s/\b7\b\|\b1[3459]\b\|\b2[36]\b\|\b3[378]\b/high/gi' | sed
#'s/\b12\b\|\b2[0245]\b\|\b3[1569]\b\|\b40\b/NA/gi' > matrix_res2_2_methylation
Methy<-read.table("matrix_res2_2_methylation", header=T, colClasses = "character", 
                  row.names=1, check.names = F)


# 1. CASE SELECTION
# ==============================================================
# Filter by race: caucasian
caucasian.patients<-rownames(clinical.BLCA[clinical.BLCA$race_list == "WHITE",])
# Filter by gender: males
male.patients<-rownames(clinical.BLCA[clinical.BLCA$gender == "MALE",])
#Select only caucasian males
caucasian.males<-intersect(caucasian.patients, male.patients)
#Select only complete cases for clinical info. 
caucasian.males.complete.cases<-intersect(caucasian.males, 
                                          rownames(clinical.BLCA.complete))
matrix_caucasian<-matrix[,names(matrix) %in% caucasian.males.complete.cases]


# 2. DATA FILTERING
# ==============================================================
# 2A. FILTER GENES WITH NA VALUES IN MORE THAN 10% OF THE SAMPLES
filter_genes_NA <- c ("PLXND1","TACC3","CLCA4","C8B","SLC4A7","PHF23","CP","PUM2","SPA17","FNDC8","IGF2BP2","TRAF4","AFM","CFHR2","HSPB11","ANKRD24","CD209","DTX2","DGCR14","GGTLC2","EIF3L","DHRS2","AMMECR1","MSLN","NME4","SGK3","CSPP1","C19orf53","TJP3","ARMC6","FKBP14","CYP3A5","SPAM1","HSD17B1","CRTAM","SLCO1B3","AICDA","MRPL2","IMPG1","C9","TFDP2","LRRC31","PPM1G","OPRD1","FBXO6","CFHR3","MPL","PHF3","RAB3GAP2","HOXB5","INSL4","CSMD2","COPA","TAF1L","GIPC1","TP53TG5","HCST","CANX","OR7A10","GFER","OR1E2","POR","FGF13","ACSS2","TPTE2","CFHR4","HRH4","UGT2A3","AKIRIN2","FAM186B","TAF5L","SETDB2","CBWD2","USP20","IL33","SLC3A1","CYP2C8","CCDC54","MTTP","SBNO1","PARN","GAS8","OSBPL1A","RPS8","GPA33","F13B","FCRL5","NVL","SCN1A","TIMD4","TAAR8","NONO","SEC16A","MSRB2","SPINT4","LPHN3","CCDC74B","CYSLTR2","UTRN","ZNF117","CPB1","TMPRSS11D","PTDSS1","FAM160B2","WFDC8","FBXL13","CYP4A22","CCDC74A","WDR82","IL31RA","CSMD3","PRKACG","C8orf34","C9orf84","GTF2A1","ALKBH3","SPIC","STXBP4","OR2D2","OR4D6","URM1","TSR1","KLK4","OR13J1","TAS2R1","PYDC1","MAP3K2","SERPINA9","DLGAP1","STAT2","OSCAR","OR7G2","OR4D5","GRIK1","PROL1","UGT2B7","OR2AT4","PDZD3","OR5T3","AGAP5","OR1S1","CBWD1","XKR3","TADA2B","OR10K1","C9orf131","MGA","OR5R1","OR10AG1","OR11H6","OR4K17","OR10H4","DMRTA1","ANKLE2","OR2T12","MAF1","CRIPAK","OR1E1","WFDC11","C12orf36","DEFB112","OR56A1","OR52E4","OR52D1","OR4C6","OR4C16","FAM153B","BACE2","OTOL1","OR10G7","IGSF5","OR56A4","PMCH","FAM167B","CLECL1","C6orf58","PIWIL3","ODF4","OR6C70","ANKS1B","ZNF479","ZBTB37","SCN10A","KCNIP4","C1orf110","CALHM1","NKAIN3","OR5D14","OR10X1","KPNA4","ZNF197","GPR141","TMPRSS11A","OR6C75","OR51F1","ZNF383","CENPP","OR6N2","SNTN","LIPI","GJB5","OR6M1","VN1R2","OR5AC2","PATE2","KRT39","CXorf40B","OR6N1","SLC28A3","DMBX1","MAP1LC3C","OR1S2","OR5H2","TMEM207","OR2M2","AKR1C4","OXCT2","SUCNR1","TRIM10","ZFP57","OR2W1","C7orf66","CRYZL1","SERPINB4","CD200R1L","GBP7","KLRK1","ZNF99","OR6B1","TAS2R41","MAGEA3","OR2F2","OR2A2","TBC1D3B","ORM2","VN1R4","DYTN","TAS2R39","OR4A47","OR7E24","STRC","CFHR1","OR2AE1","N4BP2L2","PRODH2","PYDC2","DPP3","OR5M1","TAS2R43","ABCA10")

# 2B. FILTER MCA OUTLIERS
# Filering MCA outliers in data
matrix_clean<-matrix_caucasian[! rownames(matrix_caucasian) %in% filter_genes_NA,]
matrix_t_clean <- t(matrix_clean)
result <- MCA(matrix_t_clean, graph = FALSE)
ind.mca.outliers<-names(c(which(result$ind$coord[,2]>1),
                          which(result$ind$coord[,1]>1)))
# Filering MCA outliers in SNV data
mutations_t<-t(mutations)
mutations_t<-mutations_t[row.names(matrix_t_clean),]
mutations_t<-mutations_t[,! colnames(mutations_t) %in% filter_genes_NA]
mca_mutations<-MCA(mutations_t, graph = FALSE)
ind.mut.outliers<-names(c(which(mca_mutations$ind$coord[,1]>1) , 
                          which(mca_mutations$ind$coord[,2]>1)))
# Filering MCA outliers in CNV data
CNV_t<-t(CNV)
CNV_t<-CNV_t[row.names(matrix_t_clean),]
CNV_t<-CNV_t[,! colnames(CNV_t) %in% filter_genes_NA]
mca_CNV<-MCA(CNV_t, graph = FALSE)
ind.cnv.outliers<-names(which(mca_CNV$ind$coord[,2]>2))
# Filering MCA outliers in METHYLATION data
Methy_t<-t(Methy)
Methy_t<-Methy_t[row.names(matrix_t_clean),]
Methy_t<-Methy_t[,! colnames(Methy_t) %in% filter_genes_NA]
mca_Methy<-MCA(Methy_t, graph = FALSE)
ind.met.outliers<-names(c(which(mca_Methy$ind$coord[,1]>5), 
                          which(mca_Methy$ind$coord[,2]>(5))))
# summary for outliers criteria
filter_cases <- unique(c(ind.mca.outliers, ind.mut.outliers, 
                         ind.cnv.outliers, ind.met.outliers)) 

# 2C. FILTER GENES WITH LOW AND HIGH VARIABILITY OF CONFIGURATIONS 
## Count Matrix:
count_config<-apply(matrix_t_clean, 2, table)
matrix_config_count<-matrix(nrow = ncol(matrix_t_clean), ncol = 40)
rownames(matrix_config_count)<-names(count_config)
colnames(matrix_config_count)<-c(1:40)
l<-length(count_config)
library(rgl)
for (i in 1:l){
  index.list<- melt(count_config[i])[,1]
  value.list<-melt(count_config[i])[,2]
  names(value.list)<-index.list
  total_count<-sum(value.list)
  for (q in 1:length(value.list)){
    matrix_config_count[i,index.list[q]]<-
      ((value.list[as.character(index.list[q])])/total_count)*100
  }
}
# Remove NA and fill with zeros:
matrix_config_count[is.na(matrix_config_count)] <- 0
#Measure of covariation and filtering
CV<- function(x){(mean(x)/sd(x))*100}
CV.genes<-apply(matrix_config_count, 1, CV)
# Quantiles
(cbind(cv.quantiles<-quantile(CV.genes, probs = seq(0,1, 0.01))))
# Cut-off genes with high and low coeff of variation in configurations:
selected.low<-CV.genes[CV.genes<cv.quantiles ["25%"]]
selected.high<-CV.genes[CV.genes>cv.quantiles["95%"]]
filtered_variability<-c(selected.low,selected.high)
filtered_genes_variability<-names(filtered_variability)


# 2D. SECOND MCA FILTERING TO REMOVE DIMENSIONAL OUTLIERS
matrix_t_filtered<-matrix_t_clean[! row.names(matrix_t_clean) 
                                  %in% filter_cases,]
matrix_t_filtered<-matrix_t_filtered[,! colnames(matrix_t_filtered)
                                     %in% filtered_genes_variability]
mca1.filtered<-MCA(matrix_t_filtered, graph = FALSE)
mca1.outliers<-names(c(which(mca1.filtered$ind$coord[,2]>1.5), 
                       which(mca1.filtered$ind$coord[,3]>1)))
# Summary of filter cases after second MCA filtering
filter_cases2 <- unique(c(ind.mca.outliers, ind.mut.outliers, 
                          ind.cnv.outliers, ind.met.outliers, mca1.outliers)) 


plot.MCA(result, invisible=c("var"), cex=0.7, label='ind')
abline(v=5, h=5, col="red", lwd=3, lty=2)
# Plot MCA outliers
par(mfrow=c(3,2))
plot.MCA(mca_mutations, invisible=c("var"), cex=0.7, label='ind')
abline(v=1, h=1, col="red", lwd=3, lty=2)
plot.MCA(mca_CNV, invisible=c("var"), cex=0.7, label='ind')
abline(h=2, col="red", lwd=3, lty=2)
plot.MCA(mca_Methy, invisible=c("var"), cex=0.7, label='ind')
abline(v=5, h=5, col="red", lwd=3, lty=2)
plot.MCA(mca1.filtered, invisible=c("var"), cex=0.7, label='ind')
abline(v=1.5, h=1.5, col="red", lwd=3, lty=2)
library(rgl)
open3d()
plot3d(mca1.filtered$ind$coord[,1:3],  col = "blue", size = 5, box=FALSE) 
points3d(mca1.filtered$ind$coord["GV-A3QI",1],mca1.filtered$ind$coord["GV-A3QI",2],
         mca1.filtered$ind$coord["GV-A3QI",3], col="red")


# 3. MULTICORESPONDENCE ANALYSIS ON FINAL FILTERED DATASET.
# ==============================================================
matrix_t_clean_2<-matrix_t_filtered[! row.names(matrix_t_filtered)
                                    %in% filter_cases2, ]
result_filter_var <- MCA(matrix_t_clean_2, graph = FALSE)
matrix_clean_2<-t(matrix_t_clean_2)
# 3A. Hieralchical Clustering on Principal components
cluster_var<-HCPC(result_filter_var)
plot.HCPC(cluster_var, choice="map", angle=90, ind.names=FALSE, draw.tree=FALSE )
plot( cluster_var, choice ="tree", cex = 0.6)
plot(cluster_var, choice = "3D.map")
## Individuals in each cluster
contribucion<-result_filter_var$ind$contrib
ind_in_cluster<-as.matrix(cluster_var$call$X$clust)
rownames(ind_in_cluster)<-rownames(cluster_var$call$X)
table(ind_in_cluster)

#3D MCA plotlibrary(rgl) for individuals
plot3d(result_filter_var$ind$coord[,1:3],  
       col = as.factor(ind_in_cluster[rownames(result_filter_var$ind$coord),]), 
       size = 2, type = "s")




#4. CLUSTER VALIDATION
# ==============================================================
# http://www.sigmath.es.osaka-u.ac.jp/shimo-lab/prog/pvclust/
library(pvclust)
#Pvclust function, bootstrap the data between clusters
pvclust.result <- pvclust(t(result_filter_var$ind$coord), 
                          method.dist="euclidean", method.hclust="ward.D", nboot=100, r=1)
plot(pvclust.result)
pvrect(pvclust.result, alpha=0.95)
pv.culsters<-pvpick(pvclust.result)
#pvclust provides two types of p-values: AU (Approximately Unbiased) p-value and
#BP (Bootstrap Probability) value. AU p-value, which is computed by multiscale
#bootstrap resampling, is a better approximation to unbiased p-value than BP
#value computed by normal bootstrap resampling.


# http://www.win-vector.com/blog/2015/09/bootstrap-evaluation-of-clusters/
# Clusterboot() to individuals
library(fpc)
x= result_filter_var$ind$coord[,1]
y= result_filter_var$ind$coord[,2]
z= result_filter_var$ind$coord[,3]
cboot.hclust <- clusterboot(cbind(x,y),clustermethod=hclustCBI,
                            method="ward.D", k=3)
groups<-cboot.hclust$result$partition  
table(groups)
cboot.hclust$bootmean
cboot.hclust$bootbrd 

xtable((rbind(table(groups), round(cboot.hclust$bootmean,2), round(cboot.hclust$bootbrd,2)))[2:3,])

# with all the dimensions is more inestable
#cboot.hclust <-clusterboot(result_filter_var$ind$coord,clustermethod=hclustCBI,
#                           method="ward.D", k=3)
#cboot.hclust$bootmean# this value should be neart to 1 to be strong cluster
#cboot.hclust$bootbrd # this is the number of times that the cluster is disolved


# Hopkins statistic
library(clustertend)
get_clust_tendency(result_filter_var$ind$coord, n=50)
#$hopkins_stat  [1] 0.2324831 if the value of Hopkis statistic is close to zero
#the we can reject the null hipotesis and conclude that the dataset is
#significantly a clusterable data.


# Validation by dissplot
library(seriation)
km.res<-kmeans(result_filter_var$ind$coord,3)
data_dist<- dist(result_filter_var$ind$coord)
dissplot(data_dist, labels=km.res$cluster)


# 5. CHARACTERIZING CLUSTERS
# ==============================================================
#Cluster1
ind_in_cluster_1<-as.data.frame(ind_in_cluster[which(ind_in_cluster =="1"),])
colnames(ind_in_cluster_1)<-"cluster"
index_c1<-ind_in_cluster_1

#Cluster2
ind_in_cluster_2<-as.data.frame(ind_in_cluster[which(ind_in_cluster =="2"),])
colnames(ind_in_cluster_2)<-"cluster"
index_c2<-ind_in_cluster_2

#Cluster3
ind_in_cluster_3<-as.data.frame(ind_in_cluster[which(ind_in_cluster =="3"),])
colnames(ind_in_cluster_3)<-"cluster"
index_c3<-ind_in_cluster_3


data2index<-rbind(index_c1, index_c2, index_c3)
index<-rownames(data2index)
#Analisis of configurations distribution for each gene and cluster differencies

# Total matrix of configurations
indexed_matrix<-matrix_t_clean_2[index,]
identical(rownames(indexed_matrix), rownames(data2index))
matrix_clusters<-cbind(data2index, indexed_matrix)
#Contingency table for configuration frecuencies:
xtabs(~ matrix_clusters$FGFR3 + matrix_clusters$cluster, data = matrix_clusters)
mosaicplot(t(xtabs(~ matrix_clusters$FGFR3 + matrix_clusters$cluster, data = matrix_clusters)), shade = TRUE, type = "FT")
# plot just a subset of the table
library("vcd")
assoc(head(xtabs(~ matrix_clusters$FGFR3 + matrix_clusters$cluster, data = matrix_clusters)), shade = T, las=3)
chisq.test(xtabs(~ matrix_clusters$FGFR3 + matrix_clusters$cluster, data = matrix_clusters))


# 6. CORRELATION WITH THE CLINICAL DATA
# ==============================================================

clinical.sel<-clinical.BLCA.complete[index,]

head(clinical.sel)
identical(rownames(indexed_matrix), rownames(clinical.sel))
cluster<-data2index[rownames(clinical.sel),]$cluster
data4mfa<-cbind(cluster, clinical.sel)
library(gmodels)
dim(data4mfa)
head(data4mfa)

# AGECAT INDEPENDENT p = Chi^2 =  3.302093     d.f. =  2     p =  0.191849    !!!!!
clinical.BLCA$agecat<-factor(clinical.BLCA$agecat)
struct <- structable(~ cluster + clinical.BLCA.agecat, data = data4mfa)
mosaic(struct, data = data4mfa, shade = TRUE, direction = "v",pop=F, main = "Age Cat")
labeling_cells(text = as.table(round(struct/nrow(data4mfa),2))*100, margin = 0)(as.table(round(struct/nrow(data4mfa),2)))
barplot(round(struct/nrow(data4mfa),2)*100,beside = T, col = terrain.colors(3))
legend("topleft", levels(cluster), pch = 15, col = terrain.colors(3))
CrossTable(data4mfa$cluster,data4mfa$clinical.BLCA.diagnosis_subtype, prop.r=F, prop.c=F, prop.t=F, prop.chisq=T, chisq = T)



# CLINICAL SUBTYPE INDEPENDENT p =  0.191849    !!!!!
data4mfa$clinical.BLCA.diagnosis_subtype<-factor(data4mfa$clinical.BLCA.diagnosis_subtype)
struct <- structable(~ cluster + clinical.BLCA.diagnosis_subtype, data = data4mfa)
mosaic(struct, data = data4mfa, shade = TRUE, direction = "v",pop=F, main = "Histological subtype")
labeling_cells(text = as.table(round(struct/nrow(data4mfa),2)), margin = 0)(as.table(round(struct/nrow(data4mfa),2)))
barplot(round(struct/nrow(data4mfa),2)*100,beside = T, col = terrain.colors(3))
legend("topleft", levels(cluster), pch = 15, col = terrain.colors(3))
CrossTable(data4mfa$cluster,data4mfa$clinical.BLCA.diagnosis_subtype, prop.r=F, prop.c=F, prop.t=F, prop.chisq=T, chisq = T)



# HISTOLOGIC GRADE DEPENDENT p =  0.2866076 
data4mfa$clinical.BLCA.neoplasm_histologic_grade<-factor(data4mfa$clinical.BLCA.neoplasm_histologic_grade)
struct <- structable(~ cluster + clinical.BLCA.neoplasm_histologic_grade, data = data4mfa)
mosaic(struct, data = data4mfa, shade = TRUE, direction = "v",pop=F, main = "Grade")
labeling_cells(text = as.table(round(struct/nrow(data4mfa),2)), margin = 0)(as.table(round(struct/nrow(data4mfa),2)))
barplot(round(struct/nrow(data4mfa),2)*100,beside = T, col = terrain.colors(3))
legend("topleft", levels(cluster), pch = 15, col = terrain.colors(3))
CrossTable(data4mfa$cluster,data4mfa$clinical.BLCA.neoplasm_histologic_grade, prop.r=F, prop.c=F,
           prop.t=F, prop.chisq=T, chisq = T)


# Patological stage independent  p =  0.599314
data4mfa$clinical.BLCA.stage_event_pathologic_stage<-factor(data4mfa$clinical.BLCA.stage_event_pathologic_stage)
struct <- structable(~ cluster + clinical.BLCA.stage_event_pathologic_stage, data = data4mfa)
mosaic(struct, data = data4mfa, shade = TRUE, direction = "v",pop=F, main = "Grade")
labeling_cells(text = as.table(round(struct/nrow(data4mfa),2)), margin = 0)(as.table(round(struct/nrow(data4mfa),2)))
barplot(round(struct/nrow(data4mfa),2)*100,beside = T, col = terrain.colors(3))
legend("topleft", levels(cluster), pch = 15, col = terrain.colors(3))
CrossTable(data4mfa$cluster,data4mfa$clinical.BLCA.stage_event_pathologic_stage, prop.r=F, prop.c=F,
           prop.t=F, prop.chisq=T, chisq = T)


#7. COUNT MATRIX TO EVALUATE CONFIGURATION DISTRIBUION AMONG CLUSTERS
# ==============================================================

## Count Matrix:
count_config_ind<-apply(matrix_t_filtered, 1, table)
matrix_config_ind_count<-matrix(nrow = nrow(matrix_t_filtered), ncol = 40)
rownames(matrix_config_ind_count)<-names(count_config_ind)
colnames(matrix_config_ind_count)<-c(1:40)
l<-length(count_config_ind)
for (i in 1:l){
  index.list<- melt(count_config_ind[i])[,1]
  value.list<-melt(count_config_ind[i])[,2]
  names(value.list)<-index.list
  for (q in 1:length(value.list)){
    matrix_config_ind_count[i,index.list[q]]<-value.list[as.character(index.list[q])]
  }
}
matrix_config_ind_count[is.na(matrix_config_ind_count)] <- 0
#Per cluster
countdata<-matrix_config_ind_count[rownames(ind_in_cluster),]
count.sum<-apply(countdata, 2, sum)

c1_count<-apply(countdata[ind_in_cluster == "1",1:40], 2, sum)
c2_count<-apply(countdata[ind_in_cluster == "2",1:40], 2, sum)
c3_count<-apply(countdata[ind_in_cluster == "3",1:40], 2, sum)

count_config_per_cluster<-cbind(c1_count, c2_count, c3_count)
colnames(count_config_per_cluster)<-c("cluster1","cluster2","cluster3")
rownames(count_config_per_cluster)<-c(1:40)


total_count_per_cluster<-apply(count_config_per_cluster,2,sum)
total_count_per_config<-apply(count_config_per_cluster, 1, sum)



freq_config_per_cluster<-cbind(c1_count/count.sum,c2_count/count.sum,c3_count/count.sum )
freq_config_per_cluster[is.na(freq_config_per_cluster)] <- 0.0000000001 # to avoid marascuillo error
freq_config_per_cluster<-round(freq_config_per_cluster,2)
colnames(freq_config_per_cluster)<-c("cluster1","cluster2","cluster3")
rownames(freq_config_per_cluster)<-c(1:40)
freq_config_per_cluster

# Test for proportions, Marascuillo procedure
## Set the proportions of interest.

## Compute critical values.
maras_procedure<-function(p,s){
  N = length(p)
  value = critical.range = tag = c()
  categories <- c("C1","C2", "C3")
  significntive<-c()
  for (i in 1:(N-1)){ 
    for (j in (i+1):N){
      value <- round(c(value,(abs(p[i]-p[j]))),2)
      critical.range = round(c(critical.range,
                               sqrt(qchisq(.95,N-1))*sqrt(p[i]*(1-p[i])/s + p[j]*(1-p[j])/s)),2)
      if (value > critical.range){
        significative<-"Yes"
      }else{
        significative<-"No"
      }
      tag = c(tag, paste(categories[i], categories[j], sep = " to "))
    }
  }
  df <- as.matrix(cbind(tag,critical.range,value, significative))
  print(df)
}



for (i in c(1:40)){
  cat("configuration ",i)
  print(i)
  maras_procedure(as.vector(as.matrix(freq_config_per_cluster[i,])),total_count_per_config[i]) 
}

#7A. PLOTS IN COUNTS AND FRECUENCIES OF CLUSTERS
# ==============================================================

#barplots
par(mfrow=c(1,3))
freq_config_per_cluster<-as.data.frame(freq_config_per_cluster)
barplot(freq_config_per_cluster$cluster1, col= "darkgreen")
barplot(freq_config_per_cluster$cluster2, col= "beige")
barplot(freq_config_per_cluster$cluster3, col= "gray")
# barplots beside
par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
barplot(t(as.matrix(freq_config_per_cluster)), col=terrain.colors(3), width=2, beside=TRUE)
legend("topright",inset=c(-0.25,0), fill=terrain.colors(3), legend=colnames(freq_config_per_cluster))
# boxplots
par(mfrow=c(3,1))
boxplot(countdata[ind_in_cluster == "1",])
boxplot(countdata[ind_in_cluster == "2",])
boxplot(countdata[ind_in_cluster == "3",])
# boxplots whit ggplot
library(reshape)
library(ggplot2)
c1.data<-melt(countdata[ind_in_cluster == "1",], id.var = "Config")
c1.data$cluster<-rep("c1", length(nrow(c1.data)))
c2.data<-melt(countdata[ind_in_cluster == "2",], id.var = "Config")
c2.data$cluster<-rep("c2", length(nrow(c2.data)))
c3.data<-melt(countdata[ind_in_cluster == "3",], id.var = "Config")
c3.data$cluster<-rep("c3", length(nrow(c3.data)))
df.m<-rbind(c1.data, c2.data, c3.data)
df.m$Var2<-as.character(df.m$Var2)
ggplot(data = df.m, aes(x=Var2, y=value)) + geom_boxplot(aes(fill=df.m$cluster))



# 8. PLOTS COLLAPSING THE GENOMIC DATA
# ==============================================================

# MOSAIC ########

# Mutations
c1.mut.matrix<-mutations_t[rownames(ind_in_cluster_1),]
cluster_1s=table(c1.mut.matrix)
c2.mut.matrix<-mutations_t[rownames(ind_in_cluster_2),]
cluster_2s=table(c2.mut.matrix)
c3.mut.matrix<-mutations_t[rownames(ind_in_cluster_3),]
cluster_3s=table(c3.mut.matrix)
mosaic(rbind(cluster_1s,cluster_2s,cluster_3s), shade = T)
total_s<-sum(cluster_1s, cluster_2s, cluster_3s)
round(100*(cluster_1s/total_s),2)
round(100*(cluster_2s/total_s),2)
round(100*(cluster_3s/total_s),2)

# CNVS
c1.cnv.matrix<-CNV_t[rownames(ind_in_cluster_1),]
cluster_1c<-table(c1.cnv.matrix)
c2.cnv.matrix<-CNV_t[rownames(ind_in_cluster_2),]
cluster_2c<-table(c2.cnv.matrix)
c3.cnv.matrix<-CNV_t[rownames(ind_in_cluster_3),]
cluster_3c<-table(c3.cnv.matrix)
mosaic(rbind(cluster_1c,cluster_2c,cluster_3c), shade = T)
total_c<-sum(cluster_1c, cluster_2c, cluster_3c)
round(100*(cluster_1c/total_c),2)
round(100*(cluster_2c/total_c),2)
round(100*(cluster_3c/total_c),2)

# Methylation
c1.met.matrix<-Methy_t[rownames(ind_in_cluster_1),]
cluster_1m<-table(c1.met.matrix)
c2.met.matrix<-Methy_t[rownames(ind_in_cluster_2),]
cluster_2m<-table(c2.met.matrix)
c3.met.matrix<-Methy_t[rownames(ind_in_cluster_3),]
cluster_3m<-table(c3.met.matrix)
mosaic(rbind(cluster_1m,cluster_2m,cluster_3_m), shade = T)
total_m<-sum(cluster_1m, cluster_2m, cluster_3m)
round(100*(cluster_1m/total_m),2)
round(100*(cluster_2m/total_m),2)
round(100*(cluster_3m/total_m),2)

# GGPLOT  ########

library(vcd)
library(scales)
# Normal
# as: MUT FALSE, then: low, medium and high methylathed and inside as -2, -1, 0, 1, 2 cnv:
list<-list(c("1","6","8","3","10"), c("2","11","5","4","17"), c("13","14","7","15","19"))
setwd("~/omics/code to TFM/data/Data saved objects")
for (i in 1:length(list)){
  df.1<-melt(countdata[ind_in_cluster == "1", unlist(list[i])],id.var="config")
  df.2<-melt(countdata[ind_in_cluster == "2", unlist(list[i])],id.var="config")
  df.3<- melt(countdata[ind_in_cluster == "3", unlist(list[i])],id.var="config")
  df.1$cluster<-rep("c1", length(nrow(df.1)))
  df.2$cluster<-rep("c2", length(nrow(df.2)))
  df.3$cluster<-rep("c3", length(nrow(df.3)))
  df.m<-rbind(df.1, df.2, df.3)
  df.m$Var2<-as.character(df.m$Var2)
  df.m$Var2 <- factor(df.m$Var2, levels=unique(df.m$Var2))
  p<-ggplot(data = df.m, aes(x=Var2, y=value, fill=df.m$cluster)) +
    geom_bar(position = "fill", stat = "identity") + 
    scale_y_continuous(labels= percent_format())
  #+ geom_boxplot(aes(fill=df.m$cluster)
  outname<-paste(i, ".pdf")
  ggsave(outname, plot = p)
}


# Mutated
# as: MUT TRUE, then: low, medium and high methylathed and inside as -2, -1, 0, 1, 2 cnv:
list<-list(c("28","16","21","18","34"), c("30","27","29","9","32"), c("38","26","33","23","37"))

for (i in 1:length(list)){
  df.1<-melt(countdata[ind_in_cluster == "1", unlist(list[i])],id.var="config")
  df.2<-melt(countdata[ind_in_cluster == "2", unlist(list[i])],id.var="config")
  df.3<- melt(countdata[ind_in_cluster == "3", unlist(list[i])],id.var="config")
  df.1$cluster<-rep("c1", length(nrow(df.1)))
  df.2$cluster<-rep("c2", length(nrow(df.2)))
  df.3$cluster<-rep("c3", length(nrow(df.3)))
  df.m<-rbind(df.1, df.2, df.3)
  df.m$Var2<-as.character(df.m$Var2)
  df.m$Var2 <- factor(df.m$Var2, levels=unique(df.m$Var2))
  p<-ggplot(data = df.m, aes(x=Var2, y=value, fill=df.m$cluster)) + 
    geom_bar(position = "fill", stat = "identity")  +
    scale_y_continuous(labels= percent_format())
  #+ geom_boxplot(aes(fill=df.m$cluster)
  outname<-paste(i, "M.pdf")
  ggsave(outname, plot = p)
}



# 8. COMPARISION WITH TCGA DATA
# ==============================================================
setwd("~/omics/code to TFM/data/TCGAdata")
tcga_blca<-read.table("BLCA_TCGA.csv", sep=" ")
colnames(tcga_blca)<-c("barcode", "cluster")
tcga_blca$barcode<-substr(tcga_blca$barcode, 6,12)
rownames(tcga_blca)<-tcga_blca$barcode
tcgaANDme<-intersect(rownames(tcga_blca), rownames(data2index))
mydata_tcgadata<-data.frame(tcga_blca[tcgaANDme,],data2index[tcgaANDme,])[,2:3]
table(mydata_tcgadata)


# 9. SURVIVAL DATA AND ANALYSIS
# ==============================================================
setwd("~/omics/code to TFM/data/Clinical data")
biocLite("RTCGA.clinical")
library(RTCGA.clinical)
library(survival)
library(survminer)

surv.clin<-survivalTCGA(BLCA.clinical, extract.cols= "admin.disease_code")
surv.clin<-load("~/omics/code to TFM/data/Clinical data/surv.clin.RData")
surv.clin$bcr_patient_barcode<-substr(surv.clin$bcr_patient_barcode, 6,12)
rownames(surv.clin)<-surv.clin$bcr_patient_barcode
index<-intersect(rownames(surv.clin), rownames(data2index))
surv.clin.sel<-surv.clin[index,]
dim(surv.clin.sel)
data2index[index,]

surv.data<-data.frame(surv.clin.sel[,c(1,3)], data2index[index,])
names(surv.data)<-c("times", "status", "cluster")
head(surv.data)
prop.table(table(surv.data$cluster, surv.data$status),1)


coxph(Surv(times, status)~cluster, data=surv.data)
sfit <- survfit(Surv(times, status)~cluster, data=surv.data)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table = F)


index<-intersect(rownames(surv.clin), tcgaANDme)
index
surv.clin.sel<-surv.clin[index,]
dim(surv.clin.sel)
tcga_blca[index,2]

surv.data<-data.frame(surv.clin.sel[,c(1,3)], tcga_blca[index,2])
names(surv.data)<-c("times", "status", "cluster")
head(surv.data)
prop.table(table(surv.data$cluster, surv.data$status),1)


coxph(Surv(times, status)~cluster, data=surv.data)
sfit <- survfit(Surv(times, status)~cluster, data=surv.data)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table = TRUE)





# 10. CREATE A FINGERPRINT IN THE SPACE OF GENE CONFIGURATIONS
# ==============================================================
##Create fingerprint
merge(contribucion,ind_in_cluster, by='row.names') -> merged_var
ifelse(merged_var$V1=="1",1,0) -> merged_var$cluster1
ifelse(merged_var$V1=="2",1,0) -> merged_var$cluster2
ifelse(merged_var$V1=="3",1,0) -> merged_var$cluster3
merged_var[,2:6]*merged_var$cluster1-> cluster1
rownames(cluster1)<-merged_var$Row.names
merged_var[,2:6]*merged_var$cluster2-> cluster2
rownames(cluster2)<-merged_var$Row.names
merged_var[,2:6]*merged_var$cluster3-> cluster3
rownames(cluster3)<-merged_var$Row.names
finger<-data.frame(matrix(NA, nrow=ncol(matrix_clean_2), ncol=3))
colnames(finger)<-c("cluster1","cluster2", "cluster3")
rownames(finger)<-merged_var$Row.names
as.character(ifelse(merged_var$V1=="1",8,4)) -> finger[,1]
as.character(ifelse(merged_var$V1=="2",8,4)) -> finger[,2]
as.character(ifelse(merged_var$V1=="3",8,4)) -> finger[,3]
finger2<-as.matrix(finger)
merge(matrix_t_clean_2,finger2, by='row.names') -> matrix_t_clean_3_finger
rownames(matrix_t_clean_3_finger)<-matrix_t_clean_3_finger$Row.names
matrix_t_clean_3_finger[,2:ncol(matrix_t_clean_3_finger)]->matrix_t_clean_3_finger_2#[,2:12945]


# MCA for fingerprinting
result_filter_var_finger <- MCA(matrix_t_clean_3_finger_2,
                                quali.sup=(ncol(matrix_t_clean_3_finger_2)-2):(ncol(matrix_t_clean_3_finger_2)), graph = FALSE)


# 11. MATRIX OF DISTANCES
# ==============================================================
# Generate the finger coordinates and the variable cordinates
var_coord<-result_filter_var_finger$var$coord
finger_coord<-result_filter_var_finger$quali.sup$coord
#calculate euclidean distance for each Gene_configuration to each fingerprint
matrix_distance<-data.frame(matrix(NA, nrow=nrow(var_coord), ncol=6))#(matrix(NA, nrow=121390, ncol=6))
colnames(matrix_distance)<-rownames(finger_coord)
rownames(matrix_distance)<-rownames(var_coord)
for (i in 1:nrow(var_coord)){
  for(j in 1:nrow(finger_coord)){
    matrix_distance[i,j]<-sqrt(sum((finger_coord[j,]-var_coord[i,])^2))
  }}
matrix_distance$gene<-as.factor(sapply(strsplit(rownames(matrix_distance), 
                                                split = "_", fixed = TRUE), 
                                       function(x) unlist(x)[1]))

# GENE CONFIGURATION CLOUD OF POINTS (VISUALIZATION)
# ==============================================================

config_universe<-as.data.frame(result_filter_var_finger$var$coord[,1:3])
genes_configs<-as.factor(sapply(strsplit(rownames(result_filter_var_finger$var$coord), 
                                         split = "_", fixed = TRUE), function(x) unlist(x)[2]))
gene_names<-as.factor(sapply(strsplit(rownames(result_filter_var_finger$var$coord), 
                                      split = "_", fixed = TRUE), function(x) unlist(x)[1]))
config_universe$config<-genes_configs
config_universe$gene<-gene_names


#3D MCA plotlibrary(rgl) for configurations
library(rgl)
plot3d(config_universe[,1:3], col=as.integer(config_universe$config),
       size=2, box=FALSE)

plot3d(config_universe[,1:3], col="gray", size=2, box=FALSE)
points3d(result_filter_var_finger$quali.sup$coord[2,1],result_filter_var_finger$quali.sup$coord[2,2],result_filter_var_finger$quali.sup$coord[2,3]  ,pch = 10, size= 10, col= "black")
points3d(result_filter_var_finger$quali.sup$coord[4,1],result_filter_var_finger$quali.sup$coord[4,2],result_filter_var_finger$quali.sup$coord[4,3]  ,pch = 10, size= 10, col= "black")
points3d(result_filter_var_finger$quali.sup$coord[6,1],result_filter_var_finger$quali.sup$coord[6,2],result_filter_var_finger$quali.sup$coord[6,3]  ,pch = 10, size= 10, col= "black")
text3d(result_filter_var_finger$quali.sup$coord[2,1]+0.1,result_filter_var_finger$quali.sup$coord[2,2]+0.1,result_filter_var_finger$quali.sup$coord[2,3]  ,text= c("C1"))
text3d(result_filter_var_finger$quali.sup$coord[4,1]+0.1,result_filter_var_finger$quali.sup$coord[4,2]+0.1,result_filter_var_finger$quali.sup$coord[4,3]  ,text= c("C2"))
text3d(result_filter_var_finger$quali.sup$coord[6,1]+0.1,result_filter_var_finger$quali.sup$coord[6,2]+0.1,result_filter_var_finger$quali.sup$coord[6,3]  ,text= c("C3"))
abclines3d(result_filter_var_finger$quali.sup$coord[2,1:3], a = diag(3), col = "gray")
abclines3d(result_filter_var_finger$quali.sup$coord[4,1:3], a = diag(3), col = "gray")
abclines3d(result_filter_var_finger$quali.sup$coord[6,1:3], a = diag(3), col = "gray")

#3D MCA plotlibrary(rgl) for individuals
plot3d(result_filter_var$ind$coord[,1:3],  col = as.factor(ind_in_cluster[rownames(result_filter_var$ind$coord),]), size = 2, type = "p", box=TRUE)
points3d(config_universe[,1:3], col=as.integer(config_universe$config), size=2)


# CLOUDS OF VARIABLES, CONFIGURATION SPACE CHARACTERIZATION COLLAPSING DATA
# ==============================================================


#GENE CONFIGURATIONS AS MUTATIONS
mut.false<-c(1:8,10:15,17,19:20,22)
mut.true<-c(1:40)[!c(1:40) %in% mut.false]

for (i in 1:nrow(config_universe)) {
  if(config_universe$config[i] %in% mut.false){
    config_universe$Mut[i]<-"FALSE"
  } else if(config_universe$config[i] %in% mut.true){
    config_universe$Mut[i] <-"TRUE"
  }
}

#GENE CONFIGURATIONS AS CNVS
cnv.false<-c(5,7,8,20,21,29,33,36)
cnv.true<-c(1:40)[!c(1:40) %in% cnv.false]
cnv.plus2<-c(10,17,19,25,32,34,37,40)
cnv.plus1<-c(3,4,9,12,15,18, 23, 23, 31)
cnv.less1<-c(6,11,14,16,22,26,27,35)
cnv.less2<-c(1,2,13,24,28,30,38,39)

for (i in 1:nrow(config_universe)) {
  if(config_universe$config[i] %in% cnv.false){
    config_universe$CNV[i]<-"0"
  } else if(config_universe$config[i] %in% cnv.plus2){
    config_universe$CNV[i] <-"+2"
  } else if(config_universe$config[i] %in% cnv.plus1){
    config_universe$CNV[i] <-"+1"
  } else if(config_universe$config[i] %in% cnv.less1){
    config_universe$CNV[i] <-"-1"
  } else if(config_universe$config[i] %in% cnv.less2){
    config_universe$CNV[i] <-"-2"
  }
}


# GENE CONFIGURATIONS AS Methylation
met.low<-c(1,3,6,8,10,16,18,21,24,28)
met.med<-c(2,4,5,9,11,17,27,29,30,32)
met.high<-c(7,13,14,15,19,23,26,33,37)
met.na<-c(12,20,22,24,25,31,35,36,39,40)

for (i in 1:nrow(config_universe)) {
  if(config_universe$config[i] %in% met.low){
    config_universe$MET[i]<-"Low"
  } else if(config_universe$config[i] %in% met.med){
    config_universe$MET[i] <-"Medium"
  } else if(config_universe$config[i] %in% met.high){
    config_universe$MET[i] <-"High"
  } else if(config_universe$config[i] %in% met.na){
    config_universe$MET[i] <-"NO.Val"
  }
}




# some plots3D here.




# 12. DETERMINANT FEATURE SELECTION (1% ARROUND FINGERPRINTS CLUSTERS)
# ==============================================================
identical(rownames(var_coord), rownames(matrix_distance))

#finger1
quantile(as.matrix(matrix_distance$cluster1_8), c(.01))
set.finger1<-as.data.frame(var_coord[which(matrix_distance$cluster1_8 < quantile(as.matrix(matrix_distance$cluster1_8), c(.01))),])
set.finger1$finger<-rep("finger1", nrow(set.finger1))
set.finger1$config<-as.factor(sapply(strsplit(rownames(set.finger1), split = "_"), function(x) unlist(x)[2]))
set.finger1$gene<-as.factor(sapply(strsplit(rownames(set.finger1), split = "_", fixed = TRUE), function(x) unlist(x)[1]))

#finger2
set.finger2<-as.data.frame(var_coord[which(matrix_distance$cluster2_8< quantile(as.matrix(matrix_distance$cluster2_8), c(.01))),])
set.finger2$finger<-rep("finger2", nrow(set.finger2))
dim(set.finger2)
set.finger2$config<-as.factor(sapply(strsplit(rownames(set.finger2), split = "_"), function(x) unlist(x)[2]))
set.finger2$gene<-as.factor(sapply(strsplit(rownames(set.finger2), split = "_", fixed = TRUE), function(x) unlist(x)[1]))

#finger3
set.finger3<-as.data.frame(var_coord[which(matrix_distance$cluster3_8< quantile(as.matrix(matrix_distance$cluster3_8), c(.01))),])
set.finger3$finger<-rep("finger3", nrow(set.finger3))
dim(set.finger3)
set.finger3$config<-as.factor(sapply(strsplit(rownames(set.finger3), split = "_"), function(x) unlist(x)[2]))
set.finger3$gene<-as.factor(sapply(strsplit(rownames(set.finger3), split = "_", fixed = TRUE), function(x) unlist(x)[1]))

#finger23
set.finger23<-as.data.frame(var_coord[which(matrix_distance$cluster1_4< quantile(as.matrix(matrix_distance$cluster1_4), c(.01))),])
set.finger23$finger<-rep("finger23", nrow(set.finger23))
dim(set.finger23)
set.finger23$config<-as.factor(sapply(strsplit(rownames(set.finger23), split = "_"), function(x) unlist(x)[2]))
set.finger23$gene<-as.factor(sapply(strsplit(rownames(set.finger23), split = "_", fixed = TRUE), function(x) unlist(x)[1]))

#finger12
set.finger12<-as.data.frame(var_coord[which(matrix_distance$cluster3_4< quantile(as.matrix(matrix_distance$cluster3_4), c(.01))),])
set.finger12$finger<-rep("finger12", nrow(set.finger12))
dim(set.finger12)
set.finger12$config<-as.factor(sapply(strsplit(rownames(set.finger12), split = "_"), function(x) unlist(x)[2]))
set.finger12$gene<-as.factor(sapply(strsplit(rownames(set.finger12), split = "_", fixed = TRUE), function(x) unlist(x)[1]))

#finger13
set.finger13<-as.data.frame(var_coord[which(matrix_distance$cluster2_4< quantile(as.matrix(matrix_distance$cluster2_4), c(.01))),])
set.finger13$finger<-rep("finger13", nrow(set.finger13))
dim(set.finger13)
set.finger13$config<-as.factor(sapply(strsplit(rownames(set.finger13), split = "_"), function(x) unlist(x)[2]))
set.finger13$gene<-as.factor(sapply(strsplit(rownames(set.finger13), split = "_", fixed = TRUE), function(x) unlist(x)[1]))


# 13. PLOTS3D FINGERPRINTS
# ==============================================================
set.finger<-rbind(set.finger1,set.finger12, set.finger2,
                  set.finger23, set.finger3, set.finger13)
sel.pca3d<-set.finger
set.finger.simply<-rbind(set.finger1, set.finger2,set.finger3)
sel.pca3d.simply<-set.finger.simply

open3d()
points3d(sel.pca3d[,1:3], col="dark grey", box=FALSE, size=3)
points3d(set.finger1[set.finger1$gene %in% genes1,][,1:3], col="dark blue", size=3)


#finger1 vs 23
open3d()
points3d(set.finger1[,1:3], col="blue", box=FALSE, size=3)
points3d(set.finger23[,1:3], col="light blue", box=FALSE, size=3)
#points3d(set.finger[rownames(config.nor.genes.1U23),1:3], col="red", box=FALSE, size=5)
points3d(result_filter_var_finger$quali.sup$coord[2,1],result_filter_var_finger$quali.sup$coord[2,2],result_filter_var_finger$quali.sup$coord[2,3]  ,pch = 10, size= 10, col= "black")
points3d(result_filter_var_finger$quali.sup$coord[1,1],result_filter_var_finger$quali.sup$coord[1,2],result_filter_var_finger$quali.sup$coord[1,3]  ,pch = 10, size= 10, col= "black")
text3d(result_filter_var_finger$quali.sup$coord[2,1]+0.1,result_filter_var_finger$quali.sup$coord[2,2]+0.1,result_filter_var_finger$quali.sup$coord[2,3] +0.2 ,text= c("C1"))
text3d(result_filter_var_finger$quali.sup$coord[1,1],result_filter_var_finger$quali.sup$coord[1,2],result_filter_var_finger$quali.sup$coord[1,3]  +0.2 ,text= c("C2-3"))
abclines3d(result_filter_var_finger$quali.sup$coord[2,1:3], a = diag(3), col = "gray")
abclines3d(result_filter_var_finger$quali.sup$coord[1,1:3], a = diag(3), col = "gray")

open3d()
points3d(set.finger2[,1:3], col="dark green", box=FALSE, size=3)
points3d(set.finger13[,1:3], col="light green", box=FALSE, size=3)
points3d(result_filter_var_finger$quali.sup$coord[4,1],result_filter_var_finger$quali.sup$coord[4,2],result_filter_var_finger$quali.sup$coord[4,3]  ,pch = 10, size= 10, col= "black")
points3d(result_filter_var_finger$quali.sup$coord[5,1],result_filter_var_finger$quali.sup$coord[5,2],result_filter_var_finger$quali.sup$coord[5,3]  ,pch = 10, size= 10, col= "black")
text3d(result_filter_var_finger$quali.sup$coord[4,1]+0.1,result_filter_var_finger$quali.sup$coord[4,2]+0.1,result_filter_var_finger$quali.sup$coord[4,3] +0.2 ,text= c("C2"))
text3d(result_filter_var_finger$quali.sup$coord[5,1],result_filter_var_finger$quali.sup$coord[5,2],result_filter_var_finger$quali.sup$coord[5,3]  +0.2 ,text= c("C1-3"))
abclines3d(result_filter_var_finger$quali.sup$coord[4,1:3], a = diag(3), col = "gray")
abclines3d(result_filter_var_finger$quali.sup$coord[3,1:3], a = diag(3), col = "gray")


open3d()
points3d(set.finger3[,1:3], col="dark red", box=FALSE, size=3)
points3d(set.finger12[,1:3], col="orange", box=FALSE, size=3)
points3d(result_filter_var_finger$quali.sup$coord[6,1],result_filter_var_finger$quali.sup$coord[6,2],result_filter_var_finger$quali.sup$coord[6,3]  ,pch = 10, size= 10, col= "black")
points3d(result_filter_var_finger$quali.sup$coord[5,1],result_filter_var_finger$quali.sup$coord[5,2],result_filter_var_finger$quali.sup$coord[5,3]  ,pch = 10, size= 10, col= "black")
text3d(result_filter_var_finger$quali.sup$coord[6,1]+0.1,result_filter_var_finger$quali.sup$coord[6,2]+0.1,result_filter_var_finger$quali.sup$coord[6,3] +0.2 ,text= c("C3"))
text3d(result_filter_var_finger$quali.sup$coord[5,1],result_filter_var_finger$quali.sup$coord[5,2],result_filter_var_finger$quali.sup$coord[5,3]  +0.2 ,text= c("C1-2"))
abclines3d(result_filter_var_finger$quali.sup$coord[6,1:3], a = diag(3), col = "gray")
abclines3d(result_filter_var_finger$quali.sup$coord[5,1:3], a = diag(3), col = "gray")


open3d()
points3d(set.finger3[,1:3], col="dark red", box=FALSE, size=3)
points3d(set.finger12[,1:3], col="orange", box=FALSE, size=3)
points3d(set.finger1[,1:3], col="blue", box=FALSE, size=3)
points3d(set.finger23[,1:3], col="light blue", box=FALSE, size=3)
points3d(set.finger2[,1:3], col="dark green", box=FALSE, size=3)
points3d(set.finger13[,1:3], col="light green", box=FALSE, size=3)
points3d(result_filter_var_finger$quali.sup$coord[6,1],result_filter_var_finger$quali.sup$coord[6,2],result_filter_var_finger$quali.sup$coord[6,3]  ,pch = 10, size= 10, col= "black")
points3d(result_filter_var_finger$quali.sup$coord[5,1],result_filter_var_finger$quali.sup$coord[5,2],result_filter_var_finger$quali.sup$coord[5,3]  ,pch = 10, size= 10, col= "black")
text3d(result_filter_var_finger$quali.sup$coord[6,1]+0.1,result_filter_var_finger$quali.sup$coord[6,2]+0.1,result_filter_var_finger$quali.sup$coord[6,3] +0.2 ,text= c("F3"))
text3d(result_filter_var_finger$quali.sup$coord[5,1],result_filter_var_finger$quali.sup$coord[5,2],result_filter_var_finger$quali.sup$coord[5,3]  +0.2 ,text= c("f1-3"))
points3d(result_filter_var_finger$quali.sup$coord[4,1],result_filter_var_finger$quali.sup$coord[4,2],result_filter_var_finger$quali.sup$coord[4,3]  ,pch = 10, size= 10, col= "black")
points3d(result_filter_var_finger$quali.sup$coord[3,1],result_filter_var_finger$quali.sup$coord[3,2],result_filter_var_finger$quali.sup$coord[3,3]  ,pch = 10, size= 10, col= "black")
text3d(result_filter_var_finger$quali.sup$coord[4,1]+0.1,result_filter_var_finger$quali.sup$coord[4,2]+0.1,result_filter_var_finger$quali.sup$coord[4,3] +0.2 ,text= c("F2"))
text3d(result_filter_var_finger$quali.sup$coord[3,1],result_filter_var_finger$quali.sup$coord[3,2],result_filter_var_finger$quali.sup$coord[3,3]  +0.2 ,text= c("f1-3"))
points3d(result_filter_var_finger$quali.sup$coord[2,1],result_filter_var_finger$quali.sup$coord[2,2],result_filter_var_finger$quali.sup$coord[2,3]  ,pch = 10, size= 10, col= "black")
points3d(result_filter_var_finger$quali.sup$coord[1,1],result_filter_var_finger$quali.sup$coord[1,2],result_filter_var_finger$quali.sup$coord[1,3]  ,pch = 10, size= 10, col= "black")
text3d(result_filter_var_finger$quali.sup$coord[2,1]+0.1,result_filter_var_finger$quali.sup$coord[2,2]+0.1,result_filter_var_finger$quali.sup$coord[2,3] +0.2 ,text= c("F1"))
text3d(result_filter_var_finger$quali.sup$coord[1,1],result_filter_var_finger$quali.sup$coord[1,2],result_filter_var_finger$quali.sup$coord[1,3]  +0.2 ,text= c("f2-3"))
abclines3d(result_filter_var_finger$quali.sup$coord[6,1:3], a = diag(3), col = "gray")
abclines3d(result_filter_var_finger$quali.sup$coord[4,1:3], a = diag(3), col = "gray")
abclines3d(result_filter_var_finger$quali.sup$coord[2,1:3], a = diag(3), col = "gray")







#3D plot
open3d()
par3d(windowRect = c(1000, 1000, 612, 612))
Sys.sleep(0.1) # Allow sluggish window managers to catch up
parent <- currentSubscene3d()
mfrow3d(1, 2)
points3d(sel.pca3d[,1:3], col=factor(set.finger$config), size=3)
# Points for clusters
points3d(result_filter_var_finger$quali.sup$coord[2,1],result_filter_var_finger$quali.sup$coord[2,2],result_filter_var_finger$quali.sup$coord[2,3]  ,pch = 10, size= 10, col= "black")
points3d(result_filter_var_finger$quali.sup$coord[4,1],result_filter_var_finger$quali.sup$coord[4,2],result_filter_var_finger$quali.sup$coord[4,3]  ,pch = 10, size= 10, col= "black")
points3d(result_filter_var_finger$quali.sup$coord[6,1],result_filter_var_finger$quali.sup$coord[6,2],result_filter_var_finger$quali.sup$coord[6,3]  ,pch = 10, size= 10, col= "black")
# Poitns for intermedian clusters
points3d(result_filter_var_finger$quali.sup$coord[1,1],result_filter_var_finger$quali.sup$coord[1,2], result_filter_var_finger$quali.sup$coord[1,3]  ,pch = 10, size= 10, col= "black") 
points3d(result_filter_var_finger$quali.sup$coord[3,1],result_filter_var_finger$quali.sup$coord[3,2], result_filter_var_finger$quali.sup$coord[3,3]  ,pch = 10, size= 10, col= "black") 
points3d(result_filter_var_finger$quali.sup$coord[5,1],result_filter_var_finger$quali.sup$coord[5,2],result_filter_var_finger$quali.sup$coord[5,3]  ,pch = 10, size= 10, col= "black") 
# Text for points
text3d(result_filter_var_finger$quali.sup$coord[6,1]+0.1,result_filter_var_finger$quali.sup$coord[6,2]+0.1,result_filter_var_finger$quali.sup$coord[6,3] +0.2 ,text= c("F3"))
text3d(result_filter_var_finger$quali.sup$coord[5,1],result_filter_var_finger$quali.sup$coord[5,2],result_filter_var_finger$quali.sup$coord[5,3]  +0.2 ,text= c("f1-3"))
text3d(result_filter_var_finger$quali.sup$coord[4,1]+0.1,result_filter_var_finger$quali.sup$coord[4,2]+0.1,result_filter_var_finger$quali.sup$coord[4,3] +0.2 ,text= c("F2"))
text3d(result_filter_var_finger$quali.sup$coord[3,1],result_filter_var_finger$quali.sup$coord[3,2],result_filter_var_finger$quali.sup$coord[3,3]  +0.2 ,text= c("f1-3"))
text3d(result_filter_var_finger$quali.sup$coord[2,1]+0.1,result_filter_var_finger$quali.sup$coord[2,2]+0.1,result_filter_var_finger$quali.sup$coord[2,3] +0.2 ,text= c("F1"))
text3d(result_filter_var_finger$quali.sup$coord[1,1],result_filter_var_finger$quali.sup$coord[1,2],result_filter_var_finger$quali.sup$coord[1,3]  +0.2 ,text= c("f2-3"))
# Intersect lines in 3d space
abclines3d(result_filter_var_finger$quali.sup$coord[2,1:3], a = diag(3), col = "gray")
abclines3d(result_filter_var_finger$quali.sup$coord[4,1:3], a = diag(3), col = "gray")
abclines3d(result_filter_var_finger$quali.sup$coord[6,1:3], a = diag(3), col = "gray")
# Interactive interface
#M <- r3dDefaults$userMatrix
#fn <- par3dinterp(times = (0:2)*0.75, userMatrix = list(M,
#                                                        rotate3d(M, pi/2, 1, 0, 0),
#                                                        rotate3d(M, pi/2, 0, 1, 0)))
#control <- par3dinterpControl(fn, 0, 3, steps = 15)
#control
#if (interactive())
#  rglwidget(width = 500, height = 500) %>%
#  playwidget(control,
#             step = 0.01, loop = TRUE, rate = 0.2)
next3d(reuse = FALSE)
sublegend<-configlegend[levels(as.factor(set.finger$config)),1:3]
sublegend.plot<-c()
for (i in 1:nrow(sublegend)){
  sublegend.plot[i]<-paste(unlist(sublegend[i,]), collapse =" ")
}
legend3d("center", sublegend.plot, pch = c(1, 1, 1), col =  levels(as.factor(set.finger$config)))
useSubscene3d(parent)
sublegend.plot

# 14. CHARACTERIZE FINGERPRINTS
# ==============================================================

identical((rownames(config_universe[rownames(set.finger),])), (rownames(set.finger)))
set.finger$Mut<-as.factor((config_universe[rownames(set.finger),])$Mut)
set.finger$CNV<-as.factor((config_universe[rownames(set.finger),])$CNV)
set.finger$MET<-as.factor((config_universe[rownames(set.finger),])$MET)
set.finger.simply$Mut<-as.factor((config_universe[rownames(set.finger.simply),])$Mut)
set.finger.simply$CNV<-as.factor((config_universe[rownames(set.finger.simply),])$CNV)
set.finger.simply$MET<-as.factor((config_universe[rownames(set.finger.simply),])$MET)


round(prop.table(table(set.finger$finger, set.finger$Mut),1),2)
round(prop.table(table(set.finger$finger, set.finger$CNV),1),2)
round(prop.table(table(set.finger$finger, set.finger$MET),1),2)

round(prop.table(table(set.finger.simply$finger, set.finger.simply$Mut),1),2)
round(prop.table(table(set.finger.simply$finger, set.finger.simply$CNV),1),2)
round(prop.table(table(set.finger.simply$finger, set.finger.simply$MET),1),2)

# Normal
# as: MUT FALSE, then: low, medium and high methylathed and inside as -2, -1, 0, 1, 2 cnv:
list<-list(c("1","6","8","3","10"), c("2","11","5","4","17"), c("13","14","7","15","19"))
countfinger<-table(set.finger.simply$finger, set.finger.simply$config)
library(scales)
for (i in 1:length(list)){
  df.m<-melt(countfinger[,unlist(list[i])],id.var="config")
  df.m$Var2<-ordered(factor(as.character(df.m$Var2)), unlist(list[i]))
  p<-ggplot(data=df.m, aes(x=Var2, y=value, fill=Var1)) + geom_bar(position = "fill",stat="identity") + scale_y_continuous(labels = percent_format()) 
  outname<-paste(i, "finger.pdf")
  ggsave(outname, plot = p)
}
df.m <-ordered(df.m$Var2, levels=levels(df.m$Var2)[unclass(df.m$Var2)])
x
# Mutated
# as: MUT TRUE, then: low, medium and high methylathed and inside as -2, -1, 0, 1, 2 cnv:
list<-list(c("28","16","21","18","34"), c("27","29","9","32"), c("26","33","23"))
for (i in 1:length(list)){
  df.m<-melt(countfinger[,unlist(list[i])],id.var="config")
  df.m$Var2<-ordered(factor(as.character(df.m$Var2)), unlist(list[i]))
  p<-ggplot(data=df.m, aes(x=Var2, y=value, fill=Var1)) + geom_bar(position = "fill",stat="identity") + scale_y_continuous(labels = percent_format()) 
  outname<-paste(i, "M-finger.pdf")
  ggsave(outname, plot = p)
}




#Proportions
round(prop.table(table(set.finger1$finger, set.finger1$config)),2)*100
round(prop.table(table(set.finger2$finger, set.finger2$config)),2)*100
round(prop.table(table(set.finger3$finger, set.finger3$config)),2)*100



# 15. CODE TO FURTHER FUNCIONAL ANALYSIS
# ==============================================================

gene.universe<-colnames(matrix_t_clean_2)
library(biomaRt)
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",  host="www.ensembl.org") 
go.data <- getBM(attributes=c("go_id", 'hgnc_symbol'),
                 filters = 'hgnc_symbol', values = gene.universe, mart = ensembl)
gene2GO <- split(go.data$go_id,  go.data$hgnc_symbol)
gene2GO <- lapply(gene2GO, unique)  # to remove duplicates


library(topGO)

#This is a function to retreive only those genes that are common to two fingers,
#not common to all to increase the significance of our selection.
record<-c()
for (gen in set.finger$gene){
  if (length((set.finger[set.finger$gene==gen,])$finger) ==2){
    record=unique(c(record,gen))
  }
}
length(set.finger$gene)#[1] 5577
length(unique(set.finger$gene))#[1] 4497
length(record)#[1] 816


#This is a function to extract the configuratons for each grup from a list of unique genes

#FUNCTION TO EXTRACT CONFIG LIST
config_fx<-function(genes){
  list.genes<-c()
  for(i in 1:length(genes)){list.genes[i]<-list(rownames(set.finger[set.finger$gene==genes[i],]))}
  config.genes<-set.finger[unlist(list.genes),6:7]
  config.genes$config<-factor(config.genes$config)
  return(config.genes)
}


#This function allows to get the terms and the genes for these term, in our list
list.GOgenes<-function(sig.x, GOdat, set.gene){
  rownames(sig.x)<-sig.x[,1]
  for (term in sig.x[,1]){
    gene.sig<-intersect(unlist(genesInTerm(GOdat, term)), set.gene)
    termdesc<-(sig.x[term,2])
    mygenesforterm<-paste(gene.sig, collapse = ", ")
    print(paste(term, termdesc, "Genes:", mygenesforterm))
  }
}





GOdata_fx<-function(genes, geneuniverse){
  genelist<-factor(as.integer(geneuniverse %in% genes))
  names(genelist)<-geneuniverse
  GOdata<-new("topGOdata", description="My project", ontology="BP", 
              allGenes=genelist,geneSel=as.vector(genes), annot = annFUN.gene2GO, gene2GO = gene2GO)
  return(GOdata)
}



topGO_fx<-function(GOdata){
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  allGO= usedGO(object = GOdata)
  a<-GenTable(GOdata, classicFisher = resultFisher)
  return(a)
}




# ==============================================================
# 16. LIST FOR FINGERS AND ICS
# ==============================================================


#---------------------------
# GENE LISTS:FINGER 1
#---------------------------

#FINGER1
#Select configurations in finger1
genes1<-unique(set.finger1$gene)
GOdata1<-GOdata_fx(genes1,gene.universe)
sig.finger1<-topGO_fx(GOdata1)
sig.finger1<-sig.finger1[as.numeric(sig.finger1$classicFisher)<0.001,]
list.GOgenes(sig.finger1, GOdata1, genes1)
#Bonferroni correction p=0.05/num genes to GOterms
sig.finger1


# FINGER1 NORMAL configurations
set.finger1_normal<-set.finger1[set.finger1$config=="8" | set.finger1$config=="5",]
genes1_normal<-unique(set.finger1_normal$gene)
GOdata.finger1.normal<-GOdata_fx(genes1_normal,gene.universe)
finger1.sig.normal<-topGO_fx(GOdata.finger1.normal)
finger1.sig.normal<-finger1.sig.normal[finger1.sig.normal$classicFisher<0.001,]
list.GOgenes(finger1.sig.normal, GOdata.finger1.normal, genes1.normal)
finger1.sig.normal 


#Select only FINGER1 ALTERED configurations
set.finger1_altered<-set.finger1[!rownames(set.finger1) %in% rownames(set.finger1_normal) ,]
genes1_altered<-unique(set.finger1_altered$gene)
GOdata.finger1.altered<-GOdata_fx(genes1_altered,gene.universe)
finger1.sig.altered<-topGO_fx(GOdata.finger1.altered)
finger1.sig.altered<-finger1.sig.altered[finger1.sig.altered$classicFisher<0.001,]
list.GOgenes(finger1.sig.altered, GOdata.finger1.altered, genes1_altered)
finger1.sig.altered # no significant terms

#Select only FINGER1 MUTATED configurations

set.finger1_mut<-set.finger1[set.finger1$config=="21" | set.finger1$config=="29",]
genes1.mut<-unique(set.finger1_mut$gene)
GOdata1.mut<-GOdata_fx(genes1.mut,gene.universe)
finger1.sig.mut<-topGO_fx(GOdata1.mut)
finger1.sig.mut


# FINGER23
genes23<-unique(set.finger23$gene)
GOdata23<-GOdata_fx(genes23,gene.universe)
sig.finger23<-topGO_fx(GOdata23)
rownames(sig.finger23)<-sig.finger23$Term


# GENE LISTS: ICs 1
#---------------------------

ics1<-intersect(genes1, record)
GOdata.ics1<-GOdata_fx(genes1,gene.universe)
sig.1<-topGO_fx(GOdata.ics1)
rownames(sig.1)<-sig.1$Term
sig.1

# ICS1 normal configurations
ics1_normal<-intersect(genes1_normal,record)
GOdata1_nor<-GOdata_fx(ics1_normal,gene.universe)
sig.1n<-topGO_fx(GOdata1_nor)
rownames(sig.1n)<-sig.1n$Term
sig.1n

# ICS1 altered configurations
ics1_altered<-intersect(genes1_altered,record)
GOdata1_alt<-GOdata_fx(ics1_altered,gene.universe)
sig.1a<-topGO_fx(GOdata1_alt)
rownames(sig.1a)<-sig.1a$Term
sig.1a

# GENE LISTS: NICs1 = filtered finger23
#---------------------------
#Select configurations in finger23
ics23<-intersect(genes23, record)
GOdata23<-GOdata_fx(ics23,gene.universe)
sig.23<-topGO_fx(GOdata23)
rownames(sig.23)<-sig.23$Term
sig.23




#---------------------------
# GENE LISTS:FINGER 2
#---------------------------

#Select configurations in finger2
genes2<-unique(set.finger2$gene)
GOdata2<-GOdata_fx(genes2,gene.universe)
sig.finger2<-topGO_fx(GOdata2)
#To retreive whole table of topGOdata:
#allGO1 = usedGO(object = GOdata1)
#allGO1<-GenTable(GOdata1, classicFisher = resultFisher1, topNodes = lenght(allGO1)
#Bonferroni correction p=0.05/num genes to GOterms
sig.finger2


#FINGER13
#Select configurations in finger13
genes13<-unique(set.finger13$gene)
GOdata13<-GOdata_fx(genes13,gene.universe)
sig.finger13<-topGO_fx(GOdata13)
rownames(sig.finger13)<-sig.finger13$Term





# GENE LISTS: ICs 2
#---------------------------

ics2<-intersect(genes2, record)
GOdata2<-GOdata_fx(ics2,gene.universe)
sig.2 <- topGO_fx(GOdata2)
rownames(sig.2)<-sig.2$Term

#Deleted ics2
set.finger2_main<-set.finger2[set.finger2$config=="6"| set.finger2$config=="11"| set.finger2$config=="14",]
f2_main<-unique(set.finger2_main$gene)
GOdata2c<-GOdata_fx(f2_main,gene.universe)
sig.2c<-topGO_fx(GOdata2c)
rownames(sig.2c)<-sig.2c$Term
sig.2c


#Other ics2
set.finger2_other<-set.finger2[!rownames(set.finger2) %in% rownames(set.finger2_main) ,]
f2_other<-set.finger2_other$gene
length(f2_other)
GOdata2u<-GOdata_fx(f2_other,gene.universe)
sig.2u<-topGO_fx(GOdata2u)
rownames(sig.2u)<-sig.2u$Term
sig.2u




# GENE LISTS: NICs2 = filtered finger13
#---------------------------
ics13<-intersect(genes13, record)
GOdata13<-GOdata_fx(ics13,gene.universe)
sig.13<-topGO_fx(GOdata13)
rownames(sig.13)<-sig.13$Term
sig.13






#---------------------------
# GENE LISTS:FINGER3
#---------------------------
#FINGER3
#Select configurations in finger3
genes3<-unique(set.finger3$gene)
GOdata3.finger<-GOdata_fx(genes3,gene.universe)
sig.finger3<-topGO_fx(GOdata3.finger)



#FINGER12: 
#Select configurations in finger12
genes12<-unique(set.finger12$gene)
GOdata12.finger<-GOdata_fx(genes12,gene.universe)
sig.finger12<-topGO_fx(GOdata12.finger)



# GENE LISTS: ICs 3
#---------------------------

ics.3<-intersect(genes3, record)
GOdata3<-GOdata_fx(ics.3,gene.universe)
sig.3 <- topGO_fx(GOdata3)
rownames(sig.3)<-make.names(sig.3$Term, unique=TRUE)
sig.3


#Duplicated configurations
set.finger3_main<-set.finger3[set.finger3$config=="3"| set.finger3$config=="4",]
f3_main<-unique(set.finger3_main$gene)
GOdata3c<-GOdata_fx(f3_main,gene.universe)
sig.3c <-topGO_fx(GOdata3c)
rownames(sig.3c)<-sig.3c$Term
sig.3c

#Other
set.finger3_other<-set.finger3[!rownames(set.finger3) %in% rownames(set.finger3_main) ,]
f3_other<-set.finger3_other$gene
GOdata3u<-GOdata_fx(f3_other,gene.universe)
sig.3u<-topGO_fx(GOdata3u)
rownames(sig.3u)<-sig.3u$Term
sig.3u



# GENE LISTS: NICs3 = filtered finger12
#---------------------------

ics12<-intersect(genes12, record)
GOdata12<-GOdata_fx(genes12,gene.universe)
sig.12<-topGO_fx(GOdata12)
rownames(sig.12)<-sig.12$Term
sig.12







#  CODE TO LIST GENES:

t3.gene.sig <- intersect(unlist(genesInTerm(GOdata3c, "GO:0051606")),set.finger3$gene)
t3.gene.sig
paste(t3.gene.sig, collapse=', ')
config_fx(t3.gene.sig)








# ==============================================================
# 16. LIST FOR DCS
# ==============================================================




# ICs 1 common with 23
genes.1U23<-intersect(ics1, ics23)
config.genes.1U23<-config_fx(genes.1U23)
GOdata1U23<-GOdata_fx(genes.1U23,gene.universe)
sig.1U23<-topGO_fx(GOdata1U23)
rownames(sig.1U23)<-sig.1U23[,2]



# ICs 2 common with 13
genes.2U13<-intersect(ics2, ics13)
config.genes.2U13<-config_fx(genes.2U13)
GOdata.2u13<-GOdata_fx(genes.2U13, gene.universe)
sig.2U13<-topGO_fx(GOdata.2u13)
rownames(sig.2U13)<-sig.2U13[,2]



#Finger 3 common with 12
genes.3U12<-intersect(ics3, ics12)
config.genes.3U12<-config_fx(genes.3U12)
GOdata.3U12<-GOdata_fx(genes.3U12, gene.universe)
sig.3U12<-topGO_fx(GOdata.3U12)
rownames(sig.3U12)<-make.names(sig.3U12[,2], unique=T)







#Proportions in DCS
round(prop.table(table(config.genes.1U23$finger, config.genes.1U23$config)[1,]),2)*100
round(prop.table(table(config.genes.2U13$finger, config.genes.2U13$config)[1,]),2)*100
round(prop.table(table(config.genes.3U12$finger, config.genes.3U12$config)[1,]),2)*100




# HEATMAP ALL GO TERMS
sig.1
sig.1n
sig.1a
sig.23
sig.1U23
sig.2
#sig.2U13<-sig.2U13[c(1:6),], only to plot
sig.2U13
sig.3
sig.3U12


c.GOsig.universe<-unique(c(rownames(sig.1n), rownames(sig.1a), rownames(sig.1U23), rownames(sig.2), rownames(sig.2U13), rownames(sig.3), rownames(sig.3U12)))
length(c.GOsig.universe) #63

c.GOTerm_matrix<-matrix(nrow=length(c.GOsig.universe), ncol=7)
rownames(c.GOTerm_matrix)<-c.GOsig.universe
colnames(c.GOTerm_matrix)<-as.character(seq(1:7))

for(term in c.GOsig.universe){
  c.GOTerm_matrix[term,1]<-as.numeric(sig.1U23[term,6])
  c.GOTerm_matrix[term,2]<-as.numeric(sig.1a[term,6])
  c.GOTerm_matrix[term,3]<-as.numeric(sig.1n[term,6])
  c.GOTerm_matrix[term,4]<-as.numeric(sig.2U13[term,6])
  c.GOTerm_matrix[term,5]<-as.numeric(sig.2[term,6])
  c.GOTerm_matrix[term,6]<-as.numeric(sig.3U12[term,6])
  c.GOTerm_matrix[term,7]<-as.numeric(sig.3[term,6])
}
cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(colnames(c.GOTerm_matrix))]
#hmcol<-colorRampPalette(c("white","red","red4"))(n = 299)[20:180]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
#hmcol<-colorRampPalette(c("white","red","red4"))(n = 299)[20:180]
heatmap.2(data.matrix(c.GOTerm_matrix),col=rev(hmcol), colsep=c(1:6),rowsep=(1:62),
          cellnote=ifelse(c.GOTerm_matrix==0, NA, as.matrix(c.GOTerm_matrix)),notecol="white",
          sepwidth=c(0.05,0.05), sepcolor="white", trace="none", ColSideColors=cols,
          Rowv=F,Colv=F, scale="none",dendrogram="none",key=F, lhei = c(0.05,5),margins=c(5,20))







#HEATMAP FOR ALL THE COMMON SHARED CONFIGURATIONS
GOsig.universe<-c(sig.1U23$Term, sig.2U13$Term, rownames(sig.3U12))
length(GOsig.universe) #150

GOTerm_matrix<-matrix(nrow=length(GOsig.universe), ncol=3)
rownames(GOTerm_matrix)<-GOsig.universe
colnames(GOTerm_matrix)<-c("1U23", "2U13", "3U12")

for(term in GOsig.universe){
  GOTerm_matrix[term,1]<-as.numeric(sig.1U23[term,6])
  GOTerm_matrix[term,2]<-as.numeric(sig.2U13[term,6])
  GOTerm_matrix[term,3]<-as.numeric(sig.3U12[term,6])
}

GOTerm_matrix[is.na(GOTerm_matrix)]<-0
head(GOTerm_matrix)
cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(colnames(GOTerm_matrix))]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

rownames(GOTerm_matrix)
GOTerm_matrix<-GOTerm_matrix[-c(17:20),]
heatmap.2(GOTerm_matrix, 
          trace="none", na.color="red",
          ColSideColors=cols,
          cexCol = 1,
          col=hmcol, margins=c(5,20), notecex = 0.5)


heatmap.2(data.matrix(c1.GOTerm_matrix),trace="none", na.color=hmcol[2],
          ColSideColors=cols,
          cexCol = 1,
          col=rev(hmcol), margins=c(5,20), notecex = 0.5, Rowv=F,Colv=F)


#Finger1 vs 1u23: SHOW THE DISREGULATED BIOLOGICAL PROCESES BETWEEN GROUPS 1 AND 23
c1.GOsig.universe<-unique(c(rownames(sig.1U23), rownames(sig.1n),rownames(sig.1a)))
length(c1.GOsig.universe) #63

c1.GOTerm_matrix<-matrix(nrow=length(c1.GOsig.universe), ncol=3)
rownames(c1.GOTerm_matrix)<-c1.GOsig.universe
colnames(c1.GOTerm_matrix)<-c("f1N ", "f1U23-NOR", "j")

for(term in c1.GOsig.universe){
  c1.GOTerm_matrix[term,1]<-as.numeric(sig.1U23[term,6])
  c1.GOTerm_matrix[term,2]<-as.numeric(sig.1n[term,6])
  c1.GOTerm_matrix[term,3]<-as.numeric(sig.1a[term,6])
}
cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(colnames(c1.GOTerm_matrix))]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
#hmcol<-colorRampPalette(c("white","red","red4"))(n = 299)[20:180]
heatmap.2(data.matrix(c1.GOTerm_matrix),col=rev(hmcol), colsep=c(1:6),rowsep=(1:62),
          cellnote=ifelse(c1.GOTerm_matrix==0, NA, as.matrix(c1.GOTerm_matrix)),notecol="white",
          sepwidth=c(0.05,0.05), sepcolor="white", trace="none", ColSideColors=cols,
          Rowv=F,Colv=F, scale="none",dendrogram="none",key=F, lhei = c(0.05,5),margins=c(5,20))



#Finger2 vs 2u13: SHOW THE DISREGULATED BIOLOGICAL PROCESES BETWEEN GROUPS 1 AND 23
sig.2U13<-sig.2U13[c(1:6),]
sig.2
c2.GOsig.universe<-unique(c(rownames(sig.2U13), rownames(sig.2)))
length(c2.GOsig.universe) #63

c2.GOTerm_matrix<-matrix(nrow=length(c2.GOsig.universe), ncol=2)
rownames(c2.GOTerm_matrix)<-c2.GOsig.universe
colnames(c2.GOTerm_matrix)<-c("f1N ", "f1U23-NOR")

for(term in c2.GOsig.universe){
  c2.GOTerm_matrix[term,1]<-as.numeric(sig.2U13[term,6])
  c2.GOTerm_matrix[term,2]<-as.numeric(sig.2[term,6])
}
cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(colnames(c2.GOTerm_matrix))]
hmcol<-colorRampPalette(c("white","red","red4"))(n = 299)[20:180]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
#hmcol<-colorRampPalette(c("white","red","red4"))(n = 299)[20:180]
heatmap.2(data.matrix(c2.GOTerm_matrix),col=rev(hmcol), colsep=c(1:6),rowsep=(1:62),
          cellnote=ifelse(c2.GOTerm_matrix==0, NA, as.matrix(c2.GOTerm_matrix)),notecol="white",
          sepwidth=c(0.05,0.05), sepcolor="white", trace="none", ColSideColors=cols,
          Rowv=F,Colv=F, scale="none",dendrogram="none",key=F, lhei = c(0.05,5),margins=c(5,20))




#Finger3 vs 3u12: SHOW THE DISREGULATED BIOLOGICAL PROCESES BETWEEN GROUPS 2 AND 13

c3.GOsig.universe<-unique(c(rownames(sig.3U12), rownames(sig.3)))
length(c2.GOsig.universe) #63

c3.GOTerm_matrix<-matrix(nrow=length(c3.GOsig.universe), ncol=2)
rownames(c3.GOTerm_matrix)<-c3.GOsig.universe
colnames(c3.GOTerm_matrix)<-c("f1N ", "f1U23-NOR")

for(term in c3.GOsig.universe){
  c3.GOTerm_matrix[term,1]<-as.numeric(sig.3U12[term,6])
  c3.GOTerm_matrix[term,2]<-as.numeric(sig.3[term,6])
}
cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(colnames(c3.GOTerm_matrix))]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
#hmcol<-colorRampPalette(c("white","red","red4"))(n = 299)[20:180]
heatmap.2(data.matrix(c3.GOTerm_matrix),col=rev(hmcol), colsep=c(1:6),rowsep=(1:62),
          cellnote=ifelse(c3.GOTerm_matrix==0, NA, as.matrix(c3.GOTerm_matrix)),notecol="white",
          sepwidth=c(0.05,0.05), sepcolor="white", trace="none", ColSideColors=cols,
          Rowv=F,Colv=F, scale="none",dendrogram="none",key=F, lhei = c(0.05,5),margins=c(5,20))









# To extract Terms
sig.1n
t1.gene.sig <- intersect(unlist(genesInTerm(GOdata1_nor, "GO:0009314")),set.finger1$gene)
t1.gene.sig
paste(t1.gene.sig, collapse=', ')
config_fx(t1.gene.sig)

# ==============================================================
# 15. PROPORTIONS IN GENE LISTS
# ==============================================================

# FINGER1, ICS1, DCS1


prop.table(table(set.finger1$config))*100
table(set.finger1$finger, set.finger1$config)
round(prop.table(table(set.finger1$finger, set.finger1$config)),2)*100

#Normal 
prop.table(table((set.finger[set.finger$config=="8" | set.finger$config=="5",])$finger))*100
chisq.test(table((set.finger[set.finger$config=="8" | set.finger$config=="5",])$finger))#X-squared = 10.683, df = 2, p-value = 0.004788





#Mutations
length((set.finger1[set.finger1$config=="21" | set.finger1$config=="29",])$gene)
prop.table(table((set.finger.simply[set.finger.simply$config=="21" | set.finger.simply$config=="29",])$finger))*100
chisq.test(table((set.finger[set.finger$config=="21" | set.finger$config=="29",])$finger))#X-squared = 201.79, df = 4, p-value < 2.2e-16


#Deletions
length((set.finger2[set.finger2$CNV=="-1" | set.finger2$CNV=="-2",])$gene)
prop.table(table((set.finger[set.finger$CNV=="+1" | set.finger$CNV=="+2",])$finger))*100
chisq.test(table((set.finger[set.finger$CNV=="+1" | set.finger$CNV=="+2",])$finger))#X-squared = 201.79, df = 4, p-value < 2.2e-16



#Hypermethylated configurations
(set.finger1[set.finger1$config=="7",])$gene
#Deleted configurations
(set.finger1[set.finger1$config=="11",])$gene



# Proportions ICs

c1.nor<-config_fx(ics1_normal)
round(prop.table(table(c1.nor$finger, c1.nor$config)),2)*100
c1.alt<-config_fx(ics1_altered)
round(prop.table(table(c1.alt$finger, c1.alt$config)),2)*100


# ==============================================================
# 15. CNV CHROMOSOME ANNOTATION ENRICHMENT 
# ==============================================================


ics1.mutated<-as.vector(set.finger1[set.finger1$Mut=="TRUE",]$gene)


# CNV SPLITED IDEOGRAMS
set.finger1<-set.finger[rownames(set.finger1),]
ics1.duplicated<-as.vector(set.finger1[set.finger1$CNV== "+1"|set.finger1$CNV=="+2",]$gene)
ics1.deleted<-as.vector(set.finger1[set.finger1$CNV== "-1"|set.finger1$CNV=="-2",]$gene)
set.finger2<-set.finger[rownames(set.finger2),]
ics2.duplicated<-as.vector(set.finger2[set.finger2$CNV== "+1"|set.finger2$CNV=="+2",]$gene)
ics2.deleted<-as.vector(set.finger2[set.finger2$CNV== "-1"|set.finger2$CNV=="-2",]$gene)
set.finger3<-set.finger[rownames(set.finger3),]
ics3.duplicated<-as.vector(set.finger3[set.finger3$CNV== "+1"|set.finger3$CNV=="+2",]$gene)
ics3.deleted<-as.vector(set.finger3[set.finger3$CNV== "-1"|set.finger3$CNV=="-2",]$gene)




write.table(ics1.deleted, file="ics1_deleted.txt", sep=",", append=T, col.names = F, row.names = F)
write.table(ics1.duplicated, file="ics1_duplicated.txt", sep=",", append=T, col.names = F, row.names = F)
write.table(ics2.deleted, file="ics2_deleted.txt", sep=",", append=T, col.names = F, row.names = F)
write.table(ics2.duplicated, file="ics2_duplicated.txt", sep=",", append=T, col.names = F, row.names = F)
write.table(ics3.deleted, file="ics3_deleted.txt", sep=",", append=T, col.names = F, row.names = F)
write.table(ics3.duplicated, file="ics3_dupicated.txt", sep=",", append=T, col.names = F, row.names = F)
write.table(ics1.mutated, file="ics1_mutated.txt", sep=",", append=T, col.names = F, row.names = F)

# eSTO QUE SIGUE ESTA A MEDIAS, ES PARA PLOTEAR EL KARIOGRAMA, PERO NO PUEDO CARGAR LAS LIBRERIAS.

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
all.genes<-genes(txdb)
head(all.genes)
columns(txdb)

ics3.dup.anot<-getBM(attributes=c("hgnc_symbol","entrezgene"),
filters= "hgnc_symbol",
values= ics3.duplicated,
mart=ensembl)

sel<-as.data.frame(select(txdb, keys=as.character(ics3.dup.anot$entrezgene), columns=c("TXCHROM", "TXSTART", "TXEND", "TXSTRAND"), keytype="GENEID"))
head(sel)
sel<-is.na(sel)
reference_GRange <- GRanges(seqnames= sel$TXCHROM,IRanges(start= sel$TXSTART,end= sel$TXEND),h.gene=sel$GENEID) 



#chr 1 is automatically drawn by default (subchr="chr1")
p.ideo <- Ideogram(genome = "hg19")
p.ideo

autoplot(reference_GRange, layout = "karyogram", cytoband = TRUE)


#Highlights a region on "chr2"
p.ideo + xlim(GRanges("chr2", IRanges(1e8, 1e8+10000000)))



# LIST to DAVID ENRICHMENT

setwd("~/omics/code to TFM/data/Data saved objects/DAVID list")
ics2_del.annot<-read.table("ics2_del_resDAVID.txt", header = T, sep="\t", fill = T)
ics2_del.annot<-ics2_del.annot[ics2_del.annot$Category=="CYTOBAND",]
ics2_del.sig<-ics2_del.annot[ics2_del.annot$Benjamini<0.05,]
ics2_del.sig$Chr<-sapply(ics2_del.sig$Term, function(x) substr(x,1,2))
ics2_del.sig.19<-ics2_del.sig[ics2_del.sig$Chr=="19",]


head(ics2_del.sig)
names(ics2_del.sig)
dim(ics2_del.sig)

ics3_dup.annot<-read.table("ics3_dup_resDAVID.txt", header = T, sep="\t", fill = T)

# 15. EXPRESSION DATA 
# ==============================================================


# Expression data:
setwd("~/omics/code to TFM/data/Expression data")
expression<-read.table("BLCA_rna.txt", header=TRUE, row.names=1, sep="\t", check.names = FALSE)
dim(expression) # [1] 20531   427


exp<-expression[,which(substr(colnames(expression), 14,16) == "01A")]
dim(exp) # [1] 20532   405
head(exp)
exp<-expression[c(2:20532), ]
# Patients names 
colnames(exp)<-substr(colnames(exp), 6,12)
colnames(exp)<-gsub("\\.", "-", colnames(exp))
exp.matrix<-exp[, intersect(rownames(matrix_t_clean_2), names(exp))]
dim(exp.matrix)# [1] 20531   218



# Gene names
exp.filt.names<-as.vector(sapply(rownames(exp.filt), function(x) unlist(strsplit(x, split = "|", fixed = TRUE))[1]))
exp.filt<-exp.matrix[!duplicated(exp.filt.names), ]
exp.index.genes<-which(exp.filt.names%in%colnames(matrix_t_clean_2))


# Matrix like expression set
exp.filt<-exp.filt[exp.index.genes,]
dim(exp.filt) # [1] 10244   218
rownames(exp.filt)<-as.vector(sapply(rownames(exp.filt), function(x) unlist(strsplit(x, split = "|", fixed = TRUE))[1]))
countdata<-data.matrix(exp.filt)

countdata<-countdata[,rownames(cluster.group[order(cluster.group$cluster),])]


ics<-unique(c(genes1, genes2, genes3,  genes12, genes13, genes23))
dcs<-c(genes.1U23, genes.2U13, genes.3U12)

finger.countdata<-countdata[intersect(ics,rownames(countdata)),]
dim(finger.countdata)
#[1] 815 218




# UNSUPERVISED CLUSTERING FOR MCA AND EXPRESSION IN ICS: COMPARISION
# ==============================================================



# cluster MCA
clust.hcpc <- as.numeric(cluster_var$data.clust$clust)
clust.cutree <- dendextend:::cutree(cluster_var$call$t$tree, k=3, order_clusters_as_data = FALSE)
idx <- order(as.numeric(names(clust.cutree)))
clust.cutree <- clust.cutree[idx]
( tbl <- table(clust.hcpc, clust.cutree) )
( lbls <- apply(tbl,2,which.max) )
cluster_var$call$t$tree %>% 
  color_branches(k=3, groupLabels =lbls) %>% 
  set("labels_cex", .5) %>% 
  plot(horiz=T) 
mca.dend<-as.dendrogram(cluster_var$call$t$tree)


# Hieralchical kmeans clustering for EXPRESSION DATA

d.exp<-dist(t(finger.countdata))
hc.exp<-hclust(d.exp, method="ward.D")
plot(hc.exp, hang = -1)
# Determine the optimal clusts
k.max <- 15 # Maximal number of clusters
fviz_nbclust(finger.countdata[intersect(set.finger$gene,rownames(finger.countdata)),], hcut, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)
# Expression Hkmeans for ICs
hkm.exp<-hkmeans(t(as.matrix(finger.countdata)),3)
hkmeans_tree(hkm.exp)
dend_exp<-as.dendrogram(hkm.exp$hclust)
dend_exp %>% 
  color_branches(k=3, groupLabels =lbls) %>% 
  set("labels_cex", .5) %>% 
  plot(horiz=T) 



#Compare dendograms 
dl.EXP.CLUST <- dendlist(highlight_branches(mca.dend), highlight_branches(dend_exp))
tanglegram(dl.EXP.CLUST, sort = T, common_subtrees_color_lines = T, highlight_distinct_edges  = F, highlight_branches_lwd = F)

tanglegram(mca.dend,dend_exp, color_lines=(cluster_var$call$X$clust),
           common_subtrees_color_lines = T, highlight_distinct_edges  = F, highlight_branches_lwd = F)



# CLUSTER VECTORS FOR REPRESENTATION IN BARS AND ORDERING:

# Expression type
expresion.group<-as.data.frame(hkm.exp$cluster)
dim(expresion.group)

# Cluster type
cluster.group<-data2index[rownames(expresion.group),]
dim(cluster.group)

identical(rownames(expresion.group), rownames(cluster.group))

# Clinical type
exp.clin.data<-data4mfa[rownames(cluster.group),]
identical(rownames(cluster.group), rownames(exp.clin.data))


cluster.cor<-data.frame(expresion.group$`hkm.exp$cluster`, cluster.group$cluster, exp.clin.data$clinical.BLCA.diagnosis_subtype)
rownames(cluster.cor)<-rownames(cluster.group)

# Vectors:
cluster.type<-cluster.cor$cluster.group.cluster
expresion.type<-cluster.cor$expresion.group..hkm.exp.cluster.
clinical.type<-cluster.cor$exp.clin.data.clinical.BLCA.diagnosis_subtype

# Comparision for clusters and expression clusters. 
expANDclus<-cluster.cor
colnames(expANDclus)<-c("exp", "mca", "clin")
table(expANDclus$exp, expANDclus$mca)

xtable(round(prop.table(table(expANDclus$exp, expANDclus$mca)),2)*100)


round(prop.table(table(expANDclus$mca, expANDclus$exp)[1,]),2)*100
round(prop.table(table(expANDclus$mca, expANDclus$exp)[2,]),2)*100
round(prop.table(table(expANDclus$mca, expANDclus$exp)[3,]),2)*100

#Otro tipo de representacion (expresion de ics)

# Set up colour vector for celltype variable in clinical: Papillary, non-papillary
col.cell.clin <- c("purple","orange")[clinical.type]
# Set up colour vector for celltype variable
col.cell <- c("black","red", "green")[cluster.type]


### plots horiz version:
par(mar = c(4,1,1,12))
dend_exp %>% set("labels_cex", 0.5) %>% set("labels_col", value = c("purple","orange","green"), k=3) %>% 
  plot(main = "Color labels \nper cluster", horiz =T)
abline(v = 7, lty = 2)
colored_bars(cbind(cluster.type, clinical.type), dend_exp, rowLabels = c("cluster", "clinical"),horiz = TRUE)





# 16. UNSUPERVISED MULTISCALING
# ==============================================================





groups<-factor(cluster.type)
# Countdata object sorted as cluster group to latter match with the cluster group in ColSideColoros 
finger.countdata<-finger.countdata[,rownames(cluster.cor)]
counts<-finger.countdata[,rownames(cluster.cor)]
identical(colnames(counts), rownames(cluster.cor))

# Filter COUNTS PER MILLION
y <- DGEList(counts = counts, group = as.factor(cluster.type))
#y <- DGEList(counts = counts.finger, group = as.factor(expresion.groups))
cpm <- cpm(y)
lcpm <- cpm(y, log=TRUE)




# 16. NORMALIZATION, CPM FILTERING AND HEATMAPS
# ==============================================================


# 2. UNFILTERED DATA: Heapmap expression with the top 100 more variables

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)


# We estimate the variance for each row in the logcounts matrix
finger.logcounts<-counts[rownames(y$counts)%in%ics,]   
dim(finger.logcounts)
# Variance between genes          
var_genes <- apply(finger.logcounts, 1, var)
# Get the gene names for the top 25 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:100]
highly_variable_cpm <- finger.logcounts[select_var,]
dim(highly_variable_cpm)


source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
clab<-as.matrix(cbind(col.cell, col.cell.clin))
colnames(clab)=c("Clinical subtype","Cluster subtype")
heatmap.3(highly_variable_cpm,col=rev(morecols(50)), trace="none", main="ICs Expression",ColSideColorsSize=2,ColSideColors=clab,scale="row", Colv=F, Rowv=T,  margins=c(6,12))
legend("topright",legend=c("Papillary","Not-papillary","Cluster 1","Cluster 2","Cluster 3"),
       fill=c("purple","orange", "black", "red", "green"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)



# FILTERING STEP: 
# Removing the low expressed genes

# Filtering step
keep.exprs <- rowSums(lcpm>6)>=98
y <- y [keep.exprs,]
dim(y)

nsamples <- ncol(y[,1:10])
#samplenames <- substring(colnames(y[,1:10]), 12, nchar(colnames(y[,1:10])))
samplenames <- colnames(y[,1:10])
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(y, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")


heatmap(cor(lcpm))



# 2. FILTERED DATA: Heapmap expression with the top 100 more variables

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)


# We estimate the variance for each row in the logcounts matrix
finger.logcounts<-y$counts   
dim(finger.logcounts)
# Variance between genes          
var_genes <- apply(finger.logcounts, 1, var)
# Get the gene names for the top 25 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:100]
highly_variable_cpm <- finger.logcounts[select_var,]
dim(highly_variable_cpm)


source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
clab<-as.matrix(cbind(col.cell, col.cell.clin))
colnames(clab)=c("Cluster subtype", "Clinical subtype")
heatmap.3(highly_variable_cpm,col=rev(morecols(50)), trace="none", main="ICs Expression",ColSideColorsSize=2,ColSideColors=clab,scale="row", Colv=F, Rowv=T,  margins=c(6,12))
legend("topright",legend=c("Papillary","Not-papillary","Cluster 1","Cluster 2","Cluster 3"),
       fill=c("purple","orange", "black", "red", "green"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)


# Heatmap with the 80% most variables

highly_variable_cpm<-finger.logcounts[which(var_genes> quantile(var_genes, c(.8))),]
dim(highly_variable_cpm)


source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
clab<-as.matrix(cbind(col.cell, col.cell.clin))
colnames(clab)=c("Cluster subtype", "Clinical subtype")
heatmap.3(highly_variable_cpm,col=rev(morecols(50)), trace="none", main="ICs Expression",ColSideColorsSize=2,ColSideColors=clab,scale="row", Colv=F, Rowv=T,  margins=c(6,12))
legend("topright",legend=c("Papillary","Not-papillary","Cluster 1","Cluster 2","Cluster 3"),
       fill=c("purple","orange", "black", "red", "green"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)



# Other heatmap representation

library(pheatmap)
mat<-highly_variable_cpm
mat  <- mat - rowMeans(mat)
pheatmap(mat, cluster_cols=F)



# HCLUSTERING AND DENDROGRAMS ON HIGH VARIABLE CPM FILTERED 


d.exp<-dist(t(highly_variable_cpm))
hc.exp<-hclust(d.exp, method="ward.D")
plot(hc.exp, hang = -1)
# Determine the optimal clusts
k.max <- 15 # Maximal number of clusters
fviz_nbclust(highly_variable_cpm, hcut, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)
# Expression Hkmeans for ICs
hkm.exp.hvar<-hkmeans(t(as.matrix(highly_variable_cpm)),3)
hkmeans_tree(hkm.exp.hvar)
dend_exp.hvar<-as.dendrogram(hkm.exp.hvar$hclust)
lbls<-c("1", "2", "3")
dend_exp.hvar %>% 
  color_branches(k=3, groupLabels =lbls) %>% 
  set("labels_cex", .5) %>% 
  plot(horiz=T) 



#Compare dendograms 
dl.EXP.hvar.CLUST <- dendlist(highlight_branches(mca.dend), highlight_branches(dend_exp.hvar))
tanglegram(dl.EXP.hvar.CLUST, sort = T, common_subtrees_color_lines = T, highlight_distinct_edges  = F, highlight_branches_lwd = F)

tanglegram(mca.dend,dend_exp.hvar, color_lines=(cluster_var$call$X$clust),
           common_subtrees_color_lines = T, highlight_distinct_edges  = F, highlight_branches_lwd = F)








































# 16. EXPRESSION ANALYSIS WITH LIMMA
# ==============================================================


exp.matrix<-exp.matrix[, rownames(cluster.cor[order(cluster.cor$cluster.group.cluster),])]
identical(colnames(exp.matrix), rownames(cluster.cor))
exp.matrix<-exp.matrix[complete.cases(exp.matrix), ]
dim(exp.matrix)


readES = ExpressionSet(log2(data.matrix(exp.matrix)+1))
pData(readES) = cluster.group
group= readES$cluster

design <- model.matrix(~group)
norm.factor <- calcNormFactors(data.matrix(exp.matrix))
dat.voomed <- voom(data.matrix(exp.matrix), design, plot = TRUE, lib.size = colSums(data.matrix(exp.matrix)) * norm.factor)
fit <- lmFit(dat.voomed, design)
fit <- eBayes(fit)
results <- decideTests(fit)
summary(results)
vennDiagram(results)
# Unadjusted
voom.limma.genes <- rownames(topTable(fit, p.value=0.05))


# To plot the expression of any set of genes in a heatmap.
mycol <- colorpanel(1000,"blue","white","red")

# To plot significant DEG over groups in a heatmap
selected <- p.adjust(fit$p.value[, 2], method = "fdr") <0.05 
sum(selected)
esetSel <- readES[selected, ] 
mycol <- colorpanel(1000,"blue","white","red")
i <- which(rownames(dat.voomed$E) %in% rownames(exprs(esetSel)))



heatmap.3(dat.voomed$E[i,],col=rev(morecols(50)), trace="none", main="DEG Expression",ColSideColorsSize=2,ColSideColors=clab,scale="row", Colv=F, Rowv=T,  margins=c(6,12))
legend("topright",legend=c("Papillary","Not-papillary","Cluster 1","Cluster 2","Cluster 3"),
       fill=c("purple","orange", "black", "red", "green"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)


sel.vommed<-rownames(dat.voomed$E[i,])
sel.vommed<-as.vector(sapply(sel.vommed, function(x) 
  unlist(strsplit(x, split = "|", fixed = TRUE))[1]))
cat(sel.vommed, sep="\n")
length(sel.vommed)


save(sel.vommed, file = "selvoomed.RData")
save(ics2.deleted, file = "ics2deleted.RData")
save(ics3.duplicated, file = "ics3duplicated.RData")

# 16. PLOTMEANS AND STATISTISTICS
# ==============================================================
# Plotmeans and statistics

pdf("plotmeans-voomed-summary2.pdf")
counts_t<-as.data.frame(t(data.matrix(exp.matrix)))
dim(counts_t)
counts_t$group<-cluster.group$cluster
relevant<-c()
for (gene in sel.vommed){
  gen.counts<-as.numeric(counts_t[, grep(as.character(gene),colnames(counts_t))[1]])
  aov(gen.counts~as.factor(groups))->fit
  summary(fit)
  anova.pval<-summary(fit)[[1]][["Pr(>F)"]][1]
  if (anova.pval<1) {
    par(mfrow=c(2,2))
    title = paste("Expression of ",as.character(gene)," by group")
    plot( gen.counts ~ group, data=counts_t, border="blue", col="cyan",
          main=title )
    plotmeans( gen.counts ~ groups, data=counts_t, barwidth=2, connect=FALSE,
               main="Means and 95% Confidence Intervals\nof expression by group")
    info <- sapply( split(gen.counts, counts_t$group),
                    function(x) round(c(Mean=mean(x), SD=sd(x), N=gdata::nobs(x)),2) )
    textplot( info, valign="top"  )
    title("Expression by Group")
    stat <- TukeyHSD(fit)
    textplot( capture.output(stat), valign="top")
    title("TukeyHSD Expression by Group")
    relevant<-c(relevant, gene)
  }
}

dev.off()

relevant
relevant<-as.vector(sapply(relevant, function(x) 
  unlist(strsplit(x, split = "|", fixed = TRUE))[1]))
cat(relevant, sep="\n")

table(matrix_clusters$NR2E1, matrix_clusters$cluster)
round(prop.table(table(matrix_clusters$CTDSP1, matrix_clusters$cluster)[c(1),]),2)*100
round(prop.table(table(matrix_clusters$CTDSP1, matrix_clusters$cluster)[c(6),]),2)*100
chisq.test(table(matrix_clusters$ATOH7, matrix_clusters$cluster)[c(1,7),])


table(matrix_clusters$CALCB, matrix_clusters$cluster)
cat(genes.1U23, sep="\n")




# GROUP COMPARISION with CONTRAST MATRIX

design<-model.matrix(~0+group)
colnames(design)<-c("c1","c2","c3")
fit<-lmFit(readES, design)
contrast.matrix<-makeContrasts(
  c1vsc2 = c1 - c2, 
  c1vsc3 = c1 - c3, 
  c2vsc3 = c2 - c3, 
  levels = colnames(design))
contrast.matrix

fit2  <- contrasts.fit(fit, contrast.matrix)
fit2  <- eBayes(fit2)

topTable(fit2,adjust="fdr")
results <- decideTests(fit2)
summary(results)
vennDiagram(results)



selected.c1<- p.adjust(fit2$p.value[,1], method = "fdr") <0.05 
esetSel.c1 <- readES[selected.c1, ] 
mycol <- colorpanel(1000,"blue","white","red")
i <- which(rownames(dat.voomed$E) %in% rownames(exprs(esetSel.c1)))
col.cell <- c("purple","orange", "green")[cluster.type]
heatmap.2(dat.voomed$E[i,], scale="row",
          labRow=rownames(exprs(esetSel.c1))[i], labCol=readES$group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column", ColSideColors = col.cell)

heatmap(exprs(esetSel.c1), col=topo.colors(100), ColSideColors=col.cell) 






# 9. SURVIVAL DATA FOR EXPRESSION CLUSTERS AND ANALYSIS
# ==============================================================
setwd("~/omics/code to TFM/data/Clinical data")
surv.clin<-survivalTCGA(BLCA.clinical, extract.cols= "admin.disease_code")
surv.clin<-load("~/omics/code to TFM/data/Clinical data/surv.clin.RData")
surv.clin$bcr_patient_barcode<-substr(surv.clin$bcr_patient_barcode, 6,12)
rownames(surv.clin)<-surv.clin$bcr_patient_barcode
index<-intersect(rownames(surv.clin), names(expresion.group))
surv.clin.sel<-surv.clin[index,]
dim(surv.clin.sel)
exp.index<-expresion.group[index]

surv.data<-data.frame(surv.clin.sel[,c(1,3)], exp.index)
names(surv.data)<-c("times", "status", "cluster")
head(surv.data)
prop.table(table(surv.data$cluster, surv.data$status),1)


coxph(Surv(times, status)~cluster, data=surv.data)
sfit <- survfit(Surv(times, status)~cluster, data=surv.data)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table = F)


index<-intersect(rownames(surv.clin), tcgaANDme)
index
surv.clin.sel<-surv.clin[index,]
dim(surv.clin.sel)
tcga_blca[index,2]

surv.data<-data.frame(surv.clin.sel[,c(1,3)], tcga_blca[index,2])
names(surv.data)<-c("times", "status", "cluster")
head(surv.data)
prop.table(table(surv.data$cluster, surv.data$status),1)


coxph(Surv(times, status)~cluster, data=surv.data)
sfit <- survfit(Surv(times, status)~cluster, data=surv.data)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table = TRUE)







class(surv.data$times)
fit <- survfit(Surv(times, status) ~ 1, data = surv.data)
# Drawing curves
ggsurvplot(fit, color = "#2E9FDF")

fit<- survfit(Surv(times, status) ~ cluster , data = surv.data)
# Drawing survival curves
ggsurvplot(fit)











