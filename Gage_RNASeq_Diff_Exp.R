# This code performs post alignment QC and differential gene expression (DGE) on RNASeq data.
# QC is performed by checking the
# 1-  Correlations of most expressed genes(5%)
# 2-3-  PCA and Clustering and if they are very similar
# 4- Also check number of non zero genes and number of reads aligned per sample.
# 5-MA plot and pvalue bar plots
# Differential gene expression
# 1- Performed by DESeq2

rm(list=ls())

# load libraries
library('DESeq2')
library('BiocManager')
library('limma')
library('dendextend')
library('biomaRt')

#Path to count files
path.count.files <- "/Users/muhkas/Desktop/HN/Evren/GagE_Project/Sequencing/RNASeq/countedfiles"
#Read in the counted read numbers
count.files.names <- list.files(path.count.files, pattern="*Aligned_Count$", full.names=TRUE)
#length(count.files.names)

# Extract only the relevant data from count files
# that is ensemble gene id and counts
sub_count_data<-function(x){
  return (read.table(x)[, c(1,7)])
}
count.data <- lapply(count.files.names, sub_count_data)

# Perform basic checkes
head(count.data[[2]])
summary(count.data) 
head(unlist(count.data))

#Convert list to matrix(required for dseq2)
countdatamatrix <-""
countdatamatrix <- matrix(unlist(count.data), 2 * length(count.files.names), nrow = nrow(count.data[[1]]))
head(countdatamatrix)

#Store the row names, that is ensemble gene ids.
my.rownames <- countdatamatrix[,1]

#Customized col names, that is removing the redundant parts of the file names.
count.files.name.custom <- list.files(path.count.files, pattern="*Aligned_Count$", full.names=FALSE)
custom_file_names<- function(x){
  return(unlist(strsplit(x,"_BEA"))[1] )
}

count.files.name.custom.edit  <- ""
count.files.name.custom.edit  <- unlist(lapply(count.files.name.custom, custom_file_names))


# Subset the data and only keep the columns with counts only
col.keep <- seq(from = 2, to= 2*length(count.files.name.custom.edit), by=2 )
countdatamatrix.cln <- countdatamatrix[(2:nrow(count.data[[1]])), col.keep ]
head(countdatamatrix.cln)
# Data inside the matrix is of class character. 
# Convert all data into class numeric and then remade the matrix
class(countdatamatrix.cln[1,1])
my.nrow <- nrow(countdatamatrix.cln) 
my.ncol <- ncol(countdatamatrix.cln)
countdatamatrix.qt.fr   <-  matrix( as.numeric(countdatamatrix.cln), nrow = my.nrow, ncol = my.ncol)
class(countdatamatrix.cln[1,1])

#Add col names  and row names
colnames(countdatamatrix.qt.fr) <- count.files.name.custom.edit
row.names(countdatamatrix.qt.fr) <-  my.rownames[(2:nrow(count.data[[1]])) ]
head(countdatamatrix.qt.fr)

#Arrange samples according to experimental groups
YT.group    <- c(1,11, 12)
YTS.group   <- c(13,14,15) 
KMYG1.group <- c(2, 3, 4)
NK92.group  <- c(5,6,7)
NKL.group   <- c(8,9,10)
countdatamatrix.qt.fr<- countdatamatrix.qt.fr[, c(YT.group, YTS.group, KMYG1.group,
                                                  NK92.group, NKL.group)]
# Remove low frequency genes,
# By only keeping data for the genes which have >15 reads in all samples.
countdatamatrix.qt.fr<- countdatamatrix.qt.fr[ which(rowSums(countdatamatrix.qt.fr) >15),]

#Start of diagnostics  quality controls
# Total number of reads counted per sample
  my.col.sums <- colSums(countdatamatrix.qt.fr)
  max(my.col.sums)
  min(my.col.sums)
  mean(my.col.sums)
  print( c("Number of reads aligned per sample=", my.col.sums))
  ###Number of non zero genes
  idx.nz <- apply(countdatamatrix.qt.fr, 1, function(x) { all(x > 0)})
  sum(idx.nz)
  head(countdatamatrix.qt.fr)
  #Start of clustering of data
  help(dist)
  # Calculate the distance matrix using euclidean method
  # Other methods are: 'manhattan', 'canberra', 'minkowski' etc.
  dist.euc <- dist(t(countdatamatrix.qt.fr), method = 'euclidean') 
  #Perform hierarchical clustering
  hr.clust<- hclust(dist.euc )
  dend <- as.dendrogram(hr.clust)
  #Color each cell line
  my-colors <- c('red', 'green', 'blue' , 'orange', 'yellow3')
  labels_colors(dend) <- rep(my-colors, each=3)[order.dendrogram(dend)]
  par(mar=c(15,4.1,2,2))
  plot(dend, xlab="", ylab="Distance", main="Clustering of data ")

#########################Add Metadata to create a dds objec of DSEq2
colDatamat <-0
colDatamat <- matrix(nrow = ncol( countdatamatrix.qt.fr) ,   ncol =2  )

dim(colDatamat)

colDatamat[,1]<- c( rep("YT", length(YT.group) ),  rep("YTS", length( YTS.group ) ),
                    rep("KMYG1", length(KMYG1.group)),  rep("NK92",length(NK92.group )),
                    rep("NKL",length(NKL.group ))
                  )
#Group by UT=non-transfected, Sol1=transfection method1 , Sol2= transfection method 2.

indeces.group <- lapply( as.matrix(c('*UT', '*Sol1', '*Sol2')),  function(x) {
  grep( x, colnames(countdatamatrix.qt.fr))} )

colDatamat[c(indeces.group[[1]], indeces.group[[2]], indeces.group[[3]]), 2] <-
  rep(c('UT', 'Sol1', 'Sol2'),  each=length(indeces.group[[1]]) )

rownames(colDatamat) <- colnames(countdatamatrix.qt.fr)  
colnames(colDatamat) <- c( "CellLine", "Trans_Group")
colDatamat.frm <- data.frame(colDatamat)
#Change the relevel to any control to perform that type of analysis.
colDatamat.frm$Trans_Group <- as.factor(colDatamat.frm$Trans_Group)
dds <- DESeqDataSetFromMatrix( countData=countdatamatrix.qt.fr,
                               colData= colDatamat.frm,
                               design= ~Trans_Group)
dds$Trans_Group <- relevel(dds$Trans_Group,"UT")

dds <- DESeq(dds)
summary(dds)
#Diagnostic plots
plotMA( dds, ylim = c(-5, 5) )
plotDispEsts( dds, ylim = c(1e-6, 1e1) )
hist( ut_sol1$pvalue, breaks=20, col="grey" )

resultsNames(dds)
ut_sol1 <- results(dds, contrast = c("Trans_Group", "Sol1", "UT" ), independentFiltering=FALSE  )
ut_sol2 <- results(dds, contrast = c("Trans_Group", "Sol2", "UT" ), independentFiltering=FALSE  )

#Threshold for differentially expressed gene signature
ut_sol1.sig<- ut_sol1[   which(ut_sol1$padj < 0.05),  ]
ut_sol2.sig <-  ut_sol2[   which(ut_sol2$padj < 0.05),  ]

# View summary of results
summary(ut_sol1)
summary(ut_sol2)

#Add gene names on the basis of IDs
ensembl<- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#listFilters(mart) 
#listAttributes(mart)
#listDatasets(ensembl)
dds.fl.rst.cal <-""
dds.fl.rst.cal <-     ut_sol2.sig
dds.fl.rst.cal.fr <-    data.frame( dds.fl.rst.cal )

genes <- rownames(dds.fl.rst.cal.fr)  # length(genes)
G_list <- ""
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description"), values=genes, mart= mart)


dds.fl.rst.cal.fr$ensid <- rownames(dds.fl.rst.cal.fr)

idx <- match(dds.fl.rst.cal.fr$ensid, G_list$ensembl_gene_id)  # Match ids

dds.fl.rst.cal.fr$matchedEnsID  <-  G_list$ensembl_gene_id[idx]
dds.fl.rst.cal.fr$Symbol  <-  G_list$hgnc_symbol[idx]   # Add gene symbols (names)
dds.fl.rst.cal.fr$Description  <-  G_list$description[idx] #Add gene description

head(dds.fl.rst.cal.fr)  # see first 5 rows of the results

summary(dds.fl.rst.cal.fr ) #see the summary of result
write.csv ( dds.fl.rst.cal.fr ,paste(path.count.files,"/DGE_UTvsSol2.csv" , sep="") , row.names = TRUE )





##############################
rld <- rlogTransformation(dds, blind=TRUE)
head(rld)


groupnames <- c( rep("cntrl", length(cntrl.ind.d1) + length(cntrl.ind.d2)  ),  rep("bta",length(bta.ind.d1 )+ length(bta.ind.d2 ) ),
                 rep("can", length(can.ind.d1)+ length(can.ind.d2) ),  rep("reti",length(reti.ind.d1 )+ length(reti.ind.d2 )),
                 rep("vc",length(vc.ind.d1 )+ length(vc.ind.d2 ))
)

groupnames <- c( rep("b1",20),rep("b2",20) )



colnames(rld) <- make.unique(groupnames)
#colorcodes<- c( cntrl="black", bta="blue", can="red", reti="yellow", vc="green")
colorcodes<- c( b1="black", b2="blue")
distsRL <- dist(t(assay(rld)))
hc<- hclust(distsRL )
dend <- as.dendrogram(hc)
#install.packages("dendextend")
#library(dendextend)
labels_colors(dend) <- colorcodes[groupnames][order.dendrogram(dend)]
par(mar=c(8,4.1,2,2))
plot(dend, xlab="Groups", ylab="Distance", main="Clustering of  Kashif data \n Before outlier removal",
     cex.axis=1.5, cex.lab=1.5, las=1 )
#pca

ntop = 500
Pvars <- rowVars(assay(rld))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                      length(Pvars)))]
PCA <- prcomp(t(assay(rld)[select, ]), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
dataGG = data.frame(PC1 = PCA$x[,1],PC2 = PCA$x[,2],
                    PC3 = PCA$x[,3],PC4 = PCA$x[,4],
                    condition = colData(dds))
summary(PCA)

(qplot(PC3, PC4, data = dataGG,  
       main = "PC3 vs PC4, top variable genes \n Before outlier removal of Kashif data", size = I(6))
  + labs(x = paste0("PC3, VarExp:", round(percentVar[3],4)),
         y = paste0("PC4, VarExp:", round(percentVar[4],4)))
  + scale_colour_brewer(type="qual", palette=2)
)

###End of diagnostics



#pca
library(ggplot2)
ntop = 500
Pvars <- rowVars(assay(rld))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                      length(Pvars)))]
PCA <- prcomp(t(assay(rld)[select, ]), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
dataGG = data.frame(PC1 = PCA$x[,1],PC2 = PCA$x[,2],
                    PC3 = PCA$x[,3],PC4 = PCA$x[,4],
                    condition = colData(dds)$MelanAntiOx)
summary(PCA)

(qplot(PC3, PC4, data = dataGG, color =  condition, 
       main = "PC3 vs PC4, top variable genes \n Before outlier removal of Kashif data", size = I(6))
  + labs(x = paste0("PC3, VarExp:", round(percentVar[3],4)),
         y = paste0("PC4, VarExp:", round(percentVar[4],4)))
  + scale_colour_brewer(type="qual", palette=2)
)

####################################################3





