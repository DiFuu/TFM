###############
#Preprocessing#
###############

#remove the saved variables
rm(list=ls())

#Read the data and pass it to a .txt

library(readxl)

file <- read_excel("1-Datos.xlsx", 1)
write.table(file, "1-Datos.txt",sep="\t")

#Experiment Design
file2 <- read_excel("2-Muestras.xlsx", 1)
write.table(file2, "2-Muestras.txt",sep="\t")

#Dataset
seqdata0 <- read.delim("1-Datos.txt", stringsAsFactors = FALSE, row.names = 1)

head(seqdata0)
dim(seqdata0)

#Eliminate the columns from 1 to 110
seqdata <- seqdata0[,-(1:110)]

#Experimental Design
end_layout <- read.delim("2-Muestras.txt")

#We change the name of the samples, starting from end_layout

colnames(end_layout) #In this case, we are interested in well.Sample.Cells_type.Day2 and Library
type <- substr(end_layout$Cells_type,1,1)

new_name <- paste(end_layout$well,end_layout$Sample,type,end_layout$Day2,end_layout$Library,sep=".")

end_layout$name <- new_name
colnames(seqdata) <- new_name

#The name of the genes in seqdata matches the data in the name column in endometrium layout
table(colnames(seqdata)==end_layout$name)

#We are going to filter out genes with very low expression, those that have less than 0.5 counts per million reads (CPM)
#in half of the samples. To do this, we use the "edgeR" package that provides the cpm function.

library(limma)
library(edgeR)

#expression counts
mycounts <- seqdata
dim(mycounts)

mycounts[1:5,1:3]

#Filter
#We maintain genes with at least 0.5 counts per million reads (cpm) in at least 2 samples
isexpr <- rowSums(cpm(mycounts)>0.5) >= 2
table(isexpr)

mycounts <- mycounts[isexpr,]
genes <- rownames(mycounts)

dim(mycounts)

#In our case, we want to do two different analyses, one for the stromal cells and one for the epithelial cells.
#Both analyzes will be completely independent because what we want to compare is the effect of time, not the tissue itself.

#Instead of renaming all the input data and making it double, the line of the epithelium or the stroma will be executed depending on which one is convenient

#Here we only have to switch between Stroma or Epithelium depending on the data we want to analyze

library(dplyr)

end_layout <- end_layout %>%
  filter(Cells_type == "Epithelium")
mycounts <- mycounts[,end_layout$name]
dim(mycounts)

#For normalization we use the voom function from the limma package which normalizes the read counts
#and applies a linear model before statistical analysis of differential expression

#condition group
group<-factor(end_layout$Day2)

#matrix design for limma
design <- model.matrix(~0+group)
#we substitute "group" for the layout column names
colnames(design)<- gsub("group","",colnames(design))
colnames(design) <- c("LH2","LH8")
design
#write.table(design[,1],file="Str_class.txt",row.names=F,quote=F,sep="\t")
write.table(design[,1],file="Ep_class.txt",row.names=F,quote=F,sep="\t")

#normalization factors between libraries
nf <- calcNormFactors(mycounts)

#we normalize the read counts with the 'voom' function
y <- voom(mycounts,design,lib.size=colSums(mycounts)*nf)
#we extract the normalized read counts
counts.voom <- y$E

#we save the normalized expression data in the output directory
#write.table(counts.voom,file="Str_counts.voom.txt",row.names=F,quote=F,sep="\t")
#write.table(mycounts,file="Str_mycounts.txt",row.names=F,quote=F,sep="\t")

write.table(counts.voom,file="Ep_counts.voom.txt",row.names=F,quote=F,sep="\t")
write.table(mycounts,file="Ep_mycounts.txt",row.names=F,quote=F,sep="\t")

#write.table(counts.voom,file="Ep_counts.voom.txt",row.names=F,quote=F,sep="\t")
#write.table(mycounts,file="Ep_mycounts.txt",row.names=F,quote=F,sep="\t")

#we fit the linear model for each gene
fit <- lmFit(y,design)

#we build the contrast matrix
cont.matrix <- makeContrasts(LH8-LH2,levels=design)
cont.matrix 

#we calculate estimated coefficients and standard errors for a given set of tests
fit <- contrasts.fit(fit, cont.matrix)

#moderated t-statistics of differential expression using empirical Bayesian moderation of standard errors
fit <- eBayes(fit)
options(digits=3)

#output fit
dim(fit)

#The topTable function summarizes the output of limma in a table format.
#Genes with a pvalue < 0.05 are selected, we don't filter by logFC.

#set pvalue adjusted threshold and logarithmic change threshold
mypval=0.05
#myfc=2

#name of the coefficient for the comparison of interest
colnames(fit$coefficients)

mycoef="LH8 - LH2"
#output table of the 10 most important genes for this comparison
topTable(fit,coef=mycoef)

#full table ("n = number of genes in the fit")
limma.res <- topTable(fit,coef=mycoef,n=dim(fit)[1])

#only significant genes (adjusted p-value < mypval)
limma.res.pval <- topTable(fit,coef=mycoef,n=dim(fit)[1],p.val=mypval)
dim(limma.res.pval)

#Significant genes with low adjusted p-value and high FC
#limma.res.pval.FC <- limma.res.pval[which(abs(limma.res.pval$logFC)>myfc),]

#to order the genes based on their number
limma.res.pval$Gene <- as.numeric(rownames(limma.res.pval))
limma.res.pval <- limma.res.pval[order(limma.res.pval$Gene),]

#we write the limma output table for significant genes to a tab-delimited file
#filename = paste("Genes_STROMA_",mycoef,"_pval",mypval,".txt",sep="")
#write.table(limma.res.pval,file=filename,row.names=F,quote=F,sep="\t")

filename = paste("Genes_EPITHELIUM_",mycoef,"_pval",mypval,".txt",sep="")
write.table(limma.res.pval,file=filename,row.names=F,quote=F,sep="\t")

#To output the normalized data to the data_normalized.txt file

#We need the name of the genes
dim(seqdata0)
seqdata1 <- seqdata0[,-(2:195)]

#we obtain the position of the genes that we are going to export as a result of normalization
list_genes2 <- as.integer(rownames(counts.voom))
length(list_genes2)

#We get the corresponding genes
genesnorm <- NULL
for (ngene in list_genes2){
  genesnorm[[length(genesnorm) + 1]] <- (seqdata1[(ngene)])
}

genesnorm <- as.character(genesnorm)

#We add the genes
normdata <- counts.voom #copiamos los datos
normdata <- as.data.frame(normdata) #transformamos a dataframe
normdata$Gene <- genesnorm #añadimos el nombre de los genes
normdata <- normdata %>%
  dplyr::select("Gene", everything()) #ordenamos para que la columna gen esté la primera

#We save the files
#write.table(normdata,file="Str_data_normalized.txt",row.names=F,quote=F,sep="\t")
#write.table(normdata$Gene,file="Str_id_normalized.txt",row.names=F,quote=F,sep="\t")

write.table(normdata,file="Ep_data_normalized.txt",row.names=F,quote=F,sep="\t")
write.table(normdata$Gene,file="Ep_id_normalized.txt",row.names=F,quote=F,sep="\t")


#We take a file with the final data to be able to continue with the analysis

#we get the position of the genes that we are going to export
list_genes2 <- as.integer(rownames(limma.res.pval))
length(list_genes2)

#We get the corresponding genes
genesfc <- NULL
for (ngene in list_genes2){
  genesfc[[length(genesfc) + 1]] <- (seqdata1[(ngene)])
}

genesfc <- as.character(genesfc)

#We add the genes
fcdata <- limma.res.pval #copiamos los datos
fcdata <- as.data.frame(fcdata) #transformamos a dataframe
fcdata$Gene <- genesfc #añadimos el nombre de los genes
fcdata <- fcdata %>%
  dplyr::select("Gene", everything()) #ordenamos para que la columna gen esté la primera

#We save the files
#write.table(fcdata,file="data_Str.txt",row.names=F,quote=F,sep="\t")
#write.table(fcdata$Gene,file="id_Str.txt",row.names=F,quote=F,sep="\t")

write.table(fcdata,file="data_Ep.txt",row.names=F,quote=F,sep="\t")
write.table(fcdata$Gene,file="id_Ep.txt",row.names=F,quote=F,sep="\t")

#Heatmap for general data visualization

library(DESeq2)

#we read the file
genedata <- read.table("1-Datos.txt",header=T, row.names = "Gene")
attach(genedata)

#we select the data of interest
genedata <- genedata[,111:195]
geneNames<-row.names(genedata)

#we transform the data frame to an array object
is.matrix(genedata)
data.matrix <- as.matrix(genedata)
is.matrix(data.matrix)
data.matrix[1:5,]
summary(data.matrix)

countdata <- data.matrix
head(countdata, 3)

#we open the file with the sample data and select the columns of interest
coldata<- read.table("2-Muestras.txt",header=T, row.names = "name")
attach(coldata)

#We change the name of the columns so that they are exactly the same as those of the samples

colnames(coldata) #we are interested in well.Sample.Cells_type.Day2 and Library
type <- substr(coldata$Cells_type,1,1)

name <- paste(coldata$well,coldata$Sample,type,coldata$Day2,coldata$Library,sep=".")

rownames(coldata) <- name
colnames(countdata) <- name

#Now we take only the columns that interest us for the heatmap
coldata <- coldata[,c("Cells_type","Day2")]
coldata[1:5,]
colnames(coldata) <- c("Tissue","Day")

#We separate the data according to whether they are stromal or epithelial cells
library(dplyr)
coldatast <- coldata %>%
  filter(Tissue == "Stroma")
countdatast <- countdata[,rownames(coldatast)]

coldataep <- coldata %>%
  filter(Tissue == "Epithelium")
countdataep <- countdata[,rownames(coldataep)]

#Set str + ep
ddsMat <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ Day + Tissue + Day:Tissue)

#Only str
ddsMat <- DESeqDataSetFromMatrix(countData = countdatast, colData = coldatast, design = ~ Day)

#Only ep
ddsMat <- DESeqDataSetFromMatrix(countData = countdataep, colData = coldataep, design = ~ Day)

#A warning message appears "it is better not to put symbols"
nrow(ddsMat)

#Variance stabilization for data visualization
library("vsn")
vsd <- vst(ddsMat, blind = FALSE)
head(assay(vsd), 3)

colData(vsd)

#Visualization
library("genefilter")
library("pheatmap")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 40)
topVarGenes
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Day","Tissue")])
pheatmap(mat, annotation_col = anno)


#######################
#Differential analysis#
#######################

#remove the saved variables
rm(list=ls())

require(multtest)
library(RankProd)

#We read the file with the data to analyze
#file<- read.table("Str_data_normalized.txt",header=T, row.names=1)
file<- read.table("Ep_data_normalized.txt",header=T, row.names=1)

attach(file)

geneNames<-row.names(file)

#Convert the object to an array
is.matrix(file)
data.matrix <- as.matrix(file)
is.matrix(data.matrix)
data.matrix[1:5,]

#Class
#LH2 1, LH8 0
#data.cl <- read.delim("Str_class.txt", stringsAsFactors = FALSE)
data.cl <- read.delim("Ep_class.txt", stringsAsFactors = FALSE)

#we convert it to numeric
data.cl <- as.numeric(data.cl$x)

#Statistical test
#Function mt.teststat(X,classlabel,test="t",nonpara="n")
#test="t", tests are based on two-sample Welch's t-statistics (unequal variances)
#nonpara="n", original data is used

teststat<-mt.teststat(file, data.cl, test="t",nonpara="n")

df=42 #41-2 to STR and 44-2 to EP
rawp0<-2*(1-pt(abs(teststat), df))
res.multtest <- mt.rawp2adjp(rawp0, "BH")
res.multtest <-res.multtest$adjp[order(res.multtest$index),]
res.multtest <-cbind(geneNames,res.multtest)

#write.table(res.multtest, "Str_multtest.txt")
write.table(res.multtest, "Ep_multtest.txt")

#select significants
#multtestd <- read.table("Str_multtest.txt",header=T, row.names=1)
multtestd <- read.table("Ep_multtest.txt",header=T, row.names=1)

attach(multtestd)
names(multtestd)
length.BH.significant <-length(which(multtestd[,3]<0.05))
length.BH.significant
BH.significant <-subset(multtestd, multtestd[,3]<0.05)
dim(BH.significant)

#write.table(BH.significant, "Str_BHsignificant.txt",row.names=F)
write.table(BH.significant, "Ep_BHsignificant.txt",row.names=F)


#We join the table of previous results with this new data

#Read the file
#filefc<- read.table("data_Str.txt",header=T, row.names=1)
filefc<- read.table("data_Ep.txt",header=T, row.names=1)

attach(filefc)

#Output gene names and save in a new column
geneNamesFC <- row.names(filefc)
filefc$geneNames <- geneNamesFC

#We join the two dataframes based on the gene name in common
data_FC <- merge(x = filefc, y = BH.significant, by = c("geneNames"))

library(data.table)

setnames(data_FC, "geneNames", "Gene")

#Save the data
#write.table(data_FC, "Data_Str_BH.txt",row.names=F)
write.table(data_FC, "Data_Ep_BH.txt",row.names=F)


####################
#GENOMIC ANNOTATION#
####################

#remove the saved variables
rm(list=ls())

#Annotation of EntrezGene IDs from RNAseq data can be performed using the BioMart database
#containing many species, such as humans, mice, zebrafish, chickens, and rats

library(biomaRt)
library(org.Hs.eg.db)

#We get the list of genes
#data_FC <- read.table("Data_Str_BH.txt",header=T)
data_FC <- read.table("Data_Ep_BH.txt",header=T)

attach(data_FC)

dim(data_FC)
list_genes <- data_FC[,1]

mart<- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org"))

#we use BioMart database to get gene symbol and description of these genes
detags.IDs <- getBM(
  filters= "hgnc_symbol",
  attributes= c("hgnc_symbol","description"),
  values= as.character(list_genes),
  mart= mart
)

dim(detags.IDs)

head(detags.IDs)

#We select only the names of interest
rownames(detags.IDs)<-detags.IDs$hgnc_symbol
entrez_genes.annot <- detags.IDs[as.character(list_genes),]
description <- entrez_genes.annot$description
data_FC.annot <- cbind(data_FC,description)

#We save the annotated table in a file
#write.table(data_FC.annot,file="Genes_STROMA_analyzed.txt",sep="\t",row.names=F)
write.table(data_FC.annot,file="Genes_EPITHELIUM_analyzed.txt",sep="\t",row.names=F)

###############
#GO ENRICHMENT#
###############

#remove the saved variables
rm(list=ls())

#We read the file and get the name of the genes
#all <- read.table("Genes_STROMA_analyzed.txt", row.names=NULL, quote="\"", comment.char="", header=TRUE)
all <- read.table("Genes_EPITHELIUM_analyzed.txt", row.names=NULL, quote="\"", comment.char="", header=TRUE)

query <- all[,1]

#We export the data of the genes to be able to use the GOplot
datosgenes <- all[,c("Gene","logFC")]
colnames(datosgenes) <- c("ID","logFC")

#write.table(datosgenes,"David_genes_STR.txt",row.names=F,sep="\t")
write.table(datosgenes,"David_genes_Ep.txt",row.names=F,sep="\t")

#We convert from symbol to entrez
library(org.Hs.eg.db)

#To see the available identifiers
columns(org.Hs.eg.db)

#We use mapIds to get the Entrez IDs
query2 <- mapIds(org.Hs.eg.db, query, 'ENTREZID', 'SYMBOL')

#write.table(query2, "query_entrezid_Str.txt")
write.table(query2, "query_entrezid_Ep.txt")

library(rJava)
library(RDAVIDWebService)
library(clusterProfiler)

#DAVID BP
GOBP <- enrichDAVID(query2, minGSSize = 10, idType = "ENTREZ_GENE_ID", annotation = "GOTERM_BP_DIRECT", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, david.user = "diana.fuertes@cragenomica.es")
KEGG <- enrichDAVID(query2, minGSSize = 10, idType = "ENTREZ_GENE_ID", annotation = "KEGG_PATHWAY", david.user = "diana.fuertes@cragenomica.es")

#Representation
GOBP <- setReadable(GOBP, 'org.Hs.eg.db', 'ENTREZID')
KEGG <- setReadable(KEGG, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(GOBP, foldChange=query2)
barplot(KEGG)

#Data
dataBP <- data.frame(GOBP$ID, GOBP$Description, GOBP$GeneRatio, GOBP$BgRatio, GOBP$pvalue, GOBP$p.adjust, GOBP$qvalue, GOBP$geneID, GOBP$Count, GOBP@ontology)
dataBP$GOBP.geneID <- gsub("/",", ", dataBP$GOBP.geneID)
dataBP$GOBP.ontology <- gsub("GOTERM_BP_DIRECT","BP", dataBP$GOBP.ontology)
colnames(dataBP) <- c("ID","term","GeneRatio","BgRatio","pvalue","adj_pval","qvalue","genes","count","category")

#DAVID CC
GOCC <- enrichDAVID(query2, minGSSize = 10, idType = "ENTREZ_GENE_ID", annotation = "GOTERM_CC_DIRECT", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, david.user = "diana.fuertes@cragenomica.es")

#Representation
GOCC <- setReadable(GOCC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(GOCC, foldChange=query2)

#Data
dataCC <- data.frame(GOCC$ID, GOCC$Description, GOCC$GeneRatio, GOCC$BgRatio, GOCC$pvalue, GOCC$p.adjust, GOCC$qvalue, GOCC$geneID, GOCC$Count, GOCC@ontology)
dataCC$GOCC.geneID <- gsub("/",", ", dataCC$GOCC.geneID)
dataCC$GOCC.ontology <- gsub("GOTERM_CC_DIRECT","CC", dataCC$GOCC.ontology)

#DAVID MF
#THERE ISN'T FOR EPITHELIUM, SO DO NOT RUN
GOMF <- enrichDAVID(query2, minGSSize = 10, idType = "ENTREZ_GENE_ID", annotation = "GOTERM_MF_DIRECT", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, david.user = "diana.fuertes@cragenomica.es")

#Representation
GOMF <- setReadable(GOMF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(GOMF, foldChange=query2)

#Data
dataMF <- data.frame(GOMF$ID, GOMF$Description, GOMF$GeneRatio, GOMF$BgRatio, GOMF$pvalue, GOMF$p.adjust, GOMF$qvalue, GOMF$geneID, GOMF$Count, GOMF@ontology)
dataMF$GOMF.geneID <- gsub("/",", ", dataMF$GOMF.geneID)
dataMF$GOMF.ontology <- gsub("GOTERM_MF_DIRECT","MF", dataMF$GOMF.ontology)

#Export
#write.table(dataBP,"David_results_Str.txt",row.names=F,sep="\t")
#write.table(dataCC,"David_results_Str.txt",row.names=F,col.names=FALSE,sep="\t",append=TRUE)
#write.table(dataMF,"David_results_Str.txt",row.names=F,col.names=FALSE,sep="\t",append=TRUE)

write.table(dataBP,"David_results_Ep.txt",row.names=F,sep="\t")
write.table(dataCC,"David_results_Ep.txt",row.names=F,col.names=FALSE,sep="\t",append=TRUE)
#write.table(dataMF,"David_results_Ep.txt",row.names=F,col.names=FALSE,sep="\t",append=TRUE)

#write.table(KEGG, "KEGG_Str.txt",row.names=F,sep="\t")
write.table(KEGG, "KEGG_Ep.txt",row.names=F,sep="\t")

#dotplot
#Not used
library("gridExtra")
p1 <- dotplot(GOBP, title = "Biological Process")
p2 <- dotplot(GOCC, title = "Cellular Component")
p3 <- dotplot(GOMF, title = "Molecular Function")

grid.arrange(p1, p2, p3, nrow = 1)

######################
#GOPLOT VISUALIZATION#
######################

#remove the saved variables
rm(list=ls())

library(GOplot)

#Data
#DAVIDdata <- read.table("David_results_Str.txt",header=TRUE)
#datosgenes <- read.table("David_genes_Str.txt", header=TRUE)

DAVIDdata <- read.table("David_results_Ep.txt",header=TRUE)
datosgenes <- read.table("David_genes_Ep.txt", header=TRUE)

#Data (DAVIDdata must contain ID, term, adj_pval, genes, count and category, datagenes must contain: ID and logFC)
circ <- circle_dat(DAVIDdata,datosgenes)

#Bubble plot
GOBubble(circ, labels = 3)

#GOCHORD

#We select the data from the category of Biological Process
dataGO <- circ[circ$category == "BP",]
genes <- unique(dataGO[,c("genes","logFC")])

#process <- unique(dataGO[,"term"])

#Stroma
#process <- c("cell division", "mitotic nuclear division", "SRP-dependent cotranslational protein targeting to membrane", "DNA replication", "translational initiation", "cell proliferation", "rRNA processing", "muscle contraction", "response to estradiol")

#Epithelium
process <- c("muscle contraction", "cell division", "mitotic nuclear division")

chord <- chord_dat(circ, genes, process)

head(chord) #1 is assigned to the term, 0 is not

chord <- chord_dat(data = circ, genes = genes, process = process)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 4)


##############
#VENN DIAGRAM#
##############

library(ggVennDiagram)
library(ggplot2)
library(dplyr)

#We get the genes from both the stroma and the epithelium
#We get the totals, the upregulated and the downregulated of each one

stroma <- read.table("David_genes_Str.txt", header=TRUE)
stroma_all <- stroma[,1]

stroma_up <- stroma %>%
  filter(logFC >= 0)
stroma_up <- stroma_up[,1]

stroma_down <- stroma %>%
  filter(logFC < 0)
stroma_down <- stroma_down[,1]

epithelium <- read.table("David_genes_Ep.txt", header=TRUE)
epithelium_all <- epithelium[,1]

epithelium_up <- epithelium %>%
  filter(logFC >= 0)
epithelium_up <- epithelium_up[,1]

epithelium_down <- epithelium %>%
  filter(logFC < 0)
epithelium_down <- epithelium_down[,1]

RColorBrewer::display.brewer.all() #To see all color palettes

#ALL
x <- list(A=stroma_all, B=epithelium_all)
ggVennDiagram(x, category.names = c("Stroma", "Epithelium"), color = "black", lwd = 0.8, lty = 1) + 
  scale_color_brewer(palette = "Paired")

#UPREGULATED

x <- list(A=stroma_up, B=epithelium_up)
p1 <- ggVennDiagram(x, category.names = c("Stroma", "Epithelium"), lwd = 0.8, lty = 1) + 
  scale_fill_distiller(palette = "Reds")

#DOWNREGULATED

x <- list(A=stroma_down, B=epithelium_down)
p2 <- ggVennDiagram(x, category.names = c("Stroma", "Epithelium"), lwd = 0.8, lty = 1) + 
  scale_color_brewer(palette = "Blue")

library(ggplotify)
library(cowplot)

p1 <- as.grob(p1)
p2 <- as.grob(p2)

cowplot::plot_grid(p1, p2, nrow = 1, labels = c("Genes Upregulated", "Genes Downregulated"), label_size = 12)

##################
#GENES COMPARISON#
##################

library(sqldf)

df1 <- as.data.frame(stroma_up)
df2 <- as.data.frame(epithelium_up)

#different genes

dif1 <- sqldf('SELECT * from df1 except select * FROM df2')
write.table(dif1, "only_stroma_up.txt")

dif2 <- sqldf('SELECT * from df2 except select * FROM df1')
write.table(dif2, "only_epithelium_up.txt")

#common genes

dup <- sqldf('SELECT * FROM df1 INTERSECT SELECT * FROM df2')
colnames(dup) <- "Epithelium_Stroma"
write.table(dup, "Common_genes_up.txt")

