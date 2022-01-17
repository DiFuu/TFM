setwd("~/../Documents/R")

setwd("~/../Dropbox/UPO/TFM/R/Final")

##################
#Preprocesamiento#
##################

#eliminamos las variables guardadas
rm(list=ls())

#Leemos los datos y lo pasamos a un .txt

library(readxl)

file <- read_excel("1-Datos.xlsx", 1)
write.table(file, "1-Datos.txt",sep="\t")

#Diseño del experimento
file2 <- read_excel("2-Muestras.xlsx", 1)
write.table(file2, "2-Muestras.txt",sep="\t")

#Datos
seqdata0 <- read.delim("1-Datos.txt", stringsAsFactors = FALSE, row.names = 1)

head(seqdata0)
dim(seqdata0)

#Eliminamos las columnas desde la 1 hasta la la 110
seqdata <- seqdata0[,-(1:110)]

#Diseño del experimento
end_layout <- read.delim("2-Muestras.txt")

#Cambiamos el nombre de las muestras, partiendo desde end_layout

colnames(end_layout) #En este caso, nos interesa well.Sample.Cells_type.Day2 y Library
type <- substr(end_layout$Cells_type,1,1)

new_name <- paste(end_layout$well,end_layout$Sample,type,end_layout$Day2,end_layout$Library,sep=".")

end_layout$name <- new_name
colnames(seqdata) <- new_name

#El nombre de los genes en seqdata coincide con los datos de la columna name en endometrium layout
table(colnames(seqdata)==end_layout$name)

#Vamos a filtrar los genes con expresión muy baja, aquellos que tienen menos de 0.5 recuento por millón de lecturas (CPM)
#en la mitad de las muestras. Para ello, utilizamos el paquete "edgeR" que proporciona la función cpm. 

library(limma)
library(edgeR)

#recuentos de expresiones
mycounts <- seqdata
dim(mycounts)

mycounts[1:5,1:3]

#Filtrado
#Mantenemos los genes con al menos 0.5 recuento por millón de lecturas (cpm) en al menos 2 muestras 
isexpr <- rowSums(cpm(mycounts)>0.5) >= 2
table(isexpr)

mycounts <- mycounts[isexpr,]
genes <- rownames(mycounts)

dim(mycounts)

#En nuestro caso, queremos hacer dos análisis distintos, uno para las células del estroma y otro para las células del epitelio
#Ambos análisis van a ser completamente independientes porque lo que queremos comparar es el efecto del tiempo,no del tejido en sí

#En vez de renombrar todos los datos de entrada y hacerlo doble, se va a ejecutar la línea del epitelio o del estroma en función
#de cuál convenga

#Aquí sólo hay que cambiar entre Stroma o Epithelium en función de los datos que queramos analizar

library(dplyr)

end_layout <- end_layout %>%
  filter(Cells_type == "Epithelium")
mycounts <- mycounts[,end_layout$name]
dim(mycounts)

#Para la normalización usamos la función voom del paquete limma que normaliza los recuentos de lectura
#y aplica un modelo lineal antes del análisis estadístico de expresión diferencial

#grupo de condiciones
group<-factor(end_layout$Day2)

#matriz de diseño para limma
design <- model.matrix(~0+group)
#sustituimos "grupo" de los nombres de las columnas de diseño 
colnames(design)<- gsub("group","",colnames(design))
colnames(design) <- c("LH2","LH8")
design
#write.table(design[,1],file="Str_class.txt",row.names=F,quote=F,sep="\t")
write.table(design[,1],file="Ep_class.txt",row.names=F,quote=F,sep="\t")

#factores de normalización entre bibliotecas
nf <- calcNormFactors(mycounts)

#normalizamos los recuentos de lectura con la función 'voom' 
y <- voom(mycounts,design,lib.size=colSums(mycounts)*nf)
#extraemos los recuentos de lectura normalizados 
counts.voom <- y$E

#guardamos los datos de expresión normalizados en el directorio de salida 
#write.table(counts.voom,file="Str_counts.voom.txt",row.names=F,quote=F,sep="\t")
#write.table(mycounts,file="Str_mycounts.txt",row.names=F,quote=F,sep="\t")

write.table(counts.voom,file="Ep_counts.voom.txt",row.names=F,quote=F,sep="\t")
write.table(mycounts,file="Ep_mycounts.txt",row.names=F,quote=F,sep="\t")

#write.table(counts.voom,file="Ep_counts.voom.txt",row.names=F,quote=F,sep="\t")
#write.table(mycounts,file="Ep_mycounts.txt",row.names=F,quote=F,sep="\t")

#ajustamos el modelo lineal para cada gen
fit <- lmFit(y,design)

#construimos la matriz de contraste 
cont.matrix <- makeContrasts(LH8-LH2,levels=design)
cont.matrix 

#calculamos coeficientes estimados y errores estándar para un conjunto dado de contrastes 
fit <- contrasts.fit(fit, cont.matrix)

#estadísticos t moderados de expresión diferencial mediante la moderación empírica de Bayes de los errores estándar 
fit <- eBayes(fit)
options(digits=3)

#output fit
dim(fit)

#La función topTable resume la salida de limma en un formato de tabla.
#Se seleccionan los genes con un pvalue < 0.05, no hacemos filtrado por logFC.

#establecemos el umbral ajustado de pvalue y el umbral de cambio logarítmico 
mypval=0.05
#myfc=2

#nombre del coeficiente para la comparación de interés
colnames(fit$coefficients)

mycoef="LH8 - LH2"
#tabla de salida de los 10 genes más importantes para esta comparación 
topTable(fit,coef=mycoef)

#tabla completa ("n = número de genes en el ajuste") 
limma.res <- topTable(fit,coef=mycoef,n=dim(fit)[1])

#solo genes significativos (valor p ajustado <mypval) 
limma.res.pval <- topTable(fit,coef=mycoef,n=dim(fit)[1],p.val=mypval)
dim(limma.res.pval)

#genes significativos con un valor de p ajustado bajo y un FC alto
#limma.res.pval.FC <- limma.res.pval[which(abs(limma.res.pval$logFC)>myfc),]

#para ordenar los genes en función de su número
limma.res.pval$Gene <- as.numeric(rownames(limma.res.pval))
limma.res.pval <- limma.res.pval[order(limma.res.pval$Gene),]

#escribimos la tabla de salida de limma para genes significativos en un archivo delimitado por tabulaciones 
#filename = paste("Genes_STROMA_",mycoef,"_pval",mypval,".txt",sep="")
#write.table(limma.res.pval,file=filename,row.names=F,quote=F,sep="\t")

filename = paste("Genes_EPITHELIUM_",mycoef,"_pval",mypval,".txt",sep="")
write.table(limma.res.pval,file=filename,row.names=F,quote=F,sep="\t")

#Para sacar los datos normalizados al archivo data_normalized.txt

#Necesitamos el nombre de los genes
dim(seqdata0)
seqdata1 <- seqdata0[,-(2:195)]

#obtenemos la posición de los genes que vamos a exportar fruto de la normalización
list_genes2 <- as.integer(rownames(counts.voom))
length(list_genes2)

#Obtenemos los genes correspondientes
genesnorm <- NULL
for (ngene in list_genes2){
  genesnorm[[length(genesnorm) + 1]] <- (seqdata1[(ngene)])
}

genesnorm <- as.character(genesnorm)

#Añadimos los genes
normdata <- counts.voom #copiamos los datos
normdata <- as.data.frame(normdata) #transformamos a dataframe
normdata$Gene <- genesnorm #añadimos el nombre de los genes
normdata <- normdata %>%
  dplyr::select("Gene", everything()) #ordenamos para que la columna gen esté la primera

#Guardamos los archivos
#write.table(normdata,file="Str_data_normalized.txt",row.names=F,quote=F,sep="\t")
#write.table(normdata$Gene,file="Str_id_normalized.txt",row.names=F,quote=F,sep="\t")

write.table(normdata,file="Ep_data_normalized.txt",row.names=F,quote=F,sep="\t")
write.table(normdata$Gene,file="Ep_id_normalized.txt",row.names=F,quote=F,sep="\t")


#Sacamos un archivo con los datos finales para poder continuar con el análisis

#obtenemos la posición de los genes que vamos a exportar
list_genes2 <- as.integer(rownames(limma.res.pval))
length(list_genes2)

#Obtenemos los genes correspondientes
genesfc <- NULL
for (ngene in list_genes2){
  genesfc[[length(genesfc) + 1]] <- (seqdata1[(ngene)])
}

genesfc <- as.character(genesfc)

#Añadimos los genes
fcdata <- limma.res.pval #copiamos los datos
fcdata <- as.data.frame(fcdata) #transformamos a dataframe
fcdata$Gene <- genesfc #añadimos el nombre de los genes
fcdata <- fcdata %>%
  dplyr::select("Gene", everything()) #ordenamos para que la columna gen esté la primera

#Guardamos los archivos
#write.table(fcdata,file="data_Str.txt",row.names=F,quote=F,sep="\t")
#write.table(fcdata$Gene,file="id_Str.txt",row.names=F,quote=F,sep="\t")

write.table(fcdata,file="data_Ep.txt",row.names=F,quote=F,sep="\t")
write.table(fcdata$Gene,file="id_Ep.txt",row.names=F,quote=F,sep="\t")

#Heatmap para la visualización general de los datos

library(DESeq2)

#leemos el archivo
genedata <- read.table("1-Datos.txt",header=T, row.names = "Gene")
attach(genedata)

#seleccionamos los datos de interés
genedata <- genedata[,111:195]
geneNames<-row.names(genedata)

#transformamos el data frame a objeto matriz
is.matrix(genedata)
data.matrix <- as.matrix(genedata)
is.matrix(data.matrix)
data.matrix[1:5,]
summary(data.matrix)

countdata <- data.matrix
head(countdata, 3)

#abrimos el archivo con los datos de las muestras y seleccionamos las columnas de interés
coldata<- read.table("2-Muestras.txt",header=T, row.names = "name")
attach(coldata)

#Cambiamos el nombre de las columnas para que sean exactamente igual a los de las muestras

colnames(coldata) #nos interesa well.Sample.Cells_type.Day2 y Library
type <- substr(coldata$Cells_type,1,1)

name <- paste(coldata$well,coldata$Sample,type,coldata$Day2,coldata$Library,sep=".")

rownames(coldata) <- name
colnames(countdata) <- name

#Ahora cogemos sólo las columnas que nos interesan para el heatmap
coldata <- coldata[,c("Cells_type","Day2")]
coldata[1:5,]
colnames(coldata) <- c("Tissue","Day")

#Separamos los datos en función de si se trata de células del estroma o del epitelio
library(dplyr)
coldatast <- coldata %>%
  filter(Tissue == "Stroma")
countdatast <- countdata[,rownames(coldatast)]

coldataep <- coldata %>%
  filter(Tissue == "Epithelium")
countdataep <- countdata[,rownames(coldataep)]

#Conjunto str + ep
ddsMat <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ Day + Tissue + Day:Tissue)

#Sólo str
ddsMat <- DESeqDataSetFromMatrix(countData = countdatast, colData = coldatast, design = ~ Day)

#Sólo ep
ddsMat <- DESeqDataSetFromMatrix(countData = countdataep, colData = coldataep, design = ~ Day)

#Sale un mensaje de advertencia de que mejor no poner símbolos
nrow(ddsMat)

#Estabilización de la varianza para la visualización de los datos
library("vsn")
vsd <- vst(ddsMat, blind = FALSE)
head(assay(vsd), 3)

colData(vsd)

#Visualización
library("genefilter")
library("pheatmap")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 40)
topVarGenes
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Day","Tissue")])
pheatmap(mat, annotation_col = anno)


######################
#Análisis diferencial#
######################

#eliminamos las variables guardadas
rm(list=ls())

require(multtest)
library(RankProd)

#Leemos el archivo con los datos a analizar
#file<- read.table("Str_data_normalized.txt",header=T, row.names=1)
file<- read.table("Ep_data_normalized.txt",header=T, row.names=1)

attach(file)

geneNames<-row.names(file)

#coaccionar el objeto a una matriz
is.matrix(file)
data.matrix <- as.matrix(file)
is.matrix(data.matrix)
data.matrix[1:5,]

#Clases
#LH2 1, LH8 0
#data.cl <- read.delim("Str_class.txt", stringsAsFactors = FALSE)
data.cl <- read.delim("Ep_class.txt", stringsAsFactors = FALSE)

#lo convertimos a numérico
data.cl <- as.numeric(data.cl$x)

#Test estadístico
#Función mt.teststat(X,classlabel,test="t",nonpara="n")
#test="t", los test se basan en estadísticos t de Welch de dos muestras (varianzas desiguales).
#nonpara="n", se utilizan los datos originales.

teststat<-mt.teststat(file, data.cl, test="t",nonpara="n")

df=42 #41-2 para STR y 44-2 para EP
rawp0<-2*(1-pt(abs(teststat), df))
res.multtest <- mt.rawp2adjp(rawp0, "BH")
res.multtest <-res.multtest$adjp[order(res.multtest$index),]
res.multtest <-cbind(geneNames,res.multtest)

#write.table(res.multtest, "Str_multtest.txt")
write.table(res.multtest, "Ep_multtest.txt")

#seleccionar significantes
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


#Juntamos la tabla de resultados anteriores con estos nuevos datos

#Leer el archivo
#filefc<- read.table("data_Str.txt",header=T, row.names=1)
filefc<- read.table("data_Ep.txt",header=T, row.names=1)

attach(filefc)

#Sacar los nombres de los genes y guardar en una nueva columna
geneNamesFC <- row.names(filefc)
filefc$geneNames <- geneNamesFC

#Unimos los dos dataframes en función del nombre del gen en común
data_FC <- merge(x = filefc, y = BH.significant, by = c("geneNames"))

library(data.table)

setnames(data_FC, "geneNames", "Gene")

#Guardamos
#write.table(data_FC, "Data_Str_BH.txt",row.names=F)
write.table(data_FC, "Data_Ep_BH.txt",row.names=F)


####################
#ANOTACIÓN GENÓMICA#
####################

#eliminamos las variables guardadas
rm(list=ls())

#La anotación de las ID de EntrezGene a partir de los datos de RNAseq se puede realizar utilizando la base de datos de BioMart
#que contiene muchas especies, como humanos, ratones, peces cebra, pollos y ratas

library(biomaRt)
library(org.Hs.eg.db)

#Obtenemos la lista de genes
#data_FC <- read.table("Data_Str_BH.txt",header=T)
data_FC <- read.table("Data_Ep_BH.txt",header=T)

attach(data_FC)

dim(data_FC)
list_genes <- data_FC[,1]

mart<- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org"))

#usamos la base de datos de BioMart para obtener el símbolo del gen y la descripción de estos genes 
detags.IDs <- getBM(
  filters= "hgnc_symbol",
  attributes= c("hgnc_symbol","description"),
  values= as.character(list_genes),
  mart= mart
)

dim(detags.IDs)

head(detags.IDs)

#Seleccionamos sólo los nombres de interés
rownames(detags.IDs)<-detags.IDs$hgnc_symbol
entrez_genes.annot <- detags.IDs[as.character(list_genes),]
description <- entrez_genes.annot$description
data_FC.annot <- cbind(data_FC,description)

#La guardamos la tabla anotada en un archivo
#write.table(data_FC.annot,file="Genes_STROMA_analyzed.txt",sep="\t",row.names=F)
write.table(data_FC.annot,file="Genes_EPITHELIUM_analyzed.txt",sep="\t",row.names=F)

###############
#GO ENRICHMENT#
###############

#eliminamos las variables guardadas
rm(list=ls())

#Leemos el archivo y sacamos el nombre de los genes
#all <- read.table("Genes_STROMA_analyzed.txt", row.names=NULL, quote="\"", comment.char="", header=TRUE)
all <- read.table("Genes_EPITHELIUM_analyzed.txt", row.names=NULL, quote="\"", comment.char="", header=TRUE)

query <- all[,1]

#Exportamos los datos de los genes para poder usar el GOplot
datosgenes <- all[,c("Gene","logFC")]
colnames(datosgenes) <- c("ID","logFC")

#write.table(datosgenes,"David_genes_STR.txt",row.names=F,sep="\t")
write.table(datosgenes,"David_genes_Ep.txt",row.names=F,sep="\t")

#Convertimos de symbol a entrez
library(org.Hs.eg.db)

#Para ver los identificadores disponibles
columns(org.Hs.eg.db)

#Usamos mapIds para obtener los Entrez IDs
query2 <- mapIds(org.Hs.eg.db, query, 'ENTREZID', 'SYMBOL')

#write.table(query2, "query_entrezid_Str.txt")
write.table(query2, "query_entrezid_Ep.txt")

library(rJava)
library(RDAVIDWebService)
library(clusterProfiler)

#Ejecutar DAVID BP
GOBP <- enrichDAVID(query2, minGSSize = 10, idType = "ENTREZ_GENE_ID", annotation = "GOTERM_BP_DIRECT", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, david.user = "diana.fuertes@cragenomica.es")
KEGG <- enrichDAVID(query2, minGSSize = 10, idType = "ENTREZ_GENE_ID", annotation = "KEGG_PATHWAY", david.user = "diana.fuertes@cragenomica.es")

#Representación
GOBP <- setReadable(GOBP, 'org.Hs.eg.db', 'ENTREZID')
KEGG <- setReadable(KEGG, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(GOBP, foldChange=query2)
barplot(KEGG)

#Datos del análisis
dataBP <- data.frame(GOBP$ID, GOBP$Description, GOBP$GeneRatio, GOBP$BgRatio, GOBP$pvalue, GOBP$p.adjust, GOBP$qvalue, GOBP$geneID, GOBP$Count, GOBP@ontology)
dataBP$GOBP.geneID <- gsub("/",", ", dataBP$GOBP.geneID)
dataBP$GOBP.ontology <- gsub("GOTERM_BP_DIRECT","BP", dataBP$GOBP.ontology)
colnames(dataBP) <- c("ID","term","GeneRatio","BgRatio","pvalue","adj_pval","qvalue","genes","count","category")

#Ejecutar DAVID CC
GOCC <- enrichDAVID(query2, minGSSize = 10, idType = "ENTREZ_GENE_ID", annotation = "GOTERM_CC_DIRECT", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, david.user = "diana.fuertes@cragenomica.es")

#Representación
GOCC <- setReadable(GOCC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(GOCC, foldChange=query2)

#Datos del análisis
dataCC <- data.frame(GOCC$ID, GOCC$Description, GOCC$GeneRatio, GOCC$BgRatio, GOCC$pvalue, GOCC$p.adjust, GOCC$qvalue, GOCC$geneID, GOCC$Count, GOCC@ontology)
dataCC$GOCC.geneID <- gsub("/",", ", dataCC$GOCC.geneID)
dataCC$GOCC.ontology <- gsub("GOTERM_CC_DIRECT","CC", dataCC$GOCC.ontology)

#Ejecutar DAVID MF
#PARA EPITELIO NO HAY, ASÍ QUE NO EJECUTAR
GOMF <- enrichDAVID(query2, minGSSize = 10, idType = "ENTREZ_GENE_ID", annotation = "GOTERM_MF_DIRECT", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, david.user = "diana.fuertes@cragenomica.es")

#Representación
GOMF <- setReadable(GOMF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(GOMF, foldChange=query2)

#Datos del análisis
dataMF <- data.frame(GOMF$ID, GOMF$Description, GOMF$GeneRatio, GOMF$BgRatio, GOMF$pvalue, GOMF$p.adjust, GOMF$qvalue, GOMF$geneID, GOMF$Count, GOMF@ontology)
dataMF$GOMF.geneID <- gsub("/",", ", dataMF$GOMF.geneID)
dataMF$GOMF.ontology <- gsub("GOTERM_MF_DIRECT","MF", dataMF$GOMF.ontology)

#Exportar
#write.table(dataBP,"David_results_Str.txt",row.names=F,sep="\t")
#write.table(dataCC,"David_results_Str.txt",row.names=F,col.names=FALSE,sep="\t",append=TRUE)
#write.table(dataMF,"David_results_Str.txt",row.names=F,col.names=FALSE,sep="\t",append=TRUE)

write.table(dataBP,"David_results_Ep.txt",row.names=F,sep="\t")
write.table(dataCC,"David_results_Ep.txt",row.names=F,col.names=FALSE,sep="\t",append=TRUE)
#write.table(dataMF,"David_results_Ep.txt",row.names=F,col.names=FALSE,sep="\t",append=TRUE)

#write.table(KEGG, "KEGG_Str.txt",row.names=F,sep="\t")
write.table(KEGG, "KEGG_Ep.txt",row.names=F,sep="\t")

#Visualización dotplot
#No usado
library("gridExtra")
p1 <- dotplot(GOBP, title = "Biological Process")
p2 <- dotplot(GOCC, title = "Cellular Component")
p3 <- dotplot(GOMF, title = "Molecular Function")

grid.arrange(p1, p2, p3, nrow = 1)

######################
#VISUALIZACIÓN GOPLOT#
######################

#eliminamos las variables guardadas
rm(list=ls())

library(GOplot)

#Datos
#DAVIDdata <- read.table("David_results_Str.txt",header=TRUE)
#datosgenes <- read.table("David_genes_Str.txt", header=TRUE)

DAVIDdata <- read.table("David_results_Ep.txt",header=TRUE)
datosgenes <- read.table("David_genes_Ep.txt", header=TRUE)

#Datos (DAVIDdata debe contener ID, term, adj_pval, genes, count y category, datosgenes debe conter: ID y logFC)
circ <- circle_dat(DAVIDdata,datosgenes)

#Bubble plot
GOBubble(circ, labels = 3)

#Visualización GOCHORD

#Seleccionamos los datos de la categoría de Proceso Biológico
dataGO <- circ[circ$category == "BP",]
genes <- unique(dataGO[,c("genes","logFC")])

#process <- unique(dataGO[,"term"])

#Estroma
#process <- c("cell division", "mitotic nuclear division", "SRP-dependent cotranslational protein targeting to membrane", "DNA replication", "translational initiation", "cell proliferation", "rRNA processing", "muscle contraction", "response to estradiol")

#Epitelio
process <- c("muscle contraction", "cell division", "mitotic nuclear division")

chord <- chord_dat(circ, genes, process)

head(chord) #1 está asignado al término, 0 es que no

chord <- chord_dat(data = circ, genes = genes, process = process)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 4)


##############
#VENN DIAGRAM#
##############

library(ggVennDiagram)
library(ggplot2)
library(dplyr)

#Sacamos los genes tanto del estroma como del epitelio
#Sacamos los totales, los upregulated y los downregulated de cada uno

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

RColorBrewer::display.brewer.all() #Para ver todos las paletas de colores

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

###################
#COMPARACIÓN GENES#
###################

library(sqldf)

df1 <- as.data.frame(stroma_up)
df2 <- as.data.frame(epithelium_up)

#Genes diferentes

dif1 <- sqldf('SELECT * from df1 except select * FROM df2')
write.table(dif1, "only_stroma_up.txt")

dif2 <- sqldf('SELECT * from df2 except select * FROM df1')
write.table(dif2, "only_epithelium_up.txt")

#Genes comunes

dup <- sqldf('SELECT * FROM df1 INTERSECT SELECT * FROM df2')
colnames(dup) <- "Epithelium_Stroma"
write.table(dup, "Common_genes_up.txt")

