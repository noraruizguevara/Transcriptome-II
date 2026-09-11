# Transcriptome-II
A set of R codes for DEG

# PASO 0: ¿Qué motiva el desarrollo de estos comandos?
```r
Las tablas de conteos individuales obtenidas con "Salmon" necesitan concatenarse para realizar analisis de expresion diferencial de genes, para ello la libreria "edgeR" es especialmente diseñada.
```

# PASO 1: Instalar y cargar las librerias necesarias
```r
setwd("C:/Users/HP/Documents/lepidium/metatranscriptoma")
dir()
as.data.frame(dir())

# 0: Instalar paquetes
install.packages(c("tximport", "DESeq2", "pheatmap"))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# BiocManager::install("tximport")
# BiocManager::install("DESeq2")
# BiocManager::install("edgeR")
# BiocManager::install("statmod")
# install.packages("ggplot2", type = "binary")
# install.packages("shadowtext")
# BiocManager::install("clusterProfiler", update = TRUE, ask = FALSE)

# Cargar librerías
library(tximport)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(edgeR)
library(statmod)
```

# PASO 2: cargar los archivos de conteos generados con SALMON
```r

# 1. Importar conteos (solo si no se ha creado la matriz de conteos previamente)
quant_dir <- "quant_meta"
sample_names <- c("SRR2922712", "SRR2922713", "SRR2922714", 
                  "SRR2922715", "SRR2922716", "SRR2922717",
		     "SRR2960160","SRR2960161","SRR7003712",
		     "SRR7003713","SRR7003714")
files <- file.path(quant_dir, sample_names, "quant.sf")
names(files) <- sample_names

# tximport
txi <- tximport(files, type = "salmon", txOut = TRUE)
names(txi)

# abundance
head(txi$abundance)
# raw reads counts
head(txi$counts)
# average transcript length for each gene
head(txi$length)
# if defined
head(txi$countsFromAbundance)

# 2. guardar la tabla de conteos como un nuevo archivo que
# se pueda cargar facilmente 
# write.table(cts,"metatranscriptome_counts.csv",row.names=T,sep=",")

# 3. retener la tabla de conteos como un objeto llamado "cts"

cts <- as.data.frame(txi$counts)
head(cts)
dim(cts)
```

# PASO 3: cargar la metadata en formato CSV
```r
# cts <- read.csv("metatranscriptome_counts.csv")
# head(cts)
# dim(cts)

metadata <- read.csv("metadata.csv")
```

# PASO 4: transformar la tabla de conteos al formato DGE
```r

# 3. DGEList object
y <- DGEList(cts)
# explora los objetos de "y"
names(y)
head(y$counts)
# ahora debemos configurar la columna "group" de "y$samples"
y$samples
# write.table(as.data.frame(y$samples),"samples_1.csv",row.names=T,sep=",")
# y$samples$group <- metadata$individual
# y$samples$group <- metadata$target_study
# write.table(as.data.frame(y$samples),"samples_2.csv",row.names=T,sep=",")
```

# PASO 5: remover los transcritos/genes con bajos recuentos de acuerdo con un criterio de DGE
```r
# 4. remove low counts
keep <- filterByExpr(y)
head(keep)
summary(keep)
y <- y[keep, keep.lib.sizes=FALSE]
y$samples
```

# PASO 6: estimar los factores de normalizacion
```r
# 5. normalization FACTORS
y <- calcNormFactors(y)
y$samples
# write.table(as.data.frame(y$samples),"samples_3.csv",row.names=T,sep=",")
plotMDS(y)
```

# PASO 7: obtener una "design matrix"
```r
# 6. design matrix
design <- model.matrix(~0+y$samples$group)
colnames(design) <- levels(factor(y$samples$group))
```

# PASO 8: estimar parametros de dispersion (3)
```r
# 7. disperssion
y <- estimateDisp(y,design,robust=T)
names(y)
y$common.dispersion
head(y$trended.dispersion)
head(y$tagwise.dispersion)
par(mfrow=c(1,1))
plotBCV(y)
```

# PASO 9: emplear el modelo de distribucion negativa binomial a traves de la transformacion quasi-likelihood de "edgeR"
```r
# 8. GLM model adjustment
fit <- glmQLFit(y, design, robust=T)
names(fit)
head(fit$counts)
head(fit$fitted.values)
fit$samples
plotQLDisp(fit)
fit <- na.omit(fit)
```

# PASO 1O: configurar diferentes pruebas
```r
# 9. tests
# lev_hyp <- makeContrasts(leaves-hypocotyls,levels=design)
# hyp_lev <- makeContrasts(hypocotyls-leaves,levels=design)
# stI_II <- makeContrasts(stage_I-stage_II,levels=design)
# stI_III <- makeContrasts(stage_I-stage_III,levels=design)
# stII_III <- makeContrasts(stage_II-stage_III,levels=design)
# by <- makeContrasts(black-yellow,levels=design)
# bv <- makeContrasts(black-violet,levels=design)
# yv <- makeContrasts(yellow-violet,levels=design)

# res1 <- glmQLFTest(fit,contrast=lev_hyp)
# res1 <- glmQLFTest(fit,contrast=hyp_lev)
# res1 <- glmQLFTest(fit,contrast=stI_II)
# res1 <- glmQLFTest(fit,contrast=by)

# comparison <- "leave_hypocotyl"
# comparison <- "hypocotyl_leave"
# comparison <- "stage I vs. II"
# comparison <- "black vs yellow"
```

# PASO 11: obtener valores corregidos de P > FDR
```r
names(res1)
head(res1$table)
res1_corr <- topTags(res1, n=Inf)
dim(res1_corr)
names(res1_corr)
res1a <- res1_corr$table
head(res1a)
# write.table(res1a,"DEGs.tsv",sep="\t", row.names=F, quote=F)
```

# PASO 12: Filtrar los datos y plotear
```r
res1b <- res1a[res1a$FDR <= 0.001 & abs(res1a$logFC) >= 1, ]
dim(res1b)
dim(res1a)
head(res1b)

res1a$DE <- "NO"
res1a$DE[res1a$logFC >= 1 & res1a$FDR <= 0.001] <- "UP"
res1a$DE[res1a$logFC <= -1 & res1a$FDR <= 0.001] <- "DOWN"

data <- res1a
head(data)
dim(data)

up_regulated <- length(data[data$DE %in% "UP", 1])
down_regulated <- length(data[data$DE %in% "DOWN", 1])

print("########################")
print(paste0("Numero de transcritos up-regulated: ",up_regulated))
print(paste0("Numero de transcritos down-regulated: ",down_regulated))
print("########################")


-log10(0.1)
-log10(0.01)
-log10(0.001)
-log10(0.0001)

par(mfrow=c(1,1))
p1 <- ggplot(data,aes(x=logFC,y=-log10(FDR),col=DE)) +
      geom_point(size=0.5) + theme_minimal() +
      geom_vline(xintercept=c(-1,1), col="red", linetype="dashed") +
      geom_hline(yintercept= -log10(0.001), col="black", linetype="dashed") +
      ggtitle(comparison, subtitle = "METATRANSCRIPTOME") + xlab("log2FC")
p1

head(data)

data$abs <- abs(data$logFC)
data2 <- data[order(data$FDR,-data$abs),]
data2 <- data[order(-data$abs),]
head(data2)
dim(data2)
```
