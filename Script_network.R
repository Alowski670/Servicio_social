####################################################################### Análisis con topGO

######## INSTALANDO PAQUETES
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("topGO", force = T)

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("GO.db", force = T)

#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("Rgraphviz", force = T)
########## Datos
library(pacman)

pacman::p_load(topGO, biomaRt, GO.db, dplyr)

datos <- readRDS("Datos_Costa_Castelo/ttFIRinELGANsUCtissue.rds")
datosS <- datos$Symbol

# genesDEdown conjunto de genes a comparar con el universo datos
genesDEdownS <- genesDEdown$Symbol

### convirtiendo a formato character
datosS <- as.character(datosS)
genesDEdownS <- as.character(genesDEdownS)

# Asociar los genes a terminos GO
db <- useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl', host="https://www.ensembl.org")
go_ids <- getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003', 'hgnc_symbol'), filters='external_gene_name', values= datosS, mart=db)
View(go_ids)
save(go_ids, file = "go_ids.RData")
load("go_ids.RData")
###### construyendo objeto topGO
# Lista de anotacion
gene_2_GO <- unstack(go_ids[,c(1,2)])

# Quitando genes candidatos sin anotacion GO
genesDEdownS
keep <- genesDEdownS %in% go_ids[,2] 
keep
keep <- which(keep == T)
genesDEdownS <- genesDEdownS[keep]
genesDEdownS
geneList <- factor(as.integer(datosS %in% genesDEdownS))
#geneList <- factor(as.integer(grep(datosS, genesDEdownS)))
names(geneList) <- datosS
View(geneList)

### creando el objeto topGo
GOdata_down <- new('topGOdata', ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_GO, nodeSize = 20)

### test algoritmo clasico de fisher
classic_fisher_result <- runTest(GOdata_down, algorithm='classic', statistic='fisher')

### test con algoritmo weight01
weight_fisher_result <- runTest(GOdata_down, algorithm='weight01', statistic='fisher') 

### tabla final con los resultados
allGO <- usedGO(GOdata_down)
all_res <- GenTable(GOdata_down, weightFisher=weight_fisher_result, orderBy='weightFisher', topNodes=length(allGO))

#correcion BH para los valores de p
p.adj <- round(p.adjust(all_res$weightFisher,method="BH"),digits = 4)

# creando un archivo con todos los estadisticos del analisis go
all_res_final <- cbind(all_res,p.adj)
all_res_final <- all_res_final[order(all_res_final$p.adj),]

# lista de terminos go siginificativos antes de las correcciones
results.table.p <- all_res_final[which(all_res_final$weightFisher<=0.001),]

# lista de terminos go siginificativos despues de las correcciones
results.table.bh <- all_res_final[which(all_res_final$p.adj<=0.05),]

#50 antologias por terminos de p ajustados
write.table(all_res_final[1:50,],"summary_topGO_analysis.csv",sep=",",quote=FALSE,row.names=FALSE)

View(all_res_final)
library(Rgraphviz)

pdf(file='topGOPlot_fullnames.pdf', height=12, width=12, paper='special', pointsize=18)
showSigOfNodes(GOdata_down, score(classic_fisher_result), useInfo = "all", sigForAll=FALSE, firstSigNodes=10,.NO.CHAR=50)
dev.off()


GOdata@nodeSize
genes(GOdata)
k <- (getSigGroups(GOdata, test.stat))
k@geneData

data("GOdata")
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")

resultFis <- getSigGroups(GOdata_down, test.stat)
View(resultFis@score)

################################################ Genes sobre
genesDEupS <- genesDEup$Symbol

### convirtiendo a formato character
datosS <- as.character(datosS)
genesDEupS <- as.character(genesDEupS)

# Asociar los genes a terminos GO
#(devtools::install_version("dbplyr", version = "2.3.4"))
db_up <- useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl', host="https://www.ensembl.org")
go_idsup <- getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003', 'hgnc_symbol'), filters='external_gene_name', values= datosS, mart=db_up)
save(go_idsup, file = "go_idsup.RData")
load("go_idsup.RData")
View(go_idsup)
###### construyendo objeto topGO
# Lista de anotacion
gene_2_GO_up <- unstack(go_idsup[,c(1,2)])

# Quitando genes candidatos sin anotacion GO
keep <- genesDEupS %in% go_ids[,2] 
keep
keep <- which(keep == T)
genesDEupS <- genesDEupS[keep]
genesDEupS
geneList_up <- factor(as.integer(datosS %in% genesDEupS))
#geneList <- factor(as.integer(grep(datosS, genesDEdownS)))
names(geneList_up) <- datosS

### creando el objeto topGo
GOdataup <- new('topGOdata', ontology='BP', allGenes = geneList_up, annot = annFUN.gene2GO, gene2GO = gene_2_GO_up)

### test algoritmo clasico de fisher
classic_fisher_result_up <- runTest(GOdata, algorithm='classic', statistic='fisher')

### test con algoritmo weight01
weight_fisher_result_up <- runTest(GOdata, algorithm='weight01', statistic='fisher') 

### tabla final con los resultados
allGOup <- usedGO(GOdataup)
all_res_up <-GenTable(GOdataup, weightFisher=weight_fisher_result_up, orderBy='weightFisher', topNodes=length(allGOup))

#correcion BH para los valores de p
p.adjup <- round(p.adjust(all_res_up$weightFisher,method="BH"),digits = 4)

# creando un archivo con todos los estadisticos del analisis go
all_res_final_up <- cbind(all_res_up,p.adj)
all_res_final_up <- all_res_final_up[order(all_res_final_up$p.adj),]

# lista de terminos go siginificativos antes de las correcciones
results.table.p_up <- all_res_final_up[which(all_res_final_up$weightFisher<=0.001),]

# lista de terminos go siginificativos despues de las correcciones
results.table.bh_up <- all_res_final_up[which(all_res_final_up$p.adj<=0.05),]

#50 antologias por terminos de p ajustados
write.table(all_res_final_up[1:50,],"summary_topGO_analysis.csv",sep=",",quote=FALSE,row.names=FALSE)

pdf(file='topGOPlotup_fullnames.pdf', height=12, width=12, paper='special', pointsize=18)
showSigOfNodes(GOdataup, score(classic_fisher_result_up), useInfo = "all", sigForAll=FALSE, firstSigNodes=5,.NO.CHAR=50)
dev.off()

