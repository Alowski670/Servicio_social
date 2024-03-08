library(dplyr)
library(enrichplot)
library(GOSemSim)
library(networkD3)
library(clusterProfiler)
library(kableExtra)
library(DT)
library(tools)
library(plotly)
library(htmlwidgets)
library(tidyverse)
library(networkD3)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(highcharter)
library(networkD3)
library(ggplot2)

uRead <- function(file){
  ext <- file_ext(file)
  print(paste0("Es un archivo tipo: ",ext))
  print(paste0("Links con funciones para leer archivos:"))
  print(paste0("https://bookdown.org/rdpeng/rprogdatascience/getting-data-in-and-out-of-r.html"))
  print(paste0("https://bookdown.org/rdpeng/rprogdatascience/using-the-readr-package.html"))
  print(paste0("https://rstudio-education.github.io/hopr/dataio.html"))
}
uRead("Datos_Costa_Castelo/ttFIRinELGANsUCtissue.rds")

####### Separar datos
GeneDE <- function(lista, pvalor, patron, GeneNumber){
  list_gene <-readRDS(lista)
  GeneNumber <- as.numeric(GeneNumber)
  print("Dimensiones")
  print(dim(list_gene))
  head(list_gene, n= 3)
  
  genesDE <- list_gene[list_gene$adj.P.Val < pvalor & 
                         abs(list_gene$logFC) > log2(1.5), ]
  
  patron2 <- list_gene$Symbol[grep(patron, list_gene$Symbol)]
  
  ## Tabla de contingencia para prueba de Fisher
  m <- length(patron2) 
  N <- nrow(list_gene) 
  n <- dim(genesDE)[1]
  k <- length(intersect(genesDE$Symbol, patron2)) 
  dnames <- list(CG=c("dentro","fuera"),
                 DE=c("si","no"))
  t <- matrix(c(k, n-k, m-k, N+k-n-m),
              nrow=2, ncol=2, dimnames=dnames)
  t
  head(list_gene, n = 4)
  
  ### Crear tablas interactivas con DT
  tabla_DT <- datatable(list_gene)
  htmlwidgets::saveWidget(tabla_DT, "tabla_DT.html")
  ### Tablas interactivas con kableextra
  genesDEup <- filter(genesDE, logFC > 0)
  genesDEup <- genesDEup[order(genesDEup$logFC, decreasing = F),]
  genesDEup2 <- genesDEup[1:GeneNumber, ]
  x <- c(genesDEup2$Symbol)
  y <- c(genesDEup2$logFC)
  text <- c(genesDEup2$Symbol)
  data <- data.frame(x, y, text)
  fig <- plot_ly(data, x = ~x, y = ~y, type = 'bar',
                 text = y, textposition = 'auto',
                 marker = list(color = 'rgb(58,10,232)'))
  
  fig <- fig %>% layout(title = "Sobreexpresados",
                        
                        xaxis = list(title = "Simbolo"),
                        
                        yaxis = list(title = "LogFC"))
  
  htmlwidgets::saveWidget(fig, "up.html", selfcontained = T, libdir = "lib")
  
  genesDEdown <- filter(genesDE, logFC < 0)
  genesDEdown <- genesDEdown[order(genesDEdown$logFC, decreasing = F),]
  genesDEdown2 <- genesDEdown[1:GeneNumber, ]
  x <- c(genesDEdown2$Symbol)
  y <- c(genesDEdown2$logFC)
  text <- c(genesDEdown2$Symbol)
  data <- data.frame(x, y, text)
  fig2 <- plot_ly(data, x = ~x, y = ~y, type = 'bar',
                  text = y, textposition = 'auto',
                  marker = list(color = 'rgb(58,300,23)'))
  
  fig2 <- fig2 %>% layout(title = "Subexpresados",
                          
                          xaxis = list(title = "Simbolo"),
                          
                          yaxis = list(title = "LogFC"))
  
  htmlwidgets::saveWidget(fig2, "down.html", selfcontained = T, libdir = "lib")
  
  #### Fig tree
  genesDEdown100 <- genesDEdown[1:10, ]
  tree_down <- genesDEdown100 %>%  
    hchart(type = "treemap", hcaes(x = Symbol, value = adj.P.Val, color = adj.P.Val)) |>
    hc_title(text ="Genes sub expresados organizados de acuerdo a su valor de p") |>
    hc_colorAxis(stops = color_stops(colors = viridis::inferno(10)))
  
  htmlwidgets::saveWidget(tree_down, "TreemapDE_down.html", selfcontained = T, libdir = "lib")
  
  genesDEup100 <- genesDEup[1:10, ]
  tree_up <- genesDEup100 %>%  
    hchart(type = "treemap", hcaes(x = Symbol, value = adj.P.Val, color = adj.P.Val)) |>
    hc_title(text ="Genes sobre expresados organizados de acuerdo a su valor de p") |>
    hc_colorAxis(stops = color_stops(colors = viridis::inferno(10)))
  
  htmlwidgets::saveWidget(tree_up, "TreemapDE_up.html", selfcontained = T, libdir = "lib")
  
  save(genesDE, file = "genesDE.RData")
  save(genesDEup, file = "genesDEup.RData")
  save(genesDEdown, file = "genesDEdown.RData")
}
GeneDE("Datos_Costa_Castelo/ttFIRinELGANsUCtissue.rds", pvalor = 0.01, patron = "TLR", GeneNumber = 10)

GO_analysis <- function(){
  load("genesDEup.RData")
  load("genesDEdown.RData")
  #### Analisis GO
  genesDEupSymbol <- genesDEup[ , 1]
  genesDEupSymbol
  GO_up <- enrichGO(gene = genesDEupSymbol, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
  GO_up <- as.data.frame(GO_up)
  
  genesDEdownSymbol <- genesDEdown[ , 1]
  genesDEdownSymbol
  GO_down <- enrichGO(gene = genesDEdownSymbol, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
  GO_down <- as.data.frame(GO_down)
  
  ### Gráficas con las 20 categorías más comunes, de acuerdo al análisis GO
  png("graphGO_up.png")
  GO_up10 <- GO_up[order(GO_up$Count, decreasing = T), ]
  GO_up10 <- GO_up10[1:10,]
  graphGO_up <- ggplot(GO_up10, aes(x=Description, y=Count, fill = Description)) + 
    geom_bar(stat = "identity") +
    theme(legend.position="none")+ 
    coord_flip()
  graphGO_up
  dev.off()
  
  png("graphGO_down.png")
  GO_down10 <- GO_down[order(GO_down$Count, decreasing = T), ]
  GO_down10 <- GO_down10[1:10,]
  graphGO_up <- ggplot(GO_down10, aes(x=Description, y=Count, fill = Description)) + 
    geom_bar(stat = "identity") +
    theme(legend.position="none")+ 
    coord_flip()
  graphGO_up
  dev.off()
  
  save(GO_up, file = "GO_up.RData")
  save(GO_down, file = "GO_down.RData")
  save(net_down, file = "net_down.RData")
  save(net_up, file = "net_up.RData")
}

d_up <- godata("org.Hs.eg.db", ont = "BP")
save(d_up, file = "d_up.RData")

d_down <- godata("org.Hs.eg.db", ont = "BP")
save(d_down, file = "d_down.RData")

netGO <- function(lim){
###################### A PARTIR DE AQUI EN TEORIA FUNCIONA
### para los down
lim <- as.numeric(lim)
load("d_up.RData")
load("d_down.RData")
load("GO_up.RData")
load("GO_down.RData")
load("net_down.RData")
load("net_up.RData")
###### cortar la red, para asi poder ver como interactuan entre si
GO_downT <- GO_down
GO_downT <- as.data.frame(GO_downT)
GO_downT$Description <- chartr(old = " ", new = "_", GO_downT$Description)
#GO_downT <- filter(GO_downT, Count > 20)
#GO_downT <- enrichGO(GO_downT)

invalid_down <- filter(GO_downT, Count < lim)
invalid_down <- invalid_down$Description

down <- pairwise_termsim(GO_downT, semData = d_down)
similarity_down <- down@termsim
networkDown <- data.frame()

cat_namesD <- c(colnames(similarity_down))
length(cat_namesD)

num <- 1
Source <- c()
Target <- c()
Value <- c()
for(i in 1:dim(similarity_down)[1]){
  Source <- c(Source, rep(cat_namesD[num], times = dim(similarity_down)[1], each = 1:dim(similarity_down)[1]))
  Target <- c(Target, rep(cat_namesD, times = 1, each = 1))
  Value <- c(Value, similarity_down[num, ])
  num <- num + 1
  
}

networkDown <- data.frame(Source, Target, Value)
networkDown_final <- networkDown

#### quitando los invalidos
num <- 1
for(i in 1:length(invalid_down)){
  networkDown_final <- filter(networkDown_final, Source != invalid_down[num])
  num <- num + 1
}
View(networkDown_final)
#### quitando NAs
#networkDown_final$Value[networkDown_final$Value <= 0.05] <- "NA"
networkDown_final <- filter(networkDown_final, Value != "NA")
networkDown_final

netD <- networkDown_final[order(networkDown_final$Value, decreasing = T), ]
net <- networkDown_final
netD100 <- netD[1:100, ] ### recuerda que debes correr de nuevo netD, para que se genere una vez mas la tabla
net <- net[1:300,]
netD100

Net_down <- simpleNetwork(netD100, width = 400, height = 250, fontSize = 14, nodeColour ="blue", zoom = T)
saveNetwork(Net_down, file = "Net_down.html", selfcontained = T)

#### con enrich plot
emapplot(down, showCategory = 10)
cnetplot(down)


########################## PARA LOS UP
GO_upT <- as.data.frame(GO_up)
GO_upT$Description <- chartr(old = " ", new = "_", GO_upT$Description)
#GO_upT <- GO_upT@result[order(GO_upT@result$Count, decreasing = T),]
#GO_upT <- filter(GO_upT, Count > 20)
#GO_upT <- enrichGO(GO_upT)

invalid_up <- filter(GO_upT, Count < lim)
invalid_up <- invalid_up$Description

up <- pairwise_termsim(GO_up, semData = d_up)
similarity_up <- up@termsim

network_up <- data.frame()

cat_namesD <- c(colnames(similarity_up))
length(cat_namesD)

num <- 1
Source <- c()
Target <- c()
Value <- c()
for(i in 1:dim(similarity_up)[1]){
  Source <- c(Source, rep(cat_namesD[num], times = dim(similarity_up)[1], each = 1:dim(similarity_up)[1]))
  Target <- c(Target, rep(cat_namesD, times = 1, each = 1))
  Value <- c(Value, similarity_up[num, ])
  num <- num + 1
  
}

Source
Target
Value
network_up <- data.frame(Source, Target, Value)
networkUP_final <- network_up

##### quitando los invalidos
num <- 1
for(i in 1:length(invalid_up)){
  networkDown_final <- filter(networkUP_final, Source != invalid_up[num])
  num <- num + 1
}

#### quitando NAs
networkUP_final <- filter(networkUP_final, Value != "NA")

#### red
netD_UP <- networkUP_final[order(networkUP_final$Value, decreasing = T), ]
net_UP <- networkUP_final
netD_UP <- netD_UP[1:300, ]
net_UP <- net_UP[1:300,]
#### funciona ahora añadir los conteos a las categorias

NET_UP <- simpleNetwork(netD_UP, width = 400, height = 250, fontSize = 14, nodeColour = "red", zoom = T)
saveNetwork(NET_UP, file = "Net_up.html", selfcontained = T)
emapplot(up, showCategory = 10)
cnetplot(up)
}

netGO( 20)

############################33
#############################
#   FORCE NETWORK
library(networkD3)
library(dplyr)
### Acortamos la red a 400 elementos, para que la compu no explote
netD_UP2 <- netD_UP[1:400, ]

### Asignamos los datos de la columna source de la matriz
## de adyacencia a un objeto src, lo mismo para target
src <- netD_UP2$Source 

target <- netD_UP2$Target
### Despues creamos un objeto que contenga los vectores src y target, nuestra nueva
# matriz de adyacencia
networkData <- data.frame(src, target, stringsAsFactors = FALSE)


##indicamos que no se repita ninguna categoría dentro de los vectores con unique,
# para crear un data frame que tenga solo el total de categorías
nodes <- data.frame(name = unique(c(networkData$src, networkData$target)))

### Para sacar los conteos asociados a cada categoría utilizamos %in$. Para que
# nos indique cuáles de los elementos de las categorías que tenemos se encuentran
# en el objeto del análisis GO (tiene los conteos para cada categoría)
counts <- GO_upT[GO_upT$Description %in%nodes$name , ]
### Seleccionamos solo las columnas con los nombres de las categorías y conteos (genes asociados a la categoría)
counts <- subset(counts, select = c(Description, Count))
counts <- counts[order(counts$Description),]
#View(counts)

### Asignamos los conteos a nuestro nodes2, contiene las categorías (nodos), y
# sus respectivos conteos 
nodes2 <- nodes[order(nodes$name), ]
nodes2 <- as.data.frame(nodes2)
nodes2$counts <- counts$Count
nodes2 <- as.data.frame(nodes2)
#View(nodes2)

### Creamos los grupos, donde tener 20 o más genes asociados hace que la categoría
# se considere enriquecida
nodes2$group <- ifelse(nodes2$counts >= 50, "Enriquecido", "No_enriquecido")

### Para crear las conexiones entre los nodos utilizamos "match". Es decir, que se tomarán
# aquellos elementos de la matriz de adyacencia, en su columna src, que correspondan
# con las categorías que están en nodes.
links <- data.frame(source = match(networkData$src, nodes$name) - 1,
                    target = match(networkData$target, nodes$name) - 1)

### Generar colores para los grupos
colorJS <- 'd3.scaleOrdinal().range(["forestgreen", "darkcyan"])'
### Crear RED

Net_UP_Force <- forceNetwork(Links = links, Nodes = nodes2, Source = "source",
                  Target = "target", NodeID ="nodes2", Group = "group",
                  opacity = 1, opacityNoHover = 1, colourScale = colorJS,
                  zoom = T)
saveNetwork(Net_UP_Force, file = "Net_up_Grupos.html", selfcontained = T)
Net_UP_Force
############### para los DOWN
netD_Down <- netD[1:400, ]
netD_Down$Source <- chartr(old = " ", new = "_", netD_Down$Source)
netD_Down$Target <- chartr(old = " ", new = "_", netD_Down$Target)

src <- netD_Down$Source
target <- netD_Down$Target
### Despues creamos un objeto que contenga los vectores src y target, nuestra nueva
# matriz de adyacencia
networkData <- data.frame(src, target, stringsAsFactors = FALSE)


##indicamos que no se repita ninguna categoría dentro de los vectores con unique,
# para crear un data frame que tenga solo el total de categorías
nodes <- data.frame(name = unique(c(networkData$src, networkData$target)))
View(nodes)
### Para sacar los conteos asociados a cada categoría utilizamos %in$. Para que
# nos indique cuáles de los elementos de las categorías que tenemos se encuentran
# en el objeto del análisis GO (tiene los conteos para cada categoría)
counts <- GO_downT[GO_downT$Description %in% nodes$name , ]
### Seleccionamos solo las columnas con los nombres de las categorías y conteos (genes asociados a la categoría)
counts <- subset(counts, select = c(Description, Count))
counts <- counts[order(counts$Description),]
View(counts)

### Asignamos los conteos a nuestro nodes2, contiene las categorías (nodos), y
# sus respectivos conteos 
nodes2 <- nodes[order(nodes$name), ]
nodes2 <- as.data.frame(nodes2)
nodes2$counts <- counts$Count
nodes2 <- as.data.frame(nodes2)
#View(nodes2)

### Creamos los grupos, donde tener 20 o más genes asociados hace que la categoría
# se considere enriquecida
nodes2$group <- ifelse(nodes2$counts >= 10, "Enriquecido", "No_enriquecido")

### Para crear las conexiones entre los nodos utilizamos "match". Es decir, que se tomarán
# aquellos elementos de la matriz de adyacencia, en su columna src, que correspondan
# con las categorías que están en nodes.
links <- data.frame(source = match(networkData$src, nodes$name) - 1,
                    target = match(networkData$target, nodes$name) - 1)

### Generar colores para los grupos
colorJS <- 'd3.scaleOrdinal().range(["gold", "gray0"])'
### Crear RED

Net_Down_Force <- forceNetwork(Links = links, Nodes = nodes2, Source = "source",
                      Target = "target", NodeID ="nodes2", Group = "group",
                      opacity = 1, opacityNoHover = 1, colourScale = colorJS,
                      zoom = T)
saveNetwork(Net_Down_Force, file = "Net_down_Grupos.html", selfcontained = T)
Net_Down_Force
