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
  
  ### Encontrar los genes diferencialmente expresados de acuerdo al valor de p
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
  ######### Tablas interactivas con kableextra
  
  ### Encontrar genes sobreexpresados
  genesDEup <- filter(genesDE, logFC > 0)
  ### Ordenarlos
  genesDEup <- genesDEup[order(genesDEup$logFC, decreasing = F),]
  ### Gene number representa el número de genes que quieres que se vean en las
  # gráficas de barras
  genesDEup2 <- genesDEup[1:GeneNumber, ]
  
  ###########3 Gráfica con plot_ly
  x <- c(genesDEup2$Symbol) # el eje de las x serán los símbolos de los genes
  y <- c(genesDEup2$logFC) # El logFC nos dirá qué tanto están sobreexpresados
  text <- c(genesDEup2$Symbol)
  data <- data.frame(x, y, text)
  
  fig <- plot_ly(data, x = ~x, y = ~y, type = 'bar',
                 text = y, textposition = 'auto',
                 marker = list(color = 'rgb(58,10,232)'))
  
  fig <- fig %>% layout(title = "Sobreexpresados",
                        
                        xaxis = list(title = "Simbolo"),
                        
                        yaxis = list(title = "LogFC"))
  
  ###### Guardar la gráfica
  htmlwidgets::saveWidget(fig, "up.html", selfcontained = T, libdir = "lib")
  
  ############## Lo mismo pero para los genes subexpresados
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
  
  #### Gráficas fig tree para los genes subexpresados, usando el valor de p
  genesDEdown100 <- genesDEdown[1:GeneNumber, ]
  tree_down <- genesDEdown100 %>%  
    hchart(type = "treemap", hcaes(x = Symbol, value = adj.P.Val, color = adj.P.Val)) |>
    hc_title(text ="Genes sub expresados organizados de acuerdo a su valor de p") |>
    hc_colorAxis(stops = color_stops(colors = viridis::inferno(10)))
  
  htmlwidgets::saveWidget(tree_down, "TreemapDE_down.html", selfcontained = T, libdir = "lib")
  
  ######### Fig tree para los sobreexpresados
  genesDEup100 <- genesDEup[1:GeneNumber, ]
  tree_up <- genesDEup100 %>%  
    hchart(type = "treemap", hcaes(x = Symbol, value = adj.P.Val, color = adj.P.Val)) |>
    hc_title(text ="Genes sobre expresados organizados de acuerdo a su valor de p") |>
    hc_colorAxis(stops = color_stops(colors = viridis::inferno(10)))
  
  htmlwidgets::saveWidget(tree_up, "TreemapDE_up.html", selfcontained = T, libdir = "lib")
  
  ######## Guardando los genes sobre y subexpresados
  save(genesDE, file = "genesDE.RData")
  save(genesDEup, file = "genesDEup.RData")
  save(genesDEdown, file = "genesDEdown.RData")
}
GeneDE("Datos_Costa_Castelo/ttFIRinELGANsUCtissue.rds", pvalor = 0.01, patron = "TLR", GeneNumber = 10)

######## Analisis GO
GO_analysis <- function(searchFormat, typeProcess ){
  
  serchFormat <- as.character(searchFormat)
  typeProcess <- as.character(typeProcess)
  ####### Cargar las tablas con los sobre y los subexpresados
  load("genesDEup.RData")
  load("genesDEdown.RData")
  
  ####### Tomando solo los símbolos de los genes para realizar la búsqueda en 
  # la base de datos org.Hs.eg.db
  genesDEupSymbol <- genesDEup[ , 1]
  genesDEupSymbol
  
  ###### En este caso se utiliza el formato symbol para los genes, dado por la variable searchFormat
  # En este caso se utiliza biological process, pero puede cambiarse en la variable typeProcess
  GO_up <- enrichGO(gene = genesDEupSymbol, OrgDb = "org.Hs.eg.db", keyType = searchFormat, ont = typeProcess)
  GO_up <- as.data.frame(GO_up)
  
  ###### Para los genes subexpresados
  genesDEdownSymbol <- genesDEdown[ , 1]
  genesDEdownSymbol
  GO_down <- enrichGO(gene = genesDEdownSymbol, OrgDb = "org.Hs.eg.db", keyType = searchFormat, ont = typeProcess)
  GO_down <- as.data.frame(GO_down)
  
  ### Gráficas con las 20 categorías más comunes, para los sobreexpresados, de acuerdo al análisis GO
  png("graphGO_up.png")
  GO_up10 <- GO_up[order(GO_up$Count, decreasing = T), ]
  GO_up10 <- GO_up10[1:10,]
  graphGO_up <- ggplot(GO_up10, aes(x=Description, y=Count, fill = Description)) + 
    geom_bar(stat = "identity") +
    theme(legend.position="none")+ 
    coord_flip()
  graphGO_up
  dev.off()
 
  ### Gráficas con las 20 categorías más comunes, para los subexpresados, de acuerdo al análisis GO 
  png("graphGO_down.png")
  GO_down10 <- GO_down[order(GO_down$Count, decreasing = T), ]
  GO_down10 <- GO_down10[1:10,]
  graphGO_up <- ggplot(GO_down10, aes(x=Description, y=Count, fill = Description)) + 
    geom_bar(stat = "identity") +
    theme(legend.position="none")+ 
    coord_flip()
  graphGO_up
  dev.off()
  
  #### Guardando archivos con genes sub y sobre expresados
  save(GO_up, file = "GO_up.RData")
  save(GO_down, file = "GO_down.RData")
  save(net_down, file = "net_down.RData")
  save(net_up, file = "net_up.RData")
}


############ REDES CON net con networkD3

#### Guardar las categorías encontradas en la base de datos org.Hs.eg.db
dataCategory <- godata("org.Hs.eg.db", ont = "BP")
save(dataCategory, file = "d_up.RData")

networkGO <- function(lim){
    lim <- as.numeric(lim)
    ### Guardamos el objeto GO_down de tipo enrichGO a un data frame
    GO_downT <- GO_down
    GO_downT <- as.data.frame(GO_downT)
    ### Añadimos un _ y quitamos espacios vacíos
    GO_downT$Description <- chartr(old = " ", new = "_", GO_downT$Description)
    
    ### Removemos las categorías que no cumplan con un cierto número de genes 
    # asociados a ella (lim)
    invalid_down <- filter(GO_downT, Count < lim)
    ### nos quedamos con las categorías invalidas
    invalid_down <- invalid_down$Description
    
    ### Matriz de similaridad que utilizaremos para crear nuestra red
    down <- pairwise_termsim(GO_downT, semData = d_down)
    similarity_down <- down@termsim
    networkDown <- data.frame()

    ### Obtner el nombre de cada una de las categorías en la matriz de similaridad
    cat_namesD <- c(colnames(similarity_down))

###### Ciclo para crear las 3 columnas que formaran la matriz de similaridad
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

#### Matriz de adyacencia para crear la red
    networkDown <- data.frame(Source, Target, Value)
    networkDown_final <- networkDown

#### quitando los invalidos
    num <- 1
      for(i in 1:length(invalid_down)){
        networkDown_final <- filter(networkDown_final, Source != invalid_down[num])
        num <- num + 1
}

#### quitando NAs
    networkDown_final <- filter(networkDown_final, Value != "NA")
    networkDown_final

    netD <- networkDown_final[order(networkDown_final$Value, decreasing = T), ]
    net <- networkDown_final
    netD100 <- netD[1:100, ] ### recuerda que debes correr de nuevo netD, para que se genere una vez mas la tabla


Net_down <- simpleNetwork(netD100, width = 400, height = 250, fontSize = 14, nodeColour ="blue", zoom = T)
saveNetwork(Net_down, file = "Net_down.html", selfcontained = T)

########################## PARA LOS UP
GO_upT <- as.data.frame(GO_up)
GO_upT$Description <- chartr(old = " ", new = "_", GO_upT$Description)

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


############################33
#############################
#   FORCE NETWORK
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
