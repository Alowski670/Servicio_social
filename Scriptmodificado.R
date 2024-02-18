library(tools)
library(clusterProfiler)
library(igraph)

### funcion para leer archivo
uRead <- function(file){
  ext <- file_ext(file)
  print(paste0("Es un archivo tipo: ",ext))
  print(paste0("Links con funciones para leer archivos:"))
  print(paste0("https://bookdown.org/rdpeng/rprogdatascience/getting-data-in-and-out-of-r.html"))
  print(paste0("https://bookdown.org/rdpeng/rprogdatascience/using-the-readr-package.html"))
}
uRead("Datos_Costa_Castelo/ttFIRinELGANsUCtissue.rds")
lista <- readRDS("Datos_Costa_Castelo/ttFIRinELGANsUCtissue.rds")

tablaGO2 <- function(lista, pvalor, patron){
  list_gene <-readRDS(lista)
  print("Dimensiones")
  print(dim(list_gene))
  head(list_gene, n= 3)
  
  ## Genes diferencialmente expresados
  #pvalor <- readline(prompt = "Indique el valor de p deseado: ")
  #pvalor <- as.numeric(pvalor)
  genesDE <- list_gene[list_gene$adj.P.Val < pvalor & 
                         abs(list_gene$logFC) > log2(1.5), ]
  
  ## Patrón de búsqueda
  #patron <- readline(prompt = "Indique el símbolo de búsqueda: ")
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
  tabla_DT
  ### Tablas interactivas con kableextra
  genesDEup <- filter(genesDE, logFC > 0)
  genesDEup2 <- genesDEup[1:10, ]
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
  genesDEdown2 <- genesDEdown[1:10, ]
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
  
  print(fig)
  print(fig2)
}
tablaGO2("Datos_Costa_Castelo/ttFIRinELGANsUCtissue.rds", pvalor = 0.01, patron = "TLR")

### fig tree 
View(lista)
#install.packages("highcharter")
#install.packages("htmltools")
#update.packages(ask = F)
library(highcharter)

View(genesDEdown)
help(hchart)
genesDEdown %>% 
hchart(type = "treemap", hcaes(x = Symbol, value = adj.P.Val, color = adj.P.Val)) |>
  hc_title(text ="Genes sub expresados organizados de acuerdo a valor de P")

genesDEup %>% 
  hchart(type = "treemap", hcaes(x = Symbol, value = adj.P.Val, color = adj.P.Val)) |>
  hc_title(text ="Genes sobre expresados organizados de acuerdo a valor de P")

################## Analisis GO
### instalando paquetes
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("clusterProfiler")

#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("org.Hs.eg.db", force = TRUE)

#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("AnnotationDbi")

#### para la red, asociar los simbolos que aparecen con la categoria, osea el numero de genes sobre o sub por categoria go
############ Analisis sobre y sub
genesDEupSymbol <- genesDEup[ , 1]
genesDEupSymbol
GO_sobre <- enrichGO(gene = genesDEupSymbol, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO_sobre <- as.data.frame(GO_sobre)

genesDEdownSymbol <- genesDEdown[ , 1]
genesDEdownSymbol
GO_sub <- enrichGO(gene = genesDEdownSymbol, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO_sub <- as.data.frame(GO_sub)
### grafica
fig_sobre <- plot(barplot(GO_sobre, showCategory = 20))
fig_sobre
View(GO_sobre)

fig_sub <- plot(barplot(GO_sub, showCategory = 20))
fig_sub

### ya tiene los conteos, ahora queda hacer un txt con las categorias que aparecen y los conteos para cada una, identificando a que categoria pertenece cada gen
red_sobre <- GO_sobre$Description
red_sobre
View(GO_sobre)
head(GO_sobre$geneID)

#### solo nos quedamos con los nombres de los genes para cada renglon
blank_up <- subset(GO_sobre, select = c(Description, geneID))
#View(blank_up)
#blank_up$Description_clean[2] <- chartr(old = "/", new = " ", blank_up$geneID[2])
#View(blank_up)

### Para los genes sub expresados
blank_down <- subset(GO_sub, select = c(Description, geneID))
##### funcion para limpiar los genes asociados a cada categoria, remover el "/"

blank_up <- subset(GO_sobre, select = c(Description, geneID))
#View(blank_up)

blank_down <-subset(GO_sub, select = c(Description, geneID))
View(blank_down)
### genes sobre
cont <- 0
for(i in 1:dim(blank_up[1])){
  cont <- cont + 1
    blank_up$Description_clean[cont] <- chartr(old = "/", new = " ", blank_up$geneID[cont])
    
}
View(blank_up)

### sub
contD <- 0
for(i in 1:dim(blank_down[1])){
  contD <- contD + 1
  blank_down$Description_clean[contD] <- chartr(old = "/", new = " ", blank_down$geneID[contD])
  
}
blank_down$Description_clean[1]
View(blank_down)

### Calculando la cantidad de genes asociados a cada categoria
str_length <- c()
cont2 <- 0
for(i in 1:dim(blank_up[1])){
  cont2 <- cont2 + 1
  str_length <- c(str_length, lengths(strsplit(blank_up$Description_clean[cont2], " ")))
}

head(str_length)
#suma <- sum(str_length)
#suma2 <- abs(str_length[1] - suma)
#suma - suma2 
#### Preparando la tabla para la red D3network
net_up <- c()
net_up$Source <- unlist(strsplit(blank_up$Description_clean, split = " "))
net_up <- as.data.frame(net_up)
View(net_up)

#### Asignando las categorias 
#to1 <- 1
#to2 <- str_length[1]
#net_up$Target <- "NO"
#str_length[3]

#net_up$Target[to1:to2] <- (GO_sobre$Description[1])
#to1 <- to2 + 1 #79, 149
#to2 <- str_length[2] + str_length[3] #148
#to2 <- sum(str_length[1:num2])
#View(net_up)

to1 <- 1
to2 <- str_length[1]
num1 <- 1
num2 <- 2
net_up$Target <- "NO"
net_up

for(i in 1:dim(blank_up[1])){
net_up$Target[to1:to2] <- (GO_sobre$Description[num1])
to1 <- to2 + 1 #79, 149
to2 <- sum(str_length[1:num2]) #148, 210
num1 <- num1 + 1 #2, 3
num2 <- num2 + 1 #3, 4
}
View(net_up)

##### Genes sub
str_lengthD <- c()
cont2D <- 0
for(i in 1:dim(blank_down[1])){
  cont2D <- cont2D + 1
  str_lengthD <- c(str_lengthD, lengths(strsplit(blank_down$Description_clean[cont2D], " ")))
}

head(str_lengthD)
net_down <- c()
net_down$Source <- unlist(strsplit(blank_down$Description_clean, split = " "))
net_down <- as.data.frame(net_up)
View(net_down)

to1D <- 1
to2D <- str_lengthD[1]
num1D <- 1
num2D <- 2
net_down$Target <- "NO"
net_down

for(i in 1:dim(blank_down[1])){
  net_down$Target[to1D:to2D] <- (GO_sub$Description[num1D])
  to1D <- to2D + 1 #79, 149
  to2D <- sum(str_lengthD[1:num2D]) #148, 210
  num1D <- num1D + 1 #2, 3
  num2D <- num2D + 1 #3, 4
}
View(net_down)
#################3 Red networkD3
#install.packages("networkD3")
## las 100 categorias mas populares
top100 <- GO_sobre[order(GO_sobre$Count, decreasing = T),]
top100_names <- top100$Description[1:100]
top100_names <- as.character(top100_names)
top100_names

###### las 10 mas
top10 <- GO_sobre[order(GO_sobre$Count, decreasing = T),]
top10_names <- top100$Description[1:10]
top10_names <- as.character(top10_names)
top10_names

library(dplyr)
### 100 mas
net_up100 <- net_up %>% filter(Target %in% c(top100_names))
net_up100

### 10 mas
net_up10 <- net_up %>% filter(Target %in% c(top10_names))
net_up10
dim(net_up100)
#netup2 <- net_up[1:100, ]
View(net_up)
library(networkD3)
simpleNetwork(netup2, width = 400, height = 250, fontSize = 11, nodeColour = "orange", zoom = T)

#### con igraph
#library(igraph)
#node_names <-
redi_up <- graph_from_data_frame()

######################### Red sub
top100D <- GO_sub[order(GO_sub$Count, decreasing = T),]
top100D_names <- top100D$Description[1:100]
top100D_names <- as.character(top100D_names)
top100D_names

###### las 10 mas
top10D <- GO_sub[order(GO_sub$Count, decreasing = T),]
top10D_names <- top100D$Description[1:10]
top10D_names <- as.character(top10D_names)
top10D_names

library(dplyr)
### 100 mas
net_down100 <- net_down %>% filter(Target %in% c(top100D_names))
net_down100

### 10 mas
net_down10 <- net_down %>% filter(Target %in% c(top10D_names))
net_down10
dim(net_down100)
#netup2 <- net_up[1:100, ]
View(net_down10)
library(networkD3)
simpleNetwork(net_down10, width = 400, height = 250, fontSize = 11, nodeColour = "orange", zoom = T)

#### con igraph
#library(igraph)
#node_names <-
redi_up <- graph_from_data_frame()



GO_analysis <- function(GO_ID){
  load("genesDEup.RData")
  load("genesDEdown.RData")
  #### Analisis GO
  genesDEupSymbol <- genesDEup[ , 1]
  genesDEupSymbol
  GO_up <- enrichGO(gene = genesDEupSymbol, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                    qvalueCutoff = 0.2, minGSSize = 20, maxGSSize = 500)
  GO_upT <- as.data.frame(GO_up)
  
  genesDEdownSymbol <- genesDEdown[ , 1]
  genesDEdownSymbol
  GO_down <- enrichGO(gene = genesDEdownSymbol, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                      qvalueCutoff = 0.2, minGSSize = 20, maxGSSize = 500)
  GO_downT <- as.data.frame(GO_down)
  
  ### Gráficas con las 20 categorías más comunes, de acuerdo al análisis GO
  png("graphGO_up.png")
  GO_up10 <- GO_upT[order(GO_up$Count, decreasing = T), ]
  GO_up10 <- GO_up10[1:10,]
  graphGO_up <- ggplot(GO_up10, aes(x=Description, y=Count, fill = Description)) + 
    geom_bar(stat = "identity") +
    theme(legend.position="none")+ 
    coord_flip()
  graphGO_up
  dev.off()
  
  png("graphGO_down.png")
  GO_down10 <- GO_downT[order(GO_downT$Count, decreasing = T), ]
  GO_down10 <- GO_down10[1:10,]
  graphGO_up <- ggplot(GO_down10, aes(x=Description, y=Count, fill = Description)) + 
    geom_bar(stat = "identity") +
    theme(legend.position="none")+ 
    coord_flip()
  graphGO_up
  dev.off()
  
  ##### Nuevas objetos para trabajar
  blank_up <- subset(GO_upT, select = c(Description, geneID))
  blank_down <-subset(GO_downT, select = c(Description, geneID))
  
  ### Quitando el "/" y reemplazándolo por espacios en blanco
  # SOBRE
  cont <- 0
  for(i in 1:dim(blank_up[1])){
    cont <- cont + 1
    blank_up$Description_clean[cont] <- chartr(old = "/", new = " ", blank_up$geneID[cont])
    
  }
  
  # SUB
  contD <- 0
  for(i in 1:dim(blank_down[1])){
    contD <- contD + 1
    blank_down$Description_clean[contD] <- chartr(old = "/", new = " ", blank_down$geneID[contD])
    
  }
  
  ### Calculando la cantidad de genes por cada categoría
  # SOBRE
  str_length <- c()
  cont2 <- 0
  for(i in 1:dim(blank_up[1])){
    cont2 <- cont2 + 1
    str_length <- c(str_length, lengths(strsplit(blank_up$Description_clean[cont2], " ")))
  }
  
  if(GO_ID == F){
    net_up <- c()
    net_up$Source <- unlist(strsplit(blank_up$Description_clean, split = " "))
    net_up <- as.data.frame(net_up)
  
      to1 <- 1
      to2 <- str_length[1]
      num1 <- 1
      num2 <- 2
      net_up$Target <- "NO"
  
    for(i in 1:dim(blank_up[1])){
      net_up$Target[to1:to2] <- (GO_upT$Description[num1])
      to1 <- to2 + 1 #79, 149
      to2 <- sum(str_length[1:num2]) #148, 210
      num1 <- num1 + 1 #2, 3
      num2 <- num2 + 1 #3, 4
    }
  
    # SUB
    str_lengthD <- c()
    cont2D <- 0
    for(i in 1:dim(blank_down[1])){
      cont2D <- cont2D + 1
      str_lengthD <- c(str_lengthD, lengths(strsplit(blank_down$Description_clean[cont2D], " ")))
    }
  
    head(str_lengthD)
    net_down <- c()
    net_down$Source <- unlist(strsplit(blank_down$Description_clean, split = " "))
    net_down <- as.data.frame(net_up)
  
    to1D <- 1
    to2D <- str_lengthD[1]
    num1D <- 1
    num2D <- 2
    net_down$Target <- "NO"
    net_down
  
    for(i in 1:dim(blank_down[1])){
      net_down$Target[to1D:to2D] <- (GO_downT$Description[num1D])
      to1D <- to2D + 1 #79, 149
      to2D <- sum(str_lengthD[1:num2D]) #148, 210
      num1D <- num1D + 1 #2, 3
      num2D <- num2D + 1 #3, 4
    }
  } else{
    net_up <- c()
    net_up$Source <- unlist(strsplit(blank_up$Description_clean, split = " "))
    net_up <- as.data.frame(net_up)
    
    to1 <- 1
    to2 <- str_length[1]
    num1 <- 1
    num2 <- 2
    net_up$Target <- "NO"
    
    for(i in 1:dim(blank_up[1])){
      net_up$Target[to1:to2] <- (GO_upT$ID[num1])
      to1 <- to2 + 1 #79, 149
      to2 <- sum(str_length[1:num2]) #148, 210
      num1 <- num1 + 1 #2, 3
      num2 <- num2 + 1 #3, 4
    }
    
    # SUB
    str_lengthD <- c()
    cont2D <- 0
    for(i in 1:dim(blank_down[1])){
      cont2D <- cont2D + 1
      str_lengthD <- c(str_lengthD, lengths(strsplit(blank_down$Description_clean[cont2D], " ")))
    }
    
    head(str_lengthD)
    net_down <- c()
    net_down$Source <- unlist(strsplit(blank_down$Description_clean, split = " "))
    net_down <- as.data.frame(net_up)
    
    to1D <- 1
    to2D <- str_lengthD[1]
    num1D <- 1
    num2D <- 2
    net_down$Target <- "NO"
    net_down
    
    for(i in 1:dim(blank_down[1])){
      net_down$Target[to1D:to2D] <- (GO_downT$ID[num1D])
      to1D <- to2D + 1 #79, 149
      to2D <- sum(str_lengthD[1:num2D]) #148, 210
      num1D <- num1D + 1 #2, 3
      num2D <- num2D + 1 #3, 4
    }
  }
  save(GO_up, file = "GO_up.RData")
  save(GO_down, file = "GO_down.RData")
  save(net_down, file = "net_down.RData")
  save(net_up, file = "net_up.RData")
}
GO_analysis(GO_ID = F)

View(GO_up)
View(net_up)

load("GO_up.RData")
load("GO_down.RData")
load("net_down.RData")
load("net_up.RData")

### Red con D3network
D3net <- function(GO_ID){
load("GO_up.RData")
load("GO_down.RData")
load("net_down.RData")
load("net_up.RData")
if(GO_ID == F){
  top100 <- GO_up[order(GO_up$Count, decreasing = T),]
  top100_names <- top100$Description[1:100]
  top100_names <- as.character(top100_names)
  top100_names

  ###### las 10 mas
  top10 <- GO_up[order(GO_up$Count, decreasing = T),]
  top10_names <- top100$Description[1:10]
  top10_names <- as.character(top10_names)
  top10_names

  ### 100 mas
  #net_up100 <- net_up %>% filter(Target %in% c(top100_names))
  #net_up100

  ### 10 mas
  net_up10 <- net_up %>% filter(Target %in% c(top10_names))
  net_up10
  simpleNetwork(net_up10, width = 400, height = 250, fontSize = 11, nodeColour = "orange", zoom = T)

  ######################### Red sub
  top100D <- GO_down[order(GO_down$Count, decreasing = T),]
  top100D_names <- top100D$Description[1:100]
  top100D_names <- as.character(top100D_names)
  top100D_names

  ###### las 10 menos
  top10D <- GO_down[order(GO_down$Count, decreasing = T),]
  top10D_names <- top100D$Description[1:10]
  top10D_names <- as.character(top10D_names)
  top10D_names

  ### 100 menos
  net_down100 <- net_down %>% filter(Target %in% c(top100D_names))
  #simpleNetwork(net_down100, width = 400, height = 250, fontSize = 11, nodeColour = "orange", zoom = T) 
  #saveNetwork(net_down100, file = "net_down100.html", selfcontained = T)

  ### 10 menos
  net_down10 <- net_down %>% filter(Target %in% c(top10D_names))
  simpleNetwork(net_down10, width = 400, height = 250, fontSize = 11, nodeColour = "blue", zoom = T)
  #save(net_down10)
}else{
  top100 <- GO_up[order(GO_up$Count, decreasing = T),]
  top100_names <- top100$ID[1:100]
  top100_names <- as.character(top100_names)
  top100_names
  
  ###### las 10 mas
  top10 <- GO_up[order(GO_up$Count, decreasing = T),]
  top10_names <- top100$ID[1:10]
  top10_names <- as.character(top10_names)
  top10_names
  
  ### 100 mas
  net_up100 <- net_up %>% filter(Target %in% c(top100_names))
  net_up100
  
  ### 10 mas
  net_up10 <- net_up %>% filter(Target %in% c(top10_names))
  net_up10
  simpleNetwork(net_up10, width = 400, height = 250, fontSize = 11, nodeColour = "orange", zoom = T)
  
  ######################### Red sub
  top100D <- GO_down[order(GO_down$Count, decreasing = T),]
  top100D_names <- top100D$ID[1:100]
  top100D_names <- as.character(top100D_names)
  top100D_names
  
  ###### las 10 menos
  top10D <- GO_down[order(GO_down$Count, decreasing = T),]
  top10D_names <- top100D$ID[1:10]
  top10D_names <- as.character(top10D_names)
  top10D_names
  
  ### 100 menos
  #net_down100 <- net_down %>% filter(Target %in% c(top100D_names))
  #simpleNetwork(net_down100, width = 400, height = 250, fontSize = 11, nodeColour = "orange", zoom = T) 
  #saveNetwork(net_down100, file = "net_down100.html", selfcontained = T)
  
  ### 10 menos
  net_down10 <- net_down %>% filter(Target %in% c(top10D_names))
  simpleNetwork(net_down10, width = 400, height = 250, fontSize = 11, nodeColour = "blue", zoom = T)
  #save(net_down10)
  }
}
D3net(GO_ID = T)

##############################################333
########## RED CON LAS CATEGORIAS
library(clusterProfiler)
emapplot(GO_up)

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("enrichplot", force = T)

library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)

d <- godata("org.Hs.eg.db", ont = "BP")
arriba <- pairwise_termsim(GO_up, semData = d)
View(arriba@termsim) ### para ver similaridad de terminos
emapplot(arriba)
View(emapplot)

library(networkD3)

### sacar la matriz de similaridad
simpleNetwork(arriba@termsim)

Source <- colnames(arriba@termsim)
Source
Target <- colnames(similarity)
Value <- 1:dim(similarity)[1]
networkUP <- data.frame(Source, Target, Value )
#networkUP <- ## rep los terminos de source para crear la columna source

### ciclo for para sacar el valor de veces que se tiene que repetir cada categoria para el nuevo data frame
#similarity <- arriba@termsim
#select_row <- 1 
#category_number <- c()
#for(i in 1:dim(similarity)[1]- 1){
  #category_value <-  which(similarity[select_row, ] != "NA")
  #category_number <- c(category_number, length(category_value))
  #select_row <- select_row + 1
#}
#category_number

num1 <- 1
for(i in 1:dim(similarity)[1]){
  networkUP$Source <- rep(Source[num1], times = 1, each = category_number[num1])
  num1 <- num1 + 1
}

k <- rep(Source[1], times = 1, each = category_number[1])
k

View(networkUP)

similarity <- arriba@termsim
select_row <- 1 
category_number <- c()
for(i in 1:dim(similarity)[1]- 1){
category_value <-  similarity[select_row, ]
category_number <- c(category_number, length(category_value))
select_row <- select_row + 1
}
category_number

#### otra idea
Source <- c(1:dim(similarity)[1]*dim(similarity)[1])
Source
networkUP <- data.frame(Source)
View(networkUP)

cat_names <- c(colnames(similarity))
cat_names
num <- 1
for(i in 1:dim(similarity)[1]){
  networkUP$Source <- rep(cat_names[num], times = 1, each = 1:dim(similarity)[1])
  num <- num + 1
}

############### metodo con matriz similaridad
similarity <- arriba@termsim 
networkUP <- data.frame()
cat_names <- c(colnames(similarity))
cat_names

num <- 1
Source <- c()
Target <- c()
Value <- c()
for(i in 1:dim(similarity)[1]){
  Source <- c(Source, rep(cat_names[num], times = dim(similarity)[1], each = 1:dim(similarity)[1]))
  Target <- c(Target, rep(cat_names, times = 1, each = 1))
  Value <- c(Value, similarity[num, ])
  num <- num + 1
  
}
Source
Target
Value
networkUP <- data.frame(Source, Target, Value)
View(networkUP)

rep(cat_names, times = 1, each = 1)

networkUP_final <- na.omit(networkUP)
View(networkUP_final)

library(networkD3)
simpleNetwork(networkUP_final, width = 400, height = 250, fontSize = 11, nodeColour = "orange", zoom = T)

###################### A PARTIR DE AQUI EN TEORIA FUNCIONA
### para los down
library(dplyr)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
library(GOSemSim)
###### cortar la red, para asi poder ver como interactuan entre si
GO_downT <- as.data.frame(GO_down)
GO_downT <- GO_downT@result[order(GO_downT@result$Count, decreasing = T),]
GO_downT <- filter(GO_downT, Count > 9)
GO_downT <- enrichGO(GO_downT)
View(GO_downT@result)

d <- godata("org.Hs.eg.db", ont = "BP")
abajo <- pairwise_termsim(GO_downT, semData = d)
similarity_down <- abajo@termsim
View(similarity_down)
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

#### Modificando Source y target para quedarnos con las categorias con el mayor count
Source
Target
Value
networkDown <- data.frame(Source, Target, Value)
networkDown_final <- networkDown
View(networkDown)
#### quitando los invalidos
num <- 1
  for(i in 1:length(invalid)){
    networkDown_final <- filter(networkDown_final, Source != invalid[num])
    num <- num + 1
  }
invalid

#### quitando NAs
#networkDown_final$Value[networkDown_final$Value <= 0.05] <- "NA"
networkDown_final <- filter(networkDown_final, Value != "NA")

View(networkDown_final)
netD <- networkDown_final[order(networkDown_final$Value, decreasing = T), ]
net <- networkDown_final
View(netD)
netD <- netD[1:20, ]
net <- net[1:300,]
#### funciona ahora añadir los conteos a las categorias

library(networkD3)
library(enrichplot)
simpleNetwork(netD, width = 400, height = 250, fontSize = 11, nodeColour = "orange", zoom = T)
emapplot(abajo, showCategory = 10)
cnetplot(abajo)
View(abajo@termsim)
### cambiar valores de p o q

library(igraph)
library(dplyr)
help("emapplot,compareClusterResult-method")
x <- graph_from_adjacency_matrix(similarity_down )
plot(x)

### quitando categorias que no cumplan con cierto numero de counts
View(GO_upT)
invalid <- filter(GO_downT, Count < 20)
invalid <- invalid$Description
invalid
head(invalid)
View(similarity_downFinal)
similarity_downFinal <- similarity_down
similarity_downFinal <- which(similarity_down != "striated muscle contraction")
############ enrich plot se ve interesante, revisar

####################################################
################# Red con las categorias 2do intento
