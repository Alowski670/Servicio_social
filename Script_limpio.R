
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
library(clusterProfiler)
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

GO_analysis <- function(GO_ID){
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
  
  ##### Nuevas objetos para trabajar
  blank_up <- subset(GO_up, select = c(Description, geneID))
  blank_down <-subset(GO_down, select = c(Description, geneID))
  
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
      net_up$Target[to1:to2] <- (GO_up$Description[num1])
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
      net_down$Target[to1D:to2D] <- (GO_down$Description[num1D])
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
      net_up$Target[to1:to2] <- (GO_up$ID[num1])
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
      net_down$Target[to1D:to2D] <- (GO_down$ID[num1D])
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

################ Redes con D3network
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
    saveNetwork(net_up10, file = "net_down100.html", selfcontained = T)
    
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
    saveNetwork(net_down10, file = "net_down100.html", selfcontained = T)
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
    saveNetwork(net_up10, file = "net_down100.html", selfcontained = T)
    
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
    saveNetwork(net_down10, file = "net_down100.html", selfcontained = T)
  }
}
D3net(GO_ID = F)
