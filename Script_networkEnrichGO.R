###################### A PARTIR DE AQUI EN TEORIA FUNCIONA
### para los down
library(dplyr)
library(enrichplot)
library(GOSemSim)
library(networkD3)
library(clusterProfiler)
###### cortar la red, para asi poder ver como interactuan entre si
GO_downT <- as.data.frame(GO_down)
GO_downT <- GO_downT@result[order(GO_downT@result$Count, decreasing = T),]
GO_downT <- filter(GO_downT, Count > 20)
GO_downT <- enrichGO(GO_downT)
View(GO_downT@result)

invalid_down <- filter(GO_downT, Count < 20)
invalid_down <- invalid_down$Description


d_down <- godata("org.Hs.eg.db", ont = "BP")
down <- pairwise_termsim(GO_downT, semData = d_down)
similarity_down <- down@termsim
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

Source
Target
Value
networkDown <- data.frame(Source, Target, Value)
networkDown_final <- networkDown

#### quitando los invalidos
num <- 1
for(i in 1:length(invalid_down)){
  networkDown_final <- filter(networkDown_final, Source != invalid_down[num])
  num <- num + 1
}

#### quitando NAs
#networkDown_final$Value[networkDown_final$Value <= 0.05] <- "NA"
networkDown_final <- filter(networkDown_final, Value != "NA")

View(networkDown_final)
netD <- networkDown_final[order(networkDown_final$Value, decreasing = T), ]
net <- networkDown_final
View(netD)
netD <- netD[1:100, ] ### recuerda que debes correr de nuevo netD, para que se genere una vez mas la tabla
net <- net[1:300,]

simpleNetwork(netD, width = 400, height = 250, fontSize = 11, nodeColour = "orange", zoom = T)

#### con enrich plot
emapplot(down, showCategory = 10)
cnetplot(down)


########################## PARA LOS UP
GO_upT <- as.data.frame(GO_up)
GO_upT <- GO_upT@result[order(GO_upT@result$Count, decreasing = T),]
GO_upT <- filter(GO_upT, Count > 20)
GO_upT <- enrichGO(GO_upT)

invalid_up <- filter(GO_upT, Count < 20)
invalid_up <- invalid_up$Description

d_up <- godata("org.Hs.eg.db", ont = "BP")
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
netD_UP <- netD_UP[1:20, ]
net_UP <- net_UP[1:300,]
#### funciona ahora añadir los conteos a las categorias

simpleNetwork(netD_UP, width = 400, height = 250, fontSize = 11, nodeColour = "orange", zoom = T)
emapplot(up, showCategory = 10)
cnetplot(up)