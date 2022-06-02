library(reshape2)
library(tidyr)
library(readr)
library(igraph)

setwd("../air/")
air = read_delim("China_2015/data15.csv", ";", escape_double = FALSE, trim_ws = TRUE)
cityids = unique(air$citycode)  # 338 cities
airpos = read.csv("China_2015/location15.csv", fileEncoding = "GB18030")
# distance matrix
d.mat = dist(as.matrix(airpos[,3:4]), diag = TRUE, upper = T) %>% as.matrix()
cat("max distance is", max(d.mat))
summary(c(d.mat))

# construct Amat by Wang et al. (2020)

airpos$lat = airpos$lat * 110.574  # convert to km
airpos$lon = airpos$lon * 111.320  # convert to km

d.mat = dist(as.matrix(airpos[,3:4]), diag = TRUE, upper = T) %>% as.matrix()
cat("max distance is", max(d.mat))
summary(c(d.mat))

adjmat = matrix(0, nrow = nrow(d.mat), ncol = ncol(d.mat))
adjmat[which(d.mat < 300)] <- 1
diag(adjmat) <- 0
cat("network density is", sum(adjmat)/(nrow(adjmat)^2 - nrow(adjmat)), "\n")
g0 <- graph_from_adjacency_matrix(adjmat, mode = "undirected")
table(degree(g0))
components <- components(g0, mode="weak")
components$csize
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(g0)[components$membership %in% biggest_cluster_id]
ids1 = unique(airpos$citycode)[vert_ids]
adjmat1 = adjmat[vert_ids,vert_ids]
d.mat1 = d.mat[vert_ids, vert_ids]
g1 <- induced.subgraph(g0, vert_ids)  # 310 nodes

air1 = air[which(air$citycode %in% ids1),c("citycode", "pm25", "DATES")]  # 选择不同的气体
airwide = air1 %>% 
  mutate(DATES = factor(DATES, levels = unique(DATES))) %>%
  spread(DATES, o3)

# air1$citycode = NULL

Ymat = as.matrix(airwide[,-1])
rownames(Ymat) = airwide$citycode
save(Ymat, file = "./air_data/Ymat1_o3.rda")
save(adjmat1, file = "./air_data/Adjmat1_o3.rda")
save(d.mat1, file = "./air_data/dmat1_o3.rda")
save(g1, file = "./air_data/g_air1_o3.rda")
# save(d.matrix, file = "./air_data/dmat_path.rda")



# 3. adj matrix
load("./China_2015/spmat.rda")
dijishiID <- read.csv("./China_2015/dijishiID.csv", header = T)
cityCODE <- dijishiID[which(dijishiID$市代码 %in% cityids), c("市代码" ,"市")]  # 314 cities
cityCODE_2019 <- as.character(cityCODE$市代码)
ids <- intersect(rownames(spmat2), cityCODE_2019)  # 313 cities
adjmat <- spmat2[which(rownames(spmat2) %in% ids),which(rownames(spmat2) %in% ids)]
g0 <- graph_from_adjacency_matrix(adjmat, mode = "undirected")
table(degree(g0))
g1 <- induced.subgraph(graph = g0, V(g0)[which(degree(g0) > 0)])  # 308 cities
plot(g1, vertex.size=5,vertex.color=c("orange", "steelblue","gray", "pink","dark red"), vertex.label = NA)
ids1 = names(V(g0))[which(degree(g0) > 0)]
ids2 = intersect(ids1, as.character(cityids))  # 308 cities
adjmat1 = adjmat[ids2, ids2]
d.matrix <- distances(g1, v = V(g1), to = V(g1))
cat("max distance is", max(d.mat))


air1 = air[which(as.character(air$citycode) %in% ids2),c("citycode", "pm25", "DATES")]
airwide = spread(air1, DATES, pm25)
Ymat = as.matrix(airwide[,-1])
rownames(Ymat) = airwide$citycode


save(Ymat, file = "./air_data/Ymat.rda")
save(adjmat1, file = "./air_data/Adjmat.rda")
save(d.mat, file = "./air_data/dmat.rda")
save(g1, file = "./air_data/g_air.rda")
save(d.matrix, file = "./air_data/dmat_path.rda")
