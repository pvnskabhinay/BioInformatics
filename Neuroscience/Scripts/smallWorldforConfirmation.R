install.packages("qgraph")
library(qgraph)
library(igraph)

gdata = read.csv("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\edgeListSmallWorld.csv")

colnames(gdata) <- c("Source","Destination")
head(gdata)
g1 <- qgraph(gdata)
smallworldIndex(g1)

smallworldness(g1)
