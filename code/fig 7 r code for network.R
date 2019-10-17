library(dplyr); library("psych");library(igraph);library(reshape2)

MCA=read.delim("MCA2.csv", header=T, sep=",")
cor <- cor(MCA[,-1])

cor.test.p <- function(x){
    FUN <- function(x, y) cor.test(x, y)[["p.value"]]
    z <- outer(
      colnames(x), 
      colnames(x), 
      Vectorize(function(i,j) FUN(x[,i], x[,j]))
    )
    dimnames(z) <- list(colnames(x), colnames(x))
    z
}
pval=cor.test.p(MCA[,-1])


flattenCorrMatrix <- function(cor, pval) {
  ut <- upper.tri(cor)
  data.frame(
    row = rownames(cor)[row(cor)[ut]],
    column = rownames(cor)[col(cor)[ut]],
    cor  =(cor)[ut],
    p = pval[ut]
    )
}

four=flattenCorrMatrix(cor,pval)
dim(four)
summary(four$cor)

adjlist <- four %>% filter(abs(cor) > 0.6)
dim(adjlist)
class(adjlist)

adjlistP <- adjlist %>% filter(abs(p) < 0.05 )
dim(adjlistP)
adjlistP

net <- graph.data.frame(adjlistP, directed = FALSE)
net
orig_mar <- par()$mar
#par(mar=rep(0.1, 0.1))
names(adjlistP) <- c('from', 'to', 'weight','P')
head(adjlistP)
set.seed(1)
ceb <- cluster_edge_betweenness(net)
plot(ceb, net,vertex.size=betweenness(net)/100) #Betweenness centrality quantifies the number of times a node acts as a bridge along the shortest path between two other nodes.
membership(ceb)


V(net)$size <- V(net)$audience.size*0.7
