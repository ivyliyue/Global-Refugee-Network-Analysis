install.packages("ggplot2")
install.packages("sand")
install.packages("igraph")
install.packages("intergraph")
install.packages("lpSolve")
install.packages("mcclust")
install.packages("psych")
install.packages("GGally")

library(ggplot2)
library(sand)
library(igraph)
library(intergraph)
library(lpSolve)
library(mcclust)
library(psych)
library(GGally)

library(dplyr)
library(reshape2)
library(ggplot2)
library(network)
library(ergm)
library(latticeExtra)
library(intergraph)
library(sand)

path <- "/Users/ivyli66/Desktop/CSC495FinalNew"
setwd(path)

nodes2016 <- read.csv("nodelist_2016_gdp.csv", stringsAsFactors = FALSE)
edges2016 <- read.csv("edgelist_2016_gdp.csv", stringsAsFactors = FALSE)

gr2016 <- graph_from_data_frame(edges2016, vertices = nodes2016, directed =TRUE)
summary(gr2016)



nodes2015 <- read.csv("nodelist_2015_GDP.csv", stringsAsFactors = FALSE)
edges2015 <- read.csv("edgelist_2015_gdp.csv", stringsAsFactors = FALSE)
gr2015 <- graph_from_data_frame(edges2015, vertices = nodes2015, directed =TRUE)
summary(gr2015)
is_simple(gr2015)

comp <- decompose(gr2015)
sa <- sapply(comp, vcount)
sa

#compute 2016 weighted degree

gr2016.wdeg = graph.strength(gr2016, 
                             weights=E(gr2016)$weight)
summary(gr2016.wdeg)

#filter out less than 1st quatile weighted degree 2016 (3008)
to_remove <- V(gr2016)[gr2016.wdeg < 3008]
gr2016.f <- delete_vertices(gr2016, to_remove)
summary(gr2016.f)
plot(gr2016.f)

#compute filtered weighted degree 2016
gr2016.f.wdeg = graph.strength(gr2016.f,
                               weights = E(gr2016.f)$weight)
summary(gr2016.f.wdeg)
summary(gr2016.wdeg)

#compute 2015 weighted degree

gr2015.wdeg = graph.strength(gr2015, 
                             weights = E(gr2015)$weight)
summary(gr2015.wdeg)

#filter out less than 1st quatile weighted degree 2015
to_remove <- V(gr2015)[gr2015.wdeg < 1804]
gr2015.f <- delete_vertices(gr2015, to_remove)
plot(gr2015.f)  
summary(gr2015.f)
#compare 2015 filter
summary(gr2015)
summary(gr2015.f)
gr2015.f <- delete_vertex_attr(gr2015.f, "X.1")
gr2015.f <- delete_vertex_attr(gr2015.f, "X.2")
gr2015.f <- delete_vertex_attr(gr2015.f, "X.3")
summary(gr2015.f)

#compute 2015 filtered weighted degree
gr2015.f.wdeg = graph.strength(gr2015.f,
                               weights = E(gr2015.f)$weight)
summary(gr2015.f.wdeg)
summary(gr2015.wdeg)

#visulized GDP

summary(V(gr2015.f)$GDP)
summary(V(gr2016.f)$GDP)

scatter.smooth(V(gr2015.f.GDP)$GDP)
hist(V(gr2015.f.GDP)$GDP)

#filter GDP

gr2015.f.GDP <- delete_vertices(gr2015.f, V(gr2015.f)$GDP ==0)
summary(V(gr2015.f.GDP)$GDP)

gr2016.f.GDP <- delete_vertices(gr2016.f, V(gr2016.f)$GDP ==0)
summary(V(gr2016.f.GDP)$GDP)

gr2015.f.GDP.wdeg = graph.strength(gr2015.f.GDP,
                                   weights = E(gr2015.f.GDP)$weight)
summary(gr2015.f.GDP.wdeg)
summary(gr2015.f.GDP)

gr2016.f.GDP.wdeg = graph.strength(gr2016.f.GDP,
                                   weights = E(gr2016.f.GDP)$weight)
summary(gr2016.f.GDP.wdeg)
summary(gr2016.f.GDP)

#plot final filtered networks
plot(gr2015.f.GDP)
plot(gr2016.f.GDP)

write_graph(gr2015.f.GDP, "entire2015.graphml", format = "graphml")
write_graph(gr2016.f.GDP, "entire2016.graphml", format = "graphml")
#remove singleton in 2016 netowrk
deg <- V(gr2016.f.GDP)[degree(gr2016.f.GDP)==0]
gr2016.f.GDP <- delete_vertices(gr2016.f.GDP, deg)
plot(gr2016.f.GDP)
summary(gr2016.f.GDP)


#log-log comparative

wtbin2015 <- tabulate(E(gr2015.f.GDP)$weight)
wtbin2015


bindf1 <- data.frame(bin=seq(1, length(wtbin2015)),
                     count=wtbin2015, year=c("2015"))

wtbin2016 <- tabulate(E(gr2016.f.GDP)$weight)
wtbin2016

bindf2 <- data.frame(bin=seq(1, length(wtbin2016)),
                     count=wtbin2016, year=c("2016"))

bindf1 <- bindf1[bindf1$count>0,]
bindf1$freq <- bindf1$count / sum(bindf1$count)
View(bindf1)

bindf2 <- bindf2[bindf2$count>0,]
bindf2$freq <- bindf2$count / sum(bindf2$count)

edgeweights <- rbind(bindf1, bindf2)
View(edgeweights)

p <- ggplot(data = edgeweights, aes(x=bin, y =freq, color=year))
p <- p + geom_point() + geom_smooth(se=FALSE)
P <- p + scale_y_log10("Frequency")
p <- p + scale_x_log10("weight")
print(p)

#degree distribution 2015
summary(degree(gr2015.f.GDP))
deg_2015 <- factor(degree(gr2015.f.GDP))

g <- ggplot(data.frame(Degree2015=deg_2015), aes(x=Degree2015))
g <- g+geom_bar(stat = "count")
print(g)

#degree distribution 2015
deg_2016 <- factor(degree(gr2016.f.GDP))

g <- ggplot(data.frame(Degree2016=deg_2016), aes(x=Degree2016))
g <- g+geom_bar(stat = "count")
print(g)


# Using the degree_distribution function 2015
distrib.df.2015 <- data.frame(Degree2015=seq(0,159),
                              Fraction=degree_distribution(gr2015.f.GDP))
#Show the distribution as a smoothed line: quasi-PDF.

# Create a sensible x axis
g <- ggplot(distrib.df.2015, aes(x=Degree2015, y=Fraction))

g <- g + geom_smooth()
print(g)

#2016
# Using the degree_distribution function 2016
summary(degree(gr2016.f.GDP))
distrib.df.2016 <- data.frame(Degree2016=seq(0,144),
                              Fraction=degree_distribution(gr2016.f.GDP))
#Show the distribution as a smoothed line: quasi-PDF.

# Create a sensible x axis
g <- ggplot(distrib.df.2016, aes(x=Degree2016, y=Fraction))

g <- g + geom_smooth()
print(g)
#compute density related measures
edge_density(gr2015.f.GDP)
edge_density(gr2016.f.GDP)

#compute global transitivity by year
transitivity(gr2015.f.GDP, type = "global") #0.5782171
transitivity(gr2016.f.GDP, type = "global") #0.4729305

Trans2015Global <- transitivity(gr2015.f.GDP, type = "global")
Trans2016Global <- transitivity(gr2016.f.GDP, type = "global")

#compute and plot local transitivity by year
trans2015 <- transitivity(gr2015.f.GDP, type ="local")
trans2016 <- transitivity(gr2016.f.GDP, type = "local")

trans_df <- rbind(data.frame(LocalTrans=trans2015,
                             year="2015"),
                  data.frame(LocalTrans=trans2016,
                             year="2016"))
View(trans_df)


p <- ggplot(trans_df, aes(x=LocalTrans, fill=year))
p<- p +geom_histogram(binwidth = 0.1, position = "dodge")
print(p)


#compute assotativity by GDP
assortativity_nominal(gr2015.f.GDP,
                      types = factor(V(gr2015.f.GDP)$GDP)) #-0.005493275

assortativity_nominal(gr2016.f.GDP,
                      types = factor(V(gr2016.f.GDP)$GDP)) #-0.004477029

#Comunity detection 2015
set.seed(20170527)
gr2015.wt4 <- cluster_walktrap(gr2015.f.GDP, steps = 4)
gr2015.wt5 <- cluster_walktrap(gr2015.f.GDP, steps = 5)
gr2015.wt6 <- cluster_walktrap(gr2015.f.GDP, steps = 6)
gr2015.wt7 <- cluster_walktrap(gr2015.f.GDP, steps = 7)
gr2015.wt8 <- cluster_walktrap(gr2015.f.GDP, steps = 8)
gr2015.wt9 <- cluster_walktrap(gr2015.f.GDP, steps = 9)
gr2015.wt10 <- cluster_walktrap(gr2015.f.GDP, steps = 10)

gr2015.wt4  #0.3492951(mod) 24
gr2015.wt5  #0.4840968(mod) 16
gr2015.wt6  #0.4989406(mod) 13
gr2015.wt7  #0.4780565(mod) 11
gr2015.wt8  #0.5109756(mod) 18 look at this one first!
gr2015.wt9  #0.5308325(mod) 16
gr2015.wt10 #0.5711362(mod) 13

gr2015.mod <- lapply(list(gr2015.wt4, gr2015.wt5, gr2015.wt6, gr2015.wt7, gr2015.wt8,
                          gr2015.wt9, gr2015.wt10), modularity)
gr2015.mod
gr2015.len <- lapply(list(gr2015.wt4, gr2015.wt5, gr2015.wt6, gr2015.wt7, gr2015.wt8,
                          gr2015.wt9, gr2015.wt10), length)
gr2015.len

p <- ggplot(data.frame(Modularity= as.vector(gr2015.mod, mode = "numeric"),
                       Len=as.vector(gr2015.len, mode = "numeric"),
                       Algorithm = c("wt4", "wt5", "wt6", "wt7", "wt8", "wt9", "wt10")),
            aes(y=Modularity, x=Len, label=Algorithm))
p <- p + geom_text()
print(p)


#find cluster id for Syrian 2015
which(V(gr2015.f.GDP)$label=="Syrian.Arab.Rep.") #117
which(V(gr2015.f.GDP)$label=="Afghanistan") #1
which(V(gr2015.f.GDP)$label=="Turkey") #122
membership(gr2015.wt5) [117] #cluster id 156-5
membership(gr2015.wt6) [117] #cluster id 156-3
membership(gr2015.wt7) [117] #cluster id 156-2
membership(gr2015.wt8) [117] #cluster id 156-6
membership(gr2015.wt9) [117] #cluster id 156-6
membership(gr2015.wt10) [117] #cluster id 156-8

V(gr2015.f.GDP)[membership(gr2015.wt9)==6] #16 countires
V(gr2015.f.GDP)[membership(gr2015.wt8)==6] #16 countires
V(gr2015.f.GDP)[membership(gr2015.wt10)==8] #15 countires
V(gr2015.f.GDP)[membership(gr2015.wt7)==2] #79 countires
V(gr2015.f.GDP)[membership(gr2015.wt6)==3] #64 countires
V(gr2015.f.GDP)[membership(gr2015.wt5)==5] #63 countires(with 1 and 117)!

#Comunity detection 2016
set.seed(20170527)
gr2016.wt4 <- cluster_walktrap(gr2016.f.GDP, steps = 4)
gr2016.wt5 <- cluster_walktrap(gr2016.f.GDP, steps = 5)
gr2016.wt6 <- cluster_walktrap(gr2016.f.GDP, steps = 6)
gr2016.wt7 <- cluster_walktrap(gr2016.f.GDP, steps = 7)
gr2016.wt8 <- cluster_walktrap(gr2016.f.GDP, steps = 8)
gr2016.wt9 <- cluster_walktrap(gr2016.f.GDP, steps = 9)
gr2016.wt10 <- cluster_walktrap(gr2016.f.GDP, steps = 10)

gr2016.wt4  #0.4546211(mod) 12
gr2016.wt5  #0.444855(mod) 11
gr2016.wt6  #0.4649787(mod) 10
gr2016.wt7  #0.4949791(mod) 7
gr2016.wt8  #0.4975094(mod) 8
gr2016.wt9  #0.5095541(mod) 7
gr2016.wt10 #0.4942319(mod) 10

gr2016.mod <- lapply(list(gr2016.wt4, gr2016.wt5, gr2016.wt6, gr2016.wt7, gr2016.wt8,
                          gr2016.wt9, gr2016.wt10), modularity)
gr2016.mod
gr2016.len <- lapply(list(gr2016.wt4, gr2016.wt5, gr2016.wt6, gr2016.wt7, gr2016.wt8,
                          gr2016.wt9, gr2016.wt10), length)
gr2016.len

p <- ggplot(data.frame(Modularity= as.vector(gr2016.mod, mode = "numeric"),
                       Len=as.vector(gr2016.len, mode = "numeric"),
                       Algorithm = c("wt4", "wt5", "wt6", "wt7", "wt8", "wt9", "wt10")),
            aes(y=Modularity, x=Len, label=Algorithm))
p <- p + geom_text()
print(p)

#find cluster id for Syrian 2016
which(V(gr2016.f.GDP)$label=="Syrian.Arab.Rep.") #124
which(V(gr2016.f.GDP)$label=="Afghanistan") #1
which(V(gr2016.f.GDP)$label=="Turkey") #130
membership(gr2016.wt5) [124] #cluster id 169-4
membership(gr2016.wt6) [124] #cluster id 156-5
membership(gr2016.wt7) [124] #cluster id 156-3
membership(gr2016.wt8) [124] #cluster id 156-3
membership(gr2016.wt9) [124] #cluster id 156-5
membership(gr2016.wt10) [124] #cluster id 156-4
membership(gr2016.wt10) [130]

V(gr2016.f.GDP)[membership(gr2016.wt5)==4] #56 countires(w/1)
V(gr2016.f.GDP)[membership(gr2016.wt6)==5] #58 countires(with 1)
V(gr2016.f.GDP)[membership(gr2016.wt7)==3] #46 countires(with 1)
V(gr2016.f.GDP)[membership(gr2016.wt8)==3] #46 countires(with 1)
V(gr2016.f.GDP)[membership(gr2016.wt9)==5] #42 countires(with 1)
V(gr2016.f.GDP)[membership(gr2016.wt10)==4] #42 countires(with 1)


cluster2015 <-induced.subgraph(gr2015, c(1, 2, 3, 7, 8, 10, 11, 14, 17, 22, 23, 25, 29, 30, 32, 36, 43, 46,
                                         47, 50, 54, 61, 62, 64, 65, 67, 69, 70, 73, 74, 75, 77, 79, 80, 84, 
                                         86, 87, 90, 92, 104, 107, 108, 109, 112, 114, 117, 122, 125, 126, 128, 
                                         131, 132, 133, 139, 149, 154, 156, 157, 161, 163, 164, 172, 174))

write_graph(cluster2015, "cluster2015.graphml", format = "graphml")
summary(cluster2015)


cluster2016 <- induced.subgraph(gr2016, c(1, 2, 7, 8, 9, 10, 16, 21, 25, 30, 43, 46, 50, 55, 62, 63, 65, 66,
                                          67, 69, 72, 73, 77, 80, 81, 82, 83, 85, 87, 88, 91, 93, 94, 95, 98, 
                                          101, 107, 113, 114, 120, 126, 128, 129, 131, 138, 149, 150, 155, 162,
                                          167, 168, 169, 170, 172, 178, 179, 180, 185))
write_graph(cluster2016, "cluster2016.graphml", format = "graphml")
summary(cluster2016)

#compute centrality measure 2015 community 
cluster2015.id <- degree(cluster2015, normalized = TRUE, mode = "in")
cluster2015.od <- degree(cluster2015, normalized = TRUE, mode = "out")
cluster2015.wd <- graph.strength(cluster2015) / max(graph.strength(cluster2015))
cluster2015.bet <- betweenness(cluster2015, normalized = TRUE)
cluster2015.cloi <- closeness(cluster2015, normalized = TRUE, mode = "in")
cluster2015.eig <- eigen_centrality(cluster2015, directed = TRUE)
cluster2015.pr <- page_rank(cluster2015)

cluster2015.cent <- data.frame(InDegree=cluster2015.id, 
                               OutDegree=cluster2015.od,
                               wtdegree = cluster2015.wd,
                               betw=cluster2015.bet,
                               cloi=cluster2015.cloi,
                               pr=cluster2015.pr$vector,
                               eig = cluster2015.eig$vector)
ggcorr(cluster2015.cent, label = TRUE)                    

#compute centrality measure 2016 community 
cluster2016.id <- degree(cluster2016, normalized = TRUE, mode = "in")
cluster2016.od <- degree(cluster2016, normalized = TRUE, mode = "out")
cluster2016.wd <- graph.strength(cluster2016) / max(graph.strength(cluster2016))
cluster2016.bet <- betweenness(cluster2016, normalized = TRUE)
cluster2016.cloi <- closeness(cluster2016, normalized = TRUE)
cluster2016.eig <- eigen_centrality(cluster2016, directed = TRUE)
cluster2016.pr <- page_rank(cluster2016)

cluster2016.cent <- data.frame(InDegree=cluster2016.id,
                               OutDegree=cluster2016.od,
                               wtdegree = cluster2016.wd,
                               betw=cluster2016.bet,
                               cloi=cluster2016.cloi,
                               pr=cluster2016.pr$vector,
                               eig = cluster2016.eig$vector)
ggcorr(cluster2016.cent, label = TRUE)  

#Comunity detection 2016
set.seed(20170527)
gr2016.wt4 <- cluster_walktrap(gr2016.f, steps = 4)
gr2016.wt5 <- cluster_walktrap(gr2016.f, steps = 5)
gr2016.wt6 <- cluster_walktrap(gr2016.f, steps = 6)
gr2016.wt7 <- cluster_walktrap(gr2016.f, steps = 7)
gr2016.wt8 <- cluster_walktrap(gr2016.f, steps = 8)
gr2016.wt9 <- cluster_walktrap(gr2016.f, steps = 9)
gr2016.wt10 <- cluster_walktrap(gr2016.f, steps = 10)

gr2016.wt4
gr2016.wt5
gr2016.wt6
gr2016.wt7
gr2016.wt8
gr2016.wt9
gr2016.wt10

gr2016.mod <- lapply(list(gr2016.wt4, gr2016.wt5, gr2016.wt6, gr2016.wt7, gr2016.wt8,
                          gr2016.wt9, gr2016.wt10), modularity)

gr2016.mod


gr2016.len <- lapply(list(gr2016.wt4, gr2016.wt5, gr2016.wt6, gr2016.wt7, gr2016.wt8,
                          gr2016.wt9, gr2016.wt10), length)

gr2016.len


p <- ggplot(data.frame(Modularity= as.vector(gr2016.mod, mode = "numeric"),
                       Len=as.vector(gr2016.len, mode = "numeric"),
                       Algorithm = c("wt4", "wt5", "wt6", "wt7", "wt8", "wt9", "wt10")),
            aes(y=Modularity, x=Len, label=Algorithm))
p <- p + geom_text()
print(p)

#use wt9: even if it's length is not the smallest, it has the highest modularity
#WT9 <- walktrap.community(gr2016.f)
#WT9
#modularity((WT9)) #0.5199151

#find cluster id for Syrian 2016 in wt9

which(V(gr2016.f)$label=="Syrian.Arab.Rep.")#126

membership(gr2016.wt9) #cluster id 126 in cluster 3 (total 196 countries)

V(gr2016.f)[membership(gr2016.wt9)==3]
cluster2016 <- induced.subgraph(gr2016.f, c(1,2, 7:10, 16, 21, 25, 30, 43, 46, 47, 50, 55, 62, 63, 65:69, 72, 73, 77, 
                                            80:83, 85, 87, 88, 92, 93, 95, 98, 101, 107, 112:114, 120, 126, 128, 129, 
                                            131, 138, 139, 140, 150, 152, 167:169, 172, 175, 178, 179, 185, 196))
cluster2016

#find cluster id for Syrian 2016 in wt7
membership(gr2016.wt7) #cluster id 126 in cluster 4 (total 196 countries)
V(gr2016.f)[membership(gr2016.wt7)==4]


write_graph(cluster2016, "2016clusterNew.graphml", format = "graphml")
write_graph(cluster2016, "Syrian.Arab.Rep.cluster2016.graphml", format = "graphml")


#compute centrality measure 2016 cluster 
cluster2016.dc <- degree(cluster2016, normalized = TRUE, mode = "in")
cluster2016.wd <- graph.strength(cluster2016) / max(graph.strength(cluster2016))
cluster2016.bet <- betweenness(cluster2016, normalized = TRUE)
cluster2016.cloi <- closeness(cluster2016, normalized = TRUE, mode = "in")
cluster2016.cloo <- closeness(cluster2016, normalized = TRUE, mode = "out")
cluster2016.eig <- eigen_centrality(cluster2016, directed = TRUE)
cluster2016.pr <- page_rank(cluster2016)

cluster2016.cent <- data.frame(degree=cluster2016.dc, 
                               wtdegree = cluster2016.wd,
                               betw=cluster2016.bet,
                               cloi=cluster2016.cloi,
                               cloo=cluster2016.cloo,
                               pr=cluster2016.pr$vector,
                               eig = cluster2016.eig$vector)
ggcorr(cluster2016.cent, label = TRUE)  

#compare inDegree distribution
idist.df <- rbind(
  data.frame(InDegree=degree(cluster2015, mode = "in"), Network = "2015"),
  data.frame(InDegree=degree(cluster2016, mode = "in"), Network = "2016")
)
View(idist.df)
g <- ggplot(data = idist.df, aes(x=InDegree, fill=Network))
g <- g + geom_histogram(position = "dodge", binwidth = 1.0)
print(g)


#visualize GDP didn't work 
#gdp.data <- as.data.frame(V(gr2015.f.GDP)$GDP)
#g <- ggplot(data = gdp.data, aes(x=V(gr2015.f.GDP)$GDP))
#g <- g + geom_histogram(position = "dodge", binwidth = 1.0)
#print(g)

#compara outDegree distribution
odist.df <- rbind(
  data.frame(OutDegree=degree(cluster2015, mode = "out"), Network = "2015"),
  data.frame(OutDegree=degree(cluster2016, mode = "out"), Network = "2016")
)
View(odist.df)
g <- ggplot(data = odist.df, aes(x=OutDegree, fill=Network))
g <- g + geom_histogram(position = "dodge", binwidth = 1.0)
print(g)

library(GGally)

#Comparative GDP and In-Degree

gdp <- log(V(cluster2015)$GDP)

gdp2015 <- as.data.frame(gdp)
View(gdp)

p <- ggplot(data = cluster2015, aes(x=cluster2015.id, y = gdp))
p <- p + geom_point() + geom_smooth(se=FALSE)
print(p)

logGDP2015<- log(V(cluster2015)$GDP)
logGDP

binGDP2015 <- data.frame(GDP = logGDP2015, InDegree=cluster2015.id, year=c("2015"))
binGDP2015

logGDP2016 <- log(V(cluster2016)$GDP)
binGDP2016 <- data.frame(GDP = logGDP2016, InDegree=cluster2016.id, year=c("2016"))

GDP.df <- rbind(binGDP2015, binGDP2016)
View(GDP.df)

p <- ggplot(data = GDP.df, aes(x=InDegree, y = GDP, color=year))
p <- geom_smooth(na.rm = FALSE)
print(p)

#plot GDP by in-degree for year 2015
logGDP2015 <- log(V(cluster2015)$GDP)
GDPbin2015 <- as.data.frame(logGDP2015)
GDPbin2015
bin2015 <- data.frame(GDP=GDPbin2015, InDegree=cluster2015.id)
#View(bin2015)
#bin2015
p <- ggplot(data = bin2015, aes(x=logGDP2015, y = InDegree))
p <- p +geom_point() +geom_smooth(se=FALSE)
print(p)

#plot GDP by in-degree for year 2016
logGDP2016 <- log(V(cluster2016)$GDP)
GDPbin2016 <- as.data.frame(logGDP2016)
#GDPbin2015
bin2016 <- data.frame(GDP=GDPbin2016, InDegree=cluster2016.id)
#View(bin2016)
#bin2016
p <- ggplot(data = bin2016, aes(x=logGDP2016, y = InDegree))
p <- p +geom_point() +geom_smooth(se=FALSE)
print(p)

#CUG and QAP test

source("mycugtest.R")
source("myqaptest.R")

set.seed(20170509)
#CUG 2015
#turn attribute into factor
attrGDP2015 <- as.factor(V(cluster2015)$GDP)
assortativity_nominal(cluster2015, types = attrGDP2015, directed = TRUE)
## -0.01203819
GDP2015.cug <- mycugtest(cluster2015, assortativity_nominal, cmode = "edges", types= attrGDP2015)
print.cug.test(GDP2015.cug)

plot.cug.test(GDP2015.cug)

#QAP 2015
GDP2015.qap <- myqaptest(cluster2015, assortativity_nominal, 
                      types=attrGDP2015)
summary.qaptest(GDP2015.qap)

plot.qaptest(GDP2015.qap)

#CUG 2016
attrGDP2016 <- as.factor(V(cluster2016)$GDP)
assortativity_nominal(cluster2016, types = attrGDP2016, directed = TRUE)
## -0.01203819
GDP2016.cug <- mycugtest(cluster2016, assortativity_nominal, cmode = "edges", types= attrGDP2016)
print.cug.test(GDP2016.cug)

plot.cug.test(GDP2016.cug)

#QAP 2016
GDP2016.qap <- myqaptest(cluster2016, assortativity_nominal, 
                         types=attrGDP2016)
summary.qaptest(GDP2016.qap)

plot.qaptest(GDP2016.qap)

#CUG entire network

attrGDPEntire <- as.factor(V(gr2015.f.GDP)$GDP)
assortativity_nominal(gr2015.f.GDP, types = attrGDPEntire, directed = TRUE)
## -0.01203819
GDPentire.cug <- mycugtest(gr2015.f.GDP, assortativity_nominal, cmode = "edges", types= attrGDPEntire)
print.cug.test(GDPentire.cug)

plot.cug.test(GDPentire.cug)

#QAP Entire
GDPentire.qap <- myqaptest(gr2015.f.GDP, assortativity_nominal, 
                         types=attrGDPEntire)
summary.qaptest(GDPentire.qap)

plot.qaptest(GDPentire.qap)


