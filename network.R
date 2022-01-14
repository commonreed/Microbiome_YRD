####生成网络图####
#############基于丰度相关性的微生物共发生网络
##计算微生物丰度间的相关系数
library(Hmisc)
library(eoffice)
library(patchwork)
otu_bact_1 <- read.delim('data/bacteria_f_even_depth_w_tax.xls', row.name = 1, check.names = FALSE)
otu_bact_2 <- read.delim('data/fungi_f_even_depth_w_tax.xls', row.name = 1, check.names = FALSE)
row.names(otu_bact_2) <- paste0("fungi_", row.names(otu_bact_2))
otu_bact <- rbind(otu_bact_1, otu_bact_2)
otu_bact <- otu_bact[-21]


#可选事先过滤一些低丰度或低频的类群
delOTU <- function(sp){
  sum(sp != 0) >= 5 &  sum(sp) > 100
}
otu_b_del_0 <- otu_bact[apply(otu_bact, 1, delOTU), ]
otu_b_del <- t(otu_b_del_0)

dim(otu_bact)
dim(otu_b_del_0)

names(otu_b_del) <- row.names(otu_b_del_0)
row.names(otu_b_del) <- names(otu_b_del_0)

ab <- apply(otu_b_del, 2, sum)
write.csv(ab, "output/aboundance.csv")

#计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
genus_corr <- rcorr(otu_b_del, type = 'spearman')

#阈值筛选
#将 spearman 相关系数低于 0.09 的关系剔除，即 r>=0.9
r <- genus_corr$r
r[abs(r) < 0.9] <- 0

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
p <- genus_corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0

#根据上述筛选的 r 值和 p 值保留数据
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]

#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
write.table(data.frame(z, check.names = FALSE), 'output/genus_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

##获得网络
library(igraph)

#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数 
g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
g

#自相关也可以通过该式去除
g <- simplify(g)

#孤立节点的删除（删除度为 0 的节点）
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))

#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)

#为节点（微生物属）添加属性信息（界门纲目科属水平注释）
#“genus_taxonomy.txt” 记录了微生物的属性，读入该表后根据已知网络节点匹配对应的行
# tax <- read.delim('genus_taxonomy.txt', row.name = 1, check.names = FALSE, stringsAsFactors = FALSE)
# tax <- tax[as.character(V(g)$name), ]
# 
# V(g)$kingdom <- tax$kingdom
# V(g)$phylum <- tax$phylum
# V(g)$class <- tax$class
# V(g)$order <- tax$order
# V(g)$family <- tax$family
# V(g)$genus <- tax$genus

#查看网络图
g
plot(g)

##网络文件输出，输出特定的网络文件类型，便于后续数据分析需求
#邻接矩阵，出了上述提到的在计算相关系数后，输出筛选后的相关系数矩阵外
#还可以由 igraph 的邻接列表转换
adj_matrix <- as.matrix(get.adjacency(g, attr = 'correlation'))
write.table(data.frame(adj_matrix, check.names = FALSE), 'output/network.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

#边列表
edge <- data.frame(as_edgelist(g))    #igraph 的邻接列表转为边列表

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$correlation
)
head(edge_list)

write.table(edge_list, 'output/network.edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#节点属性列表
node_list <- data.frame(
  label = names(V(g))
  # kingdom = V(g)$kingdom,
  # phylum = V(g)$phylum,
  # class = V(g)$class,
  # order = V(g)$order,
  # family = V(g)$family,
  # genus = V(g)$genus
)
head(node_list)

write.table(node_list, 'output/network.node_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#边列表节点属性列表可以导入至 gephi 或 cytoscape 等网络可视化软件中进行编辑
#此外 igraph 也提供了可以被 gephi 或 cytoscape 等直接识别的格式
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(g, 'output/network.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(g, 'output/network.gml', format = 'gml')


####计算R2和P####
OTUs <- read.delim('output/network.adj_matrix.txt', row.name = 1, check.names = FALSE)
otus1 <- OTUs
otus1 [otus1 != 0] <- 1
adjacency_unweight <-otus1 
write.table(adjacency_unweight, 'output/adjacency_unweight.txt', sep = '\t', quote = FALSE)

igraph <- graph_from_adjacency_matrix(as.matrix(otus1), mode = 'undirected', weighted = NULL, diag = FALSE)
V(igraph)$degree <- degree(igraph)
degree_dist=table(V(igraph)$degree)
degree_num = as.numeric(names(degree_dist))
degree_count = as.numeric(degree_dist)
dat = data.frame(degree=degree_num,count = degree_count)
hist(V(igraph)$degree)
head(dat)

mod <- nls(count ~ a*degree^b, 
           data = dat, 
           start = list(a = 2, b = 1.5))
summary(mod)
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2

p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value


####模块分析####
soil <- read.table("data/soil.txt", header = TRUE)
f <- function(g){
  res <- g/sum(g)
  return(res)
}
otu_b_del <- apply(otu_b_del, 2, f)
net <- read.csv("network/node.csv", 
         header = TRUE, 
         row.names = 1,
         sep = ",")
otu_b_1 <- otu_b_del[soil$EC < 1,]
otu_b_2 <- otu_b_del[soil$EC > 1,]

fM <- function(n){
  m <- net[net$modularity_class == n, ]$v_name
  res1 <- data.frame(
    module = n,
    salinity = "low",
    abundance = apply(otu_b_1[,m], 1, mean)
  )
  res1$mean <- mean(res1$abundance)
  res1$se <- sd(res1$abundance)/sqrt(11)
  res2 <- data.frame(
    module = n,
    salinity = "high",
    abundance = apply(otu_b_2[,m], 1, mean)
  )
  res2$mean <- mean(res2$abundance)
  res2$se <- sd(res2$abundance)/sqrt(11)
  print(t.test(res1$abundance, res2$abundance)$p.value)
  return(rbind(res1, res2))
}
res <- rbind(
  fM(12),
  fM(24),
  fM(28),
  fM(2),
  fM(27),
  fM(23)
  )

#13.3 11.3 10.2 9.6 8.9 8.5
p <- ggplot(res, aes(x=factor(module,
                levels = c(12,24,28,2,27,23)),
                y=mean,
                fill=salinity)) +
  geom_bar(stat="identity", 
           position=position_dodge(),
           color="black", width=.6) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.15))+
  #scale_fill_manual(values = mycol) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean +se),position=position_dodge(.6), width=.2) +
  theme_bw()
p
#eoffice::topptx(p, "figure/nework.pptx",
#                width = 6, height = 6)

####
fMM <- function(n){
  m <- net[net$modularity_class == n, ]$v_name
  res <- data.frame(
    Module = paste0("Module_", n),
    EC = soil$EC,
    Abundance = apply(otu_b_del[,m], 1, mean)
  )
  return(res)
}
res <- rbind(
  fMM(12),
  fMM(24),
  fMM(28),
  fMM(2),
  fMM(27),
  fMM(23)
)

p1 <- ggplot(res, aes(x = EC,
                y = Abundance,
                color = Module))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  labs(tag = "(B)")+
  ylab("Relative abundance")+
  xlab("Electrical conductivity (mS/cm)")+
  theme_bw()
p1

#eoffice::topptx(p, "figure/nework_lm.pptx",
#                width = 6, height = 6)

lm_res <- data.frame(matrix(nrow=0,ncol = 5))
for(i in unique(res$Module)){
  r <- summary(lm(Abundance ~ EC,data = res[res$Module == i,]))
  lm_res <- rbind(lm_res, c(i, r$r.squared, r$adj.r.squared, r$coefficients[2,1], r$coefficients[2,4]))
  
}
names(lm_res) <- c("Module", "R2", "R2Adj", "slop", "P")
lm_res
write.csv(lm_res, "output/lm_res_network.csv")

m1 <- net[net$modularity_class == 12, ]$v_name
m2 <- net[net$modularity_class == 24, ]$v_name
m3 <- net[net$modularity_class == 28, ]$v_name
m4 <- net[net$modularity_class == 2, ]$v_name
m5 <- net[net$modularity_class == 27, ]$v_name
m6 <- net[net$modularity_class == 23, ]$v_name

o1 <- otu_bact_1[m1,]$taxonomy
findPhyl <- function(n){
  m <- net[net$modularity_class == n, ]$v_name
  o <- otu_bact_1[m,]$taxonomy
  v <- strsplit(o, split = "; ")
  phyl <- c()
  for(i in 1:length(v)){
    phyl <- append(phyl, v[[i]][2])
  } 
  return(phyl)
}
  
table(findPhyl(12))
table(findPhyl(24))
table(findPhyl(28))#up
table(findPhyl(2))
table(findPhyl(27))#up
table(findPhyl(23))


phyl_b <- read.table("data/bact_phyl.xls", row.names = 1, header = TRUE)
phyl_f <- read.table("data/fungi_phyl.xls", row.names = 1, header = TRUE)
#Proteobacteria
#Actinobacteria
#Chloroflexi
#Ascomycota
#Mortierellomycota
#Basidiomycota
phyl <- as.data.frame(t(rbind(phyl_b[1:3,], phyl_f[1:3,])))
phyl$EC <- soil$EC
phyl_short <- reshape2::melt(phyl,
                       id.vars = c("EC"),
                       variable.name = "Phylum",
                       value.name = 'Abundance')
p2 <- ggplot(data = phyl_short, 
       aes(x = EC, y = Abundance, color = Phylum))+
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)+
  labs(tag = "(A)")+
  xlab("Electrical conductivity (mS/cm)")+
  ylab("Relative abundance")+
  theme_bw()
p2
#eoffice::topptx(p, "figure/phyl_lm.pptx",
#                width = 6, height = 6)
p2+p1
topptx(p2+p1, "figure/v1_network_lm.pptx", width = 10, height = 4 )

lm_res <- data.frame(matrix(nrow=0,ncol = 5))
for(i in unique(phyl_short$Phylum)){
  r <- summary(lm(Abundance ~ EC,
                  data = phyl_short[phyl_short$Phylum == i,]))
  lm_res <- rbind(lm_res, c(i, r$r.squared, r$adj.r.squared, r$coefficients[2,1], r$coefficients[2,4]))
  
}
names(lm_res) <- c("Phylum", "R2", "R2Adj", "slop", "P")
lm_res
write.csv(lm_res, "output/lm_res_phylum.csv")

phyl_1 <- phyl[phyl$EC<1,]
phyl_2 <- phyl[phyl$EC>1,]
phyl_1 <- phyl_1[,-7]
phyl_2 <- phyl_2[,-7]
mean_1 <- apply(phyl_1, 2, mean)
se_1 <- apply(phyl_1, 2, sd)/sqrt(11)
mean_2 <- apply(phyl_2, 2, mean)
se_2 <- apply(phyl_2, 2, sd)/sqrt(9)
res_1 <- data.frame(
  Salinity = "Low",
  Phylum = names(phyl_1),
  Mean = mean_1,
  SE = se_1
)
res_2 <- data.frame(
  Salinity = "High",
  Phylum = names(phyl_2),
  Mean = mean_2,
  SE = se_2
)
res <- rbind(res_1, res_2)

p <- ggplot(res, aes(x=factor(Phylum,
                                 levels = c("Proteobacteria","Actinobacteria","Chloroflexi","Ascomycota","Mortierellomycota","Basidiomycota")),
                     y=Mean,
                     fill=Salinity)) +
  geom_bar(stat="identity", 
           position=position_dodge(),
           color="black", width=.6) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.8))+
  #scale_fill_manual(values = mycol) +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean +SE),position=position_dodge(.6), width=.2) +
  theme_bw()
p
eoffice::topptx(p, "figure/phyl.pptx",
                width = 6, height = 6)

for(i in names(phyl_1)){
  print(t.test(phyl_1[i], phyl_2[i])$p.value)
}
