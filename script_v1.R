####Package and data preparation####
##Loading common packages
library(tidyverse)
library(ggplot2)
library(vegan)
library(eoffice)
library(patchwork)
library(ggcorrplot)

# Longitude Latitude
# Keep N P pH EC and remove EC_re
soil <- read_tsv("data/soil.txt") %>%
  select(-EC_re)

# Shoot height, leaf length, panicle length
# Seed weight, seed number, per seed weight
# Haplotype
reed <- read_tsv("data/reed.txt")%>%
  select(-Panicle)

# OTU (ASV)
otu_bact <- read_tsv("data/bacteria_f_even_depth_w_tax.xls")%>%
  select(-c(ID,taxonomy)) %>%
  t()
otu_fungi <- read_tsv("data/fungi_f_even_depth_w_tax.xls")%>%
  select(-c(ID,taxonomy)) %>%
  t()
fSat <- function(sp){
  sum(sp != 0)/length(sp) < 0.5
}
fCore <- function(sp){
  sum(sp != 0)/length(sp) >= 0.75
}
otu_b_sat <- otu_bact[,apply(otu_bact, 2, fSat)]
otu_b_core <- otu_bact[,apply(otu_bact, 2, fCore)]
otu_f_sat <- otu_fungi[,apply(otu_fungi, 2, fSat)]
otu_f_core <- otu_fungi[,apply(otu_fungi, 2, fCore)]

####Map#####
soil$Haplotype <- reed$Haplotype
map <- ggplot(soil) +
  coord_sf(crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
           ylim = c(min(soil$Latitude)-0.001, max(soil$Latitude+0.001)))+
  geom_point(aes(x = Longitude,
                 y = Latitude, 
                 color = EC,
                 shape = Haplotype),
             size = 4) +
  scale_color_gradient(low = "cyan",high = "red")+
  ggspatial::annotation_scale(location = "bl", width_hint = 0.4) +
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true", 
                                    pad_x = unit(0.05, "in"), pad_y = unit(0.05, "in"))+
  theme_minimal()
#topptx(map,"figure/v2_map.pptx", width = 9)



#####Alpha Diversity####
#library(BiocManager)
#BiocManager::install("microbiome")
library(microbiome)
library(knitr)

myDiv <- function(data, taxa, type){
  d <- alpha(t(data), 
             index = c("Chao1",
                       "Shannon",
                       "Observed",
                       "evenness_pielou"))
  d$taxa <- taxa
  d$type <- type
  d$Name <- row.names(d)
  return(d)
}

div_all <-rbind(
  myDiv(otu_bact, "Bacteria", "All"),
  myDiv(otu_b_core, "Bacteria", "Core"),
  myDiv(otu_b_sat, "Bacteria", "Satellite"),
  myDiv(otu_fungi, "Fungi", "All"),
  myDiv(otu_f_core, "Fungi", "Core"),
  myDiv(otu_f_sat, "Fungi", "Satellite")
)

data_all <- div_all %>%
  left_join(soil, by = "Name")%>%
  left_join(reed, by = "Name")%>%
  mutate(Haplotype.y = NULL) %>%
  rename(Haplotype = Haplotype.x)

vvDiv <- function(div){
  p <- ggplot(data_all,
              aes(x = EC, y = get(div), color = type))+
    geom_point(mapping = aes(shape = Haplotype))+
    geom_smooth(method = "lm", se = FALSE)+
    xlab("Electrical conductivity (mS/cm)")+
    ylab(div)+
    guides(color= "none",shape= "none")+
    facet_wrap(~taxa, scales = "free")+
    theme_bw()+
    theme(legend.position = "top")
  return(p)
}
t <- vvDiv("observed")
t
#topptx(t, "figure/v1_diversity_legend.pptx", width = 9, height = 6)
p <- vvDiv("observed")+vvDiv("chao1")+vvDiv("evenness_pielou")+vvDiv("diversity_shannon")
p
#topptx(p, "figure/v1_diversity.pptx", width = 9, height = 6)

## Many soil and plant parameter
data_all$group <- paste0(data_all$taxa, "_", data_all$type)
lm_res <- data.frame(matrix(nrow=0,ncol = 6))
for(i in unique(data_all$group)){
  for(j in names(data_all)[1:4]){
    r <- summary(lm(get(j) ~ EC,data = data_all[data_all$group == i,]))
    lm_res <- rbind(lm_res, c(i, j, r$r.squared, r$adj.r.squared, r$coefficients[2,1], r$coefficients[2,4]))
  }
}
names(lm_res) <- c("Group", "Index", "R2", "R2Adj", "slop", "P")
lm_res
#write_csv(lm_res, "output/v1_lm_res.csv")

data_all$Haplotype_n <- as.numeric(as.factor(div_bact$Haplotype))
div_bact <- data_all[data_all$group ==  "Bacteria_All",]
div_fungi <- data_all[data_all$group ==  "Fungi_All",]

index <- c("Richness","Chao1","Shannon","Evenness")
para <- c("Soil N", "Soil P", "Soil pH", "Soil EC", "Plant shoot", "Plant leaf", "Plant haplotype")

a1 <- Hmisc::rcorr(as.matrix(div_bact[1:4]), as.matrix(div_bact[c(10:15,22)]), type = "spearman")
r1 <- a1$r[names(div_bact[1:4]), names(div_bact[c(10:15,22)])]
p1 <- a1$P[names(div_bact[1:4]), names(div_bact[c(10:15,22)])]
row.names(r1) <- row.names(p1) <- index
colnames(r1) <- colnames(p1) <- para

a1 <- Hmisc::rcorr(as.matrix(div_bact[1:4]), as.matrix(div_bact[c(10:15,22)]), type = "spearman")
r1 <- a1$r[names(div_bact[1:4]), names(div_bact[c(10:15,22)])]
p1 <- a1$P[names(div_bact[1:4]), names(div_bact[c(10:15,22)])]
row.names(r1) <- row.names(p1) <- paste0("b_", index)
colnames(r1) <- colnames(p1) <- para

a2 <- Hmisc::rcorr(as.matrix(div_fungi[1:4]), as.matrix(div_fungi[c(10:15,22)]), type = "spearman")
r2 <- a2$r[names(div_fungi[1:4]), names(div_fungi[c(10:15,22)])]
p2 <- a2$P[names(div_fungi[1:4]), names(div_fungi[c(10:15,22)])]
row.names(r2) <- row.names(p2) <- paste0("f_", index)
colnames(r2) <- colnames(p2) <- para
r <- rbind(r1,r2)
p <- rbind(p1,p2)
pic <- ggcorrplot(r, p.mat = p,
           method = "circle",
           insig = "pch",
           legend.title = "Correlation r"
           )
pic
topptx(pic, "figure/v2_diversity_coor.pptx", width = 9, height = 6)

###计算谱系信号
library("picante") 
phy <- read.tree("data/wu01.phylo")
calK <- function(p1, p2){
  d1 <- data_all %>% 
    filter(taxa == p1, type == p2) %>% 
    select(c("Name",names(data_all)[1:4]))
  row.names(d1) <- d1$Name
  d1 <- d1[-1]
  res <- multiPhylosignal(d1, phy = phy)
  res$taxa <- p1
  res$type <- p2
  return(res)
}

res_phy <- rbind(
  calK("Bacteria", "All"),
  calK("Bacteria", "Core"),
  calK("Bacteria", "Satellite"),
  calK("Fungi", "All"),
  calK("Fungi", "Core"),
  calK("Fungi", "Satellite")
)
#write.csv(res_phy, "output/v1_phySig.csv")


#BiocManager::install("ggtree")
library(ggtree)
p <- ggtree(phy)+
  geom_tiplab()
#topptx(p, "figure/phy.pptx")


####Beta diversity####
#RDA bact
env <- div_bact[c(10:16)]
colnames(env) <- c("TN", "TP", "pH", "EC", "Haplotype","Shoot", "Leaf")
env$Haplotype <- factor(env$Haplotype)
bact.hell <- decostand(otu_bact, "hellinger")
decorana(bact.hell)
bact.rda <- rda(bact.hell ~ ., env)
#anova.cca(bact.rda,step=1000)
#anova.cca(bact.rda, by='axis',step=1000)
#vif.cca(bact.rda)

RsquareAdj(bact.rda)
#0.42 0.09

res1 <- anova(bact.rda, by='term', permutations=999)
res1$R2 <- res1$Variance/sum(res1$Variance)
res1$R2
#TN TP pH EC 0.06 0.05 0.07 0.12
#R2 0.56 0.54 0.65 0.11
res1 

#RDA fungi
fungi.hell <- decostand(otu_fungi, "hellinger")
decorana(fungi.hell)
fungi.rda <- rda(fungi.hell ~ ., env)
#anova.cca(fungi.rda,step=1000)
#anova.cca(fungi.rda, by='axis',step=1000)
#vif.cca(fungi.rda)
RsquareAdj(fungi.rda)
#0.39 0.04
res2 <- anova(fungi.rda, by='term', permutations=999)
res2$R2 <- res2$Variance/sum(res2$Variance)
res2

res1 <- as.data.frame(res1)
res1$Taxa <- "Bacteria"
res1$Factor <- row.names(res1)
res1

res2 <- as.data.frame(res2)
res2$Taxa <- "Fungi"
res2$Factor <- row.names(res2)
res2
write.csv(rbind(res1,res2), "output/v2_rda.csv")

jpeg("figure/v2_rda.jpg", width = 3200, height = 1600, res = 300)
par(mfrow = c(1,2))
plot(bact.rda,
     scaling = 1,
     display = c("lc", "cn"),
     xlab = "RDA1 (10.6%)",
     ylab = "RDA2 (4.1%)"
     
)
plot(fungi.rda,
     scaling = 1,
     display = c("lc", "cn"),
     xlab = "RDA1 (9.1%)",
     ylab = "RDA2 (5.1%)",
)
dev.off()

myRDA <- function (d){
  hell <- decostand(d, "hellinger")
  rda <- rda(hell ~ ., env)
  plot(rda,
       scaling = 1,
       display = c("lc", "cn"),
       xlab = "RDA1",
       ylab = "RDA2"
  )
  anova(rda, by='term', permutations=999)
}
sink("output/v2_rda.log")
par(mfrow = c(1,1))
myRDA(otu_b_sat) #EC pH
myRDA(otu_b_core) #EC
myRDA(otu_f_sat) #EC shoot
myRDA(otu_f_core)#EC pH
sink()

# correlation between microbial and geographics
# calculation of geograohic distance
library(geosphere)
geo_dist <- as.dist(distm(soil[c('Longitude', 'Latitude')]))
sal_dist <- dist(soil$EC)

mantel(vegdist(otu_bact), geo_dist, method = "spearman")
mantel(vegdist(otu_b_core), geo_dist, method = "spearman")
mantel(vegdist(otu_b_sat), geo_dist, method = "spearman")

mantel(vegdist(otu_fungi), geo_dist, method = "spearman")
mantel(vegdist(otu_f_core), geo_dist, method = "spearman")
mantel(vegdist(otu_f_sat), geo_dist, method = "spearman")

mantel(vegdist(otu_bact), sal_dist, method = "spearman")
mantel(vegdist(otu_b_core), sal_dist, method = "spearman")
mantel(vegdist(otu_b_sat), sal_dist, method = "spearman")

summary(lm(as.vector(vegdist(otu_bact)) ~ as.vector(sal_dist)))
summary(lm(as.vector(vegdist(otu_b_core)) ~ as.vector(sal_dist)))
summary(lm(as.vector(vegdist(otu_b_sat)) ~ as.vector(sal_dist)))


mantel(vegdist(otu_fungi), sal_dist, method = "spearman")
mantel(vegdist(otu_f_core), sal_dist, method = "spearman")
mantel(vegdist(otu_f_sat), sal_dist, method = "spearman")

summary(lm(as.vector(vegdist(otu_fungi)) ~ as.vector(sal_dist)))
summary(lm(as.vector(vegdist(otu_f_core)) ~ as.vector(sal_dist)))
summary(lm(as.vector(vegdist(otu_f_sat)) ~ as.vector(sal_dist)))



d_dist <- data.frame(
  type = rep(c("Bacteria", "Fungi"), each = 190*3),
  group = rep(c("All", "Core", "Satellite"), 
              each = 190, time = 2),
  geo_dist = as.vector(geo_dist),
  sal_dist = as.vector(sal_dist),
  mic_dist = c(as.vector(vegdist(otu_bact)),
               as.vector(vegdist(otu_b_core)),
               as.vector(vegdist(otu_b_sat)),
               as.vector(vegdist(otu_fungi)),
               as.vector(vegdist(otu_f_core)),
               as.vector(vegdist(otu_f_sat))
  )
)

p1 <- ggplot(d_dist, aes(x = geo_dist, 
                         y = mic_dist, 
                         color = group))+
  geom_point()+
  facet_wrap(~type)+
  theme_bw()+
  xlab("Geographical distance (m)")+
  ylab("Bray_Curtis distance")+
  theme(legend.position = "noe")

p2 <- ggplot(d_dist, aes(x = sal_dist, 
                         y = mic_dist, 
                         color = group))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  xlab("Delta Salinity (dS/cm)")+
  ylab("Bray_Curtis distance")+
  facet_wrap(~type)+
  theme_bw()+
  theme(legend.position = "bottom")

p1/p2
topptx(p1/p2, "figure/v1_dist.pptx", width = 6, height = 6 )

#距离模式，附加遗传
gen_dist <- as.dist(read.csv("data/gen_dist.csv", row.names = 1))
mantel(gen_dist, geo_dist, method = "spearman")
mantel(gen_dist, sal_dist, method = "spearman")
mantel(vegdist(otu_bact), gen_dist, method = "spearman")
mantel(vegdist(otu_b_core), gen_dist, method = "spearman")
mantel(vegdist(otu_b_sat), gen_dist, method = "spearman")
mantel(vegdist(otu_fungi), gen_dist, method = "spearman")
mantel(vegdist(otu_f_core), gen_dist, method = "spearman")
mantel(vegdist(otu_f_sat), gen_dist, method = "spearman")
plot(hclust(gen_dist))

gc()
sol <- metaMDS(otu_bact)
# pcoa <- pcoa(vegdist(otu_bact))
# NMDS <- data.frame(MDS1 = pcoa$vectors[,1], MDS2 = pcoa$vectors[,2])
NMDS <- data.frame(MDS1 = sol$points[,1], MDS2 = sol$points[,2])
NMDS$Haplotype <- reed$Haplotype 
NMDS$Salinity <- factor(soil$EC > 1, labels = c("Low", "High"))
NMDS$Sv <- soil$EC
#使用ggplot2画图，逐步添加

p1 <- ggplot(data = NMDS, aes(MDS1,MDS2, color = Haplotype, shape = Salinity))+
  #geom_point()+
  geom_text(aes(label = round(Sv,1), size = 1))+
  #stat_ellipse(aes(shape = NULL))+
  xlim(-0.4,0.7)+
  ylim(-0.35, 0.35) +
  labs(tag = "(A)",
       caption = paste0("Stress = ",
                        round(sol$stress, digits =  1)))+
  theme_bw()+
  theme(legend.position = "none")

p1

sol <- metaMDS(otu_fungi)
NMDS <- data.frame(MDS1 = sol$points[,1], MDS2 = sol$points[,2])
NMDS$Haplotype <- reed$Haplotype 
NMDS$Salinity <- factor(soil$EC > 1, labels = c("Low", "High"))
NMDS$Sv <- soil$EC
#使用ggplot2画图，逐步添加

p2 <- ggplot(data = NMDS, aes(MDS1,MDS2, color = Haplotype, shape = Salinity))+
  #geom_point()+
  geom_text(aes(label = round(Sv,1)))+
  #stat_ellipse(aes(shape = NULL))+
  xlim(-0.4,0.7)+
  ylim(-0.35, 0.35) +
  labs(tag = "(B)",
       caption = paste0("Stress = ",
                        round(sol$stress, digits =  2)))+
  theme_bw()
p2
topptx(p1+p2, "figure/v2_NMDS.pptx", width = 10, height = 4 )

####Community assembly####
library("spaa")
niche_b <- niche.width(otu_bact, method = "levins")
niche_width <- mean(t(niche_b[otu_bact[20,] != 0]))
niche_width <- c()
for(i in 1:20){
  niche_width <- c(niche_width, mean(t(niche_b[otu_bact[i,] != 0])))
}
niche_w_b <- data.frame(Name = data_all$Name[1:20], niche = niche_width, EC = soil$EC, Haplotype = reed$Haplotype)
summary(lm(niche ~ EC, niche_w_b)) 
ggplot(niche_w_b,aes(EC, niche))+
  geom_point(aes(color = Haplotype))+
  geom_smooth(method = "lm", se = FALSE)+
  theme_bw()
niche_f <- niche.width(otu_fungi, method = "levins")
niche_width <- mean(t(niche_b[otu_fungi[20,] != 0]))
niche_width <- c()
for(i in 1:20){
  niche_width <- c(niche_width, mean(t(niche_b[otu_fungi[i,] != 0])))
}
niche_w_f <- data.frame(Name = data_all$Name[1:20], niche = niche_width, EC = soil$EC, Haplotype = reed$Haplotype)
summary(lm(niche ~ EC, niche_w_f)) 
niche_w_b$Taxa <- "Bacteria"
niche_w_f$Taxa <- "Fungi"
niche_w <- rbind(niche_w_b, niche_w_f)

p <- ggplot(niche_w,aes(EC, niche, color = Taxa))+
  geom_point(aes(shape = Haplotype))+
  geom_smooth(method = "lm", se = FALSE)+
  xlab("Electrical conductivity (mS/cm)")+
  ylab("Habitat Niche breadth (Bcom)")+
  theme_bw()
p
topptx(p, "figure/v2_niche_breadth.pptx", width = 6, height = 4)

#Neutral community model
source("fNeutral.R")
otu_b_1 <- otu_bact[soil$EC < 1,]
otu_b_2 <- otu_bact[soil$EC > 1,]

countN <- function(a){
  sum(a != 0)
}

# otu_f_1 <- otu_f_adj[soil$EC < 1,]
# otu_f_2 <- otu_f_adj[soil$EC > 1,]
otu_f_1 <- otu_fungi[soil$EC < 1,]
otu_f_2 <- otu_fungi[soil$EC > 1,]


pdf("figure/v1_nm.pdf", width = 6, height = 8)
fNeutral(otu_bact)
fNeutral(otu_b_1)
fNeutral(otu_b_2)
fNeutral(otu_fungi)
fNeutral(otu_f_1)
fNeutral(otu_f_2)
dev.off()


?agricolae::kruskal
#library(devtools)
#install_github("GotelliLab/EcoSimR")
library(EcoSimR)

delOTU <- function(sp){
  sum(sp) > 50
}
otu_b0 <- otu_bact[,apply(otu_bact, 2, delOTU)]
otu_b1 <- otu_b_1[,apply(otu_b_1, 2, delOTU)]
otu_b2 <- otu_b_2[,apply(otu_b_2, 2, delOTU)]

n1 <- cooc_null_model(t(otu_b0) != 0,
                      algo="sim1",
                      nReps=1000,
                      burn_in = 500)
n2 <- cooc_null_model(t(otu_b1) != 0,
                      algo="sim1",
                      nReps=1000,
                      burn_in = 500)
n3 <- cooc_null_model(t(otu_b2) != 0,
                      algo="sim1",
                      nReps=1000,
                      burn_in = 500)

n4 <- cooc_null_model(t(otu_fungi) != 0,
                      algo="sim1",
                      nReps=1000,
                      burn_in = 500)
n5 <- cooc_null_model(t(otu_f_1) != 0,
                      algo="sim1",
                      nReps=1000,
                      burn_in = 500)
n6 <- cooc_null_model(t(otu_f_2) != 0,
                      algo="sim1",
                      nReps=1000,
                      burn_in = 500)
sink("output/v1_NullModel.log")
summary(n1)
summary(n2)
summary(n3)
summary(n4)
summary(n5)
summary(n6)
sink()

library(NST)
tda$group
s <- 1:20
names(s) <- soil$Name
s[soil$EC < 1] <- "Low" 
s[soil$EC > 1] <- "High" 
s <- as.data.frame(s)
res_NST <- tNST(otu_bact, s, dist.method = "bray", output.rand = T)
res_NST$index.grp
# group size ST.i.bray NST.i.bray MST.i.bray
# 1   Low   55 0.8015763  0.5125491  0.4734876
# 2  High   36 0.8253855  0.6435163  0.5723696

res_NST_f <- tNST(otu_fungi, s, dist.method = "bray", output.rand = T)
res_NST_f$index.grp
#group size ST.i.bray NST.i.bray MST.i.bray
#1   Low   55 0.6172147  0.2302436  0.1866824
#2  High   36 0.7566362  0.4512252  0.2989462
