library(tidyverse)
library(vegan)
library(pbapply)
library(plotrix)
library(ggplot2)

rm(list = ls())
setwd("C:/Users/Lenovo/Desktop/鄱阳湖干湿数据/数据分析/Fig.2/github数据上传")
bac <- read.csv("Water_ASV.csv",row.names = 1)

colSums(bac)
bac = as.data.frame(t(rrarefy(t(bac), min(colSums(bac)))))

colSums(bac)
write.csv(bac, file ="Water_ASV_rarefaction.csv")


bac <- read.csv("Water_ASV_rarefaction.csv",row.names = 1)
sam <- read.csv("Water_sample_change.csv",row.names = 1)

asv_t <- t(bac) %>% as.data.frame()
est <- estimateR(asv_t)
Richness <- est[1, ]
result <- data.frame(Richness)
sam_alp <- cbind(sam, result)

sam_Water <- filter(sam_alp,group %in% "Water") 
alp_W <- filter(sam_Water,factor %in% "W")
alp_D <- filter(sam_Water,factor %in% "D")

ct <- alp_D$Richness-alp_W$Richness

Water_ric <- data.frame(ct)
Water_ric$names <- rownames(alp_W)

rownames(Water_ric) <- Water_ric$names

Water_richness <- cbind(Water_ric, alp_W$Richness)
colnames(Water_richness)[3] <- "Richness"

ct_rs <- 1-(2*abs(Water_richness$ct))/(abs(Water_richness$ct) + Water_richness$Richness) 

Water_richness_resistance <- cbind(Water_richness, ct_rs)
write.csv(Water_richness_resistance,file = "Water_richness_resistance.csv")




#community_resistance
biodiv <- rbind(bac)
asv <- as.data.frame(t(biodiv))
sam_asv <- cbind(sam, asv) %>% as.data.frame()
ct <- filter(sam_asv,factor %in% c("W","D"))

asv_ct <- select(ct,c(-1:-2))

calc_bray <- function(data, row_pairs) {
  lapply(1:nrow(row_pairs), function(i) {
    idx <- row_pairs[i, ]
    as.vector(vegdist(data[idx[1]:idx[2], ], method = "bray"))
  })
}

row_pairs <- matrix(seq_len(nrow(asv_ct)), ncol = 2, byrow = TRUE)

datasets <- list(ct = asv_ct)
bray_list <- pblapply(datasets, function(x) {
  calc_bray(x, row_pairs)
})

df_bray <- as.data.frame(lapply(bray_list, unlist))
community_rs <- 1- df_bray
write.csv(community_rs,file = "Water_community_resistance.csv")

