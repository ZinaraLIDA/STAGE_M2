rm(list = ls())
setwd("~/ZINARA/Gnomic/1-StageM2/9-Nombre_de_copies/PureCN/Pipeline/8-Analyses")

ids_ <- list.dirs("~/ZINARA/Gnomic/1-StageM2/9-Nombre_de_copies/PureCN/Pipeline/7-sorties_purecn")
ids <- c()
for (id in ids_[2:28]) {
  ids <- c(ids, unlist(strsplit(id, "/"))[11])
}

data_all <- read.csv("~/ZINARA/Gnomic/1-StageM2/9-Nombre_de_copies/PureCN/Pipeline/7-sorties_purecn/D181208/D181208.csv",
                     header = TRUE)
for (id in ids[2:length(ids)]) {
  data_i <- read.csv(sprintf("~/ZINARA/Gnomic/1-StageM2/9-Nombre_de_copies/PureCN/Pipeline/7-sorties_purecn/%s/%s.csv", id, id))
  data_all <- rbind(data_all, data_i)
}

png(file = "~/ZINARA/Gnomic/1-StageM2/9-Nombre_de_copies/PureCN/Pipeline/8-Analyses/Pureté.png",
    width = 900, height = 800)
par(mar=c(5, 6, 5, 1))
hist(data_all$Purity, breaks = 27, freq = FALSE, ylim = c(0,100), xlim = c(0.1,0.4), las=1,
     xlab = "Pureté", ylab = "Pourcentage", main = "", cex.lab=3, cex.axis=2)
dev.off()

png(file = "~/ZINARA/Gnomic/1-StageM2/9-Nombre_de_copies/PureCN/Pipeline/8-Analyses/Ploidie.png",
    width = 900, height = 500)
par(mar=c(5, 6, 5, 1))
hist(data_all$Ploidy, breaks = 27, freq = FALSE, ylim = c(0,100), xlim = c(1.8,2.2), las=1,
     xlab = "Ploïdie", ylab = "Pourcentage", main = "", cex.lab=3, cex.axis=2)
dev.off()