rm(list = ls())
setwd("~/ZINARA/Gnomic/1-StageM2/3-Annotation")

# data_ <- data.frame( Consequences = c(rep("missense_variant", 29842), rep("stop_gained", 1276),
#                                      rep("frameshift_variant", 828), rep("start_lost", 79), 
#                                      rep("inframe_insertion", 249), rep("inframe_deletion", 343),
#                                      rep("synonymous_variant", 15055), rep("Other", 51)) )
# 
# png(file = "./3-Sortie-pieplot/Annotation.png", width = 900, height = 800)
# pie(table(data_$Consequences))
# dev.off()


library(tidyverse)
library(viridis)
library(hrbrthemes)
hrbrthemes::import_roboto_condensed()

data_2 <- data.frame(Consequences=c("missense_variant", "stop_gained", "frameshift_variant",
                                   "start_lost", "inframe_insertion", "inframe_deletion",
                                   "synonymous_variant", "Other"),
                     Effectif=c(29842, 1276, 828, 79, 249, 343, 15055, 51))
data_2["Pourcentage"] <- round(data_2$Effectif/sum(data_2$Effectif), 3)
library(scales)
png(file = "./3-Sortie-pieplot/Annotation.png", width = 900, height = 800)
ggplot(data=data_2, aes(x="", y=Pourcentage, fill=Consequences)) +
  geom_bar(width = 1, stat = "identity", color="white") +
  coord_polar("y", start=0) +
  theme_ipsum() +
  # geom_text(aes(y = Pourcentage/2 + c(0, cumsum(Pourcentage)[-length(Pourcentage)]), 
  #               label = percent(Pourcentage)), size=8) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title=element_text(size=25),
    legend.text=element_text(size=25),
    axis.text.x = element_text(size = 25)
  )
dev.off()