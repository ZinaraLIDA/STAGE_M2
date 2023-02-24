rm(list = ls())
setwd("~/ZINARA/Gnomic/1-StageM2/7-Fraction_cellules_T")
###############
data_ <- read.table( "./TCRA_out_all.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t" )
mean_ <- mean( as.numeric(data_$TCRA.tcell.fraction) )

data_["Batch"] <- unlist( lapply( 1:nrow(data_), function(i) { ifelse( substr(data_$sample[i], 1, 3)=="D21", "batch2", "batch1" ) } ) )
data_["Class"] <- ifelse( data_$TCRA.tcell.fraction < mean_, "cold", "hot" )

subset_b1 <- subset( data_, Batch == "batch1" )
subset_b2 <- subset( data_, Batch == "batch2" )

test_ <- wilcox.test(as.numeric(subset_b1$TCRA.tcell.fraction), as.numeric(subset_b2$TCRA.tcell.fraction), paired = FALSE)
p <- test_$p.value

png(file = "./images_stats/1-Comparaison_batch1_et_batch2.png", width = 900, height = 800)
par(mar=c(5, 6, 5, 1))
boxplot( data_$TCRA.tcell.fraction ~ data_$Batch, xlab = "Batch", ylab = "Fraction de cellules T", 
         main = sprintf( "Wilcoxon-test :\n p-value = %s", round(p, 6) ),
         col=c("brown1", "#00AFBB"), cex.lab=4, cex.axis=3, cex.main=3, ylim=c(0, 0.16) )
legend("right", legend=c("Batch2","Batch1"), col=c("#00AFBB","brown1"),
       pt.cex=4, pch=15, border = "white", cex = 4, bty = "n")
dev.off()

tab_ <- table( data_$Batch, data_$Class )

png(file = "./images_stats/2-Distribution_de_class_selon_batch.png", width = 900, height = 800)
mosaicplot(tab_, main="",xlab="Batch",ylab="Class",col=c("grey85", "grey55"), cex.axis = 1)
legend(x='center', fill=c("grey85", "grey55"), legend=c("cold : TCRA-T-cell-fraction < mean", "hot : TCRA-T-cell-fraction >= mean"),
       title=sprintf("Mean = %s", round(mean_, 5)), cex = 1)
dev.off()
################
data_["Sample_type"] <- unlist( lapply( 1:nrow(data_), function(i) { ifelse( substr(data_$sample[i], 1, 3)=="D18", "cytobrosse", "biopsie" ) } ) )

subset_b1_2 <- subset( data_, Batch == "batch1" )
subset_b1_2_biopsy <- subset( subset_b1_2, subset_b1_2$Sample_type=="biopsie" )
subset_b1_2_biopsy["Patient"] <- c(1:6)
subset_b1_2_cytobrosse<- subset( subset_b1_2, subset_b1_2$Sample_type=="cytobrosse" )
subset_b1_2_cytobrosse["Patient"] <- c(1:6)

test_2 <- wilcox.test(as.numeric(subset_b1_2_biopsy$TCRA.tcell.fraction), as.numeric(subset_b1_2_cytobrosse$TCRA.tcell.fraction), paired = TRUE)
p_2 <- test_2$p.value
png(file = "./images_stats/3-Comparaison_biopsie_et_cytobrosse_de_batch1.png", width = 900, height = 800)
par(mar=c(5, 6, 5, 1))
boxplot( subset_b1_2$TCRA.tcell.fraction ~ subset_b1_2$Sample_type, xlab = "Type de prélèvement", ylab = "Fraction de cellules T", 
         main = sprintf( "Wilcoxon-test :\n p-value = %s", round(p_2, 6) ),
         #col=c("brown1", "#00AFBB"),
         cex.lab=4, cex.axis=3, cex.main=3, ylim=c(0, 0.16) )
# legend("topright", legend=c("Biopsie","Cytobrosse"), col=c("#00AFBB","brown1"),
#        pt.cex=2, pch=15, border = "white", cex = 2, bty = "n")
dev.off()

mean_batch1 <- mean( as.numeric(subset_b1_2$TCRA.tcell.fraction) )
subset_b1_2["Class"] <- ifelse( subset_b1_2$TCRA.tcell.fraction < mean_batch1, "cold", "hot" )
tab_2 <- table( subset_b1_2$Sample_type, subset_b1_2$Class )
png(file = "./images_stats/4-Distribution_de_class_selon_type_echantillon_dans_batch1.png", width = 900, height = 800)
mosaicplot(tab_2, main="",xlab="Sample type",ylab="Class",col=c("grey85", "grey55"), cex.axis = 1)
legend(x='center', fill=c("grey85", "grey55"), legend=c("cold : TCRA-T-cell-fraction < Batch1 mean", "hot : TCRA-T-cell-fraction >= Batch1 mean"),
       title=sprintf("Batch1 mean = %s", round(mean_batch1, 5)), cex = 1)
dev.off()

########
# Correlation entre fraction de cellules T et et AF des HLA
variants_filtered <- read.table(file = "~/ZINARA/Gnomic/1-StageM2/5-Stats/1-text/Apres_filtre_avec_toutes_les_colonnes.txt",
                                header = TRUE, sep = "\t")
Annotation_reduced <- read.table(file = "~/ZINARA/Gnomic/1-StageM2/5-Stats/1-text/Annotation_reduit_avec_genes.txt",
                                 header = TRUE, sep = "\t")
Annotation_reduced["AF"] <- unlist( lapply( 1:nrow(Annotation_reduced), 
                                            function(i) { as.numeric(variants_filtered$AF[which(variants_filtered$ID==Annotation_reduced$Uploaded_variation[i])]) } ) )
Annotation_reduced["Sample"] <- unlist( lapply( 1:nrow(Annotation_reduced), function(i) { unlist(strsplit(Annotation_reduced$Uploaded_variation[i], "_"))[5] } ) )
Annotation_reduced_HLA_A <- subset(Annotation_reduced, Gene=="HLA-A")
Annotation_reduced_HLA_B <- subset(Annotation_reduced, Gene=="HLA-B")
Annotation_reduced_HLA_C <- subset(Annotation_reduced, Gene=="HLA-C")
data_HLA_A <- data_
data_HLA_B <- data_
data_HLA_C <- data_
data_HLA_A["mean_AF"] <- unlist( lapply(1:nrow(data_), 
                                          function(i) { mean(Annotation_reduced_HLA_A[Annotation_reduced_HLA_A$Sample==data_HLA_A$sample[i],]$AF) } ) )
data_HLA_B["mean_AF"] <- unlist( lapply(1:nrow(data_), 
                                        function(i) { mean(Annotation_reduced_HLA_B[Annotation_reduced_HLA_B$Sample==data_HLA_B$sample[i],]$AF) } ) )
data_HLA_C["mean_AF"] <- unlist( lapply(1:nrow(data_), 
                                        function(i) { mean(Annotation_reduced_HLA_C[Annotation_reduced_HLA_C$Sample==data_HLA_C$sample[i],]$AF) } ) )
data_HLA_A_withoutNA <- na.omit(data_HLA_A)
data_HLA_B_withoutNA <- na.omit(data_HLA_B)
data_HLA_C_withoutNA <- na.omit(data_HLA_C)
cor(data_HLA_A_withoutNA$TCRA.tcell.fraction, data_HLA_A_withoutNA$mean_AF)
cor(data_HLA_B_withoutNA$TCRA.tcell.fraction, data_HLA_B_withoutNA$mean_AF)
cor(data_HLA_C_withoutNA$TCRA.tcell.fraction, data_HLA_C_withoutNA$mean_AF)

library(stats)
reg1 <- lm(data_HLA_A_withoutNA$mean_AF~data_HLA_A_withoutNA$TCRA.tcell.fraction)
sum1 <- summary(reg1)
p1 <- sum1$coefficients[2,4]
png(file =  "./images_stats/5-Moyenne_AF_selon_cell_T_pour_HLA_A.png", width = 900, height = 800)
par(mar=c(5, 8.3, 5, 1))
plot(data_HLA_A_withoutNA$mean_AF~data_HLA_A_withoutNA$TCRA.tcell.fraction, pch=19,
     xlab="Fraction de cellules T", ylab="Moyenne de AF de chaque \n échantillon pour HLA-A",
     cex.lab=3, cex.axis=3, cex.main=3, cex=3, col="blue", main = sprintf("p-value = %s", round(p1, 4)))
abline(reg1, col="blue", lwd=5)
dev.off()

reg2 <- lm(data_HLA_B_withoutNA$mean_AF~data_HLA_B_withoutNA$TCRA.tcell.fraction)
sum2 <- summary(reg2)
p2 <- sum2$coefficients[2,4]
png(file =  "./images_stats/6-Moyenne_AF_selon_cell_T_pour_HLA_B.png", width = 900, height = 800)
par(mar=c(5, 8.3, 5, 1))
plot(data_HLA_B_withoutNA$mean_AF~data_HLA_B_withoutNA$TCRA.tcell.fraction, pch=19,
     xlab="Fraction de cellules T", ylab="Moyenne de AF de chaque \n échantillon pour HLA-B",
     cex.lab=3, cex.axis=3, cex.main=3, cex=3, col="blue", main = sprintf("p-value = %s", round(p2, 4)))
abline(reg2, col="blue", lwd=5)
dev.off()

reg3 <- lm(data_HLA_C_withoutNA$mean_AF~data_HLA_C_withoutNA$TCRA.tcell.fraction)
sum3 <- summary(reg3)
p3 <- sum3$coefficients[2,4]
png(file =  "./images_stats/7-Moyenne_AF_selon_cell_T_pour_HLA_C.png", width = 900, height = 800)
par(mar=c(5, 8.3, 5, 1))
plot(data_HLA_C_withoutNA$mean_AF~data_HLA_C_withoutNA$TCRA.tcell.fraction, pch=19,
     xlab="Fraction de cellules T", ylab="Moyenne de AF de chaque \n échantillon pour HLA-C",
     cex.lab=3, cex.axis=3, cex.main=3, cex=3, col="blue", main = sprintf("p-value = %s", round(p3, 4)))
abline(reg3, col="blue", lwd=5)
dev.off()

#######
# HLA-A
mutated <- data_HLA_A[!is.nan(data_HLA_A$mean_AF),]$TCRA.tcell.fraction
nonmutated <- data_HLA_A[is.nan(data_HLA_A$mean_AF),]$TCRA.tcell.fraction
HLA_A_mutated <- data.frame(Type=rep("Muté", nrow(data_HLA_A[!is.nan(data_HLA_A$mean_AF),])), Fraction=data_HLA_A[!is.nan(data_HLA_A$mean_AF),]$TCRA.tcell.fraction)
HLA_A_nonmutated <- data.frame(Type=rep("Non muté", nrow(data_HLA_A[is.nan(data_HLA_A$mean_AF),])), Fraction=data_HLA_A[is.nan(data_HLA_A$mean_AF),]$TCRA.tcell.fraction)
HLA_A <- rbind(HLA_A_mutated, HLA_A_nonmutated)
test_ <- wilcox.test(HLA_A[HLA_A$Type=="Muté",]$Fraction, HLA_A[HLA_A$Type=="Non muté",]$Fraction, paired = FALSE)
p <- test_$p.value
png(file =  "./images_stats/8-Fraction_cell_T_entre_mutE_et_non_mutE_HLA_A.png", width = 900, height = 800)
par(mar=c(5, 6, 5, 1))
boxplot(HLA_A$Fraction~HLA_A$Type, pch=19,
     xlab="", ylab="Fraction de cellule T",
     cex.lab=3, cex.axis=2, cex.main=2.5, cex=2, main = sprintf("Wilcoxon-test :\n p-value = %s", round(p, 3)))
dev.off()

# HLA-B
mutated <- data_HLA_B[!is.nan(data_HLA_B$mean_AF),]$TCRA.tcell.fraction
nonmutated <- data_HLA_B[is.nan(data_HLA_B$mean_AF),]$TCRA.tcell.fraction
HLA_B_mutated <- data.frame(Type=rep("Muté", nrow(data_HLA_B[!is.nan(data_HLA_B$mean_AF),])), Fraction=data_HLA_B[!is.nan(data_HLA_B$mean_AF),]$TCRA.tcell.fraction)
HLA_B_nonmutated <- data.frame(Type=rep("Non muté", nrow(data_HLA_B[is.nan(data_HLA_B$mean_AF),])), Fraction=data_HLA_B[is.nan(data_HLA_B$mean_AF),]$TCRA.tcell.fraction)
HLA_B <- rbind(HLA_B_mutated, HLA_B_nonmutated)
test_ <- wilcox.test(HLA_B[HLA_B$Type=="Muté",]$Fraction, HLA_B[HLA_B$Type=="Non muté",]$Fraction, paired = FALSE)
p <- test_$p.value
png(file =  "./images_stats/9-Fraction_cell_T_entre_mutE_et_non_mutE_HLA_B.png", width = 900, height = 800)
par(mar=c(5, 6, 5, 1))
boxplot(HLA_B$Fraction~HLA_B$Type, pch=19,
        xlab="", ylab="Fraction de cellule T",
        cex.lab=3, cex.axis=2, cex.main=2.5, cex=2, main = sprintf("Wilcoxon-test :\n p-value = %s", round(p, 3)))
dev.off()

# HLA-C
mutated <- data_HLA_C[!is.nan(data_HLA_C$mean_AF),]$TCRA.tcell.fraction
nonmutated <- data_HLA_C[is.nan(data_HLA_C$mean_AF),]$TCRA.tcell.fraction
HLA_C_mutated <- data.frame(Type=rep("Muté", nrow(data_HLA_C[!is.nan(data_HLA_C$mean_AF),])), Fraction=data_HLA_C[!is.nan(data_HLA_C$mean_AF),]$TCRA.tcell.fraction)
HLA_C_nonmutated <- data.frame(Type=rep("Non muté", nrow(data_HLA_C[is.nan(data_HLA_C$mean_AF),])), Fraction=data_HLA_C[is.nan(data_HLA_C$mean_AF),]$TCRA.tcell.fraction)
HLA_C <- rbind(HLA_C_mutated, HLA_C_nonmutated)
test_ <- wilcox.test(HLA_C[HLA_C$Type=="Muté",]$Fraction, HLA_C[HLA_C$Type=="Non muté",]$Fraction, paired = FALSE)
p <- test_$p.value
png(file =  "./images_stats/10-Fraction_cell_T_entre_mutE_et_non_mutE_HLA_C.png", width = 900, height = 800)
par(mar=c(5, 6, 5, 1))
boxplot(HLA_C$Fraction~HLA_C$Type, pch=19,
        xlab="", ylab="Fraction de cellule T",
        cex.lab=3, cex.axis=2, cex.main=2.5, cex=2, main = sprintf("Wilcoxon-test :\n p-value = %s", round(p, 3)))
dev.off()