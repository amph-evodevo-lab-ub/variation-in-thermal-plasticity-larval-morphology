# ===================================================
# Script Name: triturus_tail_script.R
# Author: Mihajlo Milić (mihajlo.milic155@gmail.com)
# Year Created: 2025
# Description: This script conducts geometric morphometric 
# analyses to evaluate the influence of elevated temperature 
# on larval tail morphology in Triturus newts. 
# The associated results will be reported in 
# [to be inserted upon acceptance].
# ===================================================

library(geomorph)

## importing files
tps <- readland.tps(file = "triturus_temperature_tail.tps", specID = "imageID", warnmsg = TRUE)
sliders <- read.table(file = "sliders_tail.csv", header = TRUE, sep = ",")

lm.pairs <- matrix(c(1,7,2,6,10,11,3,5,8,9), nrow = 5, ncol = 2, byrow = TRUE)
lm.pairs <- as.data.frame(lm.pairs)

links <- matrix(c(1,2,2,3,3,4,6,7,5,6,4,5,8,10,4,8,9,11,4,9), nrow = 10, ncol = 2, byrow = TRUE)
links <- as.data.frame(links)

## making classifier
classifier <- as.data.frame(do.call(rbind, strsplit(dimnames(tps)[[3]], split = "_")))
classifier <- classifier[,-2]

treatment <- substr(classifier$V1, 1, 1)
treatment <- sub("t", "E", treatment) ## E = treatment group
treatment <- sub("k", "K", treatment) ## K = control group

classifier$V3 <- sub("pre", "P", classifier$V3) ## P = beginning of the experiment
classifier$V3 <- sub("posle", "K", classifier$V3) ## K = end of the experiment
time <- classifier$V3

genotype <- substr(classifier$V1, 2, 4)
genotype <- sub("hma", "HMA", genotype) ## Hybrids with  T. macedonicus mtDNA
genotype <- sub("mac", "MAC", genotype) ## Triturus macedonicus 
genotype <- sub("iva", "IVA", genotype) ## Triturus ivanbureschi
genotype <- sub("hiv", "HIV", genotype) ## Hybrids with  T. ivanbureschi mtDNA

num <- substr(classifier$V1, 5, 6)

classifier.new <- data.frame(time, treatment, genotype, num)
classifier <- classifier.new

colnames(classifier) <- c("beg/end", "cont/treat", "genotype", "num")
classifier$ind <- interaction(classifier$genotype, classifier$num)
classifier$ind.2 <- interaction(classifier$ind, classifier$`beg/end`)
classifier$inter.sp.time <- interaction(classifier$genotype, classifier$`beg/end`)
classifier$inter.sp.exp <- interaction(classifier$genotype, classifier$`cont/treat`)
classifier$inter.time.exp <- interaction(classifier$`beg/end`, classifier$`cont/treat`)
classifier$inter.time.sp.exp <- interaction(classifier$`beg/end`, classifier$genotype, classifier$`cont/treat`)

## GPA
gpa <- gpagen(A = tps, approxBE = TRUE, curves = sliders)
plotAllSpecimens(A = gpa$coords, mean = TRUE, plot_param = list(pt.bg = as.factor(classifier$genotype)), label = TRUE)
# the orientation of tail was not the same in all individuals  i.e. their landmarks appear reflected over x axis
# next steps refer to adjusting the same landmark orientation for each individual

## LM configuration editing

dim(tps)

## plotting landmarks for each individual
# raw coords
for (i in 1:dim(tps)[3]) {
  plot(tps[,,i], pch = 19, cex = 2)
  text(tps[,1,i], tps[,2,i], labels = 1:dim(tps)[1], pos = 2)
}

# gpa coords
for (i in 1:dim(tps)[3]) {
  plot(gpa$coords[,,i], pch = 19, cex = 2)
  text(gpa$coords[,1,i], gpa$coords[,2,i], labels = 1:dim(tps)[1], pos = 2)
}

# PCA - separating right and left groups
pca <- gm.prcomp(gpa$coords) 
summary(pca)
pca_plot <- plot(pca, pch = 19)
pca_plot$PC.points[,1]
sort(pca_plot$PC.points[,1])

cbind(names(pca_plot$PC.points[,1]), dimnames(gpa$coords)[[3]])
identical(names(pca_plot$PC.points[,1]), dimnames(gpa$coords)[[3]])

sep <- pca_plot$PC.points[,1]

sep <- !sep > 0

identical(names(sep), dimnames(gpa$coords)[[3]])
dim(tps[,,sep])

tps_rotate <- rotate.coords(A = tps, type = "flipY", index = sep)

for (i in 1:dim(tps_rotate)[3]) {
  plot(tps_rotate[,,i], pch = 19, cex = 2)
  text(tps_rotate[,1,i], tps_rotate[,2,i], labels = 1:dim(tps)[1], pos = 2)
} ## GOOD!!!

## GPA again
gpa <- gpagen(A = tps_rotate, approxBE = TRUE, curves = sliders)
dimnames(gpa$coords)
summary(gpa)

plotAllSpecimens(A = gpa$coords, mean = TRUE, plot_param = list(pt.bg = as.factor(classifier$genotype)), label = TRUE)

## Asymmetry analysis
bilat.sim <- bilat.symmetry(A = gpa$coords, ind = dimnames(gpa$coords)[[3]], object.sym = TRUE, 
                             land.pairs = lm.pairs, iter = 9999, seed = 1, RRPP = TRUE)
summary(bilat.sim)

# extracting symetric shape component
shape.sym.comp <- bilat.sim$symm.shape
dim(shape.sym.comp) #GOOD

identical(dimnames(gpa$coords)[[3]], dimnames(shape.sym.comp)[[3]])

## Creating gdf
gdf <- geomorph.data.frame(shape = shape.sym.comp, 
                           logcs = log(gpa$Csize),
                           genotype = classifier$genotype,
                           time = classifier$`beg/end`, 
                           treat = classifier$`cont/treat`)


### Size comparison
# Linear model
lm.size.crossed.complex <- procD.lm(logcs ~ time * treat * genotype, iter = 9999, seed = 1, data = gdf)
summary(lm.size.crossed.complex)

## pairwise comparison for size using the pairwise() function
lm.for.pw.size <- procD.lm(logcs ~ classifier$inter.time.sp.exp, iter = 9999, seed = 1, data = gdf)
summary(lm.for.pw.size)

pw1.size <- pairwise(fit = lm.for.pw.size, groups = classifier$inter.time.sp.exp)
pw1.size.dist <- summary(pw1.size, test = "dist", stat.table = TRUE, iter = 9999)

### ggplot2 package boxplot
library(ggplot2)
par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1,1))

## UNIQUE BOXPLOT

bx.df <- data.frame(size = gdf$logcs,
                    groups = factor(classifier$inter.time.sp.exp, levels = c("P.MAC.K", "K.MAC.K", "P.MAC.E", "K.MAC.E",
                                                                             "P.IVA.K", "K.IVA.K", "P.IVA.E", "K.IVA.E",
                                                                             "P.HMA.K", "K.HMA.K", "P.HMA.E", "K.HMA.E",
                                                                             "P.HIV.K", "K.HIV.K", "P.HIV.E", "K.HIV.E")))

ggplot(bx.df, aes(x = groups, y = size, fill = groups)) +
  geom_boxplot(width = 0.5)+
  labs(title = "Size comparison among and within species and thier hybrids", x = "", y = "LogCS")+
  scale_fill_manual(values = c("P.MAC.K" = "#dfc27d", "K.MAC.K" = "#a6611a" , "P.MAC.E" = "#80cdc1", "K.MAC.E" = "#018571",
                               "P.IVA.K" = "#dfc27d", "K.IVA.K" = "#a6611a" , "P.IVA.E" = "#80cdc1", "K.IVA.E" = "#018571",
                               "P.HMA.K" = "#dfc27d", "K.HMA.K" = "#a6611a" , "P.HMA.E" = "#80cdc1", "K.HMA.E" = "#018571",
                               "P.HIV.K" = "#dfc27d", "K.HIV.K" = "#a6611a" , "P.HIV.E" = "#80cdc1", "K.HIV.E" = "#018571"))+
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "gray50", size = 1)+
  geom_vline(xintercept = 8.5, linetype = "dashed", color = "gray50", size = 1)+
  geom_vline(xintercept = 12.5, linetype = "dashed", color = "gray50", size = 1)+
  theme( plot.title = element_text(hjust = 0.5))


### SHAPE ANALYSIS
## Linear model
# exp.-time.-species interaction
lm.crossed.complex <- procD.lm(shape ~ time * treat * genotype, iter = 9999, seed = 1, RRPP = TRUE, SS.type = "I", data = gdf)
summary(lm.crossed.complex)

## pairwise comparison
# of all 16 possible groups

lm.for.pw <- procD.lm(gdf$shape ~ classifier$inter.time.sp.exp, iter = 9999, seed = 1, RRPP = TRUE)
summary(lm.for.pw)

pw1 <- pairwise(fit = lm.for.pw, groups = classifier$inter.time.sp.exp)

pw1.dist <- summary(pw1, test = "dist", stat.table = TRUE, iter = 9999)
pw1.var <- summary(pw1, test = "var", stat.table = TRUE, iter = 9999)

### Trajectory (vector) analysis
lm.ta <- procD.lm(gdf$shape ~ classifier$inter.sp.exp * classifier$`beg/end`,
                  iter = 9999, seed = 1, RRPP = TRUE)
summary(lm.ta)
ta <- trajectory.analysis(fit = lm.ta, groups = classifier$inter.sp.exp, 
                            traj.pts = classifier$`beg/end`, pca = TRUE)
ta.mag <- summary(ta, attribute = "MD", iter = 9999)
ta.angle <- summary(ta, attribute = "TC", angle.type = "deg", iter = 9999)  

# plotting genotype specific vectors
pch = c(22, 21)
names(pch) = c("K", "P") # K = experiment end; P = experiment beginning
pch = pch[classifier$`beg/end`]

col = c("cyan3", "pink2", "dodgerblue4", "deeppink3",
        "cyan1", "pink", "dodgerblue1", "deeppink")

names(col) = c("HIV.E", "HMA.E", "IVA.E", "MAC.E", 
               "HIV.K", "HMA.K", "IVA.K", "MAC.K")
cols = col[classifier$inter.sp.exp]


tp <- plot(ta, pch = pch, col = cols, type = "n", main = "Temperature - Species/Hybrid-specific Phenotypic Change Vectors")

library(scales)

points(tp$pc.points[classifier$genotype == "MAC",1:2],
       pch = pch[classifier$genotype == "MAC"], 
       bg = alpha(cols[classifier$genotype == "MAC"], 0.5), 
       col = alpha("deeppink4", 0.5),
       lwd = 2,
       cex = 1)

points(tp$pc.points[classifier$genotype == "HMA",1:2],
       pch = pch[classifier$genotype == "HMA"], 
       bg = alpha(cols[classifier$genotype == "HMA"],0.5),
       col = alpha("pink3", 0.5),
       lwd = 2,
       cex = 1)

points(tp$pc.points[classifier$genotype == "IVA",1:2],
       pch = pch[classifier$genotype == "IVA"], 
       bg = alpha(cols[classifier$genotype == "IVA"],0.5),
       col = alpha("darkblue", 0.5),
       lwd = 2,
       cex = 1)

points(tp$pc.points[classifier$genotype == "HIV",1:2],
       pch = pch[classifier$genotype == "HIV"], 
       bg = alpha(cols[classifier$genotype == "HIV"],0.5),
       col = alpha("cyan4", 0.5),
       lwd = 2,
       cex = 1)

add.trajectories(tp, 
                 traj.pch = c(22, 22, 22, 22, 21, 21, 21, 21), 
                 start.bg = c("cyan3", "pink2", "dodgerblue4", "deeppink3",
                              "cyan1", "pink", "dodgerblue1", "deeppink"),
                 end.bg = c("cyan1", "pink", "dodgerblue1", "deeppink",
                            "cyan3", "pink2", "dodgerblue4", "deeppink3"), 
                 traj.lty = 1, 
                 traj.lwd = 2.5, 
                 traj.cex = 2)


legend(x = "bottomleft", 
       title = "MAC", 
       title.cex = 1.2,
       legend = c("Treatment Beginning", "Control Beginning", "Treatment End", "Control End"), 
       col = c("deeppink4", "deeppink4", "deeppink4", "deeppink4"),
       pt.bg = c("deeppink3", "deeppink", "deeppink3", "deeppink"),
       pch = c(21,21,22,22), 
       pt.cex = 2, 
       bty = "n", 
       pt.lwd = 2.5)

legend(x = "bottomright", 
       title = "IVA", 
       title.cex = 1.2,
       legend = c("Treatment Beginning", "Control Beginning", "Treatment End", "Control End"), 
       col = c("dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4"),
       pt.bg = c("dodgerblue4", "dodgerblue1", "dodgerblue4", "dodgerblue1"),
       pch = c(21,21,22,22), pt.cex = 2, bty = "n", pt.lwd = 2.5)

legend(x = "topright", 
       title = "HIV", 
       title.cex = 1.2,
       legend = c("Treatment Beginning", "Control Beginning", "Treatment End", "Control End"), 
       col = c("cyan4", "cyan4", "cyan4", "cyan4"),
       pt.bg = c("cyan3", "cyan1", "cyan3", "cyan1"),
       pch = c(21,21,22,22), 
       pt.cex = 2, 
       bty = "n", 
       pt.lwd = 2.5)

legend(x = "topleft", 
       title = "HMA", 
       title.cex = 1.2,
       legend = c("Treatment Beginning", "Control Beginning", "Treatment End", "Control End"), 
       col = c("pink3", "pink3", "pink3", "pink3"),
       pt.bg = c("pink2", "pink", "pink2", "pink"),
       pch = c(21,21,22,22), pt.cex = 2, bty = "n", pt.lwd = 2.5)

## Plotting Shape Changes

#HMA
hma_p_k <- c(mean(tp$pc.points[classifier$inter.time.sp.exp == "P.HMA.K", 1]),
             mean(tp$pc.points[classifier$inter.time.sp.exp == "P.HMA.K", 2]))

hma_k_k <- c(mean(tp$pc.points[classifier$inter.time.sp.exp == "K.HMA.K", 1]),
             mean(tp$pc.points[classifier$inter.time.sp.exp == "K.HMA.K", 2]))

hma_p_e <- c(mean(tp$pc.points[classifier$inter.time.sp.exp == "P.HMA.E", 1]),
             mean(tp$pc.points[classifier$inter.time.sp.exp == "P.HMA.E", 2]))

hma_k_e <- c(mean(tp$pc.points[classifier$inter.time.sp.exp == "K.HMA.E", 1]),
             mean(tp$pc.points[classifier$inter.time.sp.exp == "K.HMA.E", 2]))

#MAC
mac_p_k <- c(mean(tp$pc.points[classifier$inter.time.sp.exp == "P.MAC.K", 1]),
             mean(tp$pc.points[classifier$inter.time.sp.exp == "P.MAC.K", 2]))

mac_k_k <- c(mean(tp$pc.points[classifier$inter.time.sp.exp == "K.MAC.K", 1]),
             mean(tp$pc.points[classifier$inter.time.sp.exp == "K.MAC.K", 2]))

mac_p_e <- c(mean(tp$pc.points[classifier$inter.time.sp.exp == "P.MAC.E", 1]),
             mean(tp$pc.points[classifier$inter.time.sp.exp == "P.MAC.E", 2]))

mac_k_e <- c(mean(tp$pc.points[classifier$inter.time.sp.exp == "K.MAC.E", 1]),
             mean(tp$pc.points[classifier$inter.time.sp.exp == "K.MAC.E", 2]))

#HIV
hiv_p_k <- c(mean(tp$pc.points[classifier$inter.time.sp.exp == "P.HIV.K", 1]),
             mean(tp$pc.points[classifier$inter.time.sp.exp == "P.HIV.K", 2]))

hiv_k_k <- c(mean(tp$pc.points[classifier$inter.time.sp.exp == "K.HIV.K", 1]),
             mean(tp$pc.points[classifier$inter.time.sp.exp == "K.HIV.K", 2]))

hiv_p_e <- c(mean(tp$pc.points[classifier$inter.time.sp.exp == "P.HIV.E", 1]),
             mean(tp$pc.points[classifier$inter.time.sp.exp == "P.HIV.E", 2]))

hiv_k_e <- c(mean(tp$pc.points[classifier$inter.time.sp.exp == "K.HIV.E", 1]),
             mean(tp$pc.points[classifier$inter.time.sp.exp == "K.HIV.E", 2]))

#IVA
iva_p_k <- c(mean(tp$pc.points[classifier$inter.time.sp.exp == "P.IVA.K", 1]),
             mean(tp$pc.points[classifier$inter.time.sp.exp == "P.IVA.K", 2]))

iva_k_k <- c(mean(tp$pc.points[classifier$inter.time.sp.exp == "K.IVA.K", 1]),
             mean(tp$pc.points[classifier$inter.time.sp.exp == "K.IVA.K", 2]))

iva_p_e <- c(mean(tp$pc.points[classifier$inter.time.sp.exp == "P.IVA.E", 1]),
             mean(tp$pc.points[classifier$inter.time.sp.exp == "P.IVA.E", 2]))

iva_k_e <- c(mean(tp$pc.points[classifier$inter.time.sp.exp == "K.IVA.E", 1]),
             mean(tp$pc.points[classifier$inter.time.sp.exp == "K.IVA.E", 2]))

preds_hma <- shape.predictor(A = tp$trajectoy.analysis$fit$GM$fitted,
                             x = tp$pc.points[,1:2],
                             hma_p_k = hma_p_k,
                             hma_k_k = hma_k_k,
                             hma_p_e = hma_p_e,
                             hma_k_e = hma_k_e)


preds_mac <- shape.predictor(A = tp$trajectoy.analysis$fit$GM$fitted,
                             x = tp$pc.points[,1:2],
                             mac_p_k = mac_p_k,
                             mac_k_k = mac_k_k,
                             mac_p_e = mac_p_e,
                             mac_k_e = mac_k_e)


preds_hiv <- shape.predictor(A = tp$trajectoy.analysis$fit$GM$fitted,
                             x = tp$pc.points[,1:2],
                             hiv_p_k = hiv_p_k,
                             hiv_k_k = hiv_k_k,
                             hiv_p_e = hiv_p_e,
                             hiv_k_e = hiv_k_e)


preds_iva <- shape.predictor(A = tp$trajectoy.analysis$fit$GM$fitted,
                             x = tp$pc.points[,1:2],
                             iva_p_k = iva_p_k,
                             iva_k_k = iva_k_k,
                             iva_p_e = iva_p_e,
                             iva_k_e = iva_k_e)

grid.par.treat <- gridPar(pt.bg = "gray60", 
                     tar.pt.bg = "firebrick2", 
                     pt.size = 1.1, 
                     tar.pt.size = 1.1, 
                     link.col = "gray",
                     tar.link.col = "firebrick3")

grid.par.cont <- gridPar(pt.bg = "lightgoldenrod3", 
                         tar.pt.bg = "black", 
                         pt.size = 1.1, 
                         tar.pt.size = 1.1, 
                         link.col = "lightgoldenrod2",
                         tar.link.col = "gray20")

grid.par.end <- gridPar(pt.bg = "black", 
                          tar.pt.bg = "firebrick2", 
                          pt.size = 1.1, 
                          tar.pt.size = 1.1, 
                          link.col = "gray20",
                          tar.link.col = "firebrick3")


par(mfrow = c(1,3), mar = c(0,0,0,0))

##### MAC Shape Changes
# Control
plotRefToTarget(preds_mac$mac_p_k, preds_mac$mac_k_k, method = "points", 
                links = links, mag = 4, gridPars = grid.par.cont)
# Treatment
plotRefToTarget(preds_mac$mac_p_e, preds_mac$mac_k_e, method = "points", 
                links = links, mag = 4, gridPars = grid.par.treat)
# Ending stages
plotRefToTarget(preds_mac$mac_k_k, preds_mac$mac_k_e, method = "points", 
                links = links, mag = 4, gridPars = grid.par.end)

##### IVA Shape Changes
# Control
plotRefToTarget(preds_iva$iva_p_k, preds_iva$iva_k_k, method = "points", 
                links = links, mag = 4, gridPars = grid.par.cont)
# Treatment
plotRefToTarget(preds_iva$iva_p_e, preds_iva$iva_k_e, method = "points", 
                links = links, mag = 4, gridPars = grid.par.treat)
# Ending stages
plotRefToTarget(preds_iva$iva_k_k, preds_iva$iva_k_e, method = "points", 
                links = links, mag = 4, gridPars = grid.par.end)

##### HMA Shape Changes
# Control
plotRefToTarget(preds_hma$hma_p_k, preds_hma$hma_k_k, method = "points", 
                links = links, mag = 4, gridPars = grid.par.cont)
# Treatment
plotRefToTarget(preds_hma$hma_p_e, preds_hma$hma_k_e, method = "points", 
                links = links, mag = 4, gridPars = grid.par.treat)
# Ending stages
plotRefToTarget(preds_hma$hma_k_k, preds_hma$hma_k_e, method = "points", 
                links = links, mag = 4, gridPars = grid.par.end)

##### HIV Shape Changes
# Control
plotRefToTarget(preds_hiv$hiv_p_k, preds_hiv$hiv_k_k, method = "points", 
                links = links, mag = 4, gridPars = grid.par.cont)
# Treatment
plotRefToTarget(preds_hiv$hiv_p_e, preds_hiv$hiv_k_e, method = "points", 
                links = links, mag = 4 , gridPars = grid.par.treat)
# Ending stages
plotRefToTarget(preds_hiv$hiv_k_k, preds_hiv$hiv_k_e, method = "points", 
                links = links, mag = 4, gridPars = grid.par.end)


## METAMORPHOSIS - Dependence of shape data on metamorphic categories

metamorphosis.ind <- read.csv(file = "metamorphosis_individual.csv", header = TRUE, row.names = 1)
metamorphosis.ind

# adding "0" in front of one decimal numbers
meta_ids <- do.call(rbind, (strsplit(rownames(metamorphosis.ind), split = "_")))
meta_ids[,4] <- sprintf("%02d", as.numeric(meta_ids[,4]))
probe <- apply(meta_ids, 1, function(row) paste(row, collapse = "_"))
cbind(probe, rownames(metamorphosis.ind))

rownames(metamorphosis.ind) <- probe

metamorphosis.ind.chr <- metamorphosis.ind # transforming into categorical data
metamorphosis.ind.chr <- as.character(metamorphosis.ind.chr[,1])
names(metamorphosis.ind.chr) <- rownames(metamorphosis.ind)
identical(rownames(metamorphosis.ind), names(metamorphosis.ind.chr))

# adjusting to tail shape data
new_ids <- apply(classifier.new, 1, function(row) paste(row, collapse = "_"))
dimnames(gdf$shape)[[3]] <- new_ids

ids_end <- new_ids[classifier.new$time == "K"] # K = end of the experiment
length(ids_end)

x <- names(metamorphosis.ind.chr) %in% dimnames(gdf$shape)[[3]]
which(x == "FALSE")

metamorphosis.ind.chr <- metamorphosis.ind.chr[-c(30, 47, 144)]
length(metamorphosis.ind.chr)

all(names(metamorphosis.ind.chr) %in% dimnames(gdf$shape)[[3]] == TRUE)

shape.end <- gdf$shape[,,names(metamorphosis.ind.chr)]
dimnames(shape.end)[[3]]
dim(shape.end)

classifier.end <- as.data.frame(do.call(rbind, strsplit(dimnames(shape.end)[[3]], split = "_")))

gdf.end <- geomorph.data.frame(shape = shape.end, 
                               genotype = classifier.end$V3,
                               treat = classifier.end$V2, 
                               time = classifier.end$V1, 
                               metamo.chr = metamorphosis.ind.chr)

summary(gdf.end)


# Linear Models
lm.metamo.chr <- procD.lm(shape ~ metamo.chr, iter = 9999, seed = 1, RRPP = TRUE, data = gdf.end) 
summary(lm.metamo.chr)

lm.metamo.sp.chr <- procD.lm(shape ~ metamo.chr * genotype, iter = 9999, seed = 1, RRPP = TRUE, data = gdf.end) 
summary(lm.metamo.sp.chr)

## per genotype test
lm.metamo.mac <- procD.lm(gdf.end$shape[,,gdf.end$genotype == "MAC"] ~ gdf.end$metamo.chr[gdf.end$genotype == "MAC"] * gdf.end$treat[gdf.end$genotype == "MAC"], iter = 9999, seed = 1, RRPP = TRUE) 
summary(lm.metamo.mac)

lm.metamo.iva <- procD.lm(gdf.end$shape[,,gdf.end$genotype == "IVA"] ~ gdf.end$metamo.chr[gdf.end$genotype == "IVA"] * gdf.end$treat[gdf.end$genotype == "IVA"], iter = 9999, seed = 1, RRPP = TRUE) 
summary(lm.metamo.iva)

lm.metamo.hma <- procD.lm(gdf.end$shape[,,gdf.end$genotype == "HMA"] ~ gdf.end$metamo.chr[gdf.end$genotype == "HMA"] * gdf.end$treat[gdf.end$genotype == "HMA"], iter = 9999, seed = 1, RRPP = TRUE) 
summary(lm.metamo.hma)

lm.metamo.hiv <- procD.lm(gdf.end$shape[,,gdf.end$genotype == "HIV"] ~ gdf.end$metamo.chr[gdf.end$genotype == "HIV"] * gdf.end$treat[gdf.end$genotype == "HIV"], iter = 9999, seed = 1, RRPP = TRUE) 
summary(lm.metamo.hiv)


par(mfrow = c(2,2), mar = c(2,4.5,3,2))

rs.common <- plotAllometry(fit = lm.metamo.chr, size = as.numeric(gdf.end$metamo.chr), logsz = FALSE, method = "RegScore", col = "white", main = expression(italic ("T. ivanbureschi")))

# IVA control
points(x = jitter(as.numeric(gdf.end$metamo.chr[gdf.end$genotype == "IVA" & gdf.end$treat == "K"]), 0.8),
       y = rs.common$RegScore[gdf.end$genotype == "IVA" & gdf.end$treat == "K"],
       pch = 19, col = "dodgerblue1", cex = 1.7)
# IVA treatment
points(x = jitter(as.numeric(gdf.end$metamo.chr[gdf.end$genotype == "IVA" & gdf.end$treat == "E"]), 0.8),
       y = rs.common$RegScore[gdf.end$genotype == "IVA" & gdf.end$treat == "E"],
       pch = 17, col = "dodgerblue4", cex = 1.7)


rs.common <- plotAllometry(fit = lm.metamo.chr, size = as.numeric(gdf.end$metamo.chr), logsz = FALSE, method = "RegScore", col = "white", main = expression(italic ("T. macedonicus")))

# MAC control
points(x = jitter(as.numeric(gdf.end$metamo.chr[gdf.end$genotype == "MAC" & gdf.end$treat == "K"]), 0.8),
       y = rs.common$RegScore[gdf.end$genotype == "MAC" & gdf.end$treat == "K"],
       pch = 19, col = "deeppink", cex = 1.7)
# MAC treatment
points(x = jitter(as.numeric(gdf.end$metamo.chr[gdf.end$genotype == "MAC" & gdf.end$treat == "E"]), 0.8),
       y = rs.common$RegScore[gdf.end$genotype == "MAC" & gdf.end$treat == "E"],
       pch = 17, col = "deeppink3", cex = 1.7)


rs.common <- plotAllometry(fit = lm.metamo.chr, size = as.numeric(gdf.end$metamo.chr), logsz = FALSE, method = "RegScore", col = "white", main = expression("Hybrid - " * italic ("T. ivanbureschi") * " mtDNA"))

# HIV control
points(x = jitter(as.numeric(gdf.end$metamo.chr[gdf.end$genotype == "HIV" & gdf.end$treat == "K"]), 0.8),
       y = rs.common$RegScore[gdf.end$genotype == "HIV" & gdf.end$treat == "K"],
       pch = 19, col = "cyan2", cex = 1.7)
# HIV treatment
points(x = jitter(as.numeric(gdf.end$metamo.chr[gdf.end$genotype == "HIV" & gdf.end$treat == "E"]), 0.8),
       y = rs.common$RegScore[gdf.end$genotype == "HIV" & gdf.end$treat == "E"],
       pch = 17, col = "cyan3", cex = 1.7)


rs.common <- plotAllometry(fit = lm.metamo.chr, size = as.numeric(gdf.end$metamo.chr), logsz = FALSE, method = "RegScore", col = "white", main = expression("Hybrid - " * italic ("T. macedonicus") * " mtDNA"))

# HMA control
points(x = jitter(as.numeric(gdf.end$metamo.chr[gdf.end$genotype == "HMA" & gdf.end$treat == "K"]), 0.8),
       y = rs.common$RegScore[gdf.end$genotype == "HMA" & gdf.end$treat == "K"],
       pch = 19, col = "pink", cex = 1.7)
# HMA treatment
points(x = jitter(as.numeric(gdf.end$metamo.chr[gdf.end$genotype == "HMA" & gdf.end$treat == "E"]), 0.8),
       y = rs.common$RegScore[gdf.end$genotype == "HMA" & gdf.end$treat == "E"],
       pch = 17, col = "pink2", cex = 1.7)

legend("topleft", 
       legend = c("Control", "Treatment"),
       col = "black", 
       pch = c(19, 17),
       bty = "n",
       pt.cex = 2, 
       cex = 1.3)

