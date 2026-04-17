library(geomorph)

## importing files
tps <- readland.tps(file = "triturus_temperature_head.TPS", specID = "imageID", warnmsg = TRUE, negNA = TRUE)
sliders <- read.table(file = "sliders_head.txt", header = FALSE)
links <- read.table(file = "links_head.txt", header = FALSE)
lm.pairs <- read.table(file = "lm_pairs_head.txt", header = FALSE)

## making classifier
dimnames(tps)[[3]]
length(dimnames(tps)[[3]])

classifier <- as.data.frame(do.call(rbind, strsplit(dimnames(tps)[[3]], split = "_")))
colnames(classifier) <- c("beg/end", "cont/treat", "genotype", "num")
#column 1: K = end of the experiment; P = beginning of the experiment 
#column 2: E = treatment group; K = control group
#column 3: MAC – Triturus macedonicus IVA – Triturus ivanbureschi; HMA – Hybrid with T. macedonicus mtDNA; HIV – Hybrid with T. ivanbureschi mtDNA

classifier$ind <- interaction(classifier$genotype, classifier$num)
classifier$ind.2 <- interaction(classifier$ind, classifier$`beg/end`)
classifier$inter.sp.time <- interaction(classifier$genotype, classifier$`beg/end`)
classifier$inter.sp.exp <- interaction(classifier$genotype, classifier$`cont/treat`)
classifier$inter.time.exp <- interaction(classifier$`beg/end`, classifier$`cont/treat`)
classifier$inter.time.sp.exp <- interaction(classifier$`beg/end`, classifier$genotype, classifier$`cont/treat`)

na_logical <- is.na(tps)
NA_ind <- apply(is.na(tps), 3, any)
any(NA_ind)
NA_ind <- tps[,,NA_ind]
sum(na_logical)

## Estimating missing landmarks

#dataset partiotioning

k.hiv.e.tps <- tps[,,classifier$inter.time.sp.exp == "K.HIV.E"]
dim(k.hiv.e.tps) # GOOD
dimnames(k.hiv.e.tps)[[3]]
all(substr(x = dimnames(k.hiv.e.tps)[[3]], start = 1, stop = 7) == "K_E_HIV")

p.hiv.e.tps <- tps[,,classifier$inter.time.sp.exp == "P.HIV.E"]
dim(p.hiv.e.tps) # GOOD
dimnames(p.hiv.e.tps)[[3]]
all(substr(x = dimnames(p.hiv.e.tps)[[3]], start = 1, stop = 7) == "P_E_HIV")

k.hma.e.tps <- tps[,,classifier$inter.time.sp.exp == "K.HMA.E"]
dim(k.hma.e.tps) # GOOD
dimnames(k.hma.e.tps)[[3]]
all(substr(x = dimnames(k.hma.e.tps)[[3]], start = 1, stop = 7) == "K_E_HMA")

p.hma.e.tps <- tps[,,classifier$inter.time.sp.exp == "P.HMA.E"]
dim(p.hma.e.tps) # GOOD
dimnames(p.hma.e.tps)[[3]]
all(substr(x = dimnames(p.hma.e.tps)[[3]], start = 1, stop = 7) == "P_E_HMA")

k.iva.e.tps <- tps[,,classifier$inter.time.sp.exp == "K.IVA.E"]
dim(k.iva.e.tps) # GOOD
dimnames(k.iva.e.tps)[[3]]
all(substr(x = dimnames(k.iva.e.tps)[[3]], start = 1, stop = 7) == "K_E_IVA")

p.iva.e.tps <- tps[,,classifier$inter.time.sp.exp == "P.IVA.E"]
dim(p.iva.e.tps) # GOOD
dimnames(p.iva.e.tps)[[3]]
all(substr(x = dimnames(p.iva.e.tps)[[3]], start = 1, stop = 7) == "P_E_IVA")

k.mac.e.tps <- tps[,,classifier$inter.time.sp.exp == "K.MAC.E"]
dim(k.mac.e.tps) # GOOD
dimnames(k.mac.e.tps)[[3]]
all(substr(x = dimnames(k.mac.e.tps)[[3]], start = 1, stop = 7) == "K_E_MAC")

p.mac.e.tps <- tps[,,classifier$inter.time.sp.exp == "P.MAC.E"]
dim(p.mac.e.tps) # GOOD
dimnames(p.mac.e.tps)[[3]]
all(substr(x = dimnames(p.mac.e.tps)[[3]], start = 1, stop = 7) == "P_E_MAC")

k.hiv.k.tps <- tps[,,classifier$inter.time.sp.exp == "K.HIV.K"]
dim(k.hiv.k.tps) # GOOD
dimnames(k.hiv.k.tps)[[3]]
all(substr(x = dimnames(k.hiv.k.tps)[[3]], start = 1, stop = 7) == "K_K_HIV")

p.hiv.k.tps <- tps[,,classifier$inter.time.sp.exp == "P.HIV.K"]
dim(p.hiv.k.tps) # GOOD
dimnames(p.hiv.k.tps)[[3]]
all(substr(x = dimnames(p.hiv.k.tps)[[3]], start = 1, stop = 7) == "P_K_HIV")

k.hma.k.tps <- tps[,,classifier$inter.time.sp.exp == "K.HMA.K"]
dim(k.hma.k.tps) # GOOD
dimnames(k.hma.k.tps)[[3]]
all(substr(x = dimnames(k.hma.k.tps)[[3]], start = 1, stop = 7) == "K_K_HMA")

p.hma.k.tps <- tps[,,classifier$inter.time.sp.exp == "P.HMA.K"]
dim(p.hma.k.tps) # GOOD
dimnames(p.hma.k.tps)[[3]]
all(substr(x = dimnames(p.hma.k.tps)[[3]], start = 1, stop = 7) == "P_K_HMA")

k.iva.k.tps <- tps[,,classifier$inter.time.sp.exp == "K.IVA.K"]
dim(k.iva.k.tps) # GOOD
dimnames(k.iva.k.tps)[[3]]
all(substr(x = dimnames(k.iva.k.tps)[[3]], start = 1, stop = 7) == "K_K_IVA")

p.iva.k.tps <- tps[,,classifier$inter.time.sp.exp == "P.IVA.K"]
dim(p.iva.k.tps) # GOOD
dimnames(p.iva.k.tps)[[3]]
all(substr(x = dimnames(p.iva.k.tps)[[3]], start = 1, stop = 7) == "P_K_IVA")

k.mac.k.tps <- tps[,,classifier$inter.time.sp.exp == "K.MAC.K"]
dim(k.mac.k.tps) # GOOD
dimnames(k.mac.k.tps)[[3]]
all(substr(x = dimnames(k.mac.k.tps)[[3]], start = 1, stop = 7) == "K_K_MAC")

p.mac.k.tps <- tps[,,classifier$inter.time.sp.exp == "P.MAC.K"]
dim(p.mac.k.tps) # GOOD
dimnames(p.mac.k.tps)[[3]]
all(substr(x = dimnames(p.mac.k.tps)[[3]], start = 1, stop = 7) == "P_K_MAC")


# estimating missing lm

# experimental group
k.hiv.e.tps <- estimate.missing(A = k.hiv.e.tps, method = "TPS")

p.hiv.e.tps <- estimate.missing(A = p.hiv.e.tps, method = "TPS") # no missing data

k.hma.e.tps <- estimate.missing(A = k.hma.e.tps, method = "TPS")

p.hma.e.tps <- estimate.missing(A = p.hma.e.tps, method = "TPS")

k.iva.e.tps <- estimate.missing(A = k.iva.e.tps, method = "TPS")

p.iva.e.tps <- estimate.missing(A = p.iva.e.tps, method = "TPS") # no missing data

k.mac.e.tps <- estimate.missing(A = k.mac.e.tps, method = "TPS")

p.mac.e.tps <- estimate.missing(A = p.mac.e.tps, method = "TPS") # no missing data

# control group
k.hiv.k.tps <- estimate.missing(A = k.hiv.k.tps, method = "TPS") # no missing data

p.hiv.k.tps <- estimate.missing(A = p.hiv.k.tps, method = "TPS")

k.hma.k.tps <- estimate.missing(A = k.hma.k.tps, method = "TPS") # no missing data

p.hma.k.tps <- estimate.missing(A = p.hma.k.tps, method = "TPS")

k.iva.k.tps <- estimate.missing(A = k.iva.k.tps, method = "TPS")

p.iva.k.tps <- estimate.missing(A = p.iva.k.tps, method = "TPS") # no missing data

k.mac.k.tps <- estimate.missing(A = k.mac.k.tps, method = "TPS") # no missing data

p.mac.k.tps <- estimate.missing(A = p.mac.k.tps, method = "TPS") # no missing data


#binding all together

tps.estimated <- abind::abind(k.hiv.e.tps,
                              p.hiv.e.tps,
                              k.hma.e.tps,
                              p.hma.e.tps,
                              k.iva.e.tps,
                              p.iva.e.tps,
                              k.mac.e.tps,
                              p.mac.e.tps,
                              k.hiv.k.tps,
                              p.hiv.k.tps,
                              k.hma.k.tps,
                              p.hma.k.tps,
                              k.iva.k.tps,
                              p.iva.k.tps,
                              k.mac.k.tps,
                              p.mac.k.tps,
                              along = 3)

dim(tps.estimated) #GOOD!

na_logical <- is.na(tps.estimated)
NAind <- apply(is.na(tps.estimated), 3, any)
any(NAind)
NAind <- tps.estimated[,, NAind]
sum(na_logical)

identical(dimnames(tps)[[3]], dimnames(tps.estimated)[[3]]) #FALSE
cbind(dimnames(tps)[[3]], dimnames(tps.estimated)[[3]]) 

tps.estimated1 <- tps.estimated[,, match(dimnames(tps)[[3]], dimnames(tps.estimated)[[3]])]
tps.estimated <- tps.estimated1

identical(dimnames(tps)[[3]], dimnames(tps.estimated)[[3]]) #TRUE
cbind(dimnames(tps)[[3]], dimnames(tps.estimated)[[3]])
cbind(as.vector(classifier$inter.time.sp.exp), classifier$r.br., dimnames(tps.estimated)[[3]])

## GPA
gpa <- gpagen(A = tps.estimated, approxBE = TRUE, curves = sliders)
dimnames(gpa$coords)
summary(gpa)


## Finding outliers

plotOutliers(A = gpa$coords, inspect.outliers = TRUE)
plotOutliers(A = gpa$coords, groups = as.factor(classifier$genotype), inspect.outliers = TRUE)
plotOutliers(A = gpa$coords, groups = as.factor(classifier$inter.sp.time), inspect.outliers = TRUE)
plotOutliers(A = gpa$coords, groups = as.factor(classifier$inter.time.sp.exp), inspect.outliers = TRUE)

which(x = dimnames(gpa$coords)[[3]] == "K_E_MAC_6")
dimnames(gpa$coords)[[3]] [94]
gpa$coords[,,94]

gpa$coords <- gpa$coords[,,-94]
which(x = dimnames(gpa$coords)[[3]] == "K_E_MAC_6") #GOOD

dim(gpa$coords) #final data dimension

classifier <- classifier[-94,]
dim(classifier)
table(classifier$inter.time.sp.exp)
table(classifier$inter.sp.exp)

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
                           logcs = log(gpa$Csize[-94]),
                           genotype = classifier$genotype,
                           time = classifier$`beg/end`, 
                           treat = classifier$`cont/treat`)

### Size comparison

# Linear model
lm.size.crossed.complex <- procD.lm(logcs ~ time * treat * genotype, iter = 9999, seed = 1, data = gdf)
summary(lm.size.crossed.complex)

## pairwise comparison for size using pairwsie() function
lm.for.pw.size <- procD.lm(logcs ~ classifier$inter.time.sp.exp, iter = 9999, seed = 1, data = gdf)
summary(lm.for.pw.size)

pw1.size <- pairwise(fit = lm.for.pw.size, groups = classifier$inter.time.sp.exp)
pw1.size.dist <- summary(pw1.size, test = "dist", stat.table = TRUE, iter = 9999)


## UNIQUE BOXPLOT
library(ggplot2)
par(mar = c(5.1, 4.1, 4.1, 2.1))

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


### Shape comparison

# Linear model
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

ta <- trajectory.analysis(fit = lm.ta, groups = classifier$inter.sp.exp, 
                            traj.pts = classifier$`beg/end`, pca = TRUE)
ta.mag <- summary(ta, attribute = "MD", iter = 9999)
ta.angle <- summary(ta, attribute = "TC", angle.type = "deg", iter = 9999)  


# plotting species specific vectors
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


par(mfrow = c(2,2), mar = c(0,0,0,0))

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


## Genotype-specific allometry analysis
# subsetting data per genotype

# Triturus macedonicus
mac_shape <- gdf$shape[,, gdf$genotype == "MAC"]
mac_logcs <- gdf$logcs[gdf$genotype == "MAC"]
mac_time <- gdf$time[gdf$genotype == "MAC"]
mac_treat <- gdf$treat[gdf$genotype == "MAC"]
cbind(dimnames(mac_shape)[[3]], mac_time, mac_treat)

# Triturus ivanbureschi
iva_shape <- gdf$shape[,, gdf$genotype == "IVA"]
iva_logcs <- gdf$logcs[gdf$genotype == "IVA"]
iva_time <- gdf$time[gdf$genotype == "IVA"]
iva_treat <- gdf$treat[gdf$genotype == "IVA"]
cbind(dimnames(iva_shape)[[3]], iva_time, iva_treat)

# Hybrids with T. macedonicus mtDNA
hma_shape <- gdf$shape[,, gdf$genotype == "HMA"]
hma_logcs <- gdf$logcs[gdf$genotype == "HMA"]
hma_time <- gdf$time[gdf$genotype == "HMA"]
hma_treat <- gdf$treat[gdf$genotype == "HMA"]
cbind(dimnames(hma_shape)[[3]], hma_time, hma_treat)

# Hybrids with T. ivanbureschi mtDNA
hiv_shape <- gdf$shape[,, gdf$genotype == "HIV"]
hiv_logcs <- gdf$logcs[gdf$genotype == "HIV"]
hiv_time <- gdf$time[gdf$genotype == "HIV"]
hiv_treat <- gdf$treat[gdf$genotype == "HIV"]
cbind(dimnames(hiv_shape)[[3]], hiv_time, hiv_treat)


## Allometry analysis - dependency of shape plasticity on size variation

# Triturus macedonicus
lm.allo.mac <- procD.lm(mac_shape ~ mac_logcs + mac_treat * mac_time, iter = 9999, SS.type = "I")
summary(lm.allo.mac)

# Triturus ivanbureschi
lm.allo.iva <- procD.lm(iva_shape ~ iva_logcs + iva_treat * iva_time, iter = 9999, SS.type = "I")
summary(lm.allo.iva)

# Hybrids with T. macedonicus mtDNA
lm.allo.hma <- procD.lm(hma_shape ~ hma_logcs + hma_treat * hma_time, iter = 9999, SS.type = "I")
summary(lm.allo.hma)

# Hybrids with T. ivanbureschi mtDNA
lm.allo.hiv <- procD.lm(hiv_shape ~ hiv_logcs + hiv_treat * hiv_time, iter = 9999, SS.type = "I")
summary(lm.allo.hiv)



## METAMORPHOSIS - Dependence of shape data on metamorphic categories

metamorphosis.ind <- read.csv(file = "metamorphosis_individual.csv", header = TRUE, row.names = 1)
metamorphosis.ind

which(rownames(metamorphosis.ind) == "K_E_MAC_6") # removing the outlier
data <- metamorphosis.ind[-94,]
length(data)

names <- rownames(metamorphosis.ind)[-94]
length(names)

metamorphosis.ind <- data
names(metamorphosis.ind) <- names
length(metamorphosis.ind)
which(names(metamorphosis.ind) == "K_E_MAC_6") 

metamorphosis.ind.chr <- metamorphosis.ind # transforming to categorical data
metamorphosis.ind.chr <- as.character(metamorphosis.ind.chr)
names(metamorphosis.ind.chr) <- names(metamorphosis.ind)

identical(names(metamorphosis.ind), names(metamorphosis.ind.chr))

df.end <- matrix(nrow = length(metamorphosis.ind), ncol = 4)
gdf.no.shape <- gdf[-1]

# extracting the non-shape data from the end of the experiment
for (i in 1:length(gdf.no.shape)) {
  df.end[,i] <- gdf.no.shape[[i]][gdf$time == "K"] ## K = end of the experiment
}

df.end
colnames(df.end) <- names(gdf.no.shape)
all(df.end[,"time"] == "K") ## GOOD!

# extracting the shape data from the end of the experiment
shape.end <- gdf$shape[,,gdf$time == "K"]
dimnames(shape.end)[[3]]
dim(shape.end)
identical(dimnames(shape.end)[[3]], names(metamorphosis.ind)) ## GOOD!

gdf.end <- geomorph.data.frame(shape = shape.end, 
                               logcs = df.end[,"logcs"],
                               genotype = df.end[,"genotype"],
                               time = df.end[,"time"],
                               exp = df.end[,"treat"],
                               metamo.num = metamorphosis.ind,
                               metamo.chr = metamorphosis.ind.chr)
summary(gdf.end)

# Linear Models
lm.metamo.chr <- procD.lm(shape ~ metamo.chr, iter = 9999, seed = 1, RRPP = TRUE, data = gdf.end) 
summary(lm.metamo.chr)

lm.metamo.sp.chr <- procD.lm(shape ~ metamo.chr * genotype, iter = 9999, seed = 1, RRPP = TRUE, data = gdf.end) 
summary(lm.metamo.sp.chr)

## per genotype test
lm.metamo.mac <- procD.lm(gdf.end$shape[,,gdf.end$genotype == "MAC"] ~ gdf.end$metamo.chr[gdf.end$genotype == "MAC"] * gdf.end$exp[gdf.end$genotype == "MAC"], iter = 9999, seed = 1, RRPP = TRUE) 
summary(lm.metamo.mac)

lm.metamo.iva <- procD.lm(gdf.end$shape[,,gdf.end$genotype == "IVA"] ~ gdf.end$metamo.chr[gdf.end$genotype == "IVA"] * gdf.end$exp[gdf.end$genotype == "IVA"], iter = 9999, seed = 1, RRPP = TRUE) 
summary(lm.metamo.iva)

lm.metamo.hma <- procD.lm(gdf.end$shape[,,gdf.end$genotype == "HMA"] ~ gdf.end$metamo.chr[gdf.end$genotype == "HMA"] * gdf.end$exp[gdf.end$genotype == "HMA"], iter = 9999, seed = 1, RRPP = TRUE) 
summary(lm.metamo.hma)

lm.metamo.hiv <- procD.lm(gdf.end$shape[,,gdf.end$genotype == "HIV"] ~ gdf.end$metamo.chr[gdf.end$genotype == "HIV"] * gdf.end$exp[gdf.end$genotype == "HIV"], iter = 9999, seed = 1, RRPP = TRUE) 
summary(lm.metamo.hiv)

par(mfrow = c(2,2), mar = c(2,4.5,3,2))

rs.common <- plotAllometry(fit = lm.metamo.chr, size = gdf.end$metamo.num, logsz = FALSE, method = "RegScore", col = "white", main = expression(italic ("T. ivanbureschi")))

# IVA control
points(x = jitter(gdf.end$metamo.num[gdf.end$genotype == "IVA" & gdf.end$exp == "K"], 0.8),
       y = rs.common$RegScore[gdf.end$genotype == "IVA" & gdf.end$exp == "K"],
       pch = 19, col = "dodgerblue1", cex = 1.7)
# IVA treatment
points(x = jitter(gdf.end$metamo.num[gdf.end$genotype == "IVA" & gdf.end$exp == "E"], 0.8),
       y = rs.common$RegScore[gdf.end$genotype == "IVA" & gdf.end$exp == "E"],
       pch = 17, col = "dodgerblue4", cex = 1.7)


rs.common <- plotAllometry(fit = lm.metamo.chr, size = gdf.end$metamo.num, logsz = FALSE, method = "RegScore", col = "white", main = expression(italic ("T. macedonicus")))

# MAC control
points(x = jitter(gdf.end$metamo.num[gdf.end$genotype == "MAC" & gdf.end$exp == "K"], 0.8),
       y = rs.common$RegScore[gdf.end$genotype == "MAC" & gdf.end$exp == "K"],
       pch = 19, col = "deeppink", cex = 1.7)
# MAC treatment
points(x = jitter(gdf.end$metamo.num[gdf.end$genotype == "MAC" & gdf.end$exp == "E"], 0.8),
       y = rs.common$RegScore[gdf.end$genotype == "MAC" & gdf.end$exp == "E"],
       pch = 17, col = "deeppink3", cex = 1.7)


rs.common <- plotAllometry(fit = lm.metamo.chr, size = gdf.end$metamo.num, logsz = FALSE, method = "RegScore", col = "white", main = expression("Hybrid - " * italic ("T. ivanbureschi") * " mtDNA"))

# HIV control
points(x = jitter(gdf.end$metamo.num[gdf.end$genotype == "HIV" & gdf.end$exp == "K"], 0.8),
       y = rs.common$RegScore[gdf.end$genotype == "HIV" & gdf.end$exp == "K"],
       pch = 19, col = "cyan2", cex = 1.7)
# HIV treatment
points(x = jitter(gdf.end$metamo.num[gdf.end$genotype == "HIV" & gdf.end$exp == "E"], 0.8),
       y = rs.common$RegScore[gdf.end$genotype == "HIV" & gdf.end$exp == "E"],
       pch = 17, col = "cyan3", cex = 1.7)


rs.common <- plotAllometry(fit = lm.metamo.chr, size = gdf.end$metamo.num, logsz = FALSE, method = "RegScore", col = "white", main = expression("Hybrid - " * italic ("T. macedonicus") * " mtDNA"))

# HMA control
points(x = jitter(gdf.end$metamo.num[gdf.end$genotype == "HMA" & gdf.end$exp == "K"], 0.8),
       y = rs.common$RegScore[gdf.end$genotype == "HMA" & gdf.end$exp == "K"],
       pch = 19, col = "pink", cex = 1.7)
# HMA treatment
points(x = jitter(gdf.end$metamo.num[gdf.end$genotype == "HMA" & gdf.end$exp == "E"], 0.8),
       y = rs.common$RegScore[gdf.end$genotype == "HMA" & gdf.end$exp == "E"],
       pch = 17, col = "pink2", cex = 1.7)

legend("topleft", 
       legend = c("Control", "Treatment"),
       col = "black", 
       pch = c(19, 17),
       bty = "n",
       pt.cex = 2, 
       cex = 1.3)