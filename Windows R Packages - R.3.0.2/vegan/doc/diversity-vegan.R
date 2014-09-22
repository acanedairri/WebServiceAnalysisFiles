### R code from vignette source 'diversity-vegan.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: diversity-vegan.Rnw:42-47
###################################################
par(mfrow=c(1,1))
options(width=72)
figset <- function() par(mar=c(4,4,1,1)+.1)
options(SweaveHooks = list(fig = figset))
options("prompt" = "R> ", "continue" = "+  ")


###################################################
### code chunk number 2: diversity-vegan.Rnw:67-69
###################################################
library(vegan)
data(BCI)


###################################################
### code chunk number 3: diversity-vegan.Rnw:87-88
###################################################
H <- diversity(BCI)


###################################################
### code chunk number 4: diversity-vegan.Rnw:95-96
###################################################
J <- H/log(specnumber(BCI))


###################################################
### code chunk number 5: diversity-vegan.Rnw:122-124
###################################################
k <- sample(nrow(BCI), 6)
R <- renyi(BCI[k,])


###################################################
### code chunk number 6: diversity-vegan.Rnw:130-131
###################################################
require(lattice, quietly=TRUE)


###################################################
### code chunk number 7: diversity-vegan.Rnw:134-135
###################################################
print(plot(R))


###################################################
### code chunk number 8: diversity-vegan.Rnw:146-147
###################################################
alpha <- fisher.alpha(BCI)


###################################################
### code chunk number 9: diversity-vegan.Rnw:183-184
###################################################
quantile(rowSums(BCI))


###################################################
### code chunk number 10: diversity-vegan.Rnw:187-188
###################################################
Srar <- rarefy(BCI, min(rowSums(BCI)))


###################################################
### code chunk number 11: diversity-vegan.Rnw:196-197
###################################################
S2 <- rarefy(BCI, 2)


###################################################
### code chunk number 12: diversity-vegan.Rnw:201-202
###################################################
all(rank(Srar) == rank(S2))


###################################################
### code chunk number 13: diversity-vegan.Rnw:207-208
###################################################
range(diversity(BCI, "simp") - (S2 -1))


###################################################
### code chunk number 14: diversity-vegan.Rnw:270-274
###################################################
data(dune)
data(dune.taxon)
taxdis <- taxa2dist(dune.taxon, varstep=TRUE)
mod <- taxondive(dune, taxdis)


###################################################
### code chunk number 15: diversity-vegan.Rnw:277-278
###################################################
plot(mod)


###################################################
### code chunk number 16: diversity-vegan.Rnw:303-305
###################################################
tr <- hclust(taxdis, "aver")
mod <- treedive(dune, tr)


###################################################
### code chunk number 17: diversity-vegan.Rnw:327-330
###################################################
k <- sample(nrow(BCI), 1)
fish <- fisherfit(BCI[k,])
fish


###################################################
### code chunk number 18: diversity-vegan.Rnw:333-334
###################################################
plot(fish)


###################################################
### code chunk number 19: diversity-vegan.Rnw:349-350
###################################################
confint(fish)


###################################################
### code chunk number 20: diversity-vegan.Rnw:374-375
###################################################
prestondistr(BCI[k,])


###################################################
### code chunk number 21: diversity-vegan.Rnw:406-408
###################################################
rad <- radfit(BCI[k,])
rad


###################################################
### code chunk number 22: diversity-vegan.Rnw:411-412
###################################################
print(radlattice(rad))


###################################################
### code chunk number 23: a
###################################################
sac <- specaccum(BCI)
plot(sac, ci.type="polygon", ci.col="yellow")


###################################################
### code chunk number 24: diversity-vegan.Rnw:480-481
###################################################
sac <- specaccum(BCI)
plot(sac, ci.type="polygon", ci.col="yellow")


###################################################
### code chunk number 25: diversity-vegan.Rnw:509-510
###################################################
ncol(BCI)/mean(specnumber(BCI)) - 1


###################################################
### code chunk number 26: diversity-vegan.Rnw:527-529
###################################################
beta <- vegdist(BCI, binary=TRUE)
mean(beta)


###################################################
### code chunk number 27: diversity-vegan.Rnw:536-537
###################################################
betadiver(help=TRUE)


###################################################
### code chunk number 28: diversity-vegan.Rnw:555-557
###################################################
z <- betadiver(BCI, "z")
quantile(z)


###################################################
### code chunk number 29: diversity-vegan.Rnw:567-572
###################################################
data(dune)
data(dune.env)
z <- betadiver(dune, "z")
mod <- with(dune.env, betadisper(z, Management))
mod


###################################################
### code chunk number 30: diversity-vegan.Rnw:575-576
###################################################
boxplot(mod)


###################################################
### code chunk number 31: diversity-vegan.Rnw:636-637
###################################################
specpool(BCI)


###################################################
### code chunk number 32: diversity-vegan.Rnw:642-644
###################################################
s <- sample(nrow(BCI), 25)
specpool(BCI[s,])


###################################################
### code chunk number 33: diversity-vegan.Rnw:655-656
###################################################
estimateR(BCI[k,])


###################################################
### code chunk number 34: diversity-vegan.Rnw:692-694
###################################################
veiledspec(prestondistr(BCI[k,]))
veiledspec(BCI[k,])


###################################################
### code chunk number 35: diversity-vegan.Rnw:707-708
###################################################
smo <- beals(BCI)


###################################################
### code chunk number 36: a
###################################################
j <- which(colnames(BCI) == "Ceiba.pentandra")
plot(beals(BCI, species=j, include=FALSE), BCI[,j], main="Ceiba pentandra", xlab="Probability of occurrence", ylab="Occurrence")


###################################################
### code chunk number 37: diversity-vegan.Rnw:719-720
###################################################
j <- which(colnames(BCI) == "Ceiba.pentandra")
plot(beals(BCI, species=j, include=FALSE), BCI[,j], main="Ceiba pentandra", xlab="Probability of occurrence", ylab="Occurrence")


