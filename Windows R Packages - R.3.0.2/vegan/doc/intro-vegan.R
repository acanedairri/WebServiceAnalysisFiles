### R code from vignette source 'intro-vegan.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: intro-vegan.Rnw:48-53
###################################################
par(mfrow=c(1,1))
options(width=72)
figset <- function() par(mar=c(4,4,1,1)+.1)
options(SweaveHooks = list(fig = figset))
options("prompt" = "R> ", "continue" = "+  ")


###################################################
### code chunk number 2: intro-vegan.Rnw:85-88
###################################################
library(vegan)
data(dune)
ord <- decorana(dune)


###################################################
### code chunk number 3: intro-vegan.Rnw:91-92
###################################################
ord


###################################################
### code chunk number 4: intro-vegan.Rnw:115-117
###################################################
ord <- metaMDS(dune)
ord


###################################################
### code chunk number 5: a
###################################################
plot(ord)


###################################################
### code chunk number 6: intro-vegan.Rnw:132-133
###################################################
plot(ord)


###################################################
### code chunk number 7: a
###################################################
plot(ord, type = "n")
points(ord, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow")
text(ord, display = "spec", cex=0.7, col="blue")


###################################################
### code chunk number 8: intro-vegan.Rnw:154-155
###################################################
plot(ord, type = "n")
points(ord, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow")
text(ord, display = "spec", cex=0.7, col="blue")


###################################################
### code chunk number 9: intro-vegan.Rnw:219-221
###################################################
data(dune.env)
attach(dune.env)


###################################################
### code chunk number 10: a
###################################################
plot(ord, disp="sites", type="n")
ordihull(ord, Management, col="blue")
ordiellipse(ord, Management, col=3,lwd=2)
ordispider(ord, Management, col="red", label = TRUE)
points(ord, disp="sites", pch=21, col="red", bg="yellow", cex=1.3)


###################################################
### code chunk number 11: intro-vegan.Rnw:231-232
###################################################
plot(ord, disp="sites", type="n")
ordihull(ord, Management, col="blue")
ordiellipse(ord, Management, col=3,lwd=2)
ordispider(ord, Management, col="red", label = TRUE)
points(ord, disp="sites", pch=21, col="red", bg="yellow", cex=1.3)


###################################################
### code chunk number 12: intro-vegan.Rnw:262-264
###################################################
ord.fit <- envfit(ord ~ A1 + Management, data=dune.env, perm=1000)
ord.fit


###################################################
### code chunk number 13: a
###################################################
plot(ord, dis="site")
plot(ord.fit)


###################################################
### code chunk number 14: b
###################################################
ordisurf(ord, A1, add=TRUE)


###################################################
### code chunk number 15: intro-vegan.Rnw:280-282
###################################################
plot(ord, dis="site")
plot(ord.fit)
ordisurf(ord, A1, add=TRUE)


###################################################
### code chunk number 16: intro-vegan.Rnw:302-304
###################################################
ord <- cca(dune ~ A1 + Management, data=dune.env)
ord


###################################################
### code chunk number 17: a
###################################################
plot(ord)


###################################################
### code chunk number 18: intro-vegan.Rnw:311-312
###################################################
plot(ord)


###################################################
### code chunk number 19: intro-vegan.Rnw:329-330
###################################################
cca(dune ~ ., data=dune.env)


###################################################
### code chunk number 20: intro-vegan.Rnw:339-340
###################################################
anova(ord)


###################################################
### code chunk number 21: intro-vegan.Rnw:354-355
###################################################
anova(ord, by="term", permu=200)


###################################################
### code chunk number 22: intro-vegan.Rnw:361-362
###################################################
anova(ord, by="mar")


###################################################
### code chunk number 23: a
###################################################
anova(ord, by="axis", perm=500)


###################################################
### code chunk number 24: intro-vegan.Rnw:379-381
###################################################
ord <- cca(dune ~ A1 + Management + Condition(Moisture), data=dune.env)
ord


###################################################
### code chunk number 25: intro-vegan.Rnw:386-387
###################################################
anova(ord, by="term", perm=500)


###################################################
### code chunk number 26: intro-vegan.Rnw:392-393
###################################################
anova(ord, by="term", perm=500, strata=Moisture)


###################################################
### code chunk number 27: intro-vegan.Rnw:397-398
###################################################
detach(dune.env)


