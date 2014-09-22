### R code from vignette source 'permutations.Rnw'

###################################################
### code chunk number 1: preliminary
###################################################
options("prompt" = "R> ", "continue" = "+ ")
options(useFancyQuotes="UTF-8")


###################################################
### code chunk number 2: load_jackal
###################################################
require(permute)
data(jackal)
jackal


###################################################
### code chunk number 3: ttest_jackal
###################################################
jack.t <-t.test(Length ~ Sex, data = jackal, var.equal = TRUE, alternative = "greater")
jack.t


###################################################
### code chunk number 4: ftest_jackal
###################################################
var.test(Length ~ Sex, data = jackal)
fligner.test(Length ~ Sex, data = jackal)


###################################################
### code chunk number 5: meanFun
###################################################
meanDif <- function(x, grp) {
 mean(x[grp == "Male"]) - mean(x[grp == "Female"])
}


###################################################
### code chunk number 6: randJackal
###################################################
Djackal <- numeric(length = 5000)
N <- nrow(jackal)
set.seed(42)
for(i in seq_len(length(Djackal) - 1)) {
    perm <- shuffle(N)
    Djackal[i] <- with(jackal, meanDif(Length, Sex[perm]))
}
Djackal[5000] <- with(jackal, meanDif(Length, Sex))


###################################################
### code chunk number 7: hist_jackal (eval = FALSE)
###################################################
## hist(Djackal, main = "",
##      xlab = expression("Mean difference (Male - Female) in mm"))
## rug(Djackal[5000], col = "red", lwd = 2)


###################################################
### code chunk number 8: permutations.Rnw:124-125
###################################################
(Dbig <- sum(Djackal >= Djackal[5000]))


###################################################
### code chunk number 9: permutations.Rnw:128-129
###################################################
Dbig / length(Djackal)


###################################################
### code chunk number 10: draw_hist_jackal
###################################################
hist(Djackal, main = "",
     xlab = expression("Mean difference (Male - Female) in mm"))
rug(Djackal[5000], col = "red", lwd = 2)


###################################################
### code chunk number 11: permutations.Rnw:141-142
###################################################
choose(20, 10)


###################################################
### code chunk number 12: permutations.Rnw:152-153
###################################################
args(shuffle)


###################################################
### code chunk number 13: permutations.Rnw:156-157
###################################################
str(permControl())


###################################################
### code chunk number 14: permutations.Rnw:160-165
###################################################
set.seed(2)
(r1 <- shuffle(10))
set.seed(2)
(r2 <- sample(1:10, 10, replace = FALSE))
all.equal(r1, r2)


###################################################
### code chunk number 15: permutations.Rnw:184-190
###################################################
set.seed(4)
x <- 1:10
CTRL <- permControl(within = Within(type = "series"))
perm <- shuffle(10, control = CTRL)
perm
x[perm] ## equivalent


###################################################
### code chunk number 16: permutations.Rnw:197-203
###################################################
set.seed(4)
block <- gl(3, 9)
CTRL <- permControl(strata = block,
                    within = Within(type = "grid", ncol = 3, nrow = 3))
perm <- shuffle(length(block), control = CTRL)
perm


###################################################
### code chunk number 17: permutations.Rnw:208-212
###################################################
## Original
lapply(split(1:27, block), matrix, ncol = 3)
## Shuffled
lapply(split(perm, block), matrix, ncol = 3)


###################################################
### code chunk number 18: permutations.Rnw:218-224
###################################################
set.seed(4)
CTRL <- permControl(strata = block,
                    within = Within(type = "grid", ncol = 3, nrow = 3,
                                    constant = TRUE))
perm2 <- shuffle(length(block), control = CTRL)
lapply(split(perm2, block), matrix, ncol = 3)


###################################################
### code chunk number 19: permutations.Rnw:234-238
###################################################
set.seed(4)
CTRL <- permControl(within = Within(type = "series"))
pset <- shuffleSet(10, nset = 5, control = CTRL)
pset


###################################################
### code chunk number 20: seesionInfo
###################################################
toLatex(sessionInfo())


