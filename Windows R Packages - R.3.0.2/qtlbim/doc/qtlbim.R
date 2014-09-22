### R code from vignette source 'qtlbim.Rnw'

###################################################
### code chunk number 1: qtlbim.Rnw:48-55
###################################################
# Make width of chunks 60.
options(width=80)
if(!file.exists("qtlbimPDF")) {
  dir.create("qtlbimPDF")
  warning(paste("Creating Sweave directory qtlbim"),
    call. = FALSE, immediate. = TRUE)
}


###################################################
### code chunk number 2: LoadQtlBim
###################################################
library(qtlbim)


###################################################
### code chunk number 3: LoadQtlBim
###################################################
data(qbHyper)
hyper <- qb.cross(qbHyper)


###################################################
### code chunk number 4: SummaryQbMCMC
###################################################
summary(qbHyper)


###################################################
### code chunk number 5: Scanone-HyperLPD
###################################################
temp <- qb.scanone(qbHyper,type="LPD") 
# If the pdf file hasn't already been created, then draw the plot.
if(!file.exists("qtlbimPDF/FIG-Scanone-HyperLPD.pdf"))
   { 
   par(ask=FALSE)
   pdf(file="qtlbimPDF/FIG-Scanone-HyperLPD.pdf", height = 4)
   plot(temp)
   dev.off()
   }


###################################################
### code chunk number 6: Scanone-HyperLPD
###################################################
summary(temp)


###################################################
### code chunk number 7: Scanone-HyperBF
###################################################
temp <- qb.scanone(qbHyper,type="2logBF")


###################################################
### code chunk number 8: Scanone-HyperSubset
###################################################
# If the pdf file hasn't already been created, then draw the plot.
if(!file.exists("qtlbimPDF/FIG-Scanone-HyperSubset.pdf"))
   { 
   par(ask=FALSE)
   pdf(file="qtlbimPDF/FIG-Scanone-HyperSubset.pdf", height = 4)
   plot(temp,chr=c(1,4,6,15))
   dev.off()
   }


###################################################
### code chunk number 9: ScantwoSum-HyperData
###################################################
temp <- qb.scantwo(qbHyper, chr = c(4,6,15))
summary(temp, digits = 2)


###################################################
### code chunk number 10: Scantwo-HyperData
###################################################
# If the pdf file hasn't already been created, then draw the plot.
if(!file.exists("qtlbimPDF/FIG-Scantwo-HyperData.pdf"))
   { 
   par(ask=FALSE)
   pdf(file="qtlbimPDF/FIG-Scantwo-HyperData.pdf")
   plot(temp)
   dev.off()
   }


###################################################
### code chunk number 11: qb.bf
###################################################
bf <- qb.bf(qbHyper, item = "pattern")
summary(bf)


###################################################
### code chunk number 12: qtlbim.Rnw:600-603
###################################################
bf.pat <- summary(bf)$pattern
post.pat <- row.names(bf.pat)[which.max(bf.pat$posterior)]
bf.pat <- row.names(bf.pat)[which.max(bf.pat$bf)]


###################################################
### code chunk number 13: qb.best
###################################################
best <- qb.best(qbHyper)
summary(best)


###################################################
### code chunk number 14: bestplot
###################################################
file <- paste("qtlbimPDF", "/bestmds.pdf", sep = "")
pdf(file = file, paper = "special", width = 9, height = 6)
plot(best)
invisible(dev.off())
cat("\\includegraphics{", file, "}\n\n", sep = "")


###################################################
### code chunk number 15: besthcplot
###################################################
file <- paste("qtlbimPDF", "/besthc.pdf", sep = "")
pdf(file = file, paper = "special", width = 9, height = 6)
plot(best, type = "hclust")
invisible(dev.off())
cat("\\includegraphics{", file, "}\n\n", sep = "")


###################################################
### code chunk number 16: average
###################################################
qb.best(qbHyper, include = "all")$model[[1]]
qb.best(qbHyper, include = "nested")$model[[1]]
qb.best(qbHyper, include = "exact")$model[[1]]


###################################################
### code chunk number 17: target
###################################################
target <- best$model[[1]]


###################################################
### code chunk number 18: qb.close
###################################################
close <- qb.close(qbHyper, target)
summary(close)


###################################################
### code chunk number 19: closeplot
###################################################
file <- paste("qtlbimPDF", "/closepat.pdf", sep = "")
pdf(file = file, paper = "special", width = 9, height = 6)
plot(close)
invisible(dev.off())
cat("\\includegraphics{", file, "}\n\n", sep = "")


###################################################
### code chunk number 20: closenqtlplot
###################################################
file <- paste("qtlbimPDF", "/closenqtl.pdf", sep = "")
pdf(file = file, paper = "special", width = 9, height = 6)
plot(close, category = "nqtl")
invisible(dev.off())
cat("\\includegraphics{", file, "}\n\n", sep = "")


###################################################
### code chunk number 21: Cross.Arch
###################################################
hyper.arch <- qb.arch(best)
hyper.arch


###################################################
### code chunk number 22: InitializeFitQTL
###################################################
hyper.sub <- subset(hyper, chr = hyper.arch$qtl$chr)
n.draws <- 8
hyper.sub <- sim.geno(hyper.sub, n.draws=n.draws, step=2, error=0.01)
qtl <- makeqtl(hyper.sub, as.character(hyper.arch$qtl$chr), hyper.arch$qtl$pos)


###################################################
### code chunk number 23: StepwiseFitQTL
###################################################
hyper.step <- step.fitqtl(hyper.sub, qtl, pheno.col = 1, hyper.arch)


###################################################
### code chunk number 24: SummaryFitQTL
###################################################
sum.fit <- summary(hyper.step$fit)
print(sum.fit$result.full, quote = FALSE, na.print = "")


###################################################
### code chunk number 25: qtlbim.Rnw:817-827
###################################################
if(!is.null(sum.fit$result.drop)) {
  cat("\\begin{Schunk}\n")
  cat("\\begin{Soutput}\n")
  cat("\nDrop one QTL at a time ANOVA table:\n")
  cat("----------------------------------\n")
  printCoefmat(sum.fit$result.drop[,-6], digits = 4, cs.ind = 1, P.values = TRUE, 
               has.Pvalue = TRUE, signif.legend = FALSE)
  cat("\\end{Soutput}\n\n")
  cat("\\end{Schunk}\n")
}


###################################################
### code chunk number 26: qb.multloci
###################################################
mult <- qb.multloci(qbHyper, chr = 1)


###################################################
### code chunk number 27: multplot
###################################################
file <- paste("qtlbimPDF", "/mult.pdf", sep = "")
pdf(file = file, paper = "special", width = 9, height = 9)
plot(mult)
invisible(dev.off())
cat("\\includegraphics{", file, "}\n\n", sep = "")


###################################################
### code chunk number 28: qtlbim.Rnw:868-869
###################################################
summary(mult)


###################################################
### code chunk number 29: unmerge
###################################################
summary(mult, merge = FALSE)


###################################################
### code chunk number 30: multmergeplot
###################################################
file <- paste("qtlbimPDF", "/multmerge.pdf", sep = "")
pdf(file = file, paper = "special", width = 9, height = 9)
plot(mult, merge = FALSE)
invisible(dev.off())
cat("\\includegraphics{", file, "}\n\n", sep = "")


###################################################
### code chunk number 31: qb.split
###################################################
qbHyper <- qb.split.chr(qbHyper)
qb.get(qbHyper, "split.chr")


###################################################
### code chunk number 32: bf.split
###################################################
qb.bf(qbHyper, item = "pattern") 


###################################################
### code chunk number 33: bf.split
###################################################
qb.best(qbHyper)


###################################################
### code chunk number 34: bf.split
###################################################
one <- qb.scanone(qbHyper, type = "LPD") 
summary(one)


###################################################
### code chunk number 35: oneplot
###################################################
file <- paste("qtlbimPDF", "/one.pdf", sep = "")
pdf(file = file, paper = "special", width = 9, height = 9)
plot(one, chr = 1)
invisible(dev.off())
cat("\\includegraphics{", file, "}\n\n", sep = "")


###################################################
### code chunk number 36: PlotQbCoda
###################################################
# If the pdf file hasn't already been created, then draw the plot.
if(!file.exists("qtlbimPDF/FIG-QBCODA.pdf"))
   { 
   par(ask=FALSE)
   pdf(file="qtlbimPDF/FIG-QBCODA.pdf")
   plot(qb.coda(qbHyper, variables = c("nqtl","envvar")))
   dev.off()
   }


###################################################
### code chunk number 37: PlotQbLoci34
###################################################
# If the file hasn't already been created, then draw the plot.
if(!file.exists("qtlbimPDF/FIG-QBLOCI34.pdf"))
   { 
   par(ask=FALSE)
   pdf(file="qtlbimPDF/FIG-QBLOCI34.pdf")
   plot(qb.loci(subset(qbHyper, chr=c(3,4))), labels=TRUE)
   dev.off()
   }


###################################################
### code chunk number 38: PlotQbBF
###################################################
# If the file hasn't already been created, then draw the plot.
if(!file.exists("plotQBBF.pdf"))
   { 
   par(ask=FALSE)
   pdf(file="qtlbimPDF/FIG-QBBF.pdf")
   plot(qb.BayesFactor(qbHyper,cutoff.pattern=0.5))
   dev.off()
   }


###################################################
### code chunk number 39: PlotQbHPD
###################################################
# If the file hasn't already been created, then draw the plot.
if(!file.exists("plotQBHPD.pdf"))
   { 
   par(ask=FALSE)
   pdf(file="qtlbimPDF/FIG-QBHPD.pdf")
   plot(qb.hpdone(qbHyper))
   dev.off()
   }


###################################################
### code chunk number 40: PlotQbEpi
###################################################
# If the file hasn't already been created, then draw the plot.
if(!file.exists("qtlbimPDF/FIG-QBEPI.pdf"))
   { 
   par(ask=FALSE)
   pdf(file="qtlbimPDF/FIG-QBEPI.pdf")
   plot(qb.epistasis(qbHyper))
   dev.off()
   }


###################################################
### code chunk number 41: PlotQbDiag
###################################################
# If the file hasn't already been created, then draw the plot.
if(!file.exists("plotQBDIAG.pdf"))
   { 
   par(ask=FALSE)
   pdf(file="qtlbimPDF/FIG-QBDIAG.pdf")
   plot(qb.diag(qbHyper))
   dev.off()
   }


###################################################
### code chunk number 42: SimData
###################################################
cross <- qb.sim.cross(len=rep(100,20), n.mar=11, eq.spacing=F, n.ind=100, type="f2", 
                 ordinal=c(0.3,0.3,0.2,0.2), missing.geno=0.03, missing.pheno=0.07,
                 qtl.pos=rbind(c(1,15),c(1,45),c(3,12),c(5,15),c(7,15),c(10,15),c(12,35),c(19,15)),
                 qtl.main=rbind(c(1,0.5,0),c(2,0,0.7),c(3,-0.5,0),c(4,0.5,-0.5)),
                 qtl.epis=rbind(c(4,5,-0.7,0,0,0),c(6,8,0,1.2,0,0)),
                 covariate=c(0.5,0.07), gbye=rbind(c(7,0.8,0))) 



###################################################
### code chunk number 43: qtlbim.Rnw:1220-1223
###################################################

names(cross)



###################################################
### code chunk number 44: qtlbim.Rnw:1231-1234
###################################################

summary(cross$qtl)



###################################################
### code chunk number 45: qtlbim.Rnw:1239-1240
###################################################
summary(cross)


###################################################
### code chunk number 46: qtlbim.Rnw:1392-1393
###################################################
summary(qb.scanone(qbHyper,type="variance",scan="env"))


###################################################
### code chunk number 47: qtlbim.Rnw:1403-1404
###################################################
summary(qb.scanone(qbHyper,type="heritability"))


###################################################
### code chunk number 48: qtlbim.Rnw:1459-1460
###################################################
summary(qb.scanone(qbHyper,type="2logBF"))


