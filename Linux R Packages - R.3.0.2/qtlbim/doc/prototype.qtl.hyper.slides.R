### R code from vignette source 'prototype.qtl.hyper.slides.Rnw'

###################################################
### code chunk number 1: Initialization
###################################################
## Initialization if not called from qb.sweave.
.qb.Package <- !exists(".qb.name")
if(.qb.Package) { ## Called by Sweave directly or in R build.
  .qb.name <- "hyper"
  require(qtl, quietly = TRUE)
  data(hyper)
  hyper <- subset(clean(hyper),
    chr = ("X" != unlist(lapply(hyper$geno, class))))
  .qb.cross <- hyper
  .qb.pheno <- 1
  .qb.niter <- 3000
  .qb.draws <- 8
  .qb.scan.type <- "2logBF"
  .qb.hpd.level <- 0.5
  .qb.threshold <- c(upper=2)
  .qb.remove <- TRUE
  .qb.SweaveFile <- system.file("doc", "prototype.qtl.hyper.slide.Rnw", package="qtlbim")
  .qb.SweaveExtra <- system.file("external", "hyper.slide.extra.Rnw", package="qtlbim")
  .qb.PDFDir <- paste(names(.qb.cross$pheno)[.qb.pheno], "PDF", sep = "")
}
.qb.pheno.name <- names(.qb.cross$pheno)[.qb.pheno]
if(!file.exists(.qb.PDFDir)) {
  dir.create(.qb.PDFDir)
  warning(paste("Creating PDF directory", .qb.PDFDir),
    call. = FALSE, immediate. = TRUE)
}
## Make sure Sweave.sty is locally available.
invisible(file.copy(file.path(R.home("share"), "texmf", "Sweave.sty"),"."))
## Assign visible names for script.
cross <- .qb.cross
pheno.col <- .qb.pheno
hpd.level <- .qb.hpd.level
scan.type <- .qb.scan.type
threshold <- .qb.threshold
n.iter <- .qb.niter
n.draws <- .qb.draws
remove.qb <- .qb.remove
if(!is.null(.qb.SweaveExtra)) {
  .qb.ExtraTex <- basename(.qb.SweaveExtra)
  .qb.ExtraTex <- substring(.qb.ExtraTex,1,nchar(.qb.ExtraTex)-4)
}


###################################################
### code chunk number 2: LoadQtlbim:CallSweave
###################################################
library(qtlbim)


###################################################
### code chunk number 3: SummaryCrossObject
###################################################
summary(cross)


###################################################
### code chunk number 4: MCMCSamples
###################################################
if(.qb.Package) {
  data(qbHyper)
  cross <- qb.cross(qbHyper)
  cross.qb <- qbHyper
  cross.qb$cross.name <- "cross"
  rm(hyper, qbHyper)
} else {
  cross <- qb.genoprob(cross,step=2)
}

## Create cross.qb if it does not exist.
if(exists("cross.qb")) {
  remove.qb <- FALSE
} else {
  cross.qb <- qb.mcmc(cross, pheno.col = pheno.col, genoupdate=TRUE,
    n.iter = n.iter, verbose=FALSE)
}


###################################################
### code chunk number 5: HPDSummary
###################################################
hpd.level
scan.type
cross.hpd <- qb.hpdone(cross.qb, hpd.level, scan.type)
sum.one <- summary(cross.hpd)
sum.one
chrs <- as.vector(sum.one[, "chr"])
pos <- sum.one[, "pos"]


###################################################
### code chunk number 6: HPDPlot
###################################################
file <- paste(.qb.PDFDir, "/slide1hpd.pdf", sep = "")
pdf(file = file, paper = "special", width = 9, height = 6)
plot(cross.hpd)
invisible(dev.off())
cat("\\includegraphics{", file, "}\n\n", sep = "")


###################################################
### code chunk number 7: ScanTwoSummary
###################################################
two <- qb.scantwo(cross.qb, chr = chrs, type = scan.type)
sum.two <- summary(two,sort="upper",threshold=threshold,
  refine = TRUE)
sum.two


###################################################
### code chunk number 8: InitialArchitecture
###################################################
cross.arch <- qb.arch(sum.two, chrs, pos)
cross.arch


###################################################
### code chunk number 9: InitializeFitQTL
###################################################
cross.sub <- subset(cross, chr = unique(cross.arch$qtl$chr))
n.draws
cross.sub <- sim.geno(cross.sub, n.draws=n.draws, step=2, error=0.01)
qtl <- makeqtl(cross.sub, as.character(cross.arch$qtl$chr), cross.arch$qtl$pos)


###################################################
### code chunk number 10: StepwiseFitQTL
###################################################
cross.step <- step.fitqtl(cross.sub, qtl, pheno.col, cross.arch)


###################################################
### code chunk number 11: SummaryFitQTL
###################################################
sum.fit <- summary(cross.step$fit)
print(sum.fit$result.full, quote = FALSE, na.print = "")


###################################################
### code chunk number 12: prototype.qtl.hyper.slides.Rnw:258-270
###################################################
if(!is.null(sum.fit$result.drop)) {
  cat("\\begin{frame}[fragile]\n")
  cat("\\frametitle{Stepwise Reduction}\n\n")
  cat("\\tiny\n\n")
  cat("\\begin{Schunk}\n")
  cat("\\begin{Soutput}\n")
  printCoefmat(sum.fit$result.drop[,-6], digits = 4, cs.ind = 1, P.values = TRUE, 
               has.Pvalue = TRUE, signif.legend = FALSE)
  cat("\\end{Soutput}\n\n")
  cat("\\end{Schunk}\n")
  cat("\\end{frame}\n\n")
}


###################################################
### code chunk number 13: FinalArchitecture
###################################################
cross.arch <- cross.step$arch
cross.arch


###################################################
### code chunk number 14: ScanTwoPlotByGroup
###################################################
## Note extra R overhead to produce an arbitrary number of plots.
if(!is.null(cross.arch$chr.by.set)) {
  for(i in names(cross.arch$chr.by.set)) {
    file <- paste(.qb.PDFDir, "/slide2LOD-", i, ".pdf", sep = "")
    pdf(file = file, paper = "special", width = 8, height = 6)
    plot(two, chr = cross.arch$chr.by.set[[i]], smooth = 3,
      col = "gray", contour = 3)
    invisible(dev.off())
    cat("\\begin{frame}[fragile]\n")
    cat("\\frametitle{2-D Plots: clique", i, "}\n\n")
    cat("\\includegraphics{", file, "}\n\n", sep = "")
    cat("\n\\end{frame}\n\n")
    warning(paste("writing", file), call. = FALSE, immediate. = TRUE)
  }
  cat("\n\n")
}



###################################################
### code chunk number 15: prototype.qtl.hyper.slides.Rnw:343-346
###################################################
if(!is.null(cross.arch$pair.by.chr))
  warning(paste("creating", nrow(cross.arch$pair.by.chr$chr), "epistatic pair plots"),
    call. = FALSE, immediate. = TRUE)


###################################################
### code chunk number 16: PairPlots
###################################################
## Note extra R overhead to produce an arbitrary number of plots.
if(!is.null(cross.arch$pair.by.chr)) {
  cross <- sim.geno(cross, step = qb.get(cross.qb, "step"))
  for(i in seq(nrow(cross.arch$pair.by.chr$chr))) {
    chri <- cross.arch$pair.by.chr$chr[i,]
    posi <- cross.arch$pair.by.chr$pos[i,]
    if(chri[1] != chri[2]) {
      file <- paste(.qb.PDFDir, "/slide-", chri[[1]], "-", chri[[2]], ".pdf",
        sep = "")
      pdf(file = file, paper = "special", width = 9, height = 6)
      tmp <- qb.slicetwo(cross.qb, chri, posi, scan.type)
      plot(tmp)
      invisible(dev.off())
      cat("\\begin{frame}[fragile]\n")
      cat("\\frametitle{Epistatic Pair", chri[[1]], "and", chri[[2]], "}\n\n")
      cat("\\includegraphics{", file, "}\n\n", sep = "")
      cat("\n\\end{frame}\n\n")
      warning(paste("writing", file), call. = FALSE, immediate. = TRUE)
    }
  }
  cat("\n\n")
}


###################################################
### code chunk number 17: UserExtraSweave
###################################################
if(!is.null(.qb.SweaveExtra)) {
  warning(paste("Running Sweave on Extra to create ",
      .qb.ExtraTex, ".tex", sep = ""),
    call. = FALSE, immediate. = TRUE)
  Sweave(.qb.SweaveExtra, quiet = TRUE)
  cat("\n\\input{", .qb.ExtraTex, "}\n\n", sep = "")
}


###################################################
### code chunk number 18: prototype.qtl.hyper.slides.Rnw:388-392
###################################################
.qb.sweave.tex <- basename(.qb.SweaveFile)
.qb.sweave.tex <- paste(substring(.qb.sweave.tex, 1,
   nchar(.qb.sweave.tex) - 4), "tex", sep = ".")
.qb.pheno.tex <- paste(.qb.pheno.name, "tex", sep = ".")


###################################################
### code chunk number 19: RemoveObjects
###################################################
remove.qb
if(remove.qb) {
  qb.remove(cross.qb)
  rm(cross, cross.sub, pheno.col, threshold, n.iter, n.draws, remove.qb)
}


###################################################
### code chunk number 20: prototype.qtl.hyper.slides.Rnw:416-418
###################################################
remove(list=objects(all.names = TRUE, pattern="^\\.qb\\..*"),
       pos=".GlobalEnv")


