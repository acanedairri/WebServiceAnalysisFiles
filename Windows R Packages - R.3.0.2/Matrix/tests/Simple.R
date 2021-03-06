#### Currently a collection of simple tests
##	(since 'Matrix' takes long to load, rather have fewer source files!)

##-------- *BEFORE* attaching Matrix: --------------------------------
str(Matrix::Matrix)# -> load the namespace
T <- new("ngTMatrix", i=0L, j=2L, Dim = c(2L,6L))
T
as(T, "CsparseMatrix")
## gave Error in asMethod(object) : could not find function ".M.classEnv"
## from  0.999375-23 to *-25

## another even shorter version of this:
n <- new("dgCMatrix")
n
## this:
m <- Matrix::Matrix(cbind(1,0,diag(x=2:4)))
m
##--------------------------------------------------------------------

library(Matrix)

source(system.file("test-tools.R", package = "Matrix"))# identical3() etc

if(interactive()) {
    options(error = recover)
} else options(Matrix.verbose = TRUE)# to show Matrix.msg()s

### Matrix() ''smartness''
(d4 <- Matrix(diag(4)))
(z4 <- Matrix(0*diag(4)))
(o4 <- Matrix(1+diag(4)))
(tr <- Matrix(cbind(1,0:1)))
(M4 <- Matrix(m4 <- cbind(0,rbind(6*diag(3),0))))
dM4 <- Matrix(M4, sparse = FALSE)
class(mN <-  Matrix(NA, 3,4)) # NA *is* logical
validObject(Matrix(NA))
bd4 <- bdiag(M4,dM4,M4)
stopifnot(isValid(o4, "dsyMatrix"),
          isValid(M4, "dtCMatrix"),
          validObject(dM4), validObject(mN),
          identical(bdiag(M4), bdiag(dM4)),
          identical(bd4@p, c(0L,0:3,3:6,6:9)),
          identical(bd4@i, c(0:2, 4:6, 8:10)), bd4@x == 6
          )
assert.EQ.mat(dM4, m4)
assert.EQ.mat(M4^M4, m4^m4)
assert.EQ.mat(mN, matrix(NA, 3,4))
assert.EQ.mat(bdiag(diag(4)), diag(4))
sL <- Matrix(, 3,4, sparse=TRUE)# -> "lgC"
trS <- Matrix(tr, sparse=TRUE)# failed in 0.9975-11
stopifnot(isValid(d4, "diagonalMatrix"),   isValid(z4,  "diagonalMatrix"),
          isValid(tr, "triangularMatrix"), isValid(trS, "triangularMatrix"),
          all(is.na(sL@x)), ## not yet:  all(is.na(sL)),
          !any(sL, na.rm=TRUE), all(!sL, na.rm=TRUE),
          validObject(Matrix(c(NA,0), 4, 3, byrow = TRUE)),
          validObject(Matrix(c(NA,0), 4, 4)),
          isValid(Matrix(c(NA,0,0,0), 4, 4), "sparseMatrix"))
I <- i1 <- I1 <- Diagonal(1)
I1[1,1] <- i1[1, ] <- I [ ,1] <- NA
stopifnot(identical3(I,i1,I1))
image(d4) # gave infinite recursion

assertError( Matrix(factor(letters)) )
mlp <- matrix(.leap.seconds)## 24 x 1 numeric matrix
Mlp <- Matrix(.leap.seconds)
assert.EQ.mat(Mlp, mlp)
if(FALSE) { ## TODO -- once R itself does better ...
    .Leap.seconds <- as.POSIXlt(.leap.seconds)
    mLp <- matrix(.Leap.seconds)##  24 x 1  list each of length 24 -- yuck!!!
    MLp <- Matrix(.Leap.seconds)## --> error (for now)
}

### Unit-diagonal and unitriangular  {methods need diagU2N() or similar}
I <- Diagonal(3)
(T <- as(I,"TsparseMatrix")) # unitriangular
(C <- as(I,"CsparseMatrix")) #   (ditto)
lT <- as(T,"lMatrix")
lC <- as(C,"lMatrix")
stopifnot(
    identical((n0 <- I != 0), Diagonal(3, TRUE)), I@diag == "U",
    identical(n0, I & TRUE), identical(n0, I | FALSE),
    identical(n0, TRUE & I), identical(n0, FALSE | I),
    all(n0 == !(I == 0)), all(I == n0), identical(n0 == I, I == n0)
    ,
    identical4(lT, as(Diagonal(3, x=TRUE),"TsparseMatrix"), T & TRUE, TRUE & T),
    identical4(lC, as(Diagonal(3, x=TRUE),"CsparseMatrix"), C & TRUE, TRUE & C),
    identical3(lT, T | FALSE, FALSE | T),
    identical3(lC, C | FALSE, FALSE | C),
    TRUE)
I[,1] <- NA; I[2,2] <- NA ; I[3,] <- NaN
stopifnot(isValid(I, "sparseMatrix"))
I # gave error in printSpMatrix() - because of R bug in format.info()

L <- spMatrix(9, 30, i = rep(1:9, 3), 1:27, (1:27) %% 4 != 1)
M <- drop0(crossprod(L))
diag(M) <- diag(M) + 5 # to make it pos.def.
M. <- M[1:12,1:12] # small ex
N3 <- as(Matrix(upper.tri(diag(3))), "nMatrix")
isValid(bdN <- bdiag(N3, N3),"nsparseMatrix")

stopifnot(identical(L, L == TRUE), ## used to give infinite recursion
          all(drop0((0 - L) != 0) == drop0(L)))
L[sample(length(L), 10)] <- NA
ll <- as(L,"logical")
stopifnot(all.equal(mean(L,  na.rm=TRUE),
		    mean(ll, na.rm=TRUE), tol= 1e-14),
	  all.equal(mean(L,  na.rm=TRUE, trim=1/4),# <- with a warning
		    mean(ll, na.rm=TRUE, trim=1/4), tol= 1e-14))


## Examples where  is.na(.) was wrong:
validObject(sc <- new("dsCMatrix", i=as.integer(c(0,0:1,1:2,0:1,3)), Dim=c(4L,4L),
                      p = c(0L,1L,3L,5L,8L), x = c(0,NA,NA,0:1,0,NA,1)))
validObject(gc <- as(sc, "generalMatrix"))

stopifnot(isSymmetric(M), isSymmetric(M.),
	  is(bdiag(M., M.),"symmetricMatrix"),
	  is(bdN, "triangularMatrix"),
          all(sc == gc | (is.na(sc) & is.na(gc))),
	  all.equal(N3,N3),
	  tail(all.equal(N3, t(N3)), 1) == all.equal(1,-1),# ~= "Mean relative difference: 2"
          all((bdN != t(bdN)) == (bdN + t(bdN))), # <nsparse> != <nsparse>  failed to work...
	  !any((0+bdN) > bdN), # <dsparse> o <nsparse>
	  !any(bdN != (0+bdN)), # <nsparse> o <dsparse>
	  length(grep("Length", all.equal(M., (vM <- as.vector(M.))))) > 0,
	  identical(M., (M2 <- Matrix(vM, 12,12))),
	  all.equal(M., M2, tol=0)
	  )

## large sparse ones: these now directly "go sparse":
str(m0 <- Matrix(0,     nrow=100, ncol = 1000))
str(l0 <- Matrix(FALSE, nrow=100, ncol = 200))
stopifnot(all(!l0),
          identical(FALSE, any(l0)))

if(!interactive()) warnings()
## really large {length(<dense equivalent>) is beyond R's limits}:
op <- options(warn = 2) # warnings (e.g. integer overflow!) become errors:
n <- 50000L
stopifnot(n^2 > .Machine$integer.max)
## had integer overflow in index constructions:
x <- 1:n
D <- Diagonal(n, x=x[n:1])
summary(D)# special method
summary(D != 0)
stopifnot(identical(x*D, (Dx <- D*x)),
	  identical(D != 0, as(D, "lMatrix")),
	  identical(Dx, local({d <- D; d@x <- d@x * x; d})))

Lrg <- new("dgTMatrix", Dim = c(n,n))
diag(Lrg[2:9,1:8]) <- 1:8
## ==:  Lrg[2:9,1:8] <- `diag<-`(Lrg[2:9,1:8], 1:8)
e1 <- try(Lrg == Lrg) # error message almost ok
e2 <- try(!Lrg) # error message was "bad", now perfect
ina <- is.na(Lrg)# "all FALSE"
stopifnot(grep("too large", e1) == 1,
          grep("too large", e2) == 1,
          !any(ina))# <- gave warning previously
stopifnot(suppressWarnings(any(Lrg)))# (double -> logical  warning)

## with dimnames:
m. <- matrix(c(0, 0, 2:0), 3, 5)
dimnames(m.) <- list(LETTERS[1:3], letters[1:5])
(m0 <- m <- Matrix(m.))
m@Dimnames[[2]] <- m@Dimnames[[1]]
## not valid anymore:
(val <- validObject(m, test=TRUE)); stopifnot(is.character(val))
dm <- as(m0, "denseMatrix")
stopifnot(all.equal(rcond(dm), rcond(m.), tol = 1e-14),
	  ##^^^^^^^ dm and m. are both dense, interestingly small differences
	  ## show in at least one case of optimized BLAS
	  all.equal(rcond(dm), 0.4899474520656),
          ## show(<dgRMatrix>) had revealed a bug in C:
          identical(capture.output(show(as(m0, "RsparseMatrix")))[-(1:2)],
                    gsub("0", ".",  capture.output(show(m.))[-1])))
rm(m)

###--  Sparse Triangular :

g5 <- new("dgCMatrix", Dim = c(5L, 5L),
          x = c(10, 1, 3, 10, 1, 10, 1, 10, 10),
          i = c(0L,2L,4L, 1L, 3L,2L,4L, 3L, 4L),
          p = c(0L, 3L, 5L, 7:9))
t5 <- as(g5, "triangularMatrix") # fine
stopifnot(class(t5) == "dtCMatrix",
          identical(t5, tril(g5)))
## This is really a regression test for 'methods::selectMethod()'
## Maybe move to R once 'Matrix' is recommended
sm <- selectMethod(coerce, c("dgCMatrix", "triangularMatrix"), verbose=TRUE)
stopifnot(identical(sm(g5), t5))


(t1 <- new("dtTMatrix", x= c(3,7), i= 0:1, j=3:2,
           Dim= as.integer(c(4,4))))
## Diagonal  o  Sparse
I4 <- Diagonal(4)
D4 <- Diagonal(4, x=1:4)
validObject(t2  <-   t1  + I4)
validObject(tt2 <- t(t1) + I4)
validObject(t1c <- as(t1, "CsparseMatrix"))
validObject(t2c <- as(t2, "CsparseMatrix"))
stopifnot(validObject(t1),
          isValid(2 * I4, "diagonalMatrix"),
          isValid(D4 * 3, "diagonalMatrix"),
          isValid(I4 / 5, "diagonalMatrix"),
          isValid(D4 / 2, "diagonalMatrix"),
          identical(t1, t(t(t1))),
          identical(t1c, t(t(t1c))),
          isValid(t1c + I4,"triangularMatrix"), isValid(t2c + I4,"triangularMatrix"),
          c(class(t2), class(t1c), class(t2c), class(tt2)) == "dtCMatrix",
          identical(t(tt2), t2))
assert.EQ.mat(t1, as(t1c, "matrix"))
D4. <- D4 * (A4 <- Matrix(1:4, 4,4))
D4p <- A4 + D4
Lg1 <- D4 > 0 & D4 > 1
nLg <- !Lg1
nnLg <- !nLg
D4m <- D4 * 4:1
assert.EQ.mat(D4., diag(x= (1:4)^2))
assert.EQ.mat(D4p, diag(x= (1:4)) + (1:4))
assert.EQ.mat(D4m, diag(x=c(4,6,6,4)))
assert.EQ.mat(Lg1, diag(x= c(FALSE, rep(TRUE,3))))
stopifnot(is(Lg1, "diagonalMatrix"), is(D4m, "diagonalMatrix"),
	  is(D4., "diagonalMatrix"),
          is(nLg, "symmetricMatrix"), is(nnLg, "symmetricMatrix"),
          identical3(Lg1, Matrix(nnLg), as(nnLg, "diagonalMatrix")),
          all(Lg1 != (!Lg1)))


## as(<diag>, <anything>) :
str(cls <- names(getClass("Matrix")@subclasses))# all Matrix classes
for(cl in cls)
    if(canCoerce(I4, cl)) {
	cat(cl,":")
	M  <- as(I4, cl)
	M. <- as(D4, cl)
        stopifnot(diag(4) == as(M,"matrix"),
                  if(is(cl,"dMatrix")) diag(x=1:4) == as(M.,"matrix") else TRUE)
	cat(" [Ok]\n")
    }
s4 <- as(D4,"sparseMatrix")
v <- c(11,2,2,12); s4[2:3,2:3] <- v; validObject(s4)
s4. <- D4; s4.[2:3,2:3] <- v; validObject(s4.)
stopifnot(all(s4 == s4.))
## now assign symmetrically to symmetricMatrix
s4 <- as(as(D4,"sparseMatrix"),"symmetricMatrix")
s4[2:3,2:3] <- v
validObject(s4)
stopifnot(is(s4,"symmetricMatrix"))
assert.EQ.mat(s4, as(s4.,"matrix"),tol=0)

## lower-triangular unit-diagonal
L <- new("dtCMatrix", i = 1L, p = c(0:1, 1L), Dim = c(2L, 2L),
         x = 0.5, uplo = "L", diag = "U")
stopifnot(range(L) == 0:1, all.equal(mean(L), 5/8))

## from  0-diagonal to unit-diagonal triangular {low-level step}:
tu <- t1 ; tu@diag <- "U"
tu
validObject(cu <- as(tu, "dtCMatrix"))
validObject(cnu <- diagU2N(cu))# <- testing diagU2N
validObject(tu. <- as(cu, "dtTMatrix"))
validObject(tt <- as(cu, "TsparseMatrix"))
stopifnot(## NOT: identical(tu, tu.), # since T* is not unique!
	  identical(cu, as(tu., "dtCMatrix")),
          length(cnu@i) == length(cu@i) + nrow(cu),
          identical(cu, diagN2U(cnu)),# <- testing diagN2U
	  all(cu >= 0, na.rm = TRUE), all(cu >= 0),
	  any(cu >= 7))
validObject(tcu <- t(cu))
validObject(ttu <- t(tu))
validObject(ltu <- as(ttu, "lMatrix"))
validObject(ldtu <- as(ltu, "denseMatrix"))
validObject(Cltu <- as(ltu, "CsparseMatrix"))
stopifnot(identical(asCsp(ttu > 0), asCsp(ltu)),
          all(ltu == as(ttu > 0,"denseMatrix")))
ltu - (ttu > 0) # failed
assert.EQ.mat(cu,  as(tu,"matrix"), tol=0)
assert.EQ.mat(cnu, as(tu,"matrix"), tol=0)

C <- suppressWarnings(Matrix(c(0,1,0,0), 5,5)) + Diagonal(5)
(tU <- diagN2U(tril(C))) # dtC Unitriangular
ntU <- as(tU, "nMatrix")
nT <- as(ntU, "TsparseMatrix")
R <- as(tU, "RsparseMatrix")
Tt <- diagU2N(R) # used to accidentally drop the diag.
stopifnot(R@x == c(1,1,1), diag(Tt) == 1)

lcu <- new("ltCMatrix", Dim = c(4L, 4L), i = c(0:1, 0L), p = c(0L, 0:3),
           x = c(TRUE, FALSE, FALSE), uplo = "U", diag = "U")
(lTu <- as(lcu,"TsparseMatrix"))# prints wrongly (in Matrix 0.999375-31)
stopifnot(identical3(rowSums(lcu), rowSums(lTu), rowSums(drop0(lcu))))
(ncu <- as(lcu, "nMatrix"))# -- gives the "pattern" of lcu, i.e. FALSE are *there*
ncn <- diagU2N(ncu)
(cncn <- crossprod(ncn))# works -> "nsCMatrix"
stopifnot(identical(ncu, as(lcu,"nsparseMatrix")),
	  identical(rowSums(ncu), c(3:1, 1L)),
	  Q.eq(ncn, ncu),
	  Q.eq(crossprod(drop0(lcu)), crossprod(lcu)),# crossprod works -> "dsCMatrix"
	  identical(crossprod(ncu), cncn),
	  Q.eq(cncn, t(ncu) %*% ncu)) #used to seg.fault

U <- new("dtCMatrix", Dim = c(6L, 6L),
	 i = c(0:1, 0L, 2:3, 1L, 4L),
	 p = c(0L,0L,0L, 2:3, 5L, 7L),
	 x = rep.int(-0.5, 7), diag = "U")
validObject(U)
U. <- solve(iU <- solve(U))#-> gave segmentation fault
stopifnot(validObject(U), ## had a case where solve(U) modified U !
	  validObject(iU),
	  validObject(U.),
	  ## no rounding error, since have iU@x * 8 is integer :
	  identical(U, diagN2U(drop0(U.))))

## <sparse> o <numeric> (of length > 1):
stopifnot(isValid(tm <- tu * 1:8, "sparseMatrix"),
          identical4(tm, cu * 1:8, 1:8 * cu, 1:8 * tu))

cu[1,2] <- tu[1,2] <- NA
mu <- as(tu,"matrix")
stopifnot(isValid(cu, "CsparseMatrix"), isValid(cu, "triangularMatrix"),
          isValid(tu, "TsparseMatrix"), isValid(tu, "triangularMatrix"),
          identical(cu * 1:8, tu * 1:8), # but are no longer triangular
          identical(cu > .1, as(tu > .1, "CsparseMatrix")),
          all(cu >= 0, na.rm=TRUE), !all(cu >= 1), is.na(all(tu >= 0)),
          ## Csparse_drop: preserves triangularity incl diag="U"
          identical(cu, .Call(Matrix:::Csparse_drop, cu, 0.))
          )
assert.EQ.mat(cu * 1:8, mu * 1:8)

ina <- is.na(as(cu,"matrix"))
## These 3 were each different (2008-03) !!
stopifnot(all(ina == is.na(cu)),
	  all(ina == is.na(as(cu,"generalMatrix"))),
	  all(ina == as(is.na(as(cu,"matrix")),"nMatrix")))


set.seed(7)
xx <- rpois(10, 50)
Samp <- function(n,size) sample(n, size, replace=TRUE)
Tn <- sparseMatrix(i=Samp(8, 50), j=Samp(9,50), x=xx, giveCsparse=FALSE)
Tn
stopifnot(xx == Tn@x,
	  max(xx) < max(Tn), 0 == min(Tn),
	  (sT <- sum(Tn)) == sum(colSums(Tn)), sT == sum(Tn@x),
	  range(Tn) == range(as(Tn, "CsparseMatrix")))

## tu. is diag "U", but tu2 not:
tu2 <- as(as(tu., "generalMatrix"), "triangularMatrix")
assert.EQ.mat(cu, mu, tol=0)
stopifnot(identical3(cu[cu > 1],  tu [tu > 1], mu [mu > 1]),
          identical3(cu <= 1, tu <= 1, as(mu <= 1, "lMatrix")),# all lgeMatrix
	  identical3(cu[cu <= 1], tu[tu <= 1], mu[mu <= 1]),
	  identical3(cu , triu(cu ), t(t(cu))),
	  identical3(tu , triu(tu ), t(t(tu))),
	  identical3(tu., triu(tu.), t(t(tu.))),
	  identical(tu2, triu(tu2)),
	  identical(tcu , tril(tcu)),
	  identical(ttu , tril(ttu)),
	  identical(t(tu), tril(t(tu)))
          )
assert.EQ.mat(triu(cu),   as.matrix(triu(as.matrix(cu))))
for(k in -1:1)
    assert.EQ.mat(tril(cu,k), as.matrix(tril(as.matrix(cu),k)))

(dtr <- Matrix(local({m <- diag(2); m[1,2] <- 3;m})))
identical(dtr, triu(dtr))
assert.EQ.mat(tril(dtr), diag(2))


(t4 <- new("dgTMatrix", i = 3:0, j = 0:3, x = rep(1,4), Dim = as.integer(c(4,4))))
c4 <- as(t4, "CsparseMatrix")
## the same but "dsT" (symmetric)
suppressWarnings(M <- Matrix(c(0, rep(c(0,0:1),4)), 4,4))# warning:.. length [13] is not ..multiple
tt <- as(M, "TsparseMatrix")
stopifnot(all.equal(triu(t4) + tril(t4), c4),
          all.equal(triu(tt) + tril(tt), c4))


###-- Numeric Dense: Crossprod & Solve

set.seed(123)
mm. <- mm <- Matrix(rnorm(500 * 150), nc = 150)
stopifnot(validObject(mm))
xpx <- crossprod(mm)
stopifnot(identical(mm, mm.),# once upon a time, mm was altered by crossprod()
          isValid(xpx, "dpoMatrix"))
str(mm) # 'dge*"
str(xpx)# 'dpo*"
xpy <- crossprod(mm, rnorm(500))
res <- solve(xpx, xpy)
str(xpx)# now with Cholesky factor
stopifnot(validObject(xpx),
          validObject(xpy),
          validObject(res))
stopifnot(all.equal(xpx %*% res, xpy, tol= 1e-12))
lp <- xpx >= 1
slp <- as(lp, "sparseMatrix")

ltlp  <-  lp[ lower.tri(lp) ]
sltlp <- slp[ lower.tri(slp) ]
dim(ij <- which(lower.tri(lp), arr.ind = TRUE))
ss <- slp[ij] # now fast (!)
stopifnot(identical4(lp[ij], ltlp, sltlp, as(lp, "matrix")[ij]),
          identical(ss, sltlp),
          isValid(lp, "lsyMatrix"), lp@uplo == "U")

###-- more solve() methods  {was ./solve.R }

## first for "dgeMatrix" and all kinds of RHS :
(m6 <- 1 + as(diag(0:5), "dgeMatrix"))
rcond(m6)
I6 <- as(diag(6), "dgeMatrix")
stopifnot(all.equal(I6, m6 %*% solve(m6)),
          all.equal(I6, solve(m6) %*% m6) )

(i6 <- solve(m6, Matrix(1:6)))
stopifnot(identical(i6, as(cbind(c(-4, rep(1,5))), "dgeMatrix")),
          identical(i6, solve(m6, 1:6)),
          identical(i6, solve(m6, matrix(1:6))),
          identical(i6, solve(m6, matrix(c(1,2,3,4,5,6))))
          )

## solve(<sparse>)
(m <- t1+ t(t1) + Diagonal(4))
i.m <- solve(as.mat(m))
I1 <- m %*% i.m
o4 <- diag(I1)
im <- solve(m)# is now sparse {not yet}
(I2 <- m %*% im)
(ms <- as(m, "symmetricMatrix"))
## solve(<sparse>, <sparse>):
s.mm <-  solve(m,m)
s.mms <- solve(m, ms)
## these now work "fully-sparse"
s.ms2 <- solve(ms, ms)
s.msm <- solve(ms, m)
I4c <- as(Matrix(diag(4),sparse=TRUE), "generalMatrix")
stopifnot(isValid(im, "Matrix"), isValid(I2, "Matrix"), class(I4c) == "dgCMatrix",
          all.equal(I1, as(I2,"dgeMatrix"), tol = 1e-14),
          all.equal(diag(4), as.mat(I2), tol = 1e-12),
          all.equal(s.mm,  I2, tol = 1e-14),
          all.equal(s.mms, I2, tol = 1e-14),
          all.equal(s.ms2, s.msm, tol = 4e-15),
          all.equal(s.ms2, I4c  , tol = 4e-15),
          abs(o4 - 1) < 1e-14)

image(T125 <- kronecker(kronecker(t5,t5),t5),
      main = paste("T125:",class(T125)))
dim(T3k <- kronecker(t5,kronecker(T125, t5)))
system.time(IT3 <- solve(T3k))# incredibly fast
I. <- drop0(zapsmall(IT3 %*% T3k))
I.. <- diagN2U(I.)
I <- Diagonal(5^5)
stopifnot(isValid(IT3, "dtCMatrix"),
          ## something like the equivalent of  all(I. == Diagonal(3125)) :
          identical(as(I., "diagonalMatrix"), I),
          identical(as(I..,"diagonalMatrix"), I)
          )

###-- row- and column operations  {was ./rowcolOps.R }

set.seed(321)
(m1 <- round(Matrix(rnorm(25), 5), 2))
m1k <- Matrix(round(rnorm(1000), 2), 50, 20)
m.m <- as(m1k, "matrix")
stopifnot(all.equal(colMeans(m1k), colMeans(m.m)),
          all.equal(colSums (m1k), colSums (m.m)),
          all.equal(rowMeans(m1k), rowMeans(m.m)),
          all.equal(rowSums (m1k), rowSums (m.m)))

###-- kronecker for nonsparse uses Matrix(.):
stopifnot(isValid(kr <- kronecker(m1, m6), "Matrix"))
assert.EQ.mat(kr,
              kronecker(as(m1, "matrix"),
                        as(m6, "matrix")), tol = 0)

## sparse:
(kt1 <- kronecker(t1, tu))
kt2 <- kronecker(t1c, cu)
stopifnot(identical(Matrix:::uniq(kt1), Matrix:::uniq(kt2)))
## but kt1 and kt2, both "dgT" are different since entries are not ordered!
ktf <- kronecker(as.matrix(t1), as.matrix(tu))
if(FALSE) # FIXME? our kronecker treats "0 * NA" as "0" for structural-0
assert.EQ.mat(kt2, ktf, tol= 0)
(cs1 <- colSums(kt1))
NA.or.True <- function(x) is.na(x) | x
eq <- (cs1 == colSums(as(kt1, "matrix")))
stopifnot(NA.or.True(eq), identical(is.na(eq), is.na(cs1)))
nt1 <- as(kt1, "nMatrix") # no NA's anymore
(ng1 <- as(as(nt1, "generalMatrix"),"CsparseMatrix")) # ngC
dg1 <- as(ng1, "dMatrix")# dgC
lt1 <- kt1 > 5
nt1 <- as(lt1, "nMatrix")
(colSums(nt1, sparseResult = TRUE))
(colSums(kt1, sparseResult = TRUE)) # dsparse, with NA
(colSums(lt1, sparseResult = TRUE)) # isparse, with NA
(colSums(lt1, sparseResult = TRUE, na.rm = TRUE))
(colSums(nt1, sparseResult = TRUE)) # isparse, no NA
## check correct sparseness of both:
for(M in list(kt1, nt1, ng1, dg1, lt1, nt1)) {
    m <- as(M, "matrix")
    for(na.rm in c(FALSE,TRUE)) {
	cs  <- colSums(M, na.rm = na.rm)
	cs. <- colSums(M, na.rm = na.rm, sparseResult = TRUE)
	rs  <- rowSums(M, na.rm = na.rm)
	rs. <- rowSums(M, na.rm = na.rm, sparseResult = TRUE)
	stopifnot(isValid(cs., "sparseVector"), identical(cs, as(cs., "vector")),
                  isValid(rs., "sparseVector"), identical(rs, as(rs., "vector")),
		  {eq <- cs == colSums(m, na.rm = na.rm) ; ineq <- is.na(eq)
		   all(ineq | eq) && identical(ineq, is.na(cs)) },
		  {eq <- rs == rowSums(m, na.rm = na.rm) ; ineq <- is.na(eq)
		   all(ineq | eq) && identical(ineq, is.na(rs)) } )
    }
}

i1 <- cs. == 1
cs2 <- cs.
cs2[i1] <- 0 # failed in *-31 !!
## now *index* with a NA-sparseVector :
i2 <- i1 ; i2[3] <- NA ; li2 <- as.logical(i2)
cs3 <- cs. ;	       cs3 [i2] <- 0
v3 <- as(cs.,"vector"); v3[li2] <- 0
cs4 <- cs.	     ; cs4[li2] <- 0
stopifnot(length(i1@x) == 2, identical(li2, as(i2,"vector")),
	  identical(cs3, cs4),
	  cs3 == v3, all(as(v3, "sparseVector") == cs3)
	  ## indexing simple "numeric" with sparseVector:
	  ## see 'R_FIXME' in ../R/sparseVector.R
	  ## , identical(v3[i2], v3[li2])
	  ## TODO:
	  ## sub-assigning into simple "numeric" with sparseVector index:
	  )



M <- Matrix(c(2:0,1),2); M. <- as(M, "sparseMatrix")
(N <- as(crossprod(kronecker(diag(2), M)) > 0,
         "nMatrix"))
(L. <- as(N,"lMatrix"))
stopifnot(identical(N, as(L.,"nMatrix")),
	  identical(kronecker(	  c(1,0), M),
		    kronecker(cbind(1:0), M)))
assert.EQ.mat(kronecker(M,	      c(1,0,0)),
	      kronecker(as.matrix(M), c(1,0,0)))

## coercion from "dpo" or "dsy"
xx <- as(xpx, "dsyMatrix")
stopifnot(isSymmetric(xxS  <- as(xx,  "sparseMatrix")),
          isSymmetric(xpxS <- as(xpx, "sparseMatrix")))

tm <- matrix(0, 8,8)
tm[cbind(c(1,1,2,7,8),
         c(3,6,4,8,8))] <- c(2,-30,15,20,80)
(tM <- Matrix(tm))                ## dtC
(mM <- Matrix(m <- (tm + t(tm)))) ## dsC
mT <- as(mM, "dsTMatrix")
gC <- as(as(mT, "dgTMatrix"), "dgCMatrix")
lT <- as(Matrix(TRUE, 2,2),"TsparseMatrix")
## Check that mT, lT, and gC print properly :
pr.mT <- capture.output(mT)
pr.lT <- capture.output(lT)[-(1:2)]
nn <- unlist(strsplit(gsub(" +\\.", "", sub("^....", "", pr.mT[-(1:2)])), " "))
stopifnot(as.numeric(nn[nn != ""]) == m[m != 0],
	  identical(1:2, grep("|", pr.lT, fixed=TRUE)),
	  identical(pr.lT, capture.output(as(lT, "nMatrix"))[-(1:2)]),
          capture.output(gC)[-1] == pr.mT[-1])
assert.EQ.mat(tM, tm, tol=0)
assert.EQ.mat(gC, m,  tol=0)
assert.EQ.mat(mT, m,  tol=0)
stopifnot(isValid(mM, "dsCMatrix"), isValid(tM, "dtCMatrix")
	  , identical(mT, as(mM, "TsparseMatrix"))
	  , identical(gC, as(mM, "generalMatrix"))
	  ## coercions	general <-> symmetric
	  , identical(as(as(mM, "generalMatrix"), "symmetricMatrix"), mM)
	  , identical(as(as(mM, "dgTMatrix"),     "symmetricMatrix"), mT)
	  , identical(as(as(tM, "generalMatrix"),"triangularMatrix"), tM)
          , identical(tM + Diagonal(8), tMD <- Diagonal(8) + tM)
          , isValid(tMD, "dtCMatrix")
	  )
eM <- eigen(mM) # works thanks to base::as.matrix hack in ../R/zzz.R
stopifnot(all.equal(eM$values,
                { v <- c(162.462112512353, 30.0665927567458)
                  c(v, 15, 0, 0, 160-v[1], -15, -v[2])}, tol=1e-14))

##--- symmetric -> pos.def. needs valid test:
m5 <- Matrix(diag(5) - 1)
assertError(as(m5, "dpoMatrix"))# not pos.definite!
pm5 <- as(m5, "dspMatrix") # packed
assertError(as(pm5, "dppMatrix"))# not pos.definite!
sm <- as(Matrix(diag(5) + 1),"dspMatrix")
pm <- as(sm,"dpoMatrix")## gave infinite recursion (for a day or so)
pp <- as(pm,"dppMatrix")

x <- round(100 * crossprod(Matrix(runif(25),5)))
D <- Diagonal(5, round(1000*runif(5)))
px <- pack(x)
stopifnot(is(x, "dpoMatrix"), is(px,"dppMatrix"), is(D, "ddiMatrix"))

class(x+D)#--> now "dsyMatrix"
stopifnot(is(x+D, "symmetricMatrix"),
	  is(D+px, "dspMatrix"),
	  identical(x+D, D+x), identical(px+D, D+px), identical(pack(x-D), px-D))


tx <- tril(x)
ptx <- pack(tx)
stopifnot(is(tx, "dtrMatrix"), is(ptx, "dtpMatrix"),
          is(t(tx), "dtrMatrix"), is(t(ptx), "dtpMatrix"),
          is(D + tx, "dtrMatrix"), is(tx + D, "dtrMatrix"),
          is(ptx + D, "dtpMatrix"), is(D + ptx, "dtpMatrix"))


###-- dense nonzero pattern:
class(m <- Matrix(TRUE,2,2)) # lsy
isValid(n <- as(m, "nMatrix"), "nsyMatrix")

## 1)
as(n,"CsparseMatrix") # used to give CHOLMOD error: invalid xtype...
ls2 <- as(m, "CsparseMatrix") # works fine
## and really  'm' and 'n' are interally slot identical (!!!)

as(n,"sparseMatrix")
as(m,"sparseMatrix")

### -- now when starting with nsparse :
nT <- new("ngTMatrix",
          i = as.integer(c(0, 1, 0)),
          j = as.integer(c(0, 0, 1)), Dim = as.integer(c(2,2)))
(nC <- as(nT, "ngCMatrix"))
str(nC)# of course, no 'x' slot

tt <- as(nT,"denseMatrix") # nge (was lge "wrongly")
stopifnot(is(tt,"ngeMatrix"),
	  identical(as(tt, "lMatrix"),
		    as(as(nT, "lMatrix"), "denseMatrix")))
tt
as(nC,"denseMatrix")


###-- sparse nonzero pattern : ----------

(nkt <- as(as(as(kt1, "generalMatrix"), "CsparseMatrix"), "ngCMatrix"))# ok
dkt <- as(nkt, "denseMatrix")
(clt <- crossprod(nkt))
stopifnot(isValid(nkt, "ngCMatrix"),
          isValid(clt, "nsCMatrix"))
suppressWarnings(crossprod(clt)) ## warning "crossprod() of symmetric ..."

## a Csparse with *repeated* entry is not valid!
assertError(new("ngCMatrix", p = c(0L,2L), i = c(0L,0L), Dim = 2:1))


### "d" <-> "l"  for (symmetric) sparse : ---------------------------------------
suppressWarnings( data(KNex) ) ## may warn, as 'Matrix' is recommended
                               ## and exist more than once at check-time
mm <- KNex$mm
xpx <- crossprod(mm)
## extract nonzero pattern
nxpx <- as(xpx, "nsCMatrix")
show(nxpx) ## now ok, since subsetting works
r <- nxpx[1:2,]
lmm <- as(mm, "lgCMatrix")
nmm <- as(lmm, "nMatrix")
xlx <- crossprod(lmm)
x.x <- crossprod(nmm)
## now A = lxpx and B = xlx should be close, but not quite the same
## since <x,y> = 0 is well possible when x!=0 and y!=0 .
## However,  A[i,j] != 0 ==> B[i,j] != 0:
A <- as(as(nxpx, "lMatrix"), "TsparseMatrix")
B <- as(as(xlx,  "lMatrix"), "TsparseMatrix")
ij <- function(a) a@i + ncol(a) * a@j
stopifnot(all(ij(A) %in% ij(B)))

l3 <- upper.tri(matrix(,3,3))
stopifnot(isValid(c3 <- as(l3, "CsparseMatrix"), "CsparseMatrix"),# lgC
          is(c3, "lMatrix"))
(M <- Matrix(l3))
stopifnot(isValid(M, "ltCMatrix"),
          isValid(M2 <- M %x% M, "triangularMatrix"), # is "dtT" (why not "dtC" ?)
          dim(M2) == c(9,9), identical(M2, kronecker(M,M)))
M3 <- M %x% M2 #ok
(cM3 <- colSums(M3, sparse=TRUE))
identical(as.vector(cM3),
          as(rev(rowSums(M3, sparse=TRUE)), "vector"))
M. <- M2 %x% M # gave infinite recursion

## diagonal, sparse & interactions
stopifnot(isValid(as(Diagonal(3), "TsparseMatrix"), "TsparseMatrix"),
          isValid(X <- Diagonal(7) + 1.5 * tM[1:7,1:7], "sparseMatrix"),
          isValid(X, "triangularMatrix"),
          isValid(XX <- X - chol(crossprod(X)), "triangularMatrix"))
X
XX
XX <- as(drop0(XX), "dsCMatrix")
stopifnot(identical(XX, Matrix(0, nrow(X), ncol(X))))

M <- Matrix(m., sparse = FALSE)
(sM <- Matrix(m.))
class(dlM <- M >= 1)
stopifnot(identical(dlM, !(M < 1)),
	  isValid(sM, "sparseMatrix"),
	  isValid(dlM, "denseMatrix"))
(lM  <- as(dlM, "sparseMatrix"))
lM2 <- as(dlM, "CsparseMatrix") #-> now ok
lM0 <- Matrix:::as_Csparse(dlM)
stopifnot(identical3(lM, lM2, lM0))

selectMethod("coerce",	c("lgeMatrix", "CsparseMatrix"),
	     useInherited = c(from = TRUE, to = FALSE))

ms0 <- Matrix(c(0,1,1,0), 2,2)
ms <- as(ms0, "TsparseMatrix")
cs <- as(ms, "CsparseMatrix")
ll <- as(ms, "lMatrix")
lt <- as(ll, "lgTMatrix")
nn <- as(cs, "nsparseMatrix")
l2 <- as(cs, "lsparseMatrix")
nt <- triu(nn)
n3 <- as(nt, "lsparseMatrix")
da <- nt + t(nt)
dm <- nt * t(nt) + da
##
mnt <- as(nt, "matrix")
m <- rbind(v <- 2:3)
assert.EQ.mat(nt %*% v, mnt %*% v)
assert.EQ.mat(v %*% nt, v %*% mnt)
assert.EQ.mat( crossprod(nt, v),  crossprod(mnt,v))
assert.EQ.mat( crossprod(v, nt),  crossprod(v,mnt))
assert.EQ.mat(tcrossprod(v, nt), tcrossprod(v,mnt))
assert.EQ.mat(tcrossprod(nt, m), tcrossprod(mnt, m))
##
stopifnot(isValid(ms, "dsTMatrix"),
          as(ms0,"matrix") == as(ll, "matrix"), # coercing num |-> log
	  as(lt, "matrix") == as(ll, "matrix"),
	  identical(ms, as(ll, "dMatrix")),
	  identical4(as(ll, "CsparseMatrix"), as(cs, "lMatrix"),# lsC*
		     as(nn, "lsparseMatrix"), l2),
	  identical3(da, dm, as(cs, "generalMatrix")),		# dgC*
	  identical(as(da, "lMatrix"), as(lt, "CsparseMatrix")) # lgC*
	  )
## Dense *packed* ones:
s4 <- as(D4, "symmetricMatrix")
sp <- as(as(as(D4, "symmetricMatrix"),"denseMatrix"),"dspMatrix")
tp <- as(triu(sp),"dtpMatrix")
tpL <- as(tril(sp),"dtpMatrix")
(spL <- t(sp))
stopifnot(sp @uplo=="U", tp @uplo=="U",
	  spL@uplo=="L", tpL@uplo=="L")

## band():
n <- 4 ; m <- 6
r1 <- Matrix(1:24, n,m)
validObject(M1 <- band(r1, 0,0))
(M1 <- as(M1, "sparseMatrix"))
r2 <- Matrix(1:18, 3, 6)
stopifnot(identical(M1, bandSparse(n,m, k=0, diag = list(diag(r1)))),
	  identical(band(r2, 0,4),
		    band(r2, 0,3) + band(r2, 4,4)))
s1 <- as(r1, "sparseMatrix") # such that band(s1) is sparse, too
for(k1 in (-n):m)
    for(k2 in k1:m) {
        isValid(br1 <- band(r1, k1,k2), "ddenseMatrix")
        isValid(bs1 <- band(s1, k1,k2), "CsparseMatrix")
        stopifnot(all(r1 == s1))
    }

D. <- Diagonal(x= c(-2,3:4)); D.[lower.tri(D.)] <- 1:3 ; D.
D0 <- Diagonal(x= 0:3);       D0[upper.tri(D0)] <- 1:6 ; D0
stopifnot(all.equal(list(modulus = structure(24, logarithm = FALSE), sign = -1L),
                    unclass(determinant(D.,FALSE)), tol=1e-15),
	  det(Matrix(0,1)) == 0,
          all.equal(list(modulus = structure(0, logarithm = FALSE), sign = 1L),
                    unclass(determinant(D0,FALSE)), tol=0)
          )

### More sparseVector checks: -------------------------------
validObject(new("isparseVector"))
R <- sv <- as(D4, "sparseVector")
## dim(<sparseVector>) <- (n1,n2)  --> sparse Matrix :
dim(R) <- dim(D4)
stopifnot(isValid(sv,"sparseVector"),
	  isValid(R, "sparseMatrix"),
	  identical(D4, as(R, "diagonalMatrix")))
iv <- c(rep(0, 5), 3, 0,0,7,0,0,0)
sv <- as(iv, "sparseVector")
sv. <- as(as.integer(iv), "sparseVector")
## Note: Method with signature "numeric#sparseVector" chosen ...
(sv2 <- as(sv, "isparseVector")) ## gave error
as(sv, "zsparseVector")
stopifnot(identical(sv., sv2),
	  identical(  Matrix(sv, 3,4, byrow=TRUE),
		    t(Matrix(sv, 4,3))))
options(warn = 0)# no longer error


## "Large" sparse:
n <- 100000
m <-  50000 ; nnz <- 47
M <- spMatrix(n, m,
              i = sample(n, nnz, replace = TRUE),
              j = sample(m, nnz, replace = TRUE),
              x = round(rnorm(nnz),1))
validObject(Mv <- as(M, "sparseVector"))
validObject(Dv <- as(Diagonal(60000), "sparseVector"))
validObject(LD <- Diagonal(60000, TRUE))
validObject(Lv <- as(LD, "sparseVector"))
Dm <- Dv; dim(Dm) <- c(180000L, 20000L)
stopifnot(!doExtras || isValid(Md <- M * rowSums(M, sparseResult=TRUE), "sparseMatrix"),
	  LD@diag == "U",
          isValid(Dm, "sparseMatrix"),
	  identical(Dv, as(Dm, "sparseVector")))

p. <- new("dtCMatrix", i = c(2:3, 2L), p = c(0L, 2:3, 3L, 3L),
          Dim = c(4L, 4L), x = rep(-0.5, 3), uplo = "L", diag = "U")
assert.EQ.mat(solve(solve(p.)), as(p., "matrix"))
dimnames(p.)[[1]] <- paste(1:4)
ii <- is.na(p.)
stopifnot(all(!ii), !any(as(ii, "denseMatrix")))# used to fail

lst <- ls()
table(istri <- sapply(lst, function(.) is(get(.),"triangularMatrix")))
table(triC <- sapply(lst[istri], function(.) class(get(.))))
table(uniC <- sapply(lst[istri], function(.) get(.)@diag == "U"))
lsUtr <- lst[istri][uniC]
(di <- sapply(lsUtr, function(.) dim(get(.))))
## TODO: use %*%, crossprod(), .. on all those  4 x 4 -- and check "triangular rules"
cat('Time elapsed: ', (.pt <- proc.time()),'\n') # "stats"
##

cat("doExtras:",doExtras,"\n")
if(doExtras) {
    cat("checkMatrix() of all: \n---------\n")
    Sys.setlocale("LC_COLLATE", "C")    # to keep ls() reproducible
    for(nm in ls()) if(is(.m <- get(nm), "Matrix")) {
	cat("\n", rep("-",nchar(nm)),"\n",nm, ":\n", sep='')
	checkMatrix(.m)
    }
    cat('Time elapsed: ', proc.time() - .pt,'\n') # "stats"
}

if(!interactive()) warnings()

## Platform - and other such info -- so we find it in old saved outputs
.libPaths()
SysI <- Sys.info()
structure(Sys.info()[c(4,5,1:3)], class="simple.list")
sessionInfo()
c(Matrix = packageDescription("Matrix")$Built)
if(SysI[["sysname"]] == "Linux" && require("sfsmisc")) local({
    nn <- names(.Sc <- sfsmisc::Sys.cpuinfo())
    nn <- names(.Sc <- .Sc[nn != "flags"])
    print(.Sc[grep("\\.[0-9]$", nn, invert=TRUE)])
})
