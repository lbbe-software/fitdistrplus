
#?fitdistr
#cat(paste(paste0("\\code{", attr(methods(class="fitdist"), "info")$generic, "}"), collapse=", "))

library(MASS)
set.seed(123)
x <- rgamma(100, shape = 5, rate = 0.1)
fitdistr(x, "gamma", list(shape = 1, rate = 0.1), lower = 0.001)


library(stats4)
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
nLL <- function(lambda) -sum(stats::dpois(y, lambda, log = TRUE))
(fit0 <- mle(nLL, start = list(lambda = 5), nobs = NROW(y)))


fx <- fitdist(x, "gamma")
fy <- fitdist(x, "exp")

setwd("share/prez/2026LJK/")

pdf(width=7, height=4, file="cdfcomp_gamma_exp.pdf")
par(mar=c(4,4,1,0.1))
cdfcomp(list(fx, fy), datacol = "grey40", fitlwd = 2, lines01 = TRUE, datapch = 21)
grid()
dev.off()

pdf(width=7, height=4, file="denscomp_gamma_exp.pdf")
par(mar=c(4,4,1,0.1))
denscomp(list(fx, fy), fitlwd = 2)
grid()
dev.off()

pdf(width=7, height=4, file="qqcomp_gamma_exp.pdf")
par(mar=c(4,4,1,0.1))
qqcomp(list(fx, fy))
grid()
dev.off()

pdf(width=7, height=4, file="ppcomp_gamma_exp.pdf")
par(mar=c(4,4,1,0.1))
ppcomp(list(fx, fy))
grid()
dev.off()



print2tex <- function (x, file, cap1, cap2, ...) 
{
  if (!inherits(x, "gofstat.fitdist")) 
    stop("Use only with 'gofstat.fitdist' objects")
  if (x$discrete) {
    if (!is.null(x$chisq)) {
      cat(capture.output(xtable::xtable(x$chisqtable, caption = cap1)), file=file, sep="\n")
      
      mm <- rbind(AIC = x$aic, BIC = x$bic)
      rownames(mm) <- c("Akaike's Information Criterion", 
                        "Bayesian Information Criterion")
      cat(capture.output(xtable::xtable(mm, caption = cap2)), file=file, sep="\n",  append = TRUE)
    }
    else cat("The sample is too small to automatically define cells for Chi-squared test \n")
  }
  else {
    mm <- rbind(KS = x$ks, CvM = x$cvm, AD = x$ad, 
                AIC = x$aic, BIC = x$bic)
    rownames(mm) <- c("Kolmogorov-Smirnov statistic", "Cramer-von Mises statistic", 
                      "Anderson-Darling statistic",
                      "Akaike's Information Criterion", "Bayesian Information Criterion")
    cat(capture.output(xtable::xtable(mm, caption = cap1)), file=file, sep="\n")
  }
}

print2tex(gofstat(list(fx, fy), fitnames = c("MLE gamma", "MLE exponential")), 
          file="gofstat_gamma_exp.tex", cap1="GoF stat. and criteria -- Gamma, exponential")

z <- rpois(100, lambda=3)
fz1 <- fitdist(z, "pois")
fz2 <- fitdist(z, "geom")

print2tex(gofstat(list(fz1, fz2), fitnames = c("MLE Poisson", "MLE geometric")), 
          file="gofstat_pois_geom.tex", 
          cap1="Chi'square table -- Poisson, geometric",
          cap2="GoF statistics -- Poisson, geometric")


library(RWsearch)
crandb_down()

p_deps_count("fitdistrplus")

png("essai.png")
p_graphF("fitdistrplus", reverse = TRUE)
dev.off()



library(mbbefd)


set.seed(123)
x <- rmbbefd(1000, 1/3, 1/10)
fx <- fitDR(x, "mbbefd")
bx <- bootDR(fx)

pdf(width=5, height=5, file="bootdist_mbbefd.pdf")
par(mar=c(4,4,1,0.1))
plot(bx, trueval=c(1/3, 1/10), enhance=TRUE, nbgrid=1000)
dev.off()


