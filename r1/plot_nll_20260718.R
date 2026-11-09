# Regenerated on the 2026-07-18 core (row-centering of U removed).  NEW file.
# NEW file (2026-07-16); does NOT modify plot_nll_20260626_fullEM.R.
# WHY: the previous version fitted OD-HUB with FOUR hubs at Q=4 (hb = 13,23,27,40) and
# labelled its panel "(c) OD-HUB (Q=4, non-negative)".  The paper uses THREE hubs at Q=3
# (Tokyo 13 / Osaka 27 / Fukuoka 40, P=44) and Appendix E explicitly REJECTS the four-hub
# Q=4 specification as producing a redundant basis -- so the figure in the appendix
# contradicted the text beside it.  Fixed here, and the nmfkc initialisation is now seeded
# (seed = 1) to match refit_empirical_fullEM.R and make the figure reproducible.
# The objective that the FULL procedure (inner block descent + outer variance
# updates) monotonically decreases is the MARGINAL negative log-likelihood
#   l = 0.5 sum_n [ log|2pi Sigma| + (y_n - X Theta a_n)' Sigma^{-1} (y_n - X Theta a_n) ],
#   Sigma = sigma2 I + tau2 X X^T,
# NOT the fixed-form penalized least squares at a moving lambda. Here we run pure
# ECM (no cap; sigma2 and tau2 both estimated) and plot l across all iterations.
suppressMessages({
  library(nlme); library(lavaan); library(nmfkc); library(readxl)
  library(quanteda); library(quanteda.textmodels)
})
source("nmfre_core20260718.R")

fit_nll <- function(Y, A, Q, theta_sign, x_update, tag) {
  res0 <- nmfkc(Y = Y, A = A, Q = Q, epsilon = 1e-8, method = "EU", seed = 1)
  f <- nmf.re(Y, A, Q = Q, X.ini = res0$X, Theta.ini = res0$C,
              theta_sign = theta_sign, x_update = x_update,
              dfU_control = "off",                       # no cap: pure variance estimation
              sigma2 = 1, tau2 = 1,
              sigma2_update = TRUE, tau2_update = TRUE,   # ECM variance M-steps
              inference = FALSE, maxit = 50000, epsilon = 1e-8,
              verbose = FALSE, post_message = FALSE, save_obj_trace = TRUE)
  nll <- f$nll_trace[is.finite(f$nll_trace)]
  inc <- if (length(nll) > 1) max(diff(nll)) else NA_real_
  cat(sprintf("%-10s iters=%5d  nll0=%.6g  nllN=%.6g  max step-increase=%+.3e  rel=%.2e\n",
              tag, length(nll), nll[1], nll[length(nll)], inc, inc / (abs(nll[1]) + 1e-8)))
  list(nll = nll, breaks = f$outer_break_iters, tag = tag, maxinc = inc, iters = length(nll))
}

load_ortho <- function() { d <- Orthodont; su <- unique(d$Subject); ag <- sort(unique(d$age))
  N <- length(su); P <- length(ag); Y <- matrix(d$distance, P, N)
  ml <- as.vector(matrix(ifelse(d$Sex=="Male",1,0), P, N)[1, ])
  list(Y=Y, A=rbind(intercept=rep(1,N), male=ml), Q=1, ts="free", xu="seminmf") }
load_hs <- function() { data(HolzingerSwineford1939); d <- HolzingerSwineford1939; d <- d[complete.cases(d),]
  d$age.rev<- -(d$ageyr+d$agemo/12); d$sex.2<-ifelse(d$sex==2,1,0); d$school.GW<-ifelse(d$school=="Grant-White",1,0)
  d[,c("id","sex","ageyr","agemo","school","grade")]<-NULL; d<-nmfkc::nmfkc.normalize(d)
  ex<-c("age.rev","sex.2","school.GW"); en<-setdiff(colnames(d),ex)
  list(Y=t(d[,en]), A=t(d[,ex]), Q=3, ts="free", xu="seminmf") }
load_od <- function() { dd<-as.data.frame(read_xlsx("OD_estat_R3_000040170612.xlsx",sheet=1,skip=2))
  dd<-dd[dd[,1]=="1"&dd[,3]=="02"&dd[,5]!="100"&dd[,7]!="100",]; M<-matrix(NA,47,47)
  dd[,5]<-as.numeric(dd[,5]); dd[,7]<-as.numeric(dd[,7])
  for(i in 1:47)for(j in 1:47)M[i,j]<-dd[which(dd[,5]==i&dd[,7]==j),9]; M<-log(1+M); hb<-c(13,27,40)   # Tokyo/Osaka/Fukuoka -- the paper's 3-hub spec (P=44)
  stopifnot(nrow(M[-hb,])==44, nrow(M[hb,])==3)
  list(Y=M[-hb,], A=M[hb,], Q=3, ts="nonneg", xu="mu") }
load_topic <- function() { co<-corpus(data_corpus_inaugural)
  tk<-tokens_remove(tokens(co),pattern=stopwords("en",source="snowball"))
  df<-dfm_trim(dfm_select(dfm(tk),min_nchar=3),min_termfreq=100); dm<-as.matrix(df); dm<-dm[,order(colSums(dm),decreasing=TRUE)]
  At<-rbind(pre1865=ifelse(co$Year<=1865,1,0), y1865_45=ifelse(co$Year>1865&co$Year<=1945,1,0), after1945=ifelse(co$Year>1945,1,0))
  Yt<-t(dm); colnames(At)<-colnames(Yt); list(Y=Yt, A=At, Q=3, ts="nonneg", xu="mu") }

dat <- list(Orthodont=load_ortho(), HS1939=load_hs(), `OD-HUB`=load_od(), Topic=load_topic())
cat("=== marginal negative log-likelihood (pure ECM, no cap) ===\n")
R <- list()
for (nm in names(dat)) { D <- dat[[nm]]; R[[nm]] <- fit_nll(D$Y, D$A, D$Q, D$ts, D$xu, nm) }

labs <- c(Orthodont="(a) Orthodont (Q=1, sign-free)", HS1939="(b) Holzinger--Swineford (Q=3, sign-free)",
          `OD-HUB`="(c) OD-HUB (Q=3, non-negative)", Topic="(d) Topic model (Q=3, non-negative)")
# The trace is a staircase, and the steps are not an artefact: within an inner loop
# (sigma2, tau2) -- hence Sigma = sigma2 I + tau2 X X^T -- are HELD FIXED, and l is
# dominated by the variance terms, so l is flat there (verified: within-inner change
# ~1e-7); l falls only when the OUTER EM M-steps update (sigma2, tau2).  We therefore
# mark the outer boundaries, so the staircase reads as the structure of the algorithm:
#   flat  = inner block-coordinate descent at fixed lambda
#   cliff = outer variance M-step
# (Verified on Orthodont: the 8 largest drops coincide 8/8 with outer_break_iters.)
# Full-height rules were tried and rejected: OD-HUB has 180 outer M-steps, which filled the
# panel and hid the curve.  The boundaries are marked as ticks on the x-axis instead --
# unobtrusive, and they still show WHERE l falls.  The annotation sits top-right, which is
# empty in every panel (all four traces descend to a low plateau).
draw <- function() {
  par(mfrow=c(2,2), mar=c(4.2,5.2,2.6,1.0), mgp=c(2.9,0.8,0), cex.lab=1.05)
  for (nm in names(R)) {
    nll <- R[[nm]]$nll; it <- seq_along(nll)
    # type="s", not "l": the trace IS piecewise constant -- l is flat within an inner loop
    # (verified: change ~1e-7) and falls only at the outer variance M-step.  With treads of
    # only 2-4 iterations (Orthodont), "l" joins successive points by straight segments and
    # renders the staircase as a diagonal; "s" draws what the numbers actually are.
    plot(it, nll, type="s", lwd=2.2, col="black",   # monochrome, as Figures 1-4
         xlab="Update (iteration)", ylab="Marginal negative log-likelihood",
         main=labs[[nm]])
    brk <- R[[nm]]$breaks
    brk <- brk[is.finite(brk) & brk >= 1 & brk <= length(nll)]
    if (length(brk)) rug(brk, side = 1, ticksize = 0.035, lwd = 0.8, col = "gray35")
    legend("topright", bty = "n", cex = 0.62, text.col = "gray25",
           legend = c(sprintf("%.0f to %.0f over %d iterations", nll[1], nll[length(nll)], length(nll)),
                      sprintf("%d outer variance M-steps (ticks)", length(brk))))
  }
}
pdf("Figure_nll_monotone_20260718.pdf", width=9, height=7); draw(); dev.off()
png("Figure_nll_monotone_20260718.png", width=1500, height=1150, res=150); draw(); dev.off()
cat("\nWrote Figure_nll_monotone_20260718.pdf / .png\n")
