# Empirical refits on the 2026-07-18 core (row-centering of U removed; sign-preserving
# Theta initialisation).  NEW file; refit_empirical_fullEM.R is untouched.
# Same recipe as before: nmfkc(seed=1) init, pure ECM, inference with wild_B=1000.
# Expected: Orthodont unchanged (its fixed point already had U A' = 0); HS/OD/Topic move
# slightly, and HS becomes the exact root of the profiled estimating equation.
# Re-fit ALL FOUR empirical examples under the COMPLETE-EM core
# (nmfre_core20260627_fullEM.R), PURE ECM (estimate sigma2,tau2; no cap),
# WITH wild-bootstrap inference. Emits saturation diagnostics (dfU,r,lambda,
# tau2,sigma2) + coefficient tables in the paper's
# "Variable & Basis & Estimate & SE & BSE & z & p" format.
# Mirrors emit_ecm_tables_20260626.R EXACTLY except: (1) sources the fullEM core,
# (2) OD-HUB uses the CURRENT paper spec Q=3, 3 hubs (Tokyo=13/Osaka=27/Fukuoka=40),
#     P=44 (emit_ecm_tables still had the old Q=4 / 4-hub spec).
# New file; does NOT modify any committed script.
suppressMessages({library(nlme); library(lavaan); library(nmfkc); library(readxl)
  library(quanteda); library(quanteda.textmodels)})
source("nmfre_core20260718.R")
ECM <- function(Y,A,Q,ts,xu,pfx){r0<-nmfkc(Y=Y,A=A,Q=Q,epsilon=1e-8,method="EU",seed=1)
  nmf.re(Y,A,Q=Q,X.ini=r0$X,Theta.ini=r0$C,prefix=pfx,theta_sign=ts,x_update=xu,
    dfU_control="off",sigma2=1,tau2=1,sigma2_update=TRUE,tau2_update=TRUE,
    inference=TRUE,wild_B=1000,wild_seed=123,verbose=FALSE,post_message=FALSE)}
fmtp<-function(p)ifelse(p<0.001,"$<0.001$",formatC(p,format="f",digits=4))
emit<-function(res,tag,N,P,Q){
  r<-res$dfU_frac
  cat(sprintf("\n=== %s ===  dfU=%.2f r=%.4f lambda=%.5g tau2=%.5g sigma2=%.5g conv=%s\n",
    tag,r*N*Q,r,res$sigma2/res$tau2,res$tau2,res$sigma2,res$converged))
  cf<-res$coefficients
  for(i in seq_len(nrow(cf))){v<-cf$Variable[i];b<-cf$Basis[i]
    cat(sprintf("%s & %s & %s & %s & %s & %s & %s \\\\\n",v,b,
      formatC(as.numeric(cf$Estimate[i]),format="f",digits=3),
      formatC(as.numeric(cf$SE[i]),format="f",digits=3),
      formatC(as.numeric(cf$BSE[i]),format="f",digits=3),
      formatC(as.numeric(cf$z_value[i]),format="f",digits=2),fmtp(as.numeric(cf$p_value[i]))))}
}
# Orthodont (Q=1, sign-free)
d<-Orthodont;su<-unique(d$Subject);N<-length(su);Y<-matrix(d$distance,4,N)
ml<-as.vector(matrix(ifelse(d$Sex=="Male",1,0),4,N)[1,]);A<-rbind(intercept=rep(1,N),male=ml)
emit(ECM(Y,A,1,"free","seminmf","Trend"),"ORTHODONT",27,4,1)
# HS1939 (Q=3, sign-free)
data(HolzingerSwineford1939);d<-HolzingerSwineford1939;d<-d[complete.cases(d),]
d$age.rev<- -(d$ageyr+d$agemo/12);d$sex.2<-ifelse(d$sex==2,1,0);d$school.GW<-ifelse(d$school=="Grant-White",1,0)
d[,c("id","sex","ageyr","agemo","school","grade")]<-NULL;d<-nmfkc::nmfkc.normalize(d)
ex<-c("age.rev","sex.2","school.GW");en<-setdiff(colnames(d),ex)
emit(ECM(t(d[,en]),t(d[,ex]),3,"free","seminmf","Factor"),"HS1939",300,9,3)
# OD-HUB  CURRENT paper spec: 3 hubs (Tokyo=13, Osaka=27, Fukuoka=40), Q=3, P=44
dd<-as.data.frame(read_xlsx("OD_estat_R3_000040170612.xlsx",sheet=1,skip=2))
dd<-dd[dd[,1]=="1"&dd[,3]=="02"&dd[,5]!="100"&dd[,7]!="100",];pref<-unique(dd[,6])
M<-matrix(NA,47,47);rownames(M)<-colnames(M)<-pref;dd[,5]<-as.numeric(dd[,5]);dd[,7]<-as.numeric(dd[,7])
for(i in 1:47)for(j in 1:47)M[i,j]<-dd[which(dd[,5]==i&dd[,7]==j),9];M<-log(1+M);hb<-c(13,27,40)
emit(ECM(M[-hb,],M[hb,],3,"nonneg","mu","Region"),"OD-HUB(3hub)",47,44,3)
# Topic (Q=3, nonneg) -- emits learned column order (paper relabels to content order separately)
corp<-corpus(data_corpus_inaugural);tk<-tokens_remove(tokens(corp),pattern=stopwords("en",source="snowball"))
df<-dfm_trim(dfm_select(dfm(tk),min_nchar=3),min_termfreq=100);dm<-as.matrix(df);dm<-dm[,order(colSums(dm),decreasing=TRUE)]
At<-rbind(pre1865=ifelse(corp$Year<=1865,1,0),y1865_45=ifelse(corp$Year>1865&corp$Year<=1945,1,0),after1945=ifelse(corp$Year>1945,1,0))
Yt<-t(dm);colnames(At)<-colnames(Yt)
emit(ECM(Yt,At,3,"nonneg","mu","Topic"),"TOPIC(learned order)",59,67,3)
cat("\nDONE\n")
