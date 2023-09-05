#' Identification of causal loci
#'
#' Identification of causal loci using KnockoffHybrid's feature statistics
#' @param window The result window from KnockoffHybrid. If there are multiple windows, please use rbind to combine the windows.
#' @param M A positive integer for the number of knockoffs. The default is 10.
#' @param fdr A real number in a range of (0,1) indicating the target FDR level. The default is 0.2.
#' @return A list that contains:
#' \describe{
#'   \item{window}{A data frame for an updated window that includes an extra column for KnockoffHybrid's Q-values. A locus with a Q-value <= the target FDR level, i.e., window$q<=fdr, is considered as causal.}
#'   \item{thr.w}{A positive real number indicating the significance threshold for KnockoffTrio's feature statistics. A locus with a feature statistic >= thr.w, i.e., window$w>=thr.w is considered as causal. The loci selected by window$w>=thr.w are equivalent to those by window$q<=fdr. No loci are selected at the target FDR level if thr.w=Inf.}
#' }
#' @importFrom stats as.dist cutree hclust median pcauchy pnorm princomp
#' @export
#' @examples
#' data(KnockoffHybrid.example)
#' dat.ko<-create_knockoff(KnockoffHybrid.example$dat.hap,KnockoffHybrid.example$pos,M=10)
#' weight<-calculate_weight(geno=KnockoffHybrid.example$dat.pop,y=KnockoffHybrid.example$y.pop)
#' window<-KnockoffHybrid(dat=KnockoffHybrid.example$dat,dat.ko=dat.ko,
#'         pos=KnockoffHybrid.example$pos,weight=weight)
#' result<-causal_loci(window,M=10,fdr=0.2)
causal_loci<-function(window,M=10,fdr=0.2){
  q<-MK.q.byStat(window$kappa,window$tau,M=M)
  window<-cbind(window,q)
  thr.w<-MK.threshold.byStat(window$kappa,window$tau,M=M,fdr=fdr)
  out<-list()
  out$window<-window
  out$thr.w<-thr.w
  return(out)
}
MK.q.byStat<-function (kappa,tau,M,Rej.Bound=10000){
  b<-order(tau,decreasing=T)
  c_0<-kappa[b]==0
  #calculate ratios for top Rej.Bound tau values
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    #if(i==1){temp_0=c_0[i]}
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  #calculate q values for top Rej.Bound values
  q<-rep(1,length(tau))
  for(i in 1:length(b)){
    q[b[i]]<-min(ratio[i:min(length(b),Rej.Bound)])*c_0[i]+1-c_0[i]
    if(i>Rej.Bound){break}
  }
  return(q)
}
MK.threshold.byStat<-function (kappa,tau,M,fdr = 0.2,Rej.Bound=10000){
  b<-order(tau,decreasing=T)
  c_0<-kappa[b]==0
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  ok<-which(ratio<=fdr)
  if(length(ok)>0){
    return(tau[b][ok[length(ok)]])
  }else{return(Inf)}
}

