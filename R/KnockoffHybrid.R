#' Calculate KnockoffHybrid's feature statistics
#'
#' Calculate KnockoffHybrid's feature statistics using original and knockoff genotype data.
#'
#' @param dat A 3n*p matrix for the original trio genotype data, in which n is the number of trios and p is the number of variants. Each trio must consist of father, mother, and offspring (in this order). The genotypes must be coded as 0, 1, or 2. Missing genotypes are not allowed.
#' @param dat.ko A 3n*p*M array for the knockoff trio genotype data created by function create_knockoff. M is the number of knockoffs.
#' @param pos A numeric vector of length p for the position of p variants.
#' @param allele A vector of length p for the minor allele at each position. Minor alleles for windows with multiple variants will be shown as "W" in the output.
#' @param start An integer for the first position of sliding windows. If NA, start=min(pos). Only used if you would like to use the same starting position for different cohorts/analyses.
#' @param end An integer for the last position of sliding windows. If NA, end=max(pos). Only used if you would like to use the same ending position for different cohorts/analyses.
#' @param size A numeric vector for the size(s) of sliding windows when scanning the genome
#' @param p_value_only A logical value indicating whether to perform the knockoff analysis. When p_value_only is TRUE, only the ACAT-combined p-values are to be calculated for each window. When p_value_only is FALSE, dat.ko is required and KnockoffHybrid's feature statistics are to be calculated for each window in addition to the p-values.
#' @param adjust_for_cov A logical value indicating whether to adjust for covariates. When adjust_for_cov is TRUE, y is required.
#' @param y A numeric vector of length n for the residual Y-Y_hat. Y_hat is the predicted value from the regression model in which the quantitative trait Y is regressed on the covariates. If Y is dichotomous, you may treat Y as quantitative when applying the regression model.
#' @param chr A character for the name of the chromosome, e.g., "1", "2", ..., "22", and "X".
#' @param sex A numeric vector of length n for the sex of offspring. 0s indicate females and 1s indicate males.
#' @param weight A numeric vector of length p for the weight of p variants. The weight can be obtained via the function "calculate_weight" using population genotype or summary statistics. If NULL, the weight will be calculated based on minor allele frequencies.
#' @return A data frame for the hybrid analysis results. Each row contains the p-values and, if p_value_only is FALSE, KnockoffHybrid's feature statistics for a window.
#' @importFrom stats as.dist cutree hclust median pcauchy pnorm princomp
#' @export
#' @examples
#' data(KnockoffHybrid.example)
#' dat.ko<-create_knockoff(KnockoffHybrid.example$dat.hap,KnockoffHybrid.example$pos,M=10)
#' weight<-calculate_weight(geno=KnockoffHybrid.example$dat.pop,y=KnockoffHybrid.example$y.pop)
#' window<-KnockoffHybrid(dat=KnockoffHybrid.example$dat,dat.ko=dat.ko,
#'         pos=KnockoffHybrid.example$pos,weight=weight)
KnockoffHybrid<-function (dat, dat.ko = NA, pos, allele = NA, start = NA, end = NA,
                          size = c(1,1000, 5000, 10000, 20000, 50000),
                          p_value_only = FALSE, adjust_for_cov = FALSE,
                          y = NA, chr = "1", sex = NA, weight=NULL)
{
  if (nrow(dat)%%3 != 0)
    stop("The number of rows of the original trio matrix must be a multiple of three.")
  if (is.na(start))
    start <- min(pos)
  if (is.na(end))
    end <- max(pos)
  n <- nrow(dat)/3
  nsnp <- ncol(dat)
  out_fbat <- fbat(dat, adjust_for_cov, y, sex = sex)
  z <- out_fbat$additive
  p1 <- 2 * pnorm(-abs(z))
  p1_sign <- sign(z)
  index_dad <- seq(1, (3 * n), 3)
  index_mom <- seq(2, (3 * n), 3)
  index_off <- seq(3, (3 * n), 3)
  maf <- apply(dat[-index_off,,drop=F], 2, mean)/2

  if (is.null(weight)){ #reduces to KnockoffTrio
    weight <- 1/(sqrt(n * maf * (1 - maf)))
    weight[which(weight == Inf)] <- 0
  } else{
    if (length(weight) != nsnp)
      stop("The length of the weight vector should be equal to the number of SNPs.")
  }
  Z <- sweep(out_fbat$Z, 2, weight, "*")
  Zsq <- t(Z) %*% Z
  Zsum <- apply(Z, 2, sum)
  window <- data.frame()
  index <- which(maf >= 0)
  for (i in 1:length(size)) {
    if (size[i] == 1) {
      if (length(index) > 0)
        window <- data.frame(start = pos[index], end = pos[index],
                             ind3 = index, ind4 = index, n = rep(1, length(index)))
    }
    else window <- rbind(window, getslidingwindow(start,
                                                  end, size[i], pos))
  }
  window <- window[which(window$n > 0), ]
  window <- window[!duplicated(window[, c(3, 4)]), ]
  tmp <- which(window$n == 1 & !(window$ind3 %in% index))
  if (length(tmp) > 0)
    window <- window[-tmp, ]
  nwindow <- nrow(window)

  #add allele and frq information
  alle<-rep("W",nwindow)
  frq<-rep(0,nwindow)
  for (i in 1:nwindow){
    if (window$n[i]==1){
      alle[i]<-allele[window$ind3[i]]
      frq[i]<-maf[window$ind3[i]]
    }else frq[i]<-mean(maf[window$ind3[i]:window$ind4[i]])
  }
  window<-cbind(window,allele=alle,frq=frq)

  q1 <- p.fbat <- z1 <- rep(NA, nwindow)
  direction <- rep(NA, nwindow)
  if (!p_value_only) {
    M <- dim(dat.ko)[3]
    out_fbat <- fbat(dat, adjust_for_cov, y, dosage = TRUE,
                     dat.ko, sex)
    z.ko <- out_fbat$additive
    Z2 <- out_fbat$Z
    p2 <- 2 * pnorm(-abs(z.ko))
    Z2sq <- array(dim = c(nsnp, nsnp, M))
    Z2sum <- array(dim = c(nsnp, M))
    if (dim(Z2)[2] > 1) {
      for (m in 1:M) {
        Z2[, , m] <- sweep(Z2[, , m], 2, weight, "*")
        Z2sq[, , m] <- t(Z2[, , m]) %*% Z2[, , m]
        Z2sum[, m] <- apply(Z2[, , m], 2, sum)
      }
    }
    else if (dim(Z2)[2] == 1) {
      for (m in 1:M) {
        Z2[, , m] <- Z2[, , m] * weight
        Z2sq[, , m] <- t(Z2[, , m]) %*% Z2[, , m]
        Z2sum[, m] <- sum(Z2[, , m])
      }
    }
    q2 <- z2 <- array(dim = c(M, nwindow))
  }
  for (i in 1:nwindow) {
    ind <- window$ind3[i]:window$ind4[i]
    if (length(ind) == 1) {
      if (maf[ind] >= 0.01) {
        q1[i] <- p.fbat[i] <- p1[ind]
        z1[i] <- z[ind]
        direction[i] <- p1_sign[ind]
      }
    }
    else {
      z1[i] <- fbat_set(Zsum, Zsq, ind)
      p.fbat[i] <- 2 * pnorm(-abs(z1[i]))
      ind.single <- ind
      p.single <- p1[ind.single]
      q1[i] <- ACAT(c(p.single, p.fbat[i]))
      stat <- c(z1[i], z[ind.single])
      tmp <- which.max(abs(stat))
      if (length(tmp) > 0)
        direction[i] <- sign(stat[tmp])
    }
    if (!p_value_only) {
      if (length(ind) == 1) {
        if (maf[ind] >= 0.01) {
          q2[, i] <- p2[ind, ]
          z2[, i] <- z.ko[ind, ]
        }
      }
      else {
        for (m in 1:M) {
          p.single.ko <- p2[ind.single, m]
          z2[m, i] <- fbat_set(Z2sum[, m], Z2sq[, , m],
                               ind)
          p.burden.ko <- 2 * pnorm(-abs(z2[m, i]))
          q2[m, i] <- ACAT(c(p.single.ko, p.burden.ko))
        }
      }
    }
  }
  window$ind3 <- pos[window$ind3]
  window$ind4 <- pos[window$ind4]
  colnames(window)[3] <- "actual_start"
  colnames(window)[4] <- "actual_end"
  if (!p_value_only) {
    out <- calculate_w_kappatau(q1, q2)
    w <- out$w
    kappatau <- out$kappatau
    rownames(q2) <- paste0("p_", 1:M)
    q2 <- t(q2)
    rownames(z2) <- paste0("z_", 1:M)
    z2 <- t(z2)
    window <- cbind(chr, window, dir = direction, w, p = q1,
                    z = z1, p.burden = p.fbat, kappatau, q2, z2)
  }
  else window <- cbind(chr, window, dir = direction, p = q1,
                       z = z1, p.burden = p.fbat)
  return(window)
}
fbat<-function(dat,adjust_for_covariates=FALSE,y=NA,dosage=FALSE,dat1=NA,sex=NA){
  if (is.null(ncol(dat))) dat<-as.matrix(dat)
  n<-nrow(dat)/3
  index_dad<-seq(1,(3*n),3) #index for parent1
  index_mom<-seq(2,(3*n),3) #index for parent2
  index_off<-seq(3,(3*n),3) #index for offspring
  if (dosage==FALSE) Z<-dat[index_off,]-(dat[index_dad,]+dat[index_mom,])/2 else Z<-dat1[index_off,,,drop=FALSE]-(dat1[index_dad,,,drop=FALSE]+dat1[index_mom,,,drop=FALSE])/2
  if (adjust_for_covariates){
    Z<-sweep(Z,1,y,'*')
  }
  additive<-colSums(Z)/sqrt(colSums(Z^2)) #if dosage=T, dim=c(nsnp,M)
  additive[is.na(additive)]<-0
  out<-list()
  out$additive<-additive
  out$Z<-Z
  return (out)
}
fbat_set<-function(W,V,ind){
  if (length(ind)>0){
    s1<-sum(W[ind])
    s2<-sum(V[ind,ind])
    if (s2==0) s2<-1
    return (unname(s1/sqrt(s2)))
  } else return(NA)
}
calculate_w_kappatau<-function(q1,q2){
  out<-list()
  t1<--log10(q1)
  t2<--log10(q2)
  t2_med<-apply(t2,2,median)
  t2_max<-apply(t2,2,max)
  out$w<-(t1-t2_med)*(t1>=t2_max)
  out$w.raw<-t1-t2_med
  out$kappatau<-MK.statistic(t1,t(t2),method="median")
  #out$q<-MK.q.byStat(out$kappatau[,1],out$kappatau[,2],M=M)
  return(out)
}
MK.statistic<-function (T_0,T_k,method='median'){
  T_0<-as.matrix(T_0);T_k<-as.matrix(T_k)
  T.temp<-cbind(T_0,T_k)
  T.temp[is.na(T.temp)]<-0
  kappa<-apply(T.temp,1,which.max)-1
  if(method=='median'){
    Get.OtherMedian<-function(x){median(x[-which.max(x)])}
    tau<-apply(T.temp,1,max)-apply(T.temp,1,Get.OtherMedian)
  }
  return(cbind(kappa,tau))
}
getslidingwindow<-function(left,right,size,pos){
  if (size<=0) stop("The size of sliding window must be a positive number.")
  start<-end<-n<-ind3<-ind4<-c()
  step<-ceiling(size/2) #move half of the window size each time
  l<-left-step
  r<-l+size-1
  while (r<=(right+step)){
    start<-c(start,l)
    end<-c(end,r)
    tmp<-which(pos>=l & pos<=r)
    n<-c(n,length(tmp))
    if (length(tmp)==0) tmp<-0
    ind3<-c(ind3,tmp[1])
    ind4<-c(ind4,tmp[length(tmp)])
    l<-l+step
    r<-r+step
  }
  return(data.frame(start,end,ind3,ind4,n))
}
ACAT<-function(p){
  p[p>0.99]<-0.99 #p[p>1-1e-2]<-1-1e-2
  is.small<-(p<1e-16) & !is.na(p)
  is.regular<-(p>=1e-16) & !is.na(p)
  temp<-rep(NA,length(p))
  temp[is.small]<-1/p[is.small]/pi
  temp[is.regular]<-as.numeric(tan((0.5-p[is.regular])*pi))

  cct.stat<-mean(temp,na.rm=T)
  if(is.na(cct.stat)){return(NA)}
  if(cct.stat>1e+15){return((1/cct.stat)/pi)}else{
    return(1-pcauchy(cct.stat))
  }
}
