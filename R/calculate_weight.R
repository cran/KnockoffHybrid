#' Calculate weight from population designs
#'
#' Calculate weight using population genotype or summary statistics.
#'
#' @param pval A numeric vector of length p for p-values. P-values must be between 0 and 1. If not NULL, weight will be calculated as -log10(p-value).
#' @param beta A numeric vector of length p for beta coefficients. If not NULL, weight will be calculated as the absolute value of beta coefficients.
#' @param method A character string for the name of the weight estimation method. Must not be NULL if population genotype is used to calculate weight. Weight can be calculated using "score" (i.e., single variant score test) or "lasso" (i.e., least absolute shrinkage and selection operator). The default is "score".
#' @param geno A n*p matrix for the population genotype data, in which n is the number of subjects and p is the number of variants. The genotypes must be coded as 0, 1, or 2.
#' @param y A numeric vector of length n for the phenotype data for the n subjects.
#' @param phetype A character for the variable type of the phenotype. The type can be "C" (i.e., continuous) or "D" (i.e., dichotomous). The default is "D".
#' @param PCs A n*k matrix for the principal components of population structure, in which n is the number of subjects and k is the number of (top) principal components. If not NULL, principal components will be included as covariates when calculating weight from population genotype.
#' @return A numeric vector of length p for the weight.
#' @import glmnet
#' @import SPAtest
#' @importFrom stats glm pchisq var
#' @export
#' @examples
#' data(KnockoffHybrid.example)
#' weight<-calculate_weight(geno=KnockoffHybrid.example$dat.pop,y=KnockoffHybrid.example$y.pop)

calculate_weight<-function(pval=NULL,beta=NULL,method="score",geno=NULL,y=NULL,phetype="D",PCs=NULL)
{
  if (!is.null(pval)) {
    if (any(pval<0) | any(pval>1))
      stop("The p-values must be between 0 and 1.")
    weight = -log10(pval)
    return(weight)
  }
  if (!is.null(beta)) {
    weight = abs(beta)
    return(weight)
  }
  if (is.null(method)) stop("Weight estimation method must be provided if p-values or beta coefficients are not provided.")
  if (!method %in% c("lasso","score")) stop("Please choose method from lasso or score.")
  if (is.null(geno) | is.null(y) | is.null(phetype)) stop("Genotype, phenotype, and type of phenotpye must be provided.")
  if (!phetype %in% c("C","D")) stop("Please choose phetype from C or D.")
  nsnp<-ncol(geno)
  y<-as.numeric(y)
  if (method=="lasso"){
    if (phetype=="D") family <-"binomial" else if (phetype=="C") family <- "gaussian"
    if (!is.null(PCs)) geno<-cbind(geno,PCs)
    cv.lasso.pca <- cv.glmnet(geno, y, alpha = 1, family = family, parallel = FALSE)
    l.min.pca<-cv.lasso.pca$lambda.min
    lasso.model.pca<-glmnet(geno, y, family = family, alpha=1, lambda = l.min.pca)
    beta.pca<-as.vector(lasso.model.pca$beta)
    beta.pca<-beta.pca[1:nsnp]
    weight<-abs(beta.pca)
    return(weight)
  }
  if (method=="score"){
    GeneScan.Null.Model<-function(Y, X=NULL, id=NULL, out_type="D", resampling=FALSE,B=1000){

      Y<-as.matrix(Y);n<-nrow(Y)

      if(length(X)!=0){X0<-svd(as.matrix(X))$u}else{X0<-NULL}
      X0<-cbind(rep(1,n),X0)

      if(out_type=="C"){nullglm<-glm(Y~0+X0,family="gaussian")}
      if(out_type=="D"){nullglm<-glm(Y~0+X0,family="binomial")}

      if (length(id)==0){id<-1:n}

      mu<-nullglm$fitted.values;Y.res<-Y-mu;

      re.Y.res=NULL
      if(resampling==TRUE){
        index<-sapply(1:B,function(x)sample(1:length(Y)));temp.Y.res<-Y.res[as.vector(index)]
        re.Y.res<-matrix(temp.Y.res,length(Y),B)
      }

      if(out_type=='D'){v<-mu*(1-mu)}else{v<-rep(as.numeric(var(Y.res)),length(Y))}
      inv.X0<-solve(t(X0)%*%(v*X0))

      result.null.model<-list(Y=Y,id=id,n=n,mu=mu,res=Y.res,v=v,
                              X0=X0,nullglm=nullglm,out_type=out_type,
                              re.Y.res=re.Y.res,inv.X0=inv.X0)
      return(result.null.model)
    }
    Get.p<-function(X,result.null.model){
      X<-as.matrix(X)
      mu<-result.null.model$nullglm$fitted.values;Y.res<-result.null.model$Y-mu
      outcome<-result.null.model$out_type
      if(outcome=='D'){
        p<-ScoreTest_SPA(t(X),result.null.model$Y,result.null.model$X,method=c("fastSPA"),minmac=-Inf)$p.value
      }else{
        v<-rep(as.numeric(var(Y.res)),nrow(X))
        p<-pchisq(as.vector((t(X)%*%Y.res)^2)/(apply(X*(v*X),2,sum)-apply(t(X)%*%(v*result.null.model$X0)%*%result.null.model$inv.X0*t(t(result.null.model$X0)%*%as.matrix(v*X)),1,sum)),df=1,lower.tail=F)
      }
      return(as.matrix(p))
    }
    result.null.model.pca=GeneScan.Null.Model(y, X = PCs, id = NULL, out_type=phetype, resampling=FALSE)
    p.score.pca<-as.vector(Get.p(geno,result.null.model.pca))
    weight<- -log10(p.score.pca)
    return(weight)
  }
}
