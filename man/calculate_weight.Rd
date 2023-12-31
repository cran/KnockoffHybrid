% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_weight.R
\name{calculate_weight}
\alias{calculate_weight}
\title{Calculate weight from population designs}
\usage{
calculate_weight(
  pval = NULL,
  beta = NULL,
  method = "score",
  geno = NULL,
  y = NULL,
  phetype = "D",
  PCs = NULL
)
}
\arguments{
\item{pval}{A numeric vector of length p for p-values. P-values must be between 0 and 1. If not NULL, weight will be calculated as -log10(p-value).}

\item{beta}{A numeric vector of length p for beta coefficients. If not NULL, weight will be calculated as the absolute value of beta coefficients.}

\item{method}{A character string for the name of the weight estimation method. Must not be NULL if population genotype is used to calculate weight. Weight can be calculated using "score" (i.e., single variant score test) or "lasso" (i.e., least absolute shrinkage and selection operator). The default is "score".}

\item{geno}{A n*p matrix for the population genotype data, in which n is the number of subjects and p is the number of variants. The genotypes must be coded as 0, 1, or 2.}

\item{y}{A numeric vector of length n for the phenotype data for the n subjects.}

\item{phetype}{A character for the variable type of the phenotype. The type can be "C" (i.e., continuous) or "D" (i.e., dichotomous). The default is "D".}

\item{PCs}{A n*k matrix for the principal components of population structure, in which n is the number of subjects and k is the number of (top) principal components. If not NULL, principal components will be included as covariates when calculating weight from population genotype.}
}
\value{
A numeric vector of length p for the weight.
}
\description{
Calculate weight using population genotype or summary statistics.
}
\examples{
data(KnockoffHybrid.example)
weight<-calculate_weight(geno=KnockoffHybrid.example$dat.pop,y=KnockoffHybrid.example$y.pop)
}
