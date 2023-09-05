#' Example data for KnockoffHybrid
#'
#' A toy example of the haplotype and genotype data for trio designs and the genotype and phenotype data for population designs
#'
#' @format KnockoffHybrid.example contains the following items:
#' \describe{
#'   \item{dat}{A numeric genotype matrix of 3 trios and 5 variants. Each trio contains 3 rows in the order of father, mother and offspring. Each column represents a variant.}
#'   \item{dat.hap}{A numeric haplotype matrix of 3 trios and 5 variants. Each trio contains 6 rows in the order of father, mother and offspring. Each column represents a variant.}
#'   \item{dat.pop}{A numeric genotype matrix of 8 subjects and 5 variants. Each row represents a subject and each column represents a variant.}
#'   \item{y.pop}{A numeric vector of length 8 for the dichotomous phenotypes of 8 subjects.}
#'   \item{pos}{A numeric vector of length 5 for the position of 5 variants.}
#' }
"KnockoffHybrid.example"
