#' cell renal cell carcinomas (ccRCC) dataset
#' 
#' Dataset of clear cell renal cell carcinomas (ccRCC) (Gerlinger, et al. 2014) and non-small cell lung cancer (NSCLC) (Zhang et al., 2014), combined from the multiregional sequencing variant allele frequency. The dataset contain 18 patients.
#' 
#'
#' @name ccRCC
#' @docType data
#' @format A list of 18 VAF matrices
#' @section Variables:
#' Variables:
#' \itemize{
#' \item chrom The number of chromosome
#' \item pos Genomic position
#' \item DESC Description of the position such as the gene symbol.
#' \item normal The VAF of the normal cells.
#' \item Other columns Each column indicates the VAF of each region
#' }
#' @author Yusuke Matsui & Teppei Shimamura
#' @references Gerlinger,M. et al. (2014) Genomic architecture and evolution of clear cell renal cell carcinomas defined by multiregion sequencing. Nat. Genet., 46, 225–233.
#' @references Zhang,J. et al. (2014) Intra-tumor heterogeneity in localized lung adenocarcinomas delineated by multi-region sequencing. Science, 346, 256–259.
NULL
