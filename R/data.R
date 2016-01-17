#' Variant allele frequency data (VAF) set
#' 
#' This data contains the multiregional sequencing variant allele frequency data set from https://github.com/viq854/lichee/tree/master/LICHeE/data.
#' The data set contain 14 patients VAFs consisting the two studies.
#' 
#'
#' @name vaf
#' @docType data
#' @format A list of 14 VAF matrices
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
NULL
