#' Simulation dataset
#' 
#' Simulation dataset based on BEP simulator (Niida, et al. 2015). The simulation mimics the cancer sub-clonal evolution using cellular automaton. We generated two types of cancer, that is, high and low heterogeneity.
#'
#' @name simulation
#' @docType data
#' @format A list of 60 VAF matrices
#' @section Variables:
#' Variables:
#' \itemize{
#' \item normal The VAF of the normal cells.
#' \item Other columns Each column indicates the VAF of each region
#' }
#' @author Yusuke Matsui & Teppei Shimamura
#' @references Niida,A. et al. (2015) Cancer evolution simulation identifies possible principles underlying intratumor heterogeneity. bioRxiv, 022806.
NULL
