#' Palettes for the methods of \code{GeRnika}
#'
#' A data.frame containing 3 default palettes for the parameters used in the methods of \code{GeRnika}.
#'
#' @format A data.frame containing 3 palettes.
#' \describe{
#'   \item{Lancet}{#0099B444, #AD002A77, #42B540FF}
#'   \item{NEJM}{#FFDC9177, #7876B188, #EE4C97FF}
#'   \item{Simpsons}{#FED43966, #FD744688, #197EC0FF}
#' }
#' @source Lancet, NEJM and The Simpsons palettes; inspired by the plots in Lancet journals, the plots in the New England Journal of Medicine and the colors used in the TV show The Simpsons, respectively.
"palettes"

#' A set of 10 trios of B matrices for experimenting with the methods of \code{GeRnika}
#'
#' A list of lists composed by 10 trios of B matrices; a real B matrix, a B matrix got by using the one algorithm (alg1) method and another one as a result of another algorithm (alg2). These matrices can be used as examples for the methods of \code{GeRnika}.
#'
#' @format A list of lists composed by 10 trios of B matrices.
#' \describe{
#'   \item{Trio 1}{B_real, B_alg1 and B_alg2 (matrices composed by 5 clones)}
#'   \item{Trio 2}{B_real, B_alg1 and B_alg2 (matrices composed by 5 clones)}
#'   \item{Trio 3}{B_real, B_alg1 and B_alg2 (matrices composed by 5 clones)}
#'   \item{Trio 4}{B_real, B_alg1 and B_alg2 (matrices composed by 5 clones)}
#'   \item{Trio 5}{B_real, B_alg1 and B_alg2 (matrices composed by 5 clones)}
#'   \item{Trio 6}{B_real, B_alg1 and B_alg2 (matrices composed by 10 clones)}
#'   \item{Trio 7}{B_real, B_alg1 and B_alg2 (matrices composed by 10 clones)}
#'   \item{Trio 8}{B_real, B_alg1 and B_alg2 (matrices composed by 10 clones)}
#'   \item{Trio 9}{B_real, B_alg1 and B_alg2 (matrices composed by 10 clones)}
#'   \item{Trio 10}{B_real, B_alg1 and B_alg2 (matrices composed by 10 clones)}
#' }
#' @source Local source; as a result of the Grasp and the ILS methods used for solving the Clonal Deconvolution and Evolution Problem (CDEP).
"B_mats"