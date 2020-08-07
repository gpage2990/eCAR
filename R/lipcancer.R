#' Number of recorded lip cancer cases in the 56 districts of Scotland.
#'
#' A list containing the data frame and neighborhood matrix based on 56 districts of Scotland that is with the following variables
#'
#' @format data: A data frame with 56 rows and the following 6 variables:
#' \describe{
#'   \item{observed}{observed number of cancer cases}
#'   \item{expected}{the expected number of lip cancer cases computed using indirect standardisation based on Scotland-wide disease rates}
#'   \item{pcaff}{percentage of the district’s workforce employed in agriculture, fishing and forestry}
#'   \item{latitude}{latitude coordinates}
#'   \item{longitude}{longitude coordinates}
#'   \item{name}{name}
#' }
#'
#' @format neighborhood.Matrix: A 56 x 56 matrix neighbhorhood matrix
#'
"lipcancer"
