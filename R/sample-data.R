#' COVID-19 Sample data
#'
#' Data from a an experiment investigating T cell compositions between COVID-19
#' patients and healthy control. This data has been transformed using a arcsinh transform
#' using a co-factor of 5 and randomly subsetted
#'
#' @docType data
#'
#' @usage data(COVIDSampleData)
#'
#' @format An object of class \code{"SingeCellExperiment"}
#'
#' @keywords datasets
#'
#' @references De Biasi et al. (2020) Nat Commun 11, 3434
#' (\href{https://www.nature.com/articles/s41467-020-17292-4}{Nature})
#'
#' @source \href{https://flowrepository.org/id/FR-FCM-Z2N5}{FlowRepository}
#'
#' @import SingleCellExperiment
#' @examples
#' data(COVIDSampleData)
"DeBiasi_COVID_CD8_samp"
