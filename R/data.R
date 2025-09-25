#' Spatial transcriptomics and metabolomics data of mouse brain
#'
#' Spatial transcriptomics and metabolomics data measured on the same tissue sections
#'  of mouse brains on a regular grid by \insertCite{Vicari2024}{sbivar}. Only a subset of the data, consisting of
#'  the 20 most abundant transcripts and metabolites for 6 samples are included in the package for computational and memory reasons.
#' The data consist of two lists of outcome variables and two lists of coordinates
#'
#' @format Four lists of A data matrix
#' \describe{
#'   \item{TranscriptCoords,MetaboliteCoords}{Coordinate lists}
#'   \item{TranscriptOutcomes,MetaboliteOutcomes}{Character vector with gene identities}
#' }
#' @source \doi{10.1038/s41587-023-01937-y}
#' @references
#' \insertAllCited{}
#' @usage data(Vicari)
"Vicari"
