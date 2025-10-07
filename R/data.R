#' Spatial transcriptomics and metabolomics data of mouse brain
#'
#' Spatial transcriptomics and metabolomics data measured on the same tissue sections
#'  of mouse brains on a regular grid by \insertCite{Vicari2024;nobrackets}{sbivar}. Only a subset of the data, consisting of
#'  the 10 most abundant transcripts and metabolites for 6 samples, are included in the package for computational and memory reasons.
#'  The images were pre-aligned manually with the help of MAGPIE \insertCite{Williams2025}{sbivar}.
#' The data consist of two lists of outcome variables and two lists of coordinates
#'
#' @format Four lists of data matrices:
#' \describe{
#'   \item{TranscriptCoords,MetaboliteCoords}{Coordinate lists}
#'   \item{TranscriptOutcomes,MetaboliteOutcomes}{Outcome matrices}
#' }
#' @source \doi{10.1038/s41587-023-01937-y}
#' @references
#' \insertAllCited{}
#' @usage data(Vicari)
"Vicari"
