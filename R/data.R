#' Spatial transcriptomics and metabolomics data of mouse brain
#'
#' Spatial transcriptomics and metabolomics data measured on the same tissue sections
#'  of mouse brains on a regular grid by \insertCite{Vicari2024;nobrackets}{sbivar}. Only a subset of the data, consisting of
#'  the 5 most abundant transcripts and metabolites for 6 samples, are included in the package for computational and memory reasons.
#'  The images were pre-aligned manually with the help of MAGPIE \insertCite{Williams2025}{sbivar}.
#' The data consist of two lists of outcome variables and their coordinates.
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
#' Spatial transcriptomics and proteomics data of data of human hepatocellular carcinoma (HCC)
#'
#' Single-molecule spatial transcriptomics and lattice protein immunofluorescence data measured on the same tissue section
#'  of a human hepatocellular carcinoma sample by \insertCite{Duchini2026;nobrackets}{sbivar}. Only a subset of the data, consisting of
#'  two transcripts and two proteins is included in the package for computational and memory reasons.
#'  The images were pre-aligned by the authors
#'
#' @format A point pattern, and two matrices
#' \describe{
#'   \item{TranscriptsDuchini}{A point pattern of class ppp}
#'   \item{ProteinCoordsDuchini}{Protein coordinate matrix}
#'   \item{ProteinOutcomesDuchini}{Protein outcome matrix}
#' }
#' @source \doi{10.64898/2026.08.17.742355}
#' @references
#' \insertAllCited{}
#' @usage data(Duchini)
"Duchini"
