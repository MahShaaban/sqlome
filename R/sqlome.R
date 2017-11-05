#' \code{sqlome} package
#'
#' Build SQLite tables of microRNAs and Transcription Factors-gene Correlations
#'
#' @import RTCGA.miRNASeq
#' @import RTCGA.rnaseq
#' @import RTCGA.RPPA
#' @import targetscan.Hs.eg.db
#' @import org.Hs.eg.db
#'
#' @docType package
#' @name sqlome
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
## fix by @jennybc
## source https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
if(getRversion() >= "2.15.1")  utils::globalVariables(c('BCR',
                                                        'Cohort',
                                                        'RPPA',
                                                        'ab_id',
                                                        'feature',
                                                        'mRNASeq',
                                                        'miRNA_ID',
                                                        'miRSeq',
                                                        'mirna_base',
                                                        'name',
                                                        'value',
                                                        'x'))
