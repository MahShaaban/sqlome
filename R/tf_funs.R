#' Construct a url for Cistrome Cancer file
#'
#' This function constructs a url to downlaod the tarnscription factors correlation data from Cistrome Cancer. Briefly, the function constructs a url
#' for the file using a predefined handle and a suffix.
#'
#' @param tf A character string of the official symbol of the transcription factor
#' @param all A logicall of whether to get the correlation files of all genes (default) or only targets.
#'
#' @return A character string of the cistrome cancer url for a transcription factor
#' url <- tf_url('AFF4')
#'
#' @export
tf_url <- function(tf, all = TRUE) {
    # construct a suffix depending of type of file
    # all refers to files containing correlations of all genes not only targets
    if(all == TRUE) {
        suf <- '.all.cor.csv'
    } else {
        suf <- '.cor.csv'
    }

    # construct url by adding cistrome handle and file suffix
    tf_url <- paste('http://cistrome.org/CistromeCancer/CancerTarget/examples/', tf, suf, sep = '')

    return(tf_url)
}


#' Construct a file name for a Cistrome Cancer file
#'
#' @param tf A character string of the official symbol of the transcription factor
#' @param dir A character string for a directory to download files.
#' @param all A logicall of whether to get the correlation files of all genes (default) or only targets.
#'
#' @return A character string of the cistrome cancer file name for a transcription factor
#'
#' @examples
#' tf_fname('AFF4', dir = 'tmp/tf')
#'
#' @export
tf_fname <- function(tf, dir, all = TRUE) {
    # construct a suffix depending of type of file
    # all refers to files containing correlations of all genes not only targets
    if(all == TRUE) {
        suf <- '.all.cor.csv'
    } else {
        suf <- '.cor.csv'
    }

    # construct a file name
    tf_fname <- paste(dir, tf, suf, sep = '')

    return(tf_fname)
}


#' Read transcription factor correlation data
#'
#' This is a custom function to read the local csv files obtained from Cistrome Cancer into a list of tidy data.frames.
#' First, the function constructs a vector of available files in the given directory and reads each of them to a data.frame.
#' Then, it drops the damaged files and names the items of the list with the corresponding transcription factor name.
#'
#' @param dir A character string of the directory containing the Cistrome Cancer files
#'
#' @return A named list of tidy data.frames. Each item of the list correspond to a transcription factor. The data.frame consist of
#'   a column containing the names of the features and a single column for each of the TCGA studies containing the correspoinding
#'   correlation value.
#'
#' @importFrom readr read_csv
#' @importFrom stringr str_split
#'
#' @export
tf_read <- function(dir = 'tmp/tf/') {
    # get file names and reading intro a list
    fls <- list.files(dir,
                      pattern = '*.csv',
                      full.names = TRUE)
    tf <- lapply(fls, read_csv(x))

    # subset damaged files
    ind <- unlist(lapply(tf, nrow)) > 1
    tf <- tf[ind]

    # get names of okay transcription factors' files
    nms <- str_split(list.files(dir,
                                pattern = '*.csv'),
                     pattern = '/|\\.',
                     simplify = TRUE)[, 1]
    nms <- nms[ind]
    names(tf) <- nms

    return(tf)
}
