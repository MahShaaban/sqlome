#' Extract TCGA info
#'
#' Extract TCGA info for mRNASeq, miRSeq and RPPA. The function calls the RTCGA::infoTCGA and
#' removes entries that doesn't have all three assays performed; mRNASeq, miRSeq and RPPA.
#'
#' @return A data.frame with Cohort, BCR and a columns with the number of cases in each assay
#'
#' @examples
#' sqlome_info
#'
#' @importFrom RTCGA infoTCGA
#' @importFrom magrittr %>%
#' @importFrom dplyr select filter
#'
#' @export
sqlome_info <- function() {
    infoTCGA() %>%
        select(Cohort, BCR, mRNASeq, miRSeq, RPPA) %>%
        filter(mRNASeq != 0, miRSeq != 0, RPPA != 0)
}

#' Tidy RTCGA.rnaseq data.frame
#'
#' This is a custom function to tidy up the RTCGA.rnaseq data.frames. The aim is to strip out unnessecary data
#' and get to a typical bioconductor expression matrix. To achieve that, the function drops the first column,
#' extract the patient bcr and gene symbols and add them to the transposed matrix as column and row names.
#'
#' @param mrna A RTCGA.rnaseq data.frame
#'
#' @return A tidy expression matrix with genes in the rows and patient samples in the columns.
#'
#' @examples
#' library(RTCGA.rnaseq)
#' rtcga.df <- ACC.rnaseq
#' ACC.mat <- mrna_tidy(rtcga.df)
#'
#' @importFrom magrittr %>%
#' @importFrom stringr str_split
#'
#' @export
mrna_tidy <- function(mrna) {
    # drop the first column of the data frame and transpose
    mrna_mat <- t(mrna[, -1])

    # extract bcr
    mrna_bcr <- mrna$bcr_patient_barcode
    mrna_bcr <- str_split(mrna_bcr,
                          pattern = '\\.',
                          simplify = TRUE)[ , 1] %>%
        unique

    # extract gene ids
    mrna_gene <- colnames(mrna)[-1]
    mrna_gene <- str_split(mrna_gene,
                           pattern = '\\|',
                           simplify = TRUE)[, 1]

    # confirm dimensions match the bcr and gene ids
    if(dim(mrna_mat)[2] != length(mrna_bcr)) {
        stop("BCR are not unique.")
    }

    if(dim(mrna_mat)[1] != length(mrna_gene)) {
        stop("Gene IDs are not unique.")
    }

    # make columns names with bcr
    colnames(mrna_mat) <- mrna_bcr

    # make rownames with gene ids
    ind <- mrna_gene != '?' & !duplicated(mrna_gene)
    mrna_mat <- mrna_mat[ind, ]
    rownames(mrna_mat) <- mrna_gene[ind]

    # return tidy matrix
    return(mrna_mat)
}

#' Tidy RTCGA.miRNASeq data.frame
#'
#' This is a custom function to tidy up the RTCGA.miRNASeq data.frames. The aim is to strip out unnessecary data
#' and get to a typical bioconductor expression matrix. To achieve that, the function drops unneeded columns and rows,
#' extract the patient bcr and microRNA miRBase IDs and add them to the transposed matrix as column and row names.
#'
#' @param mirna A RTCGA.miRNASeq data.frame
#'
#' @return A tidy expression matrix with microRNAs in the rows and patient samples in the columns.
#'
#' @examples
#' library(RTCGA.miRNASeq)
#' rtcga.df <- ACC.miRNASeq
#' ACC.mat <- mirna_tidy(rtcga.df)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select
#' @importFrom stringr str_split
#'
#' @export
mirna_tidy <- function(mirna) {
    # drop unneeded rows from data frame and transpose
    mirna_mat <- mirna %>%
        filter(miRNA_ID == 'read_count') %>%
        select(-1, -2) %>%
        t

    # convert to numeric
    mirna_mat <- apply(mirna_mat, 2, as.numeric)

    # extract bcr
    mirna_bcr <- rownames(mirna)
    mirna_bcr <- str_split(mirna_bcr,
                           pattern = '\\.',
                           simplify = TRUE)[, 1] %>%
        unique

    # extract miRBase ids
    mirna_mibase <- colnames(mirna)[-1:-2]

    # confirm dimensions match the bcr and miRBase ids
    if(dim(mirna_mat)[2] != length(mirna_bcr))
        stop("BCR are not unique.")

    if(dim(mirna_mat)[1] != length(mirna_mibase))
        stop("miRBase IDs are not unique.")

    # make col and row names
    colnames(mirna_mat) <- mirna_bcr
    rownames(mirna_mat) <- mirna_mibase

    # return tidy matrix
    return(mirna_mat)
}

#' Tidy RTCGA.RPPA data.frame
#'
#' This is a custom function to tidy up the RTCGA.RPPA data.frames. The aim is to strip out unnessecary data
#' and get to a typical bioconductor expression matrix. To achieve that, the function drops the first column,
#' extract the patient bcr and antibody ID and add them to the transposed matrix as column and row names.
#'
#' @param rppa A RTCGA.RPPA data.frame
#'
#' @return A tidy expression matrix with proteins in the rows and patient samples in the columns.
#'
#' @examples
#' library(RTCGA.RPPA)
#' rtcga.df <- ACC.RPPA
#' ACC.mat <- rppa_tidy(rtcga.df)
#'
#' @importFrom magrittr %>%
#' @importFrom stringr str_split
#'
#' @export
rppa_tidy <- function(rppa) {
    # drop unneeded first column of data fram
    rppa_mat <- t(rppa[, -1])

    # extract bcr
    rppa_bcr <- rppa$bcr_patient_barcode
    rppa_bcr <- str_split(rppa_bcr,
                         pattern = '\\.',
                         simplify = TRUE)[ , 1] %>%
        unique

    # confirm dimensions match the bcr and antibody ids
    if(dim(rppa_mat)[2] != length(rppa_bcr))
        stop("BCR are not unique.")

    # make colnames
    colnames(rppa_mat) <- rppa_bcr

    # return tidy matrix
    return(rppa_mat)
}

#' Make mircoRNA feature correlations
#'
#' Calculate the pearson's correlation between microRNA and feature (gene/protein) matrices.
#' The function calls the R base cor function. More importantly the function prepare the input
#' matrices before the call. First, if matches the col.names/bcr. Second, remove duplicated entries.
#' Finally, it removes the microRNAs and features with expression values less than one. In addation,
#' The function returns a tidy data.frame rather than the raw matrix output of cor.
#'
#' @param mi A tidy matrix of microRNA expression data.
#' @param m A tidy matrix of gene or protein expression data.
#' @param cohort A charachter string of the TCGA ID.
#' @param tidy A logical
#'
#' @return A tidy data.frame of four columns; miRBase ID, feature id, correlation and TCGA study name.
#'
#' @examples
#' # load required libraries
#' library(RTCGA.rnaseq)
#' library(RTCGA.miRNASeq)
#' library(RTCGA.RPPA)
#'
#' # calculate correlation for genes
#' ACC.mi <- mirna_tidy(ACC.miRNASeq)
#' ACC.m <- mrna_tidy(ACC.rnaseq)
#' corr <- cor_make(ACC.mi, ACC.m, 'ACC')
#'
#' # calculate correlation for proteins
#' ACC.mi <- mirna_tidy(ACC.miRNASeq)
#' ACC.r <- rppa_tidy(ACC.RPPA)
#' corr <- cor_make(ACC.mi, ACC.r, 'ACC')
#'
#' @importFrom magrittr %>%
#' @importFrom stringr str_split
#' @importFrom dplyr filter mutate
#' @importFrom reshape2 melt
#' @importFrom stats cor setNames
#'
#' @export
cor_make <- function(mi, m, cohort, tidy = TRUE) {
    # get patient bcr for microRNA and remove duplicates
    bcr_mi <- str_split(colnames(mi),
                        pattern = '\\-',
                        simplify = TRUE)[, 3]
    mi <- mi[, !duplicated(bcr_mi)]
    colnames(mi) <- unique(bcr_mi)
    mi <- as.matrix(mi)

    # get patient bcr for feature matrix and remove duplicates
    bcr_m <- str_split(colnames(m),
                       pattern = '\\-',
                       simplify = TRUE)[, 3]

    m <- m[, !duplicated(bcr_m)]
    colnames(m) <- unique(bcr_m)

    # get intersection of bcr and subset matrices
    ind <- intersect(colnames(m), colnames(mi))
    mi <- mi[rowSums(mi) > 0, ind]
    m <- m[rowSums(m) > 0, ind]

    # calculate correlation
    c <- cor(t(m), t(mi))

    if(tidy == TRUE) {
        # tidy correlation
        nms <- c('feature', 'mirna_base', paste(cohort))
        c <- melt(c) %>%
            mutate(value = round(value, 2) * 100) %>%
            filter(abs(value) >= 10) %>%
            setNames(nms)
    }

    # return tidy correlation data.frame
    return(c)
}

#' Get microRNA targets
#'
#' This function uses the targetscan database annotations to get the microRNA gene and protien tragets. For microRNA-gene targets
#' it is a straight forward subsetting the database and converting to ta tidy data.frame with the appropriate IDs. For microRNA-protein
#' targets, it first obtain the gene targets as explained and then change transform the IDs to the corresponding protein/antibody IDs
#' provided by the manufacturers of the used assays.
#'
#' @param mirna A charachter vector of the names of microRNAs of interes
#' @param feature_type A character string to choose the desired type of annotations. Takes one of 'gene' or 'protein'.
#'
#' @return A tidy data.frame of microRNA in one column 'mirna_base' and targets column 'feature'. When used with the 'protein' another
#'   column is added to distinguish the feature_types genes or proteins.
#'
#' @examples
#' mirna_base <- c("hsa-let-7a", "hsa-let-7b", "hsa-let-7c")
#' gene <- get_targets(mirna_base, 'gene')
#' ptn <- get_targets(mirna_base, 'protein')
#'
#' @import targetscan.Hs.eg.db
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi toTable
#' @importFrom dplyr left_join select mutate
#' @importFrom readr read_csv
#' @importFrom stats na.omit setNames
#'
#' @export
get_targets <- function(mirna, feature_type = 'gene') {
    gene_targets <- toTable(targetscan.Hs.egMIRBASE2FAMILY) %>%
        left_join(toTable(targetscan.Hs.egTARGETS)) %>%
        select(-name) %>%
        setNames(c('mirna_base', 'feature')) %>%
        mutate(mirna_base = tolower(mirna_base)) %>%
        filter(mirna_base %in% mirna) %>%
        mutate(feature = AnnotationDbi::select(org.Hs.eg.db,
                                               feature,
                                               'SYMBOL',
                                               'ENTREZID')$SYMBOL)
    switch(feature_type,
           'gene' = {
               gene_targets
           },
           'protein' = {
               readr::read_csv('https://www.dropbox.com/s/pvcf4y8x4hbkm89/ab_info.csv?raw=1') %>%
                   dplyr::left_join(gene_targets,
                                    by = c('gene_id' = 'feature')) %>%
                   dplyr::select(mirna_base, ab_id) %>%
                   setNames(c('mirna_base', 'feature')) %>%
                   na.omit()
           })

}

