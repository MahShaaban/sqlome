context("test functios for transcription factors data")

test_that("mrna_tidy returns a proper matrix", {
    library(RTCGA.rnaseq)
    rtcga.df <- ACC.rnaseq
    ACC.mat <- mrna_tidy(rtcga.df)

    expect_equal(class(ACC.mat), 'matrix')
    expect_false(all(duplicated(colnames(ACC.mat))))
    expect_false(is.null(rownames(ACC.mat)))
})

test_that("mirna_tidy returns a proper matrix", {
    library(RTCGA.miRNASeq)
    rtcga.df <- ACC.miRNASeq
    ACC.mat <- mirna_tidy(rtcga.df)

    expect_equal(class(ACC.mat), 'matrix')
    expect_false(all(duplicated(colnames(ACC.mat))))
    expect_false(is.null(rownames(ACC.mat)))
})

test_that("rppa_tidy returns a proper matrix", {
    library(RTCGA.RPPA)
    rtcga.df <- ACC.RPPA
    ACC.mat <- rppa_tidy(rtcga.df)

    expect_equal(class(ACC.mat), 'matrix')
    expect_false(all(duplicated(colnames(ACC.mat))))
    expect_false(is.null(rownames(ACC.mat)))
})

test_that("cor_make calculate correlation for genes", {
    library(RTCGA.rnaseq)
    library(RTCGA.miRNASeq)

    cohort <- 'ACC'
    ACC.mi <- mirna_tidy(ACC.miRNASeq)
    ACC.m <- mrna_tidy(ACC.rnaseq)
    corr <- cor_make(ACC.mi, ACC.m, cohort = cohort)

    expect_s3_class(corr, 'data.frame')
    expect_identical(names(corr), c('feature', 'mirna_base', cohort))
    expect_true(abs(min(corr$ACC)) >= 10)
})

test_that("cor_make calculate correlation for proteins", {
    library(RTCGA.rnaseq)
    library(RTCGA.RPPA)

    cohort <- 'ACC'
    ACC.mi <- mirna_tidy(ACC.miRNASeq)
    ACC.r <- rppa_tidy(ACC.RPPA)
    corr <- cor_make(ACC.mi, ACC.r, cohort = cohort)

    expect_s3_class(corr, 'data.frame')
    expect_identical(names(corr), c('feature', 'mirna_base', cohort))
    expect_true(abs(min(corr$ACC)) >= 10)
})

test_that('get_targets extract the right taregets for genes', {
    library(targetscan.Hs.eg.db)
    library(org.Hs.eg.db)

    mirna_base <- c("hsa-let-7a", "hsa-let-7b", "hsa-let-7c")
    gene <- get_targets(mirna_base, 'gene')

    expect_s3_class(gene, 'data.frame')
    expect_identical(names(gene), c('mirna_base', 'feature'))
    expect_true(all(unique(gene$mirna_base) %in% mirna_base))

    ptn <- get_targets(mirna_base, 'protein')
    expect_s3_class(ptn, 'data.frame')
    expect_identical(names(ptn), c('mirna_base', 'feature'))
    expect_true(all(unique(ptn$mirna_base) %in% mirna_base))
})
