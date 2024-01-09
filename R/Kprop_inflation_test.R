compute_expected_pi_for_k_NB <- function(gene_mean, k, dispersion) {
  n_para <- 1/dispersion
  p_para <- (n_para)/(n_para + gene_mean)
  expected_pi <- pnbinom(k, size = n_para, prob = p_para)
  return(expected_pi)
}

k_prop_to_one_gene <- function(data.x, cell.number){
  if(sum(data.x) == 0){
    return(c(gene_mean = 0, k = 0, k_proportion = 1))
  } else {
    gene_mean <- sum(data.x)/cell.number
    use_k <- determine_k(gene_mean)
    k_proportion <- sum(data.x <= use_k)/cell.number
    return(c(gene_mean = gene_mean, k = use_k, k_proportion = k_proportion))
  }
}

est_var_asymptotic <- function(k_proportion, gene_mean, dispersion, k){

  rrr <- 1/dispersion
  rrr_gene_mean <- rrr + gene_mean
  p_para <- rrr/rrr_gene_mean

  expected_pi <- pnbinom(k, size = rrr, prob = p_para)

  est_var_term1 <- expected_pi*(1-expected_pi) #k_proportion*(1-k_proportion)
  var_x <- gene_mean + gene_mean^2*dispersion

  one_k_function_dev <- function(k_val, rrr, p_para, gene_mean, rrr_gene_mean){
    a1 <- (k_val-gene_mean)*rrr/(gene_mean*(rrr_gene_mean))
    a1*dnbinom(k_val, size = rrr, prob = p_para)
  }

  F_dev <- sum(one_k_function_dev(0:k, rrr, p_para, gene_mean, rrr_gene_mean))

  est_var_term2 <- F_dev^2*var_x

  one_k_function_cov <- function(k_val, rrr, p_para){
    k_val*dnbinom(k_val, size = rrr, prob = p_para)
  }

  cov_X <- sum(one_k_function_cov(0:k, rrr, p_para)) - expected_pi*gene_mean # k_proportion*gene_mean

  est_var_term3 <- -2*F_dev*cov_X

  est_var <- est_var_term1 + est_var_term2 + est_var_term3
  return(list(expected_pi = expected_pi, est_var = est_var))
}

compute_zscore_asymp <- function(k_proportion, gene_mean, k, dispersion, cell.number) {

  est_var_list <- est_var_asymptotic(k_proportion, gene_mean, dispersion, k)
  expected_pi <- est_var_list$expected_pi
  est_var_asym <- est_var_list$est_var
  est_se <- sqrt(est_var_asym/cell.number)
  propdiff <- k_proportion - expected_pi
  zvalue <- propdiff/est_se
  return(list(zscore = zvalue, expected_pi = expected_pi, est_var_asym = est_var_asym))
}

determine_k <- function(gene_mean){
  possible_n <- 1:ceiling(gene_mean)
  bounds_n <- possible_n + sqrt(possible_n)
  if(gene_mean >= 6){
    diff_mean <- bounds_n - gene_mean
    k <- which(diff_mean>0)[1] - 1
  }
  if(gene_mean < 6 & gene_mean >= 4.732){ # 3+sqrt(3)
    k <- 3
  }
  if(gene_mean < 4.732 & gene_mean >= 3.414){ #sqrt(2) + 2
    k <- 2
  }
  if(gene_mean < 3.414 & gene_mean >= 2 ){
    k <- 1
  }
  if(gene_mean < 2){
    k <- 0
  }
  return(k)
}

solve_for_dispersion <- function(k_prop_mat, k_cutoff = 1){
  k_prop_mat <- k_prop_mat[k_prop_mat$k >= k_cutoff, ]
  diff_est <- function(dispersion){
    num_gene <- nrow(k_prop_mat)
    fitted_prop <- rep(0, num_gene)
    for(i in 1:num_gene) {
      fitted_prop[i] <- compute_expected_pi_for_k_NB(k_prop_mat$gene_mean[i], k_prop_mat$k[i], dispersion)
    }
    sum((fitted_prop - k_prop_mat$k_proportion)^2)
  }
  return(c(optim(0.5, diff_est, method = "L-BFGS-B", lower = 0.01, upper = 10)$par))
}

deviation_test_to_one_gene_NB_asymp <- function(data.x, cell.number, dispersion){
  if(sum(data.x) == 0){
    return(c(gene_mean = 0, k = 0, kprop_diff = 0, k_proportion = 1, expected_pi = 1, zscore = NA, est_var_asym = NA))
  } else {
    gene_mean <- sum(data.x)/cell.number
    use_k <- determine_k(gene_mean)
    k_proportion <- sum(data.x <= use_k)/cell.number
    zscore_est <- compute_zscore_asymp(k_proportion, gene_mean, use_k, dispersion, cell.number)
    return(c(gene_mean = gene_mean, k = use_k, kprop_diff = k_proportion - zscore_est$expected_pi, k_proportion = k_proportion, expected_pi = zscore_est$expected_pi, zscore = zscore_est$zscore, est_var_asym = zscore_est$est_var_asym))
  }
}


Estimate_dispersion_zscore_asymp_for_one_cluster <- function(X, k_cutoff = 1, fdr_control_method = "BH") {

  cell.number <- ncol(X)
  all.stats1 <- apply(X, 1, function(x){ k_prop_to_one_gene(x, cell.number)})
  score_mat1 <- as.data.frame(t(all.stats1))
  dispersion_est <- solve_for_dispersion(score_mat1, k_cutoff = k_cutoff)

  all.stats2 <- apply(X, 1, function(x){ deviation_test_to_one_gene_NB_asymp(x, cell.number, dispersion_est)})
  score_mat <- as.data.frame(t(all.stats2))
  score_mat$pvalue <- pnorm(score_mat$zscore, lower.tail = FALSE)
  score_mat$qvalue <- p.adjust(score_mat$pvalue, n = nrow(score_mat), method = fdr_control_method)

  return(list(score_mat = score_mat, dispersion_est = dispersion_est))
}


#' Kprop_inflation_test
#'
#' @param dat Input data matrix (preferably UMI counts) for gene \times cell.
#' @param grouping The grouping IDs of cells. If the grouping is NULL, all cells will be assessed at once. The grouping can be cell types, individuals or any phenotype.
#' @param min.cell.num.group Minimum cell number required for each group. The default is 25.
#' @param min.depth.group Minimum average depth required for each group. The default is 4000, set for 10X data.
#' @param cell.depth.ranges A two-value vector of range of total counts per cell. The default is c(2500, 60000), set for 10X data.
#' @param kprop.cutoff A pre-specified cut-off on K props. Only k >= k_cutoff will be used for estimating the dispersion of NB. The default is 1, meaning not using 0-prop for estimation.
#' @param fdr.control.method Method to adjust for multiple comparisons. The default is "BH". See function p.adjust() for more options.
#' @param qvalue.cutoff A cut-off on q-value adjusted by multiple comparisons, to call significantly inflated genes. The default is 0.05.
#' @param genemean.cutoff A cut-off on mean UMI per cell, to report in the list of inflated genes. The default is 1.
#' @param filter.gene Whether reporting IncRNA, mtRNA, or ribosomal RNA in the list of inflated genes. The default is TRUE.
#' @return A list of summary of k-proportion test for each group.
#' \itemize{
#'   \item groupID - the ID of the group. The value is NA if there is no grouping information.
#'   \item NBdispersion - estimated dispersion parameter
#'   \item Zscores - a matrix summarizing gene mean, zscore, p-value and q-value for all genes in dat. NAs are reported for genes with mean = 0.
#'   \item Inflated - a matrix of selected genes that pass qvalue.cutoff and genemean.cutoff.
#' }
#' @export

Kprop_inflation_test <- function(dat, grouping = NULL,
                                        min.cell.num.group = 25, min.depth.group = 4000, cell.depth.ranges = c(2500, 60000),
                                        kprop.cutoff = 1, fdr.control.method = "BH", qvalue.cutoff = 0.05, genemean.cutoff = 0.5, kprop.diff.cutoff = 0.05, filter.gene = TRUE){
  if(!is.null(grouping)){
    if(ncol(dat) != length(grouping)) stop("Cell numbers and grouping labels do not match.")
    else {
      cell_total_count <- apply(dat, 2, sum)
      cell_flag <- which(cell_total_count >= cell.depth.ranges[1] & cell_total_count <= cell.depth.ranges[2])
      dat <- dat[, cell_flag]
      grouping <- grouping[cell_flag]
      cell_total_count <- cell_total_count[cell_flag]

      min.depth <- log10(min.depth.group+1)
      comb_tab <- data.frame(cellnumber = table(grouping), meancount = log10(tapply(cell_total_count, grouping, mean)+1))
      colnames(comb_tab)[1:2] <- c("category", "cellnumber")
      selected_groups <- comb_tab[which(comb_tab[, "cellnumber"] >= min.cell.num.group & comb_tab[, "meancount"] >= min.depth), ]

      groupIDs <- unique(selected_groups[, 1])
      group_zscore_summary <- list(NULL)

      for(j in 1:length(groupIDs)){

        one_group_mat <- dat[, grouping==groupIDs[j]]
        res_onegroup <- Estimate_dispersion_zscore_asymp_for_one_cluster(one_group_mat, k_cutoff = kprop.cutoff, fdr_control_method = fdr.control.method)
        scores <- res_onegroup$score_mat
        rownames(scores) <- rownames(one_group_mat)
        hits <- scores[scores$qvalue <= qvalue.cutoff & !is.na(scores$zscore) & scores$gene_mean >= genemean.cutoff & scores$zscore > 0 & scores$kprop_diff >= kprop.diff.cutoff, ]
        if(filter.gene == TRUE & !is.null(rownames(hits))){
          data("otherRNA")
          filtered.id <- c(which(rownames(hits)%in%otherRNA[, "Gene"]), grep("^LINC", rownames(hits)), grep("^A[A-Z][0-9][0-9][0-9][0-9][0-9][0-9]", rownames(hits)))
          if(length(filtered.id) > 0) {
            hits <- hits[-filtered.id, ]
          }
        }
        hits <- hits[order(hits$zscore, decreasing = TRUE), ]
        group_zscore_summary[[j]] <- list(groupID = groupIDs[j], NBdispersion = res_onegroup$dispersion_est, Zscores = scores, Inflated = hits)

      }
    }
  } else {
    cell_total_count <- apply(dat, 2, sum)
    cell_flag <- which(cell_total_count >= cell.depth.ranges[1] & cell_total_count <= cell.depth.ranges[2])
    dat <- dat[, cell_flag]

    res_onegroup <- Estimate_dispersion_zscore_asymp_for_one_cluster(dat, k_cutoff = kprop.cutoff, fdr_control_method = fdr.control.method)
    scores <- res_onegroup$score_mat
    rownames(scores) <- rownames(dat)
    hits <- scores[scores$qvalue <= qvalue.cutoff & !is.na(scores$zscore) & scores$gene_mean >= genemean.cutoff & scores$zscore > 0 & scores$kprop_diff >= kprop.diff.cutoff, ]
    if(filter.gene == TRUE & !is.null(rownames(hits))){
      data(otherRNA)
      filtered.id <- c(which(rownames(hits)%in%otherRNA[, "Gene"]), grep("^LINC", rownames(hits)), grep("^A[A-Z][0-9][0-9][0-9][0-9][0-9][0-9]", rownames(hits)))
      if(length(filtered.id) > 0) {
        hits <- hits[-filtered.id, ]
      }
    }
    hits <- hits[order(hits$zscore, decreasing = TRUE), ]
    group_zscore_summary <- list(groupID = NA, NBdispersion = res_onegroup$dispersion_est, Zscores = scores, Inflated = hits)
  }

  return(group_zscore_summary)
}

