solve_for_pi_NB <- function(gene_mean, dispersion, n = 50, one_sided = TRUE){ #alpha = 0.05

  use_k <- determine_k(gene_mean)
  expected_pi <- compute_expected_pi_for_k_NB(gene_mean, use_k, dispersion)

  if(one_sided == TRUE){
    zscore_alpha <- 1.645
  } else {
    zscore_alpha <- 1.96
  }
  if(one_sided == TRUE){
    diff_est <- function(what_k_prop){
      est_se <- sqrt(expected_pi*(1-expected_pi)/n)
      propdiff <- what_k_prop - expected_pi
      zvalue <- propdiff/est_se
      abs(zvalue-zscore_alpha)
    }
    return(c(upper = optim(0.5, diff_est, method = "L-BFGS-B", lower = 0, upper = 1)$par))
  } else {

    diff_est <- function(what_k_prop){
      est_se <- sqrt(expected_pi*(1-expected_pi)/n)
      propdiff <- what_k_prop - expected_pi
      zvalue <- propdiff/est_se
      abs(zvalue-zscore_alpha)
    }
    diff_est2 <- function(what_k_prop){
      est_se <- sqrt(expected_pi*(1-expected_pi)/n)
      propdiff <- what_k_prop - expected_pi
      zvalue <- propdiff/est_se
      abs(zvalue+zscore_alpha)
    }
    return(c(upper = optim(0.5, diff_est, method = "L-BFGS-B", lower = 0, upper = 1)$par,
             lower = optim(0.5, diff_est2, method = "L-BFGS-B", lower = 0, upper = 1)$par))

  }

}


#' waveplot
#'
#' @param one.cluster.zscore.summary output from Kprop_inflation_test()
#' @param separate.RNA Default is TURE. Whether to plot IncRNA, mtRNA, or ribosomal RNA. When set to FALSE, all RNAs will be plotted.
#' @param one_sided Default is TRUE. Whether to draw the confident interval based on a one-sided test.
#' @param show_mean_cutoff Default is 1.5. Cut-off on the gene mean for inflated genes.
#' @param show_CI Default is FALSE. Whether to draw the confidence intervals.
#' @param n_CI Default is 50. Number of points to quantify the confidence intervals.
#' @param show_top Default is 25. Print names of top X genes.
#' @param xmax Default is 30. X-axis up-limit.
#' @param log_scale Default is FALSE. Whether to show X-axis (mean UMI counts) in log-scale.
#' @param add_poisson_line Whether to add expected k-proportions calculated from a Poisson distribution.
#' @return A ggplot2 object.
#' @export
#'


waveplot <- function(one.cluster.zscore.summary, separate.RNA = TRUE, one_sided = TRUE,
                     show_mean_cutoff = 1.5, show_CI = FALSE, n_CI = 50, show_top = 25,
                     xmax = 30, log_scale = FALSE, add_poisson_line = FALSE){

  dispersion <- one.cluster.zscore.summary$NBdispersion
  selected_stats <- one.cluster.zscore.summary$Inflated
  selected_stats <- selected_stats[selected_stats$gene_mean >= show_mean_cutoff, ]

  xmax <- min(xmax, max(round(selected_stats$gene_mean))+3)
  if(xmax <= 5) xmax <- 5
  hits_num <- nrow(selected_stats)
  show_top <- min(hits_num, show_top)
  showhits <- selected_stats[1:show_top, ]

  score_mat <- one.cluster.zscore.summary$Zscores
  all_genes <- rownames(score_mat)

  selected.IDs <- which(all_genes %in% rownames(selected_stats))

  if(separate.RNA == TRUE){
    data("otherRNA")
    MT_genes <- otherRNA[otherRNA[, "Family"] == "Mitochondiral", "Gene"]
    RP_genes <- otherRNA[otherRNA[, "Family"] == "Ribosomal", "Gene"]
    IncRNA_genes <- otherRNA[otherRNA[, "Family"] == "IncRNA", "Gene"]

    MT.IDs <- which(all_genes%in%MT_genes)
    RP.IDs <- which(all_genes%in%RP_genes)
    IncRNA.IDs <- which(all_genes%in%IncRNA_genes)

    MT_stats <- score_mat[MT.IDs, ]
    RP_stats <- score_mat[RP.IDs, ]
    IncRNA_stats <- score_mat[IncRNA.IDs, ]

    background_stats <- score_mat[-c(selected.IDs, MT.IDs, RP.IDs, IncRNA.IDs), ]

    stats_combined <- data.frame(category = c(rep("B", nrow(background_stats)), rep("H", nrow(selected_stats)), rep("MT", nrow(MT_stats)),
                                              rep("RP", nrow(RP_stats)), rep("IncRNA", nrow(IncRNA_stats))),
                                 rbind(background_stats, selected_stats, MT_stats, RP_stats, IncRNA_stats))
  } else {
    background_stats <- score_mat[-selected.IDs, ]
    stats_combined <- data.frame(category = c(rep("B", nrow(background_stats)), rep("H", nrow(selected_stats))),
                                 rbind(background_stats, selected_stats))
  }

  x_intervals <- seq(0.1, xmax, by = 0.2)
  p_intervals <- rep(0, length(x_intervals))
  logx_intervals <- log10(x_intervals + 0.1)

  if(show_CI == FALSE){
    upper_lower_intervals <- matrix(0, ncol = 1, nrow = length(x_intervals))
    for(i in 1:length(x_intervals)){
      p_intervals[i] <- determine_k(x_intervals[i])
      upper_lower_intervals[i, 1] <- compute_expected_pi_for_k_NB(x_intervals[i], p_intervals[i], dispersion)
    }
    expected_line <- data.frame(aa = x_intervals, bb = logx_intervals, dd = upper_lower_intervals[, 1])
  }

  if(one_sided == TRUE & show_CI == TRUE){
    upper_lower_intervals <- matrix(0, ncol = 2, nrow = length(x_intervals))
    for(i in 1:length(x_intervals)){
      p_intervals[i] <- determine_k(x_intervals[i])
      upper_lower_intervals[i, 1] <- compute_expected_pi_for_k_NB(x_intervals[i], p_intervals[i], dispersion)
      upper_lower_intervals[i, 2] <- solve_for_pi_NB(x_intervals[i], dispersion, n = n_CI, one_sided = one_sided)
    }
    expected_line <- data.frame(aa = x_intervals, bb = logx_intervals, dd = upper_lower_intervals[, 1])
    CI_line <- data.frame(aa = x_intervals, bb = logx_intervals, ee = upper_lower_intervals[, 2])
  }


  if(one_sided == FALSE & show_CI == TRUE){
    upper_lower_intervals <- matrix(0, ncol = 3, nrow = length(x_intervals))
    for(i in 1:length(x_intervals)){
      p_intervals[i] <- determine_k(x_intervals[i])
      upper_lower_intervals[i, 1] <- compute_expected_pi_for_k_NB(x_intervals[i], p_intervals[i], dispersion)
      upper_lower_intervals[i, 2:3] <- solve_for_pi_NB(x_intervals[i], dispersion, n = n_CI, one_sided = one_sided)
    }
    expected_line <- data.frame(aa = x_intervals, bb = logx_intervals, dd = upper_lower_intervals[, 1])
    CI_line <- data.frame(aa = x_intervals, bb = logx_intervals, ee = upper_lower_intervals[, 2])
    CI_line2 <- data.frame(aa = x_intervals, bb = logx_intervals, ff = upper_lower_intervals[, 3])
  }

  Poisson_line <- data.frame(aap = x_intervals, bbp = logx_intervals, eep = ppois(p_intervals, x_intervals))

  if(log_scale == FALSE){
    require(ggplot2)
    g_basic <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = gene_mean,
                                                            y = k_proportion, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
      ggplot2::ylab("K Proportion") + ggplot2::xlab("Gene Mean") + xlim(0, xmax) + ylim(0, 1) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside") + ggplot2::geom_line(data = expected_line, aes(aa, dd), color = "turquoise4", show.legend = FALSE)

    if(separate.RNA == TRUE ){

      g <- g_basic + ggplot2::scale_color_manual(values = c("B" = "black", "IncRNA" = "#0099DD", "MT" = "#C2BB00", "RP" = "#8C1F66", "H" ="#F06060"),
                                                 labels = c("B" = "Background", "IncRNA" = "IncRNA", "MT" = "Mt Genes", "RP" = "Ribo Genes", "H" = "Inflated"))

    } else {
      g <- g_basic +  ggplot2::scale_color_manual(values = c("B" = "black", "H" ="#F06060"),
                                                  labels = c("B" = "Background", "H" = "Inflated"))
    }

    if(show_CI == TRUE){
      g <- g + ggplot2::geom_line(data = CI_line, aes(aa, ee), color = "darkgreen", linetype = "dashed", show.legend = FALSE)
      if(one_sided == FALSE){
        g <- g + ggplot2::geom_line(data = CI_line2, aes(aa, ff), color = "darkgreen", linetype = "dashed", show.legend = FALSE)
      }
    }

    if(add_poisson_line == TRUE){
      g <- g + ggplot2::geom_line(data = Poisson_line, aes(aap, eep), color = "blue", linetype = "dashed", show.legend = FALSE)

    }


  } else {

    stats_combined[, "gene_mean"] <- log10(stats_combined[, "gene_mean"]+0.1)
    showhits[, "gene_mean"] <- log10(showhits[, "gene_mean"]+0.1)
    xmax <- max(showhits[, "gene_mean"]) + 0.2

    require(ggplot2)
    g_basic <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = gene_mean,
                                                            y = k_proportion, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
      ggplot2::ylab("K Proportion") + ggplot2::xlab("log10 Gene Mean") + xlim(-1, xmax) + ylim(0, 1) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside") + ggplot2::geom_line(data = expected_line, aes(bb, dd), color = "turquoise4", show.legend = FALSE)

    if(separate.RNA == TRUE ){

      g <- g_basic + ggplot2::scale_color_manual(values = c("B" = "black", "IncRNA" = "#0099DD", "MT" = "#C2BB00", "RP" = "#8C1F66", "H" ="#F06060"),
                                                 labels = c("B" = "Background", "IncRNA" = "IncRNA", "MT" = "Mt Genes", "RP" = "Ribo Genes", "H" = "Inflated"))

    } else {
      g <- g_basic +  ggplot2::scale_color_manual(values = c("B" = "black", "H" ="#F06060"),
                                                  labels = c("B" = "Background", "H" = "Inflated"))
    }



    if(show_CI == TRUE){
      g <- g + ggplot2::geom_line(data = CI_line, aes(bb, ee), color = "darkgreen", linetype = "dashed", show.legend = FALSE)
      if(one_sided == FALSE){
        g <- g + ggplot2::geom_line(data = CI_line2, aes(bb, ff), color = "darkgreen", linetype = "dashed", show.legend = FALSE)
      }
    }


    if(add_poisson_line == TRUE){
      g <- g + ggplot2::geom_line(data = Poisson_line, aes(bbp, eep), color = "darkblue", linetype = "dashed", show.legend = FALSE)

    }
  }



  return(g)

}
