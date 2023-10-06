#' Title
#'
#' @param core_threshold The true frequency above which a gene is considered a core gene (numeric)
#' @param completeness A numeric vector of completeness values for each genome sample from CheckM, must be between 0 and 1
#' @param observations A vector of the number of observations of each gene, must be between 0 and the number of genome samples
#' @param n.cores The number cores to use for parallel calculation, useful when there are >1000 genome samples
#'
#' @return
#' @export
#' @importFrom ggplot2 ggplot geom_bar labs ggtitle scale_x_continuous theme annotate theme_bw aes position_dodge element_text geom_line scale_fill_discrete
#' @import patchwork
#'
#' @examples
#' cgt(core_threshold = 0.95,
#'     completeness = runif(20, 0.95, 1),
#'     quantiles = c(0.3, 0.5),
#'     observations = rbinom(140, size = 20, prob = 0.90))
#'
cgt <- function(core_threshold = 0.95,
                rare_threshold = 0.10,
                completeness,
                observations,
                quantiles = c(0.05, 0.25, 0.5),
                plot = TRUE,
                betap1 = 1,
                betap2 = 1,
                n_posterior_samples = 10000,
                model_object) {
  n_samples <- length(completeness)

  stanfit <- model_object$sample(data = list(N = n_samples,
                                    completeness = completeness,
                                    betap1 = betap1,
                                    betap2 = betap2,
                                    core_threshold = core_threshold,
                                    rare_threshold = rare_threshold),
                                 refresh = 1000,
                        fixed_param = TRUE, iter_sampling = n_posterior_samples, chains = 1)

  res_core <- as.vector(stanfit$draws("obs_core"))
  res_notcore <- as.vector(stanfit$draws("obs_notcore"))
  res_rare <- as.vector(stanfit$draws("obs_rare"))
  res_notrare <- as.vector(stanfit$draws("obs_notrare"))

  # Calculate wrong assignment probabilities
  prob_core_as_notcore <- c()
  prob_notcore_as_core <- c()
  prob_rare_as_notrare <- c()
  prob_notrare_as_rare <- c()

  for(i in 1:n_samples){
    prob_core_as_notcore[i] <- sum(res_core < i) / n_posterior_samples
    prob_notcore_as_core[i] <- sum(res_notcore >= i) / n_posterior_samples

    prob_rare_as_notrare[i] <- sum(res_rare > i) / n_posterior_samples
    prob_notrare_as_rare[i] <- sum(res_notrare <= i) / n_posterior_samples
  }

  core_df <- data.frame(x = 1:n_samples, prob_core_as_notcore, prob_notcore_as_core)
  rare_df <- data.frame(x = 1:n_samples, prob_rare_as_notrare, prob_notrare_as_rare)

  # Create summary table for quantiles
  quant_core <- quant_rare <- c()
  for(i in 1:length(quantiles)){
    if(any(core_df$prob_core_as_notcore > quantiles[i])){
      quant_core[i] <- core_df$x[max(which(core_df$prob_core_as_notcore < quantiles[i]))]
    } else {
      quant_core[i] <- nrow(core_df)
    }

    if(any(rare_df$prob_rare_as_notrare > quantiles[i])){
      quant_rare[i] <- rare_df$x[min(which(rare_df$prob_rare_as_notrare < quantiles[i]))]
    }else{
      quant_rare[i] <- 1
    }
  }

  print(rare_df)

  quantiles_df <- data.frame("Probability of x as not x" = paste0("<", quantiles * 100, "%"),
             "core cutoff" = quant_core,
             "rare cutoff" = quant_rare)

  # Plot beta prior shape of underlying gene frequencies
  dt_prior <- data.frame(x = as.vector(stanfit$draws("prior_sample")))

  prior_plot <- ggplot2::ggplot(data = dt_prior, ggplot2::aes(x = x)) +
    ggplot2::geom_density() +
    ggplot2::ggtitle("Prior distribution of underlying population frequencies") +
    ggplot2::labs(x = "Underlying population frequency", y = "Prior probability") +
    ggplot2::theme_bw() +
    ggplot2::geom_vline(xintercept = c(rare_threshold, core_threshold), col = "red")

  # Plot assignment errors
  core_plot_df <- core_df[core_df$x %in% min(quantiles_df$core.cutoff):max(quantiles_df$core.cutoff),]
  core_plot_df <- data.frame(x = rep(core_plot_df$x, 2),
                        variable = c(rep("core as not core", nrow(core_plot_df)),
                                     rep("not core as core", nrow(core_plot_df))),
                        value = c(core_plot_df$prob_core_as_notcore,
                                  core_plot_df$prob_notcore_as_core))

  p_core <- ggplot2::ggplot(data = core_plot_df, ggplot2::aes(x = x, y = value, fill = variable)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
    ggplot2::ggtitle(paste0("Observation cut offs for a core threshold of ", core_threshold * 100, "%")) +
    ggplot2::scale_fill_discrete(labels = c("P(assign core as not core)", "P(assign not core as core)"), name = "") +
    ggplot2::scale_x_continuous(breaks = 0:n_samples) +
    ggplot2::labs(x = "", y = "Probability") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(hjust = 0.5))


  rare_plot_df <- rare_df[rare_df$x %in% min(quantiles_df$rare.cutoff):max(quantiles_df$rare.cutoff),]
  # rare_plot_df <- rare_df
  rare_plot_df <- data.frame(x = rep(rare_plot_df$x, 2),
                             variable = c(rep("rare as not rare", nrow(rare_plot_df)),
                                          rep("not rare as rare", nrow(rare_plot_df))),
                             value = c(rare_plot_df$prob_rare_as_notrare,
                                       rare_plot_df$prob_notrare_as_rare))

  p_rare <- ggplot2::ggplot(data = rare_plot_df, ggplot2::aes(x = x, y = value, fill = variable)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
    ggplot2::ggtitle(paste0("Observation cut offs for a rare threshold of ", rare_threshold * 100, "%")) +
    # ggplot2::scale_fill_discrete() +
    ggplot2::scale_fill_brewer(labels = c("P(assign not rare as rare)", "P(assign rare as not rare)"), name = "", palette = "Set2") +
    ggplot2::scale_x_continuous(breaks = 0:n_samples) +
    ggplot2::labs(x = "", y = "Probability") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(hjust = 0.5))

  # Number of genes in core for thresholds
  n_genes_core <- c()
  n_genes_rare <- c()
  for(i in 1:nrow(quantiles_df)){
    n_genes_core[i] <- sum(observations >= quantiles_df$core.cutoff[i])
    n_genes_rare[i] <- sum(observations <= quantiles_df$rare.cutoff[i])
  }

  quantiles_df$core.size <- n_genes_core
  quantiles_df$rare.size <- n_genes_rare

  # p_sizes <-
  #   ggplot2::ggplot(data = core_genome_sizes) +
  #   ggplot2::geom_bar(ggplot2::aes(x, y),
  #                     fill = "lightgrey",
  #                     stat = "identity") +
  #   ggplot2::labs(x = "Value of N used for core genome assignment (>=N considered core)",
  #                 y = "Number of genes considered core") +
  #   ggplot2::theme_bw() +
  #   ggplot2::ggtitle(paste0("Effect of cut-off choice on subsequent core genome size \n (out of ", length(observations), " genes)"))

  plot(prior_plot + p_core + p_rare)

  return(quantiles_df)
}

#' Title
#'
#' @param model
#'
#' @return
#' @export
#' @importFrom cmdstanr cmdstan_model
#'
#' @examples
stanmodel <- function(model = system.file("stan", "sampler.stan", package = "cgt")) {
  mod <- cmdstanr::cmdstan_model(model)
  return(mod)
}
