#' Title
#'
#' @param core_threshold The true frequency above which a gene is considered a core gene (numeric)
#' @param completeness A numeric vector of completeness values for each genome sample from CheckM, must be between 0 and 1
#' @param observations A vector of the number of observations of each gene, must be between 0 and the number of genome samples
#' @param n.cores The number cores to use for parallel calculation, useful when there are >1000 genome samples
#'
#' @return
#' @export
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom ggplot2 ggplot geom_bar labs ggtitle scale_x_continuous theme annotate theme_bw aes element_text
#' @import patchwork
#'
#' @examples
#' cgt(core_threshold = 0.95,
#'     completeness = runif(20, 0.95, 1),
#'     quantiles = c(0.3, 0.5),
#'     observations = rbinom(140, size = 20, prob = 0.90))
#'
cgt <- function(core_threshold = 0.95,
                        completeness, observations,
                        quantiles = c(0.05, 0.25, 0.5),
                        n.cores = 1) {
  N <- length(completeness)
  true_freq <- core_threshold

  # A function to use with the integration method
  f <- function(true_freq, y, completeness){
    out <- c()
    for(i in 1:length(true_freq)) {
      p <- true_freq[i] * completeness
      out[i] <- poibin::dpoibin(kk = y, pp = p)
    }
    out
  }

  # Calculate posterior values for a subset of possible observations
  cl <- parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl)
  resvec <- foreach::foreach(i = 1:N, .combine = "c") %dopar% {
    integrate(f = f,
              lower = true_freq,
              upper = 1, y = i,
              completeness = completeness)$value / (1 - true_freq)
  }
  parallel::stopCluster(cl)

  # Insert bit to choose where to calculate from
  # Then need to add cumulative probability from P(x <= start) to cumulative sum

  res <- data.frame(x = 1:N, prob = resvec)
  res$cum_prob <- cumsum(res$prob)

  # Create plots
  plot_data <- res[res$prob > 0.0001, ]

  if(nrow(plot_data) <= 50){
    declutter_index <- 1:nrow(plot_data)
  }else if(N <= 100){
    declutter_index <- seq(1, nrow(plot_data), 5)
  }else{
    declutter_index <- seq(1, nrow(plot_data), 10)
  }

  bar_labels <- paste0(signif(plot_data$cum_prob, 2) * 100, "%")[declutter_index]

  p1 <-
    ggplot2::ggplot(data = plot_data)+
    ggplot2::geom_bar(ggplot2::aes(x, prob),
                      fill = "lightgrey",
                      stat = "identity") +
    ggplot2::labs(y = "Probability", x = "Number of observations of gene") +
    ggplot2::ggtitle(paste0("Bars show P(N observations | true frequency >=", true_freq * 100, "%) \n",
                   "Text shows ", paste0("P(observations <= N | true frequency >=",true_freq * 100, "%)"))) +
    ggplot2::scale_x_continuous(breaks = floor(seq(min(plot_data$x),
                                                   max(plot_data$x),
                                                   length.out = 10))) +
    ggplot2::theme(plot.title = ggplot2::element_text()) +
    ggplot2::theme_bw() +
    ggplot2::annotate("text", label = bar_labels,
             x = plot_data$x[declutter_index], y = (plot_data$prob + max(plot_data$prob)/50)[declutter_index],
             fontface = "bold")

  n_genes_core <- c()
  for(i in plot_data$x){
    n_genes_core <- c(n_genes_core, sum(observations >= i))
  }

  core_genome_sizes <- data.frame(x = plot_data$x, y = n_genes_core)

  p2 <-
    ggplot2::ggplot(data = core_genome_sizes) +
    ggplot2::geom_bar(ggplot2::aes(x, y),
                      fill = "lightgrey",
                      stat = "identity") +
    ggplot2::labs(x = "Value of N used for core genome assignment (>=N considered core)",
         y = "Number of genes considered core") +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste0("Effect of cut-off choice on subsequent core genome size \n (out of ", length(observations), " genes)"))

  plot <- p1 + p2 & ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # Create summary table
  quantiles
  quant_x <- c()
  for(i in 1:length(quantiles)){
    quant_x[i] <- res$x[which(res$cum_prob > quantiles[i])[1]]
  }

  return(list(plot = plot, table = data.frame(quantile = quantiles, gene_observations = quant_x)))
}
