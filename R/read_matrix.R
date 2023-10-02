#' read_matrix
#'
#' @param path
#'
#' @return
#' @export
#' @importFrom data.table fread melt setnames
#'
#' @examples
read_matrix <- function(path){
  # Read PA matrix using data.table and covert back to matrix with named rows and columns
  dt <- data.table::fread(path)
  dt <- data.table::melt(dt, id.vars = "V1")
  dt <- dt[, .(observed = sum(value), n_samples = .N), "variable"]
  data.table::setnames(dt, "variable", "gene")
  dt[, crude_freq := observed / n_samples]

  return(dt)
}
