#' Simulated unspliced and spliced counts for a single gene
#'
#' @description
#' Unspliced and spliced counts for 200 observations (cells) are simulated for a single gene, under the model.
#'
#' @format A list containing the observed counts, mean counts from the ODE, and parameters to generate the data.
"sim.data"



#' Simulated unspliced and spliced counts for 26 genes as well as posterior summaries from individual fit
#'
#' @description
#' Unspliced and spliced counts for 400 observations (cells) are simulated for 26 single gene, under the model, and
#' the posterior summaries from individual fit, which is needed for the post-processing step.
#'
#' @format A list containing the following components:
#' \describe{
#' \enumerate{
#' \item \code{u.obs} and \code{s.obs} are generated counts for 400 cells (row) and 26 genes (column).
#' \item \code{x} is a cell-by-gene matrix where each entry is the posterior mean of log gene-wise latent time,
#' conditional on the state with a larger posterior probability \code{stored in (k.mode)}.
#' \item \code{k.mode} is a cell-by-gene matrix. Each entry stores the state (induction: 1, repression: 2)
#' with a larger posterior probability.
#' \item \code{var.logt} an array \code{[cell, gene, state]} of posterior variance of log gene-wise latent time, for each state.
#' \item \code{p.hat} an array \code{[cell, gene, state]} of posterior probability of belonging to each state for each cell in each gene.
#' }
#' }
"post.summary"
