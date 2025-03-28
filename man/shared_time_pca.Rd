% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post_process.R
\name{shared_time_pca}
\alias{shared_time_pca}
\title{Visualize shared time on PCA space}
\usage{
shared_time_pca(u.obs, s.obs, shared_time)
}
\arguments{
\item{u.obs}{a cell-by-gene matrix of unspliced counts.}

\item{s.obs}{a cell-by-gene matrix of spliced counts.}

\item{shared_time}{a vector of (point estimate) of gene-shared latent time.}
}
\value{
a plot with PC2 against PC1, and cells are colored by normalized shared latent time.
}
\description{
This function applies PCA on concatenated matrices for observed spliced and unspliced counts, and then visualize
the point estimate of shared latent time after normalizing it on \eqn{[0,1]}.
}
\examples{
shared_time_pca(u.obs = u.obs, s.obs = s.obs, shared_time = shared_time_median)
}
