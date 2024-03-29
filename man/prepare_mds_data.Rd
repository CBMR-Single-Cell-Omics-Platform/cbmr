% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_figures.R
\name{prepare_mds_data}
\alias{prepare_mds_data}
\alias{prepare_mds_data.EList}
\alias{prepare_mds_data.DGEList}
\alias{prepare_mds_data.default}
\title{Prepare data for MDS plot}
\usage{
prepare_mds_data(y, dim_plot, colour_by = NULL, metadata = NULL)

\method{prepare_mds_data}{EList}(y, dim_plot, metadata = NULL)

\method{prepare_mds_data}{DGEList}(y, dim_plot, metadata = NULL)

\method{prepare_mds_data}{default}(y, dim_plot, metadata = NULL)
}
\arguments{
\item{y}{DGEList, EList or matrix.}

\item{dim_plot}{integer vector of length two specifying which principal
components should be plotted.}

\item{metadata}{optional data.frame with additional meta data.}
}
\value{
data.frame with data necessary for an MDS plot
}
\description{
Prepare data for MDS plot
}
\section{Methods (by class)}{
\itemize{
\item \code{prepare_mds_data(EList)}: Prepare data for MDS plot DGEList

\item \code{prepare_mds_data(DGEList)}: Prepare data for MDS plot DGEList

\item \code{prepare_mds_data(default)}: Prepare data for MDS plot default

}}
