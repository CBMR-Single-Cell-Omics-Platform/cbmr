% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{volcanoplot}
\alias{volcanoplot}
\title{Volcanoplot}
\usage{
volcanoplot(
  table = NULL,
  logfc_cutoff = NULL,
  fdr_cutoff = NULL,
  extra_pval = NULL,
  extra_logfc = NULL,
  name_col = NULL,
  interactive = FALSE
)
}
\arguments{
\item{table}{data.frame with output from edgeR or limma. Must contain three
columns named "logFC", "PValue"/"P.Value" & "FDR"/"adj.P.Val"}

\item{logfc_cutoff}{numeric, absolute logFC must be higher than this to be
coloured significant}

\item{fdr_cutoff}{numeric, FDR/adj.P.Val must be less than this to be
coloured significant. Uses \link[cbmr]{find_pvalue_cutoff} to determine where
to place cutoff line.}

\item{extra_pval}{numeric, P-value to include when calculating y-limits.
Useful for setting the same limits in multiple volcanoplots}

\item{extra_logfc}{numeric, logFC to include when calculating x-limits.
Useful for setting the same limits in multiple volcanoplots}

\item{name_col}{string, column that contains the gene name, only relevant for
interactive plots}

\item{interactive}{logical, is the plot to be printed in an interactive plot?}
}
\description{
Volcanoplot
}
