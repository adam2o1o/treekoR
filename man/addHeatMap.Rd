% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualiseTree.R
\name{addHeatMap}
\alias{addHeatMap}
\title{Title}
\usage{
addHeatMap(
  p,
  cluster_medians,
  offset = 0.5,
  width = 1,
  expand_y_lim = 20,
  low = "#313695",
  mid = "ivory",
  high = "#A50026",
  colnames_angle = 90,
  metric_name = "Column z-score"
)
}
\arguments{
\item{p}{a phylogenetic tree plot created from the ggtree() function}

\item{cluster_medians}{a dataframe with the cluster medians.
The rownumbers of the clusters median data frame should correspond to the
nodes in the phylo tree. The column names should also correspond to the labels
you want to use}

\item{offset}{the distance between the tree plot and heatmap}

\item{width}{width of each tile in the heatmap}

\item{expand_y_lim}{white space below heatmap}

\item{low}{colour used for low values on heatmap}

\item{mid}{colour used for medium values on heatmap}

\item{high}{colour used for large values on heatmap}

\item{colnames_angle}{angle for x-axis label}

\item{metric_name}{legend title}
}
\value{
an interactive ggplot graph object with heatmap of median cluster
expressions plotted alongside hierarchical tree
}
\description{
a function to add a heatmap of cluster medians alongisde the
phylogenetic tree
}
