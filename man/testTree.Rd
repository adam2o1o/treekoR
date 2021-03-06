% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyseTree.R
\name{testTree}
\alias{testTree}
\title{Title}
\usage{
testTree(
  phylo,
  clusters,
  samples,
  classes,
  pos_class_name = NULL,
  subjects = NULL,
  paired = FALSE
)
}
\arguments{
\item{phylo}{a ggtree object}

\item{clusters}{a vector representing the cell type or cluster of each cell
(can be character or numeric). If numeric, cluster names need to be consecutive
starting from 1.}

\item{samples}{a vector identifying the patient each cell belongs to}

\item{classes}{a vector containing the patient outcome/class each cell belongs to}

\item{pos_class_name}{a character indicating which class is positive}

\item{subjects}{a vector containing which subject the cell belongs to, used
to identify matched samples in paired t-tests (not yet tested)}

\item{paired}{a boolean indicating whether to performed paired t-tests
(not yet tested)}
}
\value{
a ggtree object with significance testing results in embedded data
}
\description{
This function takes a hierarchical tree of the cluster
medians of a cytometry dataset, and then uses this structure to perform
t-tests between conditions of patients testing for difference using the
proportion of cluster relative to sample's n and proportion of cluster
relative to sample's n of hierarchical parent cluster.
Takes a ggtree object and returns a ggtree object with testing results
appended in the data
}
\examples{
library(SingleCellExperiment)
data(COVIDSampleData)

sce <- DeBiasi_COVID_CD8_samp
exprs <- t(assay(sce, "exprs"))
clusters <- colData(sce)$cluster_id
classes <- colData(sce)$condition
samples <- colData(sce)$sample_id

clust_tree <- getClusterTree(exprs,
                             clusters,
                             hierarchy_method="hopach")

tested_tree <- testTree(clust_tree$clust_tree,
                        clusters=clusters,
                        samples=samples,
                        classes=classes,
                        pos_class_name=NULL,
                        subjects=NULL,
                        paired = FALSE)
}
