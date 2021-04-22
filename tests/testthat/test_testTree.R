suppressPackageStartupMessages({
  library(treekoR)
  library(SingleCellExperiment)
})

data(COVIDSampleData)
sce <- DeBiasi_COVID_CD8_samp

exprs <- t(assay(sce, "exprs"))
clusters <- colData(sce)$cluster_id
classes <- colData(sce)$condition
samples <- colData(sce)$sample_id

clust_tree <- getClusterTree(exprs,
                             clusters,
                             hierarchy_method="hopach")

tested_tree <- testTree(phylo=clust_tree$clust_tree,
                        exprs=exprs,
                        clusters=clusters,
                        samples=samples,
                        classes=classes,
                        pos_class_name=NULL,
                        subjects=NULL,
                        paired = FALSE)


## Randomise clusters to be in string
cluster_names <- paste0("cluster_name", 1:50)
clusters_str <- sample(cluster_names, 5000, replace=T)

## Test framework using character as cluster names
clust_tree_str <- getClusterTree(exprs,
                             clusters_str,
                             hierarchy_method="hopach")

tested_tree_str <- testTree(phylo=clust_tree_str$clust_tree,
                        exprs=exprs,
                        clusters=clusters_str,
                        samples=samples,
                        classes=classes,
                        pos_class_name=NULL,
                        subjects=NULL,
                        paired = FALSE)

test_that("Dataframe of median cluster expression has correct dimensions", {
  expect_equal(nrow(clust_tree$median_freq), length(unique(clusters)))
  expect_equal(ncol(clust_tree$median_freq), ncol(exprs))
})

test_that("getClusterTree is a list with two values", {
  expect_output(str(clust_tree), "List of 2")
})

test_that("Cluster names are carried through pipeline", {
  expect_true(all(unique(unlist(tested_tree_str$data$clusters)) %in% unique(clusters_str)))
  expect_true(all(rownames(clust_tree_str$median_freq) %in% unique(clusters_str)))
  expect_true(all(clust_tree_str$clust_tree$tip.label %in% unique(clusters_str)))

  expect_true(all(unique(clusters_str) %in% unique(unlist(tested_tree_str$data$clusters))))
  expect_true(all(unique(clusters_str) %in% rownames(clust_tree_str$median_freq)))
  expect_true(all(unique(clusters_str) %in% clust_tree_str$clust_tree$tip.label))
})

