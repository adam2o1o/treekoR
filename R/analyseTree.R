#' Title
#' This function takes a CATALYST sce with clusters and creates a hierarchical tree
#'
#' @return
#'
#' @examples
#' @import data.table
#' @importFrom ggtree ggtree
#' @importFrom ape as.phylo
#'
#' @export
getClusterTree <- function(exprs,
                           clusters,
                           hierarchy_method="average") {
  clust_med_dt <- as.data.table(exprs)
  clust_med_dt[, cluster_id := clusters]
  # data table containing median
  res <-clust_med_dt[, lapply(.SD, median, na.rm=TRUE), by=cluster_id]
  res2 <- res[,.SD, .SDcols = !c('cluster_id')]
  rownames(res2) <- res[["cluster_id"]]

  clust_dist <- dist(scale(res2))
  hc_dend <- hclust(clust_dist, method=hierarchy_method)
  hc_phylo <- as.phylo(hc_dend)
  hc_phylo$tip.label <- as.character(res[["cluster_id"]])

  return(list(
    median_freq = res2,
    clust_tree = hc_phylo
  ))
}

#' Title
#'
#' @param tree
#'
#' @return
#' @export
#'
#' @examples
findChildren <- function(tree) {
  d <- tree$data
  d$clusters <- d$label
  d$clusters <- as.list(as.character(d$clusters))
  uNodes <- sort(unique(d$x), decreasing = TRUE)
  for(x in uNodes){
    nodes = as.matrix(d[which(d$x==x), "node"])
    for(n in nodes){
      parentNode <- as.character(d$parent[which(d$node == n)])
      childClusters <- d$clusters[[which(d$node == n)]]
      parentNodeInd <- which(d$node == parentNode)
      if (is.na(d$clusters[parentNodeInd])) {
        d$clusters[[parentNodeInd]] <- childClusters
      } else {
        d$clusters[[parentNodeInd]] <- c(d$clusters[[parentNodeInd]], childClusters)
      }
    }
  }
  d$clusters <- sapply(d$clusters, unlist)
  d$clusters <- lapply(d$clusters , unique)
  d$clusters <- lapply(d$clusters , function(x) x[!is.na(x)])

  tree$data <- d
  return(tree)
}

#' Title
#' @description This function takes a hierarchical tree of the cluster
#' medians of a cytometry dataset, and then uses this structure to perform
#' t-tests between conditions of patients testing for difference using the
#' proportion of cluster relative to sample's n and proportion of cluster
#' relative to sample's n of hierarchical parent cluster.
#' Takes a ggtree object and returns a ggtree object with testing results
#' appended in the data
#'
#' @param tree a ggtree object
#' @param cells a dataframe containing the clusters for each cell, the
#' sample id, the subject id (needed for paired tests), and the group
#' that the subject/sample belongs to
#' @param groupA
#' @param groupB
#' @param paired
#'
#' @return
#' @export
#'
#' @examples
testTree <- function(phylo,
                     exprs,
                     clusters,
                     samples,
                     classes,
                     pos_class_name=NULL,
                     subjects=NULL,
                     paired = FALSE){
  t <- ggtree(phylo, branch.length="none")
  t <- findChildren(t)
  td <- t$data

  if(paired == TRUE){
    samp2Group <- unique(data.frame(subjects, samples, classes))
    samp2Group <- samp2Group[samp2Group$group%in%c(groupA,groupB),]
    samp2Group <- samp2Group[samp2Group$subject%in%names(which(table(samp2Group$subject)==2)),]
    ### Put catch error in here
    groupA <- as.character(samp2Group[samp2Group$group==groupA,'sample'])
    groupB <- as.character(samp2Group[samp2Group$group==groupB,'sample'])

  }else{
    if (is.null(pos_class_name)) {
      pos_class_name <- unique(classes)[1]
      neg_class_name <- unique(classes)[2]
    } else {
      neg_class_name <- unique(classes)[unique(classes) != pos_class_name]
    }
    samp2Group <- unique(data.frame(samples, classes))
    groupA <- as.character(samp2Group[samp2Group$classes==pos_class_name,'samples'])
    groupB <- as.character(samp2Group[samp2Group$classes==neg_class_name,'samples'])
  }
  statParent <- statAll <- pvalParent <- pvalAll <- NULL

  for( i in 1:nrow(td)){
    childClusters <- td[i,'clusters'][[1]][[1]]
    parent <- td[i,'parent'][[1]]
    parentClusters <- td[td$node == parent,'clusters'][[1]][[1]]

    if(mean(parentClusters %in% childClusters)==1){
      statAll[i] <- 0
      pvalAll[i] <- 1
      statParent[i] <- 0
      pvalParent[i] <- 1
    } else {
      numAll <- table(samples)
      numParent <- tapply(clusters %in% parentClusters, samples , sum)
      numNode <- tapply(clusters %in% childClusters, samples, sum)

      ratioAll <- numNode/pmax(numAll[names(numNode)],1)
      ratioParent <- numNode/pmax(numParent[names(numNode)],1)

      if(paired == TRUE){
        testAll <- t.test(ratioAll[groupA],ratioAll[groupB],paired = TRUE)
        statAll[i] <- testAll$statistic
        pvalAll[i] <- testAll$p.value


        testParent <- t.test(ratioParent[groupA],ratioParent[groupB],paired = TRUE)
        statParent[i] <- testParent$statistic
        pvalParent[i] <- testParent$p.value
      } else {
        testAll <- t.test(ratioAll[groupA],ratioAll[groupB])
        statAll[i] <- testAll$statistic
        pvalAll[i] <- testAll$p.value


        testParent <- t.test(ratioParent[groupA],ratioParent[groupB])
        statParent[i] <- testParent$statistic
        pvalParent[i] <- testParent$p.value
      }

    }
  }

  td$statAll <- statAll
  td$statParent <- statParent
  td$pvalAll <- pvalAll
  td$pvalParent <- pvalParent

  t$data <- td
  return(t)

}
