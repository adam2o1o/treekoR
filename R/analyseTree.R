#' runHOPACH
#'
#' @param data dataframe containing the median expression of the clusters/cell
#' types
#' @param K  positive integer specifying the maximum number of levels in the
#' tree. Must be 15 or less, due to computational limitations (overflow)
#' @param kmax integer between 1 and 9 specifying the maximum number of children
#' at each node in the tree
#' @param dissimilarity_metric metric used to calculate dissimilarities
#' between clusters/cell types
#'
#' @import hopach
#' @return a list containing the groups each cluster belongs to at each level of
#' the hopach tree
#'
#' @examples
#' library(SingleCellExperiment)
#' library(data.table)
#' data(COVIDSampleData)
#'
#' sce <- DeBiasi_COVID_CD8_samp
#' exprs <- t(assay(sce, "exprs"))
#' clusters <- colData(sce)$cluster_id
#' classes <- colData(sce)$condition
#' samples <- colData(sce)$sample_id
#'
#' clust_med_dt <- as.data.table(exprs)
#' clust_med_dt[, cluster_id := clusters]
#' res <- clust_med_dt[, lapply(.SD, median, na.rm=TRUE), by=cluster_id]
#' res2 <- res[,.SD, .SDcols = !c('cluster_id')]
#'
#' hopach_res <- runHOPACH(as.data.frame(scale(res2)))
#' @export
runHOPACH <- function(data, K = 10, kmax = 5, dissimilarity_metric = "cor") {
    # Compute the distance matrix using correlation
    dist <- distancematrix(data, d = dissimilarity_metric)
    # Run HOPACH function
    clustresult <- hopach(data, K = K , dmat = dist, kmax = kmax)
    # Rearrange HOPACH results
    final_labels <- strsplit(as.character(format(clustresult$final$labels,
                                                 scientific = FALSE)), "")
    # Get each level seperate
    level_list <- list()
    for (i in seq_len(nchar(clustresult$final$labels[1]))) {
        if (i != 1) {
            level_list[[i]] <- paste(level_list[[i - 1]],
                                     unlist(lapply(final_labels, "[[", i)),
                                     sep = "")
        } else{
            level_list[[i]] <- unlist(lapply(final_labels, "[[", i))
        }
    }
    # Generate cutree results
    cutree_list <- list()
    # First level (All ones)
    cutree_list[[1]] <- rep(1, length(rownames(data)))
    names(cutree_list[[1]]) <- rownames(data)
    cutree_list[2:(length(level_list) + 1)] <- lapply(level_list, function(x) {
        x <- as.numeric(as.factor(x))
        names(x) <- rownames(data)
        return(x)
    })
    # Last levels will be all distinct numbers
    cutree_list[[length(cutree_list) + 1]] <-
        seq_len(length(clustresult$final$labels))
    names(cutree_list[[length(cutree_list)]]) <- rownames(data)

    return(list(cutree_list = cutree_list))
}

#' hopachToPhylo
#'
#' @param res an object returned from the runHOPACH() function
#'
#' @importFrom ape read.tree
#' @return a phylogram converted from the outputted list from the runHOPACH
#' function
#' @examples
#' library(SingleCellExperiment)
#' library(data.table)
#' data(COVIDSampleData)
#'
#' sce <- DeBiasi_COVID_CD8_samp
#' exprs <- t(assay(sce, "exprs"))
#' clusters <- colData(sce)$cluster_id
#' classes <- colData(sce)$condition
#' samples <- colData(sce)$sample_id
#'
#' clust_med_dt <- as.data.table(exprs)
#' clust_med_dt[, cluster_id := clusters]
#' res <- clust_med_dt[, lapply(.SD, median, na.rm=TRUE), by=cluster_id]
#' res2 <- res[,.SD, .SDcols = !c('cluster_id')]
#'
#' hopach_res <- runHOPACH(as.data.frame(scale(res2)))
#' phylo <- hopachToPhylo(hopach_res)
#'
#' @export
hopachToPhylo <- function(res) {
    cutree_list_df <- do.call(cbind, res$cutree_list)
    for (i in seq_len(ncol(cutree_list_df)-1)) {
        cutree_list_df[,i] <- paste0(i,"_",cutree_list_df[,i])
    }

    # We add < & > to each cluster for string matching purposes
    base_nodes <- paste0("<",as.vector(cutree_list_df[,ncol(cutree_list_df)]),">")
    # A vector containing the current grouping of nodes as the loop iterates
    current_nodes <- base_nodes

    # This loop is to iterate through the levels of the tree and
    # sequentially combine clusters in a format to be read into a
    # phylogenetic tree
    # ------------------------------------------------------------
    for (i in rev(seq_len(ncol(cutree_list_df)-1))) {
        cur_lev_nodes <- as.vector(cutree_list_df[,i])
        cur_lev_n_t <- table(cur_lev_nodes)
    # Get parent nodes with more than two children
        unique_nodes <- names(cur_lev_n_t[which(cur_lev_n_t >= 2)])
        # For each parent node, we combine the corresponding children nodes
        for (node in unique_nodes) {
            base_nodes_combo <- base_nodes[which(cur_lev_nodes== node)]
            base_nodes_rgx_s <- paste(base_nodes_combo, collapse="|")
            if (sum(grepl(base_nodes_rgx_s, current_nodes)) > 1) {
                nodes_to_combine <- current_nodes[grepl(base_nodes_rgx_s, current_nodes)]
                current_nodes <- c(current_nodes[!current_nodes %in% nodes_to_combine],
                                   paste0("(", paste(nodes_to_combine, collapse=","), ")")
                )
            }
        }
    }
    tree_string <- paste0(gsub("<|>", "", current_nodes), ";")
    # Conver tree string into a phylo object
    tree <- read.tree(text = tree_string)

    return(tree)
}


#' getClusterTree
#' This function takes a CATALYST sce with clusters and creates a hierarchical tree
#'
#' @param exprs a dataframe containing single cell expression data
#' @param clusters a vector representing the cell type or cluster of each cell
#' (can be character or numeric). If numeric, cluster names need to be consecutive
#' starting from 1.
#' @param hierarchy_method a string indicating the hierarchical tree construction
#'  method to be used
#' @param hopach_kmax integer between 1 and 9 specifying the maximum number of
#' children at each node in the tree
#' @param hopach_K positive integer specifying the maximum number of levels in the
#' tree. Must be 15 or less, due to computational limitations (overflow)
#' @param scale_exprs boolean indicating whether to scale median cluster expression
#' data before constructing hierarchical tree
#'
#' @import data.table
#' @importFrom ggtree ggtree
#' @importFrom ape as.phylo
#' @importFrom stats dist hclust median
#'
#' @return a list containing the cluster median frequencies and a phylogram of the
#' hierarchical tree
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' data(COVIDSampleData)
#'
#' sce <- DeBiasi_COVID_CD8_samp
#' exprs <- t(assay(sce, "exprs"))
#' clusters <- colData(sce)$cluster_id
#' classes <- colData(sce)$condition
#' samples <- colData(sce)$sample_id
#'
#' clust_tree <- getClusterTree(exprs,
#'                              clusters,
#'                              hierarchy_method="hopach")
getClusterTree <- function(exprs,
                           clusters,
                           hierarchy_method="hopach",
                           hopach_kmax = 5,
                           hopach_K = 10,
                           scale_exprs=TRUE) {
    clust_med_dt <- as.data.table(exprs)
    clust_med_dt[, cluster_id := clusters]
    # data table containing median
    res <- clust_med_dt[, lapply(.SD, median, na.rm=TRUE), by=cluster_id]
    res2 <- res[,.SD, .SDcols = !c('cluster_id')]
    rownames(res2) <- res[["cluster_id"]]

    if (scale_exprs) {
        res_unscaled <- res2
        res2[, (colnames(res2)) := lapply(.SD, scale), .SDcols=colnames(res2)]
    } else {
        res_unscaled <- res2
    }

    if (hierarchy_method == "hopach") {
        hp_dend <- runHOPACH(data = as.data.frame(res2),
                             kmax=hopach_kmax,
                             K=hopach_K)
        hc_phylo <- hopachToPhylo(hp_dend)
        hc_phylo$tip.label <- rownames(res2)[as.numeric(hc_phylo$tip.label)]
    } else {
        clust_dist <- dist(res2)
        hc_dend <- hclust(clust_dist, method=hierarchy_method)
        hc_phylo <- as.phylo(hc_dend)
        hc_phylo$tip.label <- as.character(res[["cluster_id"]])
    }

    return(list(
        median_freq = res_unscaled,
        clust_tree = hc_phylo
    ))
}

#' findChildren
#'
#' @param tree a ggtree object
#'
#' @return a ggtree object with the data containing a column with the clusters
#' contained in each node
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
    d$clusters <- lapply(d$clusters, unlist)
    d$clusters <- lapply(d$clusters, unique)
    d$clusters <- lapply(d$clusters, function(x) x[!is.na(x)])

    tree$data <- d
    return(tree)
}

#' getTotalProp
#' @description getCellProp helper function
#'
#' @param vars1 name of cell type, matching to column in n_cells
#' @param n_cells matrix of counts of each cell type per sample
#' @param n_cells_pat vector containing number of cells per sample
#'
#' @return a vector containing the proportions of cell type vars1 as a
#' percent of total per sample
getTotalProp <- function(vars1,
                         n_cells,
                         n_cells_pat) {
    rowSums(n_cells[,vars1])/n_cells_pat
}

#' getParentProp
#' @description getCellProp helper function
#'
#' @param vars1 name of cell type, matching to column in n_cells
#' @param vars2 name of parent cell type, matching to column in n_cells
#' @param n_cells matrix of counts of each cell type per sample
#'
#' @return a vector containing the proportions of cell type vars1 as a
#' percent of parent vars2 per sample
getParentProp <- function(vars1,
                          vars2,
                          n_cells) {
    # if parent population is null, then return NA for each sample
    if (is.null(vars2)) {return(rep(NA,length=nrow(n_cells)))}
    rowSums(n_cells[,vars1])/rowSums(n_cells[,vars2])
}

#' getCellProp
#'
#' @param phylo a phylogram with tip.labels corresponding to cell types/cluster
#' contained in 'clusters' vector
#' @param clusters a vector representing the cell type or cluster of each
#' cell (can be character or numeric). If numeric, cluster names need to be
#' consecutive starting from 1.
#' @param samples a vector identifying the patient each cell belongs to
#' @param classes a vector containing the patient outcome/class each cell belongs to
#' @param excl_top_node_parent a boolean indicating whether the %parent should be calculated
#' for cell types with the highest node as their parent
#'
#' @importFrom dplyr %>% tally group_by select distinct pull mutate n
#' @importFrom tidyr complete pivot_wider
#' @return a dataframe containing proportions calculated for each sample
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' data(COVIDSampleData)
#'
#' sce <- DeBiasi_COVID_CD8_samp
#' exprs <- t(assay(sce, "exprs"))
#' clusters <- colData(sce)$cluster_id
#' classes <- colData(sce)$condition
#' samples <- colData(sce)$sample_id
#'
#' clust_tree <- getClusterTree(exprs,
#'                              clusters,
#'                              hierarchy_method="hopach")
#'
#' prop_df <- getCellProp(clust_tree$clust_tree,
#'                         clusters=clusters,
#'                         samples=samples,
#'                         classes=classes)
getCellProp <- function(phylo,
                        clusters,
                        samples,
                        classes,
                        excl_top_node_parent=TRUE){
    t <- findChildren(ggtree(phylo, branch.length="none"))
    td <- t$data
    td$label <- ifelse(is.na(td$label),td$node, td$label)
    # Get all parent nodes
    parent_node_ind <- match(td$parent, td$node)
    child_node_ind <- td$node
    if (excl_top_node_parent) {
        iter_layer_ind <- seq_len(max(td$x)-1)
        for (lvl in iter_layer_ind) {
            td[,paste0("parent_",lvl,"_clusters")] <- list(td$clusters[parent_node_ind])
            # Don't compute proportion relative to  all cells (top node)
            td[td$parent[parent_node_ind] == parent_node_ind,paste0("parent_",lvl,"_clusters")] <- NA
            parent_node_ind <- td$parent[parent_node_ind]
        }
    } else {
        iter_layer_ind <- seq_len(max(td$x))
        for (lvl in iter_layer_ind) {
            td[,paste0("parent_",lvl,"_clusters")] <- list(td$clusters[parent_node_ind])
            td[which(child_node_ind == parent_node_ind),paste0("parent_",lvl,"_clusters")] <- NA
            child_node_ind <- parent_node_ind
            parent_node_ind <- td$parent[parent_node_ind]
        }
    }
    # Remove root node
    td <- td[-which(td$node == td$parent),]
    n_cells <- data.frame(cluster_id=factor(clusters,levels=t$data$label),
                          sample_id=samples) %>%
        group_by(cluster_id, sample_id, .drop=FALSE) %>%
        tally %>%
        tidyr::complete() %>%
        tidyr::pivot_wider(names_from=cluster_id, values_from=n, values_fill=0)

    n_cells_pat <- rowSums(n_cells %>%
        select(-sample_id))
    # Get absolute proportions of aggregated cells
    prop_total <- vapply(td %>% pull(clusters),
                         getTotalProp,
                         double(length(n_cells_pat)),
                         n_cells=n_cells,
                         n_cells_pat=n_cells_pat)
    colnames(prop_total) <- paste0("perc_total_", td$label)
    # Get parent proportions of all cells
    prop_parent <- lapply(iter_layer_ind,
                          function(lvl) {
                              # For each cluster, calculate proportion for each sample
                              res <- mapply(getParentProp,
                                            td %>% pull(clusters),
                                            td %>% pull(paste0("parent_",lvl,"_clusters")),
                                            MoreArgs = list(n_cells=n_cells))
                              colnames(res) <- paste0("perc_parent_",lvl,"_",td$label)
                              # Remove columns with all NA
                              res <- res[,colSums(is.na(res)) != nrow(res)]
                              return(res)
                              }) %>%
        do.call(cbind, .)
    # Mapping from sample to class outcome
    samp_class <- dplyr::distinct(data.frame(sample_id=samples, class=classes))
    prop_df <- cbind(data.frame(sample_id = n_cells$sample_id,
                                class = samp_class[,2][match(n_cells$sample_id,
                                                             samp_class[,1])]),
                     prop_total,
                     prop_parent)
    return(prop_df)
}

#' geometricMean
#' @description getCellGMeans helper function
#'
#' @param x vector containing numeric values
#' @param na.rm whether or not to ignore NA values
#' @return geomtric mean of vector x
geometricMean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


#' getCellGMeans
#'
#' @param phylo a phylogram with tip.labels corresponding to cell types/cluster
#' contained in 'clusters' vector
#' @param exprs a dataframe containing single cell expression data
#' @param clusters a vector representing the cell type or cluster of each
#' cell (can be character or numeric). If numeric, cluster names need to be
#' consecutive starting from 1.
#' @param samples a vector identifying the patient each cell belongs to
#' @param classes a vector containing the patient outcome/class each cell belongs to
#'
#' @importFrom dplyr %>% tally group_by select distinct pull mutate summarise_all n
#' @importFrom tidyr  pivot_wider
#' @return a dataframe containing proportions calculated for each sample
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' data(COVIDSampleData)
#'
#' sce <- DeBiasi_COVID_CD8_samp
#' exprs <- t(assay(sce, "exprs"))
#' clusters <- colData(sce)$cluster_id
#' classes <- colData(sce)$condition
#' samples <- colData(sce)$sample_id
#'
#' clust_tree <- getClusterTree(exprs,
#'                              clusters,
#'                              hierarchy_method="hopach")
#'
#' means_df <- getCellGMeans(clust_tree$clust_tree,
#'                         exprs=exprs,
#'                         clusters=clusters,
#'                         samples=samples,
#'                         classes=classes)
getCellGMeans <- function(phylo,
                          exprs,
                          clusters,
                          samples,
                          classes){
    t <- findChildren(ggtree(phylo, branch.length="none"))
    td <- t$data
    td$label <- ifelse(is.na(td$label),td$node, td$label)
    td$parentClusters <- td$clusters[match(td$parent, td$node)]
    # Remove root node
    td <- td[-which(td$node == td$parent),]

    # Calculate gmeans of base nodes
    mean_all <- as.data.frame(exprs) %>%
        mutate(
            cluster_id=factor(clusters,levels=t$data$label),
            sample_id=samples) %>%
        group_by(cluster_id, sample_id, .drop=FALSE) %>%
        summarise_all(geometricMean) %>%
        tidyr::pivot_wider(names_from=cluster_id, names_prefix="gmean_",
                           values_from=colnames(exprs))
    # Calculate number of cells per base node and sample
    n_all <- as.data.frame(exprs) %>%
        mutate(
            cluster_id=factor(clusters,levels=t$data$label),
            sample_id=samples
            ) %>%
        group_by(cluster_id, sample_id, .drop=FALSE) %>%
        tally %>%
        # tidyr::complete() %>%
        tidyr::pivot_wider(names_from=cluster_id, names_prefix="n_",
                           values_from=n, values_fill=0)
    # Calculate gmeans of each non-leaf cluster in tree
    mean_all_agg <- lapply(
        which(!td$isTip),
        function(ind) {
            x <- td$clusters[[ind]]
            mat2 <- n_all[,paste0("n_",x)]
            gmean_mat <- vapply(colnames(exprs),
                                function(y) {
                                    mat1 <- mean_all[,paste0(y,"_gmean_",x)]
                                    res <- exp(rowSums(mat2*log(mat1), na.rm=TRUE)/rowSums(mat2))
                                    return(res)
                                    },
                                double(nrow(mat2)))
            colnames(gmean_mat) <- paste0(colnames(exprs),"_gmean_",
                                          td$label[[ind]])
            return(gmean_mat)
            })

    samp_class <- dplyr::distinct(data.frame(sample_id=as.character(samples),
                                             class_id=as.character(classes)))

    mean_all_df <- cbind(
        data.frame(
            sample_id = mean_all$sample_id,
            class = samp_class[,2][match(mean_all$sample_id,
                                         samp_class[,1])]),
        # Here exclude duplicated means calculated in the initial group_by sample_id/cluster_id
        # because some parent nodes could be contained in the supplied clusters,
        # but we want to get the aggregated mean of cells belonging to that parent
        # as defiend by the hierarchical tree supplied
        mean_all[,-which(colnames(mean_all)=="sample_id" |
                             colnames(mean_all) %in%
                             colnames(cbind(mean_all[,"sample_id"], mean_all_agg)))],
        mean_all_agg)

    return(mean_all_df)
}

#' runEdgeRTests
#' @description This function runs edgeR using the treekoR inputs across all nodes
#' of the hierarchical tree of clusters, adapted from the diffcyt workflow
#'
#' @param td a dataframe of data from ggtree object
#' @param clusters a vector representing the cell type or cluster of each cell
#' (can be character or numeric). If numeric, cluster names need to be consecutive
#' starting from 1.
#' @param classes a vector containing the patient outcome/class each cell belongs to
#' @param samples a vector identifying the patient each cell belongs to
#'
#' @importFrom edgeR DGEList estimateDisp glmFit glmLRT topTags
#' @importFrom dplyr %>% distinct mutate left_join rename
#' @importFrom stats model.matrix
#'
#' @return a dataframe with pvalues, test statistic (signed -log10(p)), and FDR
runEdgeRTests <- function(td,
                          clusters,
                          samples,
                          classes,
                          pos_class_name) {
    n_parent <- NULL
    n_node <- NULL
    for (i in seq_len(nrow(td))) {
        child_clusters <- td[i, "clusters"][[1]][[1]]
        parent_node <- td[i, "parent"][[1]]
        parent_clusters <- td[td$node == parent_node, "clusters"][[1]][[1]]
        n_parent <- rbind(n_parent, tapply(clusters %in% parent_clusters,
                                           samples, sum))
        n_node <- rbind(n_node, tapply(clusters %in% child_clusters, samples,
                                       sum))
    }
    class_df <- data.frame(classes, samples) %>%
        distinct() %>%
        mutate(cl = ifelse(classes == pos_class_name,1, 0))

    ## EdgeR parent
    cl <- class_df$cl
    names(cl) <- class_df$samples
    y <- edgeR::DGEList(counts=n_node,
                        group = cl[colnames(n_node)],
                        lib.size = rep(1,ncol(n_node)))
    design <- model.matrix(~cl[colnames(n_node)])
    # log because of getOffset()
    y$offset <- log(n_parent+1)
    y <- edgeR::estimateDisp(y, design)
    fit <- edgeR::glmFit(y,design)
    lrt <- edgeR::glmLRT(fit,coef=2)
    top_parent <- edgeR::topTags(lrt, n = 10000)
    top_nodes_parent <- top_parent$table
    top_nodes_parent$node <- rownames(top_nodes_parent)

    ## EdgeR total
    y <- edgeR::DGEList(counts=n_node, group = cl[colnames(n_node)])
    y <- edgeR::estimateDisp(y)
    fit <- edgeR::glmFit(y)
    lrt <- edgeR::glmLRT(fit,coef=2)
    top_total <- edgeR::topTags(lrt, n = 10000)
    top_nodes_total <- top_total$table
    colnames(top_nodes_total) <- paste(colnames(top_nodes_total), "total", sep = "_")
    top_nodes_total$node <- rownames(top_nodes_total)

    td <- td %>%
        left_join(
            top_nodes_parent %>%
                mutate(stat_parent = -log(PValue,10)*sign(logFC),
                       node = as.numeric(node)) %>%
                rename(c("pval_parent"="PValue", "FDR_parent"="FDR")) %>%
                dplyr::select(stat_parent, pval_parent, FDR_parent, node),
            by=c("node"="node")
        ) %>%
        left_join(
            top_nodes_total %>%
                mutate(stat_total = -log(PValue_total,10)*sign(logFC_total),
                       node = as.numeric(node)) %>%
                rename(c("pval_total"="PValue_total")) %>%
                dplyr::select(stat_total, pval_total, FDR_total, node),
            by=c("node"="node")
        )

    return(td)
}

#' runGLMMTests
#' @description This function runs GLMM using the treekoR inputs across all nodes
#' of the hierarchical tree of clusters, adapted from the diffcyt workflow.
#' (Unable to get direction of test statistic currently)
#'
#' @param td a dataframe of data from ggtree object
#' @param clusters a vector representing the cell type or cluster of each cell
#' (can be character or numeric). If numeric, cluster names need to be consecutive
#' starting from 1.
#' @param classes a vector containing the patient outcome/class each cell belongs to
#' @param samples a vector identifying the patient each cell belongs to
#'
#' @importFrom diffcyt createFormula
#' @importFrom lme4 glmer
#' @importFrom multcomp glht
#' @importFrom dplyr %>% bind_cols distinct mutate left_join rename filter pull
#'
#' @return a dataframe with pvalues and test statistics
runGLMMTests <- function(td,
                         clusters,
                         samples,
                         classes,
                         pos_class_name,
                         neg_class_name) {
    n_parent <- NULL
    n_node <- NULL
    for (i in seq_len(nrow(td))) {
        child_clusters <- td[i, "clusters"][[1]][[1]]
        parent_node <- td[i, "parent"][[1]]
        parent_clusters <- td[td$node == parent_node, "clusters"][[1]][[1]]
        n_parent <- rbind(n_parent, tapply(clusters %in% parent_clusters,
                                           samples, sum))
        n_node <- rbind(n_node, tapply(clusters %in% child_clusters, samples,
                                       sum))
    }

    glmm_p_vals <- data.frame(matrix(nrow=nrow(n_node), ncol=4))
    colnames(glmm_p_vals) <- c("stat_total", "pval_total",
                               "stat_parent", "pval_parent")

    n_cells_smp <- colSums(n_node[td %>% filter(isTip) %>% pull(node),])

    experiment_info <- data.frame(sample_id=samples,
                                  condition=classes) %>%
        distinct() %>%
        left_join(data.frame(n_cells=n_cells_smp,
                             sample_id=names(n_cells_smp)),
                  by="sample_id")
    contrast <- matrix(c(0, 1), ncol=2)
    formula <- diffcyt::createFormula(experiment_info,
                                      cols_fixed="condition",
                                      cols_random = "sample_id")
    levels(formula$data$condition) <- c(pos_class_name, neg_class_name)

    for (i in seq_len(nrow(n_node))) {
        tryCatch({
            # data for cluster i
            # note: divide by total number of cells per sample (after filtering) to get
            # proportions instead of counts
            y <- n_node[i, ] / n_cells_smp
            data_i <- cbind(y, n_cells_smp, formula$data)
            # fit model
            # note: provide proportions (y) together with weights for total number of cells per
            # sample (n_cells_smp); this is equivalent to providing counts
            fit <- lme4::glmer(formula$formula, data = data_i, family = "binomial", weights = n_cells_smp)
            # test contrast
            test <- multcomp::glht(fit, contrast)
            # return p-value
            glmm_p_vals[i,"stat_total"] <- summary(test)$test$tstat
            glmm_p_vals[i,"pval_total"] <- summary(test)$test$pvalues
            # return NA as p-value if there is an error
        }, error = function(e) NA)
        tryCatch({
            # data for cluster i
            # note: divide by total number of cells per sample (after filtering) to get
            # proportions instead of counts
            y <- n_node[i, ] / n_parent[i,]
            data_i <- cbind(y, n_cells_smp, formula$data)
            # fit model
            # note: provide proportions (y) together with weights for total number of cells per
            # sample (n_cells_smp); this is equivalent to providing counts
            fit <- lme4::glmer(formula$formula, data = data_i, family = "binomial", weights = n_parent[i,])
            # test contrast
            test <- multcomp::glht(fit, contrast)
            # return p-value
            glmm_p_vals[i,"stat_parent"] <- summary(test)$test$tstat
            glmm_p_vals[i,"pval_parent"] <- summary(test)$test$pvalues
            # return NA as p-value if there is an error
        }, error = function(e) NA)
    }
    td <- td %>%
        bind_cols(glmm_p_vals)
    return(td)
}

#' testTree
#' @description This function takes a hierarchical tree of the cluster
#' medians of a cytometry dataset, and then uses this structure to perform
#' t-tests between conditions of patients testing for difference using the
#' proportion of cluster relative to sample's n and proportion of cluster
#' relative to sample's n of hierarchical parent cluster.
#' Takes a ggtree object and returns a ggtree object with testing results
#' appended in the data
#'
#' @param phylo a ggtree object
#' @param clusters a vector representing the cell type or cluster of each cell
#' (can be character or numeric). If numeric, cluster names need to be consecutive
#' starting from 1.
#' @param classes a vector containing the patient outcome/class each cell belongs to
#' @param samples a vector identifying the patient each cell belongs to
#' @param pos_class_name a character indicating which class is positive
#' @param sig_test a character, either "ttest" or "wilcox" indicating the significance test to be used
#'
#' @importFrom stats t.test wilcox.test as.formula
#'
#' @return a ggtree object with significance testing results in embedded data
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' data(COVIDSampleData)
#'
#' sce <- DeBiasi_COVID_CD8_samp
#' exprs <- t(assay(sce, "exprs"))
#' clusters <- colData(sce)$cluster_id
#' classes <- colData(sce)$condition
#' samples <- colData(sce)$sample_id
#'
#' clust_tree <- getClusterTree(exprs,
#'                              clusters,
#'                              hierarchy_method="hopach")
#'
#' tested_tree <- testTree(clust_tree$clust_tree,
#'                         clusters=clusters,
#'                         samples=samples,
#'                         classes=classes,
#'                         sig_test="ttest",
#'                         pos_class_name=NULL)
testTree <- function(phylo,
                     clusters,
                     samples,
                     classes,
                     sig_test="ttest",
                     pos_class_name=NULL){
    if (!is.null(pos_class_name)) {
        if (!pos_class_name %in% unique(classes)) {
            stop("'pos_class_name' needs to be in classes")
        }
    }
    if (length(unique(classes)) > 2) {
        stop("length(unique(classes)) > 2.
             treekoR can currently only test between two classes.")
    }
    if (!(sig_test %in% c("ttest", "wilcox", "edgeR", "GLMM"))) {
        stop("Significance test needs to be either 'ttest', 'wilcox', 'edgeR' or 'GLMM'")
    }

    t <- findChildren(ggtree(phylo, branch.length="none"))
    td <- t$data
    td$label <- ifelse(is.na(td$label),td$node, td$label)
    if (is.null(pos_class_name)) {
        pos_class_name <- names(table(classes))[1]
        neg_class_name <- names(table(classes))[2]
    } else {
        neg_class_name <- names(table(classes))[names(table(classes)) != pos_class_name]
    }
    if (sig_test == "edgeR") {
        td <- runEdgeRTests(td, clusters, samples, classes, pos_class_name)
        t$data <- td
        return(t)
    }
    if (sig_test == "GLMM") {
        td <- runGLMMTests(td, clusters, samples, classes, pos_class_name, neg_class_name)
        t$data <- td
        return(t)
    }
    prop_df <- getCellProp(phylo,
                           clusters=clusters,
                           samples=samples,
                           classes=classes,
                           excl_top_node_parent=FALSE)
    prop_df$class <- factor(prop_df$class, levels=c(pos_class_name,neg_class_name))
    stat_parent <- stat_total <- pval_parent <- pval_total <- NULL
    for( i in seq_len(nrow(td))[-which(td$node == td$parent)]){
        cell_type = td$label[i]
        if (sig_test=="ttest") {
            test_total <- t.test(formula=as.formula(paste0("`perc_total_",cell_type, "`~class")),
                                 data=prop_df)
            test_parent <- t.test(formula=as.formula(paste0("`perc_parent_1_",cell_type, "`~class")),
                                  data=prop_df)
        }
        else if (sig_test=="wilcox") {
            test_total <- wilcox.test(formula=as.formula(paste0("`perc_total_",cell_type, "`~class")),
                                      data=prop_df)
            test_parent <- wilcox.test(formula=as.formula(paste0("`perc_parent_1_",cell_type, "`~class")),
                                       data=prop_df)
        }
        stat_total[i] <- test_total$statistic
        pval_total[i] <- test_total$p.value
        stat_parent[i] <- test_parent$statistic
        pval_parent[i] <- test_parent$p.value
    }
    td[, c("stat_total", "stat_parent", "pval_total", "pval_parent")] <-
        data.frame(stat_total, stat_parent, pval_total, pval_parent)
    t$data <- td
    return(t)
}

#' getTreeResults
#'
#' @param testedTree a ggtree object outputed from testTree()
#' @param sort_by whether to sort by p-values testing via proportions to parent or
#' p-values testing via absolute proportions. Values can can be c(NA, "parent", "all")
#'
#' @return a dataframe with hierarchical tree nodes, corresponding clusters and
#' corresponding significance testing results
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' data(COVIDSampleData)
#'
#' sce <- DeBiasi_COVID_CD8_samp
#' exprs <- t(assay(sce, "exprs"))
#' clusters <- colData(sce)$cluster_id
#' classes <- colData(sce)$condition
#' samples <- colData(sce)$sample_id
#'
#' clust_tree <- getClusterTree(exprs,
#'                              clusters,
#'                              hierarchy_method="hopach")
#'
#' tested_tree <- testTree(clust_tree$clust_tree,
#'                         clusters=clusters,
#'                         samples=samples,
#'                         classes=classes,
#'                         pos_class_name=NULL)
#'
#' res_df <- getTreeResults(tested_tree)
#'
#' head(res_df, 10)
getTreeResults <- function(testedTree,
                           sort_by = "parent"){
    if (!sort_by %in% c("parent", "total")) {
        stop("'sort_by' must be 'parent' or 'total'")
    }

    res <- as.data.frame(testedTree$data[,c(
        "parent", "node", "isTip",
        "clusters", "stat_total", "stat_parent",
        "pval_total", "pval_parent"
    )])

    if (sort_by == "parent") {
        res <- res[order(res$pval_parent),]
    }

    else if (sort_by == "total") {
        res <- res[order(res$pval_total),]
    }

    return(res)
}


