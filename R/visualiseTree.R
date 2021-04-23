#' Title
#' @description a function to create a skeleton tree diagram to display significance
#' testing results on each node
#'
#' @param offset distance between leaf nodes on the tree and their labels
#' @param font_size font size of leaf labels
#' @param hjust horizontal justification as defined in ggplot2
#' @param p a phylogenetic tree plot created from the ggtree() function
#'
#' @importFrom ggtree geom_tiplab
#' @return a ggtree graph object with the hierarchical tree of clusters and corresponding labels
addTree <- function(p,
                    offset = 0.3,
                    font_size = 2.5,
                    hjust = 0) {

  p$data$label <- ifelse(is.na(p$data$label),
                         p$data$node,
                         p$data$label)

  # draw graph
  p2 <- p +
    geom_tiplab(size=font_size,
                hjust = hjust,
                offset = offset)

  return(p2)
}

#' Title
#' @description a function to add a heatmap of cluster medians alongisde the phylogenetic tree
#'
#' @param p a phylogenetic tree plot created from the ggtree() function
#' @param cluster_medians a dataframe with the cluster medians.
#' The rownumbers of the clusters median data frame should correspond to the nodes in the phylo tree.
#' The column names should also correspond to the labels you want to use
#' @param offset the distance between the tree plot and heatmap
#' @param width width of each tile in the heatmap
#' @param expand_y_lim white space below heatmap
#' @param low colour used for low values on heatmap
#' @param mid colour used for medium values on heatmap
#' @param high colour used for large values on heatmap
#' @param colnames_angle angle for x-axis label
#' @param metric_name legend title
#'
#' @import ggiraph
#' @import ggplot2
#' @return an interactive ggplot graph object with heatmap of median cluster expressions
#'  plotted alongside hierarchical tree
addHeatMap <- function(p,
                       cluster_medians,
                       offset = 0.5,
                       width = 1,
                       expand_y_lim = 20,
                       low = "#313695", #"purple",
                       mid = "ivory",#"#FFFFBF", #"white",
                       high = "#A50026",
                       colnames_angle = 90,
                       metric_name = "Column z-score") {
  tree_df <- p$data
  # Get x value for start of heatmap
  start_x <- max(tree_df$x, na.rm = TRUE) + offset

  clust_rows <- rownames(cluster_medians)
  cluster_medians <- as.data.frame(scale(cluster_medians))
  rownames(cluster_medians) <- clust_rows

  clust_med_tree_ind <- match(clust_rows, tree_df$label)

  # Create Heatmap Matrix using y values from tree and getting the x values
  heatmap_df <- cluster_medians %>%
    dplyr::mutate(y = tree_df$y[clust_med_tree_ind],
                  label = tree_df$label[clust_med_tree_ind]) %>%
    tidyr::gather(variable, value, -c(y,label)) %>%
    dplyr::mutate(x = start_x + as.numeric(as.factor(.$variable))*width)

  # Insert tooltips
  tile_tooltip <- paste0("<b>Cluster</b>: ", heatmap_df$label,
                         "\n", heatmap_df$variable," = ",
                         signif(heatmap_df$value,2))

  # Add x-axis labels
  heatmap_labels <- heatmap_df %>%
    dplyr::select(variable, x) %>%
    unique() %>%
    dplyr::mutate(y = min(heatmap_df$y) - 2)

  p2 <- p +
    geom_tile_interactive(data = heatmap_df,
                          aes(x=x, y=y,
                              fill = ifelse(abs(value)>3.5, sign(value)*3.5, value), tooltip=tile_tooltip),
                          width = width,
                          inherit.aes = FALSE) +
    scale_fill_gradient2(low = low,
                         mid = mid,
                         high = high,
                         midpoint=0,
                         na.value = NA,
                         name = metric_name, limits=c(-3.5,3.5)) +
    geom_text(data = heatmap_labels,
              aes(x=x, y=y, label=variable),
              angle = 50,
              hjust=1,
              size=3) +
    expand_limits(y = -expand_y_lim)

  return(p2)
}

#' Title
#' @description a function to add the frequency bars for each cluster
#'
#' @param clusters a vector representing the cell type or cluster of each cell (can be character or numeric)
#' @param offset distance between the heatmap and frequency bars
#' @param bar_length length of bar with max frequency
#' @param bar_width width of each frequency bar
#' @param freq_labels boolean indicated whether or not to show frequency bar labels
#' @param p a phylogenetic tree plot created from the ggtree() function
#'
#' @return an interactive ggplot graph object with frequency bars of clusters alongside
#' heatmap of cluster median expression
addFreqBars <- function(p,
                        clusters,
                        offset = 0.75,
                        bar_length = 3,
                        bar_width=0.4,
                        freq_labels = FALSE) {

  tree_df <- p$data

  # Get starting x value by getting all x values plotted in ggplot
  start_x <- max(vapply(ggplot_build(p)$data,
                        function(dat) if(!is.null(dat$x)) max(dat$x) else 0,
                        numeric(1))) + offset

  # Calculate frequency of each cluster
  cluster_freq_df <- data.frame(
    freq = as.numeric(table(clusters)),
    cluster = names(table(clusters))
  )

  clust_med_tree_ind <- match(cluster_freq_df$cluster, tree_df$label)
  clust_med_tree_y <- tree_df$y[clust_med_tree_ind]

  # Add dimensions for bars
  cluster_freq_df <- cluster_freq_df %>%
    dplyr::mutate(y = clust_med_tree_y,
                  ymin = clust_med_tree_y - bar_width,
                  ymax = clust_med_tree_y + bar_width,
                  xmin = start_x,
                  xmax = start_x + bar_length * freq / max(freq)) %>%
    dplyr::mutate(x_label = max(xmax) + 0.5) %>%
    dplyr::mutate(freq_label = paste0(round(100*freq/sum(freq), 2), "% (", freq, ")")) %>%
    dplyr::filter(!is.na(y)) # Remove data for clusters that may be in the clusters factor levels, but not present in exprs

  # Add Tooltips
  bar_tooltip <- paste0("<b>Cluster</b>: ", cluster_freq_df$cluster,
                        '\n % of Cells = ', cluster_freq_df$freq_label)

  p2 <- p +
    geom_rect_interactive(data=cluster_freq_df,
                          mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, tooltip=bar_tooltip),
                          fill='grey',
                          inherit.aes = FALSE)

  if (freq_labels) {
    p2 <- p2 +
      geom_text(data = cluster_freq_df,
                aes(x=x_label, y=y, label= freq_label),
                size=1.75,
                hjust=0)
  }

  return(p2)
}


#' colourTree
#'
#' @param tree a tree plot created from the ggtree() function
#' with p$data containing test statisic and p-
#' @param point_size size of nodes in the tree
#' @param high colour for large values
#' @param low colour for low values
#' @param mid colour for middle values
#'
#' @description Adding statistical test results onto the tree by using
#' colourful nodes and branches
#' Takes a ggtree object with test results for each node and returns
#' a ggtree graph object
#'
#' @importFrom utils modifyList
#'
#' @return an interactive ggplot graph object, plotting the hierarchical tree of clusters with
#' nodes and branches coloured by the significance testing results.
#'
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
#'                         exprs=exprs,
#'                         clusters=clusters,
#'                         samples=samples,
#'                         classes=classes,
#'                         pos_class_name=NULL,
#'                         subjects=NULL,
#'                         paired = FALSE)
#'
#' colourTree(tested_tree)
colourTree <- function(tree,
                       point_size = 1.5,
                       high = "#00c434",
                       low = "purple",
                       mid="ivory2"){

  tree$data$label <- ifelse(is.na(tree$data$label),
                            tree$data$node,
                            tree$data$label)

  tooltip <- paste("<b>Cluster</b>:", tree$data$label,
                   "\n <b>Parent</b>: stat = ", signif(tree$data$statParent,2),
                   ", p-value = ", signif(tree$data$pvalParent,2),
                   "\n <b>All</b>:    stat = ", signif(tree$data$statAll,2),
                   ", p-value = ", signif(tree$data$pvalAll,2))

  tree_df <- tree$data

  # Add colour aes mapping to the first tree layer for the branch colours
  tree$layers[[1]]$mapping <- modifyList(tree$layers[[1]]$mapping, aes(color=statParent))
  tree$layers[[2]]$mapping <- modifyList(tree$layers[[2]]$mapping, aes(color=statParent))

  tree <- tree +
    geom_point_interactive(data = tree$data,
                           aes(x,y, colour = statAll, tooltip = tooltip, data_id = label),
                           size = point_size) +
    scale_colour_gradient2(low = low,
                           mid = mid,
                           high = high,
                           midpoint = 0,
                           limits = c(min(tree$data$statParent,tree$data$statAll),
                                      max(tree$data$statParent,tree$data$statAll)),
                           name = 'Test\nStatistic')+
    guides(colour = guide_colourbar(order = 2),
           fill = guide_colourbar(order = 1)) +
    theme(legend.position = "right")

  return(tree)
}

#' Title
#'
#' @param treeDat
#'
#' @return a ggplot object, containing test statistics from testing proportions relative to parent against
#' the test statistics from testing absolute proportions.
plotSigScatter <- function(testedTree) {
  ggplot(testedTree,
         aes(x = statAll, y = statParent, shape=isTip, col=isTip,
             data_id=label,
             tooltip = scatter_tooltip))+
    geom_point_interactive() +
    # geom_point(size=1.5)+
    geom_hline(yintercept = 0, linetype="dashed")+
    geom_vline(xintercept = 0, linetype="dashed")+
    coord_equal(xlim=c(-max_val,max_val),ylim=c(-max_val,max_val))+
    labs(x = "Statistic: Relative to all",y = "Statistic: Relative to parent",
         title = "Significance between patient condition using prop. to parent vs prop. to all") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          legend.position = "top") +
    scale_color_manual(values=c("grey10", "grey50"))
}


#' Title
#' @description This function takes a hierarchical tree which has been
#' tested for proportion to all and proportion to parent cluster
#'
#' @param testedTree a ggtree object that has been run through the testTree
#' @param clust_med_df a dataframe with the cluster medians.
#' The rownumbers of the clusters median data frame should correspond to the nodes in the phylo tree.
#' The column names should also correspond to the labels you want to use
#' @param clusters a vector representing the cell type or cluster of each cell (can be character or numeric)
#' @param svg_width width of svg canvas
#' @param svg_height height of svf canvas
#' @param tr_offset distance between leaf nodes on the tree and their labels
#' @param tr_font_size font size of leaf labels
#' @param tr_point_size size of each node in the tree
#' @param tr_col_high colour used for high test statistics, coloured on the nodes and branches
#' of the tree
#' @param tr_col_low colour used for low test statistics, coloured on the nodes and branches
#' of the tree
#' @param tr_col_mid colour used for medium test statistics, coloured on the nodes and branches
#' of the tree
#' @param hm_offset distance between the tree and the heatmap
#' @param hm_tile_width width of each tile in the heatmap
#' @param hm_expand_y_lim white space below heatmap
#' @param hm_col_high colour used for large values on heatmap
#' @param hm_col_mid colour used for medium values on heatmap
#' @param hm_col_low colour used for low values on heatmap
#' @param fb_offset distance between the heatmap and frequency bars
#' @param fb_bar_length length of bar with max frequency
#' @param fb_bar_width width of each frequency bar
#' @param fb_freq_labels boolean indicated whether or not to show frequency bar labels
#'
#' @import ggiraph
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#' @importFrom dplyr select mutate
#' @importFrom patchwork plot_layout
#'
#' @return an interactive ggplot object containing the hierarchical tree of clusters
#' coloured by significance testing results, with corresponding heatmap and a scatterplot
#' comparing significance whne testing using proportions to parent vs. absolute proportions
#'
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
#'                         exprs=exprs,
#'                         clusters=clusters,
#'                         samples=samples,
#'                         classes=classes,
#'                         pos_class_name=NULL,
#'                         subjects=NULL,
#'                         paired = FALSE)
#'
#' plotInteractiveHeatmap(tested_tree,
#'                        clust_med_df = clust_tree$median_freq,
#'                        clusters=clusters)
plotInteractiveHeatmap <- function(testedTree,
                                   clust_med_df,
                                   clusters,
                                   svg_width=13,
                                   svg_height=9,
                                   tr_offset = 0.3,
                                   tr_font_size = 2,
                                   tr_point_size = 1.5,
                                   tr_col_high = "#00c434",
                                   tr_col_low = "purple",
                                   tr_col_mid="ivory2",
                                   hm_offset = 1,
                                   hm_tile_width = 1,
                                   hm_expand_y_lim = 20,
                                   hm_col_high = "#cc2010",
                                   hm_col_mid="#fff8de",
                                   hm_col_low="#66a6cc",
                                   fb_offset = 0.75,
                                   fb_bar_length = 3,
                                   fb_bar_width=0.4,
                                   fb_freq_labels = FALSE) {
  # Plot heatmap with tree results
  gTree <- testedTree %>%
    colourTree(high = tr_col_high, low = tr_col_low, mid=tr_col_mid,
               point_size = tr_point_size) %>%
    addTree(offset = tr_offset, font_size = tr_font_size) %>%
    addHeatMap(cluster_medians = clust_med_df, offset=hm_offset,
               high=hm_col_high, mid=hm_col_mid, low=hm_col_low,
               width = hm_tile_width, expand_y_lim = hm_expand_y_lim) %>%
    addFreqBars(clusters=clusters, offset = fb_offset,
                bar_length = fb_bar_length, bar_width=fb_bar_width)

  max_val <- max(abs(c(testedTree$data$statAll, testedTree$data$statParent)))
  # Make label non-null using node
  testedTreeDat <- testedTree$data %>% # Use node name where label is null
    dplyr::mutate(label = ifelse(is.na(label),node,label))
  # Insert tooltips
  scatter_tooltip <- paste0("<b>Cluster</b>: ", testedTreeDat$label,
                            "\n <b>p-value</b> (rel. to all)= ", signif(testedTreeDat$pvalAll, 2),
                            "\n <b>p-value</b> (rel. to parent)= ", signif(testedTreeDat$pvalParent, 2))
  # Plot scatterplot with parent vs. all proportions
  g1 <- plotSigScatter(testedTreeDat)

  hover_css <- "fill:cyan; stroke:darkcyan; r:4pt;"
  tooltip_css <- "border-style: solid; border-color: #c3c3c3; border-radius: 8px;
                  border-width: 0.5px; padding: 5px;box-shadow: 2px 2px 3px 0px #bbbbbb;
                  background-color: white; font: menu;"
  gf <- girafe(
    print(gTree + theme(legend.position = "top") + g1 + plot_layout(nrow = 1, widths = c(7, 4.5))),
    width_svg=svg_width,
    height_svg=svg_height)

  girafe_options(gf, opts_hover(css = hover_css), opts_tooltip(css = tooltip_css))
}
