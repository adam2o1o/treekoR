#' Title
#' @description a function to create a skeleton tree diagram to display significance
#' testing results on each node
#'
#' @param offset
#' @param width
#' @param color
#' @param cluster_angle
#' @param colnames
#' @param colnames_position
#' @param colnames_angle
#' @param colnames_level
#' @param colnames_offset_x
#' @param colnames_offset_y
#' @param font_size
#' @param hjust
#' @param legend_title
#' @param p a phylogenetic tree plot created from the ggtree() function
#'
#' @export
#' @importFrom ggtree geom_tiplab

addTree <- function(p,
                    offset = 0.3,
                    width = 0,
                    color = NULL,
                    cluster_angle = 0,
                    colnames = TRUE,
                    colnames_position = "bottom",
                    colnames_angle = NULL,
                    colnames_level = NULL,
                    colnames_offset_x = 0,
                    colnames_offset_y = -0.55,
                    font_size = 2.5,
                    hjust = 0,
                    legend_title = "Clusters") {

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
#' @param width the width of each tile
#'
#' @export

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

  # Create Heatmap Matrix using y values from tree and getting the x values
  heatmap_df <- cluster_medians %>%
    dplyr::mutate(y = tree_df$y[match(rownames(cluster_medians), tree_df$label)],
                  label = tree_df$label[match(rownames(cluster_medians), tree_df$label)]) %>%
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
#' @param p a phylogenetic tree plot created from the ggtree() function
#'
#' @export

addFreqBars <- function(p,
                        clusters,
                        offset = 0.75,
                        bar_length = 3,
                        bar_width=0.4,
                        freq_labels = FALSE) {


  tree_df <- p$data

  # Get starting x value by getting all x values plotted in ggplot
  start_x <- max(sapply(ggplot_build(p)$data, function(dat) if(!is.null(dat$x)) max(dat$x) else 0)) + offset

  # Calculate frequency of each cluster
  cluster_freq_df <- data.frame(
    freq = as.numeric(table(clusters)),
    cluster = names(table(clusters))
  )

  # Add dimensions for bars
  cluster_freq_df <- cluster_freq_df %>%
    dplyr::mutate(y = tree_df$y[match(cluster_freq_df$cluster, tree_df$label)],
                  ymin = tree_df$y[match(cluster_freq_df$cluster, tree_df$label)] - bar_width,
                  ymax = tree_df$y[match(cluster_freq_df$cluster, tree_df$label)] + bar_width,
                  xmin = start_x,
                  xmax = start_x + bar_length * cluster_freq_df$freq / max(cluster_freq_df$freq)) %>%
    dplyr::mutate(x = max(xmax) + 0.5) %>%
    dplyr::mutate(freq_label = paste0(round(100*freq/sum(freq), 2), "% (", freq, ")"))

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
                aes(x=x, y=y, label= freq_label),
                size=1.75,
                hjust=0)
  }

  return(p2)
}


#' Title
#' @description Adding statistical test results onto the tree by using
#' colourful nodes and branches
#' Takes a ggtree object with test results for each node and returns
#' a ggtree graph object
#'
#' @param p a phylogenetic tree plot created from the ggtree() function
#' with p$data containing test statisic and p-
#'
#' @export

colourTree <- function(tree,
                       point_size = 1.5,
                       high = "red",
                       mid = "lightgrey",
                       low = "blue"){


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
#' @description This function takes a hierarchical tree which has been
#' tested for proportion to all and proportion to parent cluster
#'
#' @param testedTree a ggtree object that has been run through the testTree
#' function
#' @import ggiraph
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#' @importFrom dplyr select mutate
#' @importFrom patchwork plot_layout
#'
#' @export

plotInteractiveHeatmap <- function(testedTree,
                                   clust_med_df,
                                   clusters,
                                   svg_width=13,
                                   svg_height=9) {
  # Plot heatmap with tree results
  gTree <- testedTree %>%
    colourTree(high = "#00c434", low = "purple", mid="ivory2", point_size = 1.5) %>%
    addTree(font_size = 2) %>%
    addHeatMap(cluster_medians = clust_med_df, offset=1,
               high="#cc2010", mid="#fff8de", low="#66a6cc") %>%
    addFreqBars(clusters=clusters, offset = 0.75, bar_length = 3, bar_width=0.4)

  max_val <- max(abs(c(testedTree$data$statAll, testedTree$data$statParent)))

  # Make label non-null using node
  testedTreeDat <- testedTree$data %>% # Use node name where label is null
    mutate(label = ifelse(is.na(label),node,label))

  # Insert tooltips
  scatter_tooltip <- paste0("<b>Cluster</b>: ", testedTreeDat$label,
                            "\n <b>p-value</b> (rel. to all)= ", signif(testedTreeDat$pvalAll, 2),
                            "\n <b>p-value</b> (rel. to parent)= ", signif(testedTreeDat$pvalParent, 2)
  )

  # Plot scatterplot with parent vs. all proportions
  g1 <- ggplot(testedTree$data %>% # Use node name where label is null
                 mutate(label = ifelse(is.na(label),node,label)),
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
  # + xlim(-max(abs(c(gTree$data$statAll, gTree$data$statParent))), max(abs(c(gTree$data$statAll, gTree$data$statParent))))

  hover_css <- "fill:cyan;
              stroke:darkcyan;
              r:4pt;"
  tooltip_css <- "border-style: solid;
                            border-color: #c3c3c3;
                            border-radius: 8px;
                            border-width: 0.5px;
                            padding: 5px;
                            box-shadow: 2px 2px 3px 0px #bbbbbb;
                            background-color: white;
                            font: menu;"

  gf <-girafe(
    print(gTree +
            theme(legend.position = "top")+
            g1 +
            plot_layout(nrow = 1, widths = c(7, 4.5))),
    width_svg=svg_width,
    height_svg=svg_height)

  girafe_options(gf,
                 opts_hover(css = hover_css),
                 opts_tooltip(css = tooltip_css))
}
