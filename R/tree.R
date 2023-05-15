
source("R/clone.R")
source("R/utils.R")

get_depth<-function(phylotree){
  return(max(colSums(phylotree@B)))
}

create_tree<-function(phylotree){
  root_idx<-get_root_clone(phylotree@B)
  add_children(phylotree,root_idx,phylotree@tree)
}

#' @export
#' @title Find the set of common subtrees between two phylogenetic trees
#' @description Plots the common subtrees between two phylogenetic trees and prints the information about their similarities and their differences.
#'
#' @param phylotree_1 A \code{Phylotree} class object.
#' @param phylotree_2 A \code{Phylotree} class object.
#' @param labels A boolean, if \code{TRUE} the rendered graph will be plotted with the tags of the genes in the phylogenetic trees instead of their gene index. \code{FALSE} by default.
#' @examples
#'
#' # Load the predefined B matrices of the package
#' B_mats <- GeRnika::B_mats
#'
#'
#  # Get the B matrixes from B_mats for comparing them
#' B_real <- B_mats[[2]]$B_real
#' B_opt <- B_mats[[2]]$B_opt
#' 
#' 
#' # Generate the tags for the genes of
#' # the phyogenetic tree
#' tags <- LETTERS[1:nrow(B)]
#' 
#' 
#' # Instantiate two Phylotree class objects on 
#' # the basis of the B matrices using tags
#' phylotree_real <- B_to_phylotree(
#'                     B = B_real, 
#'                     labels = tags)
#'                     
#' phylotree_opt <- B_to_phylotree(
#'                     B = B_opt, 
#'                     labels = tags)
#' 
#' 
#' # find the set of common subtrees between both 
#' # phylogenetic trees
#' find_common_subtrees(
#'   phylotree_1 = phylotree_real, 
#'   phylotree_2 = phylotree_opt)
#' 
#' 
#' # find the set of common subtrees between both
#' # phylogenetic trees using tags
#' find_common_subtrees(
#'   phylotree_1 = phylotree_real, 
#'   phylotree_2 = phylotree_opt, 
#'   labels = TRUE)
find_common_subtrees <- function(phylotree_1, phylotree_2, labels = FALSE) {
  tree1 <- (ToDiagrammeRGraph(phylotree_1@tree))
  tree2 <- (ToDiagrammeRGraph(phylotree_2@tree))
  intersection <- which(phylotree_1@parents == phylotree_2@parents)
  e1 <- length(get_edges(tree1))
  e2 <- length(get_edges(tree2))
  et <- e1+e2
  if (length(intersection) == 0) {
    cat("tree1 and tree2 have no common subtrees \n")
    cat("Independent edges of tree1: ", e1, "\n")
    cat("Independent edges of tree2: ", e2, "\n")
    cat("Common edges: 0 \n")
    cat("Distance: ", et, "\n")
  } else {
    intersection <- unique(c(intersection, 
                             phylotree_1@parents[intersection][phylotree_1@parents[intersection] > 0]))
    inter1 <- transform_to_subgraph_ws(select_nodes(tree1, conditions = (label %in% intersection)))
    inter2 <- transform_to_subgraph_ws(select_nodes(tree2, conditions = (label %in% intersection)))
    inter_edges <- get_edge_ids(tree1)[which(get_edges(tree1, return_values = "label") 
                                             %in% get_edges(tree2, return_values = "label"))]
    inter <- (transform_to_subgraph_ws(select_edges_by_edge_id(inter1, inter_edges)))
    ec <- length(inter_edges)
    e1 <- e1 - ec
    e2 <- e2 - ec
    cat("Independent edges of tree1: ", e1, "\n")
    cat("Independent edges of tree2: ", e2, "\n")
    cat("Common edges: ", ec, "\n")
    d <- e1 + e2
    cat("Distance: ", d, "\n")
    if (isTRUE(labels)) {
      inter <- set_node_attrs(inter, "label", phylotree_1@labels[as.numeric(get_node_ids(inter))])
    }
    render_graph(inter)
  }
}

#' @export
#' @title Check if two phylogenetic trees are equal
#' @description Checks wether two phylogenetic trees are equivalent or not.
#'
#' @param phylotree_1 A \code{Phylotree} class object.
#' @param phylotree_2 A \code{Phylotree} class object.
#' @return A boolean, \code{TRUE} if they are equal and \code{FALSE} if not.
#' @examples 
#' 
#' # Load the predefined B matrices of the package
#' B_mats <- GeRnika::B_mats
#'
#'
#  # Get the B matrixes from B_mats for comparing them
#' B_real <- B_mats[[2]]$B_real
#' B_opt <- B_mats[[2]]$B_opt
#' 
#' 
#' # Instantiate two \code{Phylotree} class objects on 
#' # the basis of the B matrices
#' phylotree_real <- B_to_phylotree(
#'                     B = B_real)
#'                     
#' phylotree_opt <- B_to_phylotree(
#'                     B = B_opt)
#' 
#' 
#' equals(phylotree_real, phylotree_opt)
equals <- function(phylotree_1, phylotree_2) {
  return(all(phylotree_1@parents == phylotree_2@parents))
}

#' @export
#' @title Get consensus tree between two phylogenetic trees
#' @description Returns a graph representing the consensus tree between two phylogenetic trees.
#'
#' @param phylotree_1 A \code{Phylotree} class object.
#' @param phylotree_2 A \code{Phylotree} class object.
#' @param palette A vector composed by the hexadecimal code of three colors. "The Simpsons" palette used as default.
#' @param labels A boolean, if \code{TRUE} the resulting graph will be plotted with the tags of the genes in the phylogenetic trees instead of their mutation index. \code{FALSE} by default.
#' @return a \code{dgr_graph} object representing the consensus graph between \code{phylotree_1} \code{phylotree_2}.
#' @examples 
#' 
#' # Load the predefined B matrices of the package
#' B_mats <- GeRnika::B_mats
#'
#'
#  # Get the B matrixes from B_mats for comparing them
#' B_real <- B_mats[[2]]$B_real
#' B_opt <- B_mats[[2]]$B_opt
#' 
#' 
#' # Generate the tags for the genes of
#' # the phyogenetic tree
#' tags <- LETTERS[1:nrow(B)]
#' 
#' 
#' # Instantiate two \code{Phylotree} class objects on 
#' # the basis of the B matrices
#' phylotree_real <- B_to_phylotree(
#'                     B = B_real, 
#'                     labels = tags)
#'                     
#' phylotree_opt <- B_to_phylotree(
#'                     B = B_opt, 
#'                     labels = tags)
#' 
#' 
#' # Create the consensus tree between phylotree_real
#' # and phylotree_opt
#' consensus <- combine_trees(
#'                phylotree_1 = phylotree_real,
#'                phylotree_2 = phylotree_opt)
#'                
#'                
#' # Render the consensus tree
#' render_graph(consensus)
#' 
#' 
#' # Load another palette
#' palette_1 <- GeRnika::palette["Lancet"]
#' 
#' 
#' # Create the consensus tree between phylotree_real
#' # and phylotree_opt using tags and another palette
#' consensus_tag <- combine_trees(
#'                    phylotree_1 = phylotree_real, 
#'                    phylotree_2 = phylotree_opt
#'                    palette = palette_1
#'                    labels = TRUE)
#' 
#' 
#' # Render the consensus tree using tags and the
#' # selected palette
#' render(consensus_tag)
combine_trees <- function(phylotree_1, phylotree_2, palette = GeRnika::palettes$Simpsons, labels = FALSE) {
  if (length(palette) < 3) {
    stop("\n palette argument has length<2, more than 3 elements are needed")
  }
  if (length(palette) > 3) {
    warning("\n palette argument has length>3, only the first 3 elements will be used")
  }
  color_tree1 <- adjust_transparency(palette[1], alpha = 0.25)
  color_tree2 <- adjust_transparency(palette[2], alpha = 0.25)
  color_combined <- palette[3]
  not_color <- adjust_transparency("black", alpha = 0.2)
  tree1 <- ToDiagrammeRGraph(Clone(phylotree_1@tree))
  tree2 <- ToDiagrammeRGraph(Clone(phylotree_2@tree))
  if (!(length(phylotree_1@clones)==length(phylotree_2@clones))) {
    stop("\n Phylotrees must have the same size")
  }
  if (equals(phylotree_1,phylotree_2)) {
    tree1 <- set_node_attrs(set_edge_attrs(tree1, edge_attr = color, 
                                         values = color_combined), 
                          node_attr = color, values = color_combined)
    return(tree1)
  } else {
    tree1 <- set_edge_attrs(tree1, edge_attr = color, values = color_tree1)
    intersection <- which(phylotree_1@parents == phylotree_2@parents)
    if (!length(intersection) == 0) {
      intersection <- unique(c(intersection, 
                               phylotree_1@parents[intersection][phylotree_1@parents[intersection] > 0]))
      inters1 <- select_nodes(tree1, conditions = (label %in% intersection))
      tree1 <- set_node_attrs_ws(inters1, node_attr = color, value = color_combined)
      if(length(intersection) != length(phylotree_1@genes)) {
        tree1 <- set_node_attrs_ws(invert_selection(select_nodes(tree1, conditions=(label %in% intersection))), 
                                    node_attr = color, value = not_color)
        tree1 <- set_node_attrs_ws(tree1, node_attr = fontcolor, value = "grey")
      }
      inter_edges <- get_edge_ids(tree1)[which(get_edges(tree1, return_values = "label") 
                                               %in% get_edges(tree2, return_values = "label"))]
      inters <- select_edges_by_edge_id(tree1, inter_edges)
      tree1 <- set_edge_attrs_ws(inters, edge_attr = color, value = color_combined)
    }
    else {
      tree1 <- set_node_attrs_ws(select_nodes(tree1), node_attr = color, value = not_color)
      tree1 <- set_node_attrs_ws(tree1, node_attr = fontcolor, value = "grey")
    }
    edges2 <- get_edges(tree2, return_values = "label")[which(!(get_edges(tree2, return_values = "label") 
                                                                %in% get_edges(tree1, return_values ="label")))]
    tree1 <- add_edges_w_string(tree1, edges2, use_labels = TRUE)
    tree1 <- set_edge_attrs_ws(invert_selection(select_edges(tree1, conditions = (color == color_combined | color ==  color_tree1))), 
                               edge_attr = color, value = color_tree2)
    if (isTRUE(labels)) {
      tree1 <- set_node_attrs(tree1, "label", phylotree_1@labels)
    }
    return(tree1)
  }
}

get_leaves<-function(phylotree){
  leaves<-map(1:get_depth(phylotree_1), function(x) is_leave(phylotree_1,x))
  return(which(leaves==TRUE))
}

