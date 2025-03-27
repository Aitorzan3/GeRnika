setClass("Node", slots="data.tree::Node")

#' S4 class to represent phylogenetic trees.
#'
#' @slot B A data.frame containing the square matrix that represents the ancestral relations among the clones of the phylogenetic tree.
#' @slot clones A vector representing the equivalence table of the clones in the phylogenetic tree.
#' @slot genes A vector representing the equivalence table of the genes in the phylogenetic tree.
#' @slot parents A vector representing the parents of the clones in the phylogenetic tree.
#' @slot tree A \code{Node} class object representing the phylogenetic tree.
#' @slot labels A vector representing the tags of the genes in the phylogenetic tree.
setClass("Phylotree", slots = list(
                                   B = "matrix", clones = "vector", 
                                   genes = "vector", parents = "vector", 
                                   tree = "Node", labels = "vector"))

get_genes<-function(B){
  genes<-c(unlist(purrr::map(1:nrow(B), function(x) get_mutation_idx(B,x))))
  return(genes)
}

get_clones<-function(genes){
  clones<-c(unlist(purrr::map(1:length(genes),  function(x) which(genes==x))))
}

#' @export
#' @title Create a \code{Phylotree} object
#' @description This is the general constructor of the \code{Phylotree} S4 class.
#'
#' @param B A square matrix that represents the phylogenetic tree.
#' @param clones A numeric vector representing the clones in the phylogenetic tree.
#' @param genes A numeric vector representing the genes in the phylogenetic tree.
#' @param parents A numeric vector representing the parents of the clones in the phylogenetic tree.
#' @param tree A \code{data.tree} object containing the tree structure of the phylogenetic tree.
#' @param labels An optional vector containing the tags of the genes in the phylogenetic tree. \code{NA} by default.
#' @return A \code{Phylotree} class object.
#'
#' @examples
#' # Create a B matrix instance
#' # composed by 10 subpopulations of
#' # clones
#' B <- create_instance(
#'        n = 10, 
#'        m = 4, 
#'        k = 1, 
#'        selection = "neutral")$B
#'
#'
#' # Create a new 'Phylotree' object
#' # on the basis of the B matrix
#' phylotree1 <- B_to_phylotree(B = B)
#'
#'
#' # Create a new 'Phylotree' object
#' # with the general constructor of
#' # the class
#' phylotree2 <- create_phylotree(
#'                 B = B, 
#'                 clones = phylotree1@clones, 
#'                 genes = phylotree1@genes, 
#'                 parents = phylotree1@parents, 
#'                 tree = phylotree1@tree)
#'
#'
#' # Generate the tags for the genes of
#' # the phyogenetic tree
#' tags <- LETTERS[1:nrow(B)]
#'
#'  
#' # Create a new 'Phylotree' object
#' # with the general constructor of
#' # the class using tags
#' phylotree_tags <- create_phylotree(
#'                     B = B, 
#'                     clones = phylotree1@clones, 
#'                     genes = phylotree1@genes, 
#'                     parents = phylotree1@parents, 
#'                     tree = phylotree1@tree, 
#'                     labels = tags)
create_phylotree <- function(B, clones, genes, parents, tree, labels = NA) {
  if (length(class(B)) == 1) {
    if (inherits(B, "data.frame")) {
      B <- as.matrix(B)
    }
  }
  if (class(B)[1] != "matrix") {
    stop("\n B must be a matrix or a data.frame class object")
  }
  if (nrow(B) != ncol(B)) {
    stop("\n B must be a square matrix")
  }
  phylotree <- methods::new("Phylotree", B=B, clones=clones, genes=genes, 
                   parents=parents, tree=tree, labels=NA)
  if (!is.na(labels[1])) {
    if (!length(labels) == length(genes)) {
      stop("\n \"labels\" and \"genes\" vectors must have the same size")
    }
    phylotree@labels<-labels[as.numeric(phylotree@tree$Get("name"))]
  } else {
    if (!is.null(colnames(B))) {
      phylotree@labels <- colnames(B)[as.numeric(phylotree@tree$Get("name"))]
    } else {
      mutation_names <- purrr::map(1:length(genes), function(x) paste0("mut", x))
      phylotree@labels <- mutation_names[as.numeric(phylotree@tree$Get("name"))]
    }
  }
  return(phylotree)
}

#' @export
#' @title Create a \code{Phylotree} object from a \code{B} matrix.
#' @description This function creates a \code{Phylotree} class object from a \code{B} matrix.
#'
#' @param B A square matrix that represents the phylogenetic tree.
#' @param labels An optional vector containing the tags of the genes in the phylogenetic tree. \code{NA} by default.
#' @return A \code{Phylotree} class object.
#'
#' @examples
#' # Create a B matrix instance
#' # composed by 10 subpopulations of
#' # clones
#' B <- create_instance(
#'        n = 10, 
#'        m = 4, 
#'        k = 1, 
#'        selection = "neutral")$B
#'
#' # Create a new 'Phylotree' object
#' # on the basis of the B matrix
#' phylotree <- B_to_phylotree(B = B)
#'
#' # Generate the tags for the genes of
#' # the phyogenetic tree
#' tags <- LETTERS[1:nrow(B)]
#'
#' # Create a new 'Phylotree' object
#' # on the basis of the B matrix and
#' # the list of tags
#' phylotree_tags <- B_to_phylotree(
#'                     B = B, 
#'                     labels = tags)
B_to_phylotree <- function(B, labels = NA) {
  if (length(class(B)) == 1) {
    if (inherits(B, "data.frame")) {
      B <- as.matrix(B)
    }
  }
  if (class(B)[1] != "matrix") {
    stop("\n B must be a matrix or a data.frame class object")
  }
  if (nrow(B) != ncol(B)) {
    stop("\n B must be a square matrix")
  }
  genes <- get_genes(B)
  clones <- get_clones(genes)
  new <- data.tree::Node$new(get_root_mutation(B))
  phylotree <- methods::new("Phylotree", B = B, clones = clones, genes = genes, 
                   parents = integer(nrow(B)), tree = new, labels = NA)
  phylotree@parents <- get_parents(phylotree)
  create_tree(phylotree)
  if (!is.na(labels[1])) {
    if (!length(labels) == length(genes)) {
      stop("\n \"labels\" and \"genes\" vectors must have the same size")
    }
    phylotree@labels <- labels[as.numeric(phylotree@tree$Get("name"))]
  } else {
    if (!is.null(colnames(B))) {
      phylotree@labels <- colnames(B)[as.numeric(phylotree@tree$Get("name"))]
    } else {
      mutation_names<-purrr::map(1:length(genes), function(x) paste0("mut", x))
      phylotree@labels<-unlist(mutation_names[as.numeric(phylotree@tree$Get("name"))])
    }
  }
  return(phylotree)
}


phylotree_to_B <- function(phylotree) {
  return(phylotree@B)
}


get_parents<-function(phylotree){
  parents<-c(unlist(purrr::map(1:nrow(phylotree@B), function(x) get_parent(phylotree,x))))
  parents[get_root_mutation(phylotree@B)]<--1
  return(parents)
}


phylotree_to_tree<-function(phylotree){
  return(phylotree@tree)
}


plot_phylotree<-function(phylotree, labels=FALSE){
  if(!is.logical(labels)){
    stop("\n labels must be logical")
  }
  tree<-data.tree::Clone(phylotree@tree)
  if(isTRUE(labels)){
    tree$Set(name=phylotree@labels)
  }
  DiagrammeR::render_graph(data.tree::ToDiagrammeRGraph(tree))
}


plot_p<-function(phylotree, proportions){
  if(length(phylotree@clones) != length(proportions)){
    stop("\n the proportion vectors length must be equal to the number of clones in the tree")
  }
  graph<-data.tree::ToDiagrammeRGraph(data.tree::Clone(phylotree@tree))
  order<-as.numeric(DiagrammeR::get_node_attrs_ws(DiagrammeR::select_nodes(graph), node_attr='label'))
  proportions_genes<-unlist(purrr::map(1:length(proportions), function(x) proportions[gene_to_clone(phylotree,x)]))
  ordered<-unlist(purrr::map(1:length(order), function(x) proportions_genes[order[x]]))
  circles<-unlist(purrr::map(ordered, function(x) return(x*3+0.5)))
  sizes<-unlist(purrr::map(ordered, function(x) return(x*140+10)))
  colors<-unlist(purrr::map(ordered, function(x) return(colorspace::adjust_transparency(GeRnika::palettes$Simpsons[3], alpha = x*0.95+0.05))))
  graph<-DiagrammeR::set_node_attrs(DiagrammeR::set_node_attrs(DiagrammeR::set_node_attrs(graph, node_attr = 'fontsize', values = unlist(sizes)), node_attr = 'width', values = unlist(circles)), node_attr = 'height', values = unlist(circles))
  graph<-DiagrammeR::set_node_attrs(DiagrammeR::set_node_attrs(graph, node_attr='fontcolor', values = unlist(colors)), node_attr='color', values = unlist(colors))
  return(graph)
}

#' @export
#' @title Plot a phylogenetic tree with proportional node sizes and colors
#' @description This function plots a phylogenetic tree with nodes sized and colored according to the proportions of each clone. If a matrix of proportions is provided, multiple phylogenetic trees will be plotted, each corresponding to a row of proportions.
#'
#' @param phylotree A \code{Phylotree} class object representing the phylogenetic tree to be plotted.
#' @param proportions A numeric vector or matrix representing the proportions of each clone in the phylogenetic tree. If a matrix is provided, each row should represent the proportions for a separate tree.
#' @param labels A logical value indicating whether to label the nodes with gene tags (if \code{TRUE}) or gene indices (if \code{FALSE}). Default is \code{FALSE}.
#'
#' @return A graph representing the phylogenetic tree, with node sizes and colors reflecting clone proportions.
#' @examples
#' # Create an instance
#' # composed by 5 subpopulations of clones
#' # and 4 samples
#' instance <- create_instance(
#'        n = 5, 
#'        m = 4, 
#'        k = 1, 
#'        selection = "neutral")
#'        
#' # Extract its associated B matrix
#' B <- instance$B
#'
#' # Create a new 'Phylotree' object
#' # on the basis of the B matrix
#' phylotree <- B_to_phylotree(B = B)
#'
#' # Generate the tags for the genes of
#' # the phyogenetic tree
#' tags <- LETTERS[1:nrow(B)]
#'
#' # Plot the phylogenetic tree taking
#' # into account the proportions of the
#' # previously generated instance
#' plot_proportions(phylotree, instance$U, labels=TRUE)
plot_proportions <- function(phylotree, proportions, labels = FALSE) {
  if (is.matrix(proportions)) {
    graphs <- purrr::map(1:nrow(proportions), function(x) plot_p(phylotree, proportions[x,]))
    graph <- graphs[[1]]
    if (isTRUE(labels)) {
      graph <- DiagrammeR::set_node_attrs(graph, "label", phylotree@labels)
    }
    merged <- graph
    merged <- merge_graphs(phylotree, merged, graphs, 2, labels)
  } else if (is.vector(proportions)) {
    merged <- plot_p(phylotree, proportions)
    if (isTRUE(labels)) {
      merged <- DiagrammeR::set_node_attrs(merged, "label", phylotree@labels)
    }
  } else {
    stop("\n the proportions parameter must be a matrix or a vector")
  }

  merged <- DiagrammeR::set_edge_attrs(merged, edge_attr = 'color', values = 'black')
  DiagrammeR::render_graph(merged)
}

#' Plot a Phylotree object.
#' @param object A \code{Phylotree} object.
#' @param labels A label vector.
#' @rdname plot
#' @docType methods
#' @export 
#' 
setGeneric("plot", function(object, labels = FALSE) standardGeneric("plot"))
#' @rdname plot
setMethod("plot", signature(object = "Phylotree", labels = "ANY"), 
          function(object, labels) plot_phylotree(object, labels))
