library(data.tree)
library(DiagrammeR)
library(tidyverse)

source("R/tree.R")
source("R/clone.R")
source("R/utils.R")

setClass("Node", slots="Node")
#' @exportClass
#' @name Phylotree_class
#' @title Phylotree_class
#' S4 class to represent phylogenetic trees.
#'
#' @slot B the data.frame containing the square matrix that represents the clones of the phylogenetic tree.
#' @slot clones a vector representing the equivalence table of the clones in the phylogenetic tree.
#' @slot genes a vector representing the equivalence table of the genes in the phylogenetic tree.
#' @slot parents a vector representing the parents of the clones in the phylogenetic tree.
#' @slot tree a \code{Node} class object representing the phylogenetic tree.
#' @slot labels a vector representing the tags of the genes in the phylogenetic tree.
setClass("Phylotree", slots=list(B="matrix", clones="vector", genes="vector", parents="vector", tree="Node", labels="vector"))


get_genes<-function(B){
  genes<-c(unlist(map(1:nrow(B), function(x) get_mutation_idx(B,x))))
  return(genes)
}

get_clones<-function(genes){
  clones<-c(unlist(map(1:length(genes),  function(x) which(genes==x))))
}

#' @export
#' @title Create a \code{Phylotree} object.
#' @description The general constructor of the \code{Phylotree} S4 class.
#'
#' @param B The Matrix that represents the clones in the phylogenetic tree.
#' @param clones numeric vector representing the clones in the phylogenetic tree.
#' @param genes numeric vector representing the genes in the phylogenetic tree.
#' @param parents numeric vector representing the parents the clones in the phylogenetic tree.
#' @param tree \code{data.tree} object containing the tree structure of the phylogenetic tree.
#' @param labels Optional argument that refers to the list containing the tags of the genes of the phylogenetic tree. \code{NA} by default.
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
#'                 clones = tree@clones, 
#'                 genes = tree@genes, 
#'                 parents = tree@parents, 
#'                 tree = tree@tree)
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
#'                     clones = tree@clones, 
#'                     genes = tree@genes, 
#'                     parents = tree@parents, 
#'                     tree = tree@tree, 
#'                     labels = tags)
create_phylotree<-function(B, clones, genes, parents, tree, labels=NA){
  if(length(class(B))==1){
    if(class(B)=="data.frame"){
      B<-as.matrix(B)
    }
  }
  if(class(B)[1]!="matrix"){
    stop("\n B must be a matrix or a data.frame class object")
  }
  if(nrow(B) != ncol(B)){
    stop("\n B must be a square matrix")
  }
  phylotree<-new("Phylotree", B=B, clones=clones, genes=genes, parents=parents, tree=tree, labels=NA)
  if(!is.na(labels[1])){
    if(!length(labels)==length(genes)){
      stop("\n \"labels\" and \"genes\" vectors must have the same size")
    }
    phylotree@labels<-labels[as.numeric(phylotree@tree$Get("name"))]
  }
  else{
    if(!is.null(colnames(B))){
      phylotree@labels<-colnames(B)[as.numeric(phylotree@tree$Get("name"))]
    }
    else{
      mutation_names<-map(1:length(genes), function(x) paste0("mut", x))
      phylotree@labels<-mutation_names[as.numeric(phylotree@tree$Get("name"))]
    }
  }
  return(phylotree)
}

#' @export
#' @title Create a \code{Phylotree} object from a \code{B} matrix.
#' @description Creates a \code{Phylotree} class object from a \code{B} matrix.
#'
#' @param B The Matrix that represents the clones in the phylogenetic tree.
#' @param labels Optional argument that refers to the vector containing the tags of the genes of the phylogenetic tree. \code{NA} by default.
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
B_to_phylotree<-function(B, labels=NA){
  if(length(class(B))==1){
    if(class(B)=="data.frame"){
      B<-as.matrix(B)
    }
  }
  if(class(B)[1]!="matrix"){
    stop("\n B must be a matrix or a data.frame class object")
  }
  if(nrow(B) != ncol(B)){
    stop("\n B must be a square matrix")
  }
  genes<-get_genes(B)
  clones<-get_clones(genes)
  new<-Node$new(get_root_mutation(B))
  phylotree<-new("Phylotree", B=B, clones=clones, genes=genes, parents=integer(nrow(B)), tree=new, labels=NA)
  phylotree@parents<-get_parents(phylotree)
  create_tree(phylotree)
  if(!is.na(labels[1])){
    if(!length(labels)==length(genes)){
      stop("\n \"labels\" and \"genes\" vectors must have the same size")
    }
    phylotree@labels<-labels[as.numeric(phylotree@tree$Get("name"))]
  }
  else{
    if(!is.null(colnames(B))){
      phylotree@labels<-colnames(B)[as.numeric(phylotree@tree$Get("name"))]
    }
    else{
      mutation_names<-map(1:length(genes), function(x) paste0("mut", x))
      phylotree@labels<-unlist(mutation_names[as.numeric(phylotree@tree$Get("name"))])
    }
  }
  return(phylotree)
}

#' @export
#' @title Get B from \code{Phylotree}
#' @description Returns the B matrix of a \code{Phylotree} object.
#'
#' @param phylotree a \code{phylotree} class object.
#' @return A \code{data.frame} representing the B matrix of the phylogenetic tree.
#' @examples
#' # Get the B matrix of a tumor instance
#' # composed by 10 subpopulations of
#' # clones
#' B <- create_instance(10, 4, 1, 1)$B
#'
#' # Create a new 'Phylotree' object
#' # on the basis of the B matrix
#' phylotree <- B_to_phylotree(B)
#'
#' # Get the B matrix of the phyotree
#' b1 <- phylotree_to_B(phylotree)
phylotree_to_B<-function(phylotree){
  return(phylotree@B)
}


get_parents<-function(phylotree){
  parents<-c(unlist(map(1:nrow(phylotree@B), function(x) get_parent(phylotree,x))))
  parents[get_root_mutation(phylotree@B)]<--1
  return(parents)
}

phylotree_to_tree<-function(phylotree){
  return(phylotree@tree)
}

plot_phylotree<-function(phylotree, labels=FALSE){
  if(!class(labels)=='logical'){
    stop("\n labels must be logical")
  }
  tree<-Clone(phylotree@tree)
  if(isTRUE(labels)){
    tree$Set(name=phylotree@labels)
  }
  render_graph(ToDiagrammeRGraph(tree))
}

plot_proportions<-function(phylotree, proportions){
  if(length(phylotree@clones) != length(proportions)){
    stop("\n the proportion vectors length must be equal to the number of clones in the tree")
  }
  graph<-ToDiagrammeRGraph(Clone(phylotree@tree))
  order<-as.numeric(get_node_attrs_ws(select_nodes(graph), node_attr=label))
  proportions_genes<-unlist(map(1:length(proportions), function(x) proportions[gene_to_clone(phylotree,x)]))
  ordered<-unlist(map(1:length(order), function(x) proportions_genes[order[x]]))
  circles<-unlist(map(ordered, function(x) return(x*3+0.5)))
  sizes<-unlist(map(ordered, function(x) return(x*140+10)))
  colors<-unlist(map(ordered, function(x) return(adjust_transparency(GeRnika::palettes$Simpsons[3], alpha = x*0.95+0.05)))
  graph<-set_node_attrs(set_node_attrs(set_node_attrs(graph, node_attr = fontsize, values = unlist(sizes)), node_attr = width, values = unlist(circles)), node_attr = height, values = unlist(circles))
  graph<-set_node_attrs(set_node_attrs(graph, node_attr=fontcolor, values = unlist(colors)), node_attr=color, values = unlist(colors))
  return(graph)
}

#' @exportMethod
setGeneric("plot", function(object, labels=FALSE)standardGeneric("plot"))
setMethod("plot", signature(object="Phylotree", labels="ANY"), function(object, labels) plot_phylotree(object, labels))
