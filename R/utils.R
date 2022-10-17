library(data.tree)
library(DiagrammeR)
library(tidyverse)

get_distance<-function(node_idx_1,node_idx_2, B){
  return(length(which(xor(B[node_idx_1,],B[node_idx_2,])==TRUE)))
}


get_root_mutation <- function(B) {
  return(which(colSums(B) == nrow(B)))
}


get_root_clone <- function(B) {
  return(which(rowSums(B) == 1))
}


add_children<-function(phylotree,node_idx,root){
  children<-get_children_idx(phylotree@B,node_idx)
  mutations<-unlist(map(children, function(x) clone_to_gene(phylotree,x)))
  node<-FindNode(root,clone_to_gene(phylotree,node_idx))
  map(mutations, function(x) node$AddChild(x))
  map(children, function(x) add_children(phylotree,x,root))
}

get_children_idx<-function(B, node_idx){
  distances<- map(1:nrow(B), function(x) get_distance(x,node_idx,B))
  candidates<-(which(distances==1))
  logic<-map_lgl(candidates, function(x) all(sum(B[x, ]==1) > sum(B[node_idx, ]==1)))
  return(unlist(candidates[logic]))
}

get_parent_idx<-function(B, node_idx){
  distances<-map(1:nrow(B), function(x) get_distance(x,node_idx,B))
  candidates<-which(distances==1)
  logic<-map_lgl(candidates, function(x) all(sum(B[x, ]==1) < sum(B[node_idx, ]==1)))
  ifelse(length(logic), candidates[logic],NULL)
}

clone_to_gene<-function(phylotree, idx){
  return(phylotree@genes[idx])
}

gene_to_clone<-function(phylotree, idx){
  return(phylotree@clones[idx])
}

get_mutation_idx <- function(B, node_idx) {
  parent <- get_parent_idx(B, node_idx)
  if (is.na(parent)) {
    return(which(rep(0, nrow(B)) != B[node_idx,]))
  } else {
    return(which(B[parent,] != B[node_idx,]))
  }
}

get_ascendants_idx <- function(B, node_idx) {
  parent <- get_parent_idx(B, node_idx)
  if (!is.na(parent)) {
    return(c(parent, get_ascendants_idx(B, parent)))
  } else {
    return()
  }
}

get_parent_idx <- function(B, node_idx) {
  indeces <- setdiff(1:nrow(B), node_idx)
  map_dbl(indeces, function(x) sum(B[node_idx, ] - B[x, ])) -> substraction_values
  possible_nodes <- indeces[which(substraction_values == 1)]
  map_lgl(possible_nodes, function(x) all(B[node_idx, ] >= B[x, ])) -> parent_logical
  parent <- possible_nodes[which(parent_logical)]
  ifelse(length(parent), parent, NA)
}

get_children_idx <- function(B, node_idx) {
  indeces <- setdiff(1:nrow(B), node_idx)
  map_dbl(indeces, function(x) sum(B[x, ] - B[node_idx, ])) -> substraction_values
  possible_nodes <- indeces[which(substraction_values == 1)]
  map_lgl(possible_nodes, function(x) all(B[x, ] >= B[node_idx, ])) -> children_logical
  children <- possible_nodes[which(children_logical)]
  if (length(children)) { # cannot use ifelse construction as it returns a value of the same length as the test
    children
  } else {
    return()
  }
}

get_descendants_idx <- function(B, node_idx) {
  children <- get_children_idx(B, node_idx)
  if (length(children)) {
    return(unlist(c(children, map(children, function(x) get_descendants_idx(B, x)))))
  } else {
    return()
  }
}

calc_F <- function(U, B, heterozygous = TRUE) {
  F_homozygous <- U %*% B
  if (heterozygous) {
    F_ <- 0.5 * F_homozygous
  } else {
    F_ <- F_homozygous
  }
  # wrap values of F between 0 and 1 as they represent clone proportions
  return(apply(F_, c(1,2), FUN = function(x) min(x,1)))
}

merge_graphs<-function(phylotree,merge,graphs,i,labels){
  if(i>length(graphs)){
    return(merge)
  }
  graph<-graphs[[i]]
  if(isTRUE(labels)){
    graph<-set_node_attrs(graph, "label", phylotree@labels)
  }
  merge<-combine_graphs(merge, graph)
  return(merge_graphs(phylotree, merge, graphs, i+1, labels))
}


