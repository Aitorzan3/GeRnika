get_children<-function(phylotree, node_idx){
  children_clones<-get_children_idx(phylotree@B, gene_to_clone(phylotree, node_idx))
  return(unlist(purrr::map(children_clones, function(x) clone_to_gene(phylotree,x))))
}


get_descendants <- function(phylotree, node_idx) {
  children<-get_children_idx(phylotree@B, gene_to_clone(phylotree,node_idx))
  children_mutations<-c(purrr::map(children, function(x) clone_to_gene(phylotree,x)))
  if (length(children)) {
    return(unlist(c(children_mutations, purrr::map(children_mutations, function(x) get_descendants(phylotree, x)))))
  } else {
    return()
  }
}


get_ascendants <- function(phylotree, node_idx) {
  parent <- get_parent(phylotree, node_idx)
  if (is.na(parent)) {
  } else {
    return(c(parent,get_ascendants(phylotree,parent)))
  }
}


get_parent<-function(phylotree, node_idx){
  parent_idx<-get_parent_idx(phylotree@B,gene_to_clone(phylotree,node_idx))
  ifelse(length(parent_idx), clone_to_gene(phylotree,parent_idx),NULL)
}


is_leave<-function(phylotree, node_idx){
  return(is.null(get_children(phylotree,node_idx)))
}
