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
  mutations<-unlist(purrr::map(children, function(x) clone_to_gene(phylotree,x)))
  node<-data.tree::FindNode(root,clone_to_gene(phylotree,node_idx))
  purrr::map(mutations, function(x) node$AddChild(x))
  purrr::map(children, function(x) add_children(phylotree,x,root))
}

get_children_idx<-function(B, node_idx){
  distances<- purrr::map(1:nrow(B), function(x) get_distance(x,node_idx,B))
  candidates<-(which(distances==1))
  logic<-purrr::map_lgl(candidates, function(x) all(sum(B[x, ]==1) > sum(B[node_idx, ]==1)))
  return(unlist(candidates[logic]))
}

get_parent_idx<-function(B, node_idx){
  distances<-purrr::map(1:nrow(B), function(x) get_distance(x,node_idx,B))
  candidates<-which(distances==1)
  logic<-purrr::map_lgl(candidates, function(x) all(sum(B[x, ]==1) < sum(B[node_idx, ]==1)))
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
  purrr::map_dbl(indeces, function(x) sum(B[node_idx, ] - B[x, ])) -> substraction_values
  possible_nodes <- indeces[which(substraction_values == 1)]
  purrr::map_lgl(possible_nodes, function(x) all(B[node_idx, ] >= B[x, ])) -> parent_logical
  parent <- possible_nodes[which(parent_logical)]
  ifelse(length(parent), parent, NA)
}

get_children_idx <- function(B, node_idx) {
  indeces <- setdiff(1:nrow(B), node_idx)
  purrr::map_dbl(indeces, function(x) sum(B[x, ] - B[node_idx, ])) -> substraction_values
  possible_nodes <- indeces[which(substraction_values == 1)]
  purrr::map_lgl(possible_nodes, function(x) all(B[x, ] >= B[node_idx, ])) -> children_logical
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
    return(unlist(c(children, purrr::map(children, function(x) get_descendants_idx(B, x)))))
  } else {
    return()
  }
}

merge_graphs<-function(phylotree,merge,graphs,i,labels){
  if(i>length(graphs)){
    return(merge)
  }
  graph<-graphs[[i]]
  if(isTRUE(labels)){
    graph<-DiagrammeR::set_node_attrs(graph, "label", phylotree@labels)
  }
  merge<-DiagrammeR::combine_graphs(merge, graph)
  return(merge_graphs(phylotree, merge, graphs, i+1, labels))
}

rdir <-  function(n, alpha) {
    l <- length(alpha)
    x <- matrix(stats::rgamma(l*n,alpha),ncol=l,byrow=TRUE)
    sm <- x%*%rep(1,l)
    return(x/as.vector(sm))
}

sample.vec <- function(x, ...) x[sample(length(x), ...)]

.distribute_freqs <- function(B, clone_idx, clone_proportions, selection) {
  dirich_params <- dplyr::tibble(positive = rep(0.3, 2),
                                 neutral = c(5, 10))
  children_clone_idx <- get_children_idx(B = B, node_idx = clone_idx)
  if (!is.null(children_clone_idx)) {
    params <- dplyr::pull(dirich_params, eval(selection))
    dirich_values <- rdir(1, c(params[1],
                               rep(params[2], length(children_clone_idx))))
    # Normalize the values to the parent clone's proportion
    norm_dirich_values <- dirich_values*clone_proportions[clone_idx]
    clone_proportions[c(clone_idx, children_clone_idx)] <- norm_dirich_values
  }
  return(clone_proportions)
}

calc_clone_proportions <- function(B, selection) {
  n <- nrow(B)
  clone_proportions <- rep(1, n)
  root_clone_idx <- get_root_clone(B)
  for (clone_idx in c(root_clone_idx, get_descendants_idx(B, root_clone_idx))) {
    clone_proportions <- .distribute_freqs(B, clone_idx, clone_proportions, selection)
  }
  # sum(clone_proportions)
  clone_proportion_df <- dplyr::tibble(clone_idx = paste0("clon", 1:n), proportion = clone_proportions)
  # clon_i corresponds to clone containing mutation i for the first time and not necessarility the clone in row i, but in this case it is also the clone in row i because of the way we construct the B matrix
  return(clone_proportion_df)
}

place_clones_space <- function(B) {
  n <- nrow(B)
  . <- NULL
  max_sep <- 4
  mean_diffs <- seq(from = 0, to = max_sep, by = 0.1)
  clone_order <- sample(1:n, n)
  clone_mean_diff <- sample(x = mean_diffs,
                            size = n-1,
                            prob = stats::dbeta(x = seq(0, 1, length.out = length(mean_diffs)), shape1 = 1, shape2 = 5), replace = TRUE) # change beta distribution parameters for tumor density
  clone_means <- c(0, cumsum(clone_mean_diff))
  clone_sd <- 1
  x <- seq(min(clone_means - 3*clone_sd), max(clone_means + 3*clone_sd), .01)
  y <- purrr::map(1:n, function(z) stats::dnorm(x, mean = clone_means[z], sd = clone_sd)) %>%
    do.call(cbind, .) %>%
    dplyr::as_tibble(.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)) %>%
    magrittr::set_colnames(paste0("clon", clone_order))
  spatial_coords <- y %>%
    dplyr::mutate(idx = x)
  return(list(spatial_coords = spatial_coords, x = x))
}

