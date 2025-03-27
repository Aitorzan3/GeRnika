

#' Create tumor phylogenetic tree topology
#'
#' This function generates a mutation matrix (B matrix) for a tumor phylogenetic tree with a given number of nodes. This matrix represents the topology and it is created randomly, with the probability of a node to be chosen as a parent of a new node being proportional to the number of its ascendants raised to the power of a constant `k`.
#'
#' @param n An integer representing the number of nodes in the phylogenetic tree.
#' @param k A numeric value representing the constant used to calculate the probability of a node to be chosen as a parent of a new node.
#' 
#'
#' @return  A square matrix representing the mutation relationships between the nodes in the phylogenetic tree. Each row corresponds to a node, and each column corresponds to a mutation. The value at the i-th row and j-th column is 1 if the i-th node has the j-th mutation, and 0 otherwise.
#' 
#' @examples
#' 
#' # Create a mutation matrix for a phylogenetic tree with 10 nodes and k = 2
#' B <- create_B(10, 2)
#' 
#' @export
create_B <- function(n, k) {
  B <- diag(n)
  muts <- sample(1:n, n)
  placed_muts <- rep(NA, n)
  mut_levels <- rep(NA, n)
  # Root node (first mutation)
  i <- 1
  placed_muts[i] <- muts[i]
  mut_levels[i] <- 0
  # Rest of mutations
  for (i in 2:length(muts)) {
    mut <- muts[i]
    # Calculate probability of each node
    probs <- k**mut_levels[!is.na(mut_levels)]
    parent <- sample.vec(placed_muts[!is.na(placed_muts)], 1, prob = probs)
    B[mut, ] <- B[parent, ]
    B[mut, mut] <- 1
    placed_muts[i] <- mut
    mut_levels[i] <- length(get_ascendants_idx(B, mut))
  }
  return(B)
}


#' Calculate tumor clone frequencies in samples
#'
#' This function calculates the frequencies of each clone in a set of samples, given the global clone proportions in the tumor and their spatial distribution.
#'
#' @param B A matrix representing the mutation relationships between the nodes in the phylogenetic tree (B matrix).
#' @param m An integer representing the number of samples taken from the tumor.
#' @param selection A character string representing the evolutionary mode the tumor follows. This should be either "positive" or "neutral".
#' @param n_cells An integer representing the number of cells sampled from the multinomial distribution. Default is 100.
#' @importFrom magrittr %>%
#'
#' @return A matrix where each row corresponds to a sample, and each column corresponds to a clone. The value at the i-th row and j-th column is the frequency of the j-th clone in the i-th sample.
#'
#' @examples
#' 
#' # Create random topology with 20 nodes and k = 3
#' B <- create_B(20, 3)
#' 
#' # Create U matrix with parameter m=4 and "positive" selection
#' U <- create_U(B = B, m = 4, selection = "positive")
#' @export
create_U <- function(B, m, selection, n_cells=100) {
  clone_proportions <- calc_clone_proportions(B, selection)
  # Place clones in 1D space
  clones_space <- place_clones_space(B)
  density_coords <- clones_space$spatial_coords
  x <- clones_space$x
  idx <- value <- proportion <- weights <- sampled_weight <- scaled_weight <- . <- NULL
  clone_order <- setdiff(colnames(density_coords), "idx") # for U column names
  cutpoints_idx <- seq(20, length(x)-20, length.out = m)
  cutpoints <- x[cutpoints_idx]
  U <- purrr::map(cutpoints, function(j) {
      reshape2::melt(density_coords, id = "idx") %>%
      dplyr::filter(idx == j) %>% 
      dplyr::mutate(weights = value/sum(value)) %>%
      dplyr::left_join(clone_proportions, by = c("variable" = "clone_idx")) %>% # add clone_proportions
      dplyr::mutate(scaled_weight = weights*proportion/sum(weights*proportion)) %>%
      dplyr::mutate(sampled_weight = stats::rmultinom(n = 1, size = n_cells, prob = scaled_weight)/n_cells) %>% # modelo multinomial de muestreo para conseguir clones con frecuencia 0
      dplyr::select(sampled_weight)
  }) %>%
    do.call(cbind, .) %>%
    magrittr::set_rownames(clone_order) %>%
    magrittr::set_colnames(as.character(cutpoints)) %>%
    t()
  # Sort columns from U matrix to coincide with rows in B. De todas maneras hacer esto o no hacerlo creo que da igual
  U <- U[, order(as.numeric(gsub("clon", "", colnames(U))))]
  # colSums(U)
  return(U)
}

#' Calculate the variant allele frequency (VAF) values in a set of samples
#'
#' This method generates the F matrix that contains the mutation frequency values of a series of mutations in a collection of tumor biopsies or samples.
#'
#' @param B A matrix representing the mutation relationships between the nodes in the phylogenetic tree.
#' @param U A matrix where each row corresponds to a sample, and each column corresponds to a clone. The value at the i-th row and j-th column is the frequency of the j-th clone in the i-th sample.
#' @param heterozygous A logical value indicating whether to adjust the clone proportions for heterozygous states. If `TRUE`, the clone proportions are halved. If `FALSE`, the clone proportions are not adjusted. Default is `TRUE`.
#'
#' @return A matrix containing the VAF values of a series of mutations in a set of samples.
#'
#' @examples
#' # Create random topology with 10 nodes and k = 2
#' B <- create_B(10, 2)
#' 
#' # Create U matrix with parameter m=4 and "positive" selection
#' U <- create_U(B = B, m = 4, selection = "positive")
#' 
#' # Then we compute the F matrix for a heterozygous tumor
#' F <- create_F(U = U, B = B, heterozygous = TRUE)
#'
#' @export
create_F <- function(U, B, heterozygous = TRUE) {
  F_homozygous <- U %*% B
  if (heterozygous) {
    F_ <- 0.5 * F_homozygous
  } else {
    F_ <- F_homozygous
  }
  # wrap values of F between 0 and 1 as they represent clone proportions
  return(apply(F_, c(1,2), FUN = function(x) min(x,1)))
}

#' Add noise to the VAF values in an F matrix
#'
#' This function adds noise to the variant allele frequency (VAF) values in an F matrix, simulating the effect of sequencing errors. The noise is modeled as a negative binomial distribution for the depth of the reads and a binomial distribution for both the variant allele counts and the mismatch counts.
#'
#' @param F_matrix A matrix representing the true VAF values of a series of mutations in a set of samples (F matrix).
#' @param depth A numeric value representing the mean depth of sequencing.
#' @param overdispersion A numeric value representing the overdispersion parameter for the negative binomial distribution used to simulate the depth of sequencing.
#'
#' @return A matrix containing noisy VAF values of a series of mutations in a set of samples.
#'
#' @examples
#' # Calculate the noisy VAF values of a series of mutations in a set of samples, given the true 
#' # VAF values in the F matrix F_true, a depth of 30 and an overdispersion of 5
#'
#' # Simulate the noise-free F matrix of a tumor with 50 clones,
#' # 10 samples, k = 5, following a positive selection model
#' F_true <- create_instance(
#'   n = 50,
#'   m = 10,
#'   k = 5,
#'   selection = "positive", 
#'   noisy = FALSE)$F_true
#' 
#' # Then we add the noise using a depth of 30 and an overdispersion of 5.
#' noisy_F <- add_noise(F_true, 30, 5)
#'
#' @export
add_noise <- function(F_matrix, depth, overdispersion) {
  #warn <- getOption("warn")
  #options(warn=-1)
  depth_values <- purrr::map_dbl(1:length(F_matrix), function(x) stats::rnbinom(n = 1, mu = depth, size = overdispersion))
  alt_allele_counts <- purrr::map_dbl(1:length(F_matrix), function(x) stats::rbinom(n = 1, prob = F_matrix[x], size = depth_values[x]))
  #ref_allele_counts <- map_dbl(1:length(F), function(x) rbinom(n = 1, prob = 1-F[x], size = depth_values[x]))
  alt_allele_mismatch_counts <- purrr::map_dbl(1:length(F_matrix), function(x) stats::rbinom(n = 1, prob = 0.001, size = alt_allele_counts[x]))
  ref_to_alt_allele_mismatch_counts <- purrr::map_dbl(1:length(F_matrix), function(x) stats::rbinom(n = 1, prob = 0.001/3, size = depth_values[x] - alt_allele_counts[x])) # aqui en el size no uso ref_allele_counts porque tengo en cuenta q en el sitio puede haber otro alelo y su cambio hacia el alternativo también nos estaría sumando a las reads alternativas
  noisy_vafs <- matrix(
    (alt_allele_counts - alt_allele_mismatch_counts + ref_to_alt_allele_mismatch_counts)/depth_values,
    byrow = FALSE, nrow = nrow(F_matrix))
  # For low depth values it can happen that we get depth values of 0s and so we get a NaN after dividing by 0 in the noisy VAF matrix. That would mean a VAF of 0. Hence, we'll convert those NaN values to 0s.
  noisy_vafs[is.na(noisy_vafs)] <- 0
  #options(warn=warn)
  return(noisy_vafs)
}

#' Create a tumor phylogenetic tree instance
#' @description This function generates a tumor phylogenetic tree instance, composed by a mutation matrix (B matrix), a matrix of true variant allele frequencies (F_true), a matrix of noisy variant allele frequencies (F), and a matrix of clone frequencies in samples (U).
#' @details The B matrix is a square matrix representing the mutation relationships between the clones in the tumor, or, in other words, it represents the topology of the phylogenetic tree. The F_true matrix represents the true variant allele frequencies of the mutations present in the tumor in a set of samples. The F matrix represents the noisy variant allele frequencies of the mutations in the same set of samples. The U matrix represents the frequencies of the clones in the tumor in the set of samples.
#'
#' @param n An integer representing the number of clones.
#' @param m An integer representing the number of samples.
#' @param k A numeric value that determines the linearity of the tree topology. Also referred to as the topology parameter. Increasing values of this parameter increase the linearity of the topology. When `k` is set to 1, all nodes have equal probabilities of being chosen as parents, resulting in a completely random topology.
#' @param selection A character string representing the evolutionary mode the tumor follows. This should be either "positive" or "neutral".
#' @param noisy A logical value indicating whether to add noise to the frequency matrix. If `TRUE`, noise is added to the frequency matrix. If `FALSE`, no noise is added. `TRUE` by default.
#' @param depth A numeric value representing the mean depth of sequencing. 30 by default.  
#' @param seed A numeric value used to set the seed for the random number generator. Sys.time() by default.
#' 
#' @return A list containing four elements: 'F', a matrix representing the noisy frequencies of each mutation in each sample; 'B', a matrix representing the mutation relationships between the clones in the tumor; 'U', a matrix that represents the frequencies of the clones in the tumor in the set of samples; and 'F_true', a matrix representing the true frequencies of each mutation in each sample.
#' 
#' @examples
#' # Create an instance of a tumor with 10 clones,
#' # 4 samples, k = 1, neutral evolution and
#' # added noise with depth = 500
#' I1 <- create_instance(
#'   n = 10,
#'   m = 4,
#'   k = 1,
#'   selection = "neutral",
#'   depth = 500)
#'   
#' 
#' # Create an instance of a tumor with 50 clones,
#' # 10 samples, k = 5, positive selection and
#' # added noise with depth = 500
#' I2 <- create_instance(
#'   n = 50,
#'   m = 10,
#'   k = 5,
#'   selection = "positive", 
#'   noisy = TRUE,
#'   depth = 500)
#'   
#'   
#' # Create an instance of a tumor with 100 clones,
#' # 25 samples, k = 0, positive selection without 
#' # noise
#' I3 <- create_instance(
#'   n = 100,
#'   m = 25,
#'   k = 0,
#'   selection = "positive", 
#'   noisy = FALSE)
#' @export
create_instance <- function(n, m, k, selection, noisy = TRUE, depth = 30, seed = Sys.time()) {
  if(!is.character(selection)){
    stop("\n selection must be a character")
  }
  if(selection!="neutral" & selection!="positive"){
    stop("\n selection must be neutral or positive")
  }
  if(depth < 1 || !identical(round(depth), depth)){
    stop("\n depth must be a natural number greater than or equal to 1 ")
  }
  set.seed(seed)
  overdispersion <- 5
  # Create random topology
  B <- create_B(n, k)
  # Assign proportions to each clone

  # Create U matrix
  U <- create_U(B = B, m = m, selection = selection)
  F_true <- create_F(U = U, B = B, heterozygous = FALSE)
  if (noisy) {
    F_matrix <- add_noise(F_matrix = F_true, depth = depth, overdispersion = overdispersion)
  } else {
    F_matrix <- F_true
  }
  clone_names <- purrr::map(1:n, function(x) paste0("clone", x))
  mutation_names <- purrr::map(1:n, function(x) paste0("mut", x))
  sample_names <- purrr::map(1:m, function(x) paste0("sample", x))
  rownames(F_matrix) <- sample_names
  colnames(F_matrix) <- mutation_names
  rownames(B) <- clone_names
  colnames(B) <- mutation_names
  rownames(U) <- sample_names
  colnames(U) <- clone_names
  rownames(F_true) <- sample_names
  colnames(F_true) <- mutation_names
  return(list(F_noisy = F_matrix, B = B, U = U, F_true = F_true))
}


