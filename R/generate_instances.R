source("R/utils.R")

sample.vec <- function(x, ...) x[sample(length(x), ...)]

create_B <- function(n, k) {
  B <- diag(n)
  muts <- sample(1:n, n)
  placed_muts <- rep(NA, n)
  mut_levels <- rep(NA, n)
  # Root node does not need to be touched
  i <- 1
  placed_muts[i] <- muts[i]
  mut_levels[i] <- 0
  # Rest of mutations
  for (i in 2:length(muts)) {
    mut <- muts[i]
    # Calculate probability of each node
    probs <- k**mut_levels[!is.na(mut_levels)]
    parent <- sample.vec(placed_muts[!is.na(placed_muts)], 1, prob = probs)
    B[mut, ] <- B[parent, ] # como no sabemos la topologia es necesario rellenar la B rowwise
    B[mut, mut] <- 1
    placed_muts[i] <- mut
    mut_levels[i] <- length(get_ascendants_idx(B, mut))
  }
  return(B)
}

.distribute_freqs <- function(B, clone_idx, clone_proportions, selection) {
  dirich_params <- tibble(positive = rep(0.3, 2),
                          neutral = c(5, 10))
  children_clone_idx <- get_children_idx(B = B, node_idx = clone_idx)
  if (!is.null(children_clone_idx)) {
    params <- pull(dirich_params, eval(selection))
    dirich_values <- rdirichlet(1, c(params[1],
                                     rep(params[2], length(children_clone_idx))))
    # Normalize the values to the father clone's proportion
    norm_dirich_values <- dirich_values*clone_proportions[clone_idx]
    clone_proportions[c(clone_idx, children_clone_idx)] <- norm_dirich_values # el clon i esta en la posicion i del vector clone proportions
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
  clone_proportion_df <- tibble(clone_idx = paste0("clon", 1:n), proportion = clone_proportions)
  # clon_i corresponds to clone containing mutation i for the first time and not necessarility the clone in row i, but in this case it is also the clone in row i because of the way we construct the B matrix
  return(clone_proportion_df)
}

place_clones_space <- function(B) {
  n <- nrow(B)
  max_sep <- 4
  mean_diffs <- seq(from = 0, to = max_sep, by = 0.1)
  clone_order <- sample(1:n, n)
  clone_mean_diff <- sample(x = mean_diffs,
                            size = n-1,
                            prob = dbeta(x = seq(0, 1, length.out = length(mean_diffs)), shape1 = 1, shape2 = 5), replace = TRUE) # change beta distribution parameters for tumor density
  clone_means <- c(0, cumsum(clone_mean_diff))
  clone_sd <- 1
  x <- seq(min(clone_means - 3*clone_sd), max(clone_means + 3*clone_sd), .01)
  y <- map(1:n, function(z) dnorm(x, mean = clone_means[z], sd = clone_sd)) %>%
    do.call(cbind, .) %>%
    as_tibble(.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)) %>%
    magrittr::set_colnames(paste0("clon", clone_order))
  spatial_coords <- y %>%
    mutate(idx = x)
  return(list(spatial_coords = spatial_coords, x = x))
}

create_U <- function(B, clone_proportions, density_coords, m, x) {
  n_cells <- 100
  clone_order <- setdiff(colnames(density_coords), "idx") # for U column names
  # Uniform sampling (grid). No nos vamos hasta los extremos sino que dejamos las "esquinas" fuera (20 points)
  cutpoints_idx <- seq(20, length(x)-20, length.out = m)
  cutpoints <- x[cutpoints_idx]
  U <- map(cutpoints, function(j) {
      reshape2::melt(density_coords, id = "idx") %>%
      filter(idx == j) %>% 
      mutate(weights = value/sum(value)) %>%
      left_join(clone_proportions, by = c("variable" = "clone_idx")) %>% # add clone_proportions
      mutate(scaled_weight = weights*proportion/sum(weights*proportion)) %>%
      mutate(sampled_weight = rmultinom(n = 1, size = n_cells, prob = scaled_weight)/n_cells) %>% # modelo multinomial de muestreo para conseguir clones con frecuencia 0
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

add_noise <- function(F, depth, overdispersion) {
  warn <- getOption("warn")
  #options(warn=-1)
  depth_values <- map_dbl(1:length(F), function(x) rnbinom(n = 1, mu = depth, size = overdispersion))
  alt_allele_counts <- map_dbl(1:length(F), function(x) rbinom(n = 1, prob = F[x], size = depth_values[x]))
  #ref_allele_counts <- map_dbl(1:length(F), function(x) rbinom(n = 1, prob = 1-F[x], size = depth_values[x]))
  alt_allele_mismatch_counts <- map_dbl(1:length(F), function(x) rbinom(n = 1, prob = 0.001, size = alt_allele_counts[x]))
  ref_to_alt_allele_mismatch_counts <- map_dbl(1:length(F), function(x) rbinom(n = 1, prob = 0.001/3, size = depth_values[x] - alt_allele_counts[x])) # aqui en el size no uso ref_allele_counts porque tengo en cuenta q en el sitio puede haber otro alelo y su cambio hacia el alternativo también nos estaría sumando a las reads alternativas
  noisy_vafs <- matrix(
    (alt_allele_counts - alt_allele_mismatch_counts + ref_to_alt_allele_mismatch_counts)/depth_values,
    byrow = FALSE, nrow = nrow(F))
  # For low depth values it can happen that we get depth values of 0s and so we get a NaN after dividing by 0 in the noisy VAF matrix. That would mean a VAF of 0. Hence, we'll convert those NaN values to 0s.
  noisy_vafs[is.na(noisy_vafs)] <- 0
  #options(warn=warn)
  return(noisy_vafs)
}

#' @export
#' @title Simulate tumor data
#' @description Simulates a tumor instance, composed by \code{F}, \code{F_true}, \code{B} and \code{U}.
#'
#' @param n the number of clones.
#' @param m the number of samples.
#' @param k continuous number how branchy the created topology is 
#' @param selection character that specifies the clone selection. Possible values: \code{"positive"} and \code{"neutral"}
#' @param noisy optional logical that specifies whether noise is added to values in \code{F} or not. FALSE by default.
#' @param depth optional argument representing the read sequencing depth (for noisy cases). 30 by default.   
#' @return the instance of a tumor sample, composed by \code{F}, \code{F_true}, \code{B} and \code{U} .
#' @examples
#' # Create an instance composed by 10 clones,
#' # 4 samples, k = 1, "neutral" selection and
#' # with added noise and depth = 500
#' I1 <- create_instance(
#'   n = 10,
#'   m = 4,
#'   k = 1,
#'   selection = "neutral",
#'   depth = 500)
#'   
#' 
#' # Create an instance composed by 50 clones,
#' # 10 samples, k = 5, "positive" selection with 
#' # added noise and depth = 500
#' I2 <- create_instance(
#'   n = 50,
#'   m = 10,
#'   k = 5,
#'   selection = "positive", 
#'   noisy = TRUE,
#'   depth = 500)
#'   
#'   
#' # Create an instance composed by 100 clones,
#' # 25 samples, k = 0, "positive" selection without 
#' # added noise
#' I3 <- create_instance(
#'   n = 100,
#'   m = 25,
#'   k = 0,
#'   selection = "positive", 
#'   noisy = FALSE)
create_instance <- function(n, m, k, selection, noisy = TRUE, depth = 30, seed = Sys.time()) {
  if(class(selection)!="character"){
    stop("\n selection must be a character")
  }
  if(selection!="neutral" & selection!="positive"){
    stop("\n selection must be neutral or positive")
  }
  set.seed(seed)
  # overdispersion de la binomial negativa nos da la diferencia de probabilidad del evento o su negado
  overdispersion <- 5
  # Firstly, create random topology
  B <- create_B(n, k)
  # Assign proportions to each clone
  clone_proportions <- calc_clone_proportions(B, selection)
  # Place clones in 1D space
  clones_space <- place_clones_space(B)
  density_coords <- clones_space$spatial_coords
  domain <- clones_space$x
  # Create U matrix
  U <- create_U(B = B, clone_proportions = clone_proportions, density_coords = density_coords, m = m, x = domain)
  F_true <- calc_F(U = U, B = B, heterozygous = FALSE)
  # Coarse 1s that have decimals different to 0 in the 16th position to 1
  # F_true[round(F_true, digits = 15) == 1] <- 1
  if (noisy) {
    F <- add_noise(F = F_true, depth = depth, overdispersion = overdispersion)
  } else {
    F <- F_true
  }
  clone_names<-map(1:n, function(x) paste0("clone", x))
  mutation_names<-map(1:n, function(x) paste0("mut", x))
  sample_names<-map(1:m, function(x) paste0("sample", x))
  rownames(F)<- sample_names
  colnames(F)<- mutation_names
  rownames(B)<- clone_names
  colnames(B)<- mutation_names
  rownames(U)<- sample_names
  colnames(U)<- clone_names
  rownames(F_true)<- sample_names
  colnames(F_true)<- mutation_names
  return(list(F = F, B = B, U = U, F_true = F_true))
}


