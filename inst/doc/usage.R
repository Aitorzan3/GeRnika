## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo = FALSE, message = FALSE, warning = FALSE--------------------------
library(GeRnika)
library(ggpubr)

## ----message = TRUE-----------------------------------------------------------
I <- create_instance(n = 5, m = 4, k = 0.5, selection = "neutral", seed = 1)

## ---- message = FALSE, warning = FALSE,  out.width="8.25in", fig.align ="center", fig.dim = c(6.5,6), fig.cap="On the left the plot of `tree1` (`k=0`), on the right the plot of `tree2` (`k=8`)."----
# Simulate a tumor with k=0:
I1 <- create_instance(n = 5, m = 4, k = 0, selection = "neutral", seed = 1)

# Simulate a tumor with k=1:
I2 <- create_instance(n = 5, m = 4, k = 8, selection = "neutral", seed = 1)

# Create a `Phylotree` class object for each tumor:
tree1 <- B_to_phylotree(B = I1$B)
tree2 <- B_to_phylotree(B = I2$B)

# Plot both trees
DiagrammeR::render_graph(DiagrammeR::combine_graphs(data.tree::ToDiagrammeRGraph(tree1@tree), data.tree::ToDiagrammeRGraph(tree2@tree)))

## ---- message = FALSE, echo = TRUE,  warning = FALSE, fig.align ="center", fig.dim = c(7.5,8), fig.cap="Heatmaps of the $\boldsymbol{U}$ matrices of an instance of a tumor under positive selection (top) and neutral evolution (bottom)."----
# Function to create the heatmap of the U matrix
U_to_heatmap <- function(U, values = TRUE, col_names = c("samples", "clones", "proportion")){
  Upos <-reshape2::melt(U)
  colnames(Upos) <- col_names
  Var1 <- col_names[1]
  Upos<-ggplot(Upos, aes(x = samples, y = clones, fill = proportion)) + geom_tile(col = "black") +    
    theme(
         panel.background = element_rect(fill = 'transparent'),
         plot.background = element_rect(fill = 'transparent', color = NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.background = element_rect(fill = 'transparent'),
         legend.box.background = element_rect(fill = 'transparent')
    ) + scale_fill_gradient(limits = c(0.0000000000001, 1)) + theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())
  if(values){
    Upos <- Upos + geom_text(aes(label = proportion), size = 4)  
  }
  Upos
}

# simulate a tumor with neutral selection:
Ipos <- create_instance(n = 5, m = 8, k = 0.5, selection = "positive", seed = 1)

# simulate a tumor with positive selection:
Ineu <- create_instance(n = 5, m = 8, k = 0.5, selection = "neutral", seed = 1)

Upos <- U_to_heatmap(Ipos$U)
Uneu <- U_to_heatmap(Ineu$U)

ggarrange(plotlist = list(Upos, Uneu), ncol = 1, nrow = 2)

## ---- message = FALSE, echo = TRUE, warning = FALSE, fig.align ="center", fig.dim = c(7.5, 4), fig.cap="The effect of noise."----
# Function to create a heatmap of F
F_to_heatmap <- function(U, values = TRUE, col_names = c("samples", "mutations", "VAF")){
  Upos <-reshape2::melt(U)
  colnames(Upos) <- col_names
  Var1 <- col_names[1]
  Upos<-ggplot(Upos, aes(x = samples, y = mutations, fill = VAF)) + geom_tile(col = "black") +    
    theme(
         panel.background = element_rect(fill = 'transparent'),
         plot.background = element_rect(fill = 'transparent', color = NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.background = element_rect(fill = 'transparent'),
         legend.box.background = element_rect(fill = 'transparent')
    ) + scale_fill_gradient(limits = c(0.00000000000000000001, 1)) + theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) 
  if (values) {
    Upos <- Upos + geom_text(aes(label = round(VAF, digits = 2)), size = 4)  
  }
  Upos
}

# Simulate a tumor with sequencing noise added:
Inoisy <- create_instance(m = 5, n = 8, k = 0.5, selection = "neutral", noisy = TRUE, depth = 5, seed = 1)

Fnoise <- F_to_heatmap(abs(Inoisy$F_noisy - Inoisy$F_true))

ggarrange(Fnoise)

## ---- eval = TRUE-------------------------------------------------------------
# Simulate a tumor with 5 clones, 4 samples, k=0.5:
instance <- create_instance(n = 5, m = 4, k = 0.5, selection = "neutral", noisy = FALSE, seed = 1)

# The creation of the Phylotree class object using the previously generated B matrix:
phylotree <- B_to_phylotree(B = instance$B)

## ---- out.width="8.5in", fig.align="center", fig.cap="Phylogenetic tree composed by 5 nodes."----
plot(phylotree)

## -----------------------------------------------------------------------------
# Create a vector with the tags we want to use (the length of the 
# vector must be equal to the number of clones in the phylogenetic tree):
tags <- c("mut1", "mut2", "mut3", "mut4", "mut5")

# Create the Phylotree class object and include the node labels:
phylotree <- B_to_phylotree(B = instance$B, labels = tags)

## ----  out.width="8.5in", fig.align="center",fig.cap="Phylogenetic tree of 5 clones with custom tags. The tags in the vector are assigned to the nodes in the tree according to their order. For instance, the first clone in the previous plot is now represented with the first label of the tag vector `mut1`."----

plot(phylotree, labels = TRUE)

## ---- eval=TRUE,  out.width="9.5in", fig.align='center', fig.dim = c(7.5,5), fig.cap="Phylogenetic tree of 5 clones using proportions. In this case, there are 4 plots because we have generated an instance based on 4 samples. Then, each tree represents the proportions of the clones in each sample."----
# Simulate a tumor
instance <- create_instance(n = 5, m = 4, k = 0.5, selection = "neutral", noisy = FALSE, seed = 1)

# Use the B matrix to instantiate a Phylotree class object
tree <- B_to_phylotree(B = instance$B)

# Plot the phylogenetic tree and resize the nodes according to the U matrix
plot_proportions(tree, instance$U)

## ---- eval=TRUE,  out.width="9.5in", fig.align='center', fig.dim = c(7.5,5), fig.cap="Phylogenetic tree of 5 clones using proportions and tags."----
tags <- c("A", "B", "C", "D", "E")

# Simulate a tumor
instance <- create_instance(n = 5, m = 4, k = 0.5, selection = "neutral", noisy = FALSE, seed = 1)

# Use the B matrix of the previously generated tumor for creating a Phylotree class object
tree<-B_to_phylotree(B = instance$B, labels = tags)

# Plot the phylogenetic tree, resize the nodes according to the U matrix and include the node labels
plot_proportions(tree, instance$U, labels = TRUE)

## ----fig.align ="center", out.width="9.5in", fig.align='center', echo = TRUE, fig.dim = c(8,8.75)----
# Load the predefined B matrices of the package:
B_mats <- GeRnika::B_mats

# Get one of the B matrix trios to compare them:
B_real <- B_mats[[2]]$B_real
B_alg1 <- B_mats[[2]]$B_alg1
B_alg2 <- B_mats[[2]]$B_alg2

#Create the list of the tags for the clones that compose the phylogenetic trees:
tags <- c("mut1", "mut2", "mut3", "mut4", "mut5", "mut6", 
            "mut7", "mut8", "mut9", "mut10")

#Create the Phylotree class objects using the previously loaded B matrices:
phylotree_real <- B_to_phylotree(B_real, labels = tags)
phylotree_alg1 <- B_to_phylotree(B_alg1, labels = tags)
phylotree_alg2 <- B_to_phylotree(B=B_alg2, labels=tags)


#Render both trees:
DiagrammeR::render_graph(DiagrammeR::combine_graphs(data.tree::ToDiagrammeRGraph(phylotree_real@tree),DiagrammeR::combine_graphs(data.tree::ToDiagrammeRGraph(phylotree_alg2@tree),data.tree::ToDiagrammeRGraph(phylotree_alg1@tree))))

## -----------------------------------------------------------------------------
# Checking if phylotree_real is equal to itself:
equals(phylotree_1 = phylotree_real, phylotree_2 = phylotree_real)

# Checking if phylotree_real and phylotree_alg1 are equal:
equals(phylotree_1 = phylotree_real, phylotree_2 = phylotree_alg1)

## ----fig.align ="center",  out.width="9.5in", fig.align='center'--------------
find_common_subtrees(phylotree_1 = phylotree_real, phylotree_2 = phylotree_alg2)

## ----fig.align ="center", fig.dim = c(6,4), out.width="9.5in", fig.align='center'----
find_common_subtrees(phylotree_1 = phylotree_real, phylotree_2 = phylotree_alg1)

## ----out.width="9.5in", fig.align='center',fig.dim = c(6,4)-------------------
find_common_subtrees(phylotree_1 = phylotree_real, phylotree_2 = phylotree_alg2, labels = TRUE)

## ----out.width="9.5in", fig.align='center', fig.dim = c(6,4)------------------
find_common_subtrees(phylotree_1 = phylotree_real, phylotree_2 = phylotree_alg1, labels = TRUE)

## ----fig.align ="center", fig.dim = c(8.5,9.5), out.width = "7.5in", warning=FALSE----
# Creating the consensus tree between phylotree_real and phylotree_alg2
consensus_real_alg2 <- combine_trees(phylotree_1 = phylotree_real, phylotree_2 = phylotree_alg2)


# Creating the consensus tree between phylotree_real and phylotree_alg1
consensus_real_alg1 <- combine_trees(phylotree_1 = phylotree_real, phylotree_2 = phylotree_alg1)



## ----fig.align ="center", fig.dim = c(7,8),out.width="9.5in", fig.align='center', warning=FALSE----
# Rendering the consensus between phylotree_real and phylotree_alg2
DiagrammeR::render_graph(consensus_real_alg2)


## ----  fig.dim = c(5,6), out.width="9.5in", fig.align='center', warning=FALSE----
# Rendering the consensus between phylotree_real and phylotree_alg1
DiagrammeR::render_graph(consensus_real_alg1)

## ----fig.align ="center", fig.height = 5, fig.width = 6, out.width="9.5in", fig.align='center'----
# Load one of the default palettes of the package:
palette <- GeRnika::palettes$Lancet

# Create the consensus tree between phylotree_real and phylotree_alg1 
# using clone tags and the custom palette:
consensus <- combine_trees(phylotree_1 = phylotree_real, phylotree_2 = phylotree_alg1, 
                           labels = TRUE, palette = palette)

# Render the new consensus phylogenetic tree:
DiagrammeR::render_graph(consensus)

## ---- echo = TRUE-------------------------------------------------------------
sessionInfo(package = NULL)

