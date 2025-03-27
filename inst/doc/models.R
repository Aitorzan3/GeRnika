## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo = FALSE, message = FALSE, warning = FALSE--------------------------
library(GeRnika)
library(ggpubr)
library(knitr)
set.seed(1)

## -----------------------------------------------------------------------------
I1 <- create_instance(n = 10, m = 4, k = 1, selection = "neutral", noisy = FALSE, seed = 1)
I1

## -----------------------------------------------------------------------------
I1 <- create_instance(n = 10, m = 4, k = 1, selection = "neutral", noisy = TRUE, 
                      depth = 100, seed = 2)
I1

## ----fig.align = 'center', out.width = "60%", fig.cap = " Function that determines the frequency of selection of an existing node in the tree $T$ as the parent node of new children. The function is dependant on the number of ascendants and the topology parameter $k$.", echo = FALSE----
knitr::include_graphics("figs/tumor_model_k.png")

## -----------------------------------------------------------------------------
n <- 10 
k <- 8

B_mat <- create_B(n, k)
B_mat

## ----fig.align = 'center', out.width = "60%", fig.cap = " Ternary density plots of 5,000 samples drawn from two 3-dimensional Dirichlet distributions. The parameters of the Dirichlet distribution on the left are $\\boldsymbol{\\alpha}$ = 0.3 and the distribution is used to represent positive-driven evolution. The distribution on the right has parameters $\\boldsymbol{\\alpha}$ = [5, 10, 10] and it is used to represent neutral evolution.", echo = FALSE----
knitr::include_graphics("figs/ternary_plots.png")

## ----fig.align = 'center', out.width = "80%", fig.cap = "Overlapping area of pairs of Gaussian probability density functions with $sd$=1 and mean differences ($d$) ranging between 0 and 4. The functions represent 1-dimensional projections of tumor clone masses.", echo = FALSE----
knitr::include_graphics("figs/clone_distance.png")

## -----------------------------------------------------------------------------
m <- 6  
selection <- "neutral"

U_mat <- create_U(B = B_mat, m = m, selection = selection, n_cells = 100)
U_mat

## -----------------------------------------------------------------------------
F_mat <- create_F(U = U_mat, B = B_mat, heterozygous = FALSE)

## ----fig.align = 'center', out.width = "65%", fig.cap = "Density of the mean absolute error in noisy $F$ matrices for different $\\mu_{sd}$ values that correspond to different noise levels.", echo = FALSE----
knitr::include_graphics("figs/noisy_F_error.png")

## ----fig.align = 'center', out.width = "100%", fig.cap = "Graphical summary of the models for building synthetic datasets for the CDEP problem, accompanied by an example. The simulation builds upon two core models and a third one which is optional. The tumor model first simulates a clonal tree $T$ (and an associated $\\boldsymbol{B}$ matrix) of $n$ nodes; this is done in an iterative manner and as a function of the parameter $k$, which controls for the linearity of the topology; when $k <$ 1, clones up in the tree are favoured over those deep in the topology for being the parents of the newly added nodes and when $k >$ 1, the deep ones are favoured over the others (Equation XX). In this example, $k$ was set to 0, which leads to a random topology. After that, we simulate the proportions of the clones in the tumor at the moment of sampling ($\\boldsymbol{c}$). For calculating them, a Dirichlet distribution is sampled in each multifurcation of $T$ and then the resulting values are scaled to the original proportion of the parent clone, so that at the end this results in all the values or proportions summing up to one. The proportion distribution is dependant on the evolution model we assume for the tumor (positive selection or neutral evolution) and this is controlled by using different sets of $\\alpha$ values in the Dirichlet distributions. As a result of the proportions being (scaled) samples from Dirichlet distributions, we have that each clone proportion $c_i$ follows the distributions shown in this figure. In this example we opted for a neutrally evolving tumor, which resulted in the $\\boldsymbol{c}$ values that have been overlaid over each clone in $T$. The last step in the tumor model is to simulate the tumor blend, which is modelled by a Gaussian mixture model of $n$ components, where each component represents a tumor clone. The distance between two successive components or clones $d \\in \\{0, 0.1, \\ldots,4\\}$ follows a Multinomial distribution where the probabilities are proportional to the density function of a Beta distribution. Thus, for small $d$ values, two clones are completely mixed, whereas for values clones to 4, they are physically far from each other. The sampling model draws $m$ samples from the Gaussian mixture model just described. These sampling sites correspond to the $m$ cutpoints that divide the domain of the model into $m +$ 1 equal-sized bins. In this example, we have got 4 samples. The probability density of the components of the mixture model on those sites is normalized so that their sum sums up to 1 and the resulting values are used as the probability values of a Multinomial distribution that is used to calculate the final values in the matrix $\\boldsymbol{U}$. Once we have got the $\\boldsymbol{U}$ and $\\boldsymbol{B}$ matrices, we can calculate the $\\boldsymbol{F}$ matrix using Equation \ref{eq:F}. Finally, we can optionally add noise that simulates sequencing noise to this newly created $\\boldsymbol{F}$ matrix and create a new noisy $\\boldsymbol{F}^n$ matrix. This model simulates noise at the level of the sequencing reads by sampling various distributions that simulate the total number of reads at the site, the number of reads that support each of the mutations and mismatch errors, and it then recalculates the new $f_{ij}^n$ values based on them.", echo = FALSE----
knitr::include_graphics("figs/model_summary_goodqual.jpg")

