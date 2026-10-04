# Plot the TreeMix maximum likelihood tree with migration edges for three genetic
# P. alecto populations and P. conspicillatus, rooted with P. alecto alecto.
#
# Run from this directory:  cd 4_Treemix/scripts && Rscript treemix_plot.R
#
# Inputs:  ../data/m2_reps_NAUS_NQECOAST_INDOplusWAmuseum_SFF_IBFF/<prefix>.99.2.*
#          (TreeMix output of replicate 99 of 100 with m = 2 migration edges; one of
#          the 28 replicates tied for the highest likelihood, ln(L) = 112.664)
#          plotting_funcs.R (TreeMix plotting functions by Joe Pickrell)
# Output:  ../figures/Treemix_m2_NAUS_NQECOAST_INDOplusWAmuseum_SFF_IBFF.pdf   Figure S16A
#          (population labels were renamed to manuscript codes in a graphics editor)

library(RColorBrewer)
library(R.utils)

source("plotting_funcs.R")
prefix <- file.path("..", "data", "m2_reps_NAUS_NQECOAST_INDOplusWAmuseum_SFF_IBFF",
                    "bffsff_NAUS_NQECOAST_INDOplusWAmuseum_SFF_IBFF.99.2")

dir.create("../figures", showWarnings = FALSE)
pdf("../figures/Treemix_m2_NAUS_NQECOAST_INDOplusWAmuseum_SFF_IBFF.pdf")
plot_tree(prefix)
dev.off()
