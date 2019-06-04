# Sheng's test file
source("functions.R")

# use functions to generate the plot
setwd("~/ontologyPCM/R")

taxon <- "Characiformes"
trait <- "dermal bone"

td <- Get_Tree_Data(taxon, trait)

#tdf <- filter_coverage(td, traits=0, taxa=0.1)
#tdf <- filter_coverage(tdf, traits=0.1, taxa=0)
#tdf <- filter_coverage(tdf, traits=0, taxa=0.2)
#tdf <- filter_coverage(tdf, traits=0.2, taxa=0)
#tdf <- filter_coverage(tdf, traits=0, taxa=0.5)
#tdf <- filter_coverage(tdf, traits=0.5, taxa=0)

# should not return error
td <- filter_coverage(td, traits=0.01, taxa=0)

# returns a tree with one trait
td <- filter_coverage(td, traits=0.5, taxa=0)

# should return an error
td <- filter_coverage(td, traits=0.7, taxa=0)


# make tree if traits & taxa are both > 1
if ((length(td$dat) > 1) && length(td$phy) > 1){
  njt <- makeTree(td)
}
pdf(paste("~/ontologyPCM/", taxon, "_", trait, ".pdf"), width=10, height=8)
plotData(td, njt=NULL, show.tip.label=TRUE, cex=0.25)
dev.off()
