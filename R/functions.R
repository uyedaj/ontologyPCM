library(rphenoscape)
library(treeplyr)
library(RCurl)
library(rjson)
library(readr)
library(urltools)
library(pracma)
library(httr)
library(devtools)
library(phytools)
library(phylogram)
library(geiger)

#this function takes input (a taxon name and an entity name ) and then searches through the phenoscape database to locate the relevant
# and then returns a matrix of the relevant data called td 

Get_Tree_Data <- function (taxonName, entityName)
{  
  nex <- pk_get_ontotrace_xml(taxon = c(taxonName), entity = entityName)
  
  m <- pk_get_ontotrace(nex)
  m$taxa <- gsub(" ", "_", m$taxa)
  write.csv(m, file="~/repos/ontologyPCM/data/Ontotrace_Siluriformes_AnatomicalEntity.csv")
  
  tree <- read.tree("~/repos/ontologyPCM/data/actinopt_12k_treePL.tre")
  td <- make.treedata(tree, m)
  
  #saveRDS(td, "~/repos/ontologyPCM/data/tdSiluriformesAnatomicalEntity.rds")
  
  
  
  #td <- readRDS("~/repos/ontologyPCM/data/tdSiluriformesAnatomicalEntity.rds")
  #td
  
  return (td) 
  
}





# this function takes in our matrix td from the previous function and makes a tree out the matrix data
makeTree <- function (td)
{
  
  
  traits <- colnames(td$dat)
  traits <- traits[-(1:2)] #delete otu data
  traits
  
  
  
  #Get IRI ids for each trait.
  
  
  traitDetails <- lapply(traits, function(x) pk_anatomical_detail(x, verbose=TRUE))
  
  
  
  traitDetails[1:5]
  traitIDs <- unname(do.call(c, sapply(traitDetails, function(x) x[,'@id'])))
  
  
  
  irisPhenotypes <- sapply(traitIDs, url_encode)
  
  
  
  filename <- "a.txt"
  
  cat("iris=%5B%0A%20%20", file=filename)
  irisPhenotypes <- lapply(irisPhenotypes, function(x) gsub("/", "%2F", x, fixed=TRUE))
  irisPhenotypes <- lapply(irisPhenotypes, function(x) gsub(":", "%3A", x, fixed=TRUE))
  irisPhenotypes <- lapply(irisPhenotypes, function(x) gsub("=", "%3D%0A", x, fixed=TRUE))
  dum <- lapply(irisPhenotypes[1:(length(irisPhenotypes)-1)],function(x) cat(paste0('%22', x,'%22%2C', sep=""), file=filename, append=TRUE))
  cat(paste0('%22', irisPhenotypes[[length(irisPhenotypes)]],'%22',"%5D%0A", sep=""), file=filename, append=TRUE)
  
  
  
  
  
  api.semanticSimilarity_query <- "curl -X POST -d @a.txt 'https://kb.phenoscape.org/api/similarity/jaccard'"
  semanticSimilarityAPIResults <- system(api.semanticSimilarity_query, intern=TRUE)
  
  
  
  
  
  results <- fromJSON(semanticSimilarityAPIResults)
  scores <- lapply(results$results, function(x) x$score)
  scores <- sapply(scores, function(x) if(is.null(x)) NA else(x))
  result_terms <- do.call(rbind, lapply(results$results, function(x) do.call(cbind, lapply(x$terms, curlUnescape))))
  semanticSimilarityMatrix <- matrix(NA, nrow=length(irisPhenotypes), ncol=length(irisPhenotypes))
  diag(semanticSimilarityMatrix) <- 1
  rownames(semanticSimilarityMatrix) <- colnames(semanticSimilarityMatrix) <- curlUnescape(irisPhenotypes)
  
  for(i in 1:nrow(result_terms)){
    semanticSimilarityMatrix[result_terms[i,1], result_terms[i,2]] <- semanticSimilarityMatrix[result_terms[i,2], result_terms[i,1]] <- scores[i]
  }
  
  rownames(semanticSimilarityMatrix) <- colnames(semanticSimilarityMatrix) <- traits
  
  write.csv(semanticSimilarityMatrix, file="siluriformesSemanticSimMatrix.csv")
  
  
  
  #Check to see if semantic similarity matrix makes sense.
  
  
  maxSS <- list()
  for(i in 1:ncol(semanticSimilarityMatrix)){
    ss <- semanticSimilarityMatrix[-i,]
    j <- which(ss[,i]==max(ss[,i]))[1]
    maxSS[[i]] <- cbind(colnames(semanticSimilarityMatrix)[i], rownames(ss)[j], round(ss[j, i],4))
  }
  maxSS <- do.call(rbind, maxSS)
  as.data.frame(maxSS)
  
  
  minSS <- list()
  for(i in 1:ncol(semanticSimilarityMatrix)){
    ss <- semanticSimilarityMatrix[-i,]
    j <- which(ss[,i]==min(ss[,i]))[1]
    minSS[[i]] <- cbind(colnames(semanticSimilarityMatrix)[i], rownames(ss)[j], round(ss[j, i],4))
  }
  minSS <- do.call(rbind, minSS)
  as.data.frame(minSS)
  
  
  
  #Neighbor-joining tree of SS matrix
  
  njt <- nj(1-semanticSimilarityMatrix)
  pdf("njTreeSiluriformesSemanticMatrix.pdf", height=30, width=30)
  plot(njt, type="unrooted", cex=0.35)
  dev.off()
  
  return (njt)
}



# this function plots the data from the matrix as well as the tree we made in the previous function and outputs an image with the two trees 
#and a heatmap of the data 
plotData <- function(td, njt)
{
  X <- do.call(cbind, lapply(3:ncol(td$dat), function(x) as.numeric(td$dat[[x]])))
  colnames(X) <- colnames(td$dat)[3:ncol(td$dat)]
  
  
  tree <- njt
  tree1 <- chronopl(td$phy,1)
  tree2 <- chronopl(njt, 1)
  
  #alters the length of the tips of the trees
  tree1$edge.length <- tree1$edge.length/(max(branching.times(tree1)))*20
  tree2$edge.length <- tree2$edge.length/(max(branching.times(tree2)))*200
  
  #Changes the direction of the top plot
  h1 <- plot(tree1, plot = FALSE)
  h2 <- plot(tree2, plot = FALSE, direction = "downwards")
  
  # this is all setting boundaries for the different plots and combining them into one image
  par(mar = c(0,0,0,0))
  plot(0,0, type = 'n', xlim = c(0,h1$x.lim[2]+h2$x.lim[2]), ylim=c(0,h1$y.lim[2]+h2$y.lim[2]))
  
  image(seq(h1$x.lim[2],h1$x.lim[2]+h2$x.lim[2], length.out=nrow(X)), seq(0, h1$y.lim[2], length.out=ncol(X)), X,xlim=c(h1$x.lim[2],h1$x.lim[2]+h2$x.lim[2]) ,ylim=c(0, h1$y.lim[2]), add=TRUE)
  
  par(new = TRUE)
  
  
  plot(tree1, x.lim=c(2,(h1$x.lim[2]-.5)), y.lim=c(0,h1$y.lim[2]+h2$y.lim[2]), show.tip.label=FALSE)
  
  par(new = TRUE)
  plot(tree2, direction = "downwards", x.lim=c(-h1$x.lim[2],h2$x.lim[2]), y.lim=c((-h1$y.lim[2]/2),h2$y.lim[2]), show.tip.label=FALSE)
  
  return ()
  
}




