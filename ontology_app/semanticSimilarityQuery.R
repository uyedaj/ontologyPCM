setwd("~/repos/ontologyPCM/R/")
devtools::install_github("phenoscape/rphenoscape")
library(rphenoscape)
library(RCurl)
library(rjson)
library(readr)

iri <- rphenoscape::pk_get_iri("pectoral fin spine", as="uberon")


library(RCurl)
library(urltools)
library(pracma)
library(httr)

#jb <- scan("../formdata_JB.txt", what="character")
#jbdecode <- url_decode(url_decode(url_decode(jb)))
#jbrecode <- gsub("\n", "", jbdecode, fixed=TRUE)
#jbrecode <- gsub(" ", "", jbrecode, fixed=TRUE)
#jbrecode <- curlEscape(jbrecode)
#jbrecode <- gsub("iris%3D", "iris=", jbrecode)

#jbrecode <- gsub("%3C", "%253C", jbrecode)
#jbrecode <- gsub("%3A", "%253A", jbrecode)
#jbrecode <- gsub("%2F", "%252F", jbrecode)
#jbrecode <- gsub("%3F", "%253F", jbrecode)
#rbind(substr(jb, 1, 98), substr(jbrecode, 1, 98))
#cat(jbrecode, file="../formdata.txt")

#substr(jbrecode, 301, 400)


URLencode(iri, reserved=TRUE)
api.phenotype_query <- "http://kb.phenoscape.org/api/phenotype/query?entity=IRICODE&historical_homologs=false&limit=0&offset=0&parts=false&serial_homologs=false"
queryURL <- gsub("IRICODE", URLencode(iri, reserved=TRUE), api.phenotype_query)
phenotypesJSON <- httpPOST(queryURL)
phenotypes <- rjson::fromJSON(phenotypesJSON)
phenotype_labels <- sapply(phenotypes$results, function(x) x$label)
irisPhenotypes <- lapply(phenotypes$results, function(x) curlUnescape(x$'@id'))
irisValues <- lapply(irisPhenotypes, function(x) gsub("http://purl.org/phenoscape/expression?value=<", "", x, fixed=TRUE))
#irisPhenotypes <- lapply(1:length(irisPhenotypes), function(x) paste0(irisPhenotypes[[x]], ",\n"))
#irisPhenotypes[[length(irisPhenotypes)]] <- gsub(",\n", "]", irisPhenotypes[[length(irisPhenotypes)]])
irisValues <- lapply(irisValues, curlEscape)
irisValues <- lapply(irisValues, function(x) gsub("%2B", "+", x))
cat("iris=[", file="../formdata.txt")
lapply(irisValues,function(x) cat(paste0('  "http://purl.org/phenoscape/expression?value=%3C',x,'",\n', sep=""), file="../formdata.txt", append=TRUE))
tmp <- scan("../formdata.txt", what="character")

semanticIris <- curlEscape((readr::read_file("../formdata.txt")))
#semanticIris <- gsub("%0A", "", semanticIris)
#semanticIris <- gsub("%20", "", semanticIris)
#semanticIris <- gsub("%0A", "%250A", semanticIris)
semanticIris <- sub("%2C%0A$", "%0A%5D%0A", semanticIris)
semanticIris <- gsub("iris%3D%5B", "iris=%5B%0A", semanticIris)
semanticIris
cat(semanticIris, file="../formdata2.txt")

setwd("~/repos/ontologyPCM")


api.semanticSimilarity_query <- "curl -X POST -d @formdata2.txt 'http://kb.phenoscape.org/api/similarity/jaccard'"
semanticSimilarityAPIResults <- system(api.semanticSimilarity_query, intern=TRUE)
results <- fromJSON(semanticSimilarityAPIResults)
scores <- lapply(results$results, function(x) x$score)
scores <- sapply(scores, function(x) if(is.null(x)) NA else(x))
result_terms <- do.call(rbind, lapply(results$results, function(x) do.call(cbind, lapply(x$terms, curlUnescape))))
semanticSimilarityMatrix <- matrix(NA, nrow=length(phenotype_labels), ncol=length(phenotype_labels))
diag(semanticSimilarityMatrix) <- 1
rownames(semanticSimilarityMatrix) <- colnames(semanticSimilarityMatrix) <- irisPhenotypes

for(i in 1:nrow(result_terms)){
  semanticSimilarityMatrix[result_terms[i,1], result_terms[i,2]] <- semanticSimilarityMatrix[result_terms[i,2], result_terms[i,1]] <- scores[i]
}
rownames(semanticSimilarityMatrix) <- colnames(semanticSimilarityMatrix) <- phenotype_labels
write.csv(semanticSimilarityMatrix, "semanticSimilarityMatrix.csv")


pk_get_semantic_matrix_iri <- function(iri){
  URLencode(iri, reserved=TRUE)
  api.phenotype_query <- "http://kb.phenoscape.org/api/phenotype/query?entity=IRICODE&historical_homologs=false&limit=0&offset=0&parts=false&serial_homologs=false"
  queryURL <- gsub("IRICODE", URLencode(iri, reserved=TRUE), api.phenotype_query)
  phenotypesJSON <- httpPOST(queryURL)
  phenotypes <- rjson::fromJSON(phenotypesJSON)
  phenotype_labels <- sapply(phenotypes$results, function(x) x$label)
  irisPhenotypes <- lapply(phenotypes$results, function(x) curlUnescape(x$'@id'))
  irisValues <- lapply(irisPhenotypes, function(x) gsub("http://purl.org/phenoscape/expression?value=<", "", x, fixed=TRUE))
  irisValues <- lapply(irisValues, curlEscape)
  irisValues <- lapply(irisValues, function(x) gsub("%2B", "+", x))
  cat("iris=[", file="../formdata.txt")
  dum <- lapply(irisValues,function(x) cat(paste0('  "http://purl.org/phenoscape/expression?value=%3C',x,'",\n', sep=""), file="../formdata.txt", append=TRUE))
  semanticIris <- curlEscape((readr::read_file("../formdata.txt")))
  #semanticIris <- gsub("%0A", "", semanticIris)
  #semanticIris <- gsub("%20", "", semanticIris)
  #semanticIris <- gsub("%0A", "%250A", semanticIris)
  semanticIris <- sub("%2C%0A$", "%0A%5D%0A", semanticIris)
  semanticIris <- gsub("iris%3D%5B", "iris=%5B%0A", semanticIris)
  cat(semanticIris, file="./formdata2.txt")
  api.semanticSimilarity_query <- "curl -X POST -d @formdata2.txt 'http://kb.phenoscape.org/api/similarity/jaccard'"
  semanticSimilarityAPIResults <- system(api.semanticSimilarity_query, intern=TRUE)
  results <- fromJSON(semanticSimilarityAPIResults)
  scores <- lapply(results$results, function(x) x$score)
  scores <- sapply(scores, function(x) if(is.null(x)) NA else(x))
  result_terms <- do.call(rbind, lapply(results$results, function(x) do.call(cbind, lapply(x$terms, curlUnescape))))
  semanticSimilarityMatrix <- matrix(NA, nrow=length(phenotype_labels), ncol=length(phenotype_labels))
  diag(semanticSimilarityMatrix) <- 1
  rownames(semanticSimilarityMatrix) <- colnames(semanticSimilarityMatrix) <- irisPhenotypes
  
  for(i in 1:nrow(result_terms)){
    semanticSimilarityMatrix[result_terms[i,1], result_terms[i,2]] <- semanticSimilarityMatrix[result_terms[i,2], result_terms[i,1]] <- scores[i]
  }
  rownames(semanticSimilarityMatrix) <- colnames(semanticSimilarityMatrix) <- phenotype_labels
  return(semanticSimilarityMatrix)
  #write.csv(semanticSimilarityMatrix, "semanticSimilarityMatrix.csv")
  
}