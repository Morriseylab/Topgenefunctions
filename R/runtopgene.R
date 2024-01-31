#'Run topGene API and convert it to enrichResult
#' @param genelist vector of genelist
#' @param categories GO categories to filter results by. Defaults to Molecular Function, Biological Process and
#' @import httr jsonlite dplyr tidyr multienrichjam
#' @export

#g=c("FLDB","APOE","ENSG00000113196","ENSMUSG00000020287")

runtopgene <- function(genelist,categories=c("GeneOntologyMolecularFunction","GeneOntologyBiologicalProcess","GeneOntologyCellularComponent")){
  pc_json <- list(symbols=genelist)
  res <- POST("https://toppgene.cchmc.org/API/lookup"
              , body = pc_json
              , encode = "json")
  appData <- content(res)
  genes=appData$Genes
  genes=as.data.frame(do.call(rbind, genes))
  entrez=genes$Entrez
  pc_json2 <- list(Genes = entrez)
  res2 <- POST("https://toppgene.cchmc.org/API/enrich"
               , body = pc_json2
               , encode = "json")
  appData2 <- content(res2)
  df=appData2$Annotations
  df2=as.data.frame(do.call(rbind, df)) %>% tidyr::unnest(Genes) %>% tidyr::unnest_wider(Genes) %>%
    group_by(ID) %>% mutate(geneID = paste(Entrez, collapse="/")) %>% select(-Entrez) %>% group_by(ID) %>%
    mutate(Symbol = paste(Symbol, collapse=",")) %>% distinct()

  df2$GeneRatio=paste0(df2$GenesInTermInQuery,"/",df2$GenesInTerm,sep="")
  df2$BgRatio=paste0(df2$GenesInQuery,"/",df2$TotalGenes,sep="")
  df2 = df2 %>% rename('Description'='Name','pvalue'='PValue','qvalue'='QValueBonferroni','Count'='GenesInTermInQuery')
  df3=data.frame(lapply(df2, function(x) unlist(x)))
  rownames(df3)=df3$ID
  df3 = df3[df3$Category %in% categories,]
  edo=enrichDF2enrichResult(df3,pAdjustMethod = "BH",keyColname = "ID", geneColname = "geneID", pvalueColname = "QValueFDRBH", descriptionColname = "Description", pvalueCutoff = 0.05)
  return(edo)
}
