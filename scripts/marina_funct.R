library(magrittr) # pipe operator %>%
library(igraph)
library(KEGGgraph)
library(dplyr)
library(biomaRt)
library(gdata)
library(xlsx) 
library(here)
library(drugbankR)
library(stringr)
library("org.Hs.eg.db")
library(Matching)

## KEGG to Data Base 
  # Parse KEGG metadata
  #TODO: retrieve the kegg vertex name
  #TODO: retrieve each gene/protein for each vertex, save it in a new field
  # Input: Kegg_pathway_id ex.  path_id <- "03460"
  # Output: list with all the data retrieved
  kegg2db <- function(kegg_pathway_id) {
  gnell <- tempfile() %T>%
  {
    retrieveKGML(
      kegg_pathway_id,
      organism = "hsa",
      destfile = .,
      method = "curl",
      quiet = TRUE
    )
  } %>%
  {
    
    
    parseKGML2Graph(., genesOnly = TRUE)
  }
  
  gene_db <- getKEGGnodeData(gnell)
  
  return (list("gene"=gene_db, "graph"=gnell))
}

## KEGG gene data base (gene_db) to data frame. 
  # Converts gene database obtain from kegg2db into data frame manageable
  # Input: is the outuput of kegg2db function ex. db <- kegg2db(path_id). "db"
  # Output: Data list of all info retrieve from KEGG pathway given 
  kegg2df <- function(db) {
  
  gene_db <- db$gene
  gnell <- db$graph
  
  vertex_df <- gene_db %>%
  {
    data.frame(
      kegg_id = unlist(lapply(., function(item)
        item@name[1])),
      label = unlist(lapply(., function(item)
        strsplit(item@graphics@name, ",")[[1]][1])),
      x = unlist(lapply(., function(item)
        item@graphics@x[1])),
      y = unlist(lapply(., function(item)
        -1 * item@graphics@y[1])),
      inside = unlist(lapply(., function(item)
        str_extract_all(item@link,"\\(?[0-9]+\\)?")[[1]]%>% paste("hsa:", .) %>% paste(., collapse = "; "))),
      entrzIDs = unlist(lapply(., function(item)
        str_extract_all(item@link,"\\(?[0-9]+\\)?")[[1]]%>% paste(., collapse = "; "))),
      stringsAsFactors = F
    )
  }
  
  g_init <- igraph.from.graphNEL(gnell)
  V(g_init)$name <- vertex_df$kegg_id
  vertex_df <- filter(vertex_df, !duplicated(kegg_id))
  
  edge_list <- getKEGGedgeData(gnell) %>%
    lapply(function(item) {
      if (length(item@subtype) > 0) {
        subtype_info <- item@subtype
        #TODO: use the KEGG edge hierarchy
        if (length(subtype_info) > 1) {
          return(subtype_info[[2]]@name)
        } else {
          return(subtype_info$subtype@name)
        }
      }
      NA
    }) %>%
    unlist %>%
    {
      cbind(get.edgelist(g_init), type = .)
    } %>%
    data.frame
  
  edge_list <- edge_list %>%
    as.data.frame %>%
    unique
  
  vertex_df <- vertex_df %>%
    unique %>%
    {
      filter(., !duplicated(kegg_id))
    }
  
  return(list("vertex"=vertex_df, "edges"=edge_list))
}

## Convert data frames into igraph graph
  # Input: edges_df, vertes_df ex. vertex_final, edges_final
  # Output: igraph object to read via tkplot(object)
  convert2igraph<- function(vertex_df, edge_list){
  g_final <- graph.data.frame(edge_list,
                              directed = TRUE,
                              vertices = vertex_df)
  
  V(g_final)$kegg_id <- V(g_final)$name
  V(g_final)$name <- V(g_final)$kegg_id
  V(g_final)$color <- ifelse(vertex_final$label %in% vertex_extra$label, "orange","lightgreen")
  V(g_final)$label.color <- "black"
  V(g_final)$frame.color <-ifelse(vertex_final$label %in% vertex_extra$label, "orange","Darkgreen")
  V(g_final)$label.cex <-1.5
  V(g_final)$label.family <-"Helvetica"
  V(g_final)$size <- 30
  E(g_final)$color<- ifelse(edges_final$V1 %in% edges_extra$V1, "orange","black")
  E(g_final)$width <- 3
  g_final <- g_final - V(g_final)[degree(g_final) == 0]

  
  return(g_final)
}

## Add vertex and edges if missing (therefor must be two ".xls" documents to be readed with the info for the extra vertex and edges)
  # Input: data from kegg2df  function ex. lst_df <- kegg2df(db)
  # Output: data frames with the added info (vertex_final, edges_final) and separately the df with the extra info to add (vertex_extra, edges_extra)
  adding_ve <- function (lista_datos){
  
  vertex_df <- lista_datos$vertex
  edges_df <- lista_datos$edges
  
    if(!file.exists(vertexfile)|!file.exists(edgesfile)){
      return(list("edges_final"=edges_df, "vertex_final"=vertex_df,"edges_extra"=edges_extra, "vertex_extra"=vertex_extra))
    }else{
      vertex_extra <- here("data", paste0(path_id, "_vertex.xls")) %>% read.xls(.)
      edges_extra <- here("data", paste0(path_id, "_edges.xls")) %>% read.xls(.) 
      rownames(edges_extra)<- edges_extra$X 
      edges_extra <- edges_extra[,-1]
      edges_final <- rbind(edges_df, edges_extra)
      vertex_final <- rbind(vertex_df, vertex_extra)                 
      return(list("edges_final"=edges_final,"vertex_final"=vertex_final,"edges_extra"=edges_extra, "vertex_extra"=vertex_extra))
    }
  }
    
    
 ## Checks if all the genes are in both edges and vertex df
  # Input: final vertex and edges df
  # Output: message correct or not correct
    check_genes <- function(edges, vertex){
      
      input_vertices <- edges$V1
      kegg_ids <- vertex$kegg_id
      
      checking <-input_vertices %in% kegg_ids
      
      if (all(checking == TRUE)){
        print("All vertex and edges are correct")
      }else{
        print("Some vertex names in edge list are not listed in vertex data frame or viceversa")
        print("error in vertex:") 
        print(which(!(edges_final$V1 %in% vertex_df$kegg_id)))
        print("error in edges:") 
        print(which(!(vertex_df$kegg_id %in%edges_final$V1 )))  
      }
      
    }
    
## Conversor function i.e. entrezids into UniProtKB ids
  # Input: entrez_ids i.e 2175 (vertex_final$inside without "hsa: ")
  # Output: UniProtKB ids
    conversor <- function (ids){  
      
      unlist(lapply(ids,FUN= AnnotationDbi::get, org.Hs.eg.db::org.Hs.egUNIPROT))
    }
    
## Find matched UniProtKB codes in Drugbank_dataframe$targets
  # Input: UniProtKB ids in vertex_final df and Drugbankdb targets (both are in UniProtKB ids)
  # Output: List of TRUE/FALSE /NA
    checking_targets <- function(ids, targets){
      #input uniprot-ids (take from myuni) and the drugbank data base fields targets
      hit<-ids %in% targets
      if (any(hit)) {
        res <-ids[hit]
      }else{
        res <- NA
      }
      
      return(res)
    }
    

## Find matched UniProtKB codes in Drugbank_dataframe$targets
  # Input: UniProtKB ids in vertex_final df and Drugbankdb targets (both are in UniProtKB ids)
    # Output: List of TRUE/FALSE /NA
    
  checking_targets <- function(ids, targets){
  #input uniprot-ids (take from myuni) and the drugbank data base fields targets
  hit<-ids %in% targets
  if (any(hit)) {
    res <-ids[hit]
  }else{
    res <- NA
  }
  
  return(res)
}

