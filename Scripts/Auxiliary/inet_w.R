# Creator: zchuri
# Date: 30 March 2015
# Function: iEff_w

# IEFF_W(CONNMX,TIPO)
# Computes given 'tipo' global or local efficiency given a connectivity
# matrix 'ConnMx'.
#
#   The global efficiency is the average of inverse shortest path length,
#   and is inversely related to the characteristic path length.
#   The local efficiency is the global efficiency computed on the
#   neighborhood of the node, and is related to the clustering coefficient.
#
# -------------------------------------------------------------------------
# Inputs:
# CONNMX - simmetric connectivity matrix.
# TIPO - string char, in order to compute 'global' or 'local' efficiency.
# -------------------------------------------------------------------------
# Output:
# E - vector containing efficiency measures for each node.
# -------------------------------------------------------------------------

iEff_w <- function(ConnMx,tipo="global") {
  
  # First open "igraph" package
  if(require(igraph)==0){
    install.packages("igraph")
    library(igraph)
  }
  
  ConnMx_inv <- ConnMx
  ConnMx_inv[(ConnMx!=0)] <- 1/ConnMx[(ConnMx!=0)]
  n_nodes <- dim(ConnMx_inv)[1]
  gc() # Free RAM memory
  
  # Convert connectivity matrix into a igraph
  # print("Convert ConnMx_inv into igraph...")
  # antes <- Sys.time() # Before
  ig <- graph.adjacency(ConnMx_inv,mode="undirected",weighted=T)
  # ahorita <- Sys.time() - antes # After
  # print(ahorita)
  rm(ConnMx_inv)
  gc() # Free RAM memory
  
  if(tipo=="global"){
    # Create shortest paths matrix
    ig_dist <- distances(ig,weights=E(ig)$weight,algorithm="dijkstra")
    rm(ig)
    gc() # Free RAM memory
    
    # Compute the inverse of the shortest paths
    ig_dist <- 1/ig_dist
    diag(ig_dist) <- 0
    gc() # Free RAM memory
    
    # Computes global efficiency
    E <- apply(ig_dist,1,sum)/(n_nodes-1)
  } else if(tipo=="local"){
    # Modified from code by Christopher G. Watson; cgwatson@bu.edu (http://pastebin.com/XqkEYtJS)
    
    # First get the degree (k) of each node
    k <- degree(ig)
    E <- vector("numeric",n_nodes)
    nodes <- which(k > 1)
    n_nodes <- length(nodes)
    #print(paste0("Number of nodes: ",n_nodes))
    # Computes local efficiency
    if(require(doMC)==0){
      install.packages("doMC")
      library(doMC)
    }
    registerDoMC()
    E[nodes] <- simplify2array(mclapply(nodes, function(x) {
      # Print percent
      # if(x%%trunc(n_nodes/100)==0) print(paste0(x,"/",n_nodes))
      # Get neighbours
      neighbs <- neighbors(ig, v=x)
      ig_sub <- induced.subgraph(ig, neighbs)
      ig_dist <- distances(ig_sub,weights=E(ig_sub)$weight,algorithm="dijkstra")
      # Compute the inverse of the shortest paths
      ig_dist <- 1/ig_dist
      diag(ig_dist) <- 0
      # Computes global efficiency
      sum(ig_dist)/(k[x]^2-k[x])
    }, mc.cores=getDoParWorkers()) #end mclapply
    ) #end simplify2array
  } else stop("Mi no reconocer ese 'tipo'... wey")
  
  return(E)
  
}#function iEFF_w


# Creator: Zeus Gracia
# Date: 11 August 2016
# Function: iCC_w

# ICC_W(CONNMX)
# Computes weighted clustering coefficient (per node) given a connectivity
# matrix 'ConnMx'.
#
#   The cluster coefficient is the geometric mean of triangles divided by
#   degree plus degree minus one. Here we used Barrat et al. (2004)
#   formula.
#
# -------------------------------------------------------------------------
# Inputs:
# CONNMX - simmetric (undirected) connectivity matrix.
# -------------------------------------------------------------------------
# Output:
# IG_CC - numeric vector containing clustering coefficient measure for
# each node.
# -------------------------------------------------------------------------

iCC_w <- function(ConnMx) {
  
  # First open "igraph" package
  if(require(igraph)==0){
    install.packages("igraph")
    library(igraph)
  }
  
  # Convert connectivity matrix into a igraph
  # antes <- Sys.time() # Before
  ig <- graph.adjacency(ConnMx,mode="undirected",weighted=T)
  # ahorita <- Sys.time() - antes # After
  # print(ahorita)
  ig_cc <- transitivity(ig,type="barrat",weights=E(ig)$weight)
  rm(ig); gc() # Free RAM memory
  if(sum(is.na(ig_cc))>0)  ig_cc[is.na(ig_cc)] <- 0
  
  return(ig_cc)
  
}#function iCC_w

# Creator: Zeus Gracia
# Date: June 20th 2017
# Function: iL_w

# IL_W(CONNMX)
# Computes weighted shortest path (per node) given a connectivity
# matrix 'ConnMx'.
#
#   L(i) is the average distance between node i and all other nodes.
#   Charasteristic path lenth is the average of all L(i).
#
# -------------------------------------------------------------------------
# Inputs:
# CONNMX - simmetric (undirected) connectivity matrix.
# -------------------------------------------------------------------------
# Output:
# E - numeric vector containing average shortest paths for each node.
# -------------------------------------------------------------------------

iL_w <- function(ConnMx) {
  
  # First open "igraph" package
  if(require(igraph)==0){
    install.packages("igraph")
    library(igraph)
  }
  
  ConnMx_inv <- ConnMx
  ConnMx_inv[(ConnMx!=0)] <- 1/ConnMx[(ConnMx!=0)]
  n_nodes <- dim(ConnMx_inv)[1]
  gc() # Free RAM memory
  
  # Convert connectivity matrix into a igraph
  # print("Convert ConnMx_inv into igraph...")
  # antes <- Sys.time() # Before
  ig <- graph.adjacency(ConnMx_inv,mode="undirected",weighted=T)
  # ahorita <- Sys.time() - antes # After
  # print(ahorita)
  rm(ConnMx_inv)
  gc() # Free RAM memory
  
  # Create shortest paths matrix
  ig_dist <- distances(ig,weights=E(ig)$weight,algorithm="dijkstra")
  rm(ig)
  gc() # Free RAM memory
  # Set the path length of disconnected nodes to the maximum observed of the network
  # in order to avoid infinite numbers (Fornito et al., 2010)
  if(sum(is.infinite(ig_dist))>0) ig_dist[is.infinite(ig_dist)] <- max(ig_dist[!is.infinite(ig_dist)])
    
  # Compute the inverse of the shortest paths
  diag(ig_dist) <- 0
  
  # Computes global efficiency
  E <- apply(ig_dist,1,sum)/(n_nodes-1)
  return(E)
}#function iL_w