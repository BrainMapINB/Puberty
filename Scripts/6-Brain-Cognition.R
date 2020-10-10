# First of all, clear workspace
rm(list=ls())
# Save graphical parameters
op <- par()

# Read phenotypic data
info <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Phenotypic/phenotypic.csv")
# Remove those sessions with excesive motion artifact
info <- info[which(as.logical(info$QC)),]

# Extract time-series files (those with GSR)
ts_files <- paste0("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Brain/GSR/sub-",
                   info$BIDS.ID,"/ses-",info$Session,"/pp150v_wGSR_NIHPD_MNI_P264_ts.txt")

# Remove those without neuropsychological assessment
ts_files <- ts_files[which(!is.na(info$CST_Score))]
info <- info[which(!is.na(info$CST_Score)),]

# Compute connectivity matrices
ConnMx <- sapply(ts_files,
                 function(x){
                   ts <- read.table(x)
                   cmx <- cor(ts)
                   if(sum(is.na(cmx))>0) cmx[is.na(cmx)] <- 0
                   return(cmx)
                 }
)
# Reshape object
nroi <- sqrt(dim(ConnMx)[1])
ConnMx_dim <- c(nroi, nroi, nrow(info))
ConnMx <- array(ConnMx, dim = ConnMx_dim)

# Load cost-threshold functions
source("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Scripts/Auxiliary/cost_th_ConnMx.R")
# Load 'psych' package for Fisher r-to-z transformation
if(!require(psych)) install.packages("psych"); library(psych)

# Apply thresholds to Connectivity Matrix
ConnMx_th <- apply(ConnMx,3,function(x) cost_th_ConnMx(x,0.25))
ConnMx_th <- array(ConnMx_th, dim = ConnMx_dim)
# Substract diagonal
for(ii in 1:ConnMx_dim[3]) diag(ConnMx_th[,,ii]) <- 0
# Fisher's z transformation
ConnMx_th <- fisherz(ConnMx_th)

# Extract Functional Network inferences
atlas <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Atlas/P264/power264.csv")

# Functional network connectome by weighted-degree
# Set Functional Networks combinations (connectome upper triangle)
fn_comb <- which(upper.tri(matrix(nrow = nlevels(atlas$Net), ncol = nlevels(atlas$Net)), diag = T), arr.ind = T)
# Remove 'Uncertain' combinations
fn_comb <- fn_comb[sapply(1:nrow(fn_comb), function(x) return(all(fn_comb[x,]!=12))),]
# Compute w-degree for each combination
fn_deg <- sapply(1:nrow(fn_comb), function(x){
  # Find subnetwork indices
  subidx <- unique(c(which(atlas$Net==levels(atlas$Net)[fn_comb[x,1]]),
              which(atlas$Net==levels(atlas$Net)[fn_comb[x,2]])))
  # Create subnetwork
  submx <- ConnMx_th[subidx,subidx,]
  # Compute weighted-degree
  return(rowMeans(apply(submx,3:2,sum)))
})
# Set labels
colnames(fn_deg) <- paste0(levels(atlas$Net)[fn_comb[,1]],".",
                           levels(atlas$Net)[fn_comb[,2]])

###########################################################################
# Generalized Additive Models a la Network-Based Statistics
###########################################################################

# Standardize cognitive scores
cog_labs <- c("CST_Score","CST_Persev","CST_RecPersev")
cog_idx <- match(cog_labs,names(info))
info[,cog_idx] <- scale(info[,cog_idx])

# Load 'gamm4' & 'ggplot2' packages
if(!require(gamm4)) install.packages("gamm4"); library(gamm4)
if(!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)

#Set random number generator seed
set.seed(18900217)
nPerm <- 5000
# Update the combinations to 13x13
fn_comb <- which(upper.tri(matrix(nrow = nlevels(atlas$Net)-1, ncol = nlevels(atlas$Net)-1), diag = T), arr.ind = T)

# If you prefer to generate those values follow the steps within this 'if' command
# (but be aware that some computations may be take a while)
if(0){
  # Apply GAM-NBS
  for(ii in 1:length(cog_labs)){
    
    # Create empty object to store results
    res <- matrix(as.numeric(NA), nrow = ncol(fn_deg), ncol = 4)
    colnames(res) <- c("edf","Ref.df","F","p-value")
    rownames(res) <- colnames(fn_deg)
    
    #Compute GAM edgewise
    for(kk in 1:nrow(res)){
      # Set edge
      info$y <- fn_deg[,kk]
      # Set model
      grp_form <- as.formula(paste0(cog_labs[ii]," ~ s(y, k=4) + PDSe + Sex + FDRMS + Coil"))
      fit <- tryCatch(gamm4::gamm4(grp_form, data = info, REML = F),
                      error=function(e) NA)
      # Store results
      if(!anyNA(fit)) res[kk,] <- summary(fit$gam)$s.table
    }
    
    # Apply NBS
    # Initial variables for component search
    obs_comp <- 1:(nlevels(atlas$Net)-1)
    obs_list <- vector("list",1)
    # Find F-statistics (p < 0.05) edges
    edges <- which(res[,4]<0.01)
    if(length(edges)>0){
      # Generate empty object for component label
      component <- vector("integer", length(edges))
      # Store strength
      strengh <- res[edges,3]
      for(jj in 1:length(edges)){
        # Maximum label to be changed
        max_lab <- max(obs_comp[fn_comb[edges[jj],]])
        min_lab <- min(obs_comp[fn_comb[edges[jj],]])
        if(sum(component==max_lab)>0) component[which(component==max_lab)] <- min_lab
        pre_idx <- which(obs_comp==max_lab)
        obs_comp[pre_idx] <- component[jj] <- min_lab
      }
      # Store observed components
      obs_list <- cbind(component,strengh)
      
      # Create auxiliary objects
      datPerm <- info
      permF <- matrix(0, nrow = nPerm, ncol = 2)
      colnames(permF) <- c("max_comp","max_sum")
      # Compute permutations
      for(pp in 1:nPerm){
        # Print progress
        if(pp%%100 == 0) print(pp)      
        # Permutate dependent variable
        datPerm[,cog_idx[ii]] <- datPerm[sample(1:nrow(datPerm)),cog_idx[ii]]
        
        # Create empty object to store results
        res <- matrix(as.numeric(NA), nrow = ncol(fn_deg), ncol = 4)
        
        #Compute GAM edgewise
        for(kk in 1:nrow(res)){
          # Set edge
          datPerm$y <- fn_deg[,kk]
          # Set model
          grp_form <- as.formula(paste0(cog_labs[ii]," ~ s(y, k=4) + PDSe + Sex + FDRMS + Coil"))
          fit <- tryCatch(gamm4::gamm4(grp_form, data = datPerm, REML = F),
                          error=function(e) NA)
          # Store results
          if(!anyNA(fit)) res[kk,] <- summary(fit$gam)$s.table
        }
        
        # Initial variables for component search
        perm_comp <- 1:(nlevels(atlas$Net)-1)
        perm_list <- vector("list",1)
        # Find F-statistics (p < 0.01) edges
        edges <- which(res[,4]<0.01)
        if(length(edges)>0){
          # Generate empty object for component label
          component <- vector("integer", length(edges))
          # Store strength
          strengh <- res[edges,3]
          for(jj in 1:length(edges)){
            # Maximum label to be changed
            max_lab <- max(perm_comp[fn_comb[edges[jj],]])
            min_lab <- min(perm_comp[fn_comb[edges[jj],]])
            if(sum(component==max_lab)>0) component[which(component==max_lab)] <- min_lab
            pre_idx <- which(perm_comp==max_lab)
            perm_comp[pre_idx] <- component[jj] <- min_lab
          }
          perm_list[[1]] <- cbind(component,strengh)
          permF[pp,1] <- max(table(perm_list[[1]][,1]))
          permF[pp,2] <- max(aggregate(perm_list[[1]][,2]~perm_list[[1]][,1], perm_list[[1]], sum)[,2])
        } 
      }# pp in 1:nPerm
      
      # Quantile observed values within the null distribution
      fwe_list <- vector("list",1)
      # Number of components FWE
      obs_comp <- table(obs_list[,1])
      for(cc in 1:length(obs_comp)) obs_comp[cc] <- sum(permF[,1] > obs_comp[cc])/nPerm
      # Components strength
      obs_strn <- aggregate(obs_list[,2]~obs_list[,1], obs_list, sum)
      for(cc in 1:nrow(obs_strn)) obs_strn[cc,2] <- sum(permF[,2] > obs_strn[cc,2])/nPerm
      # Save results
      obs_strn <- cbind(obs_strn, c(obs_comp))
      names(obs_strn) <- c("Component","weightFWE","ncompFWE")
      fwe_list[[1]] <- obs_strn
      # Save FWE results
      outfile <- paste0("pFWE_",cog_labs[ii],".csv")
      #write.csv(fwe_list[[1]], outfile, quote = F, row.names = F)
      
      # Save NBS results (if pFWE < 0.05)
      if(!all(fwe_list[[1]][,2]>0.05)){
        # Create output directory
        fwe_dir <- file.path("Cognition","NBS",cog_labs[ii])
        if(!dir.exists(fwe_dir)) dir.create(fwe_dir, recursive = T)
        # For each FWE component
        fwe_comp <- fwe_list[[1]][which(fwe_list[[1]][,2]<0.05),1]
        for(ff in fwe_comp){
          sig_edge <- which(obs_list[,1]==ff)
          # Save component edges
          outfile <- paste0(fwe_dir,"/",cog_labs[ii],"_comp",ff,"_edges.csv")
          write.csv(obs_list[sig_edge,], outfile, quote = F)
          # Save null distribution
          outfile <- paste0(fwe_dir,"/",cog_labs[ii],"_comp",ff,"_nudist.csv")
          write.csv(permF, outfile, quote = F, row.names = F)
          # Print individual edges
          plt_dir <- paste0(fwe_dir,"/comp",ff,"_edges")
          if(!dir.exists(plt_dir)) dir.create(plt_dir)
          for(ee in sig_edge){
            # Compute brain residuals
            edge_idx <- match(rownames(obs_list)[ee],colnames(fn_deg))
            info$y <- fn_deg[,edge_idx]
            info$res <- resid(lm(y~PDSe+Sex+FDRMS+Coil, data = info))
            # Compute GAM effect
            plt_form <- as.formula(paste0(cog_labs[ii]," ~ + s(res, k=4)"))
            info$preds <- predict(gamm4::gamm4(plt_form, data = info, REML = F)$gam)
            info$preds_se <- predict(gamm4::gamm4(plt_form, data = info, REML = F)$gam,
                                     se.fit = T)[[2]]
            gSC <- ggplot(data = info, mapping = aes(x = res)) + 
              geom_ribbon(aes(ymax = preds + 1.96 * preds_se,
                              ymin = preds - 1.96 * preds_se,
                              linetype = NA), alpha = 0.2) + 
              geom_line(aes(y = preds), size = 1.5) +
              ylab(paste0(cog_labs[ii]," (z)")) +
              xlab(paste0("FC Degree (res)")) +
              ggtitle(colnames(fn_deg)[edge_idx]) +
              theme_bw()
            # Save plot
            outfile <- paste0(plt_dir,"/",gsub("[.]","",colnames(fn_deg)[edge_idx]),".pdf")
            ggsave(outfile, gSC, "pdf")
          }#for(ee in sig_edge)
        }#for(ff in fwe_comp)
      }#if(!all(fwe_list[[1]][,2]>0.05))
    }#if(length(edges)>0)
  }#for(ii in 1:length(cog_labs))
  
}

###########################################################################
# Generate Chord diagram
###########################################################################

# Load 'circlize' package
if(!require(circlize)) install.packages("circlize"); library(circlize)
# Load 'scales' package
if(!require(scales)) install.packages("scales"); library(scales)

# List significant components (pFWE < 0.05)
nbs_files <- c("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Cognition/NBS/CST_Score/CST_Score_comp1_edges.csv",
               "https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Cognition/NBS/CST_RecPersev/CST_RecPersev_comp1_edges.csv")
# Analyse each significant NBS-GAM results
par(mar = c(1, 1, 1, 1), mfrow = c(1,2))
for(nn in 1:length(nbs_files)){
  # Read NBS-GAM result
  nbs_res <- read.csv(nbs_files[nn])
  # Extract labels
  cog_nbs <- which(sapply(cog_labs, function(x) grepl(x, nbs_files[nn])))
  # Extract edge labels
  edge_labs <- as.character(nbs_res[,1])
  # Empty objecto to store
  lin_mat <- matrix(0, nrow = length(edge_labs), ncol = 3)
  for(ee in 1:length(edge_labs)){
    # Compute brain residuals
    info$y <- fn_deg[,match(edge_labs[ee],colnames(fn_deg))]
    info$res <- resid(lm(y~PDSe+Sex+FDRMS+Coil, info))
    # Estimate fitted curve
    plt_form <- as.formula(paste0(cog_labs[cog_nbs],"~ + s(res, k=4)"))
    info$preds <- predict(gamm4::gamm4(plt_form, data = info, REML = F)$gam)
    
    # Store slope and point of inflection
    lo <- loess(preds~res, data = info)
    x <- seq(min(info$res),max(info$res),0.01)
    y <- predict(lo,x)
    # Derivative
    infl <- c(FALSE, diff(diff(y)>0)!=0)
    # Number of inflection points
    lin_mat[ee,1] <- sum(infl)
    if(lin_mat[ee,1]>0){
      # Last inflection point
      lin_mat[ee,2] <- x[which(infl==T)][lin_mat[ee,1]]
      # Compute slope after point of inflection
      y <- y[x>=lin_mat[ee,2]]
      x <- x[x>=lin_mat[ee,2]]
    }
    # Compute slope
    lin_mat[ee,3] <- coefficients(lm(y~x))[2]
  }
  rownames(lin_mat) <- edge_labs
  colnames(lin_mat) <- c("IPnum","IPx","slope")
  
  # Draw circleplot
  # Abbreviation
  net_lvl <- levels(atlas$Net)[-12]
  net_n <- length(net_lvl)
  # Find significant links
  sig_tri <- fn_comb[match(edge_labs,colnames(fn_deg)),]
  # Create network factor
  net_fac <- as.factor(1:net_n)
  levels(net_fac) <- net_lvl
  # Initialize the plot
  outfile <- file.path("Cognition", gsub("edges.csv","NBS-GAM.pdf",basename(nbs_files[nn])))
  #pdf(file = outfile, width = 4, height = 4)
  circos.clear()
  circos.initialize(factors = net_fac, xlim=c(-1,1))
  # Build the regions of track #1
  circos.trackPlotRegion(factors = net_fac, ylim = c(-1,1),
                         bg.col = "aliceblue",
                         bg.border = "black")
  # Add labels
  circos.text(2.24*(1:net_n),rep(0,net_n),
              labels = net_lvl,
              facing = "bending.inside", cex=1.2)
  # Add a links between a point and another
  # Draw links
  for(ii in 1:nrow(sig_tri)){
    # Set random position
    p1 <- runif(1,-1,1)+c(-0.1,0.1); p2 <- runif(1,-1,1)+c(-0.1,0.1)
    #p1 <- c(-0.1,0.1); p2 <- c(-0.1,0.1)
    # Set color based on positive/negative spline slope
    ifelse(lin_mat[ii,3]>0, tp_col <- "red", tp_col <- "blue")
    # Draw link
    circos.link(net_fac[sig_tri[ii,1]], p1,
                net_fac[sig_tri[ii,2]], p2,
                col = scales::alpha(tp_col,.8),
                border = scales::alpha(tp_col,.5),
                h.ratio=0.8)
  }
  #dev.off()
}#for(nn in 1:length(nbs_files))
