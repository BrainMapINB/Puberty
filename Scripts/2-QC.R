# First of all, clear workspace
rm(list=ls())
# Save graphical parameters
op <- par()

# Check if QC files exist
qc_files <- c("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/QC/edge_motion.csv",
              "https://raw.githubusercontent.com/BrainMapINB/Puberty/main/QC/edge_motion_euclidean_dist.csv",
              "https://raw.githubusercontent.com/BrainMapINB/Puberty/main/QC/graph_motion.csv")

# If you prefer to generate those values follow the steps within this 'if' command
# (but be aware that some computations may be take a while)
if(0){
  # Read phenotypic data
  info <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Phenotypic/phenotypic.csv")
  # Discard those with excesive head-motion artifact
  info <- info[which(as.logical(info$QC)),]
  # Relevel the factor ID
  info$BIDS.ID <- factor(info$BIDS.ID)
  
  # Read atlas info
  atlas <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Atlas/P264/power264.csv")
  
  # Matrix interactions (upper triangle)
  tri_pos <- which(upper.tri(matrix(nrow = nrow(atlas), ncol = nrow(atlas))), arr.ind = T)
  # Compute Euclidean distance between P264 ROIs
  tri_dist <- sapply(1:nrow(tri_pos), function(x) sqrt(sum((atlas[tri_pos[x,1],2:4]-atlas[tri_pos[x,2],2:4])^2)) )
  
  # noGSR
  # Create time-series file list
  ts_files <- paste0("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Brain/noGSR/sub-",
         info$BIDS.ID,"/ses-",info$Session,"/pp150v_NIHPD_MNI_P264_ts.txt")
  # Compute connectivity matrices
  ConnMx <- sapply(ts_files,
                   function(x){
                     ts <- read.table(x)
                     cmx <- cor(ts)
                     if(sum(is.na(cmx))>0) cmx[is.na(cmx)] <- 0
                     cmx
                   }
  )
  # Reshape object
  ConnMx_dim <- c(nrow(atlas), nrow(atlas), nrow(info))
  ConnMx <- array(ConnMx, dim = ConnMx_dim)
  
  # Substract diagonal
  for(ii in 1:ConnMx_dim[3]) diag(ConnMx[,,ii]) <- 0
  
  # Apply cost threshold
  source("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Scripts/Auxiliary/cost_th_ConnMx.R")
  # Apply thresholds to Connectivity Matrix
  ConnMx_th <- apply(ConnMx,3,function(x) cost_th_ConnMx(x,0.25))
  ConnMx_th <- array(ConnMx_th,dim=ConnMx_dim)
  
  # Load 'psych' package
  if(!require(psych)) install.packages("psych"); library(psych)
  
  # Fisher's z transformation
  ConnMx <- fisherz(ConnMx)
  ConnMx_th <- fisherz(ConnMx_th)
  
  # Compute graph theory measures
  source("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Scripts/Auxiliary/inet_w.R")
  ConnMx_str_woGSR <- apply(ConnMx,3:2,sum)
  ConnMx_deg_woGSR <- apply(ConnMx_th,3:2,sum)
  ConnMx_clus_woGSR <- t(sapply(1:ConnMx_dim[3], function(x) iCC_w(ConnMx_th[,,x])))
  ConnMx_eff_woGSR <- t(sapply(1:ConnMx_dim[3], function(x) iEff_w(ConnMx_th[,,x])))
  ConnMx_path_woGSR <- t(sapply(1:ConnMx_dim[3], function(x) iL_w(ConnMx_th[,,x])))
  
  # Load 'lme4' package
  if(!require(lme4)) install.packages("lme4"); library(lme4)
  
  # Compute head motion relationships (strength and degree)
  woGSR_str_lmeT <- sapply(1:ncol(ConnMx_str_woGSR), function(x){
    info$fc <- ConnMx_str_woGSR[,x]
    summary(lmer(fc~FDRMS + (1|BIDS.ID), data = info))$coefficients[2,3]
  })
  woGSR_deg_lmeT <- sapply(1:ncol(ConnMx_deg_woGSR), function(x){
    info$fc <- ConnMx_deg_woGSR[,x]
    summary(lmer(fc~FDRMS + (1|BIDS.ID), data = info))$coefficients[2,3]
  })
  woGSR_clus_lmeT <- sapply(1:ncol(ConnMx_clus_woGSR), function(x){
    info$fc <- ConnMx_clus_woGSR[,x]
    summary(lmer(fc~FDRMS + (1|BIDS.ID), data = info))$coefficients[2,3]
  })
  woGSR_eff_lmeT <- sapply(1:ncol(ConnMx_eff_woGSR), function(x){
    info$fc <- ConnMx_eff_woGSR[,x]
    summary(lmer(fc~FDRMS + (1|BIDS.ID), data = info))$coefficients[2,3]
  })
  woGSR_path_lmeT <- sapply(1:ncol(ConnMx_path_woGSR), function(x){
    info$fc <- ConnMx_path_woGSR[,x]
    summary(lmer(fc~FDRMS + (1|BIDS.ID), data = info))$coefficients[2,3]
  })
  
  # Compute head motion relationships (edgewise)
  tic <- Sys.time()
  woGSR_edge_lmeT <- sapply(1:nrow(tri_pos), function(x){
    info$fc <- ConnMx[tri_pos[x,1],tri_pos[x,2],]
    summary(lmer(fc~FDRMS + (1|BIDS.ID), data = info))$coefficients[2,3]
  })
  toc <- Sys.time()
  # Compare to euclidean distance
  #plot(tri_dist,woGSR_edge_lmeT); abline(lm(woGSR_edge_lmeT~tri_dist), col="red", lwd=3)
  
  #-------------------------------------------------------------------------
  # Same for Global Signal Regression data
  
  # Create time-series file list
  ts_files <- paste0("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Brain/GSR/sub-",
                     info$BIDS.ID,"/ses-",info$Session,"/pp150v_wGSR_NIHPD_MNI_P264_ts.txt")
  
  # Compute connectivity matrices
  ConnMx <- sapply(ts_files,
                   function(x){
                     ts <- read.table(x)
                     cmx <- cor(ts)
                     if(sum(is.na(cmx))>0) cmx[is.na(cmx)] <- 0
                     cmx
                   }
  )
  # Reshape object
  ConnMx <- array(ConnMx, dim = ConnMx_dim)
  
  # Substract diagonal
  for(ii in 1:ConnMx_dim[3]) diag(ConnMx[,,ii]) <- 0
  
  # Apply cost threshold
  # Apply thresholds to Connectivity Matrix
  ConnMx_th <- apply(ConnMx,3,function(x) cost_th_ConnMx(x,0.25))
  ConnMx_th <- array(ConnMx_th,dim=ConnMx_dim)
  
  # Fisher's z transformation
  ConnMx <- fisherz(ConnMx)
  ConnMx_th <- fisherz(ConnMx_th)
  
  # Compute graph theory measures
  ConnMx_str_wGSR <- apply(ConnMx,3:2,sum)
  ConnMx_deg_wGSR <- apply(ConnMx_th,3:2,sum)
  ConnMx_clus_wGSR <- t(sapply(1:ConnMx_dim[3], function(x) iCC_w(ConnMx_th[,,x])))
  ConnMx_eff_wGSR <- t(sapply(1:ConnMx_dim[3], function(x) iEff_w(ConnMx_th[,,x])))
  ConnMx_path_wGSR <- t(sapply(1:ConnMx_dim[3], function(x) iL_w(ConnMx_th[,,x])))
  
  # Compute head motion relationships (strength and degree)
  wGSR_str_lmeT <- sapply(1:ncol(ConnMx_str_wGSR), function(x){
    info$fc <- ConnMx_str_wGSR[,x]
    summary(lmer(fc~FDRMS + (1|BIDS.ID), data = info))$coefficients[2,3]
  })
  wGSR_deg_lmeT <- sapply(1:ncol(ConnMx_deg_wGSR), function(x){
    info$fc <- ConnMx_deg_wGSR[,x]
    summary(lmer(fc~FDRMS + (1|BIDS.ID), data = info))$coefficients[2,3]
  })
  wGSR_clus_lmeT <- sapply(1:ncol(ConnMx_clus_wGSR), function(x){
    info$fc <- ConnMx_clus_wGSR[,x]
    summary(lmer(fc~FDRMS + (1|BIDS.ID), data = info))$coefficients[2,3]
  })
  wGSR_eff_lmeT <- sapply(1:ncol(ConnMx_eff_wGSR), function(x){
    info$fc <- ConnMx_eff_wGSR[,x]
    summary(lmer(fc~FDRMS + (1|BIDS.ID), data = info))$coefficients[2,3]
  })
  wGSR_path_lmeT <- sapply(1:ncol(ConnMx_path_wGSR), function(x){
    info$fc <- ConnMx_path_wGSR[,x]
    summary(lmer(fc~FDRMS + (1|BIDS.ID), data = info))$coefficients[2,3]
  })
  
  # Compute head motion relationships (edgewise)
  tic <- Sys.time()
  wGSR_edge_lmeT <- sapply(1:nrow(tri_pos), function(x){
    info$fc <- ConnMx[tri_pos[x,1],tri_pos[x,2],]
    summary(lmer(fc~FDRMS + (1|BIDS.ID), data = info))$coefficients[2,3]
  })
  toc <- Sys.time()
  # Compare to euclidean distance
  #plot(tri_dist,wGSR_edge_lmeT); abline(lm(wGSR_edge_lmeT~tri_dist), col="red", lwd=3)
  
  #-------------------------------------------------------------------------
  # Set datasets
  
  # Edgewise-motion vs distance
  edge_motion_dist_data <- data.frame(edge_woGSR=woGSR_edge_lmeT,
                                      edge_wGSR=wGSR_edge_lmeT,
                                      dist=tri_dist)
  outname <- "edge_motion_euclidean_dist.csv"
  write.csv(x = edge_motion_dist_data, file = outname,
            quote = FALSE, row.names = FALSE)
  # Edgewise-motion density
  tri_n <- length(tri_dist)
  edge_motion_data <- data.frame(edge_hm=c(woGSR_edge_lmeT,wGSR_edge_lmeT),
                                 GSR=c(rep("NO",tri_n),rep("YES",tri_n)))
  outname <- "edge_motion.csv"
  write.csv(x = edge_motion_data, file = outname,
            quote = FALSE, row.names = FALSE)
  # Graph theory measures
  gt_data <- data.frame(str_hm=c(woGSR_str_lmeT,wGSR_str_lmeT),
                        deg_hm=c(woGSR_deg_lmeT,wGSR_deg_lmeT),
                        clus_hm=c(woGSR_clus_lmeT,wGSR_clus_lmeT),
                        eff_hm=c(woGSR_eff_lmeT,wGSR_eff_lmeT),
                        path_hm=c(woGSR_path_lmeT,wGSR_path_lmeT),
                        GSR=c(rep("NO",nrow(atlas)),rep("YES",nrow(atlas))))
  outname <- "graph_motion.csv"
  write.csv(x = gt_data, file = outname,
            quote = FALSE, row.names = FALSE)
}

#-------------------------------------------------------------------------
# Figures

# Load 'ggplot2', 'gridExtra', & 'hexbin' package
if(!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)
if(!require(gridExtra)) install.packages("gridExtra"); library(gridExtra)
if(!require(hexbin)) install.packages("hexbin"); library(hexbin)

# Hexbin chart (NO GSR)
edge_dist <- read.csv(qc_files[2])
edge_hm_b1 <- coef(lm(edge_woGSR~dist, data = edge_dist))[2]
g1 <- ggplot(edge_dist, aes(x=dist, y=edge_woGSR) ) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method='lm', formula= y~x, color = "red") +
  ggtitle(paste0("No GSR (slope=",round(edge_hm_b1, 4),")")) + 
  xlab("Euclidean Dist. (mm)") + 
  ylab("Edge-motion (LME-t)") +
  theme_bw() + ylim(-5, 5) + xlim(0,170)

# Hexbin chart
edge_hm_b1 <- coef(lm(edge_wGSR~dist, data = edge_dist))[2]
g2 <- ggplot(edge_dist, aes(x=dist, y=edge_wGSR) ) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method='lm', formula= y~x, color = "red") +
  ggtitle(paste0("with GSR (slope=",round(edge_hm_b1, 4),")")) + 
  xlab("Euclidean Dist. (mm)") + 
  ylab("Edge-motion (LME-t)") +
  theme_bw() + ylim(-5, 5) + xlim(0,170)

# Arrange plots
g3 <- grid.arrange(g1, g2, nrow = 1, ncol=2)
# Save
outname <- "edge_motion_euclidean_dist.pdf"
#ggsave(filename = outname, plot = g3, device = "pdf", width = 7, height = 2.5)

# Density plots
edge_data <- read.csv(qc_files[1])
d1 <- ggplot(data=edge_data, aes(x=edge_hm, group=GSR, fill=GSR)) +
  geom_density(adjust=1.5, alpha=.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("Edgewise") + xlab("Motion effect (LME-t)") +
  theme_bw()

# Auxiliary dataset
gt_data <- read.csv(qc_files[3])
# Strength density
d2 <- ggplot(data=gt_data, aes(x=str_hm, group=GSR, fill=GSR)) +
  geom_density(adjust=1.5, alpha=.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("Strength") + xlab("Motion effect (LME-t)") +
  theme_bw()
# Degree density
d3 <- ggplot(data=gt_data, aes(x=deg_hm, group=GSR, fill=GSR)) +
  geom_density(adjust=1.5, alpha=.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("Degree") + xlab("Motion effect (LME-t)") +
  theme_bw()
# Clustering coefficient density
d4 <- ggplot(data=gt_data, aes(x=clus_hm, group=GSR, fill=GSR)) +
  geom_density(adjust=1.5, alpha=.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("Clustering") + xlab("Motion effect (LME-t)") +
  theme_bw()
# Efficiency density
d5 <- ggplot(data=gt_data, aes(x=eff_hm, group=GSR, fill=GSR)) +
  geom_density(adjust=1.5, alpha=.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("Efficiency") + xlab("Motion effect (LME-t)") +
  theme_bw()
# Path-length density
d6 <- ggplot(data=gt_data, aes(x=path_hm, group=GSR, fill=GSR)) +
  geom_density(adjust=1.5, alpha=.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("Path-Length") + xlab("Motion effect (LME-t)") +
  theme_bw()
# Arrange plots
d7 <- grid.arrange(d1, d2, d3, d4, d5, d6, nrow = 2, ncol=3)
# Save
outname <- "graph_motion_density.pdf"
#ggsave(filename = outname, plot = d7, device = "pdf", width = 8, height = 4)

# Arrange all plots
grid.arrange(d7, g3, nrow = 2)
