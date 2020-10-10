# First of all, clear workspace
rm(list=ls())
# Save graphical parameters
op <- par()

# Check if whole-brain file exists
tr_file <- "https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Trajectories/connectome.csv"

# If you prefer to generate those values follow the steps within this 'if' command
# (but be aware that some computations may be take a while)
if(0){
  
  # Read phenotypic data
  info <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Phenotypic/phenotypic.csv")
  # Remove those sessions with excesive motion artifact
  info <- info[which(as.logical(info$QC)),]
  # Refactorize subjects' ID
  info$BIDS.ID <- factor(info$BIDS.ID)
  
  # Extract time-series files (those with GSR)
  ts_files <- paste0("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Brain/GSR/sub-",
                     info$BIDS.ID,"/ses-",info$Session,"/pp150v_wGSR_NIHPD_MNI_P264_ts.txt")
  
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
  
  # Load graph theory formulas
  source("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Scripts/Auxiliary/inet_w.R")
  # Load cost-threshold functions
  source("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Scripts/Auxiliary/cost_th_ConnMx.R")
  # Load 'psych' package for Fisher r-to-z transformation
  if(!require(psych)) install.packages("psych"); library(psych)
  
  # Apply thresholds to Connectivity Matrix
  ConnMx_th <- apply(ConnMx,3,function(x) cost_th_ConnMx(x,0.25))
  ConnMx_th <- array(ConnMx_th,dim=ConnMx_dim)
  # Substract diagonal
  for(ii in 1:ConnMx_dim[3]) diag(ConnMx_th[,,ii]) <- 0
  # Fisher's z transformation
  ConnMx_th <- fisherz(ConnMx_th)
  
  # GRAPH THEORY
  # Weighted-Degree
  w_deg <- apply(ConnMx_th,3:2,sum)
  # Clustering Coefficient (Barrat's formula)
  w_CC <- t(sapply(1:ConnMx_dim[3], function(x) iCC_w(ConnMx_th[,,x])))
  # Global Efficiency (inverse of Dijkstra's formula)
  w_E <- t(sapply(1:ConnMx_dim[3], function(x) iEff_w(ConnMx_th[,,x])))
  # Minimum path-length distances (Dijkstra's formula)
  w_L <- t(sapply(1:ConnMx_dim[3], function(x) iL_w(ConnMx_th[,,x])))
  
  # Compute whole-brain variables
  info$D_net <- rowMeans(w_deg)
  info$C_net <- rowMeans(w_CC)
  info$E_net <- rowMeans(w_E)
  info$L_net <- rowMeans(w_L)
  
  # Extract Functional Network inferences
  atlas <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Atlas/P264/power264.csv")
  
  # Labels
  grph_mat <- c("w_deg","w_CC","w_E","w_L")
  grph_lab <- c("D","C","E","L")
  
  # Apply Functional Network averages
  for(ii in 1:length(grph_mat)){
    fn_mat <- t(aggregate(x = t(get(grph_mat[ii])), by = list(atlas$Net), FUN = mean)[,-1])
    colnames(fn_mat) <- paste0(grph_lab[ii],"_",levels(atlas$Net))
    info <- cbind(info,fn_mat[,-12]) # Remove 'Uncertain' elements
  }
  
  # Save results
  outfile <- "connectome.csv"
  write.csv(info, outfile, quote = F, row.names = F)
}

###########################################################################
# Whole-brain
###########################################################################

# Plot result
brain_net <- read.csv(tr_file)
# Load 'gamm4', 'ggplot2' & 'gridExtra' packages
if(!require(gamm4)) install.packages("gamm4"); library(gamm4)
if(!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)
if(!require(gridExtra)) install.packages("gridExtra"); library(gridExtra)

# Graph labels
grph_var <- c("D_net", "C_net", "E_net", "L_net")
grph_lab <- c("Degree", "Clustering", "Efficiency", "Path-length")

# Create empty matrix to store results
DF <- matrix(data = as.numeric(NA), nrow = length(grph_var), ncol = 6)
rownames(DF) <- grph_var
colnames(DF) <- c("IPnum","IPval","slope","EDF","Fval","pval")

# Generate scatterplot for each graph modality
for(ii in 1:length(grph_var)){
  
  # Select column
  grph_frm <- as.formula(paste0(grph_var[ii], "~ + s(PDSe, k=4) + Sex + FDRMS + Coil"))
  # Apply GAMM model
  fit <- gamm4::gamm4(formula = grph_frm, data = brain_net,
                      random = ~ (1|BIDS.ID), REML = F)
  
  # Store GAMM results
  DF[ii,4:6] <- summary(fit$gam)$s.table[,c(1,3,4)]
  
  # Generate residuals
  res_frm <- as.formula(paste0(grph_var[ii], "~ Sex + FDRMS + Coil + (1|BIDS.ID)"))
  res_lme <- lmer(res_frm, data = brain_net, REML = F)
  res_b <- coefficients(res_lme)[[1]]
  brain_net$res <- residuals(res_lme)+res_b$`(Intercept)`[match(brain_net$BIDS.ID, rownames(res_b))]
  plt_form <- as.formula("res ~ + s(PDSe, k=4)")
  brain_net$preds <- predict(gamm4::gamm4(plt_form, data = brain_net, random = ~ (1|BIDS.ID), REML = F)$gam)
  brain_net$preds_se <- predict(gamm4::gamm4(plt_form, data = brain_net, random = ~ (1|BIDS.ID), REML = F)$gam,
                                se.fit = T)[[2]]
  # Generate plot
  gSC <- ggplot(data = brain_net, mapping = aes(x = PDSe)) + 
    geom_line(aes(y = res, group = BIDS.ID, color = Sex), size=.25) + 
    geom_point(aes(y = res, color = Sex)) + 
    geom_ribbon(aes(ymax = preds + 1.96 * preds_se,
                    ymin = preds - 1.96 * preds_se,
                    linetype = NA), alpha = 0.2) + 
    geom_line(aes(y = preds), size = 1.5) +
    ylab("Residuals") + xlab("PDS") + ggtitle(grph_lab[ii]) +
    theme_light()
  if(ii == 1) g1 <- gSC else if(ii == 2) g2 <- gSC else if(ii == 3) g3 <- gSC else if(ii == 4) g4 <- gSC
  
  # Extract standardized residuals
  brain_net$res <- scale(residuals(res_lme)+res_b$`(Intercept)`[match(brain_net$BIDS.ID, rownames(res_b))])
  plt_frm <- as.formula("res ~ + s(PDSe, k=4)")
  fit <- gamm4::gamm4(plt_form, data = brain_net, random = ~ (1|BIDS.ID), REML = F)
  
  # Store slope and point of inflection
  brain_net$preds <- predict(fit$gam)
  lo <- loess(preds~PDSe, data = brain_net)
  x <- seq(1,4,0.01)
  y <- predict(lo,x)
  #plot(x,y, type = "l")
  infl <- c(FALSE, diff(diff(y)>0)!=0)
  DF[ii,1] <- sum(infl)
  DF[ii,2] <- 0
  if(DF[ii,1]>0){
    # Store the last inflection point
    DF[ii,2] <- x[infl][DF[ii,1]]
    # Compute slope after inflection point
    y <- y[x>=DF[ii,2]]
    x <- x[x>=DF[ii,2]]
  }
  # Compute slope
  DF[ii,3] <- coefficients(lm(y~x))[2]

} #for(ii in 1:length(grph_var))

# Arrange plots
g5 <- grid.arrange(g1, g2, g3, g4, nrow = 2, ncol=2)
#ggsave("brain_trajectories.pdf", g5, device = "pdf", width = 8, height = 5)

# Save results
outfile <- "whole-brain.csv"
#write.csv(x = DF, file = outfile, quote = F)

###########################################################################
# Functional Networks
###########################################################################

# If you prefer to generate the NIfTI files follow the steps within this 'if' command
if(0){
  # Read ROI and Network labels
  atlas <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Atlas/P264/power264.csv")
  
  # Compute GAMM models each graph theory variable and functional network
  # Find Functional Network columns
  fn_idx <- grep(paste(levels(atlas$Net), collapse = "|"), names(brain_net))
  # Create empty matrix to store results
  DF <- matrix(data = as.numeric(NA), nrow = length(fn_idx), ncol = 7)
  rownames(DF) <- names(brain_net)[fn_idx]
  colnames(DF) <- c("IPnum","IPval","slope","EDF","Fval","pval","pFDR")
  
  # For each functional network
  for(nn in 1:length(fn_idx)){
    
    # Select column
    fn_frm <- as.formula(paste0(names(brain_net)[fn_idx[nn]],
                                "~ + s(PDSe, k=4) + Sex + FDRMS + Coil"))
    # Apply GAMM model
    fit <- gamm4::gamm4(formula = fn_frm, data = brain_net,
                        random = ~ (1|BIDS.ID), REML = F)
    
    # Store GAMM results
    DF[nn,4:6] <- summary(fit$gam)$s.table[,c(1,3,4)]
    
    # Extract standardized residuals
    res_frm <- as.formula(paste0(names(brain_net)[fn_idx[nn]],
                                 "~ Sex + FDRMS + Coil + (1|BIDS.ID)"))
    res_lme <- lmer(formula = res_frm, data = brain_net, REML = F)
    res_b <- coefficients(res_lme)[[1]]
    brain_net$res <- scale(residuals(res_lme)+res_b$`(Intercept)`[match(brain_net$BIDS.ID, rownames(res_b))])
    plt_frm <- as.formula("res ~ + s(PDSe, k=4)")
    fit <- gamm4::gamm4(plt_form, data = brain_net, random = ~ (1|BIDS.ID), REML = F)
    
    # Store slope and point of inflection
    brain_net$preds <- predict(fit$gam)
    lo <- loess(preds~PDSe, data = brain_net)
    x <- seq(1,4,0.01)
    y <- predict(lo,x)
    #plot(x,y, type = "l")
    infl <- c(FALSE, diff(diff(y)>0)!=0)
    DF[nn,1] <- sum(infl)
    DF[nn,2] <- 0
    if(DF[nn,1]>0){
      # Store the last inflection point
      DF[nn,2] <- x[infl][DF[nn,1]]
      # Compute slope after inflection point
      y <- y[x>=DF[nn,2]]
      x <- x[x>=DF[nn,2]]
    }
    # Compute slope
    DF[nn,3] <- coefficients(lm(y~x))[2]
    
  }
  
  # Compute FDR correction and volume assignation at the Functional Network level
  if(!require(neurobase)) install.packages("neurobase"); library(neurobase)
  # Get the number of functional networks (excluding Uncertain ROIs)
  fn_num <- nlevels(atlas$Net)-1
  for(gg in 1:length(grph_var)){
    # Get p-value indices
    p_idx <- 1:fn_num + (gg-1)*fn_num
    
    # Apply FDR
    DF[p_idx,7] <- p.adjust(p = DF[p_idx,6], method = "fdr")
    
    # Lastly, insert slope values into a brain volume
    # It is necesary to download this file first
    # https://github.com/BrainMapINB/Puberty/blob/main/Atlas/P264/FuntionalNetworks_P264.nii.gz
    # Read template
    nii <- readnii("FuntionalNetworks_P264.nii.gz")
    for(nn in 1:fn_num){
      if(nn < 12) nii[nii==nn] <- DF[p_idx[nn],3] else nii[nii==nn+1] <- DF[p_idx[nn],3]
    }
    nii[nii==12] <- 0
    outfile <- paste0("PDS_slope_",grph_lab[gg])
    #writenii(nim = nii, filename = outfile)
  }
  
  # Save results
  outfile <- "functional-networks.csv"
  #write.csv(x = DF, file = outfile, quote = F)
  
}
