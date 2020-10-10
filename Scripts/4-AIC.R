# First of all, clear workspace
rm(list=ls())
# Save graphical parameters
op <- par()

# Check if AIC file exists
aic_file <- "https://raw.githubusercontent.com/BrainMapINB/Puberty/main/AIC/AIC_results.csv"

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
  source("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Scripts/Auxiliary/cost_th.R")
  source("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Scripts/Auxiliary/cost_th_ConnMx.R")
  # Load 'gamm4' package for GAMM analyses
  if(!require(gamm4)) install.packages("gamm4"); library(gamm4)
  # Load 'psych' package for Fisher r-to-z transformation
  if(!require(psych)) install.packages("psych"); library(psych)
  
  # Graph names
  grp_var <- c("D_net","C_net","L_net","E_net")
  
  # Explore a set of thresholds
  cost <- seq(from = 0.01, to = 0.48, by = 0.01)
  
  # Store results into empty matrix
  DF <- matrix(data = 0, nrow = length(cost), ncol = 139)
  DF[,1] <- cost
  # Set column names
  colnames(DF) <- c("cost","thr_avg","thr_sd",
                    
                    "d_GAMMsAge_AIC","d_GAMMsAge_sex_b","d_GAMMsAge_sex_t","d_GAMMsAge_sex_p",
                    "d_GAMMsAge_edf","d_GAMMsAge_f","d_GAMMsAge_p",
                    "d_GAMMsAgeSex_AIC","d_GAMMsAgeSex_sex_b","d_GAMMsAgeSex_sex_t","d_GAMMsAgeSex_sex_p",
                    "d_GAMMsAgeSexF_edf","d_GAMMsAgeSexF_f","d_GAMMsAgeSexF_p",
                    "d_GAMMsAgeSexM_edf","d_GAMMsAgeSexM_f","d_GAMMsAgeSexM_p",
                    "d_GAMMsPDS_AIC","d_GAMMsPDS_sex_b","d_GAMMsPDS_sex_t","d_GAMMsPDS_sex_p",
                    "d_GAMMsPDS_edf","d_GAMMsPDS_f","d_GAMMsPDS_p",
                    "d_GAMMsPDSSex_AIC","d_GAMMsPDSSex_sex_b","d_GAMMsPDSSex_sex_t","d_GAMMsPDSSex_sex_p",
                    "d_GAMMsPDSSexF_edf","d_GAMMsPDSSexF_f","d_GAMMsPDSSexF_p",
                    "d_GAMMsPDSSexM_edf","d_GAMMsPDSSexM_f","d_GAMMsPDSSexM_p",
                    
                    "c_GAMMsAge_AIC","c_GAMMsAge_sex_b","c_GAMMsAge_sex_t","c_GAMMsAge_sex_p",
                    "c_GAMMsAge_edf","c_GAMMsAge_f","c_GAMMsAge_p",
                    "c_GAMMsAgeSex_AIC","c_GAMMsAgeSex_sex_b","c_GAMMsAgeSex_sex_t","c_GAMMsAgeSex_sex_p",
                    "c_GAMMsAgeSexF_edf","c_GAMMsAgeSexF_f","c_GAMMsAgeSexF_p",
                    "c_GAMMsAgeSexM_edf","c_GAMMsAgeSexM_f","c_GAMMsAgeSexM_p",
                    "c_GAMMsPDS_AIC","c_GAMMsPDS_sex_b","c_GAMMsPDS_sex_t","c_GAMMsPDS_sex_p",
                    "c_GAMMsPDS_edf","c_GAMMsPDS_f","c_GAMMsPDS_p",
                    "c_GAMMsPDSSex_AIC","c_GAMMsPDSSex_sex_b","c_GAMMsPDSSex_sex_t","c_GAMMsPDSSex_sex_p",
                    "c_GAMMsPDSSexF_edf","c_GAMMsPDSSexF_f","c_GAMMsPDSSexF_p",
                    "c_GAMMsPDSSexM_edf","c_GAMMsPDSSexM_f","c_GAMMsPDSSexM_p",
                    
                    "l_GAMMsAge_AIC","l_GAMMsAge_sex_b","l_GAMMsAge_sex_t","l_GAMMsAge_sex_p",
                    "l_GAMMsAge_edf","l_GAMMsAge_f","l_GAMMsAge_p",
                    "l_GAMMsAgeSex_AIC","l_GAMMsAgeSex_sex_b","l_GAMMsAgeSex_sex_t","l_GAMMsAgeSex_sex_p",
                    "l_GAMMsAgeSexF_edf","l_GAMMsAgeSexF_f","l_GAMMsAgeSexF_p",
                    "l_GAMMsAgeSexM_edf","l_GAMMsAgeSexM_f","l_GAMMsAgeSexM_p",
                    "l_GAMMsPDS_AIC","l_GAMMsPDS_sex_b","l_GAMMsPDS_sex_t","l_GAMMsPDS_sex_p",
                    "l_GAMMsPDS_edf","l_GAMMsPDS_f","l_GAMMsPDS_p",
                    "l_GAMMsPDSSex_AIC","l_GAMMsPDSSex_sex_b","l_GAMMsPDSSex_sex_t","l_GAMMsPDSSex_sex_p",
                    "l_GAMMsPDSSexF_edf","l_GAMMsPDSSexF_f","l_GAMMsPDSSexF_p",
                    "l_GAMMsPDSSexM_edf","l_GAMMsPDSSexM_f","l_GAMMsPDSSexM_p",
                    
                    "e_GAMMsAge_AIC","e_GAMMsAge_sex_b","e_GAMMsAge_sex_t","e_GAMMsAge_sex_p",
                    "e_GAMMsAge_edf","e_GAMMsAge_f","e_GAMMsAge_p",
                    "e_GAMMsAgeSex_AIC","e_GAMMsAgeSex_sex_b","e_GAMMsAgeSex_sex_t","e_GAMMsAgeSex_sex_p",
                    "e_GAMMsAgeSexF_edf","e_GAMMsAgeSexF_f","e_GAMMsAgeSexF_p",
                    "e_GAMMsAgeSexM_edf","e_GAMMsAgeSexM_f","e_GAMMsAgeSexM_p",
                    "e_GAMMsPDS_AIC","e_GAMMsPDS_sex_b","e_GAMMsPDS_sex_t","e_GAMMsPDS_sex_p",
                    "e_GAMMsPDS_edf","e_GAMMsPDS_f","e_GAMMsPDS_p",
                    "e_GAMMsPDSSex_AIC","e_GAMMsPDSSex_sex_b","e_GAMMsPDSSex_sex_t","e_GAMMsPDSSex_sex_p",
                    "e_GAMMsPDSSexF_edf","e_GAMMsPDSSexF_f","e_GAMMsPDSSexF_p",
                    "e_GAMMsPDSSexM_edf","e_GAMMsPDSSexM_f","e_GAMMsPDSSexM_p"
  )
  
  for(hh in 1:length(cost)){
    
    # Print progress
    print(paste0("Cost ",cost[hh]))
    
    # Apply cost threshold
    thr <- apply(ConnMx,3,function(x) cost_th(x,cost[hh]))
    # Save threshold average
    DF[hh,2] <- mean(thr); DF[hh,3] <- sd(thr)
    # Apply thresholds to Connectivity Matrix
    ConnMx_th <- apply(ConnMx,3,function(x) cost_th_ConnMx(x,cost[hh]))
    ConnMx_th <- array(ConnMx_th,dim=ConnMx_dim)
    
    # Apply complex networks measures.
    # Substract diagonal
    for(ii in 1:ConnMx_dim[3]) diag(ConnMx_th[,,ii]) <- 0
    # Fisher's z transformation
    ConnMx_th <- fisherz(ConnMx_th)
    if(max(ConnMx_th) == Inf) ConnMx_th[ConnMx_th==Inf] <- max(ConnMx_th[ConnMx_th!=Inf])
    
    # Weighted-Degree
    w_deg <- apply(ConnMx_th,3:2,sum)
    # Clustering Coefficiente (per node) - weighted
    tic <- Sys.time()
    w_CC <- t(sapply(1:ConnMx_dim[3], function(x) iCC_w(ConnMx_th[,,x])))
    if(sum(is.na(w_CC))>0)  w_CC[is.na(w_CC)] <- 0
    toc <- Sys.time()
    print("Clustering:"); print(toc-tic)
    # Average distances (per node) - weighted
    tic <- Sys.time()
    w_L <- t(sapply(1:ConnMx_dim[3], function(x) iL_w(ConnMx_th[,,x])))
    toc <- Sys.time()
    print("Path-length:"); print(toc-tic)
    # Inverse average distances (per node) - weighted
    tic <- Sys.time()
    w_E <- t(sapply(1:ConnMx_dim[3], function(x) iEff_w(ConnMx_th[,,x])))
    toc <- Sys.time()
    print("Efficiency:"); print(toc-tic)
    
    # Whole-brain graph-theory measures along the sample
    D_net <- rowMeans(w_deg)
    C_net <- rowMeans(w_CC)
    L_net <- rowMeans(w_L)
    E_net <- rowMeans(w_E)
    
    # Apply linear models at the group level
    info$D_net <- D_net
    info$C_net <- C_net
    info$L_net <- L_net
    info$E_net <- E_net
    
    # Test different models
    for(kk in 1:length(grp_var)){
      
      # GAMM model (age spline)
      grp_form <- as.formula(paste0(grp_var[kk], " ~ s(Age, k=4) + Sex + FDRMS + Coil"))
      fitGAMM <- tryCatch(gamm4(grp_form,
                                random = ~ (1|BIDS.ID),
                                data = info,
                                REML = F),
                          error=function(e) NA)
      # GAMM model (age*sex splines)
      grp_form <- as.formula(paste0(grp_var[kk], " ~ s(Age, k=4, by=Sex) + Sex + FDRMS + Coil"))
      fitGAMMi <- tryCatch(gamm4(grp_form,
                                 random = ~ (1|BIDS.ID),
                                 data = info,
                                 REML = F),
                           error=function(e) NA)
      # GAMM model (PDS spline)
      grp_form <- as.formula(paste0(grp_var[kk], " ~ s(PDSe, k=4) + Sex + FDRMS + Coil"))
      fitPDS <- tryCatch(gamm4(grp_form,
                               random = ~ (1|BIDS.ID),
                               data = info,
                               REML = F),
                         error=function(e) NA)
      # GAMM model (PDS*sex splines)
      grp_form <- as.formula(paste0(grp_var[kk], " ~ s(PDSe, k=4, by=Sex) + Sex + FDRMS + Coil"))
      fitPDSi <- tryCatch(gamm4(grp_form,
                                random = ~ (1|BIDS.ID),
                                data = info,
                                REML = F),
                          error=function(e) NA)# Store statistics
      # GAMM-Age
      DF[hh,4+(kk-1)*34] <- tryCatch(AIC(fitGAMM$mer), error=function(e) NA)
      DF[hh,5:7+(kk-1)*34] <- tryCatch(summary(fitGAMM$gam)$p.table[2,c(1,3,4)],
                                       error=function(e) NA)
      DF[hh,8:10+(kk-1)*34] <- tryCatch(summary(fitGAMM$gam)$s.table[,c(1,3,4)],
                                        error=function(e) NA)
      # GAMM-AgeBySex
      DF[hh,11+(kk-1)*34] <- tryCatch(AIC(fitGAMMi$mer), error=function(e) NA)
      DF[hh,12:14+(kk-1)*34] <- tryCatch(summary(fitGAMMi$gam)$p.table[2,c(1,3,4)],
                                         error=function(e) NA)
      DF[hh,15:20+(kk-1)*34] <- tryCatch(c(t(summary(fitGAMMi$gam)$s.table[,c(1,3,4)])),
                                         error=function(e) NA)
      # GAMM-PDS
      DF[hh,21+(kk-1)*34] <- tryCatch(AIC(fitPDS$mer), error=function(e) NA)
      DF[hh,22:24+(kk-1)*34] <- tryCatch(summary(fitPDS$gam)$p.table[2,c(1,3,4)],
                                         error=function(e) NA)
      DF[hh,25:27+(kk-1)*34] <- tryCatch(summary(fitPDS$gam)$s.table[,c(1,3,4)],
                                         error=function(e) NA)
      # GAMM-PDSbySex
      DF[hh,28+(kk-1)*34] <- tryCatch(AIC(fitPDSi$mer), error=function(e) NA)
      DF[hh,29:31+(kk-1)*34] <- tryCatch(summary(fitPDSi$gam)$p.table[2,c(1,3,4)],
                                         error=function(e) NA)
      DF[hh,32:37+(kk-1)*34] <- tryCatch(c(t(summary(fitPDSi$gam)$s.table[,c(1,3,4)])),
                                         error=function(e) NA)
      
    }# for(kk in 1:length(grp_var)) 
    
  }# for(hh in 1:length(cost))
  
  # Save results
  outname <- "AIC_results.csv"
  write.csv(DF, outname, quote = F, row.names = F)
}

# Read AIC results
all_fits <- read.csv(aic_file)

# Plot results
# Load 'ggplot2' package
if(!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)
# Model labels
mod_labs <- c("Age","Age-Sex","PDS","PDS-Sex")
if(!all(rowSums(abs(all_fits[,-(1:3)]))!=0)) all_fits <- all_fits[which(rowSums(abs(all_fits[,-(1:3)]))!=0),]
# Compute threshold-based FDR
all_fitsFDR <- all_fits
all_fitsFDR[,grep("_p",names(all_fits))] <- apply(all_fits[,grep("_p",names(all_fits))],
                                                      2,
                                                      function(x) p.adjust(x, method = "fdr"))
# Set initial variables
grph_names <- c("d_","c_","l_","e_")
best_aic <- matrix(0, nrow = nrow(all_fits), ncol = length(grph_names))
colnames(best_aic) <- paste0(grph_names,"bestAIC")
best_sigFDR <- best_sig <- best_aic

# Store best models
for(ee in 1:length(grph_names)){
  
  # Extract columns
  aic_colidx <- intersect(grep(paste0(grph_names[ee],"G"), names(all_fits)),grep("AIC", names(all_fits)))
  sig_colidx <- intersect(grep(paste0(grph_names[ee],"G"), names(all_fits)),grep("_p", names(all_fits)))
  sig_colidx <- sig_colidx[grep("_sex_", names(all_fits)[sig_colidx], invert = T)]
  # Get best fit
  best_aic[,ee] <- sapply(1:nrow(all_fits), function(x){
    y1 <- which.min(all_fits[x,aic_colidx]) # Get best AIC position
    if(length(y1) == 0) y1 <- NA
    return(y1)
  }
  )
  # Get best spline
  best_sig[,ee] <- sapply(1:nrow(all_fits), function(x){
    y1 <- which.min(all_fits[x,aic_colidx]) # Get best AIC position
    if(length(y1) == 0) y2 <- NA else{
      if(y1 == 1) y2 <- all_fits[x,sig_colidx[1]]
      if(y1 == 2) y2 <- min(all_fits[x,sig_colidx[2:3]])
      if(y1 == 3) y2 <- all_fits[x,sig_colidx[4]]
      if(y1 == 4) y2 <- min(all_fits[x,sig_colidx[5:6]])
    }
    return(y2)
  }
  )
  # Get best spline (FDR)
  best_sigFDR[,ee] <- sapply(1:nrow(all_fitsFDR), function(x){
    y1 <- which.min(all_fitsFDR[x,aic_colidx]) # Get best AIC position
    if(length(y1) == 0) y2 <- NA else{
      if(y1 == 1) y2 <- all_fitsFDR[x,sig_colidx[1]]
      if(y1 == 2) y2 <- min(all_fitsFDR[x,sig_colidx[2:3]])
      if(y1 == 3) y2 <- all_fitsFDR[x,sig_colidx[4]]
      if(y1 == 4) y2 <- min(all_fitsFDR[x,sig_colidx[5:6]])
    }
    return(y2)
  }
  )
}

# Plot best results
# Reshape data.frame
rsh_df <- as.data.frame(rep(all_fits$cost,length(grph_names))) # cost
rsh_df <- cbind(rsh_df,rep(1:length(grph_names), each = nrow(all_fits))) # varnames
rsh_df <- cbind(rsh_df,c(best_aic)) # best AIC
rsh_df <- cbind(rsh_df,c(best_sig)) # best p-value
rsh_df <- cbind(rsh_df,c(best_sigFDR)) # best pFDR
rsh_df <- cbind(rsh_df,rep(1,nrow(rsh_df))) # symbol size
rsh_df <- cbind(rsh_df,rep(1,nrow(rsh_df))) # symbol alpha
names(rsh_df) <- c("cost","var","AIC","p","pFDR","size","alpha")
# Set symbol sizes
rsh_df$size[rsh_df$p < 0.05] <- 1.5
rsh_df$size[rsh_df$pFDR < 0.05] <- 3
# Set symbol alphas
rsh_df$alpha[rsh_df$pFDR > 0.05] <- 0.5
rsh_df$alpha[rsh_df$p > 0.05] <- 0.1
rsh_df$var <- factor(rsh_df$var)
levels(rsh_df$var) <- c("Degree","Clustering","Path-Lengh","Efficiency")
# Relevel best-fit labels
rsh_df$AIC <- factor(rsh_df$AIC)
levels(rsh_df$AIC) <- mod_labs[as.integer(levels(rsh_df$AIC))]
# Plot best AIC
ggAIC <- ggplot(data=rsh_df, aes(x=var, y=cost)) + 
    stat_identity(aes(color=AIC, size=size)) +
    labs(color = "Model", size = "Significance") +
    ylab("Connectivity density (proportion)") + xlab("") +
    ggtitle("Lower Akaike Information Criterion (AIC)") + coord_flip()  + 
    scale_color_brewer(palette = "Dark2") +
    theme_light()
show(ggAIC)
# Significance: p > 0.05 (1); p < 0.05 (1.5); pFDR < 0.05 (3)

# Save plot
outfile <- "Lower_AIC.pdf"
#ggsave(outfile, ggAIC, device = "pdf")
