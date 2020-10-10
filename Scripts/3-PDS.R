# First of all, clear workspace
rm(list=ls())
# Save graphical parameters
op <- par()

# Read phenotypic data
info <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Phenotypic/phenotypic.csv")

# Load 'gamm4' package for GAMM implementation
if(require(gamm4)==0) install.packages("gamm4"); library(gamm4)

# Apply GAMM
fit <- gamm4::gamm4(PDS ~ + s(Age, k=4, by = Sex) + Sex,
                    data = info, random = ~ (1|BIDS.ID), REML = F)
# No NA's indices
nona <- which(!is.na(info$PDS))

# Store splines and standard error
info$PDS_pred <- info$PDS_se <- rep(as.numeric(NA), nrow(info))
info$PDS_pred[nona] <- predict.gam(fit$gam)
info$PDS_se[nona] <- predict.gam(fit$gam, se.fit = T)[[2]]
# Predict NA's
info$PDSe <- info$PDS
info$PDSe[which(is.na(info$PDS))] <- predict.gam(fit$gam, newdata = info[which(is.na(info$PDS)),])

# Plot PSD vs Age
if(require(ggplot2)==0) install.packages("ggplot2"); library(ggplot2)
ggplot(data = info, mapping = aes(x = Age)) + 
    geom_line(aes(y = PDSe, group = BIDS.ID, color = Sex), size=.25) + 
    geom_point(aes(y = PDSe, color = Sex)) + 
    geom_ribbon(data = info[!is.na(info$PDS_pred),], 
                aes(group = Sex,
                    ymax = PDS_pred + 1.96 * PDS_se,
                    ymin = PDS_pred - 1.96 * PDS_se,
                    linetype = NA), alpha = 0.2) + 
    geom_line(data = info[!is.na(info$PDS_pred),],
              aes(y = PDS_pred, color = Sex), size = 1.5) + 
    theme_light() +
    xlab("Age (years)") +
    ylab("PDS") +
    ggtitle("Pubertal Development Scale")
