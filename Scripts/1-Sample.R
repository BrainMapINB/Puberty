# First of all, clear workspace
rm(list=ls())

# Read phenotypic data
info <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Puberty/main/Phenotypic/phenotypic.csv")

# Set color
if(!require(scales)) install.packages("scales"); library(scales)
huecol <- hue_pal()(2)

# Plot histogram of age by sex
# Take initial variables
edad <- info$Age
sexo <- info$Sex
id <- factor(info$BIDS.ID)
tp <- table(table(id))
# Plot histogram
seq1y <- seq(floor(min(edad)), ceiling(max(edad)), 1)
hist(edad, breaks = seq1y, col = huecol[1], axes = F,
     main = "Histogram of age by sex", xlab = "Age (years)")
hist(edad[which(sexo=="M")], breaks = seq1y, col = huecol[2], axes = F, add = T)
seq2y <- seq(floor(min(edad)),ceiling(max(edad)),2)
axis(1,seq1y,F);axis(1,seq2y)
axis(2,seq(0,25,5),las=1)
legend("topright", c("F","M"), pch = 15, col = huecol)
text(17,20,paste("N =",length(table(id)),"\n2tp =",tp[2],"\n3tp =",tp[3]),pos=4)

# Longitudinal scatter-plot
# Sort original variables
info <- info[order(info$Age),]
edad <- info$Age
id <- as.character(info$BIDS.ID)
# Extract labels and frequencies
id_lab <- unique(id)
id_num <- sapply(id_lab,function(x)sum(id==x))
n <- length(id_lab)
sexo <- info$Sex
# Plot canvas
plot(edad, 1:length(edad), type = "n", axes = F,
     ylim = c(0,n), xlim = c(floor(min(edad)), ceiling(max(edad))),
     main = "Longitudinal plot", ylab = "Subject", xlab = "Age (years)")
axis(1, seq1y, F); axis(1,seq2y)
axis(2, seq(0, 10*ceiling(n/10), 10),las=1)
# Draw points
for(ii in 1:n){
  
  # Find sessions
  pos <- which(id==id_lab[ii])
  
  # Extract sex character
  ifelse(sexo[pos[1]] == "F", sexchar <- huecol[1], sexchar <- huecol[2])
  
  # Draw points
  if(length(pos)==1) points(edad[pos], ii, pch = 21, bg = sexchar, cex = 0.8, col = "black") else{
    lines(c(edad[min(pos)],edad[max(pos)]),c(ii,ii), col = sexchar)
    for(jj in pos) points(edad[jj], ii, pch = 21, bg = sexchar, cex = 0.8, col = "black")
  }
}
# Legend
legend(min(edad), n, c("F","M"), pch = 21, pt.bg = huecol, bty = "n")
