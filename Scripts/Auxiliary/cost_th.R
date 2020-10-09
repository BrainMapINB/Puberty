# Creator: zchuri
# Date: 14 September 2019
# Function: cost_th

# COST_TH(CONNMX,COST)
# This function applies a cost threshold to a connectivity matrix, and
# returns the corresponding connectivity threshold.
#
# -------------------------------------------------------------------------
# Inputs:
# CONNMX - square connectivity (or correlation matrix).
# COST - numeric value (proportion value)
# -------------------------------------------------------------------------

cost_th <- function(ConnMx,cost) {
  if(!is.numeric(cost)) stop("Please... insert a reasonable cost threshold")
  if(cost <= 0 || cost>=1) stop("Loc@!! Cost out of bounds!!")
  
  uptri <- ConnMx[upper.tri(ConnMx)]
  th <- quantile(uptri,probs = (1-cost))
  #print(th)
  return(th)
}