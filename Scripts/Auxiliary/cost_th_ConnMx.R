# Creator: zchuri
# Date: 10 November 2015
# Function: cost_th_ConnMx

# COST_TH_CONNMX(CONNMX,COST)
# This function applies a cost threshold to a connectivity matrix, values
# under the cost-threshold are set to zero
#
# -------------------------------------------------------------------------
# Inputs:
# CONNMX - square connectivity (or correlation matrix).
# COST - numeric value (proportion value)
# -------------------------------------------------------------------------

cost_th_ConnMx <- function(ConnMx,cost) {
  if(!is.numeric(cost)) stop("Please... insert a reasonable cost threshold")
  if(cost <= 0 || cost>=1) stop("Loc@!! Cost out of bounds!!")
  
  uptri <- ConnMx[upper.tri(ConnMx)]
  th <- quantile(uptri,probs = (1-cost))
  #print(th)
  ConnMx[which(ConnMx<th)] <- 0
  return(ConnMx)
}