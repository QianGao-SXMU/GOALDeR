
#' dual-weigt distance correlation (DWDC) for selecting optimal lambda
#'
#' @param DataM Data that include outcome, treatment and all of covariates
#' @param varlist names of potential confounders including known and unknown
#' @param trt.var names of treatment
#' @param wgt balancing weights
#' @param beta the conditional distance correlation between covariate and outcome given treatment
#'
#' @return wAMD is a vector of dual-weight distance correlation coefficients
#' @export
#'
#' @examples
DWDC_function <- function(DataM,varlist,trt.var,wgt,beta){
  diff_vec <- rep(NA,length(beta))
  names(diff_vec) <- varlist
  for(jj in 1:length(varlist)){
    diff_vec[jj]<-abs(wdcor(x=DataM[,trt.var],y=DataM[,varlist[jj]],w=DataM[,wgt]))
  }
  wdiff_vec = diff_vec * abs(beta)
  wAMD = c(sum(wdiff_vec))
  ret = list( diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD )
  return(ret)
}
