#' The Generalized Outcome-Adaptive LASSO and Doubly Robust Estimate (GOALDeR)
#'
#' @param data Data that include outcome, treatment and all of covariates
#' @param var.list names of potential unknown confounders
#' @param covar names of known confounders
#' @param Trt Trt names of treatment
#' @param out out names of outcome
#' @param lambda_vec a vector of possible lambda values and defualt is c(-10,-5,-2,-1.5,-1.25,-1, -0.75, -0.5, -0.25, 0.25, 0.49)
#' @param gamma_convergence the value of gamma converagence and defualt is 2
#'
#' @return ATE is the linear dose-response function estimated using doubly robust estimator method.
#' @export
#'
#' @examples
GOALDeR<-function(data,var.list,covar,Trt="Trt",out="Y",lambda_vec=c(-10,-5,-2,-1.5,-1.25,-1, -0.75, -0.5, -0.25, 0.25, 0.49),gamma_convergence=2){
	Data<-data
	n<-dim(Data)[1]
	temp.mean <- colMeans(Data[,var.list])
	Temp.mean <- matrix(temp.mean,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
	Data[,var.list] <- Data[,var.list] - Temp.mean
	temp.sd <- apply(Data[var.list],FUN=sd,MARGIN=2)
	Temp.sd <- matrix(temp.sd,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
	Data[,var.list] <- Data[,var.list] / Temp.sd
	rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))
	betaXY_cor<-NA
	for (i in 1:length(var.list)) {
		var.name<-var.list[i]
		betaXY_cor[i] <-as.numeric(cdcor(Data[,var.name],Data[,out],Data[,Trt])$statistic)
	}
	betaXY<-betaXY_cor/max(betaXY_cor)
	names(betaXY)<-var.list
	lambda_vec <- lambda_vec
	names(lambda_vec) <- as.character(lambda_vec)
	gamma_convergence_factor <- gamma_convergence
	gamma_vals <- 2*( gamma_convergence_factor - lambda_vec + 1 )
	names(gamma_vals) <- names(lambda_vec)
	wAMD_vec=rep(NA, length(lambda_vec))
	DR_estimated<-vector(mode = "list", length(lambda_vec))
	DR_pseudo_outcome<-matrix(NA,nrow=dim(Data)[1],ncol=length(lambda_vec))
	names(DR_estimated)=names(wAMD_vec)=names(lambda_vec)
	colnames(DR_pseudo_outcome)=names(lambda_vec)
	Dcow.var<-list()

	w.full.form <- formula(paste(Trt,"~",paste(c(covar,var.list),collapse="+")))
	for( lil in names(lambda_vec) ){
		v<-which(names(lambda_vec)==lil)
		il = lambda_vec[lil]
		ig = gamma_vals[lil]
  		oal_pen <- adaptive.lasso(lambda=n^(il),al.weights = abs(betaXY)^(-ig))
		logit_oal <- lqa.formula(w.full.form, data=Data, penalty=oal_pen,
				family=gaussian())
	# save propensity score coefficients
		coeff_XA <- coef(logit_oal)[var.list]
		Dcow.var[[v]]<-names(coeff_XA)[which(round(coeff_XA,5)!=0)]
		names(Dcow.var)[[v]]<-paste("lambda",lil,sep="_")
	# create inverse probability of treatment weights
		if (length(Dcow.var[[v]])!=0) {
			w.model<-formula(paste(Trt,"~",paste(Dcow.var[[v]],collapse="+")))
			Dcow_fit<-independence_weights(Data[,Trt], Data[,Dcow.var[[v]]])
			Data[,paste("w",lil,sep="")]<-Dcow_fit$weights
			} else {
			Data[,paste("w",lil,sep="")]<-1
			}
			wAMD_vec[lil] <- DWDC_function(DataM=Data,varlist=names(betaXY),trt.var=Trt,
				wgt=paste("w",lil,sep=""),beta=(betaXY_cor)^2)$wAMD
			set.seed(1)
			DR<-kennedy_etal(Y=Data[,out],A=Data[,Trt],X=Data[,var.list],Aseq=Data[,Trt],weights_est=Data[,paste("w",lil,sep="")])
			DR_estimated[[v]]<-DR$lm_estimated
			DR_pseudo_outcome[,v]<-DR$pseudo_out[,"pseudo_outcome"]
			rm(list=c("DR"))
	}
	ATE<-DR_estimated[[which.min(wAMD_vec)]]
	Svar<-Dcow.var[[which.min(wAMD_vec)]]
	lambda<-names(wAMD_vec)[which.min(wAMD_vec)]
	w_lil<-names(wAMD_vec)[which.min(wAMD_vec)]
	fw<-Data[,paste("w",w_lil,sep="")]
	GOALDeR_results<-list(
	ATE=ATE,
	selectedVar=Svar,
	lambda=lambda,
	fw=fw
	)
	return (GOALDeR_results)
}
