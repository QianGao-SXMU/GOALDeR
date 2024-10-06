#' Doubly robust estimator by Kennedy et al.,2017
#'
#' @param Y Outcome
#' @param A continuous treatment
#' @param X A matrix of covariates
#' @param Aseq continuous treatment
#' @param sl.lib methods combined by super learner and defualt is LASSO, XGBoost, Random Forest, and Support vector
#' @param weights_est balance weights
#'
#' @return lm_estimated is the estimated linear dose-response function
#' @export
#'
#' @examples
kennedy_etal <- function(Y, A, X, Aseq,sl.lib = c("SL.glmnet", "SL.randomForest", "SL.xgboost", "SL.svm"),weights_est)
{
	library(SuperLearner)
    library(locfit)
    library(randomForest)
    data<-cbind(Y,A,X)
    data <- as.data.frame(data)

    outcome_mod <- SuperLearner(Y = Y, X = as.data.frame(cbind(A,X)), SL.library = sl.lib)

    mu_x_a <- unname(outcome_mod$SL.predict)

    cmx <- colMeans(X)
    xm <- matrix(rep(cmx, nrow(data)), ncol = length(cmx), byrow = TRUE)

    pred_data <- cbind(A,xm)
    pred_data <- as.data.frame(pred_data)
    colnames(pred_data) <- c("A",colnames(X))

    mean_over_X <- function(t)
    {
        predict(outcome_mod, newdata = pred_data,onlySL = TRUE)$pred
    }

    if (any(weights_est > 500))
    {
        num_big <- sum(weights_est > 500)
        weights_est <- trim(weights_est, at = num_big)
    }

    pseudo_outcome <- weights_est * (Y - mu_x_a) + mean_over_X(A)


    dfx <- data.frame(Y = pseudo_outcome, TRT = A)
    locpoly_fit <- lm(Y ~ A, data = dfx)

    estimated<-summary(locpoly_fit)$coefficients
	pout<-data.frame(A, pseudo_outcome)
	kennedy_DR<-list(
		lmfit=locpoly_fit,
		lm_estimated=estimated,
		pseudo_out=pout
	)

    return(kennedy_DR)
}
