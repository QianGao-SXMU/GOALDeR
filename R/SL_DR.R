#' Doubly robust estimator using super learner to fit the GPS model and the outcome model (SL-DR)
#'
#' @param Y Outcome
#' @param A continuous treatment
#' @param X A matrix of covariates
#' @param treat_mod marginal and conditional distributions of treatment
#' @param sl.lib methods combined by super learner and defualt is LASSO, XGBoost, Random Forest, and Support vector
#'
#' @return lm_estimated is the estimated linear dose-response function
#' @export
#'
#' @examples
SL_DR <- function(Y, A, X, treat_mod = c("Normal", "Gamma"),sl.lib = c("SL.glmnet", "SL.randomForest", "SL.xgboost", "SL.svm"))
{
  library(SuperLearner)
  library(locfit)
  library(randomForest)
  library(WeightIt)
  temp.mean <- colMeans(X)
  Temp.mean <- matrix(temp.mean,ncol=dim(X)[2],nrow=nrow(X),byrow=TRUE)
  X <- X - Temp.mean
  temp.sd <- apply(X,FUN=sd,MARGIN=2)
  Temp.sd <- matrix(temp.sd,ncol=dim(X)[2],nrow=nrow(X),byrow=TRUE)
  X <- X / Temp.sd
  rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))

  treat_mod <- match.arg(treat_mod)
  data <- data.frame(Y = Y, A = A, X = X)
  trt_data <- data.frame(A = A, X = X)
  formula_t <- as.formula("A ~ .")
  if (treat_mod == "Normal") {
    result <- SuperLearner(Y = A, X = as.data.frame(X), SL.library = sl.lib)
    gps_fun_Normal <- function(tt) {
      dnorm(tt, mean = result$SL.predict, sd = sd((trt_data$A-result$SL.predict)))
    }
    avg_density_fun <- function(tt) {
      dnorm(tt, mean = mean(result$SL.predict), sd = sd(A))
    }
    gps_fun <- gps_fun_Normal
  } else if (treat_mod == "Gamma") {
    result <- glm(formula = formula_t, family = Gamma(link = "log"),
                  data = trt_data)
    marginal_mod <- glm(formula = A~1, family = Gamma(link = "log"),
                        data = trt_data)

    est_treat <- result$fitted
    shape_gamma <- as.numeric(MASS::gamma.shape(result)[1])
    theta_given_X <- result$fitted.values/shape_gamma
    theta_treat_X <- trt_data$A / shape_gamma


    shape_gamma_marginal   <- as.numeric(MASS::gamma.shape(marginal_mod)[1])
    theta_given_X_marginal <- marginal_mod$fitted.values / shape_gamma_marginal
    theta_treat_X_marginal <- trt_data$A / shape_gamma_marginal

    gps_fun_Gamma <- function(t) {
      dgamma(t, shape = shape_gamma, scale = theta_given_X)
    }
    avg_density_fun <- function(tt) {
      sapply(tt, function(t) mean(dgamma(t, shape = shape_gamma_marginal,
                                         scale = theta_given_X_marginal)) )
    }
    gps_fun <- gps_fun_Gamma
  }


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

  weight_fun <- function(A)
  {
    avg_density_fun(A) / gps_fun(A)
  }

  weights_est <- weight_fun(A)

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

  Huling_DR<-list(
    lmfit=locpoly_fit,
    lm_estimated=estimated,
    pseudo_out=pout
  )

  return(Huling_DR)
}
