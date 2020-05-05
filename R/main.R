#' Conduct two-sample Mendelian Randomization using summary statistics
#'
#' Runs the main MR-Robin algorithm: a two-sample Mendelian Randomization method
#' ROBust to correlated and some INvalid instruments.
#'
#' @param eqtl_betas matrix of coefficient estimates (betas) from eQTL study.
#' @param eqtl_se matrix of standard errors for coefficient estimates from eQTL study.
#' @param gwas_betas vector of coefficient estimates (betas) from GWAS study.
#' @param gwas_se vector of standard errors for coefficient estimates from GWAS study.
#' @param LD matrix of LD correlation coefficients (\eqn{r}, not \eqn{r^2}).
#' @param snpID vector of variant identifiers to be used as instrumental variables.
#'
#' @return An object of class \code{lmerMod}, returned from the reverse regression random slope
#' mixed model run by MR-Robin.
#'
#' To conduct inference on the returned results, use function \code{\link{MR_Robin_resample}}.
#'
#' @details The following are additional details describing the input arguments.
#' For \code{eqtl_betas} and \code{eqtl_se}, each row \eqn{i} corresponds to a SNP/variant
#' while each column \eqn{j} holds the summary statistics in condition \eqn{j} (e.g. a tissue type).
#' Note that the matrix \code{LD} should hold \emph{correlation coefficients}
#' (i.e. \eqn{r}), not their squared values (\eqn{r^2}).
#'
#' @export
#'
MR_Robin <- function(eqtl_betas,eqtl_se,gwas_betas,gwas_se,LD,snpID){

  ## determine number of tissues (or conditions/studies/etc.)
  nT <- ncol(eqtl_betas)

  ## set up coefficients for reverse regression
  beta_x <- matrix(eqtl_betas,ncol=1)
  beta_y <- rep(gwas_betas, nT)

  ## standard errors (for weights)
  se_x <- matrix(eqtl_se, ncol=1)
  se_y <- matrix(gwas_se, ncol=1) ## not used by function; used in resampling (consider returning with results)

  ##identifiers
  snpID <- rep(snpID,nT)

  ## return results from weighted regression with random slopes (and correlated errors)
  return(lme4::lmer(beta_x~(beta_y-1)+(beta_y-1|snpID),weights=1/se_x^2))
}


#' Obtain p-value from MR-Robin results
#'
#' Uses a resampling procedure to estimate a \eqn{P}-value for a MR-Robin object.
#'
#' @param MR_Robin_res an object of class \code{lmerMod}, returned by \code{MR_Robin}.
#' @param gwas_se vector of standard errors for coefficient estimates from GWAS study.
#' @param nsamp integer of the number of samples to use in estimating \eqn{P}-value
#' using a null distribution.
#' @param LD matrix of LD correlation coefficients (\eqn{r}, not \eqn{r^2}).
#'
#' @return A list of two elements:
#' \tabular{ll}{
#' \code{pvalue} \tab numeric of the estimated \eqn{P}-value.\cr
#' \code{nsamp_used} \tab integer of the number of samples used in estimating the \eqn{P}-value.\cr
#' }
#'
#' \code{nsamp_used} is returned because not all \code{nsamp} may be used (samples will be dropped
#' if the model does not converge or results in a singular fit of the random slope).
#'
#' @export
#'
MR_Robin_resample <- function(MR_Robin_res,gwas_se,nsamp=1000,LD){

  ## extract data from MR-Robin results
  eqtl_betas <- MR_Robin_res@frame$beta_x
  snpID <- MR_Robin_res@frame$snpID
  weights <- MR_Robin_res@frame$`(weights)`
  tstat_MR_Robin <- summary(MR_Robin_res)$coefficients[1,3]
  nT <- length(eqtl_betas)/length(unique(snpID))

  ## bootstrapped null distribution, accounting for LD correlations
  beta_gwas_nullMat <- mvtnorm::rmvnorm(nsamp,mean=rep(0,length(gwas_se)),sigma=diag(gwas_se) %*% LD %*% diag(gwas_se))

  ## initialize return vectors
  tstat_nulls <- NULL
  nsamp_used <- 0

  ## run MR-Robin on the null models
  for(i in 1:nsamp){
    beta_gwas_null <- rep(beta_gwas_nullMat[i,], nT)
    lme_null <- lmer(eqtl_betas~(beta_gwas_null-1)+(beta_gwas_null-1|snpID),weights=weights)
    if(is.null(summary(lme_null)$optinfo$conv$lme4$messages)){
      tstat_nulls <- c(tstat_nulls, summary(lme_null)$coefficients[1,3])
      nsamp_used <- nsamp_used + 1
    }
  }

  if(nsamp_used==0) stop("All resampled datasets failed to converge in random slope model.")

  if(nsamp_used<=(nsamp/2)) warning("More than half of resampled datasets failed to converge in random slope model.")

  pval <- mean(abs(tstat_nulls) >= abs(tstat_MR_Robin))

  return(list(pvalue=pval,nsamp_used=nsamp_used))
}
