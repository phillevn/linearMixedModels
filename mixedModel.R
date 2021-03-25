
mixedModel <- function(y, ID, X,Z, niter=1000){
  beta_coeff = solve(t(X)%*%X)%*%t(X)%*%(y) # initial value of beta coefficient
  # beta_coeff <- rep(rnorm(ncol(X)))
  D_0 <- diag(ncol(Z))
  sigm_2 = matrix(1,nrow = 1)
  error <- 1
  iter_count <- 1
  n_obsever = nrow(y)/length(ID)
  y_minus_X_beta <- matrix(y - X%*%beta_coeff)
  while ((error > 10^{-3}) & (iter_count < niter)) {
    iter_count = iter_count + 1
    D_tmp = matrix(rep(0,ncol(Z)*ncol(Z)), ncol = ncol(Z))
    sigm_tmp = 0
    mu <- NULL
    u <- NULL
    c = 0
    trace = 0
    
    D_0_inv = solve(D_0)
    for (i in unique(ID)){
      n_i_id = which(ID==i)
      y_i = y[n_i_id,]
      X_i = X[n_i_id,]
      Z_i = Z[n_i_id,]
      ZZ_i = t(Z_i)%*%Z_i
      
      Ome_i = solve(ZZ_i/sigm_2[1,1] + D_0_inv)
      
      y_minus_X_beta_i = y_minus_X_beta[n_i_id]
      mu_i = Ome_i%*%t(Z_i)%*%(y_minus_X_beta_i)/sigm_2[1,1]
      mu = c(mu, mu_i)
      mu_ii = Ome_i + mu_i%*%t(mu_i)
      u = c(u,Z_i%*%mu_i)
      D_tmp = D_tmp + mu_ii
     sigm_tmp = sigm_tmp + sum((y_i - X_i%*%beta_coeff)^2) + sum(diag(t(Z_i)%*%Z_i%*%mu_ii)) -2*t(y_i - X_i%*%beta_coeff)%*%Z_i%*%mu_i
    }
    
    beta_coeff_1 = solve(t(X)%*%X)%*%t(X)%*%(y - u)
    beta_error = abs(beta_coeff_1 - beta_coeff)
    beta_coeff = beta_coeff_1
    y_minus_X_beta <- matrix(y - X%*%beta_coeff)
    sigm_2_1 <- sigm_tmp/length(y)
    sigm_2_error = abs(sigm_2_1 - sigm_2)
    sigm_2 = sigm_2_1
    D_1 = D_tmp/length(unique(ID))
    D_error = abs(D_1 - D_0)
    D_0 = D_1
    error = sum(beta_error) + sum(sigm_2_error) + sum(D_error)
  }
  return(list(beta=beta_coeff, sigma_square = sigm_2, Covariance = D_0))
}
