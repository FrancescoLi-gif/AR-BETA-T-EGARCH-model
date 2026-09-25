##### Simulator ####

uDCS_t_model_simulator <- function(T, omega, phi, k, varsigma, nu){
  
  ###Define the Processes
  y   <- array(data = NA, dim = c(T) ) 
  
  ###Define Dynamic Location and Innovations
  mu_t <- array(data = NA, dim = c(T))
  u_t  <- array(data = NA, dim = c(T-1))
  
  ###Initial value for the recursion 
  mu_t[1]   <- omega
  
  ###Generate the first observations of the process
  y[1]   <- uSTDT_rnd(1, mu_t[1], varsigma, nu)
  
  ###Dynamics 
  for (t in 2:T) {
    
    ###Factor Innovations
    u_t[t-1] <- martingale_diff_u_t(y[t-1], mu_t[t-1], varsigma, nu)
    
    ###Updating Filters                    
    mu_t[t]   <- omega + phi * (mu_t[t-1] - omega) + k * u_t[t-1]
    
    ###Generate the observations of the processes
    y[t] <- uSTDT_rnd(1, mu_t[t], varsigma, nu)
  }
  
  ####### OUTPUT #######
  
  
  ###Make List
  out <- list(y_t_gen          = as.ts(y),
              Dynamic_Location = as.ts(mu_t),
              Innovation_u_t   = as.ts(u_t))
  
  return(out)
}


################# Univariate Student's t Random Generator ##################


uSTDT_rnd <- function(n, mu, varsigma, nu) {
  
  z <- rt(n, df = nu) 
  y <- numeric()
  for(i in 1:n){
    y[i] <- c(mu + z[i]* sqrt(varsigma) ) 
  }
  
  return(y)
}



################
##### Interprete
################
interprete_uDCS_t_model <- function(dati, param){
  
  omega    <- param[1]
  phi      <- param[2]
  k        <- param[3]
  varsigma <- param[4] # Direct scale input
  nu       <- param[5]
  
  theta_new <- c(omega, phi, k, varsigma, nu)
  
  fitness <- uDCS_t_model_filter(dati, theta_new)$Log_Likelihood
  
  if(is.na(fitness) | !is.finite(fitness)) fitness <- -1e10
  if(fitness != fitness) fitness <- -1e10
  
  return(-fitness)
}

################
##### Estimator
################
uDCS_t_model_estimator <- function(dati, param){
  
  Start <- Sys.time()
  T_obs <- length(dati)
  
  omega    <- param[1]
  phi      <- param[2]
  k        <- param[3]
  varsigma <- param[4] 
  nu       <- param[5]
  
  theta_st <- c(omega, phi, k, varsigma, nu)
  
  # Bounds: varsigma > 0 (strictly positive scale), nu > 2 (finite variance)
  lower <- c(-Inf, -0.999, -2, 1e-05, 2.001)
  upper <- c( Inf,  0.999,  2, Inf, 300)
  
  optimizer <- suppressWarnings(nlminb(start = theta_st, objective = interprete_uDCS_t_model, 
                                       dati  = dati, gradient = NULL, 
                                       control = list(trace = 0), hessian = NULL,
                                       lower = lower, upper = upper))
  
  omega_opt    <- optimizer$par[1]  
  phi_opt      <- optimizer$par[2]
  k_opt        <- optimizer$par[3]
  varsigma_opt <- optimizer$par[4]
  nu_opt       <- optimizer$par[5]
  
  theta_opt <- c(omega_opt, phi_opt, k_opt, varsigma_opt, nu_opt)
  names(theta_opt) <- c("omega", "phi", "k", "varsigma", "nu")
  
  theta_list <- list(omega = omega_opt,
                     phi = phi_opt,
                     k  = k_opt,
                     varsigma = varsigma_opt,
                     nu    = nu_opt)
  
  # ----------------------------------------------------------------------------
  # Asymptotic Covariance Matrix Evaluation (Harvey & Luati, 2014)
  # ----------------------------------------------------------------------------
  
  a_const <- phi_opt - k_opt * (nu_opt / (nu_opt + 3))
  b_const <- phi_opt^2 - 2 * phi_opt * k_opt * (nu_opt / (nu_opt + 3)) + 
    (k_opt^2 * nu_opt * (nu_opt^3 + 10 * nu_opt^2 + 35 * nu_opt + 38)) / 
    ((nu_opt + 1) * (nu_opt + 3) * (nu_opt + 5) * (nu_opt + 7))
  
  if(b_const >= 1) {
    warning("Condition b < 1 violated. Filter is non-stationary.")
    std_errors <- rep(NA, 5)
    z_stats    <- rep(NA, 5)
    p_values   <- rep(NA, 5)
  } else {
    
    sigma2_u <- varsigma_opt * (nu_opt^2) / ((nu_opt + 3) * (nu_opt + 1))
    
    D_psi <- matrix(0, nrow = 3, ncol = 3)
    
    # Elements as in paper eq.(14), with 1/(1-b) already included throughout
    D_psi[1, 1] <- sigma2_u / (1 - b_const)
    D_psi[1, 2] <- (sigma2_u * k_opt * a_const) /
      ((1 - a_const * phi_opt) * (1 - b_const))
    D_psi[2, 1] <- D_psi[1, 2]
    D_psi[2, 2] <- (sigma2_u * k_opt^2 * (1 + a_const * phi_opt)) /
      ((1 - phi_opt^2) * (1 - a_const * phi_opt) * (1 - b_const))
    D_psi[3, 3] <- ((1 - phi_opt)^2 * (1 + a_const)) /
      ((1 - a_const)^2 * (1 - b_const))
    
    # Scale to information matrix: I_psi = D_psi / varsigma
    I_psi <- D_psi / varsigma_opt
    
    I_loc <- matrix(0, nrow = 3, ncol = 3)
    I_loc[1, 1] <- I_psi[3, 3] 
    I_loc[2, 2] <- I_psi[2, 2] 
    I_loc[3, 3] <- I_psi[1, 1] 
    I_loc[2, 3] <- I_psi[2, 1] 
    I_loc[3, 2] <- I_psi[1, 2] 
    
    h_nu <- 0.5 * trigamma(nu_opt / 2) - 0.5 * trigamma((nu_opt + 1) / 2) - 
      (nu_opt + 5) / (nu_opt * (nu_opt + 3) * (nu_opt + 1))
    
    I_lambda_nu <- matrix(0, nrow = 2, ncol = 2)
    I_lambda_nu[1, 1] <- (2 * nu_opt) / (nu_opt + 3)
    I_lambda_nu[1, 2] <- -2 / ((nu_opt + 3) * (nu_opt + 1))
    I_lambda_nu[2, 1] <- I_lambda_nu[1, 2]
    I_lambda_nu[2, 2] <- 0.5 * h_nu
    
    # Transform from lambda-space to varsigma-space (In the paper we have a different parametrization)
    J_scale <- diag(c(1 / (2 * varsigma_opt), 1))
    I_varsigma_nu <- t(J_scale) %*% I_lambda_nu %*% J_scale
    
    I_full <- matrix(0, nrow = 5, ncol = 5)
    I_full[1:3, 1:3] <- I_loc
    I_full[4:5, 4:5] <- I_varsigma_nu
    
    vcv_matrix <- solve(I_full) / T_obs
    
    std_errors <- sqrt(diag(vcv_matrix))
    z_stats    <- theta_opt / std_errors
    p_values   <- 2 * (1 - pnorm(abs(z_stats)))
  }
  
  results_table <- data.frame(
    Estimate  = theta_opt,
    Std_Error = std_errors,
    Z_value   = z_stats,
    p_value   = p_values
  )
  
  if(!is.na(b_const) && b_const < 1) {
    results_table$Signif <- cut(results_table$p_value, 
                                breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, 1), 
                                labels = c("***", "**", "*", ".", " "))
  }
  
  Elapsed_Time <- Sys.time() - Start
  print(paste("Elapsed Time: ", toString(Elapsed_Time)))
  
  out <- list(results_table = results_table,
              theta_list    = theta_list,
              theta         = theta_opt,
              optimizer     = optimizer)
  
  return(out) 
}

##### Filtering #####

uDCS_t_model_filter <- function(y, theta){
  
  ###Take T
  T <- length(y)
  
  ###Define LogLikelihoods
  dloglik <- array(data = NA, dim = c(T))
  loglik  <- numeric()
  
  ###Parameter Selections Dynamic Location
  omega <- theta[1]
  phi   <- theta[2]
  k     <- theta[3]
  
  varsigma <- theta[4]
  nu       <- theta[5]
  
  ###Define Dynamic Location and Innovations
  mu_t <- array(data = NA, dim = c(T+1))
  u_t  <- array(data = NA, dim = c(T))
  v_t <-  array(data = NA, dim = c(T))
  
  ###Initialize Dynamic Location
  mu_t[1]   <- (omega)
  
  ###Initialize Likelihood
  dloglik[1] <- uSTDT_uDCS_t(y[1], mu_t[1], varsigma = varsigma, nu = nu, log = TRUE)
  loglik     <- dloglik[1]
  
  for(t in 2:(T+1)) {
    v_t[t-1] <- martingale_diff_u_t(y[t-1], mu_t[t-1], varsigma, Inf)
    ###Dynamic Location Innovations
    u_t[t-1] <- martingale_diff_u_t(y[t-1], mu_t[t-1], varsigma, nu)
    ###Updating Filter                    
    mu_t[t]   <- omega + phi * (mu_t[t-1] - omega) + k * u_t[t-1]
    
    if(t < (T+1)){
      ###Updating Likelihoods
      dloglik[t] <- uSTDT_uDCS_t(y[t], mu_t = mu_t[t], varsigma = varsigma, nu = nu, log = TRUE)
      loglik     <- loglik + dloglik[t]
    }
  }
  
  ######################
  ####### OUTPUT #######
  ######################
  mu_t <- ts(mu_t, start = start(y), frequency = frequency(y))
  u_t  <- ts(u_t, start = start(y), frequency = frequency(y))
  v_t  <- ts(v_t, start = start(y), frequency = frequency(y))
  
  ###Make List
  out <- list(Dynamic_Location = mu_t,
              Innovation_u_t   = u_t,
              Innovation_v_t   = v_t,
              Log_Densities_i  = dloglik,
              Log_Likelihood   = loglik)
  
  return(out)
}


####### ADDITIONAL FUNCTIONS #######



martingale_diff_u_t <- function(y, mu_t, varsigma, nu){
  
  u_t <- c((1 / (1 + (y - mu_t)^2/(nu*varsigma)) * (y - mu_t)))
  
  return(u_t)
}


uSTDT_uDCS_t <- function(y, mu_t, varsigma, nu, log = TRUE){
  
  ulpdf <- (lgamma((nu + 1) / 2) - lgamma(nu / 2) - (1/2) * log(varsigma) -
              (1/2)  * log(pi * nu) - ((nu + 1) / 2) * log(1 + (y - mu_t)^2 / (nu*varsigma) ))
  
  if(log != TRUE){
    ulpdf <- exp(ulpdf)
  } 
  
  return(ulpdf)
}
############# KALMAN FILTER ################

# we shall create our own functions to simulate and estimate the parameters as 
# well as recovering the predicted state by the KF for an AR(1) 
# signal plus noise model with our own created function 
# 
### warning: this file will be quite rough, you will then make it both personal and elegant 
################################################################################

# dgp stands for data generating process, i.e we are simulating an AR(1) plus noise 
# process: the input is made of the length of the series and the static parameters 

dgp <-function(n, phi, sigma_e, sigma_eta){
  
  eta = sigma_eta * rnorm(n)
  e = sigma_e * rnorm(n)
  
  mu = 0
  mu[1] = 0  # unconditional mean of an AR(1) process with no intercept
  y = 0 
  for(t in 1:(n-1)){
    mu[t+1] = phi * mu[t] + eta[t]
    y[t+1] =  mu[t+1] + e[t+1]
  }
  
  # ts.plot(y)
  # lines(mu,col='red')
  
  out <- list(y = as.ts(y), mu = as.ts(mu))
  
  return(out)
}


################################################################################
# KALMAN filter recursions with our own function
################################################################################

KF <-function(y, phi, sigma_e, sigma_eta){
  
  # we here separate the arguments of the function that will be estimated - phi, sigma_e,
  # sigma_eta by the arguments of the function that will not - mu_1|0 and P_1|0 
  
  m10 = 0
  P10 = 1
  
  n = NROW(y)
  
  # allocate space 
  mu_pred = array(data = NA, dim = c(n))  # this is mu_{t|t-1}
  P = 0                                   # this is P_{t|t-1} 
  v = array(data = NA, dim = c(n))        # this  will be the innovation error 
  K = 0                                   # the Kalman gain 
  F = 0                                   # the conditional variance of v_t 
  dllk = array(data = NA, dim = c(n))    #the log-likelihood value 
  llk = 0 
  
  # initialise the recursion 
  
  mu_pred[1] = m10
  P[1] = P10
  
  # the recursion 
  
  for(t in 1:(n-1)){
    v[t] = y[t] - mu_pred[t]; 
    F[t] = P[t] + sigma_e^2;
    K[t] = (phi * P[t])/F[t]
    P[t+1] = phi^2 * P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
    mu_pred[t+1] = phi * mu_pred[t] + K[t]*v[t]
    dllk[t] = (-0.5 * log(2 * pi)) + (-0.5 * (log(F[t]) + v[t]^2/F[t]))
    llk  = llk + dllk[t]
  }
  v[n] = y[n] - mu_pred[n]; 
  F[n] = P[n] + sigma_e^2;
  K[n] = (phi * P[n])/F[n]
  dllk[n] = (-0.5 * log(2 * pi)) + (-0.5 * (log(F[n]) + v[n]^2/F[n]))
  llk  = llk + dllk[n]
  
  #llk = sum(dllk)
  
  out <- list(mu_pred = as.ts(mu_pred), llk = llk)
  #out <- list(v, mu_pred)
  
  return(out)
}


################################################################################
# KALMAN filter likelihood 
################################################################################

loglikelihood <- function(par,y){
  
  # par should be a vector  
  
  phi       <- par[1]
  sigma_e   <- par[2]
  sigma_eta <- par[3]
  
  theta_new <- c(phi, sigma_e,sigma_eta)
  
  obj = KF(y, par[1],par[2],par[3])$llk
  
  return(-obj) 
}


################################################################################
# KALMAN filter parameter estimation
################################################################################
# input: the data and the initial values of the parameters that we aim to estimate
# objective function: the KF Gaussian likelihood given from the KF 
# (prediction error decomposition), see the function above 
# output: the estimated parameters
# the filtered mu will be recovered from the KF recursions 

estimator_KF <- function(y,par){
  
  n = NROW(y)
  
  phi       <- par[1]
  sigma_e   <- par[2]
  sigma_eta <- par[3]
  
  theta_0 <- c(phi, sigma_e, sigma_eta)
  
  # Definiamo i limiti
  lim_inf <- c(-0.999, 1e-6, 1e-6) # [phi, sigma_e, sigma_eta]
  lim_sup <- c( 0.999,  Inf,  Inf) # [phi, sigma_e, sigma_eta]
  
  # UNICA chiamata a nlminb corretta. 
  # L'obiettivo è "loglikelihood" e passiamo "y = y"
  #   optimizer <- suppressWarnings(optim(start = theta_0, 
  #                                        objective = loglikelihood, 
  #                                        y = y,
  #                                        hessian = TRUE,
  #                                        method = "L-BFGS-B",
  #                                        lower = lim_inf, 
  #                                        upper = lim_sup))
  #   
  #   # Estraiamo i parametri ottimizzati
  #   phi_opt       <- optimizer$par[1]
  #   sigma_e_opt   <- optimizer$par[2]
  #   sigma_eta_opt <- optimizer$par[3]
  #   
  #   theta_opt <- c(phi_opt, sigma_e_opt, sigma_eta_opt)
  #   
  #   ### Creiamo la lista con i parametri
  #   theta_list <- list(phi = phi_opt,
  #                      sigma_e = sigma_e_opt,
  #                      sigma_eta = sigma_eta_opt)
  #   
  #   out <- list(theta_list = theta_list,
  #               theta = theta_opt,
  #               optimizer = optimizer)
  #   
  #   return(out) 
  # }
  
  # Ottimizzazione
  optimizer <- suppressWarnings(optim(par = theta_0,            # Corretto: 'par' invece di 'start'
                                      fn = loglikelihood,       # Corretto: 'fn' invece di 'objective'
                                      y = y,
                                      hessian = TRUE,
                                      method = "L-BFGS-B",
                                      lower = lim_inf, 
                                      upper = lim_sup))
  
  # 1. Estraiamo i parametri ottimizzati
  theta_opt <- optimizer$par
  names(theta_opt) <- c("phi", "sigma_e", "sigma_eta")
  
  # 2. Calcolo degli Errori Standard
  # Invertiamo la matrice Hessiana per ottenere la matrice di varianza-covarianza
  # Usiamo tryCatch per evitare che la funzione si blocchi se l'Hessiana non è invertibile
  vcv_matrix <- tryCatch(solve(optimizer$hessian), 
                         error = function(e) {
                           warning("Hessiana non invertibile.")
                           return(matrix(NA, nrow = length(theta_opt), ncol = length(theta_opt)))
                         })
  
  # Gli standard error sono la radice quadrata della diagonale principale
  std_errors <- sqrt(diag(vcv_matrix))
  
  # 3. Z-Test (Statistica Z e P-value)
  z_stats <- theta_opt / std_errors
  p_values <- 2 * (1 - pnorm(abs(z_stats))) # P-value a due code
  
  # 4. Creazione di una tabella riassuntiva per l'output
  results_table <- data.frame(
    Estimate = theta_opt,
    Std_Error = std_errors,
    Z_value = z_stats,
    p_value = p_values
  )
  
  # Aggiungiamo le "stelline" di significatività per comodità visiva
  results_table$Signif <- cut(results_table$p_value, 
                              breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, 1), 
                              labels = c("***", "**", "*", ".", " "))
  
  # AIC-BIC
  ll_kf <- -optimizer$value
  AIC_kf <- (-2 * ll_kf) + (2 * length(theta_opt))
  BIC_kf <- (-2 * ll_kf) + (log(n) * length(theta_opt))
  
  # Prepariamo la lista di output
  out <- list(
    results_table = results_table,
    theta = theta_opt,
    AIC_kf = AIC_kf,
    BIC_kf = BIC_kf
    #vcv_matrix = vcv_matrix,
    #optimizer = optimizer
  )
  
  return(out) 
}


# ------------------------------------------------------------------------------
# ----------------- BETA T-EGARCH + AR(1) MODEL --------------------------------
# ------------------------------------------------------------------------------

Joint_AR1_BetaEGARCH_filter <- function(y, theta) {
  
  T_obs <- length(y)
  
  # 1. Extract Parameters explicitly
  omega_m <- theta[1]
  phi_m   <- theta[2]
  
  omega_v <- theta[3]
  phi_v   <- theta[4]
  k_v     <- theta[5]
  
  nu      <- theta[6]
  
  # 2. Initialize State Arrays
  lambda_t  <- rep(NA, T_obs + 1) 
  epsilon_t <- rep(NA, T_obs)
  u_v_t     <- rep(NA, T_obs)
  dloglik   <- rep(NA, T_obs)
  
  # 3. Initial Condition for Log-Volatility
  lambda_t[1] <- omega_v / (1 - phi_v)
  
  # Pre-compute constant part of the Student-t log-likelihood
  ll_const <- lgamma((nu + 1) / 2) - lgamma(nu / 2) - (0.5 * log(pi * nu))
  
  # 4. Joint Filtering Recursion
  for (t in 1:T_obs) {
    
    # --- A. Explicit AR(1) Innovation (Structural Shock) ---
    if (t == 1) {
      uncond_mean <- omega_m / (1 - phi_m)
      epsilon_t[t] <- y[t] - uncond_mean
    } else {
      epsilon_t[t] <- y[t] - omega_m - phi_m * y[t-1]
    }
    
    sigma2_t <- exp(2 * lambda_t[t])
    
    # --- B. Log-Likelihood Evaluation ---
    ll_scale  <- -lambda_t[t]
    ll_kernel <- -((nu + 1) / 2) * log(1 + (epsilon_t[t]^2) / (nu * sigma2_t))
    
    dloglik[t] <- ll_const + ll_scale + ll_kernel
    
    # --- C. Variance Score Computation ---
    u_v_t[t] <- ((nu + 1) * epsilon_t[t]^2) / (nu * sigma2_t + epsilon_t[t]^2) - 1
    
    # --- D. Variance State Update ---
    lambda_t[t+1] <- omega_v + phi_v * lambda_t[t] + k_v * u_v_t[t]
  }
  
  Total_LogLikelihood <- sum(dloglik)
  
  lambda_filtered   <- lambda_t[1:T_obs]
  cond_var_filtered <- exp(2 * lambda_filtered) * (nu / (nu - 2))
  mu_filtered       <- y - epsilon_t
  
  out <- list(
    Dynamic_Mean      = ts(mu_filtered, start = start(y), frequency = frequency(y)),
    Dynamic_LogScale  = ts(lambda_filtered, start = start(y), frequency = frequency(y)),
    Conditional_Var   = ts(cond_var_filtered, start = start(y), frequency = frequency(y)),
    Innovation_e_t    = ts(epsilon_t, start = start(y), frequency = frequency(y)),
    Log_Likelihood    = Total_LogLikelihood
  )
  
  return(out)
}

# CORRECTED: param must be the first argument for optim() compatibility
interprete_joint_model <- function(param, dati) {
  
  fitness <- Joint_AR1_BetaEGARCH_filter(dati, param)$Log_Likelihood
  
  if(is.na(fitness) | !is.finite(fitness)) fitness <- -1e10
  
  return(-fitness) 
}

estimator_AR1_BetaEGARCH <- function(dati, param_init) {
  
  Start <- Sys.time()
  T_obs <- length(dati)
  
  lim_inf <- c(-Inf, -0.999, -Inf, -0.999, -2, 2.001)
  lim_sup <- c( Inf,  0.999,  Inf,  0.999,  2, 300)
  
  optimizer <- suppressWarnings(optim(par = param_init, 
                                      fn = interprete_joint_model, 
                                      dati = dati,
                                      hessian = TRUE,
                                      method = "L-BFGS-B",
                                      lower = lim_inf, 
                                      upper = lim_sup))
  
  theta_opt <- optimizer$par
  names(theta_opt) <- c("omega_m", "phi_m", "omega_v", "phi_v", "k_v", "nu")
  
  # ----------------------------------------------------------------------------
  # NUMERICAL STANDARD ERRORS 
  # ----------------------------------------------------------------------------
  
  vcv_num <- tryCatch(solve(optimizer$hessian), 
                      error = function(e) {
                        warning("Numerical Hessian not invertible.")
                        return(matrix(NA, nrow = 6, ncol = 6))
                      })
  
  se_num <- sqrt(diag(vcv_num))
  z_num  <- theta_opt / se_num
  p_num  <- 2 * (1 - pnorm(abs(z_num)))
  
  
  results_table <- data.frame(
    Estimate = theta_opt,
    SE   = se_num,
    Z_val    = z_num,
    p_val    = p_num
  )
  
  # CORRECTED: Targeting p_Ana since p_value does not exist
  results_table$Signif <- cut(results_table$p_val, 
                              breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, 1), 
                              labels = c("***", "**", "*", ".", " "))
  
  ll_opt  <- -optimizer$value
  AIC_val <- (-2 * ll_opt) + (2 * length(theta_opt))
  BIC_val <- (-2 * ll_opt) + (log(T_obs) * length(theta_opt))
  
  Elapsed_Time <- Sys.time() - Start
  cat("Elapsed Time: ", toString(Elapsed_Time), "\n")
  
  out <- list(
    results_table = results_table,
    theta         = theta_opt,
    AIC           = AIC_val,
    BIC           = BIC_val
  )
  
  return(out) 
}

