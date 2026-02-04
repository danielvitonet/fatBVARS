
# Forecasting natural gas prices in real time

# Define the list of required packages
required_packages <- c(
  "here",
  "R.matlab",
  "dplyr",
  "fatBVARS",
  "pracma",
  "dbarts",
  "stochvol",
  "mfbvar",
  "abind",
  "mvtnorm"
)

# Check each package; if it's not installed, install it
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    renv::install(pkg)
    #renv::snapshot()
    renv::status()
  }
}
# Load the libraries
lapply(required_packages, library, character.only = TRUE)
# Source several helper functions 
set.seed(123) #set a seed to allow reproducibility

#############################################################################
# MF-BAVART SET-UP (ENVIROMENT and DEPENDENCES)                             #                                                                           
# 1) Path alla repo "vendorizzata" dentro ext/                              #
mf_dir <- here("ext", "mf-bavart")                                          #
#
# 2) Crea un environment dedicato per non sporcare il Global Env            #
#mf <- new.env(parent = baseenv())                                          #
mf <- new.env(parent = globalenv())                                         #
#
# 3) Carica i file R nell’environment mf                                    #
# sys.source() è preferibile quando vuoi specificare l’envir esplicitamente #
sys.source(file.path(mf_dir, "aux_func.R"),      envir = mf, chdir = TRUE)  #
sys.source(file.path(mf_dir, "mfbavart_func.R"), envir = mf, chdir = TRUE)  #
#############################################################################

#############################################################################
# flexBART SET-UP (ENVIROMENT and DEPENDENCES)                              #                                                                           
# 1) Path alla repo                                                         #
mix_dir <- here("flexBART")                                                 #
#
# 2) Crea un environment dedicato per non sporcare il Global Env            #
#mix <- new.env(parent = baseenv())                                         #
mix <- new.env(parent = globalenv())                                        #
#
# 3) Carica i file R nell’environment mix                                   #
# sys.source() è preferibile quando vuoi specificare l’envir esplicitamente #
sys.source(file.path(mix_dir, "aux_func.R"), envir = mix, chdir = TRUE)     #
sys.source(file.path(mix_dir, "flexBART.R"), envir = mix, chdir = TRUE)     #
#############################################################################


###################################################################################
# FUNCTIONS NEEDED:                                                               #
# 1) lagn function                                                                #
lagn <- function(data, m) {                                                       #
  # input: data matrix and lag m                                                  #
  data[(m+1):nrow(data), , drop = FALSE] - data[1:(nrow(data)-m), , drop = FALSE] #
}                                                                                 #
#
# 2) function for Quantile Score (QS)                                             #
quantile_levels <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)                       #
QS_sample <-function(true,mcmc,tau=quantile_levels){                              #
  require(pracma)                                                                 #
                                                                                  #
  tau_len <- length(tau)                                                          #
  Q.tau <- stats::quantile(mcmc,probs=tau)                                        #
  true_rep <- rep(true,tau_len)                                                   #
  QS.vec <- (true_rep-Q.tau)*(tau-((true_rep<=Q.tau)*1))                          #
                                                                                  #
  data.frame(                                                                     #
    tau       = tau,                                                              #
    tau_label = names(QS.vec),                                                    #
    quantile_score = unname(QS.vec),                                              #
    row.names = NULL                                                              #
  )                                                                               #
}                                                                                 #
#
# 3) function for quantile-weighted CRPS (qwCRPS)                                 #
qwcrps_tau <- seq(0.01, 0.99, by = 0.01)                                          #
qwcrps_weightings <- c("center", "left", "right", "tails")                        #
qwCRPS_sample <-function(true,mcmc,tau=qwcrps_tau,weighting="none"){              #
  require(pracma)                                                                 #
                                                                                  #
  tau_len <- length(tau)                                                          #
  Q.tau <- stats::quantile(mcmc,probs=tau)                                        #                                       
  true_rep <- rep(true,tau_len)                                                   #                            
  QS.vec <- (true_rep-Q.tau)*(tau-((true_rep<=Q.tau)*1))                          #                                                     
                                                                                  # 
  weights <- switch(tolower(weighting),                                           #                                      
                    "none" = 1,                                                   #                              
                    "tails" = (2*tau-1)^2,                                        #                                         
                    "right" = tau^2,                                              #                                   
                    "left" = (1-tau)^2,                                           #                                      
                    "center" = tau*(1-tau))                                       #                                        
  wghs <- QS.vec*weights                                                          #                       
  return(pracma::trapz(tau,wghs))                                                 ##############                                 
}                                                                                              #                                                                 
compute_qwcrps_scores <- function(true, mcmc, tau=qwcrps_tau, weightings=qwcrps_weightings) {  #
  data.frame(                                                                                  #
    weighting = weightings,                                                                    #           
    qwcrps = vapply(                                                                           #    
      weightings,                                                                              #    
      function(w) qwCRPS_sample(true, mcmc, tau, weighting = w),                               #                                                
      numeric(1)                                                                               #
    )                                                                                          #
  )                                                                                            #
}                                                                                              #
#
# 4) function for Directional Symmetry                                                         #
directional_symmetry <- function(actual, forecast) {                                           #               
  if (length(actual) < 2) {                                                                    #
    return(NA_real_)                                                                           #                                        
  }                                                                                            #                       
  direction_match <- (actual[-1] - actual[-length(actual)]) *                                  #                        
    (forecast[-1] - forecast[-length(forecast)]) > 0                                           #               
  100 / (length(actual) - 1) * sum(direction_match)                                            #              
}                                                                                              #                     
################################################################################################


###################################################################################

# Load real-time datasets
# Rows: T (starting in 1973M1+36)
# Columns: point in real time (1991M1 to 2024M2)

NG_HENRY <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/NG_HENRY.txt", sep="\t", header=FALSE))	         
# nominal gas price, Henry Hub
CPI_AC <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/CPI_AC.txt", sep="\t", header=FALSE))		           
# average change nowcast of US CPI

# Economic predictor variables:
IPALF <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/IPALF.txt", sep="\t", header=FALSE))		             
# industrial production index (ALFRED vintages)
CAPUALF <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/CAPUALF.txt", sep="\t", header=FALSE))		         
# US capacity utilization rate (ALFRED vintages)

PROD_DRY_SA <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/PROD_DRY_SA.txt", sep="\t", header=FALSE))	   
# marketed NG production
STORE_WORK_SA <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/STORE_WORK_SA.txt", sep="\t", header=FALSE)) 
# underground storage of working gas
CONS_TOT <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/CONS_TOT.txt", sep="\t", header=FALSE))	         
# total NG consumption
RIGS <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/RIGS.txt", sep="\t", header=FALSE))		               
# rig count

CDDdev <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/CDDdev.txt", sep="\t", header=FALSE))		           
# cooling degree days in deviation from historical average
HDDdev <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/HDDdev.txt", sep="\t", header=FALSE))		           
# heating degree days in deviation from historical average

# GAS Futures start in April 1990 (1990M4) but next year gas futures start in June 1990 (1990M6)
gas_futures <- as.matrix(read.table("C:/Users/HP/Downloads/DATA - LNG/DATA - models/gas_futures.txt", header=FALSE))
fut <- rbind(matrix(NA, nrow = 207, ncol = 9), gas_futures)

# May 2024 vintage of Henry Hub spot price and CPI (final release data)
# rows: T (1973M1-2024M5)
# columns: most recent vintage of 2024M5
# non-NaN rows start on 1997M1 (row 289)
# C:/Users/HP/Downloads/DATA - LNG/HH_CPI_May2024vintage.mat
HH_CPI <- readMat("C:/Users/HP/Downloads/DATA - LNG/DATA - models/HH_CPI_May2024vintage.mat",
                  fixNames = FALSE,          # conserva nomi MATLAB (underscore)
                  drop = "singletonLists",
                  verbose = TRUE)

CPI_May24 <- as.matrix(HH_CPI$CPI_May24)
NG_May24  <- as.matrix(HH_CPI$NG_May24)


# Basic parameters
#ind <- 288+1   # indicates length of initial real-time sample (up to 1997.1)
ind <- 437+1   # indicates length of initial real-time sample (up to 2009.6)
h <- 1         # forecast horizon
p <- 6         # fixed lag order

#eval_length <- (ncol(CPI_AC) - h) - (12 * 18 + 6) + 1
#sfe <- matrix(NA, nrow = eval_length, ncol = 2)

quantile_scores <- data.frame()
qwcrps_scores <- data.frame()
forecast_history <- data.frame()

#jx = 73   
# 12*6+1 = 73 (6 years and 1 month), is the column representing the 1997M1 data vintage for when we do not use future contracts data.
#jx = 222  
# 222 is the column representing the 2009M6 data vintage for which at least 228 observations for future contracts data is available.
# change 37 with 210 , and change 36 with 209 , when using future contracts data

for (jx in (12*18+6):(ncol(CPI_AC)-h)) {  # adjust for evaluation period 2009M6 to 2024M2
  print(jx)
  cycle_id <- jx - (12 * 18 + 6) + 1
  
  # Create VAR data in column format 
  rpg <- log(100 * NG_HENRY[210:ind, jx] / CPI_AC[210:ind, jx])  # log Real gas price (nominal deflated by US CPI)
  
  # variables in levels 
  capu <- CAPUALF[210:ind, jx]
  hd <- HDDdev[210:ind, jx]
  cd <- CDDdev[210:ind, jx]
  
  # variables in logs
  conslog <- log(CONS_TOT[210:ind, jx])
  
  # variables in growth rates
  ipg <- 100 * lagn(log(IPALF[1:ind, jx, drop = FALSE]), 1)
  dryprod <- 100 * lagn(log(PROD_DRY_SA[1:ind, jx, drop = FALSE]), 1)
  inventories <- 100 * lagn(log(STORE_WORK_SA[1:ind, jx, drop = FALSE]), 1)
  consg <- 100 * lagn(log(CONS_TOT[1:ind, jx, drop = FALSE]), 1)
  rigcount <- 100 * lagn(log(RIGS[1:ind, jx, drop = FALSE]), 1)
  
  
  # one observation lost due to differencing
  ip <- ipg[209:nrow(ipg), 1]
  prod <- dryprod[209:nrow(dryprod), 1]
  store <- inventories[209:nrow(inventories), 1]
  cons <- consg[209:nrow(consg), 1]
  rigs <- rigcount[209:nrow(rigcount), 1]
  
  
  # futures contract variable
  futgasrt <- log(fut[210:ind, ])
  spotgasrt <- log(NG_HENRY[210:ind, jx])
  #inflrt <- log(CPI_AC[(ind - 120 + 2):ind, jx]) - log(CPI_AC[(ind - 120 + 1):(ind - 1), jx]) # average inflation over the past 10 years
  inflrt <- log(CPI_AC[210:ind, jx]) - log(CPI_AC[209:(ind - 1), jx]) # monthly inflation
  
  jj <- switch(as.character(h),
               "1" = 1,
               "3" = 2,
               "6" = 3,
               "9" = 4,
               "12" = 5,
               "15" = 6,
               "18" = 7,
               "21" = 8,
               "24" = 9)
  
  #futs <- futgasrt[, jj] - spotgasrt - ((1 + mean(inflrt))^h - 1)
  infl_ma <- stats::filter(inflrt, rep(1 / 120, 120), sides = 1) # 10-year rolling average
  infl_exp <- (1 + infl_ma)^h - 1
  #futs <- ts(futgasrt[, jj] - spotgasrt - as.numeric(infl_exp), frequency = 12)
  futs <- futgasrt[, jj] - spotgasrt - as.numeric(infl_exp)
  
  # Create revised real price of natural gas (most recent vintage)
  x <- 100 * NG_May24[210:(ind + h), 1] / CPI_May24[210:(ind + h), 1]  # Real gas price (nominal deflated by US CPI)
  
  # Estimate BAVART
  # (Estimation code needed here)
  data <- list(
    log_GAS_Price = ts(rpg, frequency = 12),             # monthly
    mgr_DRY_Production = ts(prod, frequency = 12),       # monthly
    mgr_Working_Inventories = ts(store, frequency = 12), # monthly
    #CONS = ts(cons, frequency = 12),                     # monthly
    mgr_Rig_counts = ts(rigs, frequency = 12),           # monthly
    mgr_INDUSTRIAL_PRODUCTION = ts(ip, frequency = 12),  # monthly
    #CAP_UT = ts(capu, frequency = 12),                   # monthly
    HDD_dev = ts(hd, frequency = 12),                    # monthly
    #CDD_dev = ts(cd, frequency = 12),                    # monthly
    Fut_sp = ts(futs, frequency = 12)                    # monthly
  )
  
  # qgr_Real_GDP = ts(..., frequency = 4)                # quarterly
  
  
  # Estimate BVAR from fatBVARS package BVAR.SV(...)
  bvar_data <- do.call(
    cbind,
    lapply(data, function(series) as.numeric(series))
  )
  bvar_data <- bvar_data[complete.cases(bvar_data), , drop = FALSE]
  
  K <- ncol(bvar_data)
  p <- 1
  h <- h
  
  # Minnesota prior con shrinkage 0.2 e 0.5, MST + SV
  prior <- get_prior(
    y = bvar_data,
    p = p,
    priorStyle = "Minnesota",
    dist = "MST",
    SV = TRUE,
    lambda1 = 0.2,  # overall shrinkage
    lambda2 = 0.5,  # cross-variable shrinkage
    a_Vprior = 10   # a ~ N(0, 10I)
  )
  
  # Inizializzazione
  inits <- get_init(prior, samples = 15000, burnin = 5000)
  
  # ---- opzionale: forzare h0 ~ N(log sigma^2_OLS, 4)
  # (varianza 4 -> sd=2)
  #K <- ncol(bvar_data)
  #t_max <- nrow(y)
  #h0_mean <- log(prior$sigma^2)
  #inits$h <- matrix(
  #  rnorm(K * t_max, mean = rep(h0_mean, t_max), sd = 2),
  #  nrow = K,
  #  ncol = t_max
  #)
  #
  # ---- opzionale: troncare nu tra 4 e 100
  #if (!is.null(inits$nu)) {
  #  inits$nu <- pmin(pmax(inits$nu, 4), 100)
  #}
  
  # Stima MST-SV
  bvar_fit <- BVAR.SV(
    y = y,
    K = K,
    p = p,
    dist = "MST",
    y0 = NULL,
    prior = prior,
    inits = inits
  )
  
  bvar_fcst <- get_forecast(bvar_fit, t_pred = h)
  bvar_fcst_draws <- bvar_fcst$y_pred[, h, "log_GAS_Price"]
  
  
  
  # Estimate mixBART from flexBART(...)
  flex_data <- do.call(
    cbind,
    lapply(data, function(series) as.numeric(series))
  )
  flex_data <- flex_data[complete.cases(flex_data), , drop = FALSE] #Y has T= and M=
  
  # fcst.unconditional 
  mixbart_fit <- flexBART(
    Yraw = flex_data,
    nburn = 15000,
    nsave = 15000,
    thinfac = 1,
    prior = "HS",
    prior.sig = c(3, 0.9),
    pr.mean = matrix(0, ncol(flex_data), ncol(flex_data)),
    model = "mixBART",
    sv = "SV",
    #fc.approx="exact",
    fhorz = h,
    quiet = FALSE
  )
  mixbart_fcst_draws <- mixbart_fit$fcst[, 1, 1]  # This is an array where the first dimension is the number of saved draws, the second is the forecast horizons and the third refer to the different endogenous series
  
  
  # Evaluate h-step ahead forecast
  #sfe[jx - 12 * 6, 2] <- (exp(rpg[length(rpg)]) - x[t + h])^2  # Squared forecast error
  t <- length(rpg)
  actual_level <- x[t + h]
  
  rw_draws <- rep(exp(rpg[length(rpg)]), length(bvar_fcst_draws))
  bvar_draws_level <- exp(bvar_fcst_draws)
  mfbavart_draws_level <- exp(mfbavart_fcst_draws)
  mixbart_draws_level <- exp(mixbart_fcst_draws)
  
  model_draws <- list(
    BVAR = bvar_draws_level,
    MFBAVART = mfbavart_draws_level,
    mixBART = mixbart_draws_level,
    RW = rw_draws
  )
  
  for (model_name in names(model_draws)) {
    draws <- model_draws[[model_name]]
    qs_df <- QS_sample(actual_level, draws)
    qs_df$cycle <- cycle_id
    qs_df$model <- model_name
    quantile_scores <- dplyr::bind_rows(quantile_scores, qs_df)
    
    qw_df <- compute_qwcrps_scores(actual_level, draws)
    qw_df$cycle <- cycle_id
    qw_df$model <- model_name
    qwcrps_scores <- dplyr::bind_rows(qwcrps_scores, qw_df)
    
    point_fcst <- stats::median(draws)
    forecast_history <- dplyr::bind_rows(
      forecast_history,
      data.frame(
        cycle = cycle_id,
        model = model_name,
        actual = actual_level,
        forecast = point_fcst
      )
    )
  }
  
  # Update index for recursive estimation
  ind <- ind + 1
}

# Evaluate real-time recursive forecast accuracy

# Compute Average Quantile Scores
avg_quantile_scores <- quantile_scores |>
  dplyr::group_by(model, tau) |>
  dplyr::summarise(mean_score = mean(quantile_score, na.rm = TRUE), .groups = "drop")
rw_quantile_scores <- avg_quantile_scores |>
  dplyr::filter(model == "RW") |>
  dplyr::rename(rw_score = mean_score) |>
  dplyr::select(tau, rw_score)
avg_quantile_scores_rel <- avg_quantile_scores |>
  dplyr::left_join(rw_quantile_scores, by = "tau") |>
  dplyr::mutate(relative_to_rw = mean_score / rw_score)

# Compute quantile-weighted Continuous Ranked Probability Scores (qw-CRPS)
avg_qwcrps_scores <- qwcrps_scores |>
  dplyr::group_by(model, weighting) |>
  dplyr::summarise(mean_score = mean(qwcrps, na.rm = TRUE), .groups = "drop")
rw_qwcrps_scores <- avg_qwcrps_scores |>
  dplyr::filter(model == "RW") |>
  dplyr::rename(rw_score = mean_score) |>
  dplyr::select(weighting, rw_score)
avg_qwcrps_scores_rel <- avg_qwcrps_scores |>
  dplyr::left_join(rw_qwcrps_scores, by = "weighting") |>
  dplyr::mutate(relative_to_rw = mean_score / rw_score)

# Compute Directional Symmetry (DS) statistic
ds_scores <- forecast_history |>
  dplyr::group_by(model) |>
  dplyr::summarise(
    ds = directional_symmetry(actual, forecast),
    .groups = "drop"
  )
rw_ds <- ds_scores |>
  dplyr::filter(model == "RW") |>
  dplyr::pull(ds)
ds_scores_rel <- ds_scores |>
  dplyr::mutate(relative_to_rw = ds / rw_ds)