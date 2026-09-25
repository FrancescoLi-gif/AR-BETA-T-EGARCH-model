library(tseries)
library(stats)
library(e1071)
library(TSA)
library(gasmodel)
library(xts)
library(rugarch)
library(betategarch)


EXUSUK <- read.csv("C:/Users/Francesco/Desktop/uni_bo/1 anno/Advanced Time Series/Coursework/EXUSUK.csv")
source("C:/Users/Francesco/Desktop/uni_bo/1 anno/Advanced Time Series/Coursework/functions_ats_f.R")


# View(EXUSUK)

# quante sterline valgono un dollaro --> 1£ = x * $

df <- EXUSUK
dim(df)
t <- as.Date(df$observation_date)
y_plot <- xts(df$EXUSUK, order.by = t)
plot(y_plot, t)


# Considering only the complete years (1971 - 2025)
df <- df[1:660,]
dim(df)

#### Analysis of the stationarity of the series #### 
y <- ts(df[,2])
acf(y)
pacf(y)
lag.plot(y, 12)

# Clearly the series is not stationary. Apply some transformation to work with 
# it using the theory we discussed in class, as it is valid for (weakly) stationary
# time series sequences.

y_log <- log(y)
ts.plot(y_log)

# we apply the log transformation to "somehow standardize the data" (approfondire)
# The correlation does not really change and we clearly we do not have stationarity.

y_diff <- diff(y)
ts.plot(y_diff)

### !!!!!!!!!!!!!!!!!!!!! why this is necessary --> research

y_t <- diff(log(df$EXUSUK)) * 100 
y_t <- as.numeric(na.omit(y_t)) 
acf(y_t)
pacf(y_t)
lag.plot(y_t, 12)
ts.plot(y_t)

# By applying the lagged differences at lag 1 with order 1 on the log(y),
# we obtain a clear stationary sequence ready to be studied.


#### Analysis of the distribution of the residuals of the stationary series ####

# Comparison Q-Q plot of the empirical quantiles with those of a Normal distribution
qqnorm(y_t, main = "Normal Q-Q Plot of y_diff")
qqline(y_t, col = "red", lwd = 2) 

# clear heavy tailed distribution, since the extreme quantiles do not coincide
# we suppose we are dealing with a t-student distribution with low df.


# Assume (for pedagogical purposes) we want to model the time varying mean,
# despite the acf does not reveal a clear signal on the mean.

# ------------------------------------------------------------------------------
#-------------------------- KALMAN FILTER --------------------------------------
# ------------------------------------------------------------------------------


# Initialization: phi, sigma_e, sigma_eta
init <- c(0.1, 12, 4.5) # phi , sigma_e, sigma_eta

# Estimation through MLE
kf_est <- estimator_KF(y_t, init)

print("Optimized parameters obtained using Kalman Filter:")
kf_est

# 4. Renaming pf the parameters --> easier to read 
print(kf_est$theta)

phi_opt       <- kf_est$theta[1]
sigma_e_opt   <- kf_est$theta[2]
sigma_eta_opt <- kf_est$theta[3]


# 5. Using Kalman filter with optimized parameters
kf_filter <- KF(y_t, phi_opt, sigma_e_opt, sigma_eta_opt)
mu_kf_ts <- ts(kf_filter$mu_pred, start = start(y_t), frequency = 12)


y_ts <- ts(y_t, start = c(1971, 2), frequency = 12)
mu_kf_ts <- ts(kf_filter$mu_pred, start = c(1971, 2), frequency = 12)


# ------------------------------------------------------------------------------
# ---------------- SCORE DRIVEN (on the MEAN) with Student t errors ------------
# ------------------------------------------------------------------------------

## Inizialization of the static parameters
omega <- 0 # zero intercept
phi <- 0.8
k <- 0.8
varsigma <- 1
nu <- 4 # low degrees of freedom 
theta <- c(omega, phi, k, varsigma, nu) # vector of static parameters

# ESTIMATION
est  <- uDCS_t_model_estimator(y_t, theta)
est$theta_list

ts.plot(y_t)
filter <- uDCS_t_model_filter(y_t, est$theta)

################ GRAPHS ###############

# Alligning all the time series and mds
# we start from february, beacuse we lose january from the lagged sequence

y_ts <- ts(y_t, start = c(1971, 2), frequency = 12)
u_ts <- ts(filter$Innovation_u_t, start = c(1971, 2), frequency = 12)
v_ts <- ts(filter$Innovation_v_t, start = c(1971, 2), frequency = 12)
mu_ts <- ts(filter$Dynamic_Location, start = c(1971, 2), frequency = 12)

par(mfrow = c(2, 2))

# PLOT 1: Time Series VS. Scaled Error (u_t)
ts.plot(y_ts, ylab = "Log-Returns", col = "gray", 
        main = "Time Series vs Scaled Error (u_t)")
lines(u_ts, col = "green", lwd = 1.5)

# PLOT 2: Time Series VS. Raw Innovation (v_t)
ts.plot(y_ts, ylab = "Log-Returns", col = "gray", 
        main = "Time Series vs Raw Innovation (v_t)")
lines(v_ts, col = "blue", lwd = 1.5)

# PLOT 3: Time Series VS. Trend (mu_t)
ts.plot(y_ts, ylab = "Log-Returns", col = "gray", 
        main = "Time Series vs Dynamic Location (mu_t)")
lines(mu_ts, col = "red", lwd = 2)

# PLOT 4:  v_t (blue) VS. u_t (green)
ts.plot(v_ts, ylab = "Innovations", col = "blue", 
        main = "Raw Error (v_t) vs Scaled Error (u_t)")
lines(u_ts, col = "green", lwd = 1.5)

# maybe add legends and change the colors
par(mfrow = c(1, 1))

# --------------------------------------------------------
# COMPARISON BETWEEN KALMAN FILTER AND SCORE DRIVEN MODEL:
# --------------------------------------------------------

ts.plot(y_ts, ylab = "Log-Returns", col = "gray", 
        main = "Student-t DCS VS. Kalman Filter")

mu_gas_ts <- ts(filter$Dynamic_Location, start = c(1971, 2), frequency = 12)

lines(mu_kf_ts, col = "blue", lwd = 2) # KF
lines(mu_gas_ts, col = "red", lwd = 1) # score driven better


# Model comparison with AIC and BIC
# The least the best

kf_est$AIC_kf
kf_est$BIC_kf

# Here we already have - loglikelihood as output for optimizer
(2*est$optimizer$objective) + (2*length(theta)) # AIC
(2*est$optimizer$objective) + (log(length(y))*length(theta)) # BIC

# Both of the two criteria suggest using the student t error model

# ------------------------------------------------------------------------------


# Il trend non cattura l'andamento della serie stazionaria
# Inserire styliseed fact
# Studiare andamento varianza

# ------------------------------------------------------------------------------
# -------------------- TIME VARYING VARIANCE -----------------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# ---------------------------  GARCH(1,1) --------------------------------------
# ------------------------------------------------------------------------------

spec_norm <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
  mean.model = list(armaOrder = c(0, 0)), 
  distribution.model = "norm"
)

fit_norm <- ugarchfit(spec = spec_norm, data = y_t)
fit_norm@fit$matcoef

# Since the AIC and BIC are divided by "n", I multiply their default values
# in the package to make them comparable also to the location models.

# AIC for AR(1)+GARCH(1,1)
AIC_norm <- infocriteria(fit_norm)[1] * fit_norm@model$modeldata$T
AIC_norm

# BIC for AR(1)+GARCH(1,1)
BIC_norm  <- infocriteria(fit_norm)[2] * fit_norm@model$modeldata$T
BIC_norm

# ------------------------------------------------------------------------------
# Qui di sotto trovi i grafici che ho riprodotto usando lo stesso stile di quelli
# fatti per il location scale.
# Ho semplicemente rifatto i grafici inclusi nel pacchetto con lo stile grafico
# che abbiamo usato noi in precedenza.
# Lascio scelta a voi su quali vi piacciono di più.


#plot(fit_norm)

# GRAFICO 9: QQPLOT OF STANDARDIZED RESIDUALS
res_norm <- residuals(fit_norm, standardize = T)
qqnorm(res_norm, main = "Normal Q-Q Plot of y_diff")
qqline(res_norm, col = "red", lwd = 2) 


# GRAFICO 8: HISTOGRAMS OF STANDARDIZED RESIDUALS VS NORMAL DENSITY
hist(res_norm, breaks = 20, probability = TRUE, 
     main = "Histogram of Standardized Residuals VS Normal Density",
     xlab = "Standardized Residuals",
     ylim = c(0, 0.45))
curve(dnorm(x, mean = 0, sd = 1), add = TRUE, col = "red", lwd = 2)
# Solita interpretazione per entrambi (si potrebbe considerare di eliminare uno dei due)
# Il fit sembra buono, forse unico problemino è dove vediamo valori negativi
# ma still ottimo per residui


# GRAFICO 10: ACF OF RESIDUALS
acf(as.numeric(res_norm), main = "ACF of residuals")
# IMPORTANTE QUESTO!!!!!!
# ALTISSIMA CORRELAZIONE AL PRIMO LAG. QUESTO GIUSTIFICA LA COMPONENTE AR(1)
# CHE INSERIAMO NEI SUCCESSIVI MODELLI!!
# Per valori a lag più alti è OK.

# GRAFICO 11:
acf(as.numeric(res_norm^2), main = "ACF of squared residuals")
# Irrilevante quasi a ogni lag.
# Significa che per quanto riguarda la modellazione della varianza potremmo
# essere quasi a posto anche con questo modello, non servono ulteriori sforzi
# particolari per migliorarne il fit.


# FORMULA LATEX MODELLO:

# $$y_t = \mu + \epsilon_t$$
# $$\epsilon_t = \sigma_t z_t \quad \text{with} \quad z_t \sim \mathcal{N}(0, 1)$$
# $$\sigma_t^2 = \omega + \alpha_1 \epsilon_{t-1}^2 + \beta_1 \sigma_{t-1}^2$$




# ------------------------------------------------------------------------------
# ----------------------  AR(1) + GARCH(1,1) -----------------------------------
# ------------------------------------------------------------------------------

spec_norm_ar <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
  mean.model = list(armaOrder = c(1, 0)),
  distribution.model = "norm")

fit_norm_ar <- ugarchfit(spec = spec_norm_ar, data = y_t)
fit_norm_ar@fit$matcoef

# AIC for AR(1)+GARCH(1,1)
AIC_norm_ar <- infocriteria(fit_norm_ar)[1] * fit_norm_ar@model$modeldata$T
AIC_norm_ar

# BIC for AR(1)+GARCH(1,1)
BIC_norm_ar  <- infocriteria(fit_norm_ar)[2] * fit_norm_ar@model$modeldata$T
BIC_norm_ar


# ------------------------------------------------------------------------------
# Come prima, solo cha alcune interpretazioni le lascio a voi.

#plot(fit_norm_ar)

# GRAFICO 9: QQPLOT OF STANDARDIZED RESIDUALS
res_norm_ar <- residuals(fit_norm_ar, standardize = T)
qqnorm(res_norm_ar, main = "Normal Q-Q Plot of y_diff")
qqline(res_norm_ar, col = "red", lwd = 2) 


# GRAFICO 8: HISTOGRAMS OF STANDARDIZED RESIDUALS VS NORMAL DENSITY
hist(res_norm_ar, breaks = 20, probability = TRUE, 
     main = "Histogram of Standardized Residuals VS Normal Density",
     xlab = "Standardized Residuals",
     ylim = c(0, 0.45))
curve(dnorm(x, mean = 0, sd = 1), add = TRUE, col = "red", lwd = 2)


# GRAFICO 10: ACF OF RESIDUALS
acf(as.numeric(res_norm_ar), main = "ACF of residuals")
# IMPORTANTE QUESTO!!!!!!
# QUI CHIARAMENTE INDICA CHE UNA COMPONENTE AR(1) è SUFFICIENTE PER MODELLARE
# LA MEDIA!! IL MOTIVO PER CUI LO MANTENIAMO ANCHE DOPO

# GRAFICO 11:
acf(as.numeric(res_norm_ar^2), main = "ACF of squared residuals")

# INTERPRETATION:
# Big improvement in model selection by including a simple AR(1) component

# FORMULA LATEX MODELLO:

# $$y_t = \mu + \phi_1 y_{t-1} + \epsilon_t$$
# $$\epsilon_t = \sigma_t z_t \quad \text{with} \quad z_t \sim \mathcal{N}(0, 1)$$
# $$\sigma_t^2 = \omega + \alpha_1 \epsilon_{t-1}^2 + \beta_1 \sigma_{t-1}^2$$



# ------------------------------------------------------------------------------
# --------------------------- BETA-T-EGARCH + AR(1) ----------------------------
# ------------------------------------------------------------------------------

# Qui sotto trovi le formule latex del modello fittato:
# Purtroppo non c'era questo modello specifico nel pacchetto, quindi mi sono
# autinflitto del dolore e l'ho fatto a mano. Guarda function_ts per dettagli

#$$y_t = \omega_m + \phi_m y_{t-1} + \epsilon_t$$
#  $$\epsilon_t = \exp(\lambda_t) z_t \quad \text{with} \quad z_t \sim t_{\nu}(0, 1)$$
#  $$\lambda_t = \omega_v + \phi_v \lambda_{t-1} + \kappa_v u_{t-1}$$
#  $$u_{t-1} = \frac{(\nu + 1)\epsilon_{t-1}^2}{\nu \exp(2\lambda_{t-1}) + \epsilon_{t-1}^2} - 1$$


# Inizializziamo il vettore di parametri
# Il result dovrebbe essere consistente anche se cambiamo initialization
theta_init_beta_t_egarch <- c(omega_m = 0, 
                              phi_m   = 0.05, 
                              omega_v = 0.01, 
                              phi_v   = 0.90, 
                              k_v     = 0.05, 
                              nu      = 5)

# MLE estimation con loglikelihood del modello:
# ATTENTION: it may take a while to run
fit_beta_t_egarch_ar <- estimator_AR1_BetaEGARCH(y_ts, theta_init_beta_t_egarch)
filter_beta_t_ar <- Joint_AR1_BetaEGARCH_filter(y_ts, fit_beta_t_egarch_ar$theta)

# OUTPUT MODELLO
fit_beta_t_egarch_ar$results_table

# AIC and BIC
fit_beta_t_egarch_ar$AIC
fit_beta_t_egarch_ar$BIC
# INTERPRETATION:
# This model is clearly the best so far, yielding the lowest AIC and BIC.


# MANCANO LE DIAGNOSTICHE SUI RESIDUI!!!
# SIMILI A COME FATTE PRIMA


# ------------------------------------------------------------------------------
# ---------- PLOTS: COMPARSION BETWEEN TIME VARYING VARIANCE MODELS
# ------------------------------------------------------------------------------

t_start <- start(y_ts)
t_freq  <- frequency(y_ts)

# Save the Dynamic Mean
mu_norm    <- ts(as.numeric(fitted(fit_norm)), start = t_start, frequency = t_freq)
mu_norm_ar <- ts(as.numeric(fitted(fit_norm_ar)), start = t_start, frequency = t_freq)
mu_beta_t  <- filter_beta_t_ar$Dynamic_Mean

# Save Conditional Variance (sigma^2)
# Note: rugarch's sigma() extracts the standard deviation, so we square it.
var_norm    <- ts(as.numeric(sigma(fit_norm)^2), start = t_start, frequency = t_freq)
var_norm_ar <- ts(as.numeric(sigma(fit_norm_ar)^2), start = t_start, frequency = t_freq)
var_beta_t  <- filter_beta_t_ar$Conditional_Var

# Save Innovations 
e_norm    <- ts(as.numeric(residuals(fit_norm, standardize = FALSE)), start = t_start, frequency = t_freq)
e_norm_ar <- ts(as.numeric(residuals(fit_norm_ar, standardize = FALSE)), start = t_start, frequency = t_freq)
e_beta_t  <- filter_beta_t_ar$Innovation_e_t

# ------------------------------------------------------------------------------
# PLOT 1: DYNAMIC MEAN

ts.plot(y_ts, col = "lightgray", 
        main = "Conditional Mean Comparison")

lines(mu_norm,    col = "black", lwd = 1, lty = 2) # Constant mean
lines(mu_norm_ar, col = "blue",  lwd = 1, lty = 2) # Gaussian AR(1)
lines(mu_beta_t,  col = "red",   lwd = 1, lty = 2) # Beta-t-EGARCH AR(1)

legend("topleft", 
       legend = c("y_t", "GARCH", "AR(1)+GARCH", "AR(1)+Beta-t-EGARCH"),
       col = c("lightgray", "black", "blue", "red"))
# DA AGGIUSTARE LA LEGENDA

# ------------------------------------------------------------------------------
# PLOT 2: CONDITIONAL VARIANCE
# IMPORTANTE PER COMPARAZIONE!!!

ts.plot(y_ts^2, ylab = "Variance", col = "lightgray", 
        main = "Conditional Variance Comparison")

lines(var_norm,    col = "black", lwd = 1, lty = 2) 
lines(var_norm_ar, col = "blue",  lwd = 1, lty = 2)
lines(var_beta_t,  col = "red",   lwd = 1, lty = 2)

legend("topleft", 
       legend = c("y_t^2", "GARCH", "AR(1)+GARCH", "AR(1)+Beta-t-EGARCH"),
       col = c("lightgray", "black", "blue", "red"))
# Da aggiustare legenda

# Ci mostra come in ordine, il piu robusto agli shock del mercato (outliers)
# è il beta t egarch model + AR(1).
# Questo è seguito da garch + AR(1).
# Il peggiore (troppo sensibile agli shock del mercato, che interpreta le 
# variazioni storiche come "segnale") è il GARCH


# ------------------------------------------------------------------------------
# PLOT 3: INNOVATIONS 
# Background: Raw Returns (y_t). 
# This shows how much of the original return is "explained" by the mean equation.
ts.plot(y_ts, ylab = "Innovations", col = "black", 
        main = "Innovations (Residuals)", lwd = 2.25)

lines(e_norm,    col = "darkred", lwd = 1.75, lty = 2)
lines(e_norm_ar, col = "blue",  lwd = 1, lty = 2)
lines(e_beta_t,  col = "gold",   lwd = 1, lty = 2)

legend("topleft", 
       legend = c("y_t", "e_t GARCH", "e_t AR(1)+GARCH", "e_t AR(1)+Beta-t-EGARCH"),
       col = c("lightgray", "black", "blue", "red"))
# da aggiustare

# Questo si potrebbe anche cambiare, si sovrappongono i le innovazioni








###############         STYLSED FACTS          ############################

# NOTA PER NOI: QUESTI SONO DA INSERIRE NELLA TRANSIZIONE TRA LOCATION E SCALE 
# TIME VARYING MODELS (I.E. TRA SCORE DRIVEN E GARCH)


# 1.  Raw monthly returns show no obvious pattern (zero autocorrelation),
#     making markets look like a pure random walk.
plot(y_plot, t)
acf(y_ts)

# 2.  Uncorrelated does not mean independent. If they were truly independent, 
#     no mathematical transformation would reveal a pattern.
#     On the squares we see a clear auto-correlation.
#     Log returns are therefore not random walk.
acf(y_ts^2)


# 3.  Volatility Clustering: The data proves that big market moves cluster 
#     together, and quiet periods cluster together. 
#     Because this structure exists, prices do not follow a strict random walk.
#     After big growths or decays, respectively a decay or a shock are plausible

ts.plot(y_ts)

# 4. Heavy tails: if we have a high kurtosis we have a heavy tailed distribution
kurtosis(y_ts) # this is the kurtosis excess --> 1.826586

# we have a moderately heavy tailed distribution 


####################    PROVA  ########################
# PLOT 3: INNOVATIONS 
# Background: Raw Returns (y_t). 
# This shows how much of the original return is "explained" by the mean equation.

# 1. Background (y_t): Grigio chiaro per fare da vero "sfondo" senza rubare la scena
ts.plot(y_ts, ylab = "Innovations", col = "gray75", 
        main = "Innovations (Residuals)", lwd = 3.5)

# 2. Innovazioni: Colori ad alto contrasto, spessori decrescenti e linee continue
# Sulle serie storiche molto fitte, troppe linee tratteggiate (lty = 2) creano "rumore visivo".
lines(e_norm,    col = "darkred",    lwd = 2.5,  lty = 1)
lines(e_norm_ar, col = "royalblue",  lwd = 1.75, lty = 1)
lines(e_beta_t,  col = "darkorange", lwd = 1,    lty = 1) # Sostituito il "gold" (troppo chiaro)

# 3. Legenda: Sincronizzata con i colori esatti e i NUOVI SPESSORI
legend("topleft", 
       legend = c("y_t (Raw Returns)", "e_t GARCH", "e_t AR(1)+GARCH", "e_t AR(1)+Beta-t-EGARCH"),
       col = c("gray75", "darkred", "royalblue", "darkorange"),
       lty = c(1, 1, 1, 1),
       lwd = c(3.5, 2.5, 1.75, 1), # <-- Spessori allineati esattamente alle linee tracciate
       bty = "n",      # Rimuove il bordo nero attorno alla legenda per maggiore pulizia
       cex = 0.85,     # Riduce leggermente il testo per non coprire i dati
       bg = "white")   # Sfondo bianco nel caso i dati ci passino sotto