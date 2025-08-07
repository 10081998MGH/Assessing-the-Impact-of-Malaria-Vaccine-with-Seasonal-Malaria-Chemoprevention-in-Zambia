
# Assessing the Impact of Malaria Vaccine with Seasonal Malaria Chemoprevention in Zambia

# Candidate no. 1090601
# University of Oxford 
# Submitted in partial fulfilment of the requirements for the degree of Master of Science in Modelling for Global Health 
# August 2025


# Malaria Deterministic Compartmental Age-structured Model.
# Age-structured human population (0-5 months, 6-59 months, >60 months).
# Multiple interventions: Vaccination,insecticide-treated nets (ITN) data from CSV, Seasonal malaria chemoprevention (SMC). 
# Seasonal transmission patterns using rainfall data.
# Scenario analysis and sensitivity analysis.

# Clear the R environment to start fresh
rm(list=ls())


# Load required packages for the analysis

library(pacman)
p_load(deSolve,               # For solving differential equations(ODEs)
       tidyverse,             # For data manipulation and visualization
       readxl,                # For reading Excel files
       dplyr,                 # For data manipulation (part of tidyverse but explicit load)
       sf,                    # For spatial data handling (shapefiles)
       terra)                 # For raster data processing (rainfall data)

library(lhs)                  # For Latin Hypercube Sampling (LHS) for sensitivity analysis
library(triangle)             # For triangular distributions in sensitivity analysis 
library(plotly)               # For interactive visualizations
library(EnvStats)             # For environmental statistics 
library(patchwork)            # For combining ggplot2 plots
library(gridExtra)            # For arranging multiple ggplot2 plots
library(grid)                 # For grid graphics in R
library("scales")             # For scaling functions in ggplot2

getwd()
# change to your respective directory (Remember to remove this before submission as it is not need )
setwd("C:/mgh-mmid/Placement Project")
getwd()

#Read in insecticide-treated nets (ITN) data from CSV
# Note 
# 2000 - 2012 Model warm-up period for stabilization. Back-predicted data from this period were not used for analysis.
# 2013-2023 Actual net data for Zambia, sourced from the World Malaria Report.
# 2024-2035 Predicted net data, projected based on historical trends and realistic budget constraints.



itndat<-read.csv("netsdata1.csv")
itn_cov <-itndat$nets # Proportion of the population covered by nets

# Read in the observed malaria incidence(per 1000) data from Excel
observed_data <- read_excel("Malaria_incidence data.xlsx")


# Load Zambia shapefile for spatial rainfall extraction
# This defines the geographic boundaries for rainfall data extraction
Zambia_shapefile <- read_sf("C:/mgh-mmid/Placement Project/Rainfall Package/Shapefiles/zmb_admbnda_adm0_dmmu_20201124.shp")

# Get list of all rainfall raster files (CHIRPS data)
# CHIRPS = Climate Hazards Group InfraRed Precipitation with Station data
zarain_layers <- list.files(path = "C:/mgh-mmid/Placement Project/Rainfall Package/", 
                            pattern = ".tif", full.names = TRUE)

# Load all rainfall rasters into a single raster stack
zarain_raster <- rast(zarain_layers)

# Extract district names from shapefile for later joining
district_names <- Zambia_shapefile %>%  as.data.frame  %>%  dplyr::select(ADM0_EN)

# Extract mean rainfall values for each district and time period
# This calculates average rainfall across each district boundary
extract_rain <- terra::extract(zarain_raster, Zambia_shapefile, fun = mean) %>% 
  bind_cols(district_names) 

# Convert rainfall data from wide to long format for easier manipulation
# This transforms columns like "chirps-v2.0.2013.01" into separate rows
extract_rain_longer <- extract_rain %>% 
  pivot_longer(cols=starts_with("chirps-v2.0."),names_to = "variable",values_to = "value")

extract_rain_sf <- extract_rain_longer %>% left_join(Zambia_shapefile %>% dplyr::select(ADM0_EN, geometry),by = "ADM0_EN") %>% 
  st_as_sf()

extract_rain_sf <- extract_rain_longer %>%
  mutate(
    year = str_extract(variable, "\\d{4}"),
    month = str_extract(variable, "(?<=\\d{4}\\.)\\d{2}")
  ) %>%
  left_join(Zambia_shapefile %>% dplyr::select(ADM0_EN, geometry), by = "ADM0_EN") %>%
  st_as_sf() %>%
  mutate(
    # Create date from year and month
    ym_date = ymd(paste0(year, "-", month, "-01")),
    
    # Get number of days in each month
    days_in_month = days_in_month(ym_date)
  ) %>%
  arrange(ym_date) %>%
  mutate(
    # Cumulative offset starting from 0
    day_offset = lag(cumsum(days_in_month), default = 0)
  )

# Filter rainfall data for specific year (2013 used as reference)
rain_2025 <- extract_rain_sf %>%
  filter(year == "2013") 

# Normalize rainfall data using min-max scaling
# This creates values between 0 and 1 for seasonal forcing
Normalised_extract_rain_sf <- extract_rain_sf %>% 
  group_by(ADM0_EN) %>%  # Group by district
  mutate(
    # Min-max normalization
    value_normalized_mm = (value - min(value, na.rm = TRUE)) / 
      (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))
  ) %>%
  ungroup()

####################################################
# Define dynamic Human-static Vector model ####
####################################################
# Definitions ####
# S: Susceptible 
# E: Exposed 
# A: Asymptomatic infections
# C: Clinical infections that will never be treated 
# Ct: Clinical infections destined to be treated 
# Sev: Severe untreated infections
# T: Clinical infections under treatment 
# H: Severe infections under treatment in HOSPITAL
# R: Removed (temporarily immune)
# V: Vaccinated (if vaccination is implemented)
#SMC: Seasonal Malaria Chemoprevention 

# Define the main ODE model function
# This implements a compartmental model with age structure and interventions

Marks_go <- function(t, state, parms)  {
  with(as.list(c(state, parms)), {
    
    
    # Calculate total population in each age group
    # Age group 0: 0-5 months 
    P0 <-  S0 + E0 + A0 + C0 + Ct0 + Sev0 + T0 + H0 + R0
    
    # Age group 1: 6-59 months (young children, main target for interventions)
    P1 <-  S1 + E1 + A1 + C1 + Ct1 + Sev1 + T1 + H1 + R1 + E1V + V + A1V + C1V + Ct1V + Sev1V + T1V + H1V + R1V + Smcs + Smcr + VSmcs + VSmcr 
    
    # Age group 2: 60+ months (older children and adults)
    P2 <-  S2 + E2 + A2 + C2 + Ct2 + Sev2 + T2 + H2 + R2 + F2
    
    P  <- P0 + P1 + P2 # Total population
    
    
    # Control vaccination start time
    if (t < timevax) {
      cov1 = 0
    }
    
    # Calculate overall probability of receiving treatment
    # This is the product of: seeking care × getting tested × test sensitivity × receiving treatment
    ptrt <-  pseek * ptest * psens * ptreat
    
    # Extract normalized rainfall values for seasonal forcing 
    rain_cov <- Normalised_extract_rain_sf$value_normalized_mm 
    
    rain_t <-  Normalised_extract_rain_sf$day_offset 
    
    names(rain_t) <- NULL# time in days
    rain_t
    
    # Interpolate rainfall value for current time point
    # Uses constant interpolation (step function) with extrapolation
    rain<- approx(rain_t, rain_cov, max(0,t), method = "constant", rule = 2)$y
    
    # Calculate seasonal transmission multiplier
    # Combines rainfall with cosine seasonal pattern
    seas <- 1 + rain * amp* cos(2 * pi * (t / 365 - phi))^peak
    
    
    # SEASONAL MALARIA CHEMOPREVENTION (SMC) IMPLEMENTATION
    
    # SMC campaign parameters 
    rounds <- 9               # Number of SMC rounds
    sdur <- 30                # Duration of each SMC round (days)
    sprotect <- 90            # Duration of protection after SMC (days)          
    
    # Create SMC pulse pattern over time
    # This creates a time series of when SMC is administered 
    pulse <-  c(rep(0, 9490),
                rep(c(rep(1, sdur), rep(0, (365-sdur))),rounds),
                rep(0, (9490+rounds*365)))
    
    # Initialize SMC rates
    srate<-snap<-0 # Artificail fast rate that is used to simulate the remove people from the SMC compartment 
    
    # Apply SMC during specified time period(2026-2034)
    # Calculate SMC administration rate based on coverage and pulse timing
    if (t > 9490 & t < 12775) {
      srate = approx(1:length(pulse), (-log(1 - scov) / sdur) * pulse, t, method = "constant", rule = 2)$y
      
      # Calculate protection removal rate (snap back to regular compartments)
      if (scov > 0) {
        snapr_use = snapr
        snap = approx(1:length(pulse), (snapr_use * lag(pulse, default = 0, n = sdur + sprotect)), t, method = "constant", rule = 2)$y
      } else {
        snapr_use = 0
        snap = 0
      }
    }
    
    
    # Calculate ITN coverage as proportion of population covered
    itn_cov <-  itndat$nets * ppn / P # nets × persons per net / total population
    
    # Extract time points for ITN distribution
    itn_t <-  itndat$time
    
    # Calculate effective ITN coverage accounting for net attrition over time
    cum_cov1 <- itn_cov                                   # Current year distribution
    cum_cov2 <-  lag(cum_cov1, default = 0,  n = 1)       # Previous year (with attrition)
    cum_cov3 <-  lag(cum_cov1, default = 0, n = 2)        # Two years ago (with more attrition)
    cum_cov <-  att1 * cum_cov1 + att2 * cum_cov2 + att3 * cum_cov3  # Weighted average accounting for net degradation over time
    # View(cbind(cum_cov1, cum_cov2, cum_cov3, cum_cov))
    
    # Calculate effective ITN protection
    itn_effcov <-  itn_use * itn_eff * pmin(cum_cov, 1) # Combines: usage rate × effectiveness × coverage
    
    #  Interpolating ITN coverage for current time point
    itn <- approx(itn_t, itn_effcov, t, method = "constant", rule = 2)$y
    
    
    
    # Calculate total infectious individuals contributing to transmission
    Infectious <-  C0 + C1 + C2 + Ct0 + Ct1 + Ct2 + Sev0 + Sev1 + Sev2 + F2 + zeta_a * (A0+A1+A2) + zeta_t * (T0+T1+T2)  
    
    Seas <- 1 + amp * cos(2 * pi * (t / 365 - phi))^peak
   
    # Calculate force of infection (transmission rate)
    # Explanation:
    # lambda: Force of infection on humans from mosquitoes.
    # seas * (1 - itn): Scales transmission by seasonality and ITN (insecticide-treated net) usage.
    # a^2 * b * c * m * (Infectious / P): Baseline transmission intensity — 
    #   where:
    #   - a: mosquito biting rate,
    #   - b: probability of transmission from mosquito to human,
    #   - c: probability of transmission from human to mosquito,
    #   - m: mosquito-to-human ratio,
    #   - Infectious / P: proportion of humans currently infectious.
    # Denominator (a * c * Infectious / P + mu_m): Adjusts for mosquito mortality and infection saturation — 
    #   the more infectious humans, the more quickly mosquitoes get infected, but this is offset by death rate mu_m.
    # (gamma_m / (gamma_m + mu_m)): Probability that an infected mosquito survives the extrinsic incubation period 
    # (gamma_m is the rate of incubation, mu_m is mosquito death rate).
    
    lambda <- Seas * (1 - itn) * (a^2 * b * c * m * Infectious / P) / (a * c * Infectious / P + mu_m) * (gamma_m / (gamma_m + mu_m))
    
    
    # Children between 0 to 6 months old 
    
    
    # Age group 0 (0-5 months old) differential equations
    dS0 <- mu_b*P - lambda  * S0 + prho * R0 - mu_h * S0 - kappa * cov1 * S0 - kappa * (1 - cov1) * S0
    dE0 <- lambda * S0 - (gamma_h + mu_h) * E0 - kappa * E0
    dA0 <- pa * gamma_h * E0 + omega * C0 - (delta + mu_h) * A0 -  kappa * cov1 * A0 - kappa * (1 - cov1) * A0
    dC0 <- (1 - pa) * (1 - ptrt) * gamma_h * E0 - (omega + nu + mu_h) * C0 - kappa * C0
    dCt0 <- (1 - pa) * ptrt * gamma_h * E0 - (tau + mu_h) * Ct0 - kappa * Ct0
    dSev0 <- nu * C0 -  mu_h * Sev0 -  tau * mu_sev * Sev0 - tau * (1 - mu_sev) * Sev0 - kappa * Sev0
    dT0 <- tau * Ct0 - (r + mu_h) * T0 - kappa * T0
    dH0 <- tau * (1 - mu_sev) * Sev0 - (r_sev + mu_h) * H0 - kappa * H0
    dR0 <- r * T0 + r_sev * H0 + delta * A0 - (prho + mu_h) * R0 - kappa * cov1 * R0 - kappa * (1 - cov1) * R0
    
    
    # Age group 1 (6-59 months old) differential equations
    dS1 <-kappa * (1 - cov1) * S0 - lambda  * S1 + prho * R1 - mu_h * S1 - kappa1 * S1 - srate * S1 + snap * Smcs
    dE1 <-  kappa * E0 + lambda * S1 - (gamma_h + mu_h) * E1 - kappa1 
    dA1 <-  kappa * (1 - cov1) * A0 + pa * gamma_h * E1 + omega * C1 - (delta + mu_h) * A1 - kappa1 * A1 - srate * A1
    dC1 <-  kappa * C0 + (1 - pa) * (1 - ptrt) * gamma_h * E1 - (omega + nu + mu_h) * C1 - kappa1 * C1
    dCt1 <-  kappa * Ct0 + (1 - pa) * ptrt * gamma_h * E1 - (tau + mu_h) * Ct1 - kappa1 * Ct1
    dSev1 <-  kappa * Sev0 +  nu * C1 -  mu_h * Sev1 -  tau * mu_sev * Sev1 - tau * (1 - mu_sev) * Sev1 - kappa1 * Sev1
    dT1 <-  kappa * T0 + tau * Ct1 - (r + mu_h) * T1 - kappa1 * T1
    dH1 <-  kappa * H0 +  tau * (1 - mu_sev) * Sev1 - (r_sev + mu_h) * H1 - kappa1 * H1
    dR1 <-kappa * (1 - cov1) * R0 + r * T1 + r_sev * H1 + delta * A1 - (prho + mu_h) * R1 - kappa1 * R1 - srate * R1 + snap * Smcr
    
    
    # Seaonal malaria chemoprevention (SMC) compartments for age 6-59 months
    dSmcs <- srate * (S1) - (snap + mu_h) * Smcs              # children from S1 moving into SMC-protected susceptible (Smcs) 
    dSmcr <- srate * (A1 + R1) - (snap + mu_h) * Smcr         # children from A1 or R1 moving into SMC-protected recovered (Smcr)
    
    
    # Vaccinated compartment flows (for age group 6-59 months)
    dV <- (S0 + A0 + R0) * kappa * cov1 - (rho_V + mu_h) * V - lambda * V - srate * V + snap * VSmcs
    
    dE1V <-  lambda * V - (gamma_h + mu_h) * E1V + prho * R1V - rho_V 
    dA1V <-   paV * gamma_h * E1V + omega * C1V - (delta + mu_h) * A1V - rho_V * A1V - srate * A1V
    dC1V <-  (1 - paV) * (1 - ptrt) * gamma_h * E1V - (omega + nu + mu_h) * C1V - rho_V * C1V
    dCt1V <- (1 - paV) * ptrt * gamma_h * E1V - (tau + mu_h) * Ct1V - rho_V * Ct1V
    dSev1V <-  nu * C1V - mu_h * Sev1V - tau * mu_sev * Sev1V -  tau * (1 - mu_sev) * Sev1V - rho_V * Sev1V
    dT1V <-   tau * Ct1V - (r + mu_h) * T1V - rho_V * T1V
    dH1V <-  tau * (1 - mu_sev) * Sev1V - (r_sev + mu_h) * H1V - rho_V * H1V
    dR1V <-  r * T1V + r_sev * H1V + delta * A1V - (prho + mu_h) * R1V - rho_V * R1V - srate * R1V + snap * VSmcr
    
    # SMC compartments for vaccinated group (6-59 months, vaccinated)
    dVSmcs <-  srate * V - (snap + mu_h) * VSmcs
    dVSmcr <-  srate * (A1V + R1V) - (snap + mu_h) * VSmcr
    
    # Children above 60 months old 
    dS2 <-  kappa1 * S1 - lambda  * S2 + rho * R2 - mu_h * S2 +  rho_V * V 
    dE2 <-  kappa1 * E1 + lambda * S2 - (gamma_h + mu_h) * E2 + rho_V * E1V 
    dA2 <-  kappa1 * A1 + pa * gamma_h * E2 + gamma_h * pas * F2 + omega * C2 - (delta + mu_h) * A2 + rho_V * A1V
    dC2 <-  kappa1 * C1+ (1 - pa) * (1 - ptrt) * gamma_h * E2 + (1 - pas) * (1 - ptrt) * gamma_h *F2- (omega + nu + mu_h) * C2 + rho_V * C1V
    dCt2 <-  kappa1 * Ct1 + (1 - pa) * ptrt * gamma_h * E2 + (1 - pas) * ptrt * gamma_h * F2 - (tau + mu_h) * Ct2 + rho_V * Ct1V
    dSev2 <-  kappa1 * Sev1 + nu * C2 -  mu_h * Sev2 -  tau * mu_sev * Sev2 - tau * (1 - mu_sev) * Sev2 + rho_V * Sev1V
    dT2 <-  kappa1 * T1 + tau * Ct2 - (r + mu_h) * T2 + rho_V * T1V
    dH2 <- kappa1 * H1 + tau * (1 - mu_sev) * Sev2 - (r_sev + mu_h) * H2 + rho_V * H1V
    dR2 <-  kappa1 * R1 + r * T2 + r_sev * H2 + delta * A2 - (lambda + rho + mu_h) * R2 + rho_V * R1V
    dF2 <-  lambda * R2 -  (gamma_h + mu_h) * F2 
    
    
    # Counters
    
    CInc0 <-  lambda * S0 
    CInc1 <-  lambda * (S1 + V)
    CInc2 <-  lambda * (S2 + R2) 
    CIncU5 <- CInc0 + CInc1
    CInc  <-  CInc0 + CInc1 + CInc2 
    
    
    
    dSMC <-  srate*(S1+A1+R1+V+A1V+R1V)
    dCSnap <-  snap*(S1+R1+V+R1V)
    
    # Total number of cases.
    Cases0 <-   gamma_h * (1 - pa)  * E0
    Cases1 <-   gamma_h * (1 - pa)  * E1 + gamma_h * (1 - paV)  * E1V
    Cases2 <-  gamma_h * (1 - pa)  * E2 + gamma_h * (1 - pas)  * F2
    CasesU5 <-  Cases0 + Cases1 
    Cases  <- Cases0 + Cases1 + Cases2 # Total cases
    
    
    Csevere0  <-  nu * C0
    Csevere1  <-  nu * C1 + nu * C1V 
    Csevere2  <-  nu * C2
    CsevereU5 <- Csevere0 + Csevere1 
    Csevere   <-  Csevere0 + Csevere1 + Csevere2
    
    
    # Total number of Deaths 
    
    deaths0 <- tau * mu_sev * Sev0
    deaths1 <- tau * mu_sev * Sev1 + tau * mu_sev * Sev1V
    deaths2 <- tau * mu_sev * Sev2
    deathsU5 <- deaths0 + deaths1
    Tdeaths <- deaths0 + deaths1 + deaths2
    
    
    
    
    
    
    
    output <- c(dS0, dE0, dA0, dC0, dCt0, dSev0, dT0, dH0, dR0, dV,
                dS1, dE1, dA1, dC1, dCt1, dSev1, dT1, dH1, dR1,dSmcs,dSmcr,
                dE1V,dA1V, dC1V, dCt1V, dSev1V, dT1V, dH1V, dR1V,dVSmcs,dVSmcr,
                dS2, dE2, dA2, dC2, dCt2, dSev2, dT2, dH2, dR2,
                dF2,CInc0,CInc1,CInc2, CIncU5, CInc, dSMC, dCSnap,
                Cases0, Cases1, Cases2, CasesU5, Cases, Csevere0, Csevere1, Csevere2, CsevereU5,Csevere,
                deaths0, deaths1, deaths2, deathsU5, Tdeaths 
    )
    list(output,lambda= lambda, seas= seas,itn=itn)#, lambda= lambda, itn=itn, irs=irs)
  })
}

# MODEL PARAMETERS AND INITIAL CONDITIONS

# 2000
# AgeGroup0 = 408,458
# 
# AgeGroup1 = 1,628,833
# 
# AgeGroup2 = 9,058,710
# 
# Total Population = 11,096,001

# Initial state vector - all compartments start with realistic values
state <- c(S0 = 211495, E0 = 6348, A0 = 17047, C0 = 1333 , Ct0 =444, Sev0 = 0, T0 = 0, H0 = 0, R0 =170791, V=0,
           S1 = 845979, E1 = 25392, A1 = 68189, C1 = 5332, Ct1 =1776, Sev1 = 0, T1 = 0, H1 = 0, R1 = 683164,Smcs=0,Smcr=0,
           E1V = 0, A1V = 0, C1V = 0, Ct1V = 0, Sev1V = 0, T1V = 0, H1V = 0, R1V = 0,VSmcs=0,VSmcr=0,
           S2 = 4652177, E2 =139636 , A2 = 374983, C2 =29321 , Ct2 =9764, Sev2 = 0, T2 = 0, H2 = 0, R2 =3756829, 
           F2 = 0, CInc0 = 0,CInc1 = 0,CInc2 = 0,  CIncU5=0, CInc = 0, SMC=0, CSnap = 0, cases0 = 0, cases1 = 0, cases2 = 0, casesU5 = 0, Cases = 0,
           Csevere0 = 0, Csevere1 = 0, Csevere2 = 0, CsevereU5 = 0,Csevere = 0, deaths0 = 0, deaths1 = 0, deaths2 = 0, deathsU5 = 0, Tdeaths = 0
)

# Parameters
parms <- c(
  # Transmission parameters
  a = 0.33,                                   # human feeding rate per mosquito
  b = 0.098,                                  # transmission efficiency M->H
  c = 0.44,                                   # transmission efficiency H->M
  
  # Mosquito parameters
  mu_m = 1/15,                                # natural birth/death rate in mosquitoes
  gamma_m = 1/10,                             # rate of onset of infectiousness in mosquitoes
  m =0.98,                                    # mosquito-human ratio
  
  # Human demographic parameters
  mu_b =44.8/(1000*365),                      # Crude birth rate in humans
  mu_h = 11/(1000*365),                       # crude death rate in humans
  
  
  
  #  parameters  
  gamma_h = 1/14,                             # rate of onset of infectiousness in humans
  tau = 1/3,                                  # rate of treatment seeking 
  tau_sev = 1/3,                              # rate of treatment seeking to Hospital (severe infections)
  r = 1/6,                                    # treatment recovery rate 
  r_sev = 1/10,                               # treatment recovery rate 
  rho = 1/(1*365),                            # rate of loss of immunity
  prho = 1/(1*185),
  
  
  pa = 0.44,                                  # proportion of new infections from susceptible population who are asymptomatic
  pas = 0.93,                                 # likelihood of asymptomatic infection after secondary infection
  delta = 1/285,                              # natural recovery rate of asymptomatic infections
  nu = 1/7,                                   # rate of onset of severe symptoms in untreated clinical infections
  
  omega = 1/20,                               # rate of loss of symptoms in untreated clinical infections
  mu_sev = 0.00086,                           # rate of death from severe untreated illness
  zeta_a = 0.22,                              # relative infectiousness of asymptomatic infections
  zeta_t = 0.04,                              # relative infectiousness of treated infections
  
  # Treatment-seeking parameters
  
  pseek = 0.7,                                # proportion who seek treatment 
  ptest = 0.8,                                # proportion who receive a test 
  psens = 0.9,                                # sensitivity of diagnostic tool 
  ptreat = 0.8,                               # proportion who receive treatment 
  
  pi_rdt = 0.4,                               #proportion of tests that are RDT 
  rho_V= 1/(4.6*365),                         #rate of waning of immunity from the vaccinated compartment (v)
  paV = 0.77,                                  #Vaccine Efficacy 74%
  cov1 = 0.0,                                 # proportion of population vaccinated
  timevax = 9490,                             # time of vaccination in days (from start of model)
  
  scov = .0,                                  # proportion of population receiving SMC
  sdur=30,                                    # duration of SMC campaign in days
  snapr = 1,                                  # artificially fast rate of SMC to remove people from the compartment. 
  
  # Vector control parameters
  att1 = 1,                                   # proportion of nets distributed in year X circulating in year X
  att2 = 0.7,                                 # proportion of nets distributed in year X circulating in year X+1
  att3 = 0.5,                                 # proportion of nets distributed in year X circulating in year X+2
  itn_use = .6,                               # probability of sleeping under a net
  itn_eff = .6,                               # effectiveness of ITN in reducing transmissions
  ppn = 1.64,                                 # Persons per net
  
  
  # Seasonal and spatial coupling parameters
  amp = 0.5,                                  # amplitude of seasonal forcing
  peak = 1,                                   # peak of seasonal forcing
  
  phi = 1,                                    # phase of seasonal forcing
  kappa = 1/182.64,                           # aging rate from the first age group here 30.44*6months
  kappa1 = 1/1643.76                          # aging rate from the second age group here 30.44days*54 months
  
)                                   

model_start_time <- ymd("2000-01-01") # Start date of the model

# Calculate total duration in days (warm-up + simulation)
warmup_years <- 8
simulation_years <- 27
total_days <- (warmup_years + simulation_years) * 365

# Positive time steps: from 0 to total_days
times <- seq(0, total_days, by = 1) # days



# Run model
Let_go <- deSolve::ode(times = times, y = state, func = Marks_go, parms = parms)













Vod_df <- as_tibble(as.data.frame(Let_go)) %>%
  mutate( P0 = S0 + E0 + A0 + C0 + Ct0 + Sev0 + T0 + H0 + R0,
          P1 = S1 + E1 + A1 + C1 + Ct1 + Sev1 + T1 + H1 + R1 + V + E1V + A1V + C1V + Ct1V + Sev1V + T1V + H1V + R1V + Smcs + Smcr + VSmcs + VSmcr,
          P2 = S2 + E2 + A2 + C2 + Ct2 + Sev2 + T2 + H2 + R2 + F2,
          P = P0 + P1 + P2,
          
          # New cases (incidence)
          cases0 = c(0, diff(cases0)),
          cases1 = c(0, diff(cases1)),
          cases2 = c(0, diff(cases2)),
          casesU5 = c(0, diff(casesU5)),
          cases_total = c(0, diff(Cases)),
          # New severe cases
          severe0 = c(0, diff(Csevere0)),
          severe1 = c(0, diff(Csevere1)),
          severe2 = c(0, diff(Csevere2)),
          severeU5 = c(0, diff(CsevereU5)),
          severe_total = c(0, diff(Csevere)),
          # New deaths
          deaths0 = c(0, diff(deaths0)),
          deaths1 = c(0, diff(deaths1)),
          deaths2 = c(0, diff(deaths2)),
          deathsU5 = c(0, diff(deathsU5)),
          deaths_total = c(0, diff(Tdeaths)),
          
          Incr = c(0, diff((CInc))),
          Inc = c(0, diff((CInc))),
          Snap = c(0, diff(CSnap)),
          CSMC = c(0, diff(SMC)),
          cases = c(0, diff(Cases))) %>%
  pivot_longer(names_to = "variable", cols = !1) %>%
  mutate( level = case_when(
    variable %in% c("S0", "E0", "A0", "C0", "Ct0", "Sev0", "T0", "H0", "R0") ~ 0,
    variable %in% c("S1", "E1", "A1", "C1", "Ct1", "Sev1", "T1", "H1", "R1","V","E1V", "A1V", "C1V", "Ct1V", "Sev1V", "T1V" , "H1V" ,"R1V") ~ 1,
    variable %in% c("S2", "E2", "A2", "C2", "Ct2", "Sev2", "T2", "H2", "R2","F2") ~ 2 ),
    date  = time + model_start_time)

Vod_df %>%
  filter(variable %in% c("P","P0", "P1", "P2"),date >= model_start_time) %>%
  # filter(variable %in% c( "Inc", "CInc"),date >= model_start_time) %>%
  ggplot() +
  geom_line(aes(x = date, y = value,colour = as_factor(variable))) +
  theme_bw() +
  labs(title = "Populations", y = "Population") +
  facet_wrap(~variable, scales = "free", labeller = labeller(variable = c("P" = "Total Population",
                                                                          "P0" = "Population Age 0-5 months",
                                                                          "P1" = "Population Age 6-59 months",
                                                                          "P2" = "Population Age > 60 months")))




# Filter for under-five cases and date range
total_cases_u5_2025_2034 <- Vod_df %>%
  filter(variable == "casesU5", 
         year(date) >= 2026, 
         year(date) <= 2034) %>%
  summarise(total_U5 = sum(value, na.rm = TRUE)) %>%
  pull(total_U5)

# Print the result
cat("Total malaria cases in under-fives from 2025 to 2034:", format(round(total_cases_u5_2025_2034, 0), big.mark = ","), "\n")
# Human Compartments

`%nin%` = Negate(`%in%`)

Vod_df  %>%
  filter(variable %nin% c("P", "P0", "P1","P2", "Inc", "CInc", "CInc0","CInc1","CInc2","SMC", "CSnap", "Snap",
                          "S0", "E0", "A0", "C0", "Ct0", "Sev0", "T0", "H0", "R0",
                          "S2", "E2", "A2", "C2", "Ct2", "Sev2", "T2", "H2", "R2","F2",
                          "Ccases","cases","cases0","cases1","cases2","Cases","vSmcs","CSMC","vSmcr","lambda","Incr","seas")
  ) %>%
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = date, y=value, colour = as_factor(variable)))+
  theme_minimal() +
  labs(title = "Human Compartments", y =("population"))

# Snap4745 & t<5475

Vod_df %>%
  filter(variable %in% c("Snap")) %>%
  ggplot()+
  geom_line(aes(x = date, y=value, colour = as_factor(variable)))+
  theme_minimal() +
  labs(title = "Compartment", y =("population"))


# SMC
Vod_df %>%
  filter(variable %in% c("CSMC")) %>%
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = date, y=value, colour = variable))+
  theme_minimal() +
  labs(title = "Seasonal Malaria Chemoprevention", y =("population"))



# Malaria Incidence
Vod_df  %>%
  filter(variable %in% c("Inc"),) %>%
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = date, y=value, colour = variable))+
  theme_minimal() +
  labs(title = "Malaria Incidence", y =("population"))


# This right here is runing the model monthly.

Bod_df <- as_tibble(as.data.frame(Let_go)) %>%
  mutate(
    P0 = S0 + E0 + A0 + C0 + Ct0 + Sev0 + T0 + H0 + R0,
    P1 = S1 + E1 + A1 + C1 + Ct1 + Sev1 + T1 + H1 + R1 + V + E1V + A1V + C1V + Ct1V + Sev1V + T1V + H1V + R1V + Smcs + Smcr + VSmcs + VSmcr,
    P2 = S2 + E2 + A2 + C2 + Ct2 + Sev2 + T2 + H2 + R2 + F2,
    P = P0 + P1 + P2,
    Snap = c(0, diff(CSnap)),
    CSMC = c(0, diff(SMC)),
    Incr = c(0, diff(CInc)),
    Inc = c(0, diff(CInc)),
    cases = c(0, diff(Cases)),
    date = time + model_start_time,
    month_time = as.Date(format(time + model_start_time, "%Y-%m-01"))
  ) %>%
  pivot_longer(cols = -c(time, date, month_time, P), names_to = "variable", values_to = "value") %>%
  mutate(value = (value / P)*1000) %>%
  mutate( level = case_when(
    variable %in% c("S0", "E0", "A0", "C0", "Ct0", "Sev0", "T0", "H0", "R0") ~ 0,
    variable %in% c("S1", "E1", "A1", "C1", "Ct1", "Sev1", "T1", "H1", "R1","V","E1V", "A1V", "C1V", "Ct1V", "Sev1V", "T1V" , "H1V" ,"R1V","Smcs","Smcr","VSmcs","VSmcr") ~ 1,
    variable %in% c("S2", "E2", "A2", "C2", "Ct2", "Sev2", "T2", "H2", "R2","F2") ~ 2 ),
    date  = time + model_start_time) %>%
  mutate(month_time = as.Date(format(date, "%Y-%m-01"))) %>%
  group_by(month_time, variable) %>%
  summarise(value = sum(value)) #%>%
# ungroup() %>% dplyr::filter(variable %in% c("Smcs","Smcr")) %>%
# group_by(variable,month_time) %>% filter(value > 0) %>%
# summarise(value = sum(value))


Bod_df  %>%
  filter(variable %in% c("S1", "E1", "A1", "C1", "Ct1", "Sev1", "T1", "H1",
                         "R1","V","E1V", "A1V", "C1V", "Ct1V", "Sev1V", "T1V" ,
                         "H1V" ,"R1V","Smcs","Smcr","VSmcs","VSmcr"), between(month_time,
                                                                              as.Date("2025-01-01"),
                                                                              as.Date("2027-12-31")) )%>%
  group_by(variable) %>%
  ggplot()+
  geom_line(aes(x = month_time, y=value, colour = variable))+
  theme_minimal() +
  labs(title = "Human compartment", y =("Per 1000 population per month"))





# Malaria Incidence and Cases ####

annual_metrics <- as_tibble(as.data.frame(Let_go)) %>%
  mutate(
    date = time + model_start_time,
    year = lubridate::year(date),
    P0 = S0 + E0 + A0 + C0 + Ct0 + Sev0 + T0 + H0 + R0,
    P1 = S1 + E1 + A1 + C1 + Ct1 + Sev1 + T1 +  H1 + R1 + V + E1V + A1V + C1V + Ct1V + Sev1V + T1V + H1V + R1V + Smcs + Smcr,
    P2 = S2 + E2 + A2 + C2 + Ct2 + Sev2 + T2 + H2 + R2 +F2,
    PU5 = P0 + P1,
    P = P0 + P1 + P2,
    Inc0 = c(0, diff(CInc0)),
    Inc1 = c(0, diff(CInc1)),
    Inc2 = c(0, diff(CInc2)),
    IncU5  = c(0, diff(CInc5)),
    
    Inc = c(0, diff(CInc))
  ) %>%
  group_by(year) %>%
  summarise(
    total_new_cases = sum(Inc1),
    avg_population = mean(P1),
    annual_incidence_per_1000 = (total_new_cases / avg_population) * 1000
  ) %>%
  pivot_longer(
    cols = c("total_new_cases", "avg_population", "annual_incidence_per_1000"),
    names_to = "variable",
    values_to = "value"
  )

# Plot
annual_metrics %>%
  ggplot(aes(x = year, y = value, color = variable)) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Annual Malaria Metrics",
    x = "Year",
    y = "Value",
    color = "Metric"
  ) +
  facet_wrap(~variable, scales = "free")





# How to do it for the different populations 

annual_metrics <- as_tibble(as.data.frame(Let_go)) %>%
  mutate(
    date = time + model_start_time,
    year = year(date),
    
    
    P0 = S0 + E0 + A0 + C0 + Ct0 + Sev0 + T0 + H0 + R0,
    P1 = S1 + E1 + A1 + C1 + Ct1 + Sev1 + T1 +  H1 + R1 + V + E1V + A1V + C1V + Ct1V + Sev1V + T1V + H1V + R1V+ Smcs+ Smcr + VSmcs + VSmcr,
    P2 = S2 + E2 + A2 + C2 + Ct2 + Sev2 + T2 + H2 + R2 + F2,
    PU5  = P0 + P1,
    P  = P0 + P1 + P2,
    
    Inc0 = c(0, diff(CInc0)),
    Inc1 = c(0, diff(CInc1)),
    Inc2 = c(0, diff(CInc2)),
    IncU5  = c(0, diff(CIncU5)),
    Inc  = c(0, diff(CInc))
  ) %>%
  group_by(year) %>%
  summarise(
    P0 = mean(P0),
    P1 = mean(P1),
    P2 = mean(P2),
    PU5 = mean(PU5),
    P  = mean(P),
    Inc0 = sum(Inc0),
    Inc1 = sum(Inc1),
    Inc2 = sum(Inc2),
    IncU5 = sum(IncU5),
    Inc  = sum(Inc)
  ) %>%
  mutate(
    incidence0 = (Inc0 / P0) * 1000,
    incidence1 = (Inc1 / P1) * 1000,
    incidence2 = (Inc2 / P2) * 1000,
    incidenceU5 = (IncU5 / PU5) * 1000,
    total_incidence = (Inc / P) * 1000,
    
  ) %>%
  # Step 2: Pivot to long format
  pivot_longer(
    cols = c(incidence0, incidence1, incidence2,incidenceU5, total_incidence),
    names_to = "population_group",
    values_to = "incidence_per_1000"
  ) %>%
  mutate(
    population_group = recode(population_group,
                              "incidence0" = "P0",
                              "incidence1" = "P1",
                              "incidence2" = "P2",
                              "incidenceU5" = "PU5",
                              "total_incidence" = "P" ),
    
    
    population_group = factor(population_group, levels = c("P", "P0", "P1","PU5", "P2"))
  )

# Step 3: Custom labels for facets
group_labels <- labeller(
  population_group = c(
    "P" = "Total Population",
    "P0" = "Population Age 0–5 months",
    "P1" = "Population Age 6–59 months",
    "PU5"= "Population PO + P1",
    "P2" = "Population Age > 60 months"
    
  )
)

# Step 4: Plot with facets and labels
ggplot(annual_metrics, aes(x = year, y = incidence_per_1000)) +
  geom_line(color = "#009E73", size = 1) +
  facet_wrap(~population_group, labeller = group_labels, scales = "free_y") +
  theme_bw(base_size = 14) +
  labs(
    title = "Annual Malaria Incidence per 1000 by Age Group",
    x = "Year",
    y = "Incidence per 1000"
  )



# There I am ploting them on the same graph 

ggplot(annual_metrics, aes(x = year, y = incidence_per_1000, color = population_group)) +
  geom_line(size = 1.2) +
  theme_bw(base_size = 14) +
  scale_color_manual(
    values = c(
      "P" = "#000000",       # Black for Total Population
      "P0" = "#E69F00",      # Orange
      "P1" = "#56B4E9",      # Blue
      "P2" = "#009E73",
      "PU5"= "#009"# 
    ),
    labels = c(
      "P" = "Total Population",
      "P0" = "Age 0–5 months",
      "P1" = "Age 6–59 months",
      "P2" = "Age > 60 months",
      "PU5" = "Age < 60 months"
    )
  ) +
  labs(
    title = "Annual Malaria Incidence per 1000 Population",
    x = "Year",
    y = "Incidence per 1000",
    color = "Population Group"
  )


# Filter the data for the total population (P)

# Filter the data for the total population (P)
total_population_metrics <- annual_metrics %>%
  filter(population_group == "P",
         year >= 2000) # Adjust the years as needed

# Assuming you have observed data in a data frame called 'observed_data'
# with columns 'year' and 'observed_incidence'

# Plot the total population data along with observed data
ggplot() +
  # Line for modeled incidence
  geom_line(data = total_population_metrics, 
            aes(x = year, y = incidence_per_1000), 
            color = "#009E73", size = 1, linetype = "solid") +
  # Points for observed incidence
  geom_point(data = observed_data, 
             aes(x = year, y = observed_incidence), 
             color = "#D55E00", size = 2) +
  theme_bw(base_size = 14) +
  labs(
    title = "Annual Malaria Incidence per 1000: Model vs Observed Data",
    x = "Year",
    y = "Incidence per 1000"
  ) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )




# Lets do a scenario analysis and see what is missing from the plot 



##################################################################
#MODEL fitting here 
######################################################################





















##########################################################################
# Scenario Analysis for Malaria Model
##########################################################################


base_params <- parms
scenario_param_names <- c("itn_use", "cov1", "scov")

# ----------- 2. Scenario definitions -----------
scenario_analysis <- list(
  Baseline  = c(itn_use = base_params["itn_use"], cov1 = 0, scov = 0),
  Scenario1 = c(itn_use = base_params["itn_use"], cov1 =0.8, scov = 0.8),
  Scenario2 = c(itn_use = base_params["itn_use"], cov1 = 0.8, scov = 0),
  Scenario3 = c(itn_use = base_params["itn_use"], cov1 = 0.0, scov = 0.8),
  Scenario4 = c(itn_use = base_params["itn_use"], cov1 = 0.5, scov = 0.5)
  
)

# ----------- 3. Scenario runner -----------
run_scenario <- function(scenario_values, base_params, times, state) {
  params <- base_params
  params[scenario_param_names] <- scenario_values
  result <- ode(times = times, y = state, func = Marks_go, parms = params)
  as.data.frame(result)
}

# ----------- 4. Run all scenarios -----------
scenario_outputs <- bind_rows(
  lapply(names(scenario_analysis), function(name) {
    out <- run_scenario(scenario_analysis[[name]], base_params, times, state)
    out$scenario <- name
    out
  })
)

# ----------- 5. Postprocess results -----------
# Remove all the cases code from your scenario analysis pipeline


scenario_df <- scenario_outputs %>%
  mutate(
    date = time + model_start_time,
    year = lubridate::year(date),
    P0 = S0 + E0 + A0 + C0 + Ct0 + Sev0 + T0 + H0 + R0,
    P1 = S1 + E1 + A1 + C1 + Ct1 + Sev1 + T1 +  H1 + R1 + V + E1V + A1V + C1V + Ct1V + Sev1V + T1V + H1V + R1V+ Smcs+ Smcr + VSmcs + VSmcr, 
    P2 <-  S2 + E2 + A2 + C2 + Ct2 + Sev2 + T2 + H2 + R2 + F2 ,
    P2 = S2 + E2 + A2 + C2 + Ct2 + Sev2 + T2 + H2 + R2 + F2,
    PU5  = P0 + P1,
    P  = P0 + P1 + P2,
    
    Inc0 = c(0, diff(CInc0)),
    Inc1 = c(0, diff(CInc1)),
    Inc2 = c(0, diff(CInc2)),
    IncU5  = c(0, diff(CIncU5)),
    Inc  = c(0, diff(CInc))
  ) %>%
  group_by(year,scenario) %>%
  summarise(
    P0 = mean(P0),
    P1 = mean(P1),
    P2 = mean(P2),
    PU5 = mean(PU5),
    P  = mean(P),
    Inc0 = sum(Inc0),
    Inc1 = sum(Inc1),
    Inc2 = sum(Inc2),
    IncU5 = sum(IncU5),
    Inc  = sum(Inc)
  ) %>%
  mutate(
    incidence0 = (Inc0 / P0) * 1000,
    incidence1 = (Inc1 / P1) * 1000,
    incidence2 = (Inc2 / P2) * 1000,
    incidenceU5 = (IncU5 / PU5) * 1000,
    total_incidence = (Inc / P) * 1000,
    
  ) %>%
  pivot_longer(
    cols = c(incidence0, incidence1, incidence2, incidenceU5, total_incidence,
    ),
    names_to = "population_group",
    values_to = "incidence_per_1000"
  ) %>%
  mutate(
    population_group = recode(population_group,
                              "incidence0" = "P0", 
                              "incidence1" = "P1", 
                              "incidence2" = "P2", 
                              "incidenceU5" = "PU5",
                              "total_incidence" = "P"),
    population_group = factor(population_group, levels = c("P","PU5" ,"P0", "P1", "P2")),
    scenario = factor(scenario, levels = names(scenario_analysis))  # Only defined scenarios
  )

# --- 5. Confidence intervals ---
ci_summary <- scenario_df %>%
  group_by(year, scenario, population_group) %>%
  summarise(
    upper95 = quantile(incidence_per_1000, 0.95, na.rm=TRUE),
    lower5  = quantile(incidence_per_1000, 0.05, na.rm=TRUE),
    median  = median(incidence_per_1000, na.rm=TRUE),
    .groups = "drop"
  )

# --- 6. Define colors for scenarios ---
scenario_colors <- c(
  Baseline  = "gray30",
  Scenario1 = "#19c2c2",
  Scenario2 = "#1e90ff",
  Scenario3 = "#a259bc",
  Scenario4 = "#f4a259"
)

# --- 7. Plot: Faceted comparison of each Scenario vs Baseline ---

scenarios_to_plot <- setdiff(levels(scenario_df$scenario), "Baseline")

plots <- lapply(scenarios_to_plot, function(scen) {
  data_subset <- ci_summary %>%
    filter(population_group == "P", year >= 2020, scenario %in% c("Baseline", scen))
  
  
  
  ggplot(data_subset, aes(x = year)) +
    geom_ribbon(aes(ymin = lower5, ymax = upper95), alpha = 0.3) +
    geom_vline(xintercept = 2034, linetype = "dotted", color = "red", size = 1)+
    geom_line(aes(y = median, color = scenario), size = 1) +
    scale_color_manual(
      values = scenario_colors,
      labels = c(
        "Scenario1" = "Scenario1: 80% vaccination +  80% SMC",
        "Scenario2" = "Scenario2: 80% vaccination + 0% SMC",
        "Scenario3" = "Scenario3: 0% vaccination + 80% SMC",
        "Scenario4" = "Scenario4: 50% vaccination + 50% SMC"
      ),
      guide = guide_legend(ncol = 1)) +
    scale_x_continuous(
      breaks = seq(min(data_subset$year), max(data_subset$year)), 
      labels = as.character(seq(min(data_subset$year), max(data_subset$year)))
    ) +
    labs(
      title = paste("Baseline vs", scen, "(Incidence per 1000, Total Population)"),
      x = "Years",
      y = "Incidence per 1000",
      color = "Scenarios",
      linetype = "Scenarios"
    ) +
    ylim(0, NA) +  
    theme_bw(base_size = 12) +
    theme(
      axis.title.y = element_text(face = "bold", size = 12),
      axis.text.y = element_text(face = "bold", size = 12),
      axis.title.x = element_text(face = "bold", size = 12),
      axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 0.8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 8), 
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      # legend.title = element_text(size = 8),
      legend.key.width = unit(0.5, "cm"),
      legend.spacing.x = unit(0.5, "cm"))
  # legend.justification = "center")
})


# --- 8. Combine plots in a grid ---
do.call(grid.arrange, c(plots, ncol = 2))






#######################################
#Same plot
####################################
data_subset <- ci_summary %>%
  filter(population_group == "P", year >= 2020)


ggplot(data_subset, aes(x = year, y = median, color = scenario)) +
  geom_ribbon(aes(ymin = lower5, ymax = upper95, fill = scenario), alpha = 0.2) +
  geom_vline(xintercept = 2034, linetype = "dotted", color = "red", size = 1) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = scenario_colors,
    labels = c(
      "Baseline"  = "Baseline: 0% vaccination",
      "Scenario1" = "80% vaccination + 80% SMC",
      "Scenario2" = "80% vaccination + 0% SMC",
      "Scenario3" = "0% vaccination + 80% SMC",
      "Scenario4" = "50% vaccination + 50% SMC"
    ),
    guide = guide_legend(ncol = 2)
  ) +
  scale_fill_manual(
    values = scenario_colors,
    labels = c(
      "Baseline"  = "Baseline: 0% vaccination",
      "Scenario1" = "80% vaccination + 80% SMC",
      "Scenario2" = "80% vaccination + 0% SMC",
      "Scenario3" = "0% vaccination + 80% SMC",
      "Scenario4" = "50% vaccination + 50% SMC"
    ),
    guide = guide_legend(ncol = 2)
  ) +
  scale_x_continuous(
    breaks = seq(min(data_subset$year), max(data_subset$year)), 
    labels = as.character(seq(min(data_subset$year), max(data_subset$year)))
  ) +
  labs( x = "Years",
        y = "Incidence per 1000 in the total population",
        color = "Scenarios",
        fill = "Scenarios"
  ) +
  ylim(0, NA) +
  theme_bw(base_size = 14) +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.key.width = unit(1.5, "cm"),
    legend.spacing.x = unit(0.5, "cm")
  )


##################################################################
# the following code is the bar chart that simulates our scenario analysis
#################################################################
# 
# 
# 
# ci_summary_filtered <- ci_summary %>% 
#   filter(population_group == "P", year >= 2025 & year <= 3034)
# 
# # Aggregate across years per scenario
# total_summary <- ci_summary_filtered %>%
#   group_by(scenario) %>%
#   summarise(
#     avg_median = mean(median),
#     avg_upper95 = mean(upper95),
#     avg_lower5 = mean(lower5),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     error_upper = avg_upper95 - avg_median,
#     error_lower = avg_median - avg_lower5
#   )
# 
# # 1. Set the scenario factor levels in the desired order
# total_summary$scenario <- factor(total_summary$scenario, levels = c("Baseline", "Scenario1", "Scenario2", "Scenario3", "Scenario4"))
# 
# # 2. Create the caption text
# # scenario_caption <- paste(
# #   "Scenario Key:",
# #   "Baseline: cov1 = 0, scov = 0",
# #   "Scenario1: cov1 = 0.8, scov = 0.80",
# #   "Scenario2: cov1 = 0.8, scov = 0",
# #   "Scenario3: cov1 = 0.0, scov = 0.80",
# #   "Scenario4: cov1 = 0.5, scov = 0.5",
# #   sep = "\n"
# # )
# 
# # 3. Plot
# ggplot(total_summary, aes(x = scenario, y = avg_median, fill = scenario)) +
#   geom_bar(stat = "identity") +
#   geom_errorbar(aes(ymin = avg_median - error_lower, ymax = avg_median + error_upper), width = 0.2) +
#   scale_fill_manual(values = scenario_colors, name = "Scenario") +
#   labs(
#     title = "Total Incidence per 1000",
#     x = "Scenario",
#     y = "Average Incidence per 1000",
#     caption = scenario_caption
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "right",
#     plot.caption = element_text(hjust = 0, size = 10, face = "italic")
#   )


##############################################################
# This is the bar chat that simulates our scenario analysis
##############################################################


# 1. Filter for total population (P) and relevant years
ci_summary_filtered <- ci_summary %>% 
  filter(population_group == "P", year >= 2026 & year <= 2034)

# 2. Aggregate incidence by scenario across years
total_summary <- ci_summary_filtered %>%
  group_by(scenario) %>%
  summarise(
    avg_median = sum(median, na.rm = TRUE),
    avg_upper95 = sum(upper95, na.rm = TRUE),
    avg_lower5 = sum(lower5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    error_upper = avg_upper95 - avg_median,
    error_lower = avg_median - avg_lower5
  )

# 3. Set scenario order for consistent plotting
total_summary$scenario <- factor(total_summary$scenario, levels = c("Baseline", "Scenario1", "Scenario2", "Scenario3", "Scenario4"))

# 4. Calculate differences from Baseline
baseline_value <- total_summary %>%
  filter(scenario == "Baseline") %>%
  pull(avg_median)

total_summary <- total_summary %>%
  mutate(
    diff_from_baseline = (avg_median - baseline_value)/ baseline_value * 100, # Calculate percentage difference
    label_text = ifelse(
      scenario != "Baseline",
      sprintf("Δ = %+0.1f", diff_from_baseline),
      ""
    )
  )


# 6. Plot
ggplot(total_summary, aes(x = scenario, y = avg_median, fill = scenario)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", aes(alpha = scenario != "Baseline")) +
  geom_errorbar(aes(ymin = avg_lower5, ymax = avg_upper95), width = 0.2, size = 1) +
  geom_text(aes(label = label_text), vjust = -0.8, size = 4.2, fontface = "bold") +
  scale_fill_manual(values = scenario_colors, name = "Scenario") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5), guide = "none") +
  labs(
    title = "Total Malaria Incidence per 1000 (2026–2034)",
    subtitle = "Δ = difference vs Baseline",
    x = "Scenario",
    y = "total Incidence per 1000",
    caption = scenario_caption
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.caption = element_text(hjust = 0, size = 10, face = "italic"),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(face = "italic", size = 12)
  )


########################################################################
# this is 6-60 months age group
#########################################################################
plots <- lapply(scenarios_to_plot, function(scen) {
  data_subset <- ci_summary %>%
    filter(population_group == "PU5", year >= 2020, scenario %in% c("Baseline", scen))
  
  
  
  ggplot(data_subset, aes(x = year)) +
    geom_ribbon(aes(ymin = lower5, ymax = upper95), alpha = 0.3) +
    geom_vline(xintercept = 2034, linetype = "dotted", color = "red", size = 1)+
    geom_line(aes(y = median, color = scenario), size = 1) +
    scale_color_manual(
      values = scenario_colors,
      labels = c(
        "Scenario1" = "80% vaccination +  80% SMC",
        "Scenario2" = "80% vaccination + 0% SMC",
        "Scenario3" = "0% vaccination + 80% SMC",
        "Scenario4" = "50% vaccination + 50% SMC"
      ),
      guide = guide_legend(ncol = 1)) +
    scale_x_continuous(
      breaks = seq(min(data_subset$year), max(data_subset$year)), 
      labels = as.character(seq(min(data_subset$year), max(data_subset$year)))
    ) +
    labs(
      title = paste("Baseline vs", scen, "(Incidence per 1000, 6-60 months Population)"),
      x = "Years",
      y = "Incidence per 1000",
      color = "Scenarios",
      linetype = "Scenarios"
    ) +
    ylim(0, NA) +  
    theme_bw(base_size = 12) +
    theme(
      axis.title.y = element_text(face = "bold", size = 12),
      axis.text.y = element_text(face = "bold", size = 12),
      axis.title.x = element_text(face = "bold", size = 12),
      axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 0.8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 8), 
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      # legend.title = element_text(size = 8),
      legend.key.width = unit(0.5, "cm"),
      legend.spacing.x = unit(0.5, "cm"))
  # legend.justification = "center")
})


# --- 8. Combine plots in a grid ---
do.call(grid.arrange, c(plots, ncol = 2))



####################################

#################################

data_subset <- ci_summary %>%
  filter(population_group == "PU5", year >= 2020)


ggplot(data_subset, aes(x = year, y = median, color = scenario)) +
  geom_ribbon(aes(ymin = lower5, ymax = upper95, fill = scenario), alpha = 0.2) +
  geom_vline(xintercept = 2034, linetype = "dotted", color = "red", size = 1) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = scenario_colors,
    labels = c(
      "Baseline"  = "Baseline: 0% vaccination&SMC",
      "Scenario1" = "80% vaccination + 80% SMC",
      "Scenario2" = "80% vaccination + 0% SMC",
      "Scenario3" = "0% vaccination + 80% SMC",
      "Scenario4" = "50% vaccination + 50% SMC"
    ),
    guide = guide_legend(ncol = 2)
  ) +
  scale_fill_manual(
    values = scenario_colors,
    labels = c(
      "Baseline"  = "Baseline: 0% vaccination&SMC",
      "Scenario1" = "80% vaccination + 80% SMC",
      "Scenario2" = "80% vaccination + 0% SMC",
      "Scenario3" = "0% vaccination + 80% SMC",
      "Scenario4" = "50% vaccination + 50% SMC"
    ),
    guide = guide_legend(ncol = 2)
  ) +
  scale_x_continuous(
    breaks = seq(min(data_subset$year), max(data_subset$year)), 
    labels = as.character(seq(min(data_subset$year), max(data_subset$year)))
  ) +
  labs(
    title = "Malaria Incidence per 1000, < 5 years Population",
    x = "Years",
    y = "Incidence per 1000",
    color = "Scenarios",
    fill = "Scenarios"
  ) +
  ylim(0, NA) +
  theme_bw(base_size = 14) +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.key.width = unit(1.5, "cm"),
    legend.spacing.x = unit(0.5, "cm")
  )



##################################################################


#####################################################################

ci_summary_filtered <- ci_summary %>% 
  filter(population_group == "PU5", year >= 2026 & year <= 2034)

# 2. Aggregate incidence by scenario across years
total_summary <- ci_summary_filtered %>%
  group_by(scenario) %>%
  summarise(
    avg_median = sum(median, na.rm = TRUE),
    avg_upper95 = sum(upper95, na.rm = TRUE),
    avg_lower5 = sum(lower5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    error_upper = avg_upper95 - avg_median,
    error_lower = avg_median - avg_lower5
  )

# 3. Set scenario order for consistent plotting
total_summary$scenario <- factor(total_summary$scenario, levels = c("Baseline", "Scenario1", "Scenario2", "Scenario3", "Scenario4"))

# 4. Calculate differences from Baseline
baseline_value <- total_summary %>%
  filter(scenario == "Baseline") %>%
  pull(avg_median)

total_summary <- total_summary %>%
  mutate(
    diff_from_baseline = (avg_median - baseline_value)/ baseline_value * 100, # Calculate percentage difference
    label_text = ifelse(
      scenario != "Baseline",
      sprintf("Δ = %+0.1f", diff_from_baseline),
      ""
    )
  )

# 5. Create scenario explanation caption
scenario_caption <- paste(
  "Scenario Key:",
  "Baseline: cov1 = 0, scov = 0",
  "Scenario1: cov1 = 0.8, scov = 0.80",
  "Scenario2: cov1 = 0.8, scov = 0",
  "Scenario3: cov1 = 0.0, scov = 0.80",
  "Scenario4: cov1 = 0.5, scov = 0.5",
  sep = "\n"
)

# 6. Plot
ggplot(total_summary, aes(x = scenario, y = avg_median, fill = scenario)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", aes(alpha = scenario != "Baseline")) +
  geom_errorbar(aes(ymin = avg_lower5, ymax = avg_upper95), width = 0.2, size = 1) +
  geom_text(aes(label = label_text), vjust = -0.8, size = 4.2, fontface = "bold") +
  scale_fill_manual(values = scenario_colors, name = "Scenario") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5), guide = "none") +
  labs(
    subtitle = "Δ = difference vs Baseline",
    x = "Scenario",
    y = "total Incidence per 1000 in children under five",
    caption = scenario_caption
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.caption = element_text(hjust = 0, size = 10, face = "italic"),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(face = "italic", size = 12)
  )

##############################################################################
# Now here we want to look at he cases only 
##############################################################################

scenario_Cases <- scenario_outputs %>%
  mutate(
    date = time + model_start_time,
    year = lubridate::year(date),
    # New cases (incidence)
    cases0 = c(0, diff(cases0)),
    cases1 = c(0, diff(cases1)),
    cases2 = c(0, diff(cases2)),
    casesU5 = c(0, diff(casesU5)),
    cases_total = c(0, diff(Cases)),
    # New severe cases
    severe0 = c(0, diff(Csevere0)),
    severe1 = c(0, diff(Csevere1)),
    severe2 = c(0, diff(Csevere2)),
    severeU5 = c(0, diff(CsevereU5)),
    severe_total = c(0, diff(Csevere)),
    # New deaths
    deaths0 = c(0, diff(deaths0)),
    deaths1 = c(0, diff(deaths1)),
    deaths2 = c(0, diff(deaths2)),
    deathsU5 = c(0, diff(deathsU5)),
    deaths_total = c(0, diff(Tdeaths))
  ) %>%
  group_by(year, scenario) %>%
  summarise(
    cases0 = sum(cases0),
    cases1 = sum(cases1),
    cases2 = sum(cases2),
    casesU5 = sum(casesU5),
    cases_total = sum(cases_total),
    severe0 = sum(severe0),
    severe1 = sum(severe1),
    severe2 = sum(severe2),
    severeU5 = sum(severeU5),
    severe_total = sum(severe_total),
    deaths0 = sum(deaths0),
    deaths1 = sum(deaths1),
    deaths2 = sum(deaths2),
    deathsU5 = sum(deathsU5),
    deaths_total = sum(deaths_total),
    .groups = "drop"
  )

# Set scenario order if you want
cases_summary <- scenario_Cases %>%
  filter(year >= 2023 & year <= 2034) %>%
  group_by(scenario) %>%
  summarise(
    cases_0_5mo   = sum(cases0),
    cases_6_60mo  = sum(cases1),
    cases_5plus   = sum(cases2),
    casesU5_total = sum(casesU5),
    cases_total   = sum(cases_total),
    severe_0_5mo   = sum(severe0),
    severe_6_60mo  = sum(severe1),
    severe_5plus   = sum(severe2),
    severeU5_total = sum(severeU5),
    severe_total   = sum(severe_total),
    deaths_0_5mo   = sum(deaths0),
    deaths_6_60mo  = sum(deaths1),
    deaths_5plus   = sum(deaths2),
    deathsU5_total = sum(deathsU5),
    deaths_total   = sum(deaths_total),
    .groups = "drop"
  )

cases_long <- scenario_Cases %>%
  filter(year >= 2023 & year <= 2034) %>%
  pivot_longer(
    cols = c(cases0, cases1, cases2, casesU5, cases_total,
             severe0, severe1, severe2, severeU5, severe_total,
             deaths0, deaths1, deaths2, deathsU5, deaths_total),
    names_to = "metric",
    values_to = "value"
  )



##################################
#We want to use the code for printing the results. This generates the number we want to see
#################################
cases_summary_2026_2030 <- scenario_Cases %>%
  filter(year >= 2026 & year <= 2034) %>%
  group_by(scenario) %>%
  summarise(
    total_deathsU5 = sum(casesU5),
    .groups = "drop"
  )

# Arrange with Baseline first, then others
cases_summary_2026_2030$scenario <- factor(
  cases_summary_2026_2030$scenario,
  levels = c("Baseline", setdiff(unique(cases_summary_2026_2030$scenario), "Baseline"))
)

cases_summary_2026_2030 <- cases_summary_2026_2030 %>%
  arrange(scenario)

# View result
print(cases_summary_2026_2030)

#####################################################



#################################################
# Here we look at the total cases. 
#################################################
Pcases_total1 <- ggplot(cases_long %>% filter(metric == "cases_total"),
                        aes(x = year, y = value,color = scenario)) +
  geom_vline(xintercept = 2034, linetype = "dotted", color = "red", size = 1) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = scenario_colors,
    labels = c(
      "Scenario1" = "Scenario1: 80% vaccination +  80% SMC",
      "Scenario2" = "Scenario2: 80% vaccination + 0% SMC",
      "Scenario3" = "Scenario3: 0% vaccination + 80% SMC",
      "Scenario4" = "Scenario4: 50% vaccination + 50% SMC"
    ),
    guide = guide_legend(ncol = 2)
  ) +
  scale_x_continuous(
    breaks = seq(min(cases_long$year), max(cases_long$year)),
    labels = as.character(seq(min(cases_long$year), max(cases_long$year)))
  ) +
  labs(
    x = "Years",
    y = "Annual clinical cases in total population",
    color = "Scenarios"
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    labels = scales::label_comma()) +
  theme_bw(base_size = 14) +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.key.width = unit(1.5, "cm"),
    legend.spacing.x = unit(0.5, "cm")
  )

# ggsave("Pcases_total1.tiff", plot = Pcases_total1, width=12, height=6, dpi = 600)


########################################
# Bar plot for the Total cases 
#######################################
cases_summary <- cases_long %>%
  filter(metric == "cases_total", year >= 2026, year <= 2034) %>%
  group_by(scenario) %>%
  summarise(
    total_cases = sum(value),
    .groups = "drop"
  )

# Calculate percent difference from Baseline
baseline_cases <- cases_summary %>% filter(scenario == "Baseline") %>% pull(total_cases)

cases_summary <- cases_summary %>%
  mutate(
    diff_from_baseline = (total_cases - baseline_cases) / baseline_cases * 100,
    label_text = ifelse(
      scenario != "Baseline",
      sprintf("Δ = %+0.1f", diff_from_baseline),
      ""
    )
  )


# Bar plot code
Pcases_bartotal1 <- ggplot(cases_summary, aes(x = scenario, y = total_cases, fill = scenario)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", aes(alpha = scenario != "Baseline")) +
  geom_text(aes(label = label_text), vjust = -0.8, size = 4.2, fontface = "bold") +
  scale_fill_manual(
    values = scenario_colors,
    name = "Scenarios",
    labels = c(
      "Baseline"   = "Baseline",
      "Scenario1"  = "Scenario1:\n80% vaccination +\n80% SMC",  # Added \n for line breaks
      "Scenario2"  = "Scenario2:\n80% vaccination +\n0% SMC",
      "Scenario3"  = "Scenario3:\n0% vaccination +\n80% SMC",
      "Scenario4"  = "Scenario4:\n50% vaccination +\n50% SMC"
    )
  ) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5), guide = "none") +
  labs(
    subtitle = "Δ = Difference vs Baseline",
    x = "Scenario",
    y = "Total clinical cases in the overall population"
  ) +
  scale_y_continuous(
    labels = label_comma(),
    expand = expansion(mult = c(0, 0.1))  # Adds 10% space above top bar
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.8, "cm"),
    legend.spacing.x = unit(0.5, "cm"),
    plot.caption = element_text(hjust = 0, size = 10, face = "italic"),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(face = "italic", size = 12),
    legend.key.height = unit(0.8, "cm")  # Add this to control key height
  )


# ggsave("Pcases_bartotal1.tiff", plot =Pcases_bartotal1, width=12, height=6, dpi = 600)



###############################################
# The plot below is for the age group < 5
###############################################


Pcases_U5 <-  ggplot(cases_long %>% filter(metric == "casesU5"),
                     aes(x = year, y = value,color = scenario)) +
  geom_vline(xintercept = 2034, linetype = "dotted", color = "red", size = 1) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = scenario_colors,
    labels = c(
      "Scenario1" = "Scenario1: 80% vaccination +  80% SMC",
      "Scenario2" = "Scenario2: 80% vaccination + 0% SMC",
      "Scenario3" = "Scenario3: 0% vaccination + 80% SMC",
      "Scenario4" = "Scenario4: 50% vaccination + 50% SMC"
    ),
    guide = guide_legend(ncol = 2)
  ) +
  scale_x_continuous(
    breaks = seq(min(cases_long$year), max(cases_long$year)),
    labels = as.character(seq(min(cases_long$year), max(cases_long$year)))
  ) +
  labs(
    x = "Years",
    y = "Annual clinical cases in children under 5",
    color = "Scenarios"
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    labels = scales::label_comma()) +
  theme_bw(base_size = 14) +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.key.width = unit(1.5, "cm"),
    legend.spacing.x = unit(0.5, "cm")
  )

ggsave("Pcases_U5.tiff", plot = Pcases_U5, width=12, height=6, dpi = 600)


########################################
# Bar plot for the unders five  cases 
#######################################
cases_summary <- cases_long %>%
  filter(metric == "casesU5", year >= 2026, year <= 2034) %>%
  group_by(scenario) %>%
  summarise(
    total_u5cases = sum(value),
    .groups = "drop"
  )

# Calculate percent difference from Baseline
baseline_cases <- cases_summary %>% filter(scenario == "Baseline") %>% pull(total_u5cases)

cases_summary <- cases_summary %>%
  mutate(
    diff_from_baseline = (total_u5cases - baseline_cases) / baseline_cases * 100,
    label_text = ifelse(
      scenario != "Baseline",
      sprintf("Δ = %+0.1f", diff_from_baseline),
      ""
    )
  )



# Bar plot code
Pcases_barU5 <- ggplot(cases_summary, aes(x = scenario, y = total_u5cases, fill = scenario)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", aes(alpha = scenario != "Baseline")) +
  geom_text(aes(label = label_text), vjust = -0.8, size = 4.2, fontface = "bold") +
  scale_fill_manual(
    values = scenario_colors,
    name = "Scenarios",
    labels = c(
      "Baseline"   = "Baseline",
      "Scenario1"  = "Scenario1:\n80% vaccination +\n80% SMC",  # Added \n for line breaks
      "Scenario2"  = "Scenario2:\n80% vaccination +\n0% SMC",
      "Scenario3"  = "Scenario3:\n0% vaccination +\n80% SMC",
      "Scenario4"  = "Scenario4:\n50% vaccination +\n50% SMC"
    )
  ) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5), guide = "none") +
  labs(
    subtitle = "Δ = Difference vs Baseline",
    x = "Scenario",
    y = "Total clinical cases in children under 5"
  ) +
  scale_y_continuous(
    labels = label_comma(),
    expand = expansion(mult = c(0, 0.1))  # Adds 10% space above top bar
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.8, "cm"),
    legend.spacing.x = unit(0.5, "cm"),
    plot.caption = element_text(hjust = 0, size = 10, face = "italic"),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(face = "italic", size = 12),
    legend.key.height = unit(0.8, "cm")  # Add this to control key height
  )


ggsave("Pcasesbar_U5.tiff", plot = Pcases_barU5, width=12, height=6, dpi = 600)


##############################################################
# We now at the Total Severe compartment 
######################b#######################################

PSevere_total <- ggplot(cases_long %>% filter(metric == "severe_total"),
                        aes(x = year, y = value,color = scenario)) +
  geom_vline(xintercept = 2034, linetype = "dotted", color = "red", size = 1) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = scenario_colors,
    labels = c(
      "Scenario1" = "Scenario1: 80% vaccination +  80% SMC",
      "Scenario2" = "Scenario2: 80% vaccination + 0% SMC",
      "Scenario3" = "Scenario3: 0% vaccination + 80% SMC",
      "Scenario4" = "Scenario4: 50% vaccination + 50% SMC"
    ),
    guide = guide_legend(ncol = 2)
  ) +
  scale_x_continuous(
    breaks = seq(min(cases_long$year), max(cases_long$year)),
    labels = as.character(seq(min(cases_long$year), max(cases_long$year)))
  ) +
  labs(
    x = "Years",
    y = "Annual severe cases in the total population",
    color = "Scenarios"
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    labels = scales::label_comma()) +
  theme_bw(base_size = 14) +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.key.width = unit(1.5, "cm"),
    legend.spacing.x = unit(0.5, "cm")
  )

ggsave("PSevere_total.tiff", plot = PSevere_total, width=12, height=6, dpi = 600)

########################################
# Bar plot for the Total cases 
#######################################
cases_summary <- cases_long %>%
  filter(metric == "severe_total", year >= 2026, year <= 2034) %>%
  group_by(scenario) %>%
  summarise(
    total_severe_cases = sum(value),
    .groups = "drop"
  )

# Calculate percent difference from Baseline
baseline_cases <- cases_summary %>% filter(scenario == "Baseline") %>% pull(total_severe_cases)

cases_summary <- cases_summary %>%
  mutate(
    diff_from_baseline = (total_severe_cases - baseline_cases) / baseline_cases * 100,
    label_text = ifelse(
      scenario != "Baseline",
      sprintf("Δ = %+0.1f", diff_from_baseline),
      ""
    )
  )


# Bar plot code
PSeverebar_total <-  ggplot(cases_summary, aes(x = scenario, y = total_severe_cases, fill = scenario)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", aes(alpha = scenario != "Baseline")) +
  geom_text(aes(label = label_text), vjust = -0.8, size = 4.2, fontface = "bold") +
  scale_fill_manual(
    values = scenario_colors,
    name = "Scenarios",
    labels = c(
      "Baseline"   = "Baseline",
      "Scenario1"  = "Scenario1:\n80% vaccination +\n80% SMC",  # Added \n for line breaks
      "Scenario2"  = "Scenario2:\n80% vaccination +\n0% SMC",
      "Scenario3"  = "Scenario3:\n0% vaccination +\n80% SMC",
      "Scenario4"  = "Scenario4:\n50% vaccination +\n50% SMC"
    )
  ) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5), guide = "none") +
  labs(
    subtitle = "Δ = Difference vs Baseline",
    x = "Scenario",
    y = "Total severe cases in the overall population"
  ) +
  scale_y_continuous(
    labels = label_comma(),
    expand = expansion(mult = c(0, 0.1))  # Adds 10% space above top bar
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.8, "cm"),
    legend.spacing.x = unit(0.5, "cm"),
    plot.caption = element_text(hjust = 0, size = 10, face = "italic"),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(face = "italic", size = 12),
    legend.key.height = unit(0.8, "cm")  # Add this to control key height
  )


ggsave("PSeverebar_total.tiff", plot = PSeverebar_total, width=12, height=6, dpi = 600)

###########################################################
# We also have the Severe cases in under five 
##########################################################
PSevere_U5 <- ggplot(cases_long %>% filter(metric == "severeU5"),
                     aes(x = year, y = value,color = scenario)) +
  geom_vline(xintercept = 2034, linetype = "dotted", color = "red", size = 1) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = scenario_colors,
    labels = c(
      "Scenario1" = "Scenario1: 80% vaccination +  80% SMC",
      "Scenario2" = "Scenario2: 80% vaccination + 0% SMC",
      "Scenario3" = "Scenario3: 0% vaccination + 80% SMC",
      "Scenario4" = "Scenario4: 50% vaccination + 50% SMC"
    ),
    guide = guide_legend(ncol = 2)
  ) +
  scale_x_continuous(
    breaks = seq(min(cases_long$year), max(cases_long$year)),
    labels = as.character(seq(min(cases_long$year), max(cases_long$year)))
  ) +
  labs(
    x = "Years",
    y = "Annual severe cases in children under 5",
    color = "Scenarios"
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    labels = scales::label_comma()) +
  theme_bw(base_size = 14) +
  theme_bw(base_size = 14) +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.key.width = unit(1.5, "cm"),
    legend.spacing.x = unit(0.5, "cm")
  )

ggsave("PSevere_U5.tiff", plot = PSevere_U5, width=12, height=6, dpi = 600)

###########################################
# Bar plot for the under five severe cases 
###########################################
cases_summary <- cases_long %>%
  filter(metric == "severeU5", year >= 2026, year <= 2034) %>%
  group_by(scenario) %>%
  summarise(
    total_u5_sevcases = sum(value),
    .groups = "drop"
  )

# Calculate percent difference from Baseline
baseline_cases <- cases_summary %>% filter(scenario == "Baseline") %>% pull(total_u5_sevcases)

cases_summary <- cases_summary %>%
  mutate(
    diff_from_baseline = (total_u5_sevcases - baseline_cases) / baseline_cases * 100,
    label_text = ifelse(
      scenario != "Baseline",
      sprintf("Δ = %+0.1f", diff_from_baseline),
      ""
    )
  )



# Bar plot code
PSeverebar_U5 <- ggplot(cases_summary, aes(x = scenario, y = total_u5_sevcases, fill = scenario)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", aes(alpha = scenario != "Baseline")) +
  geom_text(aes(label = label_text), vjust = -0.8, size = 4.2, fontface = "bold") +
  scale_fill_manual(
    values = scenario_colors,
    name = "Scenarios",
    labels = c(
      "Baseline"   = "Baseline",
      "Scenario1"  = "Scenario1:\n80% vaccination +\n80% SMC",  # Added \n for line breaks
      "Scenario2"  = "Scenario2:\n80% vaccination +\n0% SMC",
      "Scenario3"  = "Scenario3:\n0% vaccination +\n80% SMC",
      "Scenario4"  = "Scenario4:\n50% vaccination +\n50% SMC"
    )
  ) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5), guide = "none") +
  labs(
    subtitle = "Δ = Difference vs Baseline",
    x = "Scenario",
    y = "Total death in children under 5"
  ) +
  scale_y_continuous(
    labels = label_comma(),
    expand = expansion(mult = c(0, 0.1))  # Adds 10% space above top bar
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.8, "cm"),
    legend.spacing.x = unit(0.5, "cm"),
    plot.caption = element_text(hjust = 0, size = 10, face = "italic"),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(face = "italic", size = 12),
    legend.key.height = unit(0.8, "cm")  # Add this to control key height
  )



ggsave("PSeverebar_U5.tiff", plot = PSeverebar_U5, width=12, height=6, dpi = 600)

#################################################
# Here we look at the total deaths. 
#################################################
Pdeaths_total <- ggplot(cases_long %>% filter(metric == "deaths_total"),
                        aes(x = year, y = value,color = scenario)) +
  geom_vline(xintercept = 2034, linetype = "dotted", color = "red", size = 1) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = scenario_colors,
    labels = c(
      "Scenario1" = "80% vaccination +  80% SMC",
      "Scenario2" = "80% vaccination + 0% SMC",
      "Scenario3" = "0% vaccination + 80% SMC",
      "Scenario4" = "50% vaccination + 50% SMC"
    ),
    guide = guide_legend(ncol = 2)
  ) +
  scale_x_continuous(
    breaks = seq(min(cases_long$year), max(cases_long$year)),
    labels = as.character(seq(min(cases_long$year), max(cases_long$year)))
  ) +
  labs(
    x = "Years",
    y = "Annual deaths in total population",
    color = "Scenarios"
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    labels = scales::label_comma()) +
  theme_bw(base_size = 14) +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.key.width = unit(1.5, "cm"),
    legend.spacing.x = unit(0.5, "cm")
  )

ggsave("Pdeaths_total.tiff", plot = Pdeaths_total, width=12, height=6, dpi = 600)

########################################
# Bar plot for the Total deaths  
#######################################
cases_summary <- cases_long %>%
  filter(metric == "deaths_total", year >= 2026, year <= 2034) %>%
  group_by(scenario) %>%
  summarise(
    total_deaths = sum(value),
    .groups = "drop"
  )

# Calculate percent difference from Baseline
baseline_cases <- cases_summary %>% filter(scenario == "Baseline") %>% pull(total_deaths)

cases_summary <- cases_summary %>%
  mutate(
    diff_from_baseline = (total_deaths - baseline_cases) / baseline_cases * 100,
    label_text = ifelse(
      scenario != "Baseline",
      sprintf("Δ = %+0.1f", diff_from_baseline),
      ""
    )
  )

# Set scenario factor levels for consistent plotting

# Bar plot code
Pdeathsbar_total <- ggplot(cases_summary, aes(x = scenario, y = total_deaths, fill = scenario)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", aes(alpha = scenario != "Baseline")) +
  geom_text(aes(label = label_text), vjust = -0.8, size = 4.2, fontface = "bold") +
  scale_fill_manual(
    values = scenario_colors,
    name = "Scenarios",
    labels = c(
      "Baseline"   = "Baseline",
      "Scenario1"  = "Scenario1:\n80% vaccination +\n80% SMC",  # Added \n for line breaks
      "Scenario2"  = "Scenario2:\n80% vaccination +\n0% SMC",
      "Scenario3"  = "Scenario3:\n0% vaccination +\n80% SMC",
      "Scenario4"  = "Scenario4:\n50% vaccination +\n50% SMC"
    )
  ) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5), guide = "none") +
  labs(
    subtitle = "Δ = Difference vs Baseline",
    x = "Scenario",
    y = "Total deaths in the overall population"
  ) +
  scale_y_continuous(
    labels = label_comma(),
    expand = expansion(mult = c(0, 0.1))  # Adds 10% space above top bar
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.8, "cm"),
    legend.spacing.x = unit(0.5, "cm"),
    plot.caption = element_text(hjust = 0, size = 10, face = "italic"),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(face = "italic", size = 12),
    legend.key.height = unit(0.8, "cm")  # Add this to control key height
  )

ggsave("Pdeathsbar_total.tiff", plot = Pdeathsbar_total, width=12, height=6, dpi = 600)


#################################################
# Here we look at under five malaria deaths. 
#################################################
Pdeaths_U5 <- ggplot(cases_long %>% filter(metric == "deaths2"),
                     aes(x = year, y = value,color = scenario)) +
  geom_vline(xintercept = 2034, linetype = "dotted", color = "red", size = 1) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = scenario_colors,
    labels = c(
      "Scenario1" = "80% vaccination +  80% SMC",
      "Scenario2" = "80% vaccination + 0% SMC",
      "Scenario3" = "0% vaccination + 80% SMC",
      "Scenario4" = "50% vaccination + 50% SMC"
    ),
    guide = guide_legend(ncol = 2)
  ) +
  scale_x_continuous(
    breaks = seq(min(cases_long$year), max(cases_long$year)),
    labels = as.character(seq(min(cases_long$year), max(cases_long$year)))
  ) +
  labs(
    x = "Years",
    y = "Annual deaths in children under 5",
    color = "Scenarios"
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    labels = scales::label_comma()) +
  theme_bw(base_size = 14) +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.key.width = unit(1.5, "cm"),
    legend.spacing.x = unit(0.5, "cm")
  )

ggsave("Pdeaths_U5.tiff", plot = Pdeaths_U5, width=12, height=6, dpi = 600)

########################################
# Bar plot for under five malaria deaths 
#######################################
cases_summary <- cases_long %>%
  filter(metric == "deathsU5", year >= 2026, year <= 2034) %>%
  group_by(scenario) %>%
  summarise(
    total_u5deaths = sum(value),
    .groups = "drop"
  )

# Calculate percent difference from Baseline
baseline_cases <- cases_summary %>% filter(scenario == "Baseline") %>% pull(total_u5deaths)

cases_summary <- cases_summary %>%
  mutate(
    diff_from_baseline = (total_u5deaths - baseline_cases) / baseline_cases * 100,
    label_text = ifelse(
      scenario != "Baseline",
      sprintf("Δ = %+0.1f", diff_from_baseline),
      ""
    )
  )

# Set scenario factor levels for consistent plotting

# Bar plot code
Pdeathsbar_U5 <- ggplot(cases_summary, aes(x = scenario, y = total_u5deaths, fill = scenario)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", aes(alpha = scenario != "Baseline")) +
  geom_text(aes(label = label_text), vjust = -0.8, size = 4.2, fontface = "bold") +
  scale_fill_manual(
    values = scenario_colors,
    name = "Scenarios",
    labels = c(
      "Baseline"   = "Baseline",
      "Scenario1"  = "Scenario1:\n80% vaccination +\n80% SMC",  # Added \n for line breaks
      "Scenario2"  = "Scenario2:\n80% vaccination +\n0% SMC",
      "Scenario3"  = "Scenario3:\n0% vaccination +\n80% SMC",
      "Scenario4"  = "Scenario4:\n50% vaccination +\n50% SMC"
    )
  ) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5), guide = "none") +
  labs(
    subtitle = "Δ = Difference vs Baseline",
    x = "Scenario",
    y = "Total deaths in children under 5"
  ) +
  scale_y_continuous(
    labels = label_comma(),
    expand = expansion(mult = c(0, 0.1))  # Adds 10% space above top bar
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.8, "cm"),
    legend.spacing.x = unit(0.5, "cm"),
    plot.caption = element_text(hjust = 0, size = 10, face = "italic"),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(face = "italic", size = 12),
    legend.key.height = unit(0.8, "cm")  # Add this to control key height
  )


ggsave("Pdeathsbar_U5.tiff", plot = Pdeathsbar_U5, width=12, height=6, dpi = 600)




















# Let us look at sensitivity analysis
n_samples <- 2

# runs 50


parms_sense <- read_excel("Sensitivity Analysis.xlsx", sheet = "Sheet1") # or the actual sheet name

# OPTIONAL: Ensure only relevant columns are kept and named correctly
parms_sense <- parms_sense[, c("Name", "Distribution", "mean", "low", "high")]

# Function to draw samples for a single parameter
sense_samples <- function(dist, mn, lo, hi, n) {
  if (dist == "tri") {
    return(rtriangle(n, a = lo, b = hi, c = mn))
  } else if (dist == "trirate") {
    days <- rtriangle(n, a = lo, b = hi, c = mn)
    return(1/days)
  } else if (dist == "unif") {
    return(runif(n, min = lo, max = hi))
  } else if (dist == "unifrate") {
    days <- runif(n, min = lo, max = hi)
    return(1/days)
  } else if (dist == "const") {
    return(rep(mn, n))
  } else {
    stop(paste("Unknown distribution:", dist))
  }
}

sense_list <- lapply(seq_len(nrow(parms_sense)), function(i) {
  with(parms_dist[i,], sense_samples (Distribution, mean, low, high, n_samples))
})

# Combine into a data.frame with columns for each parameter
sense_df <- as.data.frame(do.call(cbind, sense_list))
names(sense_df) <- parms_sense$Name
sense_df$sense_id  <- 1:n_samples

# Preview the sampled parameter sets
head(sense_df)

data_Fsen <- NULL
n_runs <- nrow(sense_df)

for (i in 1:n_runs) {
  # Convert the i-th row of sense_df to a named vector or list
  parms <- as.list(sense_df[i, names(sense_df) != "sense_id"])
  
  # Run your model
  Let_go <- deSolve::ode(times = times, y = state, func = Marks_go, parms = parms)
  
  # Collect total incidence (or other output, e.g. last time point, column 7)
  TInc = get_out[length(get_out[,"CInc"]),"CInc"]  # adapt column if needed
  
  # Store all parameter values and TInc for this run
  data_Fsen <- rbind(data_Fsen, as.data.frame(as.list(c(parms, TInc))))
  
}

lm1 <- lm(CInc ~ ., data = data_Fsen[, c(which(names(data_Fsen) == "a") : which(names(data_Fsen) == "kappa1"), which(names(data_Fsen) == "CInc"))])
library(QuantPsyc)
lm2 <- lm.beta(lm1)
stdcoeff <- cbind(sort(abs(lm2), decreasing = TRUE))



##########################################################
# Calculate Spearman correlations for each parameter with CInc
##########################################################

cors <- sapply(data_Fsen[, param_names], function(x) {
  if (sd(x, na.rm = TRUE) == 0) return(NA)  # skip if no variation
  cor(x, data_Fsen$CInc, method = "spearman", use = "complete.obs")
})

# Remove parameters with NA correlations
cors_clean <- cors[!is.na(cors)]


cor_df <- data.frame(
  Parameter = names(cors_clean),
  Correlation = cors_clean
)


ggplot(cor_df, aes(x = "", y = reorder(Parameter, Correlation), fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "purple", mid = "white", high = "red", midpoint = 0) +
  labs(
    title = "Heatmap of Parameter Correlations with CInc",
    x = NULL,
    y = "Parameter"
  ) +
  theme_minimal()






####################################################
# Sensitivity Analysis for the Vaccine 
###################################################







# 1. Define sensitivity grid
# Define grid
cov1_vals  <- seq(0.5, 0.9, by = 0.1)
paV_vals   <- seq(0.5, 0.9, by = 0.1)
rho_V_vals <- c(1/(5*365), 1/(4*365), 1/(3*365))
sens_grid <- expand.grid(cov1 = cov1_vals, paV = paV_vals, rho_V = rho_V_vals)

# Sensitivity loop
annual_metrics_list <- lapply(seq_len(nrow(sens_grid)), function(i) {
  sens_parms <- parms
  sens_parms["cov1"]  <- sens_grid$cov1[i]
  sens_parms["paV"]   <- sens_grid$paV[i]
  sens_parms["rho_V"] <- sens_grid$rho_V[i]
  
  out <- deSolve::ode(times = times, y = state, func = Marks_go, parms = sens_parms)
  df <- as_tibble(as.data.frame(out))
  
  # Calculate annual metrics (as in your code)
  metrics <- df %>%
    mutate(
      date = time + model_start_time,
      year = lubridate::year(date),
      P0 = S0 + E0 + A0 + C0 + Ct0 + Sev0 + T0 + H0 + R0,
      P1 = S1 + E1 + A1 + C1 + Ct1 + Sev1 + T1 + H1 + R1 + V + E1V + A1V + C1V + Ct1V + Sev1V + T1V + H1V + R1V+ Smcs+ Smcr + VSmcs + VSmcr,
      P2 = S2 + E2 + A2 + C2 + Ct2 + Sev2 + T2 + H2 + R2 + F2,
      P  = P0 + P1 + P2,
      Inc0 = c(0, diff(CInc0)),
      Inc1 = c(0, diff(CInc1)),
      Inc2 = c(0, diff(CInc2)),
      Inc  = c(0, diff(CInc))
    ) %>%
    group_by(year) %>%
    summarise(
      P0 = mean(P0),
      P1 = mean(P1),
      P2 = mean(P2),
      P  = mean(P),
      Inc0 = sum(Inc0),
      Inc1 = sum(Inc1),
      Inc2 = sum(Inc2),
      Inc  = sum(Inc)
    ) %>%
    mutate(
      incidence0 = (Inc0 / P0) * 1000,
      incidence1 = (Inc1 / P1) * 1000,
      incidence2 = (Inc2 / P2) * 1000,
      total_incidence = (Inc / P) * 1000,
      cov1 = sens_grid$cov1[i],
      paV = sens_grid$paV[i],
      rho_V = sens_grid$rho_V[i]
    )
  
  metrics
})

# Combine all runs
annual_metrics_all <- bind_rows(annual_metrics_list)

# Choose target year and metric
target_year <- 2034

heatmap_data <- annual_metrics_all %>%
  filter(year == target_year) %>%
  select(cov1, paV, rho_V, total_incidence) %>%
  mutate(
    cov1 = factor(cov1),
    paV = factor(paV),
    rho_V_year = factor(rho_V, labels = c("1/5 years", "1/4 years", "1/3 years"))
  )

# Plot heatmap of total annual incidence
ggplot(heatmap_data, aes(x = cov1, y = paV, fill = total_incidence)) +
  geom_tile() +
  facet_wrap(~ rho_V_year) +
  scale_fill_gradient2(
    low = "#a1d76a", mid = "#f7f7f7", high = "#e9a3c9",
    midpoint = median(heatmap_data$total_incidence, na.rm = TRUE),
    name = "Incidence per 1000"
  ) +
  labs(
    title = paste("Malaria Incidence per 1000 (Total Population), Year", target_year),
    x = "Vaccine coverage", y = "Vaccine effectiveness", fill = "Incidence per 1000"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    legend.text = element_text(size = 12)
  )


target_year <-2034

heatmap_data <- annual_metrics_all %>%
  filter(year == target_year) %>%
  select(cov1, paV, rho_V, total_incidence) %>%
  mutate(
    cov1 = factor(cov1),
    paV = factor(paV),
    rho_V_year = factor(rho_V, labels = c("1/5 years", "1/4 years", "1/3 years"))
  )

ggplot(heatmap_data, aes(x = cov1, y = paV, fill = total_incidence)) +
  geom_tile() +
  facet_wrap(~ rho_V_year) +
  scale_fill_gradient2(
    low = "#a1d76a", 
    mid = "#f7f7f7", 
    high = "#e9a3c9",
    midpoint = median(heatmap_data$total_incidence, na.rm = TRUE),
    name = "Incidence per 1000",
    breaks = c(-10, 0, 10, 20, 30),         # manual legend breaks
    labels = c("-10", "0", "10", "20", "30") # manual legend labels
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(6, "cm"),
      barheight = unit(0.5, "cm")
    )
  ) +
  labs(
    title = paste("Malaria Incidence per 1000 (Total Population), Year", target_year),
    x = "Vaccine coverage", y = "Vaccine effectiveness"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold")
  )


inc_breaks <- pretty(range(heatmap_data$total_incidence, na.rm = TRUE), n = 5)

ggplot(heatmap_data, aes(x = cov1, y = paV, fill = total_incidence)) +
  geom_tile() +
  facet_wrap(~ rho_V_year) +
  scale_fill_gradient2(
    low = "#a1d76a", mid = "#f7f7f7", high = "#e9a3c9",
    midpoint = median(heatmap_data$total_incidence, na.rm = TRUE),
    name = "Incidence per 1000",
    breaks = c(-10, 0, 10, 20, 30),
    labels = c("-10", "0", "10", "20", "30")
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(6, "cm"),
      barheight = unit(0.5, "cm")
    )
  ) +
  labs(
    x = "Vaccine coverage", 
    y = "Vaccine effectiveness"
    # Removed plot title
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # plot.title = element_text(face = "bold", size = 12), # <-- Removed
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold")
  )











# ---- Sensitivity Parameters ----
parameter_ranges <- list(
  cov1  = c(0, 1),
  scov  = c(0, 1),
  rho_V = c(1/(5*365), 1/(3*365)),
  a     = c(0.2, 0.5),
  mu_m  = c(1/20, 1/10),
  r     = c(1/8, 1/4),
  pa    = c(0.3, 0.6)
)

run_parameters <- list(
  cov1  = c(0, 0.25, 0.5, 0.75, 1),
  scov  = c(0, 0.25, 0.5, 0.75, 1),
  rho_V = c(1/(5*365), 1/(4*365), 1/(3*365)),
  a     = c(0.25, 0.33, 0.4, 0.5),
  mu_m  = c(1/20, 1/15, 1/10),
  r     = c(1/8, 1/6, 1/4),
  pa    = c(0.3, 0.44, 0.5, 0.6)
)

# ---- Model wrapper functions ----
solveMalariaModel <- function(init_state, changed_params = NULL, run_time = max(times)) {
  model_parms <- parms
  if (!is.null(changed_params)) {
    for (param in names(changed_params)) {
      model_parms[[param]] <- changed_params[[param]]
    }
  }
  out <- as.data.frame(ode(
    y = init_state,
    times = times,
    func = Marks_go,   # <-- your model function name
    parms = model_parms
  ))
  return(out)
}

getAnnualIncidence <- function(solution, target_year = 2025, start_date = model_start_time) {
  solution <- solution %>%
    mutate(date = start_date + time,
           year = year(date))
  sum(solution$CInc[solution$year == target_year], na.rm = TRUE)
}

latin_hypercube_sampling <- function(n_samples, parameter_ranges) {
  n_parameters <- length(parameter_ranges)
  samples_unit <- randomLHS(n_samples, n_parameters)
  samples <- samples_unit
  for (i in seq_along(parameter_ranges)) {
    min_val <- parameter_ranges[[i]][1]
    max_val <- parameter_ranges[[i]][2]
    samples[, i] <- min_val + (max_val - min_val) * samples_unit[, i]
  }
  colnames(samples) <- names(parameter_ranges)
  return(samples)
}

runLHSForFixAndRangesParallel <- function(init_state, n_samples, variable_ranges, fixed_param, fixed_value) {
  samples <- latin_hypercube_sampling(n_samples, variable_ranges)
  n_cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, varlist = c("init_state", "samples", "solveMalariaModel",
                                          "getAnnualIncidence", "parms", "Marks_go", "model_start_time", "times"),
                          envir = environment())
  outputs <- parallel::parLapply(cl, 1:n_samples, function(i) {
    changed_params <- as.list(samples[i, ])
    changed_params[[fixed_param]] <- fixed_value
    sol <- solveMalariaModel(init_state, changed_params)
    getAnnualIncidence(sol, target_year = 2025)
  })
  parallel::stopCluster(cl)
  return(unlist(outputs))
}

sensAnalysisGetLHS <- function(init_state, n_samples_each, runs, ranges) {
  results_list <- list()
  for (run in names(runs)) {
    for (fix in runs[[run]]) {
      fixed_name <- run
      fixed_val <- fix
      cat("Running: ", fixed_name, "=", fixed_val, "\n")
      # Special run: all parameters at middle of range, except fixed parameter
      default_params <- sapply(ranges, function(r) mean(r))
      default_params[[fixed_name]] <- fixed_val
      special_sol <- solveMalariaModel(init_state, default_params)
      special_result <- data.frame(
        outcome = getAnnualIncidence(special_sol, target_year = 2025),
        fixed_param = fixed_name,
        fixed_value = fix,
        type = "special"
      )
      # LHS varying other parameters
      variable_ranges <- ranges[names(ranges) != run]
      lhs_outcomes <- runLHSForFixAndRangesParallel(init_state, n_samples_each, variable_ranges, fixed_name, fixed_val)
      lhs_result <- data.frame(
        outcome = lhs_outcomes,
        fixed_param = fixed_name,
        fixed_value = fix,
        type = "lhs"
      )
      results_list[[length(results_list) + 1]] <- rbind(special_result, lhs_result)
    }
  }
  final_df <- do.call(rbind, results_list)
  return(final_df)
}

# ---- Run Sensitivity Analysis ----
n_samples_each <- 200 # Number of LHS samples per fixed value
final_df <- sensAnalysisGetLHS(state, n_samples_each, run_parameters, parameter_ranges)

# ---- Plot Results ----
ggplot(final_df, aes(x = factor(fixed_value), y = outcome)) +
  geom_boxplot(data = final_df[final_df$type == "lhs", ], outlier.shape = NA) +
  geom_point(
    data = final_df[final_df$type == "special", ],
    aes(x = factor(fixed_value), y = outcome),
    color = "red", size = 3, shape = 18
  ) +
  facet_wrap(~ fixed_param, scales = "free_x") +
  labs(
    x = "Fixed Parameter Value",
    y = "Annual Incidence (2025)",
    title = "Sensitivity Analysis on Malaria Model"
  ) +
  theme_minimal()

# ---- Save results ----
write.csv(final_df, "malaria_sensitivity_results.csv", row.names = FALSE)

# ---- Optionally: Density plot for a single parameter ----
lhs_subset <- subset(final_df, fixed_param == "cov1" & type == "lhs")
special_subset <- subset(final_df, fixed_param == "cov1" & type == "special")

ggplot(lhs_subset, aes(x = outcome)) +
  geom_density(fill = "skyblue", alpha = 0.6, bw = 0.5) +
  geom_vline(
    data = special_subset,
    aes(xintercept = outcome),
    color = "red", linetype = "dashed", size = 1
  ) +
  facet_wrap(~ fixed_value, ncol = 1, scales = "fixed") +
  labs(
    x = "Annual Incidence (2025)",
    y = "Density",
    title = "Outcome Densities for Different Vaccine Coverage Values (LHS)",
    subtitle = "Red dashed line = special run (only cov1 changed)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
