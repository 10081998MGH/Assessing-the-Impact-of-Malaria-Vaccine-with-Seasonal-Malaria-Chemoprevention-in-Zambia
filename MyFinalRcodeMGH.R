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



getwd()
# change to your respective directory (Remember to remove this before submission as it is not need )
setwd("")
getwd()

#Read in insecticide-treated nets (ITN) data from CSV
# Note 
# 2000-2012 Model warm-up period for stabilization. Back-predicted data from this period were not used for analysis.
# 2013-2023 Actual net data for Zambia, sourced from the World Malaria Report.
# 2024-2035 Predicted net data, projected based on historical trends and realistic budget constraints.

itndat<-read.csv("netsdata1.csv")
itn_cov <-itndat$nets # Proportion of the population covered by nets

# Read in the observed malaria incidence(per 1000) data from Excel
observed_data <- read_excel("Malaria_incidence data.xlsx")

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
    dSmcs <- srate * (S1) - (snap + mu_h) * Smcs              # children from S1 moving into SMC compartment (Smcs) 
    dSmcr <- srate * (A1 + R1) - (snap + mu_h) * Smcr         # children from A1 or R1 moving into SMC-protected compartment (Smcr)
    
    
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
    dVSmcs <-  srate * V - (snap + mu_h) * VSmcs               # vaccinated susceptible moving into SMC compartment
    dVSmcr <-  srate * (A1V + R1V) - (snap + mu_h) * VSmcr     # vaccinated recovered/asymptomatic moving into SMC compartment
    
    # Age group 2 (>= 60 months old) differential equations
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
    
    
    # Cumulative incidence counters (for tracking total infections over time)
    CInc0 <-  lambda * S0                       # new infections in 0-5 months
    CInc1 <-  lambda * (S1 + V)                 # new infections in 6-59 months 
    CInc2 <-  lambda * (S2 + R2)                # new infections in >=60 months  
    CIncU5 <- CInc0 + CInc1                     # new infections in under-five (0-59 months)
    CInc  <-  CInc0 + CInc1 + CInc2             # total new infections in all ages
    
    
    # SMC dosing counters (to track SMC actions, if needed)
    dSMC <-  srate*(S1+A1+R1+V+A1V+R1V)         # total individuals put on SMC
    dCSnap <-  snap*(S1+R1+V+R1V)               # total individuals leaving SMC protection 
    
    
    # Case counts (incidence of clinical cases) 
    Cases0 <-   gamma_h * (1 - pa)  * E0                                 # clinical cases in 0-5 months
    Cases1 <-   gamma_h * (1 - pa)  * E1 + gamma_h * (1 - paV)  * E1V    # clinical cases in 6-59 months
    Cases2 <-  gamma_h * (1 - pa)  * E2 + gamma_h * (1 - pas)  * F2      # clinical cases in >=60 months
    CasesU5 <-  Cases0 + Cases1                                          # clinical cases in under-fives
    Cases  <- Cases0 + Cases1 + Cases2                                   # clinical cases in total population
    
    # Severe case counts (cumulative, untreated turning severe)
    Csevere0  <-  nu * C0                             # Severe cases in 0-5 months             
    Csevere1  <-  nu * C1 + nu * C1V                  # severe cases in 6-59 months 
    Csevere2  <-  nu * C2                             # severe cases in >=60 months
    CsevereU5 <- Csevere0 + Csevere1                  # severe cases in under-fives
    Csevere   <-  Csevere0 + Csevere1 + Csevere2      # total severe cases
    
    
    # Death counts (from severe cases)
    deaths0 <- tau * mu_sev * Sev0                              # deaths in 0-5 months
    deaths1 <- tau * mu_sev * Sev1 + tau * mu_sev * Sev1V       # deaths in 6-59 months 
    deaths2 <- tau * mu_sev * Sev2                              # deaths in >=60 months
    deathsU5 <- deaths0 + deaths1                               # deaths in under-fives
    Tdeaths <- deaths0 + deaths1 + deaths2                      # Total deaths
    
    
    
    
    # Return the rate of change (derivatives) for all state variables
    output <- c(dS0, dE0, dA0, dC0, dCt0, dSev0, dT0, dH0, dR0, dV,
                dS1, dE1, dA1, dC1, dCt1, dSev1, dT1, dH1, dR1,dSmcs,dSmcr,
                dE1V,dA1V, dC1V, dCt1V, dSev1V, dT1V, dH1V, dR1V,dVSmcs,dVSmcr,
                dS2, dE2, dA2, dC2, dCt2, dSev2, dT2, dH2, dR2,
                dF2,CInc0,CInc1,CInc2, CIncU5, CInc, dSMC, dCSnap,
                Cases0, Cases1, Cases2, CasesU5, Cases, Csevere0, Csevere1, Csevere2, CsevereU5,Csevere,
                deaths0, deaths1, deaths2, deathsU5, Tdeaths 
    )
    list(output,lambda= lambda, Seas= Seas,itn=itn)#, lambda= lambda, itn=itn, irs=irs)
  })
}

# MODEL PARAMETERS AND INITIAL CONDITIONS

# An apporximate population for Zambia following the warm period in 2000

# AgeGroup0 = 408,458
# 
# AgeGroup1 = 1,628,833
# 
# AgeGroup2 = 9,058,710
# 
# Total Population = 11,096,001

# Initial state vector - all compartments start with approximate realistic values
state <- c(S0 = 211495, E0 = 6348, A0 = 17047, C0 = 1333 , Ct0 =444, Sev0 = 0, T0 = 0, H0 = 0, R0 =170791, V=0,
           S1 = 845979, E1 = 25392, A1 = 68189, C1 = 5332, Ct1 =1776, Sev1 = 0, T1 = 0, H1 = 0, R1 = 683164,Smcs=0,Smcr=0,
           E1V = 0, A1V = 0, C1V = 0, Ct1V = 0, Sev1V = 0, T1V = 0, H1V = 0, R1V = 0,VSmcs=0,VSmcr=0,
           S2 = 4652177, E2 =139636 , A2 = 374983, C2 =29321 , Ct2 =9764, Sev2 = 0, T2 = 0, H2 = 0, R2 =3756829, 
           F2 = 0, CInc0 = 0,CInc1 = 0,CInc2 = 0,  CIncU5=0, CInc = 0, SMC=0, CSnap = 0, cases0 = 0, cases1 = 0, cases2 = 0, casesU5 = 0, Cases = 0,
           Csevere0 = 0, Csevere1 = 0, Csevere2 = 0, CsevereU5 = 0,Csevere = 0, deaths0 = 0, deaths1 = 0, deaths2 = 0, deathsU5 = 0, Tdeaths = 0
)

#Model Parameters
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
  paV = 0.77,                                 #Vaccine Efficacy 
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

# Simulation time settings
model_start_time <- ymd("2000-01-01") # Start date of the model

# Calculate total duration in days (warm-up + simulation)
warmup_years <- 8                        # years of warm up (just reach equilibrium)
simulation_years <- 27                   # years of simulation after warm up
total_days <- (warmup_years + simulation_years) * 365 # total days to simulate

times <- seq(0, total_days, by = 1) # time steps in days from 0 to total_days


# Run the ODE model simulation
# Use deSolve to integrate the ODEs over time
Let_go <- deSolve::ode(times = times, y = state, func = Marks_go, parms = parms)







##########################################################################
# Model calibration and Uncertainty intervals in the model 
##########################################################################


# This script performs uncertainty analysis by sampling from parameter distributions
# and running multiple model simulations to assess prediction confidence intervals

# Note: this takes some time to run, depending on the number of samples.(200 samples = 6-8hours)
n_samples <- 200 # Set the number of samples to draw from each parameter distribution

# Load parameter distributions from Excel file containing model parameters
# The "universal" sheet contains parameter names, distribution types, and bounds
parms_dist <- read_excel("ParametersSetwork.xlsx", sheet = "universal") # or the actual sheet name

#Extract only the essential columns needed for parameter sampling
# Name: parameter identifier, Distribution: type (tri/unif/const), mean/low/high: distribution parameters
parms_dist <- parms_dist[, c("Name", "Distribution", "mean", "low", "high")]

# Function to sample from a given distribution type
# dist: distribution type (e.g., tri, unif, const)
# mn: mean or mode value
# lo, hi: lower and upper bounds
# n: number of samples to draw
draw_samples <- function(dist, mn, lo, hi, n) {
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
  } else if (dist == "constrate") {
    return(rep(1 / mn, n))  
  } else {
    stop(paste("Unknown distribution:", dist))
  }
}

# Apply sampling function to each row in the parameter table
sampled_list <- lapply(seq_len(nrow(parms_dist)), function(i) {
  with(parms_dist[i,], draw_samples(Distribution, mean, low, high, n_samples))
})

# Combine sampled parameters into a data frame
samples_df <- as.data.frame(do.call(cbind, sampled_list))
names(samples_df) <- parms_dist$Name
samples_df$sample_id <- 1:n_samples

# Preview the sampled parameter sets

head(samples_df)

# Sample_imp <- saveRDS(samples_df, "samples_df.rds") # Save the sampled parameters for later use



# Convert sampled parameters into a list format
# Each element in the list corresponds to one parameter set
Param_list <- split(samples_df, samples_df$sample_id) %>%
  purrr::map(~ setNames(as.numeric(.x), names(.x)))

param_try <- function(model,Param_l) {
  deSolve::ode(times = times, y = state, func = Marks_go, parms = Param_l)
}

# Define function to run ODE model using one parameter set
output_list <- purrr::map(Param_list, ~ param_try(model, Param_l = .x))

# saveRDS(output_list, "output_list.rds") # Save the output list for later use
# results3 <- readRDS("output_list.rds") # Load the output list

#Combine all simulation results into a single data frame for analysis
# .id = "sample_id" preserves which parameter set generated each result
Combined_daf1 <-  output_list %>%
  purrr::map_df(~ as.data.frame(.x), .id = "sample_id") %>%
  bind_rows()

# saveRDS(Combined_df, "Combined_df.rds") # Save the combined data frame for later use
# load34 <- readRDS("Combined_df.rds")
# 
# Combined_df<- readRDS("Combined_df_0.rds")


#####################################################
## Data Processing and Visualization of Results
#####################################################
# Calculate annual metrics and population totals for each age group and simulation run
annual_metrics_combined <- as_tibble(as.data.frame(Combined_daf1)) %>%
  mutate(
    date = time + model_start_time,
    year = year(date),
    
    P0 = S0 + E0 + A0 + C0 + Ct0 + Sev0 + T0 + H0 + R0,
    P1 = S1 + E1 + A1 + C1 + Ct1 + Sev1 + T1 +  H1 + R1 + V + E1V + A1V + C1V + Ct1V + Sev1V + T1V + H1V + R1V+ Smcs+ Smcr + VSmcs + VSmcr,
    P2 = S2 + E2 + A2 + C2 + Ct2 + Sev2 + T2 + H2 + R2 + F2,
    P  = P0 + P1 + P2,
    
    Inc0 = c(0, diff(CInc0)),
    Inc1 = c(0, diff(CInc1)),
    Inc2 = c(0, diff(CInc2)),
    Inc  = c(0, diff(CInc))
  ) %>%
  group_by(year,sample_id) %>%
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
    total_incidence = (Inc / P) * 1000
  ) %>%
  # Step 2: Pivot to long format
  pivot_longer(
    cols = c(incidence0, incidence1, incidence2, total_incidence),
    names_to = "population_group",
    values_to = "incidence_per_1000"
  ) %>%
  mutate(
    population_group = recode(population_group,
                              "incidence0" = "P0",
                              "incidence1" = "P1",
                              "incidence2" = "P2",
                              "total_incidence" = "P"),
    population_group = factor(population_group, levels = c("P", "P0", "P1", "P2"))
  )

# Create custom labels for facet panels to improve plot readability
group_labels <- labeller(
  population_group = c(
    "P" = "Total Population",
    "P0" = "Population Age 0–5 months",
    "P1" = "Population Age 6–59 months",
    "P2" = "Population Age > 60 months"
  )
)


# Create focused plot for total population comparing model predictions to observed data
malaria_sample <- annual_metrics_combined %>%
  filter(population_group == "P") %>%
  ggplot(aes(x = year, y = incidence_per_1000, color = sample_id)) +
  geom_line(size = 1) +
  geom_point(data = observed_data, 
             aes(x = year, y = observed_incidence,), 
             color = "#D55E00", size = 2, inherit.aes = FALSE) +
  coord_cartesian(ylim = c(0, 900)) + 
  theme_bw(base_size = 14) +
  labs(
    title = "Annual Malaria Incidence per 1000 (Total Population)",
    x = "Year",
    y = "Incidence per 1000"
  )
# Save high-resolution plot for thesis inclusion
# ggsave("malaria_incidence_total.tiff", plot = malaria_sample, width = 12, height = 6, dpi = 600)

# Prepare data for model fitting by matching simulation results with observed data
# Only include years where both simulated and observed data are available
annual_metrics_combined_fit <- annual_metrics_combined %>% 
  filter(population_group == "P") %>% 
  dplyr::select(year, incidence_per_1000, sample_id) %>% 
  inner_join(observed_data, by = "year")


# Calculate log-likelihood for each simulation run using Poisson distribution
# This quantifies how well each parameter set reproduces the observed data
loglikelihoods <- annual_metrics_combined_fit %>%  # Poisson log-likelihood: assumes observed incidence follows Poisson distribution
  mutate(loglik = dpois(round(observed_incidence, 0), lambda = incidence_per_1000, log = TRUE)) %>%  # Round observed values since Poisson requires integer counts
  group_by(sample_id) %>%
  summarise(loglik = sum(loglik, na.rm = TRUE), .groups = "drop")


# Select best-fitting parameter sets based on log-likelihood criteria
# Only retain runs with reasonable fit to observed data
bestruns <- loglikelihoods %>%  
  arrange(desc(loglik)) %>% 
  filter(loglik > max(loglik)*5)   # Filter criterion: keep runs with log-likelihood > 5 * maximum log-likelihood
# This removes poorly fitting parameter sets while retaining reasonable fits
# slice_max(loglik, n = floor(n_samples*0.7))# Get the top 40% runs with the highest log likelihoods


# View(bestruns)
# hist(bestruns$loglik, breaks = 100) # Histogram of log likelihoods
# saveRDS(bestruns, "bestruns.rds") # Save the best runs
# Read the best runs back in
# we now collect the best runs for our analysis 


# Filter the annual metrics to include only the best-fitting parameter sets
# This reduces uncertainty by excluding parameter sets that poorly match observed data
best_metrics_combined_fit <- annual_metrics_combined_fit %>% 
  filter(sample_id %in% bestruns$sample_id)


# Create diagnostic plot comparing best-fitting simulations to observed data
# Shows individual simulation trajectories along with observed data points
best_metrics_combined_fit %>%
  ggplot(aes(x = year)) +
  geom_line(aes( y = incidence_per_1000, color = sample_id)) +
  geom_point(aes(y = observed_incidence), 
             color = "#D55E00", size = 2) +
  coord_cartesian(ylim = c(0, 900)) + 
  theme_bw(base_size = 14) +
  labs(
    title = "Annual Malaria Incidence per 1000 (Total Population)",
    x = "Year",
    y = "Incidence per 1000"
  )

# Calculate uncertainty bounds (confidence intervals) from selected simulation runs
# Uses quantiles to define prediction intervals around the median
selected_vars_ci <-  best_metrics_combined_fit %>% 
  group_by(year) %>% 
  summarise(upper95 = quantile(incidence_per_1000, 0.95),
            lower5 = quantile(incidence_per_1000 , 0.05),
            median = median(incidence_per_1000))

# Create uncertainty visualization showing prediction intervals
# This plot shows model uncertainty as a shaded region around the median prediction
best_metrics_combined_fit %>%
  inner_join(selected_vars_ci, by = "year") %>% 
  ggplot(aes(x = year)) +
  #geom_line(aes( y = incidence_per_1000, color = sample_id)) +
  geom_point(aes(y = observed_incidence), 
             color = "#D55E00", size = 2) +
  geom_ribbon(aes(ymin = lower5, ymax = upper95), 
              fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = median), color = "black", show.legend = F) +
  coord_cartesian(ylim = c(0, 900)) + 
  theme_bw(base_size = 14) +
  labs(
    title = "Annual Malaria Incidence per 1000 (Total Population)",
    x = "Year",
    y = "Incidence per 1000"
  )


# Model prediction data with source label for legend
model_data <- best_metrics_combined_fit %>%
  inner_join(selected_vars_ci, by = "year") %>%
  mutate(Source = "Model Prediction")

# Observed data with source label for legend 
obs_data <- observed_data %>%
  mutate(Source = "Observed Malaria Cases")

# Create final publication-ready plot with proper legend and the observed data
P1 <- ggplot() +
  # Ribbon for model uncertainty (no legend)
  geom_ribbon(data = model_data, 
              aes(x = as.factor(year), ymin = lower5, ymax = upper95, group = 1),
              fill = "lightblue", alpha = 0.3, inherit.aes = FALSE, show.legend = FALSE) +
  
  # Median modelled line
  geom_line(data = model_data,
            aes(x = as.factor(year), y = median, color = Source, group = 1),
            size = 1.2) +
  
  # Observed data: dashed red line + points
  geom_line(data = obs_data,
            aes(x = as.factor(year), y = observed_incidence, color = Source, group = 1),
            linetype = "dashed", size = 0.9) +
  geom_point(data = obs_data,
             aes(x = as.factor(year), y = observed_incidence, color = Source),
             size = 2) +
  
  # Define color mapping
  scale_color_manual(
    values = c("Model Prediction" = "black", "Observed Malaria Cases" = "#D55E00")
  ) +
  
  # Aesthetics
  coord_cartesian(ylim = c(0, 900)) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    title = "Modelled and Observed Annual Malaria Incidence per 1000 (Total Population)",
    x = "Year",
    y = "Incidence per 1000",
    color = ""
  )

# ggsave("Modelled and Observed.tiff", plot = P1, width = 12, height = 6, dpi = 600)




##########################################################################
# Scenario Analysis : Vaccine and SMC intervention
##########################################################################
# Define baseline and alternative scenarios for interventions

base_params <- parms
scenario_param_names <- c("itn_use", "cov1", "scov")

#Scenario definitions (Only varry SMC and Vaccination coverage)
scenario_analysis <- list(
  Baseline  = c(itn_use = base_params["itn_use"], cov1 = 0, scov = 0),
  Scenario1 = c(itn_use = base_params["itn_use"], cov1 =0.8, scov = 0.8),
  Scenario2 = c(itn_use = base_params["itn_use"], cov1 = 0.8, scov = 0),
  Scenario3 = c(itn_use = base_params["itn_use"], cov1 = 0.0, scov = 0.8),
  Scenario4 = c(itn_use = base_params["itn_use"], cov1 = 0.5, scov = 0.5)
  
)

# Scenario runner
# Function to run a scenario by updating parameters and solving the model

run_scenario <- function(scenario_values, base_params, times, state) {
  params <- base_params
  params[scenario_param_names] <- scenario_values
  result <- ode(times = times, y = state, func = Marks_go, parms = params)
  as.data.frame(result)
}

# Run all scenarios 
# Run all scenarios and combine results
scenario_outputs <- bind_rows(
  lapply(names(scenario_analysis), function(name) {
    out <- run_scenario(scenario_analysis[[name]], base_params, times, state)
    out$scenario <- name
    out
  })
)


# Postprocess results
# Process scenario outputs to calculate annual incidence per 1000 for each scenario
# We use the cases, severe cases and deaths for our write up. All the plots can be found in the suplementary folder

scenario_df <- scenario_outputs %>%
  mutate(
    date = time + model_start_time,
    year = lubridate::year(date),
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

# Compute median and 5th/95th percentiles
ci_summary <- scenario_df %>%
  group_by(year, scenario, population_group) %>%
  summarise(
    upper95 = quantile(incidence_per_1000, 0.95, na.rm=TRUE),
    lower5  = quantile(incidence_per_1000, 0.05, na.rm=TRUE),
    median  = median(incidence_per_1000, na.rm=TRUE),
    .groups = "drop"
  )

# Define colors for scenarios in plots
scenario_colors <- c(
  Baseline  = "gray30",
  Scenario1 = "#19c2c2",
  Scenario2 = "#1e90ff",
  Scenario3 = "#a259bc",
  Scenario4 = "#f4a259"
)

# Plot Baseline vs each Scenario for total population incidence

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

# Arrange the scenario comparison plots in a grid (2 columns)
do.call(grid.arrange, c(plots, ncol = 2))


# Combined plot of all scenarios including Baseline for total population
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




# Bar chart comparing total incidence (2026–2034) across scenarios

ci_summary_filtered <- ci_summary %>% 
  filter(population_group == "P", year >= 2026 & year <= 2034)

# Sum incidence over 2026-2034 for each scenario
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

# Set scenario order for consistent plotting
total_summary$scenario <- factor(total_summary$scenario, levels = c("Baseline", "Scenario1", "Scenario2", "Scenario3", "Scenario4"))

# Calculate differences from Baseline
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

# Plot bar chart and percentage difference labels
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


# Repeat similar plots for under five population (PU5)
# Line plots Baseline vs scenarios for U5 incidence
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

# Combine plots in a grid
do.call(grid.arrange, c(plots, ncol = 2))



# Combined line plot for Under 5 incidence, all scenarios
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


# Bar chart for total under-five incidence (2026-2034) across scenarios
ci_summary_filtered <- ci_summary %>% 
  filter(population_group == "PU5", year >= 2026 & year <= 2034)

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

# Set scenario order for consistent plotting
total_summary$scenario <- factor(total_summary$scenario, levels = c("Baseline", "Scenario1", "Scenario2", "Scenario3", "Scenario4"))

# Calculate differences from Baseline
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




# These are the scenario results that we used in our report for analysi(clinacal cases, severe cases and deaths)
##############################################################################
# Scenario Outcomes: Cases, Severe cases, Deaths
##############################################################################

scenario_Cases <- scenario_outputs %>%
  mutate(
    date = time + model_start_time,
    year = lubridate::year(date),
    # New cases (incidence)
    cases0 = c(0, diff(cases0)),cases1 = c(0, diff(cases1)),cases2 = c(0, diff(cases2)),
    casesU5 = c(0, diff(casesU5)),cases_total = c(0, diff(Cases)),
    # New severe cases
    severe0 = c(0, diff(Csevere0)),severe1 = c(0, diff(Csevere1)),severe2 = c(0, diff(Csevere2)),
    severeU5 = c(0, diff(CsevereU5)),severe_total = c(0, diff(Csevere)),
    # New deaths
    deaths0 = c(0, diff(deaths0)),deaths1 = c(0, diff(deaths1)),deaths2 = c(0, diff(deaths2)),
    deathsU5 = c(0, diff(deathsU5)), deaths_total = c(0, diff(Tdeaths))
  ) %>%
  group_by(year, scenario) %>%
  summarise(
    cases0 = sum(cases0),cases1 = sum(cases1),cases2 = sum(cases2),
    casesU5 = sum(casesU5),cases_total = sum(cases_total),
    severe0 = sum(severe0), severe1 = sum(severe1),severe2 = sum(severe2),
    severeU5 = sum(severeU5),severe_total = sum(severe_total),
    deaths0 = sum(deaths0), deaths1 = sum(deaths1),deaths2 = sum(deaths2),
    deathsU5 = sum(deathsU5),deaths_total = sum(deaths_total),
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

#################################################
#  Plot annual clinical cases in total population for each scenario
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
# Save image for analysis 
# ggsave("Pcases_total1.tiff", plot = Pcases_total1, width=12, height=6, dpi = 600)



# Bar chart for total clinical cases (2026-2034) by scenario (overall population)

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

Pcases_bartotal1 <- ggplot(cases_summary, aes(x = scenario, y = total_cases, fill = scenario)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", aes(alpha = scenario != "Baseline")) +
  geom_text(aes(label = label_text), vjust = -0.8, size = 4.2, fontface = "bold") +
  scale_fill_manual(
    values = scenario_colors,
    name = "Scenarios",
    labels = c(
      "Baseline"   = "Baseline",
      "Scenario1"  = "Scenario1:\n80% vaccination +\n80% SMC",  
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
    legend.key.height = unit(0.8, "cm")  
  )

# Save image for analysis 
# ggsave("Pcases_bartotal1.tiff", plot =Pcases_bartotal1, width=12, height=6, dpi = 600)




###############################################
#  Plot annual clinical cases in under-5 population for each scenario
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
# Save image for analysis 
# ggsave("Pcases_U5.tiff", plot = Pcases_U5, width=12, height=6, dpi = 600)



# Bar chart for total under-5 clinical cases (2026-2034) by scenario

cases_summary <- cases_long %>%
  filter(metric == "casesU5", year >= 2026, year <= 2034) %>%
  group_by(scenario) %>%
  summarise(
    total_u5cases = sum(value),
    .groups = "drop"
  )

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

Pcases_barU5 <- ggplot(cases_summary, aes(x = scenario, y = total_u5cases, fill = scenario)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", aes(alpha = scenario != "Baseline")) +
  geom_text(aes(label = label_text), vjust = -0.8, size = 4.2, fontface = "bold") +
  scale_fill_manual(
    values = scenario_colors,
    name = "Scenarios",
    labels = c(
      "Baseline"   = "Baseline",
      "Scenario1"  = "Scenario1:\n80% vaccination +\n80% SMC",  
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
    expand = expansion(mult = c(0, 0.1))  # 
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
    legend.key.height = unit(0.8, "cm")  
  )

# Save image for analysis 
# ggsave("Pcasesbar_U5.tiff", plot = Pcases_barU5, width=12, height=6, dpi = 600)



##############################################################
# Malaria severe cases plots for both under five and total population 
######################b#######################################

# Plot annual severe cases in total population for each scenario
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

# Bar chart for total severe cases (2026-2034) by scenario (total pop)

cases_summary <- cases_long %>%
  filter(metric == "severe_total", year >= 2026, year <= 2034) %>%
  group_by(scenario) %>%
  summarise(
    total_severe_cases = sum(value),
    .groups = "drop"
  )

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

PSeverebar_total <-  ggplot(cases_summary, aes(x = scenario, y = total_severe_cases, fill = scenario)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", aes(alpha = scenario != "Baseline")) +
  geom_text(aes(label = label_text), vjust = -0.8, size = 4.2, fontface = "bold") +
  scale_fill_manual(
    values = scenario_colors,
    name = "Scenarios",
    labels = c(
      "Baseline"   = "Baseline",
      "Scenario1"  = "Scenario1:\n80% vaccination +\n80% SMC",  
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
    expand = expansion(mult = c(0, 0.1)) 
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
    legend.key.height = unit(0.8, "cm")  
  )

# Save image for analysis 
# ggsave("PSeverebar_total.tiff", plot = PSeverebar_total, width=12, height=6, dpi = 600)


###########################################################
# Severe cases in under-5 population for each scenario
##########################################################

# Plot annual severe cases in under-5 population for each scenario
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
# Save image for analysis 
# ggsave("PSevere_U5.tiff", plot = PSevere_U5, width=12, height=6, dpi = 600)


# Bar chart for total under-5 severe cases (2026-2034) by scenario

cases_summary <- cases_long %>%
  filter(metric == "severeU5", year >= 2026, year <= 2034) %>%
  group_by(scenario) %>%
  summarise(
    total_u5_sevcases = sum(value),
    .groups = "drop"
  )

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


PSeverebar_U5 <- ggplot(cases_summary, aes(x = scenario, y = total_u5_sevcases, fill = scenario)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", aes(alpha = scenario != "Baseline")) +
  geom_text(aes(label = label_text), vjust = -0.8, size = 4.2, fontface = "bold") +
  scale_fill_manual(
    values = scenario_colors,
    name = "Scenarios",
    labels = c(
      "Baseline"   = "Baseline",
      "Scenario1"  = "Scenario1:\n80% vaccination +\n80% SMC",  
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


# Save image for analysis 
# ggsave("PSeverebar_U5.tiff", plot = PSeverebar_U5, width=12, height=6, dpi = 600)



#################################################
# Malaria death in both under five and Total population 
#################################################


# Plot annual deaths in total population for each scenario

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

# Save image for analysis 
# ggsave("Pdeaths_total.tiff", plot = Pdeaths_total, width=12, height=6, dpi = 600)


# Bar chart for total deaths (2026-2034) by scenario (total population)

cases_summary <- cases_long %>%
  filter(metric == "deaths_total", year >= 2026, year <= 2034) %>%
  group_by(scenario) %>%
  summarise(
    total_deaths = sum(value),
    .groups = "drop"
  )

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

# Save image for analysis 
# ggsave("Pdeathsbar_total.tiff", plot = Pdeathsbar_total, width=12, height=6, dpi = 600)


#################################################
# Under five malaria deaths. 
#################################################

# Plot annual deaths in under-5 population for each scenario
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
# Save image for analysis 
# ggsave("Pdeaths_U5.tiff", plot = Pdeaths_U5, width=12, height=6, dpi = 600)


# Bar chart for total under-5 deaths (2026-2034) by scenario

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

# Save image for analysis 
# ggsave("Pdeathsbar_U5.tiff", plot = Pdeathsbar_U5, width=12, height=6, dpi = 600)



