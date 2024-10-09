# Coupled photosynthesis, stomatal and energy balance model
# Steady-state solved with iterative procedure (Newton-Rhapson method)

# Author: Killian Dupont
# Date: 08/10/2024
# 
# This R script uses data from Vilà-Guerau de Arellano, J. et al. (2020). 
# CloudRoots: integration of advanced instrumental techniques and process modelling of sub-hourly and sub-kilometre land–atmosphere interactions. Biogeosciences, 17(17), 4375-4404
# https://doi.org/10.5194/bg-17-4375-2020
# Data taken from https://www.tr32db.uni-koeln.de/site/index.php (last accessed 08/10/2024), search for "cloudroots"
# Vegetation profiles of CO2, H2O, temperature and wind speed during the CloudRoots campaign 2018


library(tidyverse)
library(readr)

###############################
#Read in data
###############################
setwd("C:/R/new_phyto_review/Cloudroots")

root_dir <- "C:/R/new_phyto_review/Cloudroots"
file_paths <- list.files(root_dir, pattern = "prof.csv", full.names = TRUE, recursive = TRUE)

read_prof_file <- function(file_path) {
  path_parts <- strsplit(file_path, "/")[[1]]
  date_time <- path_parts[length(path_parts) - 1] 
  date <- substr(date_time, 1, 8)
  time <- substr(date_time, 9, 12)
  
  df <- read.csv(file_path)
  df <- df[-c(1, 2), ]
  df <- df %>%
    mutate(
      date = as.Date(date, format = "%Y%m%d"),
      time = sprintf("%02d:%02d", as.integer(substr(time, 1, 2)), as.integer(substr(time, 3, 4)))
    )
  
  return(df)
}

all_data <- file_paths %>%
  map_df(read_prof_file) %>%
  group_by(date, time) %>%
  mutate(rep = cur_group_id()) %>%
  ungroup() %>%
  mutate_at(vars(1:28), as.numeric)

cleaned_data <- all_data %>%
  filter(!rep %in% c(65|66|67)) #empty datasets 

cleaned_data %>%
  group_by(date) %>%
  summarise(rep_count = n_distinct(rep)) 

times <- cleaned_data %>%
  group_by(date, time) %>%
  summarise(rep_count = n_distinct(rep)) 

plot_wind_speed <- cleaned_data %>%
  ggplot() +
  geom_point(aes(Aux2, X, color=rep)) +
  xlim(0, 4.5)

###############################
#Selected day
###############################
# 7 may = IOP 1, LAI 4.5, height 0.45
# 15 june = IOP 2, LAI 5.5, height 0.80
LAI4.5 <- 4.5
LAI5.5 <- 5.5
H0.45 <- 0.45
H0.80 <- 0.80

#Averages for an afternoon for one day
df_sum_iop1 <- cleaned_data %>%
  filter(date %in%as.Date(c("2018-07-02")),
         X <= H0.80) %>% 
  mutate(time = lubridate::hms(paste0(time, ":00"))) %>%  
  filter(time >= lubridate::hms("12:00:00") & time <= lubridate::hms("16:00:00")) %>%  
  group_by(X) %>%
  summarize(
    across(c(Aux2, CO2B, H2OB, FW_DIFF2), list(
      mean = ~mean(.x, na.rm = TRUE),
      n = ~sum(!is.na(.x)),
      se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
    ), .names = "{.col}_{.fn}")
  )

###############################
#PAR profiles
###############################
# Light extinction in wheat canopy (Burgess, 2015)
# https://doi.org/10.1104/pp.15.00722
# Note: LAI  is > 5 but height < 0.50
H47 <- 0.47 # height of the canopy in the paper
calculate_PAR_canopy <- function(PAR_values, z, h) {
  d <- h - z # depth in canopy (m)
  k <- 0.49 # For line 1 in Burgess
  
  clai_max_original <- 11.824 * H47 - 0.0105  # Original maximum cumulative lai
  clai <- clai_max_original * (d / h)  # Proportionally distribute LAI with depth
  
  return(PAR_values * exp(-k * clai)) 
}

df_sum_iop1 <- df_sum_iop1 %>%
  mutate(PAR = calculate_PAR_canopy(1600, X, H0.80))


###############################
# Gb profiles
###############################
Pa <- 101325 # pressure
R <- 8.31446261815324 # gas constant
Dw20 <- 24.2 # diffusion coefficient water, Jones, H. G. (1992). Plants and microclimate: a quantitative approach to environmental plant physiology.
Dh20 <- 21.5 # diffusion coefficient heat
v20 <- 15.1 #k inematic viscosity 
d_wheat <- 0.015# leaf width (Moeller, Evers, Rebetzke, 2014 (https://doi.org/10.3389/fpls.2014.00617) and Evers et al, 2004 (https://doi.org/10.1016/j.ecolmodel.2006.07.042)

temperature_correction <- function (T_C, D){ #Jones 2012
  return (D * ((273.12+T_C)/293.2)^1.75)
}

calculate_gb <- function(d, T_C, u){ # Brenner and Jarvis, 1995 (https://doi.org/10.1016/0168-1923(94)02160-L)
  Dw <- temperature_correction(T_C, Dw20)
  v <- temperature_correction(T_C, v20)
  T_K <- T_C+273.15
  gb_mm_w <- 1*(0.664*Dw^(2/3)*u^0.5)/(d^0.5*v^0.17) # one sided
  gb_w <- Pa/(R*T_K)*gb_mm_w/1000
  return(gb_w)
}

df_sum_iop1 <- df_sum_iop1 %>% 
  mutate(gb = calculate_gb(d_wheat, FW_DIFF2_mean, Aux2_mean)) 

###############################
# Photosynthesis, stomata and energy balance models
###############################
LatentHeat <- function(Tk) { #J kg -1
  return(1.91846e6 * (Tk / (Tk - 33.91))^2) # latent heat of vaporization (Henderson-Sellers, 1984)
}
SatVap <- function(Tc) {
  return(613.65 * exp(17.502 * Tc / (240.97 + Tc))) # Buck's equation
}
Ta <- function(T, c, dHa) { # Harley et al. 1992 (https://doi.org/10.1111/j.1365-3040.1992.tb00974.x)
  return(exp(c - dHa / (8.314e-3 * (273.15 + T)))) # Dependence of reaction rate on temperature
}
Tad <- function(T, c, dHa, dHd, dS) {
  return(exp(c - dHa / (8.314e-3 * (273.15 + T))) / (1.0 + exp((dS * (T + 273.15) - dHd) / (8.314e-3 * (T + 273.15))))) # Note: there is a mistake in the formula in sharkeys paper (https://doi.org/10.1111/j.1365-3040.2007.01710.x), which is corrected here
} 
NRH <- function(PAR, Abs, Jmax, Curv) {
Qabs <- Abs*PAR
return((Qabs+Jmax-sqrt((Qabs+Jmax)^2-4*Curv*Qabs*Jmax))/(2*Curv))
}

#FvCB
FvCB <- function(PAR, Ci, T, Pa, Vcmax25, Jmax25, gm25, Abs, Curv, Rd25) {
  Pi <- Ci * Pa * 1e-6 # Atmospheric pressure
  
  # Sharkey 2007
  Kc <- Ta(T, 35.9774, 93.72) # Michaelis constant rubisco, Silva-Perez 2017 (https://doi.org/10.1111/pce.12953)
  Ko <- Ta(T, 12.3772, 33.6) # Inhibition constant, Silva-Perez 2017 (https://doi.org/10.1111/pce.12953)
  O <- 21.0 # Oxygen
  Km <- Kc * (1.0 + O / Ko) # Diffusion resistance
  GS <- Ta(T, 11.187, 24.46) # GS
  
  Vcmax <- Vcmax25 * Ta(T, 26.355, 65.33)  # Bernacchi 2001
  Jmax <- Jmax25 * Ta(T, 17.71, 43.9) # Bernacchi 2003
  gm <- gm25 * Tad(T, 20.01, 49.6, 437.4, 1.4) # Bernacchi 2002
  Rd <- Rd25 * Ta(T, 18.715, 46.39) # Bernacchi 2001
  
  J <- NRH(PAR, 0.5 * Abs, Jmax, Curv)
  
  # Ethier 2004
  a <- -1.0 / gm
  b <- (Vcmax - Rd) / gm + Pi + Km
  c <- Rd * (Pi + Km) - Vcmax * (Pi - GS)
  d <- b^2 - 4.0 * a * c
  Ac <- (-b + sqrt(max(0.0, d))) / (2.0 * a)
  
  b <- (J / 4.0 - Rd) / gm + Pi + 2.0 * GS
  c <- Rd * (Pi + 2.0 * GS) - J / 4.0 * (Pi - GS)
  d <- b^2 - 4.0 * a * c
  Aj <- (-b + sqrt(max(0.0, d))) / (2.0 * a)
  
  # CO2 compensation point
  G <- (Rd * Km + Vcmax * GS) / (Vcmax / Rd)
  
  return(list(min(Ac, Aj), G, Rd))
}

Leuning <- function(A, Ds, g0, g1, Cs, G, Ds0) {
  return(g0 + g1 * A / ((Cs - G) * (1.0 + Ds / Ds0)))
}

energy_balance <- function(PAR, Tair, E, gb_leaf, RH) { #taken from LI6800 (https://www.licor.com/env/support/LI-6800/topics/equation-summary.html)
  ea <- SatVap(Tair) * RH
  SH <- 0.622 * ea / (Pa - ea) # Specific humidity
  Cp <- 1005.0 + 1820.0 * SH  # Specific heat capacity air J kg-1 K-1
  Tw <- Tair 
  
  Rabs <- PAR * Abssw * kSun 
  gbm <- gb_leaf * 8.314 * (Tair + 273.15) / Pa # mol m-2 s-1 * J/(mol K) * K / Pa = m s-1
  
  net_rad <- Rabs + 2 * 0.95 * 5.67e-8 * ((Tw + 273.15)^4 - (Tair + 273.15)^4)  # Net radiation balance
  tran <- LatentHeat(Tair + 273.15) * E * 0.01801528 # Transpiration in W m-2

  dT <- (net_rad - tran) / (1.84 * Cp * gbm + 8 * 0.95 * 5.67e-8 * (Tair + 273.15)^3)
  Tl <- Tair + dT
  
  return(Tl) 
}

leaf_transpiration <- function(gs, gb, RH, Tair, Tl){
  ea <- SatVap(Tair) * RH
  es <- SatVap(Tl) #Pa
  gt <- 1.0 / (1.0 / gb  + 1.0 / gs) 
  E <- gt * (es - ea) / 101325  #mol m-2 s-1 
  return(E)
}

ci_intercellular <- function(Ca, A, gs, gb) {
  gtc <- 1.0 / (1.6 / gs + 1.37 / gb)
  Ci <- Ca - A / gtc
}

ds_surface <- function(Ca, Tl, Tair, A, gb, gs, RH) {
  # Aphalo (1993), Eq 8 (https://doi.org/10.1111/j.1365-3040.1993.tb00499.x)
  
  ea <- SatVap(Tair) * RH
  vpd <- SatVap(Tl) - ea
  
  Csurface <- Ca - 1.37 * A / gb
  gt <- 1.0 / (1.0 / gb  + 1.0 / gs)
  Ds <- vpd * (1.0 - gt / gb) * 1e-3
  return(list(Ds, Csurface))
}

###############################
#Parameters
###############################
Pa <- 101325 #pressure
Vcmax25 <- 140 # Mcausland 2020
Jmax25 <- 220 # Mcausland 2020
Abs <- 0.86 # Mcausland 2020
Abssw <- 0.5 
Rd25 <- 0.9 # https://doi.org/10.1111/nph.17538
gm25 <- 5 # Barbour & Kaiser 2016; Tazoe et al 2009
g0 <- 0.02 # corresponds to lower values in Leuning et al for a variety of species (https://doi.org/10.1111/j.1365-3040.1995.tb00370.x)
g1 <- 7.38 # yu et al, 2004
Curv <- 0.7
kSun <- 1 / (500 / 120 * 0.47)
Ds0 <- 2 

###############################
#Optimization
##############################
photosynthesis_iteration <- function(Ca, Tair, PPFD, gb, RH){
  #Initial parameters for estimation
  Tl <- Tair
  Ci <- 0.7 * Ca
  gs <- 0.5
  
  threshold_Ci <- 0.001
  threshold_Tl <- 0.001
  threshold_gs <- 0.001
  iteration_count <- 0
  max_iterations <- 1000
  
  repeat {
    iteration_count <- iteration_count + 1  
    print(paste("Iteration:", iteration_count))
    
    results <- FvCB(PPFD, Ci, Tl, Pa, Vcmax25, Jmax25, gm25, Abs, Curv, Rd25)
    A <- results[[1]]
    G <- results[[2]]
    Rd <- results[[3]]
      
    
    results_surface <- ds_surface(Ca, Tl, Tair, A, gb, gs, RH)
    Ds_new <- results_surface[[1]]
    Cs_new <- results_surface[[2]]
    gs_new <- Leuning(A, Ds_new, g0, g1, Cs_new, G, Ds0)
    
    Ci_new <- ci_intercellular(Ca, A, gs_new, gb)
    print(paste("Current Ci:", Ci, "New Ci:", Ci_new))
    
    # Damping to prevent oscillation due to big changes
    Ci_new <- (Ci_new + Ci) / 2
    
  
    #If Ci converges, proceed to check Tl and gs 
    if (abs(Ci - Ci_new) < threshold_Ci) {
      Tl_new <- energy_balance(PPFD, Tair, leaf_transpiration(gs_new, gb, RH, Tair, Tl), gb, RH)
      print(paste("Current Tl:", Tl, "New Tl:", Tl_new))
      print(paste("Current gs:", gs, "New gs:", gs_new))
      
      #Check if both Tl and gs have converged
      if (abs(Tl - Tl_new) < threshold_Tl && abs(gs - gs_new) < threshold_gs) {
        break
        }
    }
    else {
      Tl_new <- Tl #Ensure Tl is always defined even if Ci does not converge
    }
    
    if (iteration_count >= max_iterations) {
      warning("No convergence.")
      break
      
      }
    #Update values for next iteration
    Ci <- Ci_new 
    Tl <- Tl_new
    gs <- gs_new
  }
  
  return(tibble(Ci, Tl, gs))
}

df_final <- df_sum_iop1 %>%
  mutate(RH = (H2OB_mean*Pa/1000)/SatVap(FW_DIFF2_mean)) %>%
  rowwise() %>%
  mutate(results = list(photosynthesis_iteration(CO2B_mean, FW_DIFF2_mean, PAR, gb, RH))) %>% #
  unnest(cols = c(results)) %>%
  rowwise() %>%
  mutate(photosynthesis = FvCB(PAR, Ci, Tl, Pa, Vcmax25, Jmax25, gm25, Abs, Curv, Rd25)[[1]]) %>%
  ungroup()

write.csv(df_final, "0702_afternoon_PAR1600_d015_gbvar.csv")


###############################
#No gb limitation (i.e. complete coupling to air)
###############################
df_final_infinite_gb <- df_sum_iop1 %>%
  mutate(RH = (H2OB_mean*Pa/1000)/SatVap(FW_DIFF2_mean)) %>%
  rowwise() %>%
  mutate(results = list(photosynthesis_iteration(CO2B_mean, FW_DIFF2_mean, PAR, Inf, RH))) %>% 
  unnest(cols = c(results)) %>%
  rowwise() %>%
  mutate(photosynthesis = FvCB(PAR, Ci, Tl, Pa, Vcmax25, Jmax25, gm25, Abs, Curv, Rd25)[[1]]) %>%
  ungroup()      

write.csv(df_final_infinite_gb, "0702_afternoon_PAR1600_d015_gbinf.csv")

ggplot() +
  geom_point(data = df_final, aes(x = X, y = photosynthesis)) +  
  geom_line(data = df_final_infinite_gb, aes(x = X, y = photosynthesis), color="red") +  
  theme_classic() 
