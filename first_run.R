library(devtools)
#install_github("mrc-ide/malariasimulation")
install_github("mrc-ide/malariasimulation", ref="vivax_vaccine", force = TRUE)
#install_github("mrc-ide/malariasimulation", ref="vivax_main")

# main is installed

# TO DO: In human_infection, need a separate source_humans for vivax!!
# Might be different for infections and patent infections too


library(malariasimulation)
library(tidyverse)
library(ggplot2)

year <- 365
month <- year/12
warmup <- 25*365 
sim_length <- 10*365 
starting_EIR <- 1     # EIR = 1 ~ 12.7% PCR prevalence, EIR = 0.05 ~ 1.9% PCR prevalence
population <- 10000

# Run vivax

# get starting parameters ----------
params <- get_parameters(list(
  human_population = population,
  model_seasonality = FALSE,
  individual_mosquitoes = FALSE), parasite = "vivax")

params$clinical_incidence_rendering_min_ages = 0
params$clinical_incidence_rendering_max_ages = 100*year
params$severe_incidence_rendering_min_ages = 0
params$severe_incidence_rendering_max_ages = 100*year  

#params$clinical_incidence_rendering_min_ages = c(0, 0.25, seq(5,100,5))*year
#params$clinical_incidence_rendering_max_ages = c(0.25, seq(5,100,5),200)*year
#params$severe_incidence_rendering_min_ages = c(0, 0.25, seq(5,100,5))*year
#params$severe_incidence_rendering_max_ages = c(0.25, seq(5,100,5),200)*year
# params$clinical_incidence_rendering_min_ages = c(0, 0.25, 1,2,3,4, seq(5,100,5))*year
# params$clinical_incidence_rendering_max_ages = c(0.25, 1,2,3,4, seq(5,100,5),200)*year
# params$severe_incidence_rendering_min_ages = c(0, 0.25, 1,2,3,4, seq(5,100,5))*year
# params$severe_incidence_rendering_max_ages = c(0.25, 1,2,3,4, seq(5,100,5),200)*year
params$prevalence_rendering_min_ages = c(0, 2 * year)
params$prevalence_rendering_max_ages = c(100*year, 10 * year)

# Treatment

# params <- set_clinical_treatment(
#   parameters = params,
#   drug = 1,
#   timesteps = c(1),
#   coverages = c(0.45)
# ) 

# params <- set_pev_epi(
#   params,
#   profile = rtss_profile, # We will model implementation of the RTSS vaccine.
#   timesteps = 30 * year, # Vaccination will begin at 8 years into the simulation.
#   coverages = 0.8, 
#   min_wait = 0, # There is no minimum wait since the last vaccination.
#   age = 5 * month, # Individuals will be vaccinated once they reach 5 months of age.
#   booster_timestep = 12 * month, # The booster is administered 12 months following the third dose. 
#   booster_coverage = 0.8, 
#   booster_profile = list(rtss_booster_profile) # We will model implementation of the RTSS booster.
# )

params <- set_mass_pev(
  params,
  profile = vivax_pev_profile, #rtss_profile,
  timesteps = 30 * year, # implemented in year 30
  coverages = 0.8,
  min_ages = 5*month,
  max_ages = 60*year,
  min_wait = 0,
  booster_timestep = 12*month,
  booster_coverage = 0,
  booster_profile = list(vivax_pev_booster_profile)
)

# EIR equilibrium ----------
set.seed(123)
params <- set_equilibrium(params, as.numeric(starting_EIR))

# run simulation ----------
set.seed(123)

output <- run_simulation(
  timesteps = warmup + sim_length,
  parameters = params)

output2 <- output %>% 
  # add vars to output
  mutate(EIR = starting_EIR,
         warmup = warmup,
         sim_length = sim_length,
         population = population,
         timestep = timestep - warmup) %>%
  ungroup() %>%
  filter(timestep > 0) %>% # remove warmup period
  # statistics by month
  mutate(year = ceiling(timestep/year),
          month = ceiling(timestep/month)) 

# Vivax processing
output_pv <- output2 %>%
   
  # only necessary variables
  dplyr::select(EIR, warmup, sim_length, population, month, year,
                starts_with("n_inc_severe"), starts_with("p_inc_severe"),
                starts_with("n_inc"), starts_with("p_inc"),
                starts_with("n_detect"), starts_with("p_detect"),
                starts_with("n_"), -n_bitten,  n_new_bite_infections, n_new_relapse_infections) %>%
  
  # take means of populations and sums of cases by month
  group_by(EIR, warmup, sim_length, population, year) %>%

  # why is there pcr and lm with falciparum?? Also no severe...
  mutate_at(vars(n_0_36500, n_730_3650,
                 n_detect_pcr_730_3650,  n_detect_lm_730_3650), mean, na.rm = TRUE) %>%
  mutate_at(vars(n_inc_clinical_0_36500:p_inc_clinical_0_36500,
                 n_new_bite_infections, n_new_relapse_infections), sum, na.rm = TRUE) %>%
   
   dplyr::select(EIR:year, 
                 n_0_36500, n_730_3650,
                 n_detect_pcr_730_3650,  n_detect_lm_730_3650,
                 n_inc_clinical_0_36500:p_inc_clinical_0_36500,
                 n_new_bite_infections, n_new_relapse_infections) %>%
   
   distinct()


# Run falciparum

# get starting parameters ----------
params <- get_parameters(list(
  human_population = population,
  model_seasonality = FALSE,
  individual_mosquitoes = FALSE), parasite = "falciparum")

params$clinical_incidence_rendering_min_ages = 0
params$clinical_incidence_rendering_max_ages = 100*year
params$severe_incidence_rendering_min_ages = 0
params$severe_incidence_rendering_max_ages = 100*year  

#params$clinical_incidence_rendering_min_ages = c(0, 0.25, seq(5,100,5))*year
#params$clinical_incidence_rendering_max_ages = c(0.25, seq(5,100,5),200)*year
#params$severe_incidence_rendering_min_ages = c(0, 0.25, seq(5,100,5))*year
#params$severe_incidence_rendering_max_ages = c(0.25, seq(5,100,5),200)*year
# params$clinical_incidence_rendering_min_ages = c(0, 0.25, 1,2,3,4, seq(5,100,5))*year
# params$clinical_incidence_rendering_max_ages = c(0.25, 1,2,3,4, seq(5,100,5),200)*year
# params$severe_incidence_rendering_min_ages = c(0, 0.25, 1,2,3,4, seq(5,100,5))*year
# params$severe_incidence_rendering_max_ages = c(0.25, 1,2,3,4, seq(5,100,5),200)*year
params$prevalence_rendering_min_ages = c(0, 2 * year)
params$prevalence_rendering_max_ages = c(100*year, 10 * year)

# Treatment

# params <- set_clinical_treatment(
#   parameters = params,
#   drug = 1,
#   timesteps = c(1),
#   coverages = c(0.45)
# ) 
 
params <- set_pev_epi(
  params,
  profile = rtss_profile, # We will model implementation of the RTSS vaccine.
  timesteps = 9 * year, # Vaccination will begin at 9 years into the simulation.
  coverages = 0.8, 
  min_wait = 0, # There is no minimum wait since the last vaccination.
  age = 5 * month, # Individuals will be vaccinated once they reach 5 months of age.
  booster_timestep = 12 * month, # The booster is administered 12 months following the third dose. 
  booster_coverage = 0.8, 
  booster_profile = list(rtss_booster_profile) # We will model implementation of the RTSS booster.
)

# EIR equilibrium ----------
params <- set_equilibrium(params, as.numeric(starting_EIR))

# run simulation ----------
set.seed(123)

sim_pf <- run_simulation(
  timesteps = warmup + sim_length,
  parameters = params)

output2 <- sim_pf %>% 
  # add vars to output
  mutate(EIR = starting_EIR,
         warmup = warmup,
         sim_length = sim_length,
         population = population,
         timestep = timestep - warmup) %>%
  ungroup() %>%
  filter(timestep > 0) %>% # remove warmup period
  # statistics by month
  mutate(year = ceiling(timestep/year),
         month = ceiling(timestep/month)) 

output_pf <- output2 %>%
  
  # only necessary variables
  dplyr::select(EIR, warmup, sim_length, population, month, year,
                starts_with("n_inc"), starts_with("p_inc"),
                starts_with("n_detect"), starts_with("p_detect"),
                starts_with("n_"), -n_bitten,  n_infections) %>%
  
  # take means of populations and sums of cases by month
  group_by(EIR, warmup, sim_length, population,month, year) %>%
  
  # why is there pcr and lm with falciparum?? Also no severe...
  mutate_at(vars(n_0_36500, n_730_3650,
                 n_detect_pcr_730_3650,  n_detect_lm_730_3650, p_detect_lm_730_3650), 
            mean, na.rm = TRUE) %>%
  mutate_at(vars(n_inc_clinical_0_36500:p_inc_clinical_0_36500,
                 n_infections), sum, na.rm = TRUE) %>%
 
 dplyr::select(EIR:year, 
                n_0_36500, n_730_3650,
                n_detect_pcr_730_3650,  n_detect_lm_730_3650,
                n_inc_clinical_0_36500:p_inc_clinical_0_36500, n_infections) %>%

 distinct()
 

# Plot output

#output_pv_vacc <- output_pv

ggplot() +
  #geom_point(data= output_pv_vacc, aes(x=year, y = n_inc_clinical_0_36500*1000/n_0_36500, col = "Pre-erythrocytic\nvaccine")) +
  geom_point(data= output_pv, aes(x=year, y = n_inc_clinical_0_36500*1000/n_0_36500, col = "Baseline")) +
  labs(title = "P. vivax, EIR = 1", x = "Years", y = "Clinical cases per 1000 person-years (all ages)",
       col = "Scenario") +
  ylim(0,400) +
  theme_bw()
#ggsave("preerythrocytic_eir1.png")

ggplot() +
  geom_point(data= output_pv_vacc, aes(x=year, y =  n_inc_clinical_0_36500*1000/n_0_36500, col = "Pre-erythrocytic\nvaccine")) +
  geom_point(data= output_pv, aes(x=year, y =  n_inc_clinical_0_36500*1000/n_0_36500, col = "Baseline")) +
  labs(title = "P. vivax, EIR = 0.05", x = "Years", y = "Clinical cases per 1000 person-years (all ages)",
       col = "Scenario") +
  ylim(0,150) +
  theme_bw()
#ggsave("preerythrocytic_eir0.05.png")

###

ggplot(output_pv) +
  geom_point(aes(x=month/12, y = n_inc_clinical_0_36500)) +
  labs(title = "P. vivax", x = "Years", y = "Number of cases (all ages)") +
  ylim(0,1200)

ggplot(output_pf) +
  geom_point(aes(x=month/12, y = n_inc_clinical_0_36500)) +
  labs(title = "P. falciparum", x = "Years", y = "Number of cases (all ages)") +
  ylim(0,1200)

ggplot(output_pv) +
  geom_point(aes(x=month/12, y = n_detect_lm_730_3650/n_730_3650)) +
  labs(title = "P. vivax", x = "Years", y = "Parasite prevalence") +
  ylim(0,0.2)

ggplot(output_pf) +
  geom_point(aes(x=month/12, y = n_detect_lm_730_3650/n_730_3650)) +
  labs(title = "P. falciparum", x = "Years", y = "Parasite prevalence") +
  ylim(0,1)





ggplot(output_pv_vacc) +
  geom_point(aes(x=month/12, y = n_inc_clinical_0_36500)) +
  labs(title = "P. vivax", x = "Years", y = "Number of cases (all ages)") +
  ylim(0,1200)

ggplot(output_pv_vacc) +
  geom_point(aes(x=month/12, y = n_detect_lm_730_3650/n_730_3650)) +
  labs(title = "P. vivax", x = "Years", y = "Parasite prevalence") +
  ylim(0,1)


ggplot(output_pf_vacc) +
  geom_point(aes(x=month/12, y = n_inc_clinical_0_36500)) +
  labs(title = "P. falciparum", x = "Years", y = "Number of cases (all ages)") +
  ylim(0,1200)

ggplot(output_pf_vacc) +
  geom_point(aes(x=month/12, y = n_detect_lm_730_3650/n_730_3650)) +
  labs(title = "P. falciparum", x = "Years", y = "Parasite prevalence") +
  ylim(0,1)


t <- seq(0,365*3,1)
plot(t/365, 0.6*exp(-log(2)/365*t))



