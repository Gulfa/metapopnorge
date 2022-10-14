
library(dplyr)
library(metapop)
library(metapopnorge)
library(ggplot2)


## Basic usage of metapop
L <- 200
n_vac <- 1
n_strain <- 1
N <- 5

basic_params <- list(
  N_steps=L,
  n_vac=n_vac,
  n_strain=n_strain,
  dt=1,## 
  T_waning=array(1e10, dim=c(N, n_vac)),
  vaccinations=array(0,dim=c(L, N, n_vac)),
  beta_day=matrix(0.05, ncol=N, nrow=L),
  mixing_matrix=matrix(1.0, nrow=N, ncol=N),
  tot_vac_ini=array(0, dim=c(N,n_vac)),
  migration_matrix=matrix(0.0, nrow=N, ncol=N),
  latent_period=3.0,
  beta_strain=rep(1, n_strain),
  cross_protection=matrix(1, ncol=n_strain, nrow=n_strain),
  infectious_period=4.0,
  import_vec=array(0,dim=c(L,N,n_vac,n_strain)),
  length_hosp=array(4,dim=c(N,n_vac,n_strain)),
  length_icu=array(15,dim=c(N,n_vac,n_strain)),
  hosp_prob=array(0.04, dim=c(N,n_vac,n_strain)),
  icu_prob=array(0.5, dim=c(N,n_vac,n_strain)),
  pre_icu=array(2,dim=c(N,n_vac,n_strain)),
  post_icu=array(2,dim=c(N,n_vac,n_strain)),
  time_before_death=4,
  time_before_death_hosp=3,
  time_before_death_icu=2,
  susceptibility=array(1,dim=c(N,n_vac,n_strain)),
  transmisibility=array(1,dim=c(N,n_vac,n_strain)),
  susceptibility_asymp=array(1,dim=c(N,n_vac,n_strain)),
  susceptibility_symp=array(1,dim=c(N,n_vac,n_strain)),
  symp_trans=array(1.0,dim=c(N,n_vac,n_strain)),
  prob_death_non_hosp=array(0.03, dim=c(N,n_vac,n_strain)),
  prob_death_hosp=array(0.1, dim=c(N,n_vac,n_strain)),
  prob_death_icu=array(0.5,dim=c(N,n_vac,n_strain)),
  waning_immunity_vax = array(10000, dim=c(N,n_vac,n_strain)),
  sympt_frac=array(0.7,dim=c(N,n_vac,n_strain)),
  asympt_frac=array(0.3,dim=c(N,n_vac,n_strain)),
  pre_sympt_infect=1.2,
  asympt_infect=0.2,
  pre_sympt_period=2,
  n=N,
  S_ini=matrix(1e5, nrow=N, ncol=n_vac),
  I_ini=array(20, dim=c(N,n_vac,n_strain)),
  I_imp_ini=array(0, dim=c(N,n_vac,n_strain)),
  Ea_ini=array(0, dim=c(N,n_vac,n_strain)),
  Es_ini=array(0, dim=c(N,n_vac,n_strain)),
  A_ini=array(0, dim=c(N,n_vac,n_strain)),
  P_ini=array(0, dim=c(N,n_vac,n_strain)),
  H_ini=array(0, dim=c(N,n_vac,n_strain)),
  ICU_H_ini=array(0, dim=c(N,n_vac,n_strain)),
  
  ICU_R_ini=array(0, dim=c(N,n_vac,n_strain)),
  ICU_P_ini=array(0, dim=c(N,n_vac,n_strain)),
  B_D_ini=array(0, dim=c(N,n_vac,n_strain)),
  B_D_H_ini=array(0, dim=c(N,n_vac,n_strain)),
  B_D_ICU_ini=array(0, dim=c(N,n_vac,n_strain)),
  R_ini=array(0, dim=c(N,n_vac,n_strain)),
  D_ini=array(0, dim=c(N,n_vac,n_strain)),
  tot_infected_ini=array(0, dim=c(N,n_vac,n_strain)),
  tot_hosp_ini=array(0, dim=c(N,n_vac,n_strain)),
  tot_resp_ini=array(0, dim=c(N,n_vac,n_strain)),
  beta_norm=rep(1e5,N),
  reg_pop_long=rep(1e5,N),
  waning_inf=5e10,
  N_regions=1,
  rand_beta_sd=0.1,
  rand_beta_factors=rep(0.05,N),
  beta_cut_peak_param=c(0,0,0,0),
  age_groups=N,
  change_factor=c(0,0,0,0),
  threshold=c(100,0),
  beta_mode=1,
  spont_behav_change_params=c(0,0,0,0),
  expected_health_loss=array(1, dim=c(N,n_vac))
)

results <- run_params(basic_params, L=200, N_particles = 30, N_threads=3)

names(results)

ggplot(results) + geom_line(aes(x=t, y=incidence, group=sim))




## Using convenience functions in metapopnorge


param_file <- "parameters_example.xlsx"
#Vaccination RRs for unvaccinated and vaccinated
vac_pars <- list(rr_inf = c(1,0.7),
              rr_hosp = c(1, 0.6),
              rr_death = c(1, 0.6),
              rr_icu = c(1, 0.6),
              rr_los_hosp = c(1, 0.85),
              rr_inf_asymp = c(1,0.8),
              rr_trans = c(1,1))


new_params <- get_variant_params(param_file, vac_pars=vac_pars)

R <- 1.3
dt <- 0.5

# beta_1 normalises beta such that we get R=1
beta_1 <- fix_beta_large(new_params, new_params$S_ini, new_params$I_ini, 1, beta=new_params$beta_day[1,], use_eig=TRUE)
beta_0 <- 0.2*beta_1*R
new_params$beta_day <- matrix(beta_0, ncol=9, nrow=500)
new_params <- change_dt(new_params, dt)
#new_params <- update_severity(new_params, severity, change_icu_prob=change_icu_prob, change_los_hosp=change_los_hosp)
res <- run_params(new_params, L=168, N_particles = 30, N_threads=10)

colnames(res)
ggplot(res) + geom_line(aes(x=t, y=incidence, group=sim))

ggplot(res %>% select(t, sim, paste0("tot_infected_age_", 1:9)) %>% tidyr::pivot_longer(starts_with("tot_"))) +
  geom_line(aes(x=t, y=value, color=name, group=paste(sim, name)))




params <- get_params(param_file = "parameter_files/parameters_example_large.xlsx",
                     N_age=18,
                     L=212,
                     n_vac=3,
                     n_strain=1,
                     dose_file = "",
                     vaccine_file = "parameter_files/vaccine_profile_example.csv",
                     vac_priority=NULL,
                     regional_file = glue::glue("parameter_files/geo_0.csv"),
                     initial_conditions_file = "parameter_files/epi_scenario_municip_ibm.txt",
                     import_file = "parameter_files/import_M.csv",
                     import_age_dist_file = "parameter_files/import_age_dist.txt",
                     mobility_file = "parameter_files/municip_matrix_fake_data.RDS",
                     regions = TRUE,
                     add_plus_minus_neutral=TRUE)
dt <- 1

params <- change_dt(params, dt)



beta_1 <- fix_beta_large(params, params$S_ini, params$I_ini, 1, beta=params$beta_day[1,], use_eig=TRUE)
params$beta_day <- params$beta_day*beta_1*R
res <- run_params(params, L=168, N_particles = 30, N_threads=10)

colnames(res)
ggplot(res) + geom_line(aes(x=t, y=tot_infected, group=sim))

ggplot(res %>% select(t, sim, paste0("tot_infected_", c("plus", "minus", "neutral"))) %>% tidyr::pivot_longer(starts_with("tot_"))) +
  geom_line(aes(x=t, y=value, color=name, group=paste(sim, name)))


