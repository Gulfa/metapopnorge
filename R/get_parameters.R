#' @import abind
#' @import data.table

NULL

#' @export
get_variant_params <- function(param_file, daily_import=0,
                               vac_pars=list(rr_inf = c(1,0.7),
                                             rr_hosp = c(1, 0.6),
                                             rr_death = c(1, 0.6),
                                             rr_icu = c(1, 0.6),
                                             rr_los_hosp = c(1, 0.85),
                                             rr_inf_asymp = c(1,0.8),
                                             rr_trans = c(1,1)), L=500) {

  params <- read_param_file_simplified(param_file)
  n_vac <- 2
  n_strain <- 1
  N <- 9
  params <- c(params,
              list(
                N_steps=L,
                n_vac=n_vac,
                n_strain=n_strain,
                dt=0.5,
                T_waning=array(1e10, dim=c(N, n_vac)),
                vaccinations=array(0,dim=c(L, N, n_vac)),
                                        #    import_vec=rep(0, L),
                                        #    import_age_prio=rep(1,dim=(N,n_vac,1)),
                beta_day=matrix(0.2, ncol=N, nrow=L),
                                        #    vac_time_full_effect=array(14.0, N),
                beta_strain=rep(1, n_strain),
                cross_protection=matrix(1, ncol=n_strain, nrow=n_strain),
                n=9,
                S_ini=matrix(c(get_age_groups(), get_age_groups()*0), nrow=N, ncol=n_vac),
                import_vec=array(0, dim=c(L,N, n_vac, n_strain)),
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
                tot_vac_ini=array(0, dim=c(N, n_vac)),
                tot_vac_adm_ini=array(0, dim=c(N, n_vac)),
                beta_norm=get_age_groups(),
                reg_pop_long=get_age_groups(),
                N_regions=1,
                waning_immunity_vax = array(100000, dim=c(N,n_vac,n_strain)),
                waning_inf = 100000,
                age_groups=N,
                beta_mode=1
              )
              )
  params$I_ini[,2,] <- 0
  
  params$import_vec[1:L, 1:7, 1, 1] <- round(daily_import/7)
  basic_params <- fix_params(params, N, n_vac, n_strain, vac_pars)
}

#' @export
get_age_groups <- function(){
  ## spldata::norway_population_by_age_cats(cats = list(c(1:10), c(11:20), c(21:30), c(31:40),
  ##                                                    c(41:50), c(51:60), c(61:70), c(71:80), c(80:120))) %>% dplyr::filter(year==2022 & location_code=="norge") %>% pull(pop)
  #spldata::nor_population_by_age_cats(cats=list("1"=-1:9, "2"=10:19, "3"=20:29, "4"=30:39, "5"=40:49, "6"=50:59, "7"=60:69, "8"=70:79, "9"=80:120),
                                        #                                    border=2020, include_total=FALSE) %>% filter(calyear==2022 & location_code=="norge")%>%pull(pop_jan1_n)
  return(c(587714,647020,701583,749122,714233,724859,596118,464328,240293))

}

get_vac_hist_sysvak <- function(filename, reg_data, start_date, L, vac_adherence=NULL){
  sysvak_data <- fread(filename)
  sysvak_data$date_vax <- as.Date(sysvak_data$date_vax)
  dat <- sysvak_data %>% left_join(reg_data$reg_spec_municip, by=c("municip_code"="fhidata.municip_code")) %>%
    mutate(risikogruppe=replace(risikogruppe, risikogruppe > 1, 1)) %>%
    group_by(date_vax, aldersgruppe, risikogruppe, name) %>% summarize(n=sum(n)) %>% filter(date_vax >= start_date & date_vax < start_date + L & !is.na(name))

  skeleton <- expand.grid(date_vax=seq(min(dat$date_vax), max(dat$date_vax), by=1), name=unique(dat$name), aldersgruppe=unique(dat$aldersgruppe), risikogruppe=unique(dat$risikogruppe))

  final_df <- skeleton %>% left_join(dat, on=c("date_vax"="date_vax", "name"="name", "aldersgruppe"="aldersgruppe", "risikogruppe"="risikogrupp")) %>% mutate(n=tidyr::replace_na(n, 0)) %>% arrange(date_vax, name, risikogruppe,  aldersgruppe)

  m <- matrix(final_df %>% pull(n), nrow=L, byrow=TRUE)
  
  if(!is.null(vac_adherence)){
    max_doses_to_give <- rep(vac_adherence$adherence, reg_data$N_regions) * reg_data$pop
    max_doses_matrix <- matrix(max_doses_to_give, nrow=L, ncol=length(max_doses_to_give), byrow=T)
    doses_given <- matrixStats::colCumsums(m)
    mask <- doses_given > max_doses_matrix
    m[mask] <- 0
  }
  
  return(m)
}


get_S_ini <- function(filename, start_date, n_strian){
  sysvak_data <- fread(filename)
  
  sysvak_data <- tidyr::complete(sysvak_data, age_group, dosenummer, date_vax, fill=list(n=0))
  
  tmp <- sysvak_data %>% filter(date_vax<=start_date-7) %>%group_by(age_group,dosenummer) %>% summarize(n= sum(n))# %>%#, by=.(age_group, dosenummer)][order(age_group)]
  
  third <- as.data.table(tmp %>% dplyr::filter(dosenummer=="tredje")) %>% pull(n) #%>% arrange(c(1,9,2,3,4,5,6,7,8))
  second <- as.data.table(tmp %>% dplyr::filter(dosenummer=="andre")) %>% pull(n) - third#%>% arrange(c(1,9,2,3,4,5,6,7,8))
  first <- as.data.table(tmp %>% dplyr::filter(dosenummer=="forste")) %>% pull(n) - third - second#%>% arrange(c(1,9,2,3,4,5,6,7,8))
  unvax <- get_age_groups() - third - second - first
  unvax[unvax <0] <- 0

  return(array(c(unvax, first, second, third), dim=c(9,4)))
  
}
  
get_seas <- function(index_date, L){
  seasonality <- fread("seasonality_35.csv")
  i <- lubridate::yday(index_date)
  seas <- c(seasonality[day >= i, beta], rep(seasonality[, beta],ceiling(L/365)))
  seas <- seas/seas[1]
 return(seas[1:L]) 
}

#' @export
read_param_file_simplified <- function(param_file){
  par <- readxl::read_excel(param_file)
  
  di <- par$Value
  names(di) <- par$variable_name
  
  
  hosp_prob <- c(di[["p_hosp_1"]],di[["p_hosp_2"]],di[["p_hosp_3"]],di[["p_hosp_4"]],
                 di[["p_hosp_5"]],di[["p_hosp_6"]],di[["p_hosp_7"]],di[["p_hosp_8"]],
                 di[["p_hosp_9"]])
  ifr <- c(di[["p_death_1"]],di[["p_death_2"]],di[["p_death_3"]],di[["p_death_4"]],
           di[["p_death_5"]],di[["p_death_6"]],di[["p_death_7"]],di[["p_death_8"]],
           di[["p_death_9"]])
  icu_stay <- c(di[["icu_stay_1"]],di[["icu_stay_2"]],di[["icu_stay_3"]],
                di[["icu_stay_4"]],di[["icu_stay_5"]],di[["icu_stay_6"]],
                di[["icu_stay_7"]],di[["icu_stay_8"]],di[["icu_stay_9"]])
  pre_icu <- c(di[["pre_icu_1"]], di[["pre_icu_2"]], di[["pre_icu_3"]],
               di[["pre_icu_4"]], di[["pre_icu_5"]], di[["pre_icu_6"]],
               di[["pre_icu_7"]], di[["pre_icu_8"]], di[["pre_icu_9"]])
  asympt_prob <- c(di[["asympt_prob_1"]], di[["asympt_prob_2"]], di[["asympt_prob_3"]],
                   di[["asympt_prob_4"]], di[["asympt_prob_5"]], di[["asympt_prob_6"]],
                   di[["asympt_prob_7"]], di[["asympt_prob_8"]], di[["asympt_prob_9"]])
  health_workers <- c(di[["n_health_care_1"]], di[["n_health_care_2"]], di[["n_health_care_3"]], di[["n_health_care_4"]],
                      di[["n_health_care_5"]], di[["n_health_care_6"]], di[["n_health_care_7"]], di[["n_health_care_8"]],
                      di[["n_health_care_9"]])
  params <- list(
    latent_period = di[["latent_periode"]],
    pre_sympt_period = di[["pre_sympt_periode"]],
    infectious_period = di[["time_infectious"]],
    sympt_frac=1-asympt_prob,
    asympt_frac=asympt_prob,
    pre_sympt_infect=di[["infectiousness_pre_sympt"]],
    asympt_infect=di[["infectiousness_asympt"]],
    susceptibility=c(di[["sus_1"]], di[["sus_2"]], di[["sus_3"]],
                         di[["sus_4"]], di[["sus_5"]], di[["sus_6"]],
                         di[["sus_7"]], di[["sus_8"]], di[["sus_9"]]),
    transmisibility=c(di[["trans_1"]], di[["trans_2"]], di[["trans_3"]],
                         di[["trans_4"]], di[["trans_5"]], di[["trans_6"]],
                         di[["trans_7"]], di[["trans_8"]], di[["trans_9"]]),
    mixing_matrix=mixing_matrix,
    hosp_prob=hosp_prob,
    prob_death_non_hosp=ifr,
    prob_death_hosp=ifr,
    prob_death_icu=ifr,
    time_before_death=di[["time_to_death_I"]],
    time_before_death_hosp=di[["time_to_death_H"]],
    time_before_death_icu=di[["time_to_death_ICU"]],
    icu_prob=c(di[["icu_prob_1"]], di[["icu_prob_2"]], di[["icu_prob_3"]], di[["icu_prob_4"]],
               di[["icu_prob_5"]], di[["icu_prob_6"]], di[["icu_prob_7"]], di[["icu_prob_8"]],
               di[["icu_prob_9"]]),
    n_hosp_workers=health_workers,
    length_hosp=c(di[["hosp_stay_1"]], di[["hosp_stay_2"]], di[["hosp_stay_3"]],
                      di[["hosp_stay_4"]], di[["hosp_stay_5"]], di[["hosp_stay_6"]],
                  di[["hosp_stay_7"]], di[["hosp_stay_8"]], di[["hosp_stay_9"]]),
    pre_icu=pre_icu,
    length_icu=icu_stay,
    post_icu=c(di[["icu_hosp_1"]], di[["icu_hosp_2"]], di[["icu_hosp_3"]],
      di[["icu_hosp_4"]], di[["icu_hosp_5"]], di[["icu_hosp_6"]],
      di[["icu_hosp_7"]], di[["icu_hosp_8"]], di[["icu_hosp_9"]]) -pre_icu - icu_stay
    #waning_immunity=di["waning_immunity"]
  )

  params$post_icu[params$post_icu < 0] <- 1
  return(params)



}


read_param_file <- function(param_file){
  
  par <- readxl::read_excel(param_file)
  
  di <- par$Value
  names(di) <- par$variable_name
  

  if("risk_group_1" %in% names(di)){
    frac_risk_group <- c(di[["risk_group_1"]],di[["risk_group_2"]],di[["risk_group_3"]],
                         di[["risk_group_4"]],di[["risk_group_5"]],di[["risk_group_6"]],
                         di[["risk_group_7"]],di[["risk_group_8"]],di[["risk_group_9"]])
  }
  
  hosp_prob <- c(di[["p_hosp_1"]],di[["p_hosp_2"]],di[["p_hosp_3"]],di[["p_hosp_4"]],
                 di[["p_hosp_5"]],di[["p_hosp_6"]],di[["p_hosp_7"]],di[["p_hosp_8"]],
                 di[["p_hosp_9"]])
  
  
  ifr <- c(di[["p_death_1"]],di[["p_death_2"]],di[["p_death_3"]],di[["p_death_4"]],
           di[["p_death_5"]],di[["p_death_6"]],di[["p_death_7"]],di[["p_death_8"]],
           di[["p_death_9"]])
  
  


  icu_stay <- c(di[["icu_stay_1"]],di[["icu_stay_2"]],di[["icu_stay_3"]],
                di[["icu_stay_4"]],di[["icu_stay_5"]],di[["icu_stay_6"]],
                di[["icu_stay_7"]],di[["icu_stay_8"]],di[["icu_stay_9"]])

  pre_icu <- c(di[["pre_icu_1"]], di[["pre_icu_2"]], di[["pre_icu_3"]],
               di[["pre_icu_4"]], di[["pre_icu_5"]], di[["pre_icu_6"]],
               di[["pre_icu_7"]], di[["pre_icu_8"]], di[["pre_icu_9"]])
  
  
  asympt_prob <- c(di[["asympt_prob_1"]], di[["asympt_prob_2"]], di[["asympt_prob_3"]],
                   di[["asympt_prob_4"]], di[["asympt_prob_5"]], di[["asympt_prob_6"]],
                   di[["asympt_prob_7"]], di[["asympt_prob_8"]], di[["asympt_prob_9"]])
  health_workers <- c(di[["n_health_care_1"]], di[["n_health_care_2"]], di[["n_health_care_3"]], di[["n_health_care_4"]],
                      di[["n_health_care_5"]], di[["n_health_care_6"]], di[["n_health_care_7"]], di[["n_health_care_8"]],
                      di[["n_health_care_9"]])
  
  if("risk_group_1" %in% names(di)){
    health_workers <- round(rep(health_workers,2) *c(1-frac_risk_group, frac_risk_group))
  }
  params <- list(
    latent_period = di[["latent_periode"]],
    pre_sympt_period = di[["pre_sympt_periode"]],
    infectious_period = di[["time_infectious"]],
    sympt_frac=1-asympt_prob,
    asympt_frac=asympt_prob,
    pre_sympt_infect=di[["infectiousness_pre_sympt"]],
    asympt_infect=di[["infectiousness_asympt"]],
    susceptibility=c(di[["sus_1"]], di[["sus_2"]], di[["sus_3"]],
                         di[["sus_4"]], di[["sus_5"]], di[["sus_6"]],
                         di[["sus_7"]], di[["sus_8"]], di[["sus_9"]]),
    transmisibility=c(di[["trans_1"]], di[["trans_2"]], di[["trans_3"]],
                         di[["trans_4"]], di[["trans_5"]], di[["trans_6"]],
                         di[["trans_7"]], di[["trans_8"]], di[["trans_9"]]),
    mixing_matrix=mixing_matrix,
    hosp_prob=hosp_prob,
    prob_death_non_hosp=ifr,
    prob_death_hosp=ifr,
    prob_death_icu=ifr,
    time_before_death=di[["time_to_death_I"]],
    time_before_death_hosp=di[["time_to_death_H"]],
    time_before_death_icu=di[["time_to_death_ICU"]],
    ## prob_misc=rep(c(di[["p_misc_1"]], di[["p_misc_2"]], di[["p_misc_3"]], di[["p_misc_4"]],
    ##            di[["p_misc_5"]], di[["p_misc_6"]], di[["p_misc_7"]], di[["p_misc_8"]],
    ##            di[["p_misc_9"]]),2),
    icu_prob=c(di[["icu_prob_1"]], di[["icu_prob_2"]], di[["icu_prob_3"]], di[["icu_prob_4"]],
               di[["icu_prob_5"]], di[["icu_prob_6"]], di[["icu_prob_7"]], di[["icu_prob_8"]],
               di[["icu_prob_9"]]),
    n_hosp_workers=health_workers,
    length_hosp=c(di[["hosp_stay_1"]], di[["hosp_stay_2"]], di[["hosp_stay_3"]],
                      di[["hosp_stay_4"]], di[["hosp_stay_5"]], di[["hosp_stay_6"]],
                  di[["hosp_stay_7"]], di[["hosp_stay_8"]], di[["hosp_stay_9"]]),
    pre_icu=pre_icu,
    length_icu=icu_stay,
    frac_risk_group=frac_risk_group,
    post_icu=c(di[["icu_hosp_1"]], di[["icu_hosp_2"]], di[["icu_hosp_3"]],
      di[["icu_hosp_4"]], di[["icu_hosp_5"]], di[["icu_hosp_6"]],
      di[["icu_hosp_7"]], di[["icu_hosp_8"]], di[["icu_hosp_9"]]) -pre_icu - icu_stay
    #waning_immunity=di["waning_immunity"]
  )

  params$post_icu[params$post_icu < 0] <- 1
  return(params)
}



get_vac_params_old <- function(vaccine_file, dose_file){
  v_par <- fread( vaccine_file)
  di <- v_par$x
  names(di) <- v_par$V1
  

  VE_hosp_conditional_1 <- (1 - di[["VACCINE_EFF_1ST_1_S"]])/(1  -di[["VACCINE_EFF_1ST_1_I"]])
  VE_hosp_conditional_2 <- (1 - di[["VACCINE_EFF_2ND_1_S"]])/(1  -di[["VACCINE_EFF_2ND_1_I"]] )

  VE_death_conditional_1 <- (1 - di[["VACCINE_EFF_1ST_1_D"]])/(1  -di[["VACCINE_EFF_1ST_1_I"]] )
  VE_death_conditional_2 <- (1 - di[["VACCINE_EFF_2ND_1_D"]])/(1  -di[["VACCINE_EFF_2ND_1_I"]] )

  if(!is.null(dose_file) & ! dose_file==""){
    v_doses <- fread(dose_file)
  }else{
    v_doses = list(vac_1=rep(0, 500))
  }
  vac_pars <- list(rr_inf = c(1, 1-di[["VACCINE_EFF_1ST_1_I"]], 1-di[["VACCINE_EFF_2ND_1_I"]]),
                   rr_inf_asymp = c(1,1-di[["VACCINE_EFF_1ST_1_Ia"]],
                                    1- di[["VACCINE_EFF_2ND_1_Ia"]]),
                   rr_hosp = c(1, VE_hosp_conditional_1, VE_hosp_conditional_2),
                   rr_death = c(1,VE_death_conditional_1, VE_death_conditional_2),
                   rr_icu = c(1, 1,1),
                   rr_los_hosp = c(1, 1,1),
                   rr_trans = c(1,di[["VACCINE_EFF_1ST_1_V"]],di[["VACCINE_EFF_2ND_1_V"]]),
                   ramp_up_time=di[["RAMP_UP_1ST_1"]],
                   time_dose_2=di[["RAMP_UP_1ST_1"]] + di[["DELAY_EFFECT_1"]],
                   vac_doses= v_doses$vac1
                   )

  
  return(vac_pars)


}





create_regional_contact_matrix <- function(N_age, reg_data,params, file="mobility.RDS"){
  mob <- readRDS(file)
                                        # add municip1856 <-> number 80
  mob <- rbind(mob[1:79,], municip1856=0,mob[80:355,])
  mob <- cbind(mob[,1:79], municip1856=0,mob[,80:355])
  df <- reshape2::melt(mob, value.name="mobility")
  setDT(df)
  regs <- as.data.table(reg_data$reg_spec_municip)
  a <- regs %>% left_join(regs %>% left_join(df, by=c("fhidata.municip_code"="Var1")), by=c("fhidata.municip_code"="Var2"))
  new <- a %>% group_by(name.x, name.y) %>% summarize(V1=sum(mobility))
  setDT(new)
  wide <- data.table::dcast(new, name.x ~ name.y, value.var = "V1")
  names <- wide %>% pull(name.x)
  wide <- wide %>% select(!name.x)
  d <- as.matrix(wide)
  rownames(d) <- names
  diag(d) <- 0
  mob <- d
  #Leak factor - For each person who travels, how much of their "infectiousness leaks"
  leak_factor <- 1
  cont_matrix <- mob*leak_factor
  diag(cont_matrix) <- reg_data$reg_pop - colSums(cont_matrix)
  reg_contact_matrix <- sweep(cont_matrix, MARGIN=2, 1/reg_data$reg_pop, `*`)
  ## cont_matrix <- 
  
  ## regs <- unique(pop$location_code)
  ## ags <-  unique(levels(pop$group))

  ## if(N_age == 9){
  ##   mobility_age_factors <- c(0.009, 0.1, 0.149, 0.126, 0.153, 0.184, 0.145, 0.099, 0.036)
  ## }else{
  ##   mobility_age_factors <- rep(c(0.009, 0.1, 0.149, 0.126, 0.153, 0.184, 0.145, 0.099, 0.036), 2)
  ## }


  very_large_mixing_matrix <-matrix(0, ncol=reg_data$N_regions*18, nrow=reg_data$N_regions*18)
  for(a in 1:N_age){
    for(b in 1:N_age){
      for(i in 1:reg_data$N_regions){
        for(j in 1:reg_data$N_regions){
          very_large_mixing_matrix[(i-1)*N_age+a,(j-1)*N_age+b] <- reg_contact_matrix[i,j]*params$mixing_matrix[a,b]
        }
      }
    }
  }
  return(very_large_mixing_matrix)
}











get_regional_data <- function(regional_file, regions="prior_0"){

  reg_dat <- fread(regional_file)
  colnames(reg_dat) <- c("fhidata.municip_code", "fhidata.municip_name","value", "prior_0")

  regions = list()
  reg_dat <- reg_dat %>% mutate(fylke=substr(fhidata.municip_code, 8,9),
                                region=recode(fylke, "03"= "Ost",
                                              "30"= "Ost",
                                              "34"= "Ost",
                                              "38"= "Ost",
                                              "54"= "Nord",
                                              "18"= "Nord",
                                              "50"= "Trond",
                                              "15"= "Vest",
                                              "46"= "Vest",
                                              "11"= "Vest",
                                              "42"= "Agder"),
                                name:=paste(fylke, prior_0))

  pop_data <- spldata::nor_population_by_age_cats(cats=list("1"=-1:9, "2"=10:19, "3"=20:29, "4"=30:39,
                                                            "5"=40:49, "6"=50:59, "7"=60:69, "8"=70:79, "9"=80:120),
                                                  border=2020, include_total=FALSE) %>%
    filter(calyear==2022 & granularity_geo=="municip") %>%mutate(pop=pop_jan1_n)

  pop_data18p <- spldata::nor_population_by_age_cats(cats=list("1"=18:120),
                                                  border=2020, include_total=FALSE) %>%
    filter(calyear==2022 & granularity_geo=="municip") %>%mutate(pop=pop_jan1_n)

  
  pop <- reg_dat %>%
    right_join(pop_data, by=c("fhidata.municip_code"="location_code")) %>%
    group_by(name, age) %>% summarize(pop=sum(pop))


  pop18p <- reg_dat %>%
    right_join(pop_data18p,
               by=c("fhidata.municip_code"="location_code")) %>%
    group_by(name) %>% summarize(pop=sum(pop))

  
  f <- reg_dat %>% right_join(pop_data %>% group_by(location_code) %>% summarize(pop=sum(pop)), by=c("fhidata.municip_code"="location_code")) %>%
    group_by(name, prior_0) %>% summarize(municips=paste(fhidata.municip_code, collapse=","), R_value=sum(pop*value)/sum(pop)) %>%ungroup()
  f <- f %>% mutate(reg_number=1:nrow(f))
  

  return(list(reg_spec=f,
              reg_spec_municip=reg_dat,
              pop=pop$pop,
              pop_18p=pop18p %>%pull(pop),
              reg_pop=pop%>%group_by(name) %>% summarize(pop=sum(pop)) %>% pull(pop),
              age_pop=pop%>%group_by(age) %>% summarize(pop=sum(pop)) %>% pull(pop),
              N_regions=nrow(f)))
}


read_initial_conditions <- function(filename, params, reg_data){
  
  data <- fread(filename)
  data[grepl("ward", county), county:="municip0301"]

  comps <- c("E1"="Ea", "E1"="Es", "Ia"="A", "E2"="P", "I"="I", "prev_H"="H", "prev_ICU"="ICU_R", "R"="R")
  N <- rep(0, params$n)
  for(var in names(comps)){
    
    d1 <- data[compartment %in% paste0(var, "_", 1:9)] %>% left_join(reg_data$reg_spec_municip, by=c("county"="fhidata.municip_code")) %>% group_by(name, compartment) %>% summarize(value=round(sum(value.x)))
    
    vals <- rbind(d1 %>% mutate(risk_group=0, frac=1 - params$frac_risk_group), d1 %>% mutate(risk_group=1, frac=params$frac_risk_group)) %>% arrange(name, risk_group, compartment) %>% mutate(n=round(value*frac)) %>% pull(n)
    N <- N + vals
    params[[paste0(comps[var], "_ini")]][,1,1] <- vals
  }
  params$S_ini[,1] <-  params$reg_pop_long - N
  return(params)
}

add_import <- function(import_file,
                       age_dist_file,
                       N_age,
                       N_regions,
                       n_vac,
                       n_strain
                       ){
  imp_vec <- fread(import_file)$x
  age_dist <- fread(age_dist_file)$V1
  age_dist <- age_dist[1:9] + c(rep(0,8), age_dist[10])
  if(N_age==18){
    age_dist <- rep(age_dist, 2)/2
  }
  age_dist <- rep(age_dist, N_regions)/N_regions
  import_array <- matrix(NA, nrow=0, ncol=N_age*N_regions)
  for(i in 1:length(imp_vec)){
    import_array <- rbind(import_array, t(rmultinom(1, round(imp_vec[i]), age_dist)))
  }
  import_array <- outer(import_array, c(1, rep(0, n_vac-1)))
  import_array <- array(import_array, dim=c(length(imp_vec), N_age*N_regions, n_vac, n_strain))
  return(import_array)
}
                       

#' Create a parameter object based on a parameter file
#'
#' @keywords parameters
#' @export
get_params <- function(
                       param_file="input_files/parameters_vaccination.xlsx",
                       vaccine_file=NULL,
                       dose_file=NULL,
                       N_age=9,
                       L=500,
                       n_vac=3,
                       n_strain=1,
                       regions=NULL,
                       regional_file="input_files/regions.csv",
                       mobility_file="input_files/contact_matrix.csv",
                       vac_priority=NULL,
                       initial_conditions_file=NULL,
                       import_file=NULL,
                       import_age_dist_file=NULL,
                       add_plus_minus_neutral=FALSE
                       ){

  params <- read_param_file(param_file)

  if(!is.null(vaccine_file)){
    vac_pars <- get_vac_params_old(vaccine_file, dose_file)
  }

  if(!is.null(regions)){
    reg_data <- get_regional_data(regional_file)
    params$N_regions <- reg_data$N_regions
    
  }else{
    reg_data <- list(pop=get_age_groups(),
                     N_regions=1)
  }
  

  if(! N_age %in% c(9,18)){
    stop("N_age needs to be 9 or 18")
  }

  if(N_age==18){
    new_pop_wrong_order <- round(rep(reg_data$pop, 2)*c(rep(1-params$frac_risk_group,reg_data$N_regions), rep( params$frac_risk_group, reg_data$N_regions)))
    new_pop <- new_pop_wrong_order[rep((1:reg_data$N_regions)-1, each=18)*9 + rep(1:9, 2*reg_data$N_regions) + rep(c(rep(0, 9), rep(length(new_pop_wrong_order)/2, 9)), reg_data$N_regions)]
    reg_data$pop <- new_pop
    reg_data$age_pop <- round(rep(reg_data$age_pop, 2)*c(1-params$frac_risk_group, params$frac_risk_group))
    params$mixing_matrix <- cbind(rbind(mixing_matrix, mixing_matrix), rbind(mixing_matrix, mixing_matrix))
    params$merge_half_age <- TRUE
    params$transmisibility <- c(params$transmisibility, params$transmisibility)
    params$susceptibility <- c(params$susceptibility, params$susceptibility)
    params$hosp_prob <- rep(params$hosp_prob, 2)
    params$icu_prob <- rep(params$icu_prob, 2)
    params$length_hosp <- rep(params$length_hosp, 2)
    params$length_icu <- rep(params$length_icu, 2)
    params$pre_icu <- rep(params$pre_icu, 2)
    params$post_icu <- rep(params$post_icu, 2)
    params$asympt_frac <- rep(params$asympt_frac,2)
    params$sympt_frac <- rep(params$sympt_frac,2)
    params$prob_death_non_hosp <- rep(params$prob_death_non_hosp, 2)
    params$prob_death_hosp <- rep(params$prob_death_hosp, 2)
    params$prob_death_icu <- rep(params$prob_death_icu, 2)
  }
  
  N <- N_age*reg_data$N_regions

    params$transmisibility <- rep(params$transmisibility, reg_data$N_regions)
    params$susceptibility <- rep(params$susceptibility, reg_data$N_regions)
    params$hosp_prob <- rep(params$hosp_prob, reg_data$N_regions)
    params$icu_prob <- rep(params$icu_prob, reg_data$N_regions)
    params$length_hosp <- rep(params$length_hosp, reg_data$N_regions)
    params$length_icu <- rep(params$length_icu, reg_data$N_regions)
    params$pre_icu <- rep(params$pre_icu, reg_data$N_regions)
    params$post_icu <- rep(params$post_icu, reg_data$N_regions)
    params$asympt_frac <- rep(params$asympt_frac,reg_data$N_regions)
    params$sympt_frac <- rep(params$sympt_frac,reg_data$N_regions)
    params$prob_death_non_hosp <- rep(params$prob_death_non_hosp, reg_data$N_regions)
    params$prob_death_hosp <- rep(params$prob_death_hosp, reg_data$N_regions)
    params$prob_death_icu <- rep(params$prob_death_icu, reg_data$N_regions)

  
  if(reg_data$N_regions >1 ){
    params$mixing_matrix <- create_regional_contact_matrix(N_age, reg_data, params, file=mobility_file)
  }

  
  
  
  params <- c(params,
              list(
                N_steps=L,
                n_vac=n_vac,
                n_strain=n_strain,
                dt=1,
                T_waning=array(1e10, dim=c(N, n_vac)),
                vaccinations=array(0,dim=c(L, N, n_vac)),
                #    import_vec=rep(0, L),
                #    import_age_prio=rep(1,dim=(N,n_vac,1)),
                beta_day=matrix(rep(reg_data$reg_spec$R_value, each=N_age), ncol=N, nrow=L, byrow=TRUE),
                #    vac_time_full_effect=array(14.0, N),
                beta_strain=rep(1, n_strain),
                cross_protection=matrix(1, ncol=n_strain, nrow=n_strain),
                import_vec=array(0, dim=c(L,N, n_vac, n_strain)),
                n=N,
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
                S_ini=matrix(0, nrow=N, ncol=n_vac),
                I_ini=array(0, dim=c(N,n_vac,n_strain)),
                tot_hosp_ini=array(0, dim=c(N,n_vac,n_strain)),
                tot_resp_ini=array(0, dim=c(N,n_vac,n_strain)),
                tot_vac_ini=array(0, dim=c(N, n_vac)),
                tot_vac_adm_ini=array(0, dim=c(N, n_vac)),
                beta_norm=reg_data$pop,
                reg_pop_long=reg_data$pop,
                N_regions=reg_data$N_regions,
                waning_immunity_vax = array(1000, dim=c(N,n_vac,n_strain)),
                waning_inf = 10000,
                age_groups=N_age,
                beta_mode=1
              )
              )
  params$reg_pop18 <- reg_data$pop_18p
  

  params$age_pop <- reg_data$age_pop
  params$reg_pop <- reg_data$reg_pop
  
  if(!is.null(vac_priority)){
    adherence <- get_vac_adherence(vac_priority$adherence)
    if(vac_priority$type == "history"){
      vac <- get_vac_hist_sysvak(vac_priority$filename, reg_data, vac_priority$start_date, L, vac_adherence = adherence)
      params$vaccinations  <- vac_mat_to_vac_par(vac_pars, vac)[1:L,,]
    }else if(vac_priority$type == "scenario"){
      params$reg_prio = reg_data$reg_spec %>% filter(prior_0==1) %>%pull(reg_number)
      params$reg_prio_neutral = reg_data$reg_spec %>% filter(prior_0==0) %>%pull(reg_number)
      params$reg_prio_minus = reg_data$reg_spec %>% filter(prior_0==-1) %>%pull(reg_number)
      params$reg_pri_until <- 10
      params$reg_prio_amount <- vac_priority$amount
      params$start_day_prio <- vac_priority$start_day
      params$doses_per_day <- vac_pars$vac_doses
      params$already_vaccinated <- rep(0, N_age)
     
      params$adherence <- adherence$adherence
      params$adherence_hospital <- adherence$adherence_hosp
      vac <- create_vaccination_strategy_reg(params, vac_priority$file)$vac_1
      params$vaccinations  <- vac_mat_to_vac_par(vac_pars, vac)[1:L,,]
    }
  }
  if(add_plus_minus_neutral){
    params$reg_prio = reg_data$reg_spec %>% filter(prior_0==1) %>%pull(reg_number)
    params$reg_prio_neutral = reg_data$reg_spec %>% filter(prior_0==0) %>%pull(reg_number)
    params$reg_prio_minus = reg_data$reg_spec %>% filter(prior_0==-1) %>%pull(reg_number)
  }

  if(!is.null(import_file)){

    params$import_vec <- add_import(import_file, import_age_dist_file, N_age, params$N_regions, n_vac, n_strain)
  }
  
  params <- read_initial_conditions(initial_conditions_file, params, reg_data)
  
  params <- fix_params(params, N, n_vac, n_strain, vac_pars)

}


get_vac_adherence <- function(file){
  v_par <- fread(file)
  di <- v_par$x
  names(di) <- v_par$V1
  adherence <- rep(c(
    di[["ADHERENCE_1_0_11"]],
    (di[["ADHERENCE_1_12_15"]] +  di[["ADHERENCE_1_16_17"]]+  di[["ADHERENCE_1_18_24"]])/3,
    (di[["ADHERENCE_1_18_24"]] +di[["ADHERENCE_1_25_39"]])/2,
    di[["ADHERENCE_1_25_39"]],
    (di[["ADHERENCE_1_40_44"]] + di[["ADHERENCE_1_45_54"]])/2,
    (di[["ADHERENCE_1_45_54"]] + di[["ADHERENCE_1_55_64"]])/2,
    (di[["ADHERENCE_1_55_64"]] + di[["ADHERENCE_1_65_74"]])/2,
    (di[["ADHERENCE_1_65_74"]] + di[["ADHERENCE_1_75_84"]])/2,
    di[["ADHERENCE_1_75_84"]]),2)

    adherence_hosp <- c(
      di[["ADHERENCE_2_0_11"]],
      (di[["ADHERENCE_2_12_15"]] +  di[["ADHERENCE_2_16_17"]]+  di[["ADHERENCE_2_18_24"]])/3,
      (di[["ADHERENCE_2_18_24"]] +di[["ADHERENCE_2_25_39"]])/2,
      di[["ADHERENCE_2_25_39"]],
      (di[["ADHERENCE_2_40_44"]] + di[["ADHERENCE_2_45_54"]])/2,
      (di[["ADHERENCE_2_45_54"]] + di[["ADHERENCE_2_55_64"]])/2,
      (di[["ADHERENCE_2_55_64"]] + di[["ADHERENCE_2_65_74"]])/2,
      (di[["ADHERENCE_2_65_74"]] + di[["ADHERENCE_2_75_84"]])/2,
      di[["ADHERENCE_2_75_84"]])
  return(list(adherence=adherence,
         adherence_hosp=adherence_hosp))
}

vac_mat_to_vac_par <- function(vac_pars, vac){
  vac_delay <- round(rbind(matrix(0, nrow=round(vac_pars$ramp_up_time/2), ncol=dim(vac)[2]), vac)[1:dim(vac)[1],])
  vac_dose_2 <- round(rbind(matrix(0, nrow=vac_pars$ramp_up_time + vac_pars$time_dose_2, ncol=dim(vac)[2]), vac)[1:dim(vac)[1],])
  vaccinations=array(0,dim=c(dim(vac)[1], dim(vac)[2], 3))
  vaccinations[,,1] <- - vac_delay
  vaccinations[,,2] <-  vac_delay - vac_dose_2
  vaccinations[,,3] <- vac_dose_2
  return(vaccinations)
}
