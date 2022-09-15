get_vaccines_per_day <- function(params, tot_doses, key="vac_1_doses"){
  return(params$doses_per_day)
}

get_pop <- function(params, loc=NULL, offset=1){
  if(is.null(loc)) return(params$age_pop)
  return(params$reg_pop_long[((loc-1)*params$age_groups+1):(loc*params$age_groups)])
    
}
  

create_prio_groups <- function(priority_file){
  prio <- fread(priority_file)

  priority_list <- list()
  priorities <- c()
  groups <- c()
  ordered_groups <- list()
  population_fraction_list <- c()
  N_groups <- max(as.numeric(unique(unlist(prio))[5:length(unique(unlist(prio)))])) +1

  if(length(unique(unlist(prio)))==4) N_groups <- 1
  prio[, V1:=1:3]
  long_prio <- prio %>% tidyr::pivot_longer(!V1) %>% mutate(name2=as.numeric(lapply(name, function(X) as.numeric(strsplit(X, "_")[[1]][[1]])))) %>% mutate(age_group = if_else(name2 < 90, name2%/% 10 +1, 9)) %>% mutate(group_number=9*(V1-1) + age_group)
  
  for(i in 0:(N_groups-1)){
    group_numbers <- long_prio %>% filter(value==i) %>%pull(group_number)
    ordered_groups[[length(ordered_groups) +1]] <- list()
    pop_frac <- rep(0, 36)
    for(g in unique(group_numbers)){
      ordered_groups[[length(ordered_groups)]][[length(ordered_groups[[length(ordered_groups)]])+1]] <- g

      if(g == 9 | g==18 | g==27){
        pop_frac[g] <- sum(group_numbers==g) / 30
      }else{
        pop_frac[g] <- sum(group_numbers==g) / 10
      }
      
      if(g > 18){
        ordered_groups[[length(ordered_groups)]][[length(ordered_groups[[length(ordered_groups)]])+1]] <- g + 9
        if(g==27){
          pop_frac[g + 9] <- sum(group_numbers==g) / 30
        }else{
          pop_frac[g + 9] <- sum(group_numbers==g) / 10
        }
      }
      
    }
    population_fraction_list[[length(population_fraction_list) + 1]] <- pop_frac
    
  }
  l <- Reduce(`+`, population_fraction_list)
  population_fraction_list_full <- lapply(population_fraction_list, function(x) x/(l+1e-10))
  return(list(ordered_groups=ordered_groups, population_fraction_list=population_fraction_list,
              population_fraction_list_full=population_fraction_list_full))
}
  
calc_doses <- function(ordered_groups, population_fraction_list, pops, params, locs=NULL){
  if(is.null(locs))locs=1:params$N_regions
  tot_doses <- 0
  group_doses <- list()
  group_fractions <- list()
  doses_per_subgroup <- numeric(18)


  already_vaccinated <- params$already_vaccinated
  s <- rep(0,18)
  if(length(already_vaccinated)==params$age_groups*params$N_regions){
    for(i in locs){
      s <- s + already_vaccinated[ ((i-1)*18 +1):(i*18)]
    }
  }else{
    s <- already_vaccinated
  }
  reg_hosp_workers <- params$n_hosp_workers*sum(pops)/sum(params$reg_pop)
  wax_non_hc <- s - reg_hosp_workers
  wax_non_hc[wax_non_hc < 0] <- 0
  wax_hc <-  s - wax_non_hc
  already_vaccinated <- c(wax_non_hc, wax_hc)
                                        # print(glue::glue("Already vaccinated:{sum(already_vaccinated)}"))

  for(i in 1:length(ordered_groups)){
    t <- (c(pops - reg_hosp_workers, reg_hosp_workers) *
          c(params$adherence, params$adherence_hospital, params$adherence[10:18]) -already_vaccinated)*population_fraction_list[[i]]
    t[t < 0] <- 0
    group_pop <- sum(t)
    tot_doses <- tot_doses + group_pop

    group_doses[[length(group_doses) + 1]] <- group_pop

    doses_per_subgroup <- doses_per_subgroup + t

    if(group_pop <=0) group_pop <-1
    group_fractions[[length(group_fractions) +1]] <-t / group_pop
  }

  return(list(tot_doses=tot_doses,
              group_doses=group_doses,
              group_fractions=group_fractions,
              doses_per_subgroup=doses_per_subgroup
              ))

}



update_tot <- function(tot, subtract){

  new_tot <- tot - subtract
  new_tot[new_tot<0] <- 0
  return(new_tot)
}






vac_strat <- function(ordered_groups_1, group_doses_1, group_fractions_1,
                      population_fraction_list_1,
                      daily_doses_1, doses_subgroup_1, 
                      date_finished_group=9){

### Assume that doses_subgroup_1 includes everyone that should be vaccinated

  current_priority_1 <- 1
  m_1 <- matrix(0, ncol=0, nrow=18)
  day <- 0

  tot_doses <- doses_subgroup_1
  tot_doses_current <- tot_doses

  date_finished <- 1
  setEnd<- TRUE
  prev_pop_fracs_1 <- population_fraction_list_1[[1]]
  while( ncol(m_1) < length(daily_doses_1)){
    cp <- sum(tot_doses)
    skip_1 <- FALSE
    if(current_priority_1 <=length(ordered_groups_1)){
      current_group_1 <- ordered_groups_1[[current_priority_1]]
      if(current_priority_1 == 1){
        current_fraction_1 <- tot_doses*prev_pop_fracs_1
        current_fraction_1[! 1:36 %in% unlist(current_group_1)] <- 0
        current_fraction_1 <- current_fraction_1/(sum(current_fraction_1) + 1e-10)
        tot_doses_current <- tot_doses*prev_pop_fracs_1
      }
      
      current_pop_1 <- sum(tot_doses_current[unlist(current_group_1)])
    }else{
      current_group_1 <- 0
      current_fraction_1 <- 0
      current_pop_1 <-  0
      skip_1 <- TRUE
    }
    if(current_pop_1 <= 0){
       current_group_1 <- 0
       current_fraction_1 <- 0
       current_pop_1 <-  0
       skip_1 <- TRUE
    }

    
    vac_1 <- rep(0, 18)

    doses_1 <- daily_doses_1[day +1]

    subtract <- rep(0, 36)
    cur_subtract <- rep(0,36)
    while(current_pop_1 < doses_1 & current_priority_1 <= length(ordered_groups_1)){
      for( g in current_group_1){
        if(g < 19){
          vac_1[g] <- vac_1[g] + current_pop_1*current_fraction_1[g]
          subtract[g] <- subtract[g] + current_pop_1*current_fraction_1[g]
          cur_subtract[g] <- current_pop_1*current_fraction_1[g]
        }else{
          i <- g - 18
          vac_1[i] <- vac_1[i] + current_pop_1*current_fraction_1[g]
          subtract[g] <- subtract[g] + current_pop_1*current_fraction_1[g]
          cur_subtract[g] <- current_pop_1*current_fraction_1[g]
        }
      }
      current_priority_1 <- current_priority_1  + 1
      doses_1  <- doses_1 - current_pop_1
      if(current_priority_1 <= length(ordered_groups_1)){
        prev_pop_fracs_1 <- prev_pop_fracs_1 + population_fraction_list_1[[current_priority_1]]
        current_group_1 <- ordered_groups_1[[current_priority_1]]
        current_pop_1 <- sum((tot_doses*prev_pop_fracs_1)[unlist(current_group_1)])

        tot_doses_current <- tot_doses*prev_pop_fracs_1        
        current_fraction_1 <- (tot_doses - cur_subtract)*prev_pop_fracs_1
        current_fraction_1[! 1:36 %in% unlist(current_group_1)] <- 0
        current_fraction_1 <- current_fraction_1/(sum(current_fraction_1) + 1e-10)
        cur_subtract <- rep(0, 36)


      }else{
        skip_1 <- TRUE
      }
    }

    if(! skip_1){
      for(g in current_group_1){
        if(g < 19){
          vac_1[g] <- vac_1[g] + doses_1*current_fraction_1[g]
          subtract[g] <- subtract[g] + doses_1*current_fraction_1[g]
        }else{
          i <- g - 18
          vac_1[i] <-vac_1[i] + doses_1*current_fraction_1[g]
          subtract[g] <- subtract[g] + doses_1*current_fraction_1[g]
          
        }
      }
    }
    tot_doses <- update_tot(tot_doses, subtract)
    tot_doses_current <- update_tot(tot_doses_current, subtract)
    m_1 <- cbind(m_1, vac_1)

    if(current_priority_1 > date_finished_group & setEnd){
      date_finished <- day
      setEnd <- FALSE
    }

    day <- day + 1
  }
  return(list("vac_1"=m_1, "day_finished"=date_finished))
}




create_vaccination_strategy_3_vax  <- function(params, priority_file){


  grps <- create_prio_groups(priority_file_1)
  ordered_groups_1 <- grps$ordered_groups
  population_fraction_list_1 <- grps$population_fraction_list
  pops <- get_pop(params)
  doses <- calc_doses(ordered_groups_1, population_fraction_list_1, pops, params)
  tot_doses_1 <- doses$tot_doses
  group_doses_1 <- doses$group_doses
  group_fractions_1 <- doses$group_fractions
  doses_subgroup_1 <- doses$doses_per_subgroup


 
  
  daily_doses_1 <- get_vaccines_per_day(params, tot_doses_1, key="vac_1_doses")



  r_list <- vac_strat(ordered_groups_1, group_doses_1, group_fractions_1,
                      population_fraction_list_1,
                      daily_doses_1,
                      doses_subgroup_1)
                            
  m_1 <- r_list$vac_1
  if(ncol(m_1) == 0){
    m_1 <- matrix(0, ncol=2, nrow=nrow(m_1))
  }
  l <- list()
  for(i in 1:length(params$reg_pop)){
    l[[i]] <- m_1*params$reg_pop[i]/sum(params$reg_pop)
  }
  return(vac_1)
  
}



create_vaccination_strategy_reg  <- function(params, priority_file){
  grps <- create_prio_groups(priority_file)
  ordered_groups_1 <- grps$ordered_groups
  population_fraction_list_1 <- grps$population_fraction_list
  population_fraction_list_1_full <- grps$population_fraction_list_full
  pops <- get_pop(params)
  doses <- calc_doses(ordered_groups_1, population_fraction_list_1, pops, params)
  tot_doses_1 <- doses$tot_doses
  group_doses_1 <- doses$group_doses
  group_fractions_1 <- doses$group_fractions
  doses_subgroup_1 <- doses$doses_per_subgroup

  daily_doses_1 <- get_vaccines_per_day(params)
  should_be_prio <- 0

  for(i in 1:length(group_doses_1)){
    should_be_prio <- should_be_prio + group_doses_1[[i]]
    if(i == params$reg_pri_until){
      break
    }
  }
  
  
  
  prio_frac <- sum(params$reg_pop18[params$reg_prio])/sum(params$reg_pop18)
  prio_reg_pop18 <- prio_frac*should_be_prio
  L <- length(daily_doses_1)
  start_day_prio <- params$start_day_prio
  if(params$reg_prio_amount == 0){
    end_day_prio <- start_day_prio
  }else{
    before_prio_start <- cumsum(daily_doses_1)[start_day_prio]*prio_frac
    dose_c <-cumsum((daily_doses_1)[start_day_prio:L])*((1+params$reg_prio_amount)*prio_frac) + before_prio_start
    end_day_prio <- min((start_day_prio:L)[dose_c> prio_reg_pop18])
  }
  
  final_days <- c()
  for(reg in 1:params$N_regions){
    pops <- get_pop(params, loc=reg)
    doses <- calc_doses(ordered_groups_1, population_fraction_list_1, pops, params, locs=reg)
    tot_doses_1 <- doses$tot_doses
    f <- params$reg_pop18[reg]/sum(params$reg_pop18)
    if(reg %in% params$reg_prio){
      f_prio <- (1+params$reg_prio_amount)*params$reg_pop18[reg]/sum(params$reg_pop18)
    }else if(reg %in% params$reg_prio_neutral){
      f_prio <- f
    }else{
      f_prio <- (sum(params$reg_pop18[params$reg_prio_minus]) - params$reg_prio_amount*sum(params$reg_pop18[params$reg_prio]))/sum(params$reg_pop18[params$reg_prio_minus])*f
    }
    
    f_ar <- c(rep(f, start_day_prio -1),
              rep(f_prio, end_day_prio - start_day_prio + 1),
              rep(f, L - end_day_prio)
              )
    
    daily_doses_1_reg <-daily_doses_1*f_ar
    
    final_day <- min((1:L)[cumsum(daily_doses_1_reg)>tot_doses_1], na.rm=T)
    final_days <- c(final_days, final_day)
  }
  
  
  l1 <- list()
  d <- 0
  params$end_day_prio <- end_day_prio
  regions_remaining <- list()
  for(t in end_day_prio:L){
    regions_remaining[[t]] <- 1:params$N_regions
  }
  
  day_finished <- -1
  for(reg in order(final_days)){
    pops <- get_pop(params, loc=reg)
    doses <- calc_doses(ordered_groups_1, population_fraction_list_1, pops, params, locs=reg)
    tot_doses_1 <- doses$tot_doses
    group_doses_1 <- doses$group_doses
    group_fractions_1 <- doses$group_fractions
    doses_subgroup_1 <- doses$doses_per_subgroup
    f <- params$reg_pop18[reg]/sum(params$reg_pop18)
    if(reg %in% params$reg_prio){
      f_prio <- (1+params$reg_prio_amount)*params$reg_pop18[reg]/sum(params$reg_pop18)
    }else if(reg %in% params$reg_prio_neutral){
      f_prio <- f
    }else{
      f_prio <- (sum(params$reg_pop18[params$reg_prio_minus]) - params$reg_prio_amount*sum(params$reg_pop18[params$reg_prio]))/sum(params$reg_pop18[params$reg_prio_minus])*f
    }
    f_after <- rep(params$reg_pop18[reg]/sum(params$reg_pop18), L-end_day_prio)
    for(t in end_day_prio:(L-1)){
      f_after[t-end_day_prio + 1] <- params$reg_pop18[reg]/sum(params$reg_pop18[regions_remaining[[t]]])
    }
                                        #      print(f_after)
    f_ar <- c(rep(f, start_day_prio -1),
              rep(f_prio, end_day_prio - start_day_prio + 1),
              f_after
              )
    
    daily_doses_1_reg <-daily_doses_1*f_ar
    
    r_list <- vac_strat(ordered_groups_1, group_doses_1,
                              group_fractions_1, population_fraction_list_1_full,
                              daily_doses_1_reg, doses_subgroup_1)
    if(r_list$day_finished > day_finished) day_finished=r_list$day_finished
    
    
    m_1 <- r_list$vac_1
    if(ncol(m_1) == 0){
      m_1 <- matrix(0, ncol=2, nrow=nrow(m_1))
    }
    l1[[reg]] <- m_1
    
    possible_days <- (1:L)[colSums(m_1) == 0]
    end_day <- possible_days[1]
    if(is.na(end_day)) end_day <- end_day_prio
    if(end_day < end_day_prio){
      end_day <- possible_days[2]
    }
    if (end_day < end_day_prio | is.na(end_day)){
      end_day <- end_day_prio
    }
    for(t in end_day:(L-1)){
      regions_remaining[[t]] <- regions_remaining[[t]][regions_remaining[[t]] != reg]
    }
  }
  
  vac_1 <- aperm(abind(l1, along=3), c(3,1,2))
  
  vac1_mat <- matrix(0, nrow=dim(vac_1)[3], ncol=params$age_groups*params$N_regions)

  for(i in 1:dim(vac_1)[3]){
    vac1_mat[i, ] <- as.vector(t(vac_1[,,i]))
  }
  
  return(list("vac_1"=vac1_mat, day_finished=day_finished
              ))
  
}


create_vaccination_strategy_3_vax_reg_specify_prio <- function(params, priority_file_1,
                                                               priority_file_2, priority_file_3
                                                               ){

  grps <- create_prio_groups(priority_file_1)
  ordered_groups_1 <- grps$ordered_groups
  population_fraction_list_1 <- grps$population_fraction_list
  pops <- get_pop(params)
  doses <- calc_doses(ordered_groups_1, population_fraction_list_1, pops, params)
  tot_doses_1 <- doses$tot_doses
  group_doses_1 <- doses$group_doses
  group_fractions_1 <- doses$group_fractions
  doses_subgroup_1 <- doses$doses_per_subgroup
  grps <- create_prio_groups(priority_file_2)
  ordered_groups_2 <- grps$ordered_groups
  population_fraction_list_2 <- grps$population_fraction_list
  grps <- create_prio_groups(priority_file_3)
  ordered_groups_3 <- grps$ordered_groups
  population_fraction_list_3 <- grps$population_fraction_list
  pops <- get_pop(params)


  
  daily_doses_1 <- get_vaccines_per_day(params, tot_doses_1, key="vac_1_doses")
  daily_doses_2 <- get_vaccines_per_day(params, tot_doses_2, key="vac_2_doses")
  daily_doses_3 <- get_vaccines_per_day(params, tot_doses_3, key="vac_3_doses")

  
  l1 <- list()
  l2 <- list()
  l3 <- list()
  d <- 0

  day_finished <- -1
  for(reg in 1:params$N_regions){

    f <- params$vaccination_priorities[,reg]
  #  print(paste("vac", "reg"))
   # print(f)
    pops <- get_pop(params, loc=reg)
    doses <- calc_doses(ordered_groups_1, population_fraction_list_1, pops, params, locs=reg)
    tot_doses_1 <- doses$tot_doses
    group_doses_1 <- doses$group_doses
    group_fractions_1 <- doses$group_fractions
    doses_subgroup_1 <- doses$doses_per_subgroup
    doses <- calc_doses(ordered_groups_2, population_fraction_list_2, pops, params, locs=reg)
    tot_doses_2 <- doses$tot_doses
    group_doses_2 <- doses$group_doses
    group_fractions_2 <- doses$group_fractions
    doses_subgroup_2 <- doses$doses_per_subgroup
    doses <- calc_doses(ordered_groups_3, population_fraction_list_3, pops, params, locs=reg)
    tot_doses_3 <- doses$tot_doses
    group_doses_3 <- doses$group_doses
    group_fractions_3 <- doses$group_fractions
    doses_subgroup_3 <- doses$doses_per_subgroup
    daily_doses_1_reg <- f*daily_doses_1
    daily_doses_2_reg <- f*daily_doses_2
    daily_doses_3_reg <- f*daily_doses_3
    r_list <- vac_strat_3_vac(ordered_groups_1, group_doses_1,
                              group_fractions_1, population_fraction_list_1,
                              daily_doses_1_reg, doses_subgroup_1,
                              ordered_groups_2, group_doses_2,
                              group_fractions_2,population_fraction_list_2,
                              daily_doses_2_reg, doses_subgroup_2,
                              ordered_groups_3, group_doses_3,
                              group_fractions_3, population_fraction_list_3,
                              daily_doses_3_reg, doses_subgroup_3,
                              date_finished_group=params$day_finished_group)
    if(r_list$day_finished > day_finished) day_finished=r_list$day_finished
    m_1 <- r_list$vac_1
    if(ncol(m_1) == 0){
      m_1 <- matrix(0, ncol=2, nrow=nrow(m_1))
    }
    l1[[reg]] <- m_1
    m_2 <- r_list$vac_2
    if(ncol(m_2) == 0){
      m_2 <- matrix(0, ncol=2, nrow=nrow(m_2))
    }
    l2[[reg]] <- m_2
    m_3 <- r_list$vac_3
    if(ncol(m_3) == 0){
      m_3 <- matrix(0, ncol=3, nrow=nrow(m_3))
    }
    l3[[reg]] <- m_3
    #print(sum(m_1) + sum(m_2) + sum(m_3))
  }
  vac_1 <- aperm(abind(l1, along=3), c(3,1,2))
  vac_2 <- aperm(abind(l2, along=3), c(3,1,2))
  vac_3 <- aperm(abind(l3, along=3), c(3,1,2))


  
  return(list("vac_1"=vac_1, "vac_2"=vac_2, "vac_3"=vac_3, day_finished=day_finished
              ))
  
}



## test <- function(){


##   out <- list()
##   for(i in c(0,0.2, 0.4, 0.8, 1)){
##     should_be_prio <- 0
##     params$reg_prio_amount <- i
##     for(i in 1:length(group_doses_1)){
##       should_be_prio <- should_be_prio + group_doses_1[[i]]
##       if(i == params$reg_pri_until){
##         break
##       }
##     }
##     prio_reg_pop <- params$reg_pop[params$reg_prio]/sum(params$reg_pop)*should_be_prio
##     print(prio_reg_pop)
##     L <- length(daily_doses_1)
##     if(params$reg_prio_amount == 0){
##       end_day_prio <- 1
##     }else{
##       end_day_prio <- min((1:L)[cumsum(daily_doses_1 + daily_doses_2 + daily_doses_3)*params$reg_prio_amount> prio_reg_pop])
##     }
    
##     l1 <- list()
##     l2 <- list()
##     l3 <- list()
##     d <- 0
##     for(reg in 1:11){
##       f <- params$reg_pop[reg]/sum(params$reg_pop)
##       if(reg %in% params$reg_prio){
##         f_prio <- params$reg_prio_amount*params$reg_pop[reg]/sum(params$reg_pop[params$reg_prio]) + f*(1-params$reg_prio_amount)
##         print(f_prio)
##       }else{
##         f_prio <- f*(1-params$reg_prio_amount)
##       }
      
##       daily_doses_1_reg <- c(
##         f_prio * daily_doses_1[1:end_day_prio],
##         f*daily_doses_1[(end_day_prio+1):L])
##       daily_doses_2_reg <- c(
##         f_prio * daily_doses_2[1:end_day_prio],
##         f*daily_doses_2[(end_day_prio+1):L])
##       daily_doses_3_reg <- c(
##         f_prio * daily_doses_3[1:end_day_prio],
##         f*daily_doses_3[(end_day_prio+1):L])
##       r_list <- vac_strat_3_vac(ordered_groups_1, as.list(unlist(group_doses_1)*f),
##                                 group_fractions_1, population_fraction_list_1,
##                                 daily_doses_1_reg, doses_subgroup_1*f,
##                                 ordered_groups_2, as.list(unlist(group_doses_2)*f),
##                                 group_fractions_2,population_fraction_list_2,
##                                 daily_doses_2_reg, doses_subgroup_2*f,
##                                 ordered_groups_3, as.list(unlist(group_doses_3)*f),
##                                 group_fractions_3, population_fraction_list_3,
##                                 daily_doses_3_reg, doses_subgroup_3*f)
      
##       m_1 <- r_list$vac_1
##       if(ncol(m_1) == 0){
##         m_1 <- matrix(0, ncol=2, nrow=nrow(m_1))
##       }
##       l1[[reg]] <- m_1
##       m_2 <- r_list$vac_2
##       if(ncol(m_2) == 0){
##         m_2 <- matrix(0, ncol=2, nrow=nrow(m_2))
##       }
##       l2[[reg]] <- m_2
##       m_3 <- r_list$vac_3
##       if(ncol(m_3) == 0){
##         m_3 <- matrix(0, ncol=3, nrow=nrow(m_3))
##       }
##       l3[[reg]] <- m_3
##     }
##     vac_1 <- aperm(abind(l1, along=3), c(3,1,2))
##     vac_2 <- aperm(abind(l2, along=3), c(3,1,2))
##     vac_3 <- aperm(abind(l3, along=3), c(3,1,2))

##     v <- vac_1 + vac_2 + vac_3
##     out[[length(out) +1]] <- data.table(
##       amount =params$reg_prio_amount,
##       oslo_8 = v[1,8,],
##       oslo_8_c = cumsum(v[1,8,]),
##       next_8_c = cumsum(v[2,8,]),
##       oslo=cumsum(colSums(v[1, ,])),
##       next_s=cumsum(colSums(v[2, ,])),
##       t=1:L
##       )
##   }

##   df <- rbindlist(out)
##   df[, amount_f :=factor(amount)]
##   x <- df[amount==1]
##   plot(spline(x$t, y=x$oslo_8, method="natural"))
##   x$oslo_8)
##   ggplot(df) + geom_line(aes(x=t, y=oslo_8, color=amount_f))
## }
find_regions_to_vaccinate <- function(scenarios){

  prios <- list()
  i <- 1
  for(scenario in scenarios){
      params <- read_params(scenario$param_file, "Mixingsym10year.csv")
      params <- add_vaccine_params(params, scenario$profile_file, scenario$doses)
      params$mobility_scale <- scenario$mobility_scale
      params$equal_regions <- scenario$equal_regions
      params <- read_epi_scenario(params, scenario$epi_file, scenario$R)
      half_60_69 <- FALSE
      if(length(grep("elderly", scenario$name))> 0){
        half_60_69<- TRUE
        
      }
      params$vaccinated1 <- create_vaccination_strategy(params, scenario$prio_file, regional=scenario$geo=="regional", half_60_69=half_60_69)
      
      regs <- apply(params$vaccinated1, 1 ,sum) > 1

      print(paste(scenario$name, scenario$doses, sum(params$vaccinated1)))
      d <- fhidata::norway_locations_long_b2020[granularity_geo=="county"][regs, location_code]
      prios[[i]] <- list(name=scenario$name,
                         regions=paste(d, collapse=","))
      i <- i+1
                         
  }
  fwrite(rbindlist(prios), "prio_regions.csv")

}
