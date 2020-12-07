
get_rite_p <- function(dat = madat,
                       country = "all",
                       depvar = NULL,
                       cov = "",
                       with_N = "FALSE",
                       news = "both",
                       Q_good = FALSE,
                       exclude_councilors = "TRUE",
                       interactvar = NULL,
                       treat1and2 = FALSE,
                       weights = "TRUE",
                       contested_seats = "TRUE",
                       sims = NULL){
  
  require(doParallel)
  n_cores <- detectCores()
  ifelse(n_cores > 4, corei <- makeCluster(6), corei <- makeCluster(as.integer(n_cores/2)))
  registerDoParallel(corei)
  
  # from <- which(names(dat)== "1")
  # to <- which(names(dat)== as.character(sims))
  # cov <- cov
  
  dat$row_n <- 1:nrow(dat)
  perms <- dat[,c(paste0("iter", 1:sims), "row_n")]
  dat <- dat[,!names(dat) %in% grep("^iter", names(madat_iter), value = TRUE)] #this makes it faster for estimates function
  
  df <- estimates(dat = dat, exclude_councilors = exclude_councilors,
                  cov = cov, weights = weights, with_N = with_N,
                  news = news, Q_good = Q_good, depvar = depvar,
                  country = country, interactvar = interactvar,
                  treat1and2 = treat1and2, contested_seats = contested_seats,
                  return_df = TRUE, skip_prep = FALSE)
  
  b <- estimates(dat = df, exclude_councilors = exclude_councilors,
                 cov = cov, weights = weights, with_N = with_N,
                 news = news, Q_good = Q_good, depvar = depvar,
                 country = country, interactvar = interactvar,
                 treat1and2 = treat1and2, contested_seats = contested_seats,
                 return_df = FALSE, skip_prep = TRUE)
  
  perms <- perms[perms$row_n %in% df$row_n,] 
  
  b_ri <- NA
  
  #b_dist <- sapply(1:sims, function(j) {
  #print(j)
  
  b_ri <- foreach(i = 1:sims, .combine = 'cbind', .packages = c('lfe','dplyr'), .export = 'estimates') %dopar% {
    
    df$treat <- perms[[i]]
    
    #reassign inverse propensity weights for brazil based on sector treatment proportion
    #maybe unnecessary
    df <- df %>%
      group_by(fe) %>%
      mutate(n = n(), treat_prop = sum(treat==1)/n) %>%
      ungroup()
    
    if("brz" %in% unique(df$ctry)){
      
      df$inv_wts[df$ctry == "brz"] <- with(df[df$ctry == "brz",], ifelse(treat == 1, 1/treat_prop, 1/(1-treat_prop)))
      
      sum_ipw <- sum(df$inv_wts[df$ctry == "brz"], na.rm = T)
      sum_metawt_brz <- sum(df$wt[df$ctry=="brz"], na.rm = T)
      scale_wt <- sum_metawt_brz/sum_ipw
      df$wt[df$ctry=="brz"] <- df$inv_wts[df$ctry=="brz"] * scale_wt
    }
    
    # reassign inverse propensity weights for ug2 based on location sample size
    # i.e., villages with fewer than 15 subjects had assignment probability of .5, and otherwise paired with assignment prob .2 and .8
    # we infer the propensity of treatment based on share of treated within a village,
    # which due to sample size are between .75 and .84 for .8 treatment density and between .157 and .238 for .2 density.
    
    if("ug2" %in% unique(df$ctry)){
      
      df$inv_wts[df$ctry == "ug2"] <- with(df[df$ctry == "ug2",], ifelse(n<15, 1/.5,ifelse(treat_prop > 2/3,
                                                                                           treat*1/.8 + (1-treat)*1/.2,
                                                                                           treat*1/.2 + (1-treat)*1/.8)))
      
      sum_ipw2 <- sum(df$inv_wts[df$ctry=="ug2"], na.rm = T)
      sum_metawt_ug2 <- sum(df$wt[df$ctry=="ug2"], na.rm = T)
      scale_wt2 <- sum_metawt_ug2/sum_ipw2
      df$wt[df$ctry=="ug2"] <- df$inv_wts[df$ctry=="ug2"] * scale_wt2
    }
    
    # brazil weights if country-specific estimation
    if(country == "brz") df$wt <- df$inv_wts
    
    # uganda 2 weights if country-specific estimation
    if(country == "ug2") df$wt <- df$inv_wts2
    
    estimates(dat = df, country = country, depvar = depvar, cov = cov, with_N = with_N, news = news, exclude_councilors = exclude_councilors, interactvar = interactvar, weights = weights, skip_prep = TRUE)$coef[1]
  }

return(list(est = b,
            p = mean(abs(b_ri) >= abs(b$coefficients[1])),
            rite = b_ri))

on.exit(stopCluster(corei))
invisible(gc())
}


estimates <- function(dat = madat,
                      exclude_councilors = "TRUE", # ug2 chair sample if TRUE, pooled sample of chairs and councilors if FALSE
                      cov = "",                # list covariates in string united by '+' sign
                      weights = "TRUE",            # weight estimates
                      with_N = "FALSE",            # regress on Nij
                      news = "both",             # restrict sample to good/bad news
                      Q_good = FALSE,            # use median Q rule for good news
                      depvar = NULL,             # dependent variable
                      country = "all",            # restrict data set to country
                      interactvar = NULL,
                      treat1and2 = FALSE,        # whether to regress on treat1 and treat2 (private and public treatments)
                      contested_seats = "TRUE",    #if TRUE, exclude non-competitive races in uganda 2
                      return_df = FALSE,         # return prepared data.frame instead of model estimates (to run with function 'intercept')
                      skip_prep = FALSE
){
  df <- dat
  
  if(!skip_prep){
    
    if(depvar == "m1" & !treat1and2){
      df <- df %>% subset(!is.na(m1))
    }
    
    if(depvar == "m3" & !treat1and2){
      df <- df %>% subset(!is.na(m3))
    }
    
    #subset to contested seats in Uganda 2
    if(contested_seats=="TRUE" & "ug2" %in% unique(df$ctry)){
      df <- filter(df, ctry != "ug2" | (ctry == "ug2" & ug2_contested == 1))
    }
    
    # add council and chair fixed effects from ug2 where individual observations are counted twice
    df$councilor_dummy <- 0 # if individual obs were not duplicated for analysis to include both councilor and chair, create constant fixed effect variable to avoid error in felm functions below.
    df$councilor_dummy[df$ug2_councilor_dummy == 1] <- 1
    
    if(exclude_councilors=="TRUE"){
      df <- filter(df, councilor_dummy != 1)
    }
    
    #subset to good or bad news sample
    df$goodnews <- NA
    if(Q_good){
      df$goodnews <- df$Q_good
    }else{
      df$goodnews <- df$N_good
    }
    
    if(news == "good") df <- filter(df, goodnews==1)
    if(news == "bad")  df <- filter(df, goodnews==0)
    
    #prepare for modeling
    c <- length(unique(df$ctry))
    # imputing
    if(nchar(cov)>0){
      df <- df %>%
        group_by(ctry, fe) %>%
        mutate(m14i = ifelse(is.na(m14), mean(m14, na.rm=T), m14)) %>%
        mutate(m17i = ifelse(is.na(m17), mean(m17, na.rm=T), m17)) %>%
        mutate(m18i = ifelse(is.na(m18), mean(m18, na.rm=T), m18)) %>%
        mutate(m20i = ifelse(is.na(m20), mean(m20, na.rm=T), m20)) %>%
        mutate(m21i = ifelse(is.na(m21), mean(m21, na.rm=T), m21)) %>%
        mutate(m22i = ifelse(is.na(m22), mean(m22, na.rm=T), m22)) %>%
        mutate(m24i = ifelse(is.na(m24), mean(m24, na.rm=T), m24)) %>%
        mutate(m26i = ifelse(is.na(m26), mean(m26, na.rm=T), m26)) %>%
        mutate(m27i = ifelse(is.na(m27), mean(m27, na.rm=T), m27)) %>%
        ungroup(df)
    }
    # Mean centering/stdising
    if("data.frame" %in% class(df[,depvar])){
      df$y <- df[,depvar][[1]]
    }
    if("numeric" %in% class(df[,depvar]) | "integer" %in% class(df[,depvar])){
      df$y <- df[,depvar]
    }
    
    df <- df %>% drop_na(y, fe,cl)
    
    if(with_N=="TRUE"){
      df <- df %>% drop_na(N)
    }
    if(nchar(cov)>0){
      df <- df %>% drop_na(m14i, m17i, m18i, m20i, m21i, m22i, m24i, m26i, m27i)
    }
    
    if(!is.null(interactvar)){
      df <- df %>% drop_na(interactvar)
      if(interactvar == "m11"){
        df <- df %>% subset(ctry != "mex")
      }
      if(interactvar == "m23"){
        df <- df %>% subset(ctry %in% c("ben", "brz", "ug2", "ug1"))
      }
    }
    
    if(depvar == "m8"){ # restrict to mexico and benin sample
      df <- df %>% subset(ctry %in% c("mex", "ben"))
    }
    
    df <- df %>%
      # mutate(obs = n()) %>%
      group_by(ctry) %>%
      dplyr::mutate(ctryobs = n()) %>%
      dplyr::mutate(wt = 1/(c*ctryobs)) %>% # create weight
      ungroup()
    
    if(with_N=="TRUE"){
      df <- df %>% group_by(ctry) %>%
        mutate(N_temp = (N - mean(N, na.rm = T))/(sd(N, na.rm =T))) %>%
        ungroup(df)
      df$N_temp[df$ctry=="mex"] <- 0 # because when subsetting to good/bad news, N has constant value
      df$N <- df$N_temp
    }
    
    if(depvar == "m6"){ #standardize if regressing m6 because of scaling differences
      df <- df %>% group_by(ctry) %>%
        mutate(y = (y - mean(y, na.rm = T))/(sd(y, na.rm =T))) %>%
        ungroup(df)
    }
    
    if(nchar(cov)>0){
      df <- df %>% group_by(ctry) %>%
        mutate(m14i = m14i - mean(m14i, na.rm = T)) %>%
        mutate(m17i = m17i - mean(m17i, na.rm = T)) %>%
        mutate(m18i = m18i - mean(m18i, na.rm = T)) %>%
        mutate(m20i = m20i - mean(m20i, na.rm = T)) %>%
        mutate(m21i = m21i - mean(m21i, na.rm = T)) %>%
        mutate(m22i = m22i - mean(m22i, na.rm = T)) %>%
        mutate(m24i = m24i - mean(m24i, na.rm = T)) %>%
        mutate(m26i = m26i - mean(m26i, na.rm = T)) %>%
        mutate(m27i = m27i - mean(m27i, na.rm = T)) %>%
        ungroup(df)
    }
    
    if(!is.null(interactvar)){
      if("data.frame" %in% class(df[,interactvar])){
        df$interactvar <- df[,interactvar][[1]]
      }
      if("numeric" %in% class(df[,interactvar])){
        df$interactvar <- df[,interactvar]
      }
      df <- df %>% group_by(ctry) %>%
        mutate(interactvar = interactvar - mean(interactvar, na.rm = T)) %>%
        ungroup(df)
    }
    
    if("brz" %in% unique(df$ctry)){
      
      sum_ipw <- sum(df$inv_wts, na.rm = T)
      sum_metawt_brz <- sum(df$wt[df$ctry=="brz"], na.rm = T)
      scale_wt <- sum_metawt_brz/sum_ipw
      df$wt[df$ctry=="brz"] <- df$inv_wts[df$ctry=="brz"] * scale_wt
    }
    
    if("ug2" %in% unique(df$ctry)){
      
      sum_ipw2 <- sum(df$inv_wts2[df$ctry=="ug2"], na.rm = T)
      sum_metawt_ug2 <- sum(df$wt[df$ctry=="ug2"], na.rm = T)
      scale_wt2 <- sum_metawt_ug2/sum_ipw2
      df$wt[df$ctry=="ug2"] <- df$inv_wts2[df$ctry=="ug2"] * scale_wt2
    }
    
    if (news == "good") df$N[df$ctry == "brz"] <- 0
    
    #restrict to country sample if applicable
    if(country != "all"){
      df <- df %>% subset(ctry == country)
    }
    
    # brazil weights if country-specific estimation
    if(country != "all"){
      if(country == "brz"){
        df$wt <- df$inv_wts
      }}
    
    # uganda 2 weights if country-specific estimation
    if(country != "all"){
      if(country == "ug2"){
        df$wt <- df$inv_wts2
      }}
    
    if(return_df){
      return(df)
    }
  }
  
  if(!return_df){
    
    # specify model and whether is weighted
    
    if(with_N=="FALSE" & nchar(cov)==0){
      if(weights=="TRUE" & !is.null(interactvar)){
        model <- felm(as.formula(paste0("y ~ treat*", interactvar, "|fe + councilor_dummy|0|cl")), weights = df$wt, data = df) 
      }
      if(weights=="TRUE" & is.null(interactvar)){
        if(!treat1and2){
          model <- felm(y ~ treat|fe + councilor_dummy|0|cl, weights = df$wt, data = df)
        }else{
          model <- felm(y ~ treat1 + treat2|fe + councilor_dummy|0|cl, weights = df$wt, data = df)   
        }}
      if(weights=="FALSE" & is.null(interactvar)){
        if(!treat1and2){
          model <- felm(y ~ treat|fe + councilor_dummy|0|cl, data = df)
        }else{
          model <- felm(y ~ treat1 + treat2|fe + councilor_dummy|0|cl, data = df)
        }}
    }
    
    if(with_N=="TRUE" & nchar(cov)==0){
      if(weights=="TRUE"){
        model <- felm(y ~ treat*N|fe + councilor_dummy|0|cl, weights = df$wt, data = df)
      }else{
        model <- felm(y ~ treat*N|fe + councilor_dummy|0|cl, data = df)
      }
    }
    
    if(with_N=="FALSE" & nchar(cov)>0){
      if(weights=="TRUE"){
        model <- felm(as.formula(paste0("y ~ treat*(", cov, ")|fe + councilor_dummy|0|cl")), weights = df$wt, data = df)
      }else{
        model <- felm(as.formula(paste0("y ~ treat*(", cov, ")|fe + councilor_dummy|0|cl")), data = df)
      }
    }
    if(with_N=="TRUE" & nchar(cov)>0){
      if(weights=="TRUE"){
        model <- felm(as.formula(paste0("y ~ treat*(N+", cov, ")|fe + councilor_dummy|0|cl")), weights = df$wt, data = df)
      }else{
        model <- felm(as.formula(paste0("y ~ treat*(N+", cov, ")|fe + councilor_dummy|0|cl")), data = df)
      }
    }
    if(!return_df){
      return(model)
    }
  }
}

# Main results (returns estimates, RI p-value and control means) ----------

results <- function(ri_p = "generate",         # 'generate' runs RI code, 'ignore' generates NAs, 'load' loads existing RITE from file in 'fetch_rite_obj'
                    file_rite_obj = NULL,      # path to RITE object if 'ri_p = "load"'
                    # save_rite_obj = NULL,      # path to save RITE object
                    dat = madat,               # dataset (with permutations for randomization inference)
                    exclude_councilors = "TRUE", # ug2 chair sample if TRUE, pooled sample of chairs and councilors if FALSE
                    cov = "",                # list covariates in string united by '+' sign
                    weights = "TRUE",             # weight estimates
                    with_N = "FALSE",            # regress on Nij
                    news = "both",             # restrict sample to good/bad news
                    Q_good = FALSE,            # use median Q rule for good news
                    depvar = NULL,             # dependent variable
                    country = "all",           # restrict data set to country
                    interactvar = NULL,        # string name of variable interacted with treatment
                    treat1and2 = FALSE,        # TRUE if public vs private analysis
                    contested_seats = "TRUE",    # if TRUE, exclude non-competitive races in uganda 2
                    return_df = FALSE,         # return prepared data.frame instead of model estimates (to run with function 'intercept')
                    skip_prep = FALSE,         # perform estimations only (skips subsetting and weight calculation)
                    sims = NULL){             # number of simulations of treatment assignment
  
  
  if(ri_p == "ignore") ignore_perms <- TRUE
  if(ri_p == "generate"){
    
    if(sum(names(dat)== "iter1")>0){
      ignore_perms <- FALSE
      
      ri <- get_rite_p(dat = dat, exclude_councilors = exclude_councilors,
                       cov = cov, weights = weights, with_N = with_N,
                       news = news, Q_good = Q_good, depvar = depvar,
                       country = country, interactvar = interactvar,
                       treat1and2 = treat1and2, contested_seats = contested_seats,
                       sims = sims)
      
      est <- ri$est
      b_ri <- ri$rite
      p <- ri$p
      
      if(!is.null(file_rite_obj)){
        save(ri, file = file_rite_obj)
      }
      
    }else{
      ignore_perms <- TRUE
      warning("data object does not have treatment permutation columns,
              will return NA for randomization inference treatment effects and p-value.")
    }
    }
  
  if(ri_p == "load"){
    ignore_perms <- FALSE
    
    load_obj <- function(obj){
      env <- new.env()
      ri <- load(obj, env)[1]
      env[[ri]]
    }
    
    rite <- load_obj(file_rite_obj)
    
    est <- rite$est
    b_ri <- rite$rite
    p <- rite$p
    
  }
  
  if(!skip_prep){
    dat <- dat[,!names(dat) %in% grep("^iter", names(dat), value = TRUE)]
    df <- estimates(dat = dat, exclude_councilors = exclude_councilors,
                    cov = cov, weights = weights, with_N = with_N,
                    news = news, Q_good = Q_good, depvar = depvar,
                    country = country, interactvar = interactvar,
                    treat1and2 = treat1and2, contested_seats = contested_seats,
                    return_df = TRUE, skip_prep = skip_prep)
  }else{
    df <- dat
  }
  
  
  df$y <- df[[depvar]]
  if(!treat1and2){
    mean.control <- round(lm(y ~ treat, weights = df$wt, data = df)$coefficients[1], 3)
  }else{
    mean.control <- round(lm(y ~ treat1 + treat2, weights = df$wt, data = df)$coefficients[1], 3)
  }
  
  if(ri_p == "ignore" | ignore_perms){
    
    est <- estimates(dat = df, exclude_councilors = exclude_councilors,
                     cov = cov, weights = weights, with_N = with_N,
                     news = news, Q_good = Q_good, depvar = depvar,
                     country = country, interactvar = interactvar,
                     treat1and2 = treat1and2, contested_seats = contested_seats,
                     return_df = FALSE, skip_prep = TRUE)
    b_ri <- NA
    p <- NA
  }
  
  return(list(estimates = est,
              mean.control = mean.control,
              rite = b_ri,
              p = p
  ))
  
  gc()
  }