# ------------------------------
# Libraries
source("used_libraries.R")

# ------------------------------
# Superlearning procedures

# Learning algorithms employed to learn treatment/outcome mechanisms -----------

# Mean model
lrnr_me = make_learner(Lrnr_mean) 
# GLM
lrnr_lm = make_learner(Lrnr_glm_fast)
# Spline regression
lrnr_sp = make_learner(Lrnr_earth)  

# Meta-learners: to stack together predictions from the learners ---------------

# Combine continuous predictions with non-negative least squares
meta_C = make_learner(Lrnr_nnls)   

# Combine binary predictions with logit likelihood (augmented Lagrange optimizer)
meta_B = make_learner(Lrnr_solnp,                     
  loss_function = loss_loglik_binomial,              
  learner_function = metalearner_logistic_binomial)   

# Super-learners: learners + meta-learners together ----------------------------

# Continuous super-learning
super_C = Lrnr_sl$new(learners = list(lrnr_me, lrnr_lm, lrnr_sp), 
                      metalearner = meta_C)
# Binary super-learning
super_B = Lrnr_sl$new(learners = list(lrnr_me, lrnr_lm, lrnr_sp), 
                      metalearner = meta_B)
# Super-learners put together
super_list = list(A = super_B,
                  Y = super_C)

# -----------------------------------
# Estimators

################################################################################
# TMLE-CC
################################################################################
tmlecc_estimator = function(Q1pre='',PSmis=''){
  
  if(Q1pre=='' & PSmis==''){
    # Default TMLE package
    tmle.pret = tmle3(
      tmle_spec = tmle_ATE(1,0),                     # Targeting the ATE
      node_list = list(W = 'W', A = 'A', Y = 'Y'),   # Variables involved                         
      data = full.data[sel,],                        # Data                                  
      learner_list = super_list)                     # Super-learners  
    
    # Save estimate and CI
    tmle.ptest = c(tmle.pret$summary$lower,
                   tmle.pret$summary$upper,
                   tmle.pret$summary$tmle_est)
  } 
  else if(Q1pre!='' & PSmis==''){
    # Model for treatment assignment
    train.A = make_sl3_Task(
      data = full.data[sel,], outcome = 'A', covariates = 'W')
    
    pred.A = make_sl3_Task(
      data = full.data, outcome = 'A', covariates = 'W')
    
    A_fit = super_B$train(task = train.A)
    full.data = full.data[, A_pre := A_fit$predict(task = pred.A)]
    
    # Model for Q1: MISSPECIFIED
    train.Q1 = lm(as.formula(Q1pre), data = full.data[sel,])
    
    full.data = full.data[, Q1 := predict(train.Q1, newdata = full.data)] %>%
      .[, clever.H1 := ((A/A_pre)-(1-A)/(1-A_pre))] 
    
    # Define fluctuation model
    fluct.model.1 = lm(Y ~ -1 + offset(Q1) + clever.H1, data=full.data[sel,])
    
    # Auxiliary data frame for prediction
    temp.dt = copy(full.data)
    
    # Using estimated fluctuation parameter, update Q1.1
    temp.dt$A = 1
    
    full.data = full.data[, Q1.1 := predict(train.Q1, newdata = temp.dt)] %>%
      .[, up.Q1.1 := predict(fluct.model.1, 
                            newdata = data.frame(Q1=full.data$Q1.1,
                                                 clever.H1=(1/A_pre)))]
    
    # Using estimated fluctuation parameter, update Q1.0
    temp.dt$A = 0
    
    full.data = full.data[, Q1.0 := predict(train.Q1, newdata = temp.dt)] %>%
      .[, up.Q1.0 := predict(fluct.model.1, 
                             newdata = data.frame(Q1=full.data$Q1.0,
                                                  clever.H1=(-1/(1-A_pre))))]
    
    # Compute the updated difference Q1.1 - Q1.0
    # Compute the observed value for up.Q1
    # Compute the value of the efficient influence function
    full.data = full.data[, delta.up.Q1 := up.Q1.1-up.Q1.0] %>%
      .[, up.Q1.A := A*up.Q1.1 + (1-A)*up.Q1.0] %>%
      .[, EIF := delta.up.Q1 + clever.H1*(Y - up.Q1.A)]
    
    # Using the EIF, compute the asymptotic error
    asymp.sd = sd(full.data[sel,EIF])/sqrt(nrow(full.data[sel,]))
    
    # Save estimate and CI
    tmle.ptest = c(mean(full.data[sel,delta.up.Q1])-qnorm(0.975)*asymp.sd,
                   mean(full.data[sel,delta.up.Q1])+qnorm(0.975)*asymp.sd,
                   mean(full.data[sel,delta.up.Q1]))
  }
  else if(Q1pre=='' & PSmis!=''){
    # Model for treatment assignment
    train.A = glm(as.formula(PSmis), full.data[sel,], family=binomial('logit'))
      
    full.data = full.data[, A_pre := predict(train.A, type='response', newdata = full.data)]
    
    # Model for Q1
    train.Q1 = make_sl3_Task(
      data = full.data[sel,], outcome = 'Y', covariates = c('W','A'))
    
    Q1_fit = super_C$train(task = train.Q1)
    
    # Prediction task for Q1
    pred.Q1 = make_sl3_Task(
      data = full.data, outcome = 'Y', covariates = c('W','A'))
    
    full.data = full.data[, Q1 := Q1_fit$predict(task = pred.Q1)] %>%
      .[, clever.H1 := ((A/A_pre)-(1-A)/(1-A_pre))] 
    
    # Define fluctuation model
    fluct.model.1 = lm(Y ~ -1 + offset(Q1) + clever.H1, data=full.data[sel,])
    
    # Auxiliary data frame for prediction
    temp.dt = copy(full.data)
    
    
    temp.dt$A = 1
    pred.Q1 = make_sl3_Task(
      data = temp.dt, outcome = 'Y', covariates = c('W','A'))
    
    full.data = full.data[, Q1.1 := Q1_fit$predict(task = pred.Q1)] %>%
      .[, up.Q1.1 := predict(fluct.model.1, 
                            newdata = data.frame(Q1=full.data$Q1.1,
                                                 clever.H1=(1/A_pre)))]
    
    # Using estimated fluctuation parameter, update Q1.0
    temp.dt$A = 0
    pred.Q1 = make_sl3_Task(
      data = temp.dt, outcome = 'Y', covariates = c('W','A'))
    
    full.data = full.data[, Q1.0 := Q1_fit$predict(task = pred.Q1)] %>%
      .[, up.Q1.0 := predict(fluct.model.1, 
                            newdata = data.frame(Q1=full.data$Q1.0,
                                                 clever.H1=(-1/(1-A_pre))))]
    
    # Compute the updated difference Q1.1 - Q1.0
    # Compute the observed value for up.Q1
    # Compute the value of the efficient influence function
    full.data = full.data[, delta.up.Q1 := up.Q1.1-up.Q1.0] %>%
      .[, up.Q1.A := A*up.Q1.1 + (1-A)*up.Q1.0] %>%
      .[, EIF := delta.up.Q1 + clever.H1*(Y - up.Q1.A)]
    
    # Using the EIF, compute the asymptotic error
    asymp.sd = sd(full.data[sel,EIF])/sqrt(nrow(full.data[sel,]))
    
    # Save estimate and CI
    tmle.ptest = c(mean(full.data[sel,delta.up.Q1])-qnorm(0.975)*asymp.sd,
                   mean(full.data[sel,delta.up.Q1])+qnorm(0.975)*asymp.sd,
                   mean(full.data[sel,delta.up.Q1]))
  }
  return(tmle.ptest)
}

################################################################################
# DW
################################################################################
tmledw_estimator = function(Rpos='',PSmis=''){
  if(PSmis==''){
    # Model for treatment assignment
    train.A = make_sl3_Task(
      data = full.data, outcome = 'A', covariates = 'W')
    
    A_fit = super_B$train(task = train.A)
    full.data = full.data[, A_pre := A_fit$predict(task = train.A)]
  }
  else{
    # Model for treatment assignment
    train.A = glm(as.formula(PSmis), full.data, family=binomial('logit'))
      
    full.data = full.data[, A_pre := predict(train.A, type='response')]
  }
  
  if(Rpos==''){
    # Model for missingness mechanism
    train.R = make_sl3_Task(
      data = full.data, outcome = 'R0', covariates = c('W','A','M','Z'))
    
    R_fit = super_B$train(task = train.R)
    full.data = full.data[, R_pre := R_fit$predict(task = train.R)]
  }
  else{
    train.R = glm(as.formula(Rpos), full.data, family=binomial('logit'))
    
    full.data = full.data[, R_pre := predict(train.R, type='response')]
  }

  full.data = full.data[, DIW := (1/R_pre)*((A/A_pre)+(1-A)/(1-A_pre))]
  
  diw.mod = lm_robust(Y~A, full.data[sel,], weights = DIW)
  
  # Save estimate and CI
  diw.est = unname(c(diw.mod$conf.low[2],
                     diw.mod$conf.high[2],
                     diw.mod$coefficients[2]))
  return(diw.est)
}

################################################################################
# SR
################################################################################
sr_estimator = function(Q1pos='',BS=30,Q2mis=''){
  if(Q1pos==''){
    # Model for Q1
    train.Q1 = make_sl3_Task(
      data = full.data[sel,], outcome = 'Y', covariates = c('W','A','M','Z'))
    
    Q1_fit = super_C$train(task = train.Q1)
    
    # Prediction task for Q1
    pred.Q1 = make_sl3_Task(
      data = full.data, outcome = 'Y', covariates = c('W','A','M','Z'))
    
    full.data = full.data[, Q1 := Q1_fit$predict(task = pred.Q1)]
  } 
  else {
    # Model for Q1: MISSPECIFIED
    train.Q1 = lm(as.formula(Q1pos), data = full.data[sel,])
    
    full.data = full.data[, Q1 := predict(train.Q1, newdata = full.data)]
  }
  
  if(Q2mis==''){
    # Model for Q2.1 
    train.Q2.1 = make_sl3_Task(
      data = full.data[A==1,], outcome = 'Q1', covariates = c('W'))
    
    Q2.1_fit = super_C$train(task = train.Q2.1)
    
    pred.Q2.1 = make_sl3_Task(
      data = full.data, outcome = 'Q1', covariates = c('W'))
    
    full.data = full.data[, Q2.1 := Q2.1_fit$predict(task = pred.Q2.1)]
    
    # Model for Q2.0
    train.Q2.0 = make_sl3_Task(
      data = full.data[A==0,], outcome = 'Q1', covariates = c('W'))
    
    Q2.0_fit = super_C$train(task = train.Q2.0)
    
    pred.Q2.0 = make_sl3_Task(
      data = full.data, outcome = 'Q1', covariates = c('W'))
    
    full.data = full.data[, Q2.0 := Q2.0_fit$predict(task = pred.Q2.0)] %>%
      .[, delta := full.data$Q2.1-full.data$Q2.0]
  }
  else {
    # Model for Q2 misspecified
    train.Q2 = lm(as.formula(Q2mis), full.data)
    
    full.data = full.data[, Q2.1 := predict(train.Q2, newdata=data.frame(W=full.data$W,A=1))] %>%
      .[, Q2.0 := predict(train.Q2, newdata=data.frame(W=full.data$W,A=0))] %>%
      .[, delta := full.data$Q2.1-full.data$Q2.0]
  }
  
  # Bootstrap procedure to compute the standard deviation
  bsamples = c()
  for(j in 1:BS){
    
    # Boostraped data
    ind = sample(1:N,N,replace=T)
    bs.data = copy(full.data[ind,])
    
    if(Q1pos==''){
      # Model for Q1
      train.bs.Q1 = make_sl3_Task(
        data = bs.data[sel,], outcome = 'Y', covariates = c('W','A','M','Z'))
      
      bs.Q1_fit = super_C$train(task = train.bs.Q1)
      
      # Prediction task for Q1
      pred.bs.Q1 = make_sl3_Task(
        data = bs.data, outcome = 'Y', covariates = c('W','A','M','Z'))
      
      bs.data$Q1 = bs.Q1_fit$predict(task = pred.bs.Q1)
    }
    else {
      # Model for Q1: MISSPECIFIED
      train.bs.Q1 = lm(as.formula(Q1pos), data = bs.data[sel,])
      
      bs.data$Q1 = predict(train.bs.Q1, newdata = bs.data)
    }
    
    if(Q2mis==''){
      # Model for Q2.1.
      train.bs.Q2.1 = make_sl3_Task(
        data = bs.data[A==1,], outcome = 'Q1', covariates = c('W'))
      
      bs.Q2.1_fit = super_C$train(task = train.bs.Q2.1)
      
      pred.bs.Q2.1 = make_sl3_Task(
        data = bs.data, outcome = 'Q1', covariates = c('W'))
      
      bs.data$Q2.1 = bs.Q2.1_fit$predict(task = pred.bs.Q2.1)
      
      # Model for Q2.0
      train.bs.Q2.0 = make_sl3_Task(
        data = bs.data[A==0,], outcome = 'Q1', covariates = c('W'))
      
      bsQ2.0_fit = super_C$train(task = train.bs.Q2.0)
      
      pred.bs.Q2.0 = make_sl3_Task(
        data = bs.data, outcome = 'Q1', covariates = c('W'))
      
      bs.data$Q2.0 = Q2.0_fit$predict(task = pred.bs.Q2.0)
    }
    else {
      # Model for Q2
      train.bs.Q2 = lm(as.formula(Q2mis), bs.data)
      
      full.data$Q2.1 = predict(train.bs.Q2, newdata=data.frame(W=bs.data$W,A=1))
      
      full.data$Q2.0 = predict(train.bs.Q2, newdata=data.frame(W=bs.data$W,A=0))
    }
    
    # Add predicted difference Q2.1 - Q2.0 to boostrap vessel
    bsamples = c(bsamples, mean(bs.data$Q2.1-bs.data$Q2.0))
  }
  
  # Save estimate and CI
  nesreg.T = c(as.numeric(mean(full.data$delta) 
                          - (quantile(bsamples, 0.975)-quantile(bsamples, 0.025))/2),
               as.numeric(mean(full.data$delta) 
                          + (quantile(bsamples, 0.975)-quantile(bsamples, 0.025))/2),
               mean(full.data$delta))
}

################################################################################
# TSR
################################################################################
tsr_estimator = function(Q1pos='',Q2mis=''){
  
  # STEP 1 ------------------
  if(Q1pos==''){
    # Model for Q1
    train.Q1 = make_sl3_Task(
      data = full.data[sel,], outcome = 'Y', covariates = c('W','A','M','Z'))
    
    Q1_fit = super_C$train(task = train.Q1)
    
    # Prediction task for Q1
    pred.Q1 = make_sl3_Task(
      data = full.data, outcome = 'Y', covariates = c('W','A','M','Z'))
    
    full.data = full.data[, Q1 := Q1_fit$predict(task = pred.Q1)]
  } 
  else {
    # Model for Q1: MISSPECIFIED
    train.Q1 = lm(as.formula(Q1pos), data = full.data[sel,])
    
    full.data = full.data[, Q1 := predict(train.Q1, newdata = full.data)]
  }
  
  # Define clever variable
  full.data = full.data[, clever.H1 := (1/R_pre)*((A/A_pre)-(1-A)/(1-A_pre))] 
  
  # Define fluctuation model
  fluct.model.1 = lm(Y ~ -1 + offset(Q1) + clever.H1, data=full.data[sel,])
  
  # Auxiliary data frame for prediction
  temp.dt = copy(full.data)
  
  if(Q1pos==''){
    # Using estimated fluctuation parameter, update Q1.1
    temp.dt$A = 1
    pred.Q1 = make_sl3_Task(
      data = temp.dt, outcome = 'Y', covariates = c('W','A','M','Z'))
    
    full.data = full.data[, Q1.1 := Q1_fit$predict(task = pred.Q1)] %>%
      .[, up.Q1.1 := predict(fluct.model.1, 
                            newdata = data.frame(Q1=full.data$Q1.1,
                                                 clever.H1=(1/R_pre)*(1/A_pre)))]
    
    # Using estimated fluctuation parameter, update Q1.0
    temp.dt$A = 0
    pred.Q1 = make_sl3_Task(
      data = temp.dt, outcome = 'Y', covariates = c('W','A','M','Z'))
    
    full.data = full.data[, Q1.0 := Q1_fit$predict(task = pred.Q1)] %>%
      .[, up.Q1.0 := predict(fluct.model.1, 
                            newdata = data.frame(Q1=full.data$Q1.0,
                                                 clever.H1=(1/R_pre)*(-1/(1-A_pre))))]
  }
  else {
    # Using estimated fluctuation parameter, update Q1.1
    temp.dt$A = 1
    
    full.data = full.data[, Q1.1 := predict(train.Q1, newdata = temp.dt)] %>%
      .[, up.Q1.1 := predict(fluct.model.1, 
                            newdata = data.frame(Q1=full.data$Q1.1,
                                                 clever.H1=(1/R_pre)*(1/A_pre)))]
    
    # Using estimated fluctuation parameter, update Q1.0
    temp.dt$A = 0
    
    full.data = full.data[, Q1.0 := predict(train.Q1, newdata = temp.dt)] %>%
      .[, up.Q1.0 := predict(fluct.model.1, 
                            newdata = data.frame(Q1=full.data$Q1.0,
                                                 clever.H1=(1/R_pre)*(-1/(1-A_pre))))]
  }
  
  if(Q2mis==''){
    # Learn Q2.1 from up.Q1.1, using A=1 cases
    temp.dt = copy(full.data[A==1,])
    
    train.Q2 = make_sl3_Task(
      data = temp.dt, outcome = 'up.Q1.1', covariates = c('W'))
    
    Q2_fit = super_C$train(task = train.Q2)
    
    pred.Q2 = make_sl3_Task(
      data = full.data, outcome = 'up.Q1.1', covariates = c('W'))
    
    full.data = full.data[, Q2.1 := Q2_fit$predict(task = pred.Q2)]
    
    # Learn Q2.0 from up.Q1.0, using A=0 cases
    temp.dt = copy(full.data[A==0,])
    
    train.Q2 = make_sl3_Task(
      data = temp.dt, outcome = 'up.Q1.0', covariates = c('W'))
    
    Q2_fit = super_C$train(task = train.Q2)
    
    pred.Q2 = make_sl3_Task(
      data = full.data, outcome = 'up.Q1.0', covariates = c('W'))
    
    full.data = full.data[, Q2.0 := Q2_fit$predict(task = pred.Q2)]
    
    # STEP 2 ------------------
    
    # Compute the observed values for up.Q1 and Q2
    # Define clever variable
    full.data = full.data[, up.Q1.A := A*up.Q1.1 + (1-A)*up.Q1.0 ] %>% 
      .[, Q2.A := A*Q2.1 + (1-A)*Q2.0 ] %>% 
      .[, clever.H2 := ((A/A_pre)-(1-A)/(1-A_pre))] 
    
    # Define fluctuation model
    fluct.model.2 = lm(up.Q1.A ~ -1 + offset(Q2.A) + clever.H2, data=full.data)
    
    # Using estimated fluctuation parameter, update Q2.1
    full.data = full.data[, up.Q2.1 := predict(fluct.model.2, 
                                              newdata = data.frame(Q2.A=full.data$Q2.1,
                                                                   clever.H2=1/A_pre))] %>%
      .[, up.Q2.0 := predict(fluct.model.2,
                             newdata = data.frame(Q2.A=full.data$Q2.0,
                                                  clever.H2=-1/(1-A_pre)))] 
  }
  else {
    full.data = full.data[, up.Q1.A := A*up.Q1.1 + (1-A)*up.Q1.0 ]
    
    Q2mismod = gsub("Q1", "up.Q1.A", Q2mis)
    train.Q2 = lm(as.formula(Q2mismod), full.data)
    
    full.data = full.data[, Q2.1 := predict(train.Q2, newdata=data.frame(W=full.data$W,A=1))] %>%
      .[, Q2.0 := predict(train.Q2, newdata=data.frame(W=full.data$W,A=0))]
    
    # STEP 2 ------------------
    
    # Compute the observed values for up.Q1 and Q2
    # Define clever variable
    full.data = full.data[, Q2.A := A*Q2.1 + (1-A)*Q2.0 ] %>% 
      .[, clever.H2 := ((A/A_pre)-(1-A)/(1-A_pre))] 
    
    # Define fluctuation model
    fluct.model.2 = lm(up.Q1.A ~ -1 + offset(Q2.A) + clever.H2, data=full.data)
    
    # Using estimated fluctuation parameter, update Q2.1
    full.data = full.data[, up.Q2.1 := predict(fluct.model.2, 
                                              newdata = data.frame(Q2.A=full.data$Q2.1,
                                                                   clever.H2=1/A_pre))] %>%
      .[, up.Q2.0 := predict(fluct.model.2, 
                            newdata = data.frame(Q2.A=full.data$Q2.0,
                                                 clever.H2=-1/(1-A_pre)))]
  }
  
  # Compute the updated difference Q2.1 - Q2.0
  # Compute the observed value for up.Q2
  # Compute the value of the efficient influence function
  full.data = full.data[, delta.up.Q2 := up.Q2.1-up.Q2.0] %>%
    .[, up.Q2.A := A*up.Q2.1 + (1-A)*up.Q2.0] %>%
    .[, EIF := delta.up.Q2 + clever.H2*(up.Q1.A - up.Q2.A) + clever.H1*(Y - up.Q1.A)*R0]
  
  # Using the EIF, compute the asymptotic error
  asymp.sd = sd(full.data$EIF)/sqrt(N)
  
  # Save estimate and CI
  tmle.2step = c(mean(full.data$delta.up.Q2)-qnorm(0.975)*asymp.sd,
                 mean(full.data$delta.up.Q2)+qnorm(0.975)*asymp.sd,
                 mean(full.data$delta.up.Q2))
  
  return(tmle.2step)
}

################################################################################
# TMLE-1R
################################################################################
tmle1r_estimator = function(Q1pre='',Rpre=''){
  # Update R-predictions
  
  if(Rpre==''){
    # Model for missingness mechanism
    train.R = make_sl3_Task(
      data = full.data, outcome = 'R0', covariates = c('W','A'))
    
    R_fit = super_B$train(task = train.R)
    full.data = full.data[, R_pre := R_fit$predict(task = train.R)]
  }
  else{
    train.R = glm(as.formula(Rpre), full.data, family=binomial('logit'))
    
    full.data = full.data[, R_pre := predict(train.R, type='response')]
  }
  
  # Update Q-predictions
  if(Q1pre==''){
    # Model for Q1
    train.Q1 = make_sl3_Task(
      data = full.data[sel,], outcome = 'Y', covariates = c('W','A'))
    
    Q1_fit = super_C$train(task = train.Q1)
    
    # Prediction task for Q1
    pred.Q1 = make_sl3_Task(
      data = full.data, outcome = 'Y', covariates = c('W','A'))
    
    full.data = full.data[, Q1 := Q1_fit$predict(task = pred.Q1)]
  }
  else {
    # Model for Q1: MISSPECIFIED
    train.Q1 = lm(as.formula(Q1pre), data = full.data[sel,])
    
    full.data = full.data[, Q1 := predict(train.Q1, newdata = full.data)]
  }
  
  # Define clever variable
  full.data = full.data[, clever.H1 := (1/R_pre)*((A/A_pre)-(1-A)/(1-A_pre))] 
  
  # Define fluctuation model
  fluct.model.1 = lm(Y ~ -1 + offset(Q1) + clever.H1, data=full.data[sel,])
  
  # Auxiliary data frame for prediction
  temp.dt = copy(full.data)
  
  if(Q1pre==''){
    # Using estimated fluctuation parameter, update Q1.1
    temp.dt$A = 1
    pred.Q1 = make_sl3_Task(
      data = temp.dt, outcome = 'Y', covariates = c('W','A'))
    
    full.data = full.data[, Q1.1 := Q1_fit$predict(task = pred.Q1)] %>%
      .[, up.Q1.1 := predict(fluct.model.1, 
                            newdata = data.frame(Q1=full.data$Q1.1,
                                                 clever.H1=(1/R_pre)*(1/A_pre)))]
    
    # Using estimated fluctuation parameter, update Q1.0
    temp.dt$A = 0
    pred.Q1 = make_sl3_Task(
      data = temp.dt, outcome = 'Y', covariates = c('W','A'))
    
    full.data = full.data[, Q1.0 := Q1_fit$predict(task = pred.Q1)] %>%
      .[, up.Q1.0 := predict(fluct.model.1, 
                            newdata = data.frame(Q1=full.data$Q1.0,
                                                 clever.H1=(1/R_pre)*(-1/(1-A_pre))))]
  }
  else{
    # Using estimated fluctuation parameter, update Q1.1
    temp.dt$A = 1
    
    full.data = full.data[, Q1.1 := predict(train.Q1, newdata = temp.dt)] %>%
      .[, up.Q1.1 := predict(fluct.model.1, 
                           newdata = data.frame(Q1=full.data$Q1.1,
                                                clever.H1=(1/R_pre)*(1/A_pre)))]
    
    
    # Using estimated fluctuation parameter, update Q1.0
    temp.dt$A = 0
    
    full.data = full.data[, Q1.0 := predict(train.Q1, newdata = temp.dt)] %>%
      .[, up.Q1.0 := predict(fluct.model.1, 
                            newdata = data.frame(Q1=full.data$Q1.0,
                                                 clever.H1=(1/R_pre)*(-1/(1-A_pre))))]
  }
  
  # Compute the updated difference Q1.1 - Q1.0
  # Compute the observed value for up.Q1
  # Compute the value of the efficient influence function
  full.data = full.data[, delta.up.Q1 := up.Q1.1-up.Q1.0] %>%
    .[, up.Q1.A := A*up.Q1.1 + (1-A)*up.Q1.0] %>%
    .[, EIF := delta.up.Q1 + clever.H1*(Y - up.Q1.A)*R0]
  
  # Using the EIF, compute the asymptotic error
  asymp.sd = sd(full.data$EIF)/sqrt(N)
  
  # Save estimate and CI
  tmle.1step = c(mean(full.data$delta.up.Q1)-qnorm(0.975)*asymp.sd,
                 mean(full.data$delta.up.Q1)+qnorm(0.975)*asymp.sd,
                 mean(full.data$delta.up.Q1))
  
  return(tmle.1step)
}