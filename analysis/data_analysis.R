# ------------------------------
# Libraries
source("used_libraries.R")
source("estimators.R")

# ------------------------------
# Get the path to the directories
parent_dir = here::here()
sister_folder = file.path(parent_dir, "plots")
# Create the folder if it doesn't exist
if (!dir.exists(sister_folder)) {
  dir.create(sister_folder)
}

# ------------------------------
# Print graph
plot_1 = dagify(
  A ~ W,
  M ~ A,
  Z ~ A + M,
  Y ~ W + A + M + Z,
  R ~ M + Z
) %>% tidy_dagitty(layout = "kk") %>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_point(color='white',size=0.5) +
  geom_dag_edges() +
  geom_dag_text(color='black') +
  theme_dag()

# Save the plot
ggsave(filename = file.path(sister_folder, "plot_1_graph.png"),
    plot = plot_1, width = 6, height = 4, dpi = 300, bg = "white")

# ------------------------------
# Read generated data
full.data <- readRDS("gen_data.rds")
N = nrow(full.data)

###############################################################
#---------------------------------------------------------------
# Case 1: severe missingness/selection and no misspecification
#---------------------------------------------------------------
################################################################

# SELECTION MECHANISM
full.data$R0 = full.data$R1

#-------------------------------------------------------------------------------
# Estimator A: Oracle PATE -----------------------------------------------------

oracle.ttest = t.test(full.data$ITE)
pate = unname(c(oracle.ttest$conf.int[1],
                oracle.ttest$conf.int[2],
                oracle.ttest$estimate))

#-------------------------------------------------------------------------------
# Estimator B: Oracle SATE -----------------------------------------------------

sel = (full.data$R0==1)

oracle.ttest = t.test(full.data[sel,ITE])
sate = unname(c(oracle.ttest$conf.int[1],
                oracle.ttest$conf.int[2],
                oracle.ttest$estimate))

#-------------------------------------------------------------------------------
# Estimator 1: Unadjusted t-test -----------------------------------------------

unadj.mod = lm_robust(Y~A,full.data[sel,])

# Save estimate and CI
unadj.est = unname(c(unadj.mod$conf.low[2],
                     unadj.mod$conf.high[2],
                     unadj.mod$coefficients[2]))

#-------------------------------------------------------------------------------
# Estimator 2: TMLE CC ---------------------------------------------------------

tmle.ptest = tmlecc_estimator()

#-------------------------------------------------------------------------------
# Estimator 3: Doubly-inverse weighted estimator -------------------------------

diw.est = tmledw_estimator()

#-------------------------------------------------------------------------------
# Estimator 4: Nested regressions, T-learner -----------------------------------

nesreg.T = sr_estimator()

#-------------------------------------------------------------------------------
# Estimator 5: 2-step TMLE  ----------------------------------------------------

tmle.2step = tsr_estimator()

#-------------------------------------------------------------------------------
# Estimator 6: pre-exp TMLE  ---------------------------------------------------

tmle.1step = tmle1r_estimator()

#-------------------------------------------------------------------------------
# PUT ALL TOGETHER -------------------------------------------------------------

estimators = data.table(rbind(pate,sate,unadj.est,tmle.ptest,
                              tmle.1step,diw.est,nesreg.T,tmle.2step))
colnames(estimators) = c('lower','upper','point.est')
estimators$type = c('Oracle PATE','Oracle SATE',"Unadjusted",'TMLE CC',
                    'TMLE 1R','DW','SR','TSR')

estimators$type = factor(estimators$type,levels = rev(estimators$type))

# Bias in scale of percentage points over the outcome's standard deviation

usd = sd(full.data$Y)
bia = as.numeric(estimators[1,3])

estimators = estimators[,lower := (lower-bia)/usd] %>%
  .[,upper := (upper-bia)/usd] %>%
  .[,point.est := (point.est-bia)/usd] 

# Estimate and CI plot
plot_2 = ggplot(estimators, aes(y=type, x=point.est, group=type)) +
  geom_point(position=position_dodge(0.78)) +
  geom_errorbar(aes(xmin=lower, xmax=upper, color=type),
                width=0.5, position=position_dodge(0.78)) + 
  guides(color="none") + labs(x='Estimate',y='Estimator') + 
  theme_linedraw() + xlim(c(-0.2,0.32)) + geom_vline(xintercept = 0, linetype="dashed")

ggsave(filename = file.path(sister_folder, "plot_2_case_1.png"),
    plot = plot_2, width = 6, height = 4, dpi = 300, bg = "white")


###############################################################
#---------------------------------------------------------------
# Case 2: moderate missingness/selection and propensity score misspecification
#---------------------------------------------------------------
################################################################

# SELECTION MECHANISM
full.data$R0 = full.data$R2

# Misspecified propensity score
PSmis = 'A ~ I(W^-2)+I(W^2)'

#-------------------------------------------------------------------------------
# Estimator A: Oracle PATE -----------------------------------------------------

oracle.ttest = t.test(full.data$ITE)
pate = unname(c(oracle.ttest$conf.int[1],
                oracle.ttest$conf.int[2],
                oracle.ttest$estimate))

#-------------------------------------------------------------------------------
# Estimator B: Oracle SATE -----------------------------------------------------

sel = (full.data$R0==1)

oracle.ttest = t.test(full.data[sel,ITE])
sate = unname(c(oracle.ttest$conf.int[1],
                oracle.ttest$conf.int[2],
                oracle.ttest$estimate))

#-------------------------------------------------------------------------------
# Estimator 1: Unadjusted t-test -----------------------------------------------

unadj.mod = lm_robust(Y~A,full.data[sel,])

# Save estimate and CI
unadj.est = unname(c(unadj.mod$conf.low[2],
                     unadj.mod$conf.high[2],
                     unadj.mod$coefficients[2]))

#-------------------------------------------------------------------------------
# Estimator 2: TMLE CC ---------------------------------------------------------

tmle.ptest = tmlecc_estimator(PSmis = PSmis)

#-------------------------------------------------------------------------------
# Estimator 3: Doubly-inverse weighted estimator -------------------------------

diw.est = tmledw_estimator(PSmis = PSmis)

#-------------------------------------------------------------------------------
# Estimator 4: Nested regressions, T-learner -----------------------------------

nesreg.T = sr_estimator()

#-------------------------------------------------------------------------------
# Estimator 5: 2-step TMLE  ----------------------------------------------------

tmle.2step = tsr_estimator()

#-------------------------------------------------------------------------------
# Estimator 6: pre-exp TMLE  ---------------------------------------------------

tmle.1step = tmle1r_estimator()

#-------------------------------------------------------------------------------
# PUT ALL TOGETHER -------------------------------------------------------------

estimators = data.table(rbind(pate,sate,unadj.est,tmle.ptest,
                              tmle.1step,diw.est,nesreg.T,tmle.2step))
colnames(estimators) = c('lower','upper','point.est')
estimators$type = c('Oracle PATE','Oracle SATE',"Unadjusted",'TMLE CC',
                    'TMLE 1R','DW','SR','TSR')

estimators$type = factor(estimators$type,levels = rev(estimators$type))

# Bias in scale of percentage points over the outcome's standard deviation

usd = sd(full.data$Y)
bia = as.numeric(estimators[1,3])

estimators = estimators[,lower := (lower-bia)/usd] %>%
  .[,upper := (upper-bia)/usd] %>%
  .[,point.est := (point.est-bia)/usd] 

# Estimate and CI plot
plot_3 = ggplot(estimators, aes(y=type, x=point.est, group=type)) +
  geom_point(position=position_dodge(0.78)) +
  geom_errorbar(aes(xmin=lower, xmax=upper, color=type),
                width=0.5, position=position_dodge(0.78)) + 
  guides(color="none") + labs(x='Estimate',y='Estimator') + 
  theme_linedraw() + xlim(c(-0.2,0.32)) + geom_vline(xintercept = 0, linetype="dashed")

ggsave(filename = file.path(sister_folder, "plot_3_case_2.png"),
       plot = plot_3, width = 6, height = 4, dpi = 300, bg = "white")

###############################################################
#---------------------------------------------------------------
# Case 3: low missingness/selection and Q1 misspecification
#---------------------------------------------------------------
################################################################

# SELECTION MECHANISM
full.data$R0 = full.data$R3

# MISSPECIFIED MODELS
Q1pre = 'Y ~ I(W^(A+1))'
Q1pos = 'Y ~ I(W^(A+1))+M+I(log(Z))'

#-------------------------------------------------------------------------------
# Estimator A: Oracle PATE -----------------------------------------------------

oracle.ttest = t.test(full.data$ITE)
pate = unname(c(oracle.ttest$conf.int[1],
                oracle.ttest$conf.int[2],
                oracle.ttest$estimate))

#-------------------------------------------------------------------------------
# Estimator B: Oracle SATE -----------------------------------------------------

sel = (full.data$R0==1)

oracle.ttest = t.test(full.data[sel,ITE])
sate = unname(c(oracle.ttest$conf.int[1],
                oracle.ttest$conf.int[2],
                oracle.ttest$estimate))

#-------------------------------------------------------------------------------
# Estimator 1: Unadjusted t-test -----------------------------------------------

unadj.mod = lm_robust(Y~A,full.data[sel,])

# Save estimate and CI
unadj.est = unname(c(unadj.mod$conf.low[2],
                     unadj.mod$conf.high[2],
                     unadj.mod$coefficients[2]))

#-------------------------------------------------------------------------------
# Estimator 2: TMLE CC ---------------------------------------------------------

tmle.ptest = tmlecc_estimator(Q1pre = Q1pre)

#-------------------------------------------------------------------------------
# Estimator 3: Doubly-inverse weighted estimator -------------------------------

diw.est = tmledw_estimator()

#-------------------------------------------------------------------------------
# Estimator 4: Nested regressions, T-learner -----------------------------------

nesreg.T = sr_estimator(Q1pos = Q1pos)

#-------------------------------------------------------------------------------
# Estimator 5: 2-step TMLE  ----------------------------------------------------

tmle.2step = tsr_estimator(Q1pos = Q1pos)

#-------------------------------------------------------------------------------
# Estimator 6: pre-exp TMLE  ---------------------------------------------------

tmle.1step = tmle1r_estimator(Q1pre = Q1pre)

#-------------------------------------------------------------------------------
# PUT ALL TOGETHER -------------------------------------------------------------

estimators = data.table(rbind(pate,sate,unadj.est,tmle.ptest,
                              tmle.1step,diw.est,nesreg.T,tmle.2step))
colnames(estimators) = c('lower','upper','point.est')
estimators$type = c('Oracle PATE','Oracle SATE',"Unadjusted",'TMLE CC',
                    'TMLE 1R','DW','SR','TSR')

estimators$type = factor(estimators$type,levels = rev(estimators$type))

# Bias in scale of percentage points over the outcome's standard deviation

usd = sd(full.data$Y)
bia = as.numeric(estimators[1,3])

estimators = estimators[,lower := (lower-bia)/usd] %>%
  .[,upper := (upper-bia)/usd] %>%
  .[,point.est := (point.est-bia)/usd] 

# Estimate and CI plot
plot_4 = ggplot(estimators, aes(y=type, x=point.est, group=type)) +
  geom_point(position=position_dodge(0.78)) +
  geom_errorbar(aes(xmin=lower, xmax=upper, color=type),
                width=0.5, position=position_dodge(0.78)) + 
  guides(color="none") + labs(x='Estimate',y='Estimator') + 
  theme_linedraw() + xlim(c(-0.2,0.32)) + geom_vline(xintercept = 0, linetype="dashed")

ggsave(filename = file.path(sister_folder, "plot_4_case_3.png"),
       plot = plot_4, width = 6, height = 4, dpi = 300, bg = "white")

###############################################################
#---------------------------------------------------------------
# Case 4: low missingness/selection and  + selection probability misspecification
#---------------------------------------------------------------
################################################################

# SELECTION MECHANISM
full.data$R0 = full.data$R3

# MISSPECIFIED MODELS
Rpre = 'R0 ~ W+A'
Rpos = 'R0 ~ W+A+poly(M,2)+I(log(Z))'
Q2mis = 'Q1 ~ A+I(log(4+W))+A:I(W^2)'

#-------------------------------------------------------------------------------
# Estimator A: Oracle PATE -----------------------------------------------------

oracle.ttest = t.test(full.data$ITE)
pate = unname(c(oracle.ttest$conf.int[1],
                oracle.ttest$conf.int[2],
                oracle.ttest$estimate))

#-------------------------------------------------------------------------------
# Estimator B: Oracle SATE -----------------------------------------------------

sel = (full.data$R0==1)

oracle.ttest = t.test(full.data[sel,ITE])
sate = unname(c(oracle.ttest$conf.int[1],
                oracle.ttest$conf.int[2],
                oracle.ttest$estimate))

#-------------------------------------------------------------------------------
# Estimator 1: Unadjusted t-test -----------------------------------------------

unadj.mod = lm_robust(Y~A,full.data[sel,])

# Save estimate and CI
unadj.est = unname(c(unadj.mod$conf.low[2],
                     unadj.mod$conf.high[2],
                     unadj.mod$coefficients[2]))

#-------------------------------------------------------------------------------
# Estimator 2: TMLE CC ---------------------------------------------------------

tmle.ptest = tmlecc_estimator()

#-------------------------------------------------------------------------------
# Estimator 3: Doubly-inverse weighted estimator -------------------------------

diw.est = tmledw_estimator(Rpos = Rpos)

#-------------------------------------------------------------------------------
# Estimator 4: Nested regressions, T-learner -----------------------------------

nesreg.T = sr_estimator(Q2mis = Q2mis)

#-------------------------------------------------------------------------------
# Estimator 5: 2-step TMLE  ----------------------------------------------------

tmle.2step = tsr_estimator(Q2mis = Q2mis)

#-------------------------------------------------------------------------------
# Estimator 6: pre-exp TMLE  ---------------------------------------------------

tmle.1step = tmle1r_estimator(Rpre = Rpre)

#-------------------------------------------------------------------------------
# PUT ALL TOGETHER -------------------------------------------------------------

estimators = data.table(rbind(pate,sate,unadj.est,tmle.ptest,
                              tmle.1step,diw.est,nesreg.T,tmle.2step))
colnames(estimators) = c('lower','upper','point.est')
estimators$type = c('Oracle PATE','Oracle SATE',"Unadjusted",'TMLE CC',
                    'TMLE 1R','DW','SR','TSR')

estimators$type = factor(estimators$type,levels = rev(estimators$type))

# Bias in scale of percentage points over the outcome's standard deviation

usd = sd(full.data$Y)
bia = as.numeric(estimators[1,3])

estimators = estimators[,lower := (lower-bia)/usd] %>%
  .[,upper := (upper-bia)/usd] %>%
  .[,point.est := (point.est-bia)/usd] 

# Estimate and CI plot
plot_5 = ggplot(estimators, aes(y=type, x=point.est, group=type)) +
  geom_point(position=position_dodge(0.78)) +
  geom_errorbar(aes(xmin=lower, xmax=upper, color=type),
                width=0.5, position=position_dodge(0.78)) + 
  guides(color=FALSE) + labs(x='Estimate',y='Estimator') + 
  theme_linedraw() + xlim(c(-0.2,0.32)) + geom_vline(xintercept = 0, linetype="dashed")

ggsave(filename = file.path(sister_folder, "plot_5_case_4.png"),
       plot = plot_5, width = 6, height = 4, dpi = 300, bg = "white")
