# ------------------------------
# Libraries
source("used_libraries.R")

# ------------------------------
# Causal mechanisms
# Causal mechanism for treatment assignment
#coef.A = c(0.00, 0.90, -0.09)
fun.A = function(w, u){
  dat = as.numeric(c(1, w, sign(w)*w^2))
  lo = dot(coef.A, dat) + u
  return(as.numeric(lo>0))
}

# Causal mechanism for mediator 1
#coef.M = c(-0.50, 1.00)
fun.M = function(a, u){
  dat = as.numeric(c(1, a))
  li = dot(coef.M, dat) + u
  return(li)
}

# Causal mechanism for mediator 2
#coef.Z = c(4.20, 0.25, 0.30, 0.05)
fun.Z = function(a, m, u){
  dat = as.numeric(c(1, 2*a-1, m, (2*a-1)*m))
  li = 0.12*(dot(coef.Z, dat) + u)^2
  return(li)
}

# Causal mechanism for outcome
#coef.Y = c(0.00, 1.80, 0.20, 0.75, 0.50, 2.00, 0.50, 0.80)
fun.Y = function(w, a, m, z, u){
  dat = as.numeric(c(1, w, w^3, 2*a-1, 
                    (2*a-1)*w, m, (2*a-1)*m, (2*a-1)*z))
  li = dot(coef.Y, dat) + u
  return(li)
}

# Causal mechanism for missingness mechanisms. R
#coef.R = c(0.29, 0.54) 
#r.bias = 0
fun.R = function(m, z, u){
  dat = as.numeric(c(m, z))
  lo = r.bias + dot(coef.R, dat) + u
  return(as.numeric(lo>0))
}

# ------------------------------
# Generate data function
generate_data = function(N = 1e4,        # sample size
                         used_seed = 77, # seed
                         coefA = c(0.00, 0.90, -0.09),
                         coefM = c(-0.50, 1.00),
                         coefZ = c(4.20, 0.25, 0.30, 0.05),
                         coefY = c(0.00, 1.80, 0.20, 0.75, 0.50, 2.00, 0.50, 0.80),
                         coefR = c(0.29, 0.54),
                         biasR = c(-1.20, -0.35, 0.40)){
    # Seed 
    set.seed(used_seed)

    # Set coefficients
    coef.A <<- coefA
    coef.M <<- coefM
    coef.Z <<- coefZ
    coef.Y <<- coefY
    coef.R <<- coefR

    # Noises and exogenous variables
    # Generate exogenous variables: independent confounders and noises, 
    # plus fixed treatment assignments A.1 and A.0
    full.data = data.table(noise.A = rnorm(N, 0, 1),
                            noise.M = rnorm(N, 0, 1),
                            noise.Z = rnorm(N, 0, 1),
                            noise.Y = rnorm(N, 0, 11),
                            noise.R1 = rnorm(N, 0, 1),
                            noise.R2 = rnorm(N, 0, 1),
                            noise.R3 = rnorm(N, 0, 1),
                            W = rnorm(N, 0, 1),
                            A.1 = 1, A.0 = 0)

    # Generate observations
    full.data = full.data[, A:=mapply(fun.A, W, noise.A)] %>%
    .[, M:=mapply(fun.M, A, noise.M)] %>%
    .[, M.1:=mapply(fun.M, A.1, noise.M)] %>%
    .[, M.0:=mapply(fun.M, A.0, noise.M)] %>%
    .[, Z:=mapply(fun.Z, A, M, noise.Z)] %>%
    .[, Z.1:=mapply(fun.Z, A.1, M.1, noise.Z)] %>%
    .[, Z.0:=mapply(fun.Z, A.0, M.0, noise.Z)] %>%
    .[, Y:=mapply(fun.Y, W, A, M, Z, noise.Y)] %>%
    .[, Y.1:=mapply(fun.Y, W, A.1, M.1, Z.1, noise.Y)] %>%
    .[, Y.0:=mapply(fun.Y, W, A.0, M.0, Z.0, noise.Y)] %>%
    .[, ITE:=Y.1-Y.0] 

    # Generate missingness indicator, case 1: R1
    r.bias <<- biasR[1]
    full.data = full.data[, R1:=mapply(fun.R, M, Z, noise.R1)]

    # Generate missingness indicator, case 2: R2
    r.bias <<- biasR[2]
    full.data = full.data[, R2:=mapply(fun.R, M, Z, noise.R2)]

    # Generate missingness indicator, case 2: R3 y R4
    r.bias <<- biasR[3]
    full.data = full.data[, R3:=mapply(fun.R, M, Z, noise.R3)]

    # Return
    return(full.data)
}

# ------------------------------
# Save data
gen_data = generate_data()
saveRDS(gen_data, file = "gen_data.rds")