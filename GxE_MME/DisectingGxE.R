# R codes associated with Tutorial
 # Load the library 
  library(lme4)

  rm(list=ls())  # Remove previous work
# Load data
  demo.data <- read.csv("./Data/demo.data1.csv", header = TRUE)
# Factor conversion
  demo.data$Environment <- as.factor(demo.data$Environment)
  demo.data$Rep <- as.factor(demo.data$Rep)
  demo.data$Genotype <- as.factor(demo.data$Genotype)
# Interaction term for GxE
  demo.data$gen_env <- interaction(demo.data$Genotype, demo.data$Environment)
# str(demo.data) makes sure everything’s in order—no stray variables 
  str(demo.data)

kable(head(demo.data))

# Incidence/design matrix for fixed effects: Rep
  X <- model.matrix(~ Rep, demo.data)
# Random effect design matrix for Genotype
  Z_g  <- model.matrix(~ -1 + Genotype, demo.data)
# Random effect design matrix for Environment   
  Z_e  <- model.matrix(~ -1 + Environment, demo.data)   
# Random effect design matrice for G x E
  Z_ge <- model.matrix(~ -1 + gen_env, demo.data)


# Enter your estimated variance components, usually from lme4::VarCorr(mod_gxe_env)
# Genotype variance
 sigmau   <- 114010 
# Residual variance
  sigmae   <- 331603  
# Environment variance 
  sigmaenv <- 16218   
# GxE variance
  sigmage  <- 229771  
# Build Diagonal Matrix for each:assumption that all random effects are independent and have equal variance.
  I_g <- diag(ncol(Z_g))
  I_e <- diag(ncol(Z_e))
  I_ge <- diag(ncol(Z_ge))
# Calculate lambda (variance ratios): Shrinkage factors.
  lambda_g  <- sigmae / sigmau # For genotypes
  lambda_e  <- sigmae / sigmaenv # For environments
  lambda_ge <- sigmae / sigmage # For g x e Interactions

  XtX <- t(X) %*% X

  XtZg   <- t(X) %*% Z_g   # (X'Z_g) of MME
  XtZe   <- t(X) %*% Z_e   # (X'Z_e) of MME
  XtZge  <- t(X) %*% Z_ge  # (X'Z_{ge}) of MME

  ZgtX  <- t(Z_g) %*% X    # (Z_g'X) of MME
  ZetX  <- t(Z_e) %*% X    # (Z_e'X) of MME
  ZgetX <- t(Z_ge) %*% X   # (Z_{ge}'X) of MME

  ZgtZg   <- t(Z_g) %*% Z_g    # (Z_g'Z_g) of MME
  ZgtZe   <- t(Z_g) %*% Z_e    # (Z_g'Z_e) of MME
  ZgtZge  <- t(Z_g) %*% Z_ge   # (Z_g'Z_{ge}) of MME
  
  ZetZg   <- t(Z_e) %*% Z_g    # (Z_e'Z_g) of MME
  ZetZe   <- t(Z_e) %*% Z_e    # (Z_e'Z_e) of MME
  ZetZge  <- t(Z_e) %*% Z_ge   # (Z_e'Z_{ge}) of MME
  
  ZgetZg   <- t(Z_ge) %*% Z_g   # (Z_{ge}'Z_g) of MME
  ZgetZe   <- t(Z_ge) %*% Z_e   # (Z_{ge}'Z_e) of MME
  ZgetZge  <- t(Z_ge) %*% Z_ge  # (Z_{ge}'Z_{ge}) of MME

  LHS_top   <- cbind(XtX, XtZg, XtZe, XtZge)

  LHS_mid1  <- cbind(ZgtX, ZgtZg + lambda_g * I_g, ZgtZe, ZgtZge)

  LHS_mid2  <- cbind(ZetX, ZetZg, ZetZe + lambda_e * I_e, ZetZge)

  LHS_bot   <- cbind(ZgetX, ZgetZg, ZgetZe, ZgetZge + lambda_ge * I_ge)

  LHS <- rbind(LHS_top, LHS_mid1, LHS_mid2, LHS_bot)

# Replication's view of Yield
  Xty    <- t(X) %*% demo.data$Yield 
# Genotype’s take
  Zgty   <- t(Z_g) %*% demo.data$Yield   
  # Environment’s feedback
  Zety   <- t(Z_e) %*% demo.data$Yield   
# G × E interaction's drama
  Zgety  <- t(Z_ge) %*% demo.data$Yield  
# RHS
  RHS <- rbind(Xty, Zgty, Zety, Zgety) 

sol <- solve(LHS, RHS)

# Saving the sol in Outputs as dataframe
  Outputs<-data.frame(sol)
  kable(head(Outputs))

# Number of fixed effects (Replication)
    n_fixed <- ncol(X) 
# Number of genotypes is 200
    n_gen   <- ncol(Z_g)  
 # Number of environments is 2
   n_env <- ncol(Z_e)           
# Number of genotype × environment 400
  n_ge    <- ncol(Z_ge)  
# BLUEs for Replication effects (Second)
  BLUEs    <- sol[1:n_fixed, ]
# BLUPs for genotype effects (starts after Blues of replications)
  BLUPs_gen <- sol[(n_fixed+1):(n_fixed+n_gen), ]
# BLUPs for environment (starts after replication and genotypes effects)
  BLUPs_env <- sol[(n_fixed+n_gen+1):(n_fixed+n_gen+n_env), ]
# BLUPs for g x e interaction effects (starts after fixed effect plus genotype plus environment)
  BLUPs_ge  <- sol[(n_fixed+n_gen+n_env+1):(n_fixed+n_gen+n_env+n_ge), ]
# Add intercept for comparison (assume first BLUE is intercept)

  intercept <- BLUEs[1]
  BLUPs_G_MAN <- data.frame(BLUPsG_MAN=BLUPs_gen + intercept) # BLUps for Genotypes
  kable(head(BLUPs_G_MAN))
  BLUPs_Env_MAN <- data.frame(BLUPsE_MAN=BLUPs_env + intercept) # BLUps for Environment
  kable(head( BLUPs_Env_MAN))
  BLUPs_GE_MAN  <- data.frame(BLUPsGE_MAN=BLUPs_ge + intercept)# BLUps for G X e
  kable(head(BLUPs_GE_MAN ))

BLUPs_GE_MAN$Genotype <- sub("\\.Env[12]$", "", rownames(BLUPs_GE_MAN))
BLUPs_GE_MAN$Environment <- sub(".*\\.(Env[12])$", "\\1", rownames(BLUPs_GE_MAN))
BLUPs_GE_MAN_ordered <- BLUPs_GE_MAN[order(BLUPs_GE_MAN$Genotype, BLUPs_GE_MAN$Environment), ]
BLUPs_GE_MAN_ordered$Genotype<- sub("^gen_env", "", BLUPs_GE_MAN_ordered$Genotype)
BLUPs_GE_MAN_ordered <- BLUPs_GE_MAN_ordered[order(as.numeric(BLUPs_GE_MAN_ordered$Genotype)), ]

# G x E Model in lme4
  mod_gxe <- lmer(Yield ~ Rep + (1 | Genotype) + (1 | Environment) + (1 | Genotype:Environment), 
                  data = demo.data)
# Get Summary of model
  summary(mod_gxe)

# For Replications
  fixef_vals <- fixef(mod_gxe)
  fixef_vals

# Get BLUPs for Genotypes and add intercept
  BLUP_G_lme4 <- data.frame(BLUPsG_LME4=ranef(mod_gxe)$Genotype[,1] + fixef_vals[1])
# Add names similar to manually generated file
 row.names(BLUP_G_lme4)<-paste0("Genotype", row.names(BLUP_G_lme4))
# Get BLUPs for Environment and add intercept
  BLUP_E_lme4 <- data.frame(BLUPsE_LME4=ranef(mod_gxe)$Environment[,1] + fixef_vals[1])
  # Add names similar to manually generated file
  row.names(BLUP_E_lme4)<-paste0("EnvironmentEnv", row.names(BLUP_E_lme4))
# Get BLUPs for G x E  and add intercept
  BLUP_GE_lme4 <- data.frame(BLUPsGE_LME4=ranef(mod_gxe)$`Genotype:Environment`[,1] + fixef_vals[1])
  

# Extract the fixed effects (includes the intercept)
fixef_vals <- fixef(mod_gxe)
# BLUPs for genotypes
blup_gen <- ranef(mod_gxe)$gen[,1] + fixef_vals[1]
# BLUPs for environments
blup_env <- ranef(mod_gxe)$env[,1] + fixef_vals[1]
# BLUPs for G × E interaction
blup_gxe <- ranef(mod_gxe)$`gen:env`[,1] + fixef_vals[1]

# BLUPs of Genotypes
    table(round(BLUPs_G_MAN$BLUPsG_MAN,2)==round(BLUP_G_lme4$BLUPsG_LME4,2))
# BLUPs of G x E
    table(round(BLUPs_GE_MAN_ordered$BLUPsGE_MAN, 1)==round(BLUP_GE_lme4$BLUPsGE_LME4, 1))

BLUPs_G_COMB <- cbind(BLUPs_G_MAN, BLUP_G_lme4)
kable(head(BLUPs_G_COMB))

row.names(BLUP_GE_lme4)<-row.names(BLUPs_GE_MAN_ordered)
BLUPs_GE_COMB <- cbind(BLUPs_GE_MAN_ordered, BLUP_GE_lme4)
kable(head(BLUPs_GE_COMB))
