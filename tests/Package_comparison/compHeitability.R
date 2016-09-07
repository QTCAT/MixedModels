library("coxme")
library("rrBLUP")
library("relMM")
require("heritability")

data(LDV)
data(K_atwell)

################################################################################
LDV1 <- LDV[as.character(LDV$replicate) == 1, ]
K1 <- K_atwell[as.character(LDV1$genotype), as.character(LDV1$genotype)]
hist(LDV1$LDV)

# total variance for comparison
var(LDV1$LDV)

# coxme (nlme like with "BFGS" optimizer)
lmmy1 <- lmekin(LDV ~ 1 +(1 | genotype), data = LDV1, 
                method = "REML", varlist = list(genotype = K1))
c(lmmy1$vcoef$genotype, lmmy1$sigma^2)

# rrBLUP (EMMA like with "Brent" optimizer)
lmmz1 <- mixed.solve(y = LDV1$LDV, K = K1, 
                    method = "REML")
c(lmmz1$Vu, lmmz1$Ve)

# relMM (lme4 with "Nelder-Mead" optimizer)
checks <- lmerControl(check.nobs.vs.nlev = "ignore", 
                      check.nobs.vs.nRE = "ignore")
lmmx1 <- relMM(LDV ~ 1 + (1 | genotype), data = LDV1, 
               REML = TRUE, covarrel = list(genotype = K1), 
               control = checks)
c(unlist(VarCorr(lmmx1)), sigma(lmmx1)^2)

# heritability
lmm1 <- marker_h2(data.vector = LDV1$LDV, geno.vector = LDV1$genotype, K = K1)
c(lmm1$va, lmm1$ve)

################################################################################
gt <- unique(as.character(LDV1$genotype))
K2 <- K_atwell[gt, gt]
hist(LDV$LDV)

# total variance for comparison
var(LDV$LDV) # or
var(residuals(lm(LDV ~ 1 + replicate, data = LDV)))

# coxme (nlme like with "BFGS" optimizer)
lmmy2 <- lmekin(LDV ~ 1 + replicate +(1 | genotype), data = LDV, 
                method = "REML", varlist = list(genotype = K2))
c(lmmy2$vcoef$genotype, lmmy2$sigma^2)

# rrBLUP (EMMA like with "Brent" optimizer)
Zrrblup <- model.matrix(~ 0 + genotype, LDV)
colnames(Zrrblup) <- gsub("genotype",  "", colnames(Zrrblup))
lmmz2 <- mixed.solve(y = LDV$LDV, 
                     X = model.matrix(~ 1 + replicate, LDV)[, -1],
                     Z = Zrrblup, 
                     K = K2[colnames(Zrrblup), colnames(Zrrblup)], 
                    method = "REML")
c(lmmz2$Vu, lmmz2$Ve)

# relMM (lme4 with "Nelder-Mead" optimizer)
checks <- lmerControl(check.nobs.vs.nlev = "ignore", 
                      check.nobs.vs.nRE = "ignore")
lmmx2 <- relMM(LDV ~ 1 + replicate + (1 | genotype), data = LDV, 
               REML = TRUE, covarrel = list(genotype = K2), 
               control = checks)
c(unlist(VarCorr(lmmx2)), sigma(lmmx2)^2)

# heritability
lmm2 <- marker_h2(data.vector = LDV$LDV, geno.vector = LDV$genotype, K = K2)
c(lmm2$va, lmm2$ve)

################################################################################

# total variance for comparison
var(LDV$LDV)
var(residuals(lm(LDV ~ 1 + replicate, data = LDV)))

# coxme (nlme like with "BFGS" optimizer)
Klmekin2 <- diag(nrow(K_atwell))
colnames(Klmekin2) <- rownames(Klmekin2) <- rownames(K_atwell)                
lmmy3 <- lmekin(LDV ~ 1 + replicate +(1 | genotype), data = LDV, 
                method = "REML", varlist = list(K2, Klmekin2))
c(lmmy3$vcoef$genotype, lmmy3$sigma^2)

# rrBLUP (EMMA like with "Brent" optimizer)
#   not aplicable

# relMM (lme4 with "Nelder-Mead" optimizer)
LDV$genotype2 <- LDV$genotype
checks <- lmerControl(check.nobs.vs.nlev = "ignore", 
                      check.nobs.vs.nRE = "ignore")
lmmx3 <- relMM(LDV ~ 1 + replicate + (1 | genotype) + (1 | genotype2), data = LDV, 
               REML = TRUE, covarrel = list(genotype = K2), 
               control = checks)
c(unlist(VarCorr(lmmx3)), sigma(lmmx3)^2)

# heritability
#   not aplicable

################################################################################



