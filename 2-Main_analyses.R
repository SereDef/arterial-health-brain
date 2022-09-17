# ==============================================================================
# =========================== 2. Main Analyses =================================
# ==============================================================================

# Required packages
invisible(lapply(c('mice','miceadds','openxlsx'), require, character.only = T));

# Load mids object
if (exists("imp_smri") == F) { 
  # Define path
  genrpath <- dirname(file.choose()) # project folder
  # Load imputed datasets
  mriset <- readRDS(file.path(genrpath,'results','imputation_list_smri.rds'))
  dtiset <- readRDS(file.path(genrpath,'results','imputation_list_dti.rds'))
  # Transform into mids
  imp_smri <- miceadds::datlist2mids(mriset)
  imp_dti <- miceadds::datlist2mids(dtiset)
  # clean up
  rm(mriset,dtiset)
}

# ------------------------------------------------------------------------------
# Subcortical volumes
subc_vol <- paste(c('accumbens','amygdala','caudate','hippocampus','pallidum','putamen','thalamus'),'z',sep='_')

# White matter tracts
tracts <- c('cgc','cst','unc','ilf','slf','fma','fmi')
trcts_FA <- paste(tracts,'FA','z', sep='_'); trcts_MD <- paste(tracts,'MD','z', sep='_')    

# ------------------------------------------------------------------------------
pool_mod <- function(outc, exp, sex='') {
  # some basic values
  sex_cov <- '+ sex'
  exp_name <- toupper(substr(outc,1,3))
  
  # Define the correct sample
  if (outc %in% c('tbv_13_z','gmv_13_z',subc_vol)) { 
    if (sex=='') { impset <- imp_smri } else { 
      impset <- miceadds::subset_datlist(imp_smri, subset = imp_smri$data$sex == sex,  toclass = 'mids') 
      sex_cov <- ''}
  } else if (any(sapply(c('FA','MD'), grepl, outc, ignore.case=T))) { 
    if (sex=='') { impset <- imp_dti } else { 
      impset <- miceadds::subset_datlist(imp_dti, subset = imp_dti$data$sex == sex,  toclass = 'mids') 
      sex_cov <- '' }
  }
  
  # Add total intracranial volume to covariates for subcortical outcomes
  if (outc %in% c(subc_vol,trcts_FA,trcts_MD)) { 
    tiv_cov <- '+ tiv_13'; exp_name <- toupper(sub("_z$", "", outc))
  } else { tiv_cov <- '' }
  
  # Define covariates
  covs1 <- paste(sex_cov,'+ age_13 + height_9 + ethnicity_dich',tiv_cov)
  covs2 <- '+ bmi_9_z + m_educ_6 + m_age'
  
  # Fit model and pool estimates
  pool_fit <- function(adj) {
    if (adj=='base') {
      fit <- with(impset, lm(as.formula(paste(outc,'~',exp,covs1))))
    } else {
      fit <- with(impset, lm(as.formula(paste(outc,'~',exp,covs1,covs2))))
    }
    p_fit <- mice::pool(fit) # pool results 
    mod <- summary(p_fit) # extract relevant information
    
    # Add confidence intervals
    mod$lci <- (mod$estimate - 1.96*mod$std.error)
    mod$uci <- (mod$estimate + 1.96*mod$std.error)
    # Add R squarred
    mod$rsq <- c(pool.r.squared(fit)[1], rep(NA, nrow(mod)-1)) # add a column for R2
    mod$rsq_adj <- c(pool.r.squared(fit, adjusted = T)[1], rep(NA, nrow(mod)-1)) # adjusted R2
    # Incorrect but just to check differences
    # mod$FDR_mod <- p.adjust(mod$p.value, method = 'fdr') # FDR correction per model
    
    # Round everything
    mod[,-1] <-round(mod[,-1], 3)
    
    # add model name as first column
    mod_name <- paste(exp_name,'-',toupper(substr(exp,1,3)), adj,'model')
    mod <- cbind(data.frame("model" = rep(mod_name, nrow(mod))), mod)
    # And one space in between models
    mod[nrow(mod)+1,] <- NA
    
    return(mod)
  }
  # Bid minimally and fully adjusted models
  mods <- rbind(pool_fit('base'), pool_fit('full'))
  
  return(mods)
}

# ------------------------------------------------------------------------------
# for each exposure loop through outcomes and pool models, also add proper FDR calculation
exp_mod <- function(exp, outcomes = c('tbv_13_z', 'gmv_13_z', 'mfa_13_z', 'mmd_13_z'), sex='') {
  mod <- data.frame()
  for (out in outcomes) { p <- pool_mod(out, exp, sex=sex) 
                          mod <- rbind(mod, p) }
  # Calculate FDR per predictor
  mod[,'FDR'] <- NA
  terms <- as.vector(unique(mod[,'term'])) # get predictors
  terms <- terms[!is.na(terms)] # clean NAs
  for (t in terms) {
    fdrs <- p.adjust(mod[!is.na(mod[,'term']) & mod[,'term']==t, 'p.value'], method = 'fdr')
    mod[!is.na(mod[,'term']) & mod[,'term']==t, 'FDR'] <- fdrs
  }
  mod[,'FDR'] <- round(mod[,'FDR'], 3)
  # add a column to highlight significant terms
  mod[,'sign_raw'] <- ifelse(mod[,'p.value'] < 0.05, '*', '') 
  mod[,'sign_fdr'] <- ifelse(mod[,'FDR'] < 0.05, '*', '')
  return(mod)
}

# ==============================================================================

for (e in c('imt','dis','sbp','dbp')) {
  # MAIN ANALYSES
  assign(e, exp_mod(paste0(e,'_9_z')))
  # SUBCORTICAL VOLUMES
  assign(paste0(e,'_subc_vol'), exp_mod(paste0(e,'_9_z'), subc_vol))
  # WHITE MATTER TRACTS
  assign(paste0(e,'_tract_fa'), exp_mod(paste0(e,'_9_z'), trcts_FA))
  assign(paste0(e,'_tract_md'), exp_mod(paste0(e,'_9_z'), trcts_MD))
  # SEX-STRATIFIED
  assign(paste0(e,'_f'), exp_mod(paste0(e,'_9_z'), sex='girl'))
  assign(paste0(e,'_m'), exp_mod(paste0(e,'_9_z'), sex='boy'))
}

# INTERACTION
for (e in c('sbp','dbp')) {
  assign(paste0(e,'_inter'), exp_mod(paste0('imt_9_z * ',e,'_9_z')))
}
inter <- rbind(sbp_inter, dbp_inter)

# ==============================================================================

main_modls <- list("IMT" = imt, "Dis" = dis, "SBP" = sbp, "DBP" = dbp)

openxlsx::write.xlsx(main_modls, file = file.path(genrpath,'results',
                                             paste0(Sys.Date(), "_Results.xlsx")), overwrite = T)

# ==============================================================================

supp_modls <- list("IMT_subc_vol" = imt_subc_vol, "Dis_subc_vol" = dis_subc_vol,
                   "SBP_subc_vol" = sbp_subc_vol, "DBP_subc_vol" = dbp_subc_vol,
                   "IMT_tract_FA" = imt_tract_fa, "Dis_tract_FA" = dis_tract_fa, 
                   "SBP_tract_FA" = sbp_tract_fa, "DBP_tract_FA" = dbp_tract_fa,
                   "IMT_tract_MD" = imt_tract_md, "Dis_tract_MD" = dis_tract_md, 
                   "SBP_tract_MD" = sbp_tract_md, "DBP_tract_MD" = dbp_tract_md,
                   "IMT_f" = imt_f, "Dis_f" = dis_f, "SBP_f" = sbp_f, "DBP_f" = dbp_f, 
                   "IMT_m" = imt_m, "Dis_m" = dis_m, "SBP_m" = sbp_m, "DBP_m" = dbp_m, 
                   "Inter" = inter)

openxlsx::write.xlsx(supp_modls, file = file.path(genrpath, 'results',
                                             paste0(Sys.Date(), "_SuppResults.xlsx")), overwrite = T)

# ==============================================================================
