# ==============================================================================
# =========================== 2. Main Analyses =================================
# ==============================================================================

# Required packages
invisible(lapply(c('mice','miceadds','splines','effects','openxlsx'), require, character.only = T));

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
  # Exclude one outlier
  imp_smri <- miceadds::subset_datlist(imp_smri, 
              subset = imp_smri$data$IDC != imp_smri$data$IDC[903],  toclass = 'mids') 
  # clean up
  rm(mriset,dtiset)
}

# ------------------------------------------------------------------------------
# Subcortical volumes
rois <- c('Accumbens', 'Amygdala', 'Caudate', 'Hippocampus', 'Pallidum', 'Putamen', 'Thalamus')
subc_vol <- paste(tolower(rois),'z',sep='_')
names(subc_vol) <- rois

# White matter tracts
tr_names <- c('Cingulate gyrus', 'Cortico-spinal tract', 'Uncinate fasciculus', 
              'Inferior longitudinal fasciculus', 'Superior longitudinal fasciculus', 
              'Major forceps', 'Minor forceps')
tracts <- c('cgc','cst','unc','ilf','slf','fma','fmi')
trcts_FA <- paste(tracts,'FA','z', sep='_'); trcts_MD <- paste(tracts,'MD','z', sep='_')    
names(trcts_FA) <- paste(tr_names,'(FA)'); names(trcts_MD) <- paste(tr_names,'(MD)'); 
# ------------------------------------------------------------------------------
pool_mod <- function(outc, exp, sex='') {
  # some basic values
  sex_cov <- '+ sex'
  exp_name <- toupper(substr(outc,1,3))
  
  # Define the correct sample
  if (any(sapply(c('tbv','gmv',subc_vol), grepl, outc))) { 
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
    tiv_cov <- '+ tiv_13'
  } else { tiv_cov <- '' }
  
  # Specify the names of the supplementary outcomes 
  if (outc %in% subc_vol) { exp_name <- names(subc_vol)[which(subc_vol==outc)] 
  } else if (outc %in% trcts_FA) { exp_name <- names(trcts_FA)[which(trcts_FA==outc)] 
  } else if (outc %in% trcts_MD) { exp_name <- names(trcts_MD)[which(trcts_MD==outc)]  }
  
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
    # cat(t,'\t', length(fdrs),'\n')
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
  # non standardized (original scale)
  assign(paste0(e,'_unscaled'), 
         exp_mod(paste0(e,'_9'), 
                 paste0(c('tbv', 'gmv', 'mfa', 'mmd'),'_13')))
}

# INTERACTION
for (e in c('sbp','dbp')) {
  assign(paste0(e,'_inter'), exp_mod(paste0('imt_9_z * ',e,'_9_z')))
}
inter <- rbind(sbp_inter, dbp_inter)

# ==============================================================================

main_modls <- list("IMT" = imt, "Dis" = dis, "SBP" = sbp, "DBP" = dbp, 
                   "IMT_orig" = imt_unscaled, "Dis_orig" = dis_unscaled, 
                   "SBP_orig" = sbp_unscaled, "DBP_orig" = dbp_unscaled)

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
# non linear terms 
imt_mmd_e3 <- exp_mod('imt_9_z + I(imt_9_z^3)', c('mmd_13_z'))
imt_mmd_e5 <- exp_mod('imt_9_z + I(imt_9_z^5)', c('mmd_13_z'))
# imt_mmd <- exp_mod('imt_9_z + I(imt_9_z^3) + I(imt_9_z^5)', c('mmd_13_z'))
imt_nlin <- rbind(imt_mmd_e3,imt_mmd_e5)

sbp_mfa_p3 <- exp_mod('poly(sbp_9_z, 3, raw=T)', c('mfa_13_z'))
sbp_mfa_s3 <- exp_mod('ns(sbp_9_z, 3)', c('mfa_13_z'))
sbp_mfa_s4 <- exp_mod('ns(sbp_9_z, 4)', c('mfa_13_z'))
sbp_mfa_s5 <- exp_mod('ns(sbp_9_z, 5)', c('mfa_13_z'))

sbp_mmd_e2 <- exp_mod('sbp_9_z + I(sbp_9_z^2)', c('mmd_13_z'))
sbp_mmd_e4 <- exp_mod('sbp_9_z + I(sbp_9_z^4)', c('mmd_13_z'))
sbp_mmd_s2 <- exp_mod('ns(sbp_9_z, 2)', c('mmd_13_z'))
sbp_mmd_s4 <- exp_mod('ns(sbp_9_z, 4)', c('mmd_13_z'))
sbp_mmd_s5 <- exp_mod('ns(sbp_9_z, 5)', c('mmd_13_z'))

sbp_nlin <- rbind(sbp_mfa_p3,sbp_mfa_s3,sbp_mfa_s4,sbp_mfa_s5,
                  sbp_mmd_e2,sbp_mmd_e4,sbp_mmd_s2,sbp_mmd_s4,sbp_mmd_s5)

dbp_tbv_e2 <- exp_mod('dbp_9_z + I(dbp_9_z^2)', c('tbv_13_z'))
dbp_tbv_s2 <- exp_mod('ns(dbp_9_z, 2)', c('tbv_13_z'))

dbp_gmv_e2 <- exp_mod('dbp_9_z + I(dbp_9_z^2)', c('gmv_13_z'))
dbp_gmv_s2 <- exp_mod('ns(dbp_9_z, 2)', c('gmv_13_z'))

dbp_mmd_e4 <- exp_mod('dbp_9_z + I(dbp_9_z^4)', c('mmd_13_z'))
dbp_mmd_p4 <- exp_mod('poly(dbp_9_z, 4, raw=T)', c('mmd_13_z'))

dbp_nlin <- rbind(dbp_tbv_e2,dbp_tbv_s2,dbp_gmv_e2,dbp_gmv_s2,dbp_mmd_e4,dbp_mmd_p4)

nlin_modls <- list(imt_nlin, sbp_nlin, dbp_nlin)

openxlsx::write.xlsx(nlin_modls, file = file.path(genrpath, 'results',
                                                  paste0(Sys.Date(), "_NonlinResults.xlsx")), overwrite = T)

# PLOT SPLINES -----------------------------------------------------------------
plot_nln <- function(outc, exp, m=15) {
  covs <- '+ sex + age_13 + height_9 + ethnicity_dich + bmi_9_z + m_educ_6 + m_age'
  fit <- lm(as.formula(paste0(outc,' ~ ns(',exp,', 5)',covs)), data=complete(imp_dti,m))
  eff <- effects::Effect(exp, fit) # transformation=list(inverse=exp))
  plot(eff, main=paste(toupper(substr(outc,1,3)),'~',toupper(substr(exp,1,3)),'non-linear relationship'), 
       xlab=paste(toupper(substr(exp,1,3)),'(z-score)'), ylab=paste(toupper(substr(outc,1,3)),'(z-score)'))
}

# for (e in paste0(c('imt','dis','sbp','dbp'),'_9_z')) {
#   for (o in paste0(c('tbv', 'gmv', 'mfa', 'mmd'),'_13_z')) {
#     cat('plot_nln("',o,'", "',e,'")\n', sep='')
#   }
# }

pdf(file.path(genrpath, 'results','nln_imp_fit.pdf'))
plot_nln("tbv_13_z", "imt_9_z")
plot_nln("gmv_13_z", "imt_9_z")
plot_nln("mfa_13_z", "imt_9_z")
plot_nln("mmd_13_z", "imt_9_z") # *
plot_nln("tbv_13_z", "dis_9_z")
plot_nln("gmv_13_z", "dis_9_z")
plot_nln("mfa_13_z", "dis_9_z")
plot_nln("mmd_13_z", "dis_9_z")
plot_nln("tbv_13_z", "sbp_9_z")
plot_nln("gmv_13_z", "sbp_9_z")
plot_nln("mfa_13_z", "sbp_9_z") # *
plot_nln("mmd_13_z", "sbp_9_z") # *
plot_nln("tbv_13_z", "dbp_9_z") # *
plot_nln("gmv_13_z", "dbp_9_z") # *
plot_nln("mfa_13_z", "dbp_9_z")
plot_nln("mmd_13_z", "dbp_9_z") # *
dev.off()

# ==============================================================================