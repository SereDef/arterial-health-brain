# ==============================================================================
# =========================== 2. Main Analyses =================================
# ==============================================================================

# Required packages
invisible(lapply(c('mice','miceadds','splines','effects','openxlsx'), require, character.only = T));
invisible(lapply(c('broom.mixed','lme4','lmerTest'), require, character.only = T)); # 'jtools','reghelper','nlme'

datedat <- '240523'
daterun <- '190923' # datedat # format(Sys.Date(), '%d%m%y')

# Load mids object
if (exists('genrpath') == F) { 
  # Define paths
  genrpath <- dirname(file.choose()) # <== choose project folder
  datpath  <- file.path(genrpath, paste0('results_',datedat)) 
  # results folder
  if (daterun != datedat) {
    respath  <- file.path(genrpath, paste0('results_',daterun))
    dir.create(respath) } else {  respath <- datpath }
}

# Exclude one outliers
# impset <- miceadds::subset_datlist(impset, 
#             subset = impset$data$IDC != impset$data$IDC[903],  toclass = 'mids')

# ------------------------------------------------------------------------------
# Sub-cortical volumes
rois <- c('Accumbens', 'Amygdala', 'Caudate', 'Hippocampus', 'Pallidum', 'Putamen', 'Thalamus')
subc_vol <- paste(tolower(rois),'13_z',sep='_')
names(subc_vol) <- rois

# White matter tracts
tr_names <- c('Cingulate gyrus', 'Cortico-spinal tract', 'Uncinate fasciculus', 
              'Inferior longitudinal fasciculus', 'Superior longitudinal fasciculus', 
              'Major forceps', 'Minor forceps')
tracts <- c('cgc','cst','unc','ilf','slf','fma','fmi')
trcts_FA <- paste(tracts,'FA','13','z', sep='_'); names(trcts_FA) <- paste(tr_names,'(FA)')
trcts_MD <- paste(tracts,'MD','13','z', sep='_'); names(trcts_MD) <- paste(tr_names,'(MD)')

# ==============================================================================
# INSTRUCTIONS =================================================================

samples <- c('base9') # default = c('base9', 'full', 'base5', 'smri', 'dti')
outcs <- c('other') # ,'unscaled','regional','main_10y','other') # default = c('main','unscaled','regional','main_10y','other') 
expos <- c('imt','dis','sbp','dbp') # default = c('imt','dis','sbp','dbp')

NAME_ANALYSIS <- 'tiv_sex' # default = ''
supp <- '' # ' + gest_weight + gest_age_birth' # default = ''

sex_interaction <- TRUE # default = T
BPIMT_interaction <- FALSE # default = T
longit <- FALSE # default = T
run_fulladj = TRUE # default = T

# ------------------------------------------------------------------------------
# for each exposure loop through outcomes and pool models, also add proper FDR calculation
exp_mod <- function(exp, outcomes = c('tbv_13_z', 'gmv_13_z', 'mfa_13_z', 'mmd_13_z'), impset=imput, sex='', add_covs='') {
  
  pool_mod <- function(outc, exp, impset, sexsub, add_covs) {
    # Initialize
    out_name <- toupper(substr(outc,1,3)) # first three letters
    sex_cov <- '+ sex'; tiv_cov <- ''
    
    # Stratify by sex
    if (sexsub!='') { impset <- mice::filter(impset, sex==sexsub); sex_cov <- ''
      # Print out sample size
      cat('Sample: ', nrow(impset$data),'\n')
    }
    
    # Add total intracranial volume to covariates for subcortical outcomes
    if (outc %in% c(subc_vol, trcts_FA, trcts_MD)) { tiv_cov <- '+ tiv_13' }
    
    # Specify the names of the supplementary outcomes 
    if (outc %in% subc_vol) {        out_name <- names(subc_vol)[which(subc_vol==outc)] 
    } else if (outc %in% trcts_FA) { out_name <- names(trcts_FA)[which(trcts_FA==outc)] 
    } else if (outc %in% trcts_MD) { out_name <- names(trcts_MD)[which(trcts_MD==outc)]  }
    
    # Define covariates
    covs1 <- paste(sex_cov,'+ age_mri_13 + age_gap_13 + height_10',tiv_cov)
    covs2 <- paste0('+ ethn_dich + bmi_10_z + tot_educ_6 + m_age', add_covs)
    
    # Fit model and pool estimates
    pool_fit <- function(adj) {
      if (adj=='base') { covs2 <- '' } # reduce covariates to base only 
      # Fit the model and pool estimates
      fit <- with(impset, lm(as.formula(paste(outc,'~',exp, covs1, covs2))))
      p_fit <- mice::pool(fit)
      mod <- summary(p_fit, 'all', conf.int = 0.95) # extract relevant information
      
      # Add confidence intervals
      names(mod)[which(names(mod) %in% c('2.5 %','97.5 %'))] <- c('lci','uci')
      
      # Add R squarred
      mod$rsq <- c(pool.r.squared(fit)[1], rep(NA, nrow(mod)-1)) # add a column for R2
      mod$rsq_adj <- c(pool.r.squared(fit, adjusted = T)[1], rep(NA, nrow(mod)-1)) # adjusted R2
      
      # Round everything
      mod[,-1] <-round(mod[,-1], 3)
      
      # add model name as first column
      mod_name <- paste(out_name,'-',toupper(substr(exp,1,3)), adj,'model')
      mod <- cbind(data.frame('model' = rep(mod_name, nrow(mod))), mod)
      
      return(mod)
    }
    # Bid minimally and fully adjusted models
    if (run_fulladj) { mods <- rbind(pool_fit('base'), pool_fit('full')) } else { mods <- pool_fit('base') }
    
    # And one space in between models
    mods[nrow(mods)+1,] <- NA
    
    return(mods)
  }
  # ----------------------------------------------------------------------------
  mod <- data.frame()
  for (out in outcomes) { p <- pool_mod(out, exp, impset=impset, sexsub=sex, add_covs=add_covs) 
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
  
  col_order <- c('model','term','estimate','std.error','statistic','df','p.value','FDR',
                 'lci','uci','rsq','rsq_adj','sign_raw','sign_fdr',
                 'm','riv','lambda','fmi','ubar','b','t','dfcom')
  mod <- mod[, col_order]
  
  return(mod)
}

# Let's do this ================================================================
for (sample in samples) {
  # Load imputed dataset
  imput <- readRDS(file.path(datpath, paste0('imp_',sample,'_',datedat,'.rds')))
  
  # RM LATER !!! ===============================================================
  imput2 =  readRDS(paste0('/Users/Serena/Desktop/arterial brain/results_190923/imp_',sample,'_190923.rds'))
  il2 = complete(imput2, 'long',include=T)[,c('.imp','.id','IDC','p_educ_cont','ventrix_10','ventrix_13')]
  
  l = merge(complete(imput, 'long',include=T), il2, by=c('.imp','.id','IDC'))
  
  l$tot_educ_6 = l$m_educ_cont + l$p_educ_cont
  l$age_gap_13 = l$age_mri_13 - l$age_10
  
  imput = as.mids(l)
  
  # Standardize
  mains <- paste0(c('tiv','cortex','subcort','wmv','csf','ventrix'),'_13') # c('gest_weight', 'gest_age_birth')
  sampz <- miceadds::scale_datlist(miceadds::mids2datlist(imput), orig_var = mains, trafo_var = paste0(mains, '_z'))
  imput <- miceadds::datlist2mids(sampz)
  
  # Print out sample size
  cat('\nSample: ',sample,' N =', nrow(imput$data),'\n')
  
  for (e in expos) { cat(e,'... ')

    ze <- paste0(e,'_10_z')
     
    # MAIN ANALYSES
    if ('main' %in% outcs) { assign(e, exp_mod(ze, add_covs=supp)) }
    
    # non standardized (original scale)
    if ('unscaled' %in% outcs) { assign(paste0(e,'_unscaled'), 
                                        exp_mod(paste0(e,'_10'),
                                                paste0(c('tbv', 'gmv', 'mfa', 'mmd'),'_13'), 
                                                add_covs=supp)) }
    
    if (sex_interaction) {
      # sex interaction
      assign(paste0(e,'_sex_inter'), exp_mod(paste0('sex * ',e,'_10_z'), add_covs=supp))
      # SEX-STRATIFIED
      assign(paste0(e,'_f'), exp_mod(ze, sex='girl', add_covs=supp))
      assign(paste0(e,'_m'), exp_mod(ze, sex='boy',  add_covs=supp))
    }
    
    # Regional analyses
    if ('regional' %in% outcs) {
      # SUBCORTICAL VOLUMES
      assign(paste0(e,'_subc_vol'), exp_mod(ze, subc_vol, add_covs=supp))
      # WHITE MATTER TRACTS
      assign(paste0(e,'_tract_fa'), exp_mod(ze, trcts_FA, add_covs=supp))
      assign(paste0(e,'_tract_md'), exp_mod(ze, trcts_MD, add_covs=supp))
    }
    
    # Cross-sectional brain (MRI at 10 years)
    if ('main_10y' %in% outcs) { assign(paste0(e,'_brain10'), 
                                        exp_mod(ze, paste0(c('tbv','gmv','mfa','mmd'),'_10_z'), 
                                                add_covs=supp)) }
    # Other MRI outcomes
    if ('other' %in% outcs) { assign(paste0(e,'_other'),
                                     exp_mod(ze, paste0(c('tiv','cortex','subcort','wmv','csf','ventrix'),'_13_z'), 
                                             add_covs=supp)) }
  }
  
  # IMT * BP INTERACTION
  if (BPIMT_interaction) {
    for (e in c('sbp','dbp')) { assign(paste0(e,'_inter'), exp_mod(paste0('imt_10_z * ',e,'_10_z'), add_covs=supp)) }
    inter <- rbind(sbp_inter, dbp_inter)
  }
  
  # SAVE =======================================================================
  modls = list()
  if ('main' %in% outcs) { modls <- append(modls,
             list('IMT' = imt, 'Dis' = dis, 'SBP' = sbp, 'DBP' = dbp)) }
 
  if ('unscaled' %in% outcs) { modls <- append(modls,
             list('IMT_orig' = imt_unscaled, 'Dis_orig' = dis_unscaled,
                  'SBP_orig' = sbp_unscaled, 'DBP_orig' = dbp_unscaled)) }
  
  if (sex_interaction) { modls <- append(modls,
             list('IMT_sexint' = imt_sex_inter, 'Dis_sexint' = dis_sex_inter, 'SBP_sexint' = sbp_sex_inter, 'DBP_sexint' = dbp_sex_inter,
                  'IMT_f' = imt_f, 'Dis_f' = dis_f, 'SBP_f' = sbp_f, 'DBP_f' = dbp_f,
                  'IMT_m' = imt_m, 'Dis_m' = dis_m, 'SBP_m' = sbp_m, 'DBP_m' = dbp_m)) }
  
  if ('regional' %in% outcs) { modls <- append(modls,
             list('IMT_subc_vol' = imt_subc_vol, 'Dis_subc_vol' = dis_subc_vol, 'SBP_subc_vol' = sbp_subc_vol, 'DBP_subc_vol' = dbp_subc_vol,
                  'IMT_tract_FA' = imt_tract_fa, 'Dis_tract_FA' = dis_tract_fa, 'SBP_tract_FA' = sbp_tract_fa, 'DBP_tract_FA' = dbp_tract_fa,
                  'IMT_tract_MD' = imt_tract_md, 'Dis_tract_MD' = dis_tract_md, 'SBP_tract_MD' = sbp_tract_md, 'DBP_tract_MD' = dbp_tract_md)) }
  
  if ('main_10y' %in% outcs) { modls <- append(modls,
              list('IMT_brain10' = imt_brain10, 'Dis_brain10' = dis_brain10, 'SBP_brain10' = sbp_brain10, 'DBP_brain10' = dbp_brain10)) }
  
  if ('other' %in% outcs) { modls <- append(modls,
              list('IMT_other' = imt_other, 'Dis_other' = dis_other, 'SBP_other' = sbp_other, 'DBP_other' = dbp_other)) }
  
  if (BPIMT_interaction) { modls <- append(modls, list('Inter_IMT-BP' = inter)) }
  
  openxlsx::write.xlsx(modls,
                       file = file.path(respath, paste0('Results_',NAME_ANALYSIS,sample,'_',daterun,'.xlsx')),
                       overwrite = T)

  # ============================================================================
  # LONGITUDINAL BRAIN ANALYSIS 
  # ============================================================================
  if (longit) {
    # Transform mice object to long format 
    long_imp <- data.frame()
    for (i in 0:20) {
      d <- mice::complete(imput, i)
      long_d <- reshape(data=d, idvar='IDC',
                        varying = list(c('age_mri_10','age_mri_13'),
                                       c('tbv_10','tbv_13'), c('gmv_10','gmv_13'),
                                       c('mfa_10','mfa_13'), c('mmd_10','mmd_13'),
                                       c('tiv_10','tiv_13'), c('wmv_10','wmv_13'),
                                       c('cortex_10','cortex_13'), c('subcort_10','subcort_13'),
                                       c('ventrix_10','ventrix_13'),
                                       c('accumbens_10','accumbens_13'),
                                       c('amygdala_10', 'amygdala_13'),
                                       c('caudate_10', 'caudate_13'),
                                       c('hippocampus_10','hippocampus_13'),
                                       c('pallidum_10', 'pallidum_13'),
                                       c('putamen_10', 'putamen_13'),
                                       c('thalamus_10','thalamus_13')),
                        v.name=c('AGE','TBV','GMV','MFA','MMD','TIV','WMV','CRT','SBC','VNT',
                                 'ACC', 'AMY', 'CAU', 'HIP', 'PAL', 'PUT', 'THA'),
                        timevar='Visit',
                        direction='long')
      long_d['.imp'] <- i
      long_imp <- rbind(long_imp, long_d)
    }
    
    # Transform back to mm (?)
    # tomm = c('TBV','GMV','TIV','WMV','CRT','SBC','ACC','AMY','CAU','HIP','PAL','PUT','THA')
    # long_imp[,tomm] <- long_imp[,tomm]*1000
    
    implong <- as.mids(long_imp)
    
    # imp0 <- complete(implong, 0)
    # write.csv(imp0, file.path(datpath, paste0('Data_long_',sample,'.csv')))

  # Fit the models and extract model paramenters as well as simple slopes
  fit_long <- function(e, o) {
    covs1 <- '+ sex + height_10'

    mod <- data.frame()
    for(adj in c('base','full')) {
      cat('\n-----------------------------------------------------------------------\n',
          o,' ~ ',e,' (',adj,' model)','\n-----------------------------------------------------------------------\n')
      # Define covariates
      if (adj=='base') { covs2 <- '' # reduce covariates to base only
      } else { covs2 <- '+ ethn_dich + bmi_10_z + tot_educ_6 + m_age' }
      # Construct formula
      f = paste(o,'~ AGE *',e,covs1,covs2,'+ (1|IDC)')
      # fit the model
      # fit <- lmerTest::lmer(f, data = implong, REML = F) #in complete data
      fit <- with(implong, lmerTest::lmer(as.formula(f), data = data.frame(mget(ls())), REML = F))
      sum <- summary(pool(fit), conf.int = TRUE)
      # sum = as.data.frame(coef(summary(fit)))
      extract_REsd <- function(model){
        var <- lapply(model$analyses, VarCorr)
        out <- mean(unlist(lapply(var, function(x) attr(x$IDC, 'stddev'))))
        return(out)
      }
      # Add AIC and SD of random effect
      sum$REsd <- extract_REsd(fit)
      sum$AIC <- mean(sapply(fit$analyses, AIC))

      # display the model
      # print(sum)
      cat('\n-----------------------------------------------------------------------\n')
      # Round everything
      sum[,-1] <- round(sum[,-1], 3)
      # Row names as column
      sum <- data.frame(model = o, sum)
      # And one space in between models
      sum[nrow(sum)+1,] <- NA
      # Stack base and fully adjusted models
      mod <- rbind(mod, sum)
      rownames(mod) <- NULL
    }
    return(mod)
  }

  exp_long <- function(exp, outcomes = c('TBV', 'GMV', 'MFA', 'MMD','WMV','CRT','SBC',
                                         'ACC', 'AMY', 'CAU', 'HIP', 'PAL', 'PUT', 'THA')) {
    mod <- data.frame()
    for (out in outcomes) {
      p <- fit_long(exp, out)
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

    #mod$lci <- (mod$Estimate - 1.96*mod$Std..Error)
    #mod$uci <- (mod$Estimate + 1.96*mod$Std..Error)

    return(mod)
  }

  for (e in expos) {
    ze <- paste0(e,'_10_z')
    assign(paste0(e,'_long'), exp_long(ze))
  }
  
  modls_long <- list('imt_long' = imt_long, 'dis_long' = dis_long,'sbp_long' = sbp_long,'dbp_long' = dbp_long)
  openxlsx::write.xlsx(modls_long,
                       file = file.path(respath, paste0('LongRes_',sample,'_',daterun,'.xlsx')),
                       overwrite = T)
  }
}