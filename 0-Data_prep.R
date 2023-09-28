# ==============================================================================
# ======================= 0. Data prep & imputation ============================
# ==============================================================================

# Required packages
invisible(lapply(c('foreign','dplyr','mice','miceadds'), require, character.only = T));

# Define paths 
genrpath <- dirname(file.choose()) # <== choose project folder
datapath <- file.path(genrpath,'DATA') # data folder
respath  <- file.path(genrpath,'results_190923') # results folder

date <- format(Sys.Date(), "%d%m%y")

readsav <- function(file, summary=F) {
  d <- foreign::read.spss(file.path(datapath, file), use.value.labels = T, to.data.frame = T)
  # make sure missing values are read in correctly
  for (col in d) { if (max(as.numeric(col), na.rm=T) == 999) { d[!is.na(col) & col == 999,] <- NA } }
  if (summary==T) { print(summary(d)) } 
  return(d)
}

# LOAD DATA ====================================================================
# bp2 <- readsav('CHILDBLOODPRESSURE_05102010.sav')
bp5 <- readsav('CHILDBLOODPRESSURE5_10022012.sav') # ageChildYF5, MeanSBP_ex1_5child, MeanDBP_ex1_5child
bp9 <- readsav('CHILDBLOODPRESSURE9_21042016.sav') # agechild9_visit1, MEANSBP_ex1_child9, MEANDBP_ex1_child9

dis <- readsav('CHILDF9_IMT_Distensibility_07122020.sav') # agechild9_visit1, DIS_mean_combined
imt <- readsav('CHILDIMTF9_23112020.sav') # agechild9_visit1, IMT_mean_combined

bmi5  <- readsav('CHILDGROWTH5_10122014.sav') # agey5child, sdsbmiforage
bmi9  <- readsav('CHILDGROWTH9_06072021.sav') # agechild9, heightchild9, bmichild9, sdsheight9childF, sdsbmiforage9childF
bmi13 <- readsav('CHILDGROWTH13_10122020.sav') # AGECHILD13, heightchild13, bmichild13, sdsbmiforage13childT, sdsheight13childT

gen   <- readsav('CHILD-ALLGENERALDATA_07072020.sav',summary=F)
ethn  <- read.csv(file.path(datapath, 'Ethn_rec.csv'), na.strings='')[,-1]
m_bmi <- readsav('MOTHERANTHROPOMETRY_18022013.sav')

brain <- read.csv(file.path(datapath, 'Brain.csv'))[,-1] # see brain merging script
brain[,sapply(brain, class)=='integer'] <- sapply(brain[,sapply(brain, class)=='integer'], as.numeric); # transform integers to numeric
names(brain)[1] <- 'IDC' # fix lower case idc

# Merge
dataset <- merge(gen, m_bmi, by='MOTHER')
dataset <- Reduce(function(x,y) merge(x = x, y = y, by = 'IDC', all.x = T),
                      list(dataset, ethn, bp5, bp9, dis, imt, bmi5, bmi9, bmi13, brain) ) 
rm(gen, ethn, bp5, bp9, dis, imt, bmi5, bmi9, bmi13, brain)

# Select only children alive at birth 
dataset <- dataset[dataset$OUTCOMECHILD=='live birth',]

# Select, order and transform needed variables 
data <- data.frame('IDC' = dataset$IDC,
                   'sex' = as.factor(dataset$GENDER), # 1 = boy; 2 = girl.
             'ethn_cont' = as.factor(dataset$ethn_genr), # recoded into 11 categories
                 'm_age' = dataset$AGE_M_v2, # maternal age at intake
              'm_educ_6' = dataset$EDUCM5,   # maternal education @6
             'height_10' = dataset$heightchild9, # in cm
                'bmi_10' = dataset$bmichild9,
              'bmi_10_z' = dataset$sdsbmiforage9childF,
                'age_10' = dataset$agechild9_visit1,
             # ------ EXPOSURES ------- #
                'imt_10' = dataset$IMT_mean_combined,
                'dis_10' = dataset$DIS_mean_combined,
                'sbp_10' = dataset$MEANSBP_ex1_child9, 
                'dbp_10' = dataset$MEANDBP_ex1_child9,
             # -------- BRAIN -------- #
            'age_mri_10' = dataset$age_child_mri_f09,
                'tbv_10' = dataset$genr_tbv_f09/1000, # total brain volume, in cm3
                'gmv_10' = dataset$TotalGrayVol_f09/1000, # total grey matter volume, in cm3
                'mfa_10' = dataset$mean_FA_genr_f09, # FA
                'mmd_10' = dataset$mean_MD_genr_f09*1000, # MD scaled for convergence 
                'tiv_10' = dataset$eTIV_f09/1000, # total intracranial volume, in cm3
            
            'age_mri_13' = dataset$age_child_mri_f13,
                'tbv_13' = dataset$genr_tbv_f13/1000, # total brain volume, in cm3
                'gmv_13' = dataset$TotalGrayVol_f13/1000, # total grey matter volume, in cm3
                'mfa_13' = dataset$mean_FA_genr_f13, # FA
                'mmd_13' = dataset$mean_MD_genr_f13*1000, # MD scaled for convergence 
                'tiv_13' = dataset$eTIV_f13/1000, # total intracranial volume, in cm3
            
              # ----- BRAIN (other) ----- #
          'accumbens_10' = dataset$Accumbens_area_vol_f09/1000, # in cm3
           'amygdala_10' = dataset$Amygdala_vol_f09/1000, # in cm3
            'caudate_10' = dataset$Caudate_vol_f09/1000, # in cm3
        'hippocampus_10' = dataset$Hippocampus_vol_f09/1000, # in cm3
           'pallidum_10' = dataset$Pallidum_vol_f09/1000, # in cm3
            'putamen_10' = dataset$Putamen_vol_f09/1000, # in cm3
           'thalamus_10' = dataset$Thalamus_Proper_vol_f09/1000, # in cm3
          'accumbens_13' = dataset$Accumbens_area_vol_f13/1000, # in cm3
           'amygdala_13' = dataset$Amygdala_vol_f13/1000, # in cm3
            'caudate_13' = dataset$Caudate_vol_f13/1000, # in cm3
        'hippocampus_13' = dataset$Hippocampus_vol_f13/1000, # in cm3
           'pallidum_13' = dataset$Pallidum_vol_f13/1000, # in cm3
            'putamen_13' = dataset$Putamen_vol_f13/1000, # in cm3
           'thalamus_13' = dataset$Thalamus_Proper_vol_f13/1000, # in cm3
        
             'cgc_FA_10' = dataset$cgc_FA_f09, # Cincgulate gyrus 
             'cgc_MD_10' = dataset$cgc_MD_f09*1000, 
             'cst_FA_10' = dataset$cst_FA_f09, # Cortico-spinal tract
             'cst_MD_10' = dataset$cst_MD_f09*1000, 
             'unc_FA_10' = dataset$unc_FA_f09, # uncinate fasciculus
             'unc_MD_10' = dataset$unc_MD_f09*1000, 
             'ilf_FA_10' = dataset$ilf_FA_f09, # inferior longitudinal fasciculus
             'ilf_MD_10' = dataset$ilf_MD_f09*1000, 
             'slf_FA_10' = dataset$slf_FA_f09, # superior longitudinal fasciculus
             'slf_MD_10' = dataset$slf_MD_f09*1000, 
             'fma_FA_10' = dataset$fma_FA_f09, # major forceps 
             'fma_MD_10' = dataset$fma_MD_f09*1000,
             'fmi_FA_10' = dataset$fmi_FA_f09, # minor forceps 
             'fmi_MD_10' = dataset$fmi_MD_f09*1000,
             'cgc_FA_13' = dataset$cgc_FA_f13, # Cincgulate gyrus 
             'cgc_MD_13' = dataset$cgc_MD_f13*1000, 
             'cst_FA_13' = dataset$cst_FA_f13, # Cortico-spinal tract
             'cst_MD_13' = dataset$cst_MD_f13*1000, 
             'unc_FA_13' = dataset$unc_FA_f13, # uncinate fasciculus
             'unc_MD_13' = dataset$unc_MD_f13*1000, 
             'ilf_FA_13' = dataset$ilf_FA_f13, # inferior longitudinal fasciculus
             'ilf_MD_13' = dataset$ilf_MD_f13*1000, 
             'slf_FA_13' = dataset$slf_FA_f13, # superior longitudinal fasciculus
             'slf_MD_13' = dataset$slf_MD_f13*1000, 
             'fma_FA_13' = dataset$fma_FA_f13, # major forceps 
             'fma_MD_13' = dataset$fma_MD_f13*1000,
             'fmi_FA_13' = dataset$fmi_FA_f13, # minor forceps 
             'fmi_MD_13' = dataset$fmi_MD_f13*1000,
        
             'cortex_10' = dataset$CortexVol_f09/1000, # in cm3
            'subcort_10' = dataset$SubCortGrayVol_f09/1000, # in cm3
                'wmv_10' = dataset$CerebralWhiteMatterVol_f09/1000, # in cm3
                'csf_10' = dataset$CSF_vol_f09/1000, # in cm3
            'ventrix_10' = rowSums(dataset[, grepl('Ventricle_vol_f09', names(dataset))], na.rm=FALSE)/1000, # in cm3
              'bstem_10' = dataset$Brain_Stem_vol_f09/1000, # in cm3
      'crbellum_cort_10' = dataset$Cerebellum_Cortex_vol_f09/1000, # in cm3
       'crbellum_wmv_10' = dataset$Cerebellum_White_Matter_vol_f09/1000, # in cm3
      #   'mean_thick_10' = dataset$MeanThickness_f09/1000, # in cm3
             'cortex_13' = dataset$CortexVol_f13/1000, # in cm3,
            'subcort_13' = dataset$SubCortGrayVol_f13/1000, # in cm3
                'wmv_13' = dataset$CerebralWhiteMatterVol_f13/1000, # in cm3
                'csf_13' = dataset$CSF_vol_f13/1000, # in cm3
            'ventrix_13' = rowSums(dataset[, grepl('Ventricle_vol_f13', names(dataset))], na.rm=FALSE)/1000, # in cm3
              'bstem_13' = dataset$Brain_Stem_vol_f13/1000, # in cm3
      'crbellum_cort_13' = dataset$Cerebellum_Cortex_vol_f13/1000, # in cm3
       'crbellum_wmv_13' = dataset$Cerebellum_White_Matter_vol_f13/1000, # in cm3
       #  'mean_thick_13' = dataset$MeanThickness_f13/1000, # in cm3
          
              # ------ OTHER ------- #
       'm_bmi_prepregn' = dataset$BMI_0,   # self-reported Maternal BMI (used for imputation)
               'parity' = dataset$PARITY,  # parity (used for imputation)
       'gest_age_birth' = dataset$GESTBIR, # gestational age at birth (used for imputation)
          'gest_weight' = dataset$WEIGHT,  # gestational weight (used for imputation)
         'm_educ_pregn' = dataset$EDUCM,   # maternal education pregnancy (used for imputation)
             'm_educ_3' = dataset$EDUCM3,  # maternal education @3 (used for imputation)
             'p_educ_6' = dataset$EDUCP5,  # maternal education @3 (used for imputation)
             'income_6' = dataset$INCOME5, # household income @6 (used for imputation)
              'm_bmi_6' = dataset$BMIMotherF5, # measured (used for imputation)
                'sbp_6' = dataset$MeanSBP_ex1_5child, # (used for imputation)
                'dbp_6' = dataset$MeanDBP_ex1_5child, # (used for imputation)
             'height_6' = dataset$length5child, # (used for imputation)
              'bmi_6_z' = dataset$sdsbmiforage, # (used for imputation)
            
                # ------ SELECTION ------- #
              'visit_5' = dataset$visitChildF5,
             'visit_10' = dataset$visit1.F9,
                 'm_ID' = dataset$MOTHER, # mother id used to identify siblings (for exclusion)
                 'twin' = as.factor(dataset$TWIN), # (for exclusion)
       'mri_consent_10' = as.factor(dataset$mri_consent_f09), # consent status
        'mri_braces_10' = as.factor(dataset$has_braces_mri_f09), # has braces
      'mri_inc_find_10' = as.factor(dataset$exclude_incidental_f09), # incidental findings
       't1_scantype_10' = as.factor(dataset$t1_has_nii_f09), # scan type
             't1_qc_10' = as.factor(dataset$freesurfer_qc_f09), # QC
      'dti_scantype_10' = as.factor(dataset$dti_has_nii_f09), # scan type
            'dti_qc_10' = as.factor(ifelse(is.na(dataset$dti_overall_qc_f09), 'usable','unusable')), # QC
       'mri_consent_13' = as.factor(dataset$mri_consent_f13), # consent status
        'mri_braces_13' = as.factor(dataset$has_braces_mri_f13), # has braces
      'mri_inc_find_13' = as.factor(dataset$exclude_incidental_f13), # incidental findings
       't1_scantype_13' = as.factor(dataset$t1_has_nii_f13), # scan type
             't1_qc_13' = as.factor(dataset$freesurfer_qc_f13), # QC
      'dti_scantype_13' = as.factor(dataset$dti_has_nii_f13), # scan type
            'dti_qc_13' = as.factor(ifelse(is.na(dataset$dti_overall_qc_f13), 'usable','unusable'))) # QC
        
# summary(data)
# ==============================================================================
# Organize variable names if groups
comb_timepoints <- function(names) { vec <- c(paste0(names, '_10'), paste0(names, '_13')); return(vec) }

exposures  <- paste0(c('imt', 'dis', 'sbp', 'dbp'),'_10')
outcomes   <- comb_timepoints(c('tbv', 'gmv', 'mfa', 'mmd'))
covariates <- c('sex', 'age_mri_13', 'height_10', 'ethnicity', 'bmi_10_z', 'm_educ_6', 'm_age')
subc_vol   <- comb_timepoints(c('accumbens', 'amygdala', 'caudate', 'hippocampus', 'pallidum', 'putamen', 'thalamus'))
wm_trcts   <- comb_timepoints(c('cgc_FA', 'cgc_MD', 'cst_FA', 'cst_MD', 'unc_FA', 'unc_MD', 'ilf_FA', 'ilf_MD', 'slf_FA', 'slf_MD', 'fma_FA', 'fma_MD', 'fmi_FA', 'fmi_MD'))
othr_brain <- comb_timepoints(c('tiv','cortex','subcort','wmv','ventrix','csf','bstem','crbellum_cort','crbellum_wmv')) # 'mean_thick'
selection  <- c('visit_5','visit_10', comb_timepoints(c('mri_consent', 'mri_braces', 'mri_inc_find', 't1_scantype', 't1_qc', 'dti_scantype', 'dti_qc')))
auxiliary  <- c('age_10','gest_weight','gest_age_birth','parity', 'm_bmi_prepregn', 'm_bmi_6','income_6','sbp_6','dbp_6','height_6', 'bmi_6_z',
                'm_educ_pregn','m_educ_3','p_educ_6','age_mri_10','m_ID', 'twin', othr_brain)

write.csv(data, file.path(datapath, 'Data_dirty_brain.csv'))
# ==============================================================================
# create clean brain variables for imputation
clean_brain <- function(set) {
  if (substring(set[1],1,3)=='tbv') { scan = 't1_' } else { scan = 'dti_'}
  for (im in set) {
    if (substring(im, nchar(im)-1)=='_10') { age = '_10'} else { age = '_13'}
    message(im)
    cat(nrow(data[!is.na(data[,im]),]))
    data[,im] <- ifelse(data[,paste0('mri_consent',age)]!='yes' |  
                        data[,paste0('mri_braces',age)]!='no' | 
                        data[,paste0('mri_inc_find',age)]!='include' | 
                        is.na(data[,paste0(scan,'scantype',age)]) | 
                        data[,paste0(scan,'qc',age)]!='usable', NA, data[,im])
    cat(' -->', nrow(data[!is.na(data[,im]),]), '\n')
    # cat(im, ' - ', nrow(data[!is.na(data[,im]),]), '\n')
  }
  return(data)
}
data <- clean_brain(c('tbv_10','gmv_10','tbv_13','gmv_13', subc_vol, othr_brain))
data <- clean_brain(c('mfa_10','mmd_10','mfa_13','mmd_13', wm_trcts))

# summary(data)
# ==============================================================================
# Save histograms of all variables in the full sample
pdf(file.path(respath, paste0('0-hist_orig_sample_',date,'.pdf')))
for (var in names(data)) { hist(as.numeric(data[,var]), col = 'blue', main = toupper(var), xlab = var) }
dev.off() 

# ==============================================================================
# Save dataset
write.csv(data, file.path(datapath, 'Data.csv'))

# Check selection and missing rate =============================================
'%notin%' <- Negate('%in%')

chop <- function(dset, var, oper=`==`, value=NULL, 
                 mult=NULL, var2=NULL, oper2=NULL, value2=NULL) {
  # Print summaries for checking 
  cat(var,'\n'); if(var!='IDC') { print(summary(as.factor(dset[,var]))) }
  if (!is.null(var2)) { cat('\n',var2,'\n'); print(summary(as.factor(dset[,var2]))) }
  # Simple (one variable) selection condition
  # Note: as above, not all operations require a value argument (e.g. is.na)
  if (!is.null(value)) { cond <- oper(dset[,var], value) } else { cond <- oper(dset[,var]) }
  # Multiple selection conditions 
  if (!is.null(mult)) {
    if (!is.null(value2)) { cond2 <- oper2(dset[,var2], value2) } else { cond2 <- oper2(dset[,var2]) }
    cond <- mult(cond, cond2) # paste conditions together
  }
  dat = subset(dset, cond) 
  
  # Display log of selection
  nsel <- nrow(dat); loss <- nrow(dset) - nsel
  cat(nrow(dset),' - ', loss, ' = ', nsel, '\n\n')
  
  return(dat)
}

select_sibling <- function(dt, column_selection = c(), random = F, seed = 31081996, mother_id = 'm_ID', child_id = 'IDC') {
  # if no selection is specified, missingness in the entire dataframe is used
  if (length(column_selection) > 0) { dt <- dt[, c(child_id, mother_id, column_selection)] } 
  # First randomly shuffle the dataset 
  set.seed(seed)
  dt <- dt[sample(nrow(dt)),]
  # Get rid of empty NA values for mother
  dt <- dt[!is.na(dt[,mother_id]),]
  # Determine a list of children that have a sibling in the set
  if (random==T) { 
    sibling_ids <- dt[duplicated(dt[, mother_id]), child_id] # i.e.  which mother IDs recur more than once
  } else {
    dt$missing <- rowSums(is.na(dt)) # compute how many missing values in the columns of interest 
    dt_ord <- dt[order(dt$missing),] # order based on number of missing 
    sibling_ids <- dt_ord[duplicated(dt_ord[, mother_id]), child_id] # selection
  }
  # message(length(sibling_ids), ' siblings identified.')
  return(sibling_ids)
}

apply_selection <- function(set, age='13', start=data) {
  cat('\n--------------------------------------------\n',set,
      '\n--------------------------------------------\n')
  step <- chop(start,paste0('mri_consent_', age), `==`,'yes')
  step <- chop(step, paste0('mri_braces_', age),  `==`,'no')
  step <- chop(step, paste0('mri_inc_find_', age),`==`,'include')
  step <- chop(step, paste0(set,'_scantype_', age),`==`,'yes')
  step <- chop(step, paste0(set,'_qc_', age),`==`,'usable')
  step <- chop(step,'twin',`==`,'No')
  if (set=='t1') { outc = paste0(c('tbv_','gmv_'), age) } else { outc = paste0(c('mfa_','mmd_'), age) }
  outp <- chop(step,'IDC',`%notin%`,select_sibling(step, column_selection=c(exposures,outc))) 
  return(outp)
}

sink(file.path(respath, paste0('0-flowchart_selection_',date,'.txt')))
# select sample for fully imputed set ------------------------------------------
notwins <- chop(data,'twin',`==`,'No')
fullsmp <- chop(notwins,'IDC',`%notin%`,select_sibling(notwins, column_selection=c(exposures,outcomes))) 

include_full <- fullsmp$IDC

basesmp5 <- chop(notwins,'visit_5',`==`,'Yes')
basesmp5 <- chop(basesmp5,'IDC',`%notin%`,select_sibling(basesmp5, column_selection=c(exposures,outcomes))) 

include_base5 <- basesmp5$IDC

basesmp9 <- chop(notwins,'visit_10',`==`,'Yes')
basesmp9 <- chop(basesmp9,'IDC',`%notin%`,select_sibling(basesmp9, column_selection=c(exposures,outcomes))) 

include_base9 <- basesmp9$IDC

# select sample for sensitivity analyses (complete outcome) --------------------
mri <- apply_selection('t1')
dti <- apply_selection('dti')

# NOTE: Get rid of additional missing in dti
check <- dti[is.na(dti$mfa_13),]
write.csv(check, 'check_missing_dti.csv')
dti <- dti[!is.na(dti$mfa_13),]; cat('\nAdditional missing exclusion: ',nrow(dti))

# Save list of selected IDs for later selection
include_mri <- mri$IDC
include_dti <- dti$IDC

sink()

# ------------------------------------------------------------------------------
# compare full sample with complete outcome samples. 
compare_samps <- function(var, subsample) {
  subsmp_name <- deparse(substitute(subsample))
  if (!is.null(levels(fullsmp[,var]))) { # this is a factor, need to compare proportions
    for (l in levels(fullsmp[,var])) {
      t = prop.test(x = c(nrow(fullsmp[fullsmp[var]==l,]), nrow(subsample[subsample[var]==l,])), 
                    n = c(nrow(fullsmp), nrow(subsample))) # n = c(nrow(fullsmp[!is.na(fullsmp[var]),]), nrow(subsample[!is.na(subsample[var]),]))) 
      if (t$estimate[1] > t$estimate[2]) { rel = 'full > selected' } else { rel = 'full < selected' }
      if (!is.na(t$p.value) & t$p.value < 0.05) { cat('\n',subsmp_name,'-',var,'-', l,'=',rel, round(t$p.value,3)) }
      }
    } else {
    t = t.test(fullsmp[var], subsample[var], var.equal=TRUE)
    if (t$estimate[1] > t$estimate[2]) { rel = 'full > selected' } else { rel = 'full < selected' }
    if (!is.na(t$p.value) & t$p.value < 0.05) { cat('\n',subsmp_name,'-',var, '=',rel, round(t$p.value,3)) }
    }
}
for (v in names(fullsmp)[-which(names(fullsmp)%in%c('IDC','m_ID','twin',selection))]) {
  compare_samps(v, basesmp5)
  #compare_samps(v, basesmp9)
  #compare_samps(v, mri)
  #compare_samps(v, dti)
}

# ------------------------------------------------------------------------------

# Calculate the percentage missing data ----------------------------------------
calc_miss <- function(ds) {  return(paste0(colSums(is.na(ds)),' (',round((colSums(is.na(ds))/nrow(ds))*100, 1), '%)'))
}
calc_name <- function(ds, name) { return(paste0(name, ' (N = ',nrow(ds), ')')) }

miss <- data.frame(calc_miss(fullsmp),
                   calc_miss(basesmp5),
                   calc_miss(basesmp9),
                   calc_miss(mri),
                   calc_miss(dti), row.names = names(data))
names(miss) = c(calc_name(fullsmp, 'full sample'),
                calc_name(basesmp5,'baseline 5y'),
                calc_name(basesmp9,'baseline 9y'),
                calc_name(mri, 'MRI sample'),
                calc_name(dti, 'DTI sample'))
# View(miss)
write.csv(miss, file.path(respath, paste0('0-miss_pattern_',date,'.csv')))

# Create correlation matrix (full sample) --------------------------------------
c <- round(cor(basesmp9[,-which(names(basesmp9) %in% c('IDC', selection, # binary and categorical
                                                     'sex','ethn_cont','twin','income_6',
                                                     'm_educ_pregn','m_educ_3','m_educ_6','p_educ_6'))], 
               use='pairwise.complete.obs'), 2)
write.csv(c, file.path(respath, paste0('0-corr_matrix_base9_',date,'.csv')))

# IMPUTATION ===================================================================
data <- data[ , -which(names(data) %in% selection)] # get rid of selection variables 

# Use random forest for continuous data, logistic regression for binomial (sex) and polyreg for categorical(ethnicity, m_educ_6)
meth <- make.method(data, defaultMethod = c("rf", "logreg", "polyreg"))

# Random forest imputation ran in parallel 
imp_rf <- futuremice(data, method = meth, m = 20, maxit = 40, ntree = 10, rfPackage = 'ranger',
                    n.core = 5, n.imp.core = 4, parallelseed = 310896, print=T) 

# POST-PROCESSING & QUALITY CONTROL ============================================ 

# Add dichotomized ethnicity and factor/numeric versions of categocal variables
long.impdata <- complete(imp_rf, 'long', include = TRUE) %>%
  mutate(ethn_dich = if_else(ethn_cont %in% c('Dutch','Europe'), 'European', 'non-European')) %>%
  mutate(m_educ_cont = as.numeric(m_educ_6)) %>%
  mutate(p_educ_cont = as.numeric(p_educ_6)) %>%
  mutate(educ_6 = m_educ_cont + p_educ_cont) %>%
  mutate(income_cat = as.factor(income_6)) %>%
  mutate(age_gap_13 = age_mri_13 - age_10) %>%
  mutate(age_gap_10 = age_mri_10 - age_10)
# Convert back to mids
imp_rf <- as.mids(long.impdata)

# Save original set 
saveRDS(imp_rf, file.path(respath, paste0('imp_orig_',date,'.rds')))

# Post-processing: selection and scaling 
postprocess <- function(inclusion, filename) {
  # Subset
  samp <- miceadds::subset_datlist(imp_rf, subset = imp_rf$data$IDC %in% inclusion)
  # Standardize
  mains <- c(exposures, outcomes, subc_vol, wm_trcts, othr_brain)
  sampz <- miceadds::scale_datlist(samp, orig_var = mains, trafo_var = paste0(mains, '_z'))
  sampimp <- miceadds::datlist2mids(sampz)
  # Save dataset
  saveRDS(sampimp, file.path(respath, paste0(filename,'_',date,'.rds')))
  
  # Save QC pdf 
  impqc <- function(dset, name) {
    # dir.create(file.path(respath,'QC-imputation'))
    pdf(file.path(respath,'QC-imputation', paste0(name,'_',date,'.pdf')))
    for (v in names(data)) { if (nrow(dset$imp[[v]]) > 1) {
      print(densityplot(dset, as.formula(paste('~',v)))) }
    }
    dev.off()
  }
  impqc(sampimp, filename)
  
  return(sampimp)
}

fullset <- postprocess(include_full, 'imp_full')
bas5set <- postprocess(include_base5,'imp_base5')
bas9set <- postprocess(include_base9,'imp_base9')
mriset  <- postprocess(include_mri,  'imp_smri')
dtiset  <- postprocess(include_dti,  'imp_dti')

# Save list of datasets to use as input for QDECR ------------------------------
main_sample = bas9set

dl = list()
dir.create(file.path(respath,'df_by_imp'))

for (n in 0:20) {
  d <- mice::complete(main_sample, action = n)
  write.csv(d, file.path(respath, 'df_by_imp', paste0('Data_imp',n,'.csv')))
  
  if (n>0) {
    ds <- data.frame(complete(mriset, n))
    dl[[n]] <- ds
  }
}
saveRDS(dl, file.path(respath, paste0('implist_qdecr_',date,'.rds')))

# IMPUTATION (old method) ======================================================

# # Quick prediction -------------------------------------------------------------
# qp <- mice::quickpred(data, mincor=0.2, exclude ='IDC')
# # Do not impute the outcome variables or complete ones
# # qp[c(outcomes, subc_vol, wm_trcts, 'tiv_13', 'sex', 'age_mri_13', 'm_age'), ] <- 0
# 
# write.csv(qp, file.path(respath,'QC-imputation','Qpred.csv'))
# 
# # Create an empty predictor matrix to fill with instructions -------------------
# pm <- matrix(0, ncol(data), ncol(data), dimnames=list(names(data),names(data)))

# pm[, c(exposures, outcomes, wm_trcts, subc_vol)]<- 1

# pm['imt_10', c(exposures, covariates[!covariates=='age_mri_13'], 'age_10','gest_weight')]<- 1
# pm['dis_10', c(exposures, covariates[!covariates=='age_mri_13'], 'sbp_6','gest_weight')] <- 1
# pm['sbp_10', c(exposures, covariates[!covariates=='age_mri_13'], 'sbp_6','dbp_6','m_bmi_prepregn','gest_weight')] <- 1
# pm['dbp_10', c(exposures, covariates[!covariates=='age_mri_13'], 'sbp_6','dbp_6','m_bmi_prepregn','income_6', 'gest_weight')] <- 1
# pm['height_10', c(exposures, covariates[!covariates%in%c('age_mri_13','bmi_10_z','tiv_13')], auxiliary[!auxiliary%in%c('bmi_6_z','m_bmi_6')])] <- 1
# pm['ethnicity',c(exposures, covariates[!covariates%in%c('age_mri_13','bmi_10_z')], auxiliary[!auxiliary%in%c('age_10','m_bmi_prepregn')])] <- 1
# pm['bmi_10_z',  c(exposures, covariates[!covariates%in%c('age_mri_13','height_10','tiv_13')],auxiliary[!auxiliary%in%c('height_6','m_bmi_6')])] <- 1
# pm['m_educ_6', c(exposures, covariates[!covariates%in%c('height_10','tiv_13')], auxiliary[!auxiliary%in%c('age_10','m_bmi_prepregn')])] <- 1
# pm[auxiliary, auxiliary] <- 1
# 
# diag(pm) <- 0 # Make sure no circularity arises
# 
# write.csv(pm, file.path(respath,'QC-imputation','PredM.csv'))
# 
# # Check differences between quickpred and our matrix ---------------------------
# 
# # check_pm <- function(var) {
# #   vq <- names(which(qp[var,]==1)); # print(vq)
# #   vs <- names(which(pm[var,]==1)); # print(vs)
# #   
# #   print(setdiff(vq, vs))     
# #   print(setdiff(vs, vq))   
# #   
# #   vs = vs[!vs %in% c('sex','ethnicity')]
# #   m = round(cor(data[,vs], use='pairwise.complete.obs'),2)
# # 
# #   maxes <- apply(matrix(m[row(m)!=col(m)], ncol=ncol(m)), 2, max)
# #   mx <- data.frame(colnames(m), maxes)
# #   # View(m); View(mx)
# # }
# 
# # Run imputation ---------------------------------------------------------------
# imputation <- mice::mice(data, m = 30, # nr of imputed datasets
#                          maxit = 60, # nr of iterations
#                          seed = 310896, # set a seed for the random number generation
#                          predictorMatrix = pm, 
#                          method = meth)
# ==============================================================================
# END