# ==============================================================================
# ======================= 0. Data prep & imputation ============================
# ==============================================================================

# Required packages
invisible(lapply(c('foreign','dplyr','mice','miceadds'), require, character.only = T));

# Define paths 
genrpath <- dirname(file.choose())
datapath <- file.path(genrpath,'DATA') # data folder
respath  <- file.path(genrpath,'results') # results folder

readsav <- function(file, summary=F) {
  d <- foreign::read.spss(file.path(datapath, file), use.value.labels = T, to.data.frame = T)
  # make sure missing are read in correctly
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

bmi5 <- readsav('CHILDGROWTH5_10122014.sav') # agey5child, sdsbmiforage
bmi9 <- readsav('CHILDGROWTH9_06072021.sav') # agechild9, heightchild9, bmichild9, sdsheight9childF, sdsbmiforage9childF
bmi13 <- readsav('CHILDGROWTH13_10122020.sav') # AGECHILD13, heightchild13, bmichild13, sdsbmiforage13childT, sdsheight13childT

gen <- readsav("CHILD-ALLGENERALDATA_07072020.sav",summary=F)
m_bmi <- readsav('MOTHERANTHROPOMETRY_18022013.sav')

brain <- read.csv(file.path(datapath, 'Brain.csv')) # see brain merging script
names(brain)[2] <- 'IDC' # fix lower case idc

# Merge
dataset <- merge(gen, m_bmi, by='MOTHER')
dataset <- Reduce(function(x,y) merge(x = x, y = y, by = 'IDC', all.x = T),
                      list(dataset, bp5, bp9, dis, imt, bmi5, bmi9, bmi13, brain) ) 
rm(gen, bp5, bp9, dis, imt, bmi5, bmi9, bmi13, brain)

# Select only children alive at birth 
dataset <- dataset[dataset$OUTCOMECHILD=='live birth',]

# Select and transform needed variables 
data <- data.frame('IDC' = dataset$IDC, 
                 'imt_9' = dataset$IMT_mean_combined,
                 'dis_9' = dataset$DIS_mean_combined,
                 'sbp_9' = dataset$MEANSBP_ex1_child9, 
                 'dbp_9' = dataset$MEANDBP_ex1_child9,
                'tbv_13' = dataset$genr_tbv_f13,
                'gmv_13' = as.numeric(dataset$TotalGrayVol_f13),
                'mfa_13' = dataset$mean_FA_genr_f13, # or all
                'mmd_13' = dataset$mean_MD_genr_f13*1000,
              # ------ COVARIATES ------- #                
                   'sex' = as.factor(dataset$GENDER), # 1 = boy; 2 = girl.
                'age_13' = dataset$age_child_mri_f13,
              'height_9' = dataset$heightchild9,
             'ethnicity' = as.factor(ifelse(dataset$ETHNINFv2 == 'Dutch', 'dutch', 
                                            ifelse(dataset$ETHNINFv2 %in% c('European','American,western','Oceanie'), 'european_descent',
                                                   ifelse(dataset$ETHNINFv2 %in% c('Moroccan','Turkish'), 'turkish_moroccan', 
                                                          ifelse(dataset$ETHNINFv2 %in% c('Dutch Antilles','Surinamese'), 'surinames_antillian',
                                                                 ifelse(dataset$ETHNINFv2 %in% c('Indonesian','Cape Verdian','African',
                                                                                                 'American, non western','Asian, western', 
                                                                                                 'Asian, non western'), 'other', NA)))))),
               'bmi_9_z' = dataset$sdsbmiforage9childF,
              'm_educ_6' = dataset$EDUCM5, # maternal education @5
                 'm_age' = dataset$AGE_M_v2, # maternal age at intake
                'tiv_13' = as.numeric(dataset$eTIV_f13), # total intracranial volume
               # ------ SECONDARY ------- #
             'accumbens' = dataset$Accumbens_area_vol_f13,
              'amygdala' = dataset$Amygdala_vol_f13,
               'caudate' = dataset$Caudate_vol_f13,
           'hippocampus' = dataset$Hippocampus_vol_f13,
              'pallidum' = dataset$Pallidum_vol_f13,
               'putamen' = dataset$Putamen_vol_f13,
              'thalamus' = dataset$Thalamus_Proper_vol_f13,
                'cgc_FA' = dataset$cgc_FA_f13, # Cincgulate gyrus 
                'cgc_MD' = dataset$cgc_MD_f13*1000, 
                'cst_FA' = dataset$cst_FA_f13, # Cortico-spinal tract
                'cst_MD' = dataset$cst_MD_f13*1000, 
                'unc_FA' = dataset$unc_FA_f13, # uncinate fasciculus
                'unc_MD' = dataset$unc_MD_f13*1000, 
                'ilf_FA' = dataset$ilf_FA_f13, # inferior longitudinal fasciculus
                'ilf_MD' = dataset$ilf_MD_f13*1000, 
                'slf_FA' = dataset$slf_FA_f13, # superior longitudinal fasciculus
                'slf_MD' = dataset$slf_MD_f13*1000, 
                'fma_FA' = dataset$fma_FA_f13, # major forceps 
                'fma_MD' = dataset$fma_MD_f13*1000,
                'fmi_FA' = dataset$fmi_FA_f13, # minor forceps 
                'fmi_MD' = dataset$fmi_MD_f13*1000,
               # ------ BRAIN @10 ------- #
                 'tbv_9' = dataset$genr_tbv_f09,
                 'gmv_9' = as.numeric(dataset$TotalGrayVol_f09),
                 'mfa_9' = dataset$mean_FA_genr_f09,
                 'mmd_9' = dataset$mean_MD_genr_f09*1000,
             'age_mri_9' = dataset$age_child_mri_f09,
                 'tiv_9' = as.numeric(dataset$eTIV_f09),
              'cortex_9' = as.numeric(dataset$CortexVol_f09),
             'subcort_9' = as.numeric(dataset$SubCortGrayVol_f09),
                 'wmv_9' = as.numeric(dataset$CerebralWhiteMatterVol_f09),
               # ------ OTHER ------- #
             'cortex_13' = as.numeric(dataset$CortexVol_f13),
            'subcort_13' = as.numeric(dataset$SubCortGrayVol_f13),
                'wmv_13' = as.numeric(dataset$CerebralWhiteMatterVol_f13),
        'gest_age_birth' = dataset$GESTBIR,   # gestational age at birth (used for imputation)
           'gest_weight' = dataset$WEIGHT,    # gestational weight (used for imputation)
                'parity' = dataset$PARITY,    # parity (used for imputation)
        'm_bmi_prepregn' = dataset$BMI_0, # self-reported Maternal BMI
               'm_bmi_6' = dataset$BMIMotherF5,
              'income_6' = as.numeric(dataset$INCOME5), # income @5
                 'sbp_6' = dataset$MeanSBP_ex1_5child, 
                 'dbp_6' = dataset$MeanDBP_ex1_5child,
              'height_6' = dataset$length5child,
               'bmi_6_z' = dataset$sdsbmiforage,
          'm_educ_pregn' = dataset$EDUCM,  # maternal education pregnancy
              'm_educ_3' = dataset$EDUCM3, # maternal education @3
                 'age_9' = dataset$agechild9_visit1,
                # ------ SELECTION ------- #
        'mri_consent_13' = as.factor(dataset$mri_consent_f13), # consent status
         'mri_braces_13' = as.factor(dataset$has_braces_mri_f13), # has braces
       'mri_inc_find_13' = as.factor(dataset$exclude_incidental_f13), # incidental findings
        't1_scantype_13' = as.factor(dataset$t1_has_nii_f13), # scan type
              't1_qc_13' = as.factor(dataset$freesurfer_qc_f13), # QC
       'dti_scantype_13' = as.factor(dataset$dti_has_nii_f13), # scan type
             'dti_qc_13' = as.factor(ifelse(is.na(dataset$dti_overall_qc_f13), 'usable','unusable')), # QC
         'mri_consent_9' = as.factor(dataset$mri_consent_f09), # consent status
          'mri_braces_9' = as.factor(dataset$has_braces_mri_f09), # has braces
        'mri_inc_find_9' = as.factor(dataset$exclude_incidental_f09), # incidental findings
         't1_scantype_9' = as.factor(dataset$t1_has_nii_f09), # scan type
               't1_qc_9' = as.factor(dataset$freesurfer_qc_f09), # QC
        'dti_scantype_9' = as.factor(dataset$dti_has_nii_f09), # scan type
              'dti_qc_9' = as.factor(ifelse(is.na(dataset$dti_overall_qc_f09), 'usable','unusable')), # QC
                  'm_ID' = dataset$MOTHER, # mother id used to identify siblings (for exclusion)
                  'twin' = as.factor(dataset$TWIN))
# summary(data)
# ==============================================================================
# Organize variable names if groups
exposures  <- paste0(c('imt', 'dis', 'sbp', 'dbp'),'_9')
outcomes   <- paste0(c('tbv', 'gmv', 'mfa', 'mmd'),'_13')
covariates <- c('sex', 'age_13', 'height_9', 'ethnicity', 'bmi_9_z', 'm_educ_6', 'm_age','tiv_13')
subc_vol   <- c('accumbens', 'amygdala', 'caudate', 'hippocampus', 'pallidum', 'putamen', 'thalamus')
tracts     <- c('cgc_FA', 'cgc_MD', 'cst_FA', 'cst_MD', 'unc_FA', 'unc_MD', 'ilf_FA', 'ilf_MD', 'slf_FA', 'slf_MD', 'fma_FA', 'fma_MD', 'fmi_FA', 'fmi_MD')       
outcomes_9 <- paste0(c('tbv', 'gmv', 'mfa', 'mmd'),'_9')
auxiliary  <- c('age_9','gest_weight','gest_age_birth','parity', 'm_bmi_prepregn', 'm_bmi_6','income_6','sbp_6','dbp_6','height_6', 'bmi_6_z',
                'm_educ_pregn','m_educ_3','age_mri_9','tiv_9','cortex_9','subcort_9','wmv_9','cortex_13','subcort_13','wmv_13','m_ID', 'twin')
selection  <- c('mri_consent_13', 'mri_braces_13', 'mri_inc_find_13', 't1_scantype_13', 't1_qc_13', 'dti_scantype_13', 'dti_qc_13',
                'mri_consent_9', 'mri_braces_9', 'mri_inc_find_9', 't1_scantype_9', 't1_qc_9', 'dti_scantype_9', 'dti_qc_9')

# ==============================================================================
# create clean brain variables for imputation
clean_brain <- function(set) {
  age =  # = _13 or _9
  if (substring(set[1],1,3)=='tbv') { scan = 't1_' } else { scan = 'dti_'}
  for (im in set) {
    if (substring(im, nchar(im)-1)=='_9') { age = '_9'} else { age = '_13'}
    data[,im] <- ifelse(data[,paste0('mri_consent',age)]!='yes' |  
                        data[,paste0('mri_braces',age)]!='no' | 
                        data[,paste0('mri_inc_find',age)]!='include' | 
                        is.na(data[,paste0(scan,'scantype',age)]) | 
                        data[,paste0(scan,'qc',age)]!='usable', NA, data[,im])
    
    # cat(im, ' - ', nrow(data[!is.na(data[,im]),]), '\n')
  }
  return(data)
}
data <- clean_brain(c('tbv_13','gmv_13','tbv_9','gmv_9',subc_vol,
                      'tiv_13','cortex_13','subcort_13','wmv_13',
                      'tiv_9', 'cortex_9', 'subcort_9', 'wmv_9'))
data <- clean_brain(c('mfa_13','mmd_13','mfa_9','mmd_9',tracts))

# summary(data)
# ==============================================================================
# Save histograms of all variables in the full sample
pdf(file.path(respath,'hist_full_sample.pdf'))
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

mri <- apply_selection('t1')
dti <- apply_selection('dti')

# Get rid of additional missing in dti 
check <- dti[is.na(dti$mfa_13),]
write.csv(check, 'check_missing_dti.csv')
dti <- dti[!is.na(dti$mfa_13),]; nrow(dti)

# Save list of selected IDs for later selection
include_mri <- mri$IDC
include_dti <- dti$IDC

# select sample for cross-sectional analyses -----------------------------------
mri_long <- apply_selection('t1', age='9', start = mri)
dti_long <- apply_selection('dti',age='9', start = dti)

# Get rid of additional missing in dti 
dti_long <- dti_long[!is.na(dti_long$mfa_9),]; nrow(dti_long)
# Save list of selected IDs for later selection
include_mri_long <- mri_long$IDC
include_dti_long <- dti_long$IDC

# ------------------------------------------------------------------------------

# Calculate the percentage missing data ----------------------------------------
miss <- data.frame('miss_full' = paste0(colSums(is.na(data)),' (',round((colSums(is.na(data))/nrow(data))*100, 1), '%)'),
                 'miss_mri_13' = paste0(colSums(is.na(mri)),' (',round((colSums(is.na(mri))/nrow(mri))*100, 1), '%)'),
                 'miss_dti_13' = paste0(colSums(is.na(dti)),' (',round((colSums(is.na(dti))/nrow(dti))*100, 1), '%)'),
               'miss_mri_9&13' = paste0(colSums(is.na(mri_long)),' (',round((colSums(is.na(mri_long))/nrow(mri_long))*100, 1), '%)'),
               'miss_dti_9&13' = paste0(colSums(is.na(dti_long)),' (',round((colSums(is.na(dti_long))/nrow(dti_long))*100, 1), '%)'),
                    row.names = names(data))
# View(miss)
write.csv(miss, file.path(respath, 'MissPattern.csv'))

# Create correlation matrix (full sample) --------------------------------------
data <- data[ , -which(names(data) %in% selection)] # get rid of selection variables 

c <- round(cor(data[,-which(names(data) %in% c('sex','ethnicity','twin', # binary and categorical
                                               'm_educ_pregn','m_educ_3','m_educ_6'))], 
               use='complete.obs'), 2)
write.csv(c, file.path(respath,'CorrMatrix.csv'))

# IMPUTATION ===================================================================

# Use random forest for continuous data, logistic regression for binomial (sex) and polyreg for categorical(ethnicity, m_educ_6)
meth <- make.method(data, defaultMethod = c("rf", "logreg", "polyreg"))

# Random fortest imputation 
imp_rf <- mice(data, method = meth, m = 20, maxit = 40, ntree = 100) 

# IMPUTATION (old method) ======================================================
# # RF single dataset ------------------------------------------------------------
# library(missForest)
# imp <- missForest(data, maxiter = 10, ntree = 100, # takes about 5 min per iteration
#                  variablewise = T, decreasing = F, verbose = T)
# impd <-imp$ximp
# write.csv(impd, file.path(respath,'missForest_imp.csv'))
# # ------------------------------------------------------------------------------
# imp0 <- mice::mice(data, maxit = 0)
# 
# # Quick prediction -------------------------------------------------------------
# qp <- mice::quickpred(data, mincor=0.2, exclude ='IDC')
# # Do not impute the outcome variables or complete ones
# # qp[c(outcomes, subc_vol, tracts, 'tiv_13', 'sex', 'age_13', 'm_age'), ] <- 0
# 
# write.csv(qp, file.path(respath,'impQC','Qpred.csv'))
# 
# # Create an empty predictor matrix to fill with instructions -------------------
# pm <- matrix(0, ncol(data), ncol(data), dimnames=list(names(data),names(data)))

# pm[, c(exposures, outcomes, outcomes_9, tracts, subc_vol)]<- 1

# pm['imt_9', c(exposures, covariates[!covariates=='age_13'], 'age_9','gest_weight')]<- 1
# pm['dis_9', c(exposures, covariates[!covariates=='age_13'], 'sbp_6','gest_weight')] <- 1
# pm['sbp_9', c(exposures, covariates[!covariates=='age_13'], 'sbp_6','dbp_6','m_bmi_prepregn','gest_weight')] <- 1
# pm['dbp_9', c(exposures, covariates[!covariates=='age_13'], 'sbp_6','dbp_6','m_bmi_prepregn','income_6', 'gest_weight')] <- 1
# pm['height_9', c(exposures, covariates[!covariates%in%c('age_13','bmi_9_z','tiv_13')], auxiliary[!auxiliary%in%c('bmi_6_z','m_bmi_6')])] <- 1
# pm['ethnicity',c(exposures, covariates[!covariates%in%c('age_13','bmi_9_z')], auxiliary[!auxiliary%in%c('age_9','m_bmi_prepregn')])] <- 1
# pm['bmi_9_z',  c(exposures, covariates[!covariates%in%c('age_13','height_9','tiv_13')],auxiliary[!auxiliary%in%c('height_6','m_bmi_6')])] <- 1
# pm['m_educ_6', c(exposures, covariates[!covariates%in%c('height_9','tiv_13')], auxiliary[!auxiliary%in%c('age_9','m_bmi_prepregn')])] <- 1
# pm[auxiliary, auxiliary] <- 1
# 
# diag(pm) <- 0 # Make sure no circularity arises
# 
# write.csv(pm, file.path(respath,'impQC','PredM.csv'))
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
# # Imputation method ------------------------------------------------------------
# meth <- imp0$method
# 
# # Run imputation ---------------------------------------------------------------
# imputation <- mice::mice(data, m = 30, # nr of imputed datasets
#                          maxit = 60, # nr of iterations
#                          seed = 310896, # set a seed for the random number generation
#                          predictorMatrix = pm, 
#                          method = meth)

# POST-PROCESSING ==============================================================

# Add dichotomized ethnicity 
long.impdata <- complete(imp_rf, 'long', include = TRUE) %>%
  mutate(ethnicity_dich = if_else(ethnicity %in% c('dutch','european_descent'), 0, 1))
# Convert back to mids
imp_rf <- as.mids(long.impdata)

# Save
saveRDS(imp_rf, file.path(respath,'imputation_list_full.rds'))

postprocess <- function(inclusion, filename) {
  # Subset
  samp <- miceadds::subset_datlist(imp_rf, subset = imp_rf$data$IDC %in% inclusion)
  # Standardize
  mains <- c(exposures, outcomes, outcomes_9, subc_vol, tracts)
  sampz <- miceadds::scale_datlist(samp, orig_var = mains, trafo_var = paste0(mains, '_z'))
  sampimp <- miceadds::datlist2mids(sampz)
  # Save
  saveRDS(sampimp, file.path(respath, filename))
  return(sampimp)
}

mriset <- postprocess(include_mri, 'imputation_list_smri.rds')
dtiset <- postprocess(include_dti, 'imputation_list_dti.rds')

mriset_long <- postprocess(include_mri_long, 'imputation_list_smri_long.rds')
dtiset_long <- postprocess(include_dti_long, 'imputation_list_dti_long.rds')

# Save list of datasets to use as input for QDECR ------------------------------

dl = list()
for (n in 1:20) {
  d <- data.frame(complete(mriset, n))
  dl[[n]] <- d
}
saveRDS(dl, file.path(respath,'implist_qdecr.rds'))

# QUALITY CONTROL ==============================================================

# pdf(file.path(respath,'impQC','covergence-plots-full.pdf'))
# plot(imp_rf)
# dev.off()

# for (var in colnames(data)) { cat(paste('densityplot(imp_rf, ~', var, ')', '\n')) }
pdf(file.path(respath,'impQC', 'imp-vs-obs-full_RF100.pdf'))
densityplot(imp_rf, ~ imt_9 ) 
densityplot(imp_rf, ~ dis_9 ) 
densityplot(imp_rf, ~ sbp_9 ) 
densityplot(imp_rf, ~ dbp_9 ) 
densityplot(imp_rf, ~ tbv_13 ) 
densityplot(imp_rf, ~ gmv_13 ) 
densityplot(imp_rf, ~ mfa_13 ) 
densityplot(imp_rf, ~ mmd_13 ) 
densityplot(imp_rf, ~ sex ) 
densityplot(imp_rf, ~ age_13 ) 
densityplot(imp_rf, ~ height_9 ) 
densityplot(imp_rf, ~ ethnicity ) 
densityplot(imp_rf, ~ bmi_9_z ) 
densityplot(imp_rf, ~ m_educ_6 ) 
densityplot(imp_rf, ~ m_age ) 
densityplot(imp_rf, ~ tiv_13 ) 
densityplot(imp_rf, ~ accumbens ) 
densityplot(imp_rf, ~ amygdala ) 
densityplot(imp_rf, ~ caudate ) 
densityplot(imp_rf, ~ hippocampus ) 
densityplot(imp_rf, ~ pallidum ) 
densityplot(imp_rf, ~ putamen ) 
densityplot(imp_rf, ~ thalamus ) 
densityplot(imp_rf, ~ cgc_FA ) 
densityplot(imp_rf, ~ cgc_MD ) 
densityplot(imp_rf, ~ cst_FA ) 
densityplot(imp_rf, ~ cst_MD ) 
densityplot(imp_rf, ~ unc_FA ) 
densityplot(imp_rf, ~ unc_MD ) 
densityplot(imp_rf, ~ ilf_FA ) 
densityplot(imp_rf, ~ ilf_MD ) 
densityplot(imp_rf, ~ slf_FA ) 
densityplot(imp_rf, ~ slf_MD ) 
densityplot(imp_rf, ~ fma_FA ) 
densityplot(imp_rf, ~ fma_MD ) 
densityplot(imp_rf, ~ fmi_FA ) 
densityplot(imp_rf, ~ fmi_MD ) 
densityplot(imp_rf, ~ tbv_9 ) 
densityplot(imp_rf, ~ gmv_9 ) 
densityplot(imp_rf, ~ mfa_9 ) 
densityplot(imp_rf, ~ mmd_9 ) 
densityplot(imp_rf, ~ age_mri_9 ) 
densityplot(imp_rf, ~ tiv_9 ) 
densityplot(imp_rf, ~ cortex_9 ) 
densityplot(imp_rf, ~ subcort_9 ) 
densityplot(imp_rf, ~ wmv_9 ) 
densityplot(imp_rf, ~ cortex_13 ) 
densityplot(imp_rf, ~ subcort_13 ) 
densityplot(imp_rf, ~ wmv_13 ) 
densityplot(imp_rf, ~ gest_age_birth ) 
densityplot(imp_rf, ~ gest_weight ) 
densityplot(imp_rf, ~ parity ) 
densityplot(imp_rf, ~ m_bmi_prepregn ) 
densityplot(imp_rf, ~ m_bmi_6 ) 
densityplot(imp_rf, ~ income_6 ) 
densityplot(imp_rf, ~ sbp_6 ) 
densityplot(imp_rf, ~ dbp_6 ) 
densityplot(imp_rf, ~ height_6 ) 
densityplot(imp_rf, ~ bmi_6_z ) 
densityplot(imp_rf, ~ m_educ_pregn ) 
densityplot(imp_rf, ~ m_educ_3 ) 
densityplot(imp_rf, ~ age_9 ) 
dev.off()

pdf(file.path(respath,'impQC', 'imp-vs-obs-smri_RF100.pdf'))
densityplot(mriset, ~ imt_9 ) 
densityplot(mriset, ~ dis_9 ) 
densityplot(mriset, ~ sbp_9 ) 
densityplot(mriset, ~ dbp_9 ) 
densityplot(mriset, ~ mfa_13 ) 
densityplot(mriset, ~ mmd_13 ) 
densityplot(mriset, ~ height_9 ) 
densityplot(mriset, ~ ethnicity ) 
densityplot(mriset, ~ bmi_9_z ) 
densityplot(mriset, ~ m_educ_6 ) 
densityplot(mriset, ~ cgc_FA ) 
densityplot(mriset, ~ cgc_MD ) 
densityplot(mriset, ~ cst_FA ) 
densityplot(mriset, ~ cst_MD ) 
densityplot(mriset, ~ unc_FA ) 
densityplot(mriset, ~ unc_MD ) 
densityplot(mriset, ~ ilf_FA ) 
densityplot(mriset, ~ ilf_MD ) 
densityplot(mriset, ~ slf_FA ) 
densityplot(mriset, ~ slf_MD ) 
densityplot(mriset, ~ fma_FA ) 
densityplot(mriset, ~ fma_MD ) 
densityplot(mriset, ~ fmi_FA ) 
densityplot(mriset, ~ fmi_MD ) 
densityplot(mriset, ~ tbv_9 ) 
densityplot(mriset, ~ gmv_9 ) 
densityplot(mriset, ~ mfa_9 ) 
densityplot(mriset, ~ mmd_9 ) 
densityplot(mriset, ~ age_mri_9 ) 
densityplot(mriset, ~ tiv_9 ) 
densityplot(mriset, ~ cortex_9 ) 
densityplot(mriset, ~ subcort_9 ) 
densityplot(mriset, ~ wmv_9 ) 
densityplot(mriset, ~ gest_age_birth ) 
densityplot(mriset, ~ gest_weight ) 
densityplot(mriset, ~ parity ) 
densityplot(mriset, ~ m_bmi_prepregn ) 
densityplot(mriset, ~ m_bmi_6 ) 
densityplot(mriset, ~ income_6 ) 
densityplot(mriset, ~ sbp_6 ) 
densityplot(mriset, ~ dbp_6 ) 
densityplot(mriset, ~ height_6 ) 
densityplot(mriset, ~ bmi_6_z ) 
densityplot(mriset, ~ m_educ_pregn ) 
densityplot(mriset, ~ m_educ_3 ) 
densityplot(mriset, ~ age_9 ) 
dev.off()

pdf(file.path(respath,'impQC', 'imp-vs-obs-dti_RF100.pdf'))
densityplot(dtiset, ~ imt_9 ) 
densityplot(dtiset, ~ dis_9 ) 
densityplot(dtiset, ~ sbp_9 ) 
densityplot(dtiset, ~ dbp_9 ) 
densityplot(dtiset, ~ tbv_13 ) 
densityplot(dtiset, ~ gmv_13 ) 
densityplot(dtiset, ~ height_9 ) 
densityplot(dtiset, ~ ethnicity ) 
densityplot(dtiset, ~ bmi_9_z ) 
densityplot(dtiset, ~ m_educ_6 ) 
densityplot(dtiset, ~ tiv_13 ) 
densityplot(dtiset, ~ accumbens ) 
densityplot(dtiset, ~ amygdala ) 
densityplot(dtiset, ~ caudate ) 
densityplot(dtiset, ~ hippocampus ) 
densityplot(dtiset, ~ pallidum ) 
densityplot(dtiset, ~ putamen ) 
densityplot(dtiset, ~ thalamus ) 
densityplot(dtiset, ~ tbv_9 ) 
densityplot(dtiset, ~ gmv_9 ) 
densityplot(dtiset, ~ mfa_9 ) 
densityplot(dtiset, ~ mmd_9 ) 
densityplot(dtiset, ~ age_mri_9 ) 
densityplot(dtiset, ~ tiv_9 ) 
densityplot(dtiset, ~ cortex_9 ) 
densityplot(dtiset, ~ subcort_9 ) 
densityplot(dtiset, ~ wmv_9 ) 
densityplot(dtiset, ~ cortex_13 ) 
densityplot(dtiset, ~ subcort_13 ) 
densityplot(dtiset, ~ wmv_13 ) 
densityplot(dtiset, ~ gest_age_birth ) 
densityplot(dtiset, ~ gest_weight ) 
densityplot(dtiset, ~ parity ) 
densityplot(dtiset, ~ m_bmi_prepregn ) 
densityplot(dtiset, ~ m_bmi_6 ) 
densityplot(dtiset, ~ income_6 ) 
densityplot(dtiset, ~ sbp_6 ) 
densityplot(dtiset, ~ dbp_6 ) 
densityplot(dtiset, ~ height_6 ) 
densityplot(dtiset, ~ bmi_6_z ) 
densityplot(dtiset, ~ m_educ_pregn ) 
densityplot(dtiset, ~ m_educ_3 ) 
densityplot(dtiset, ~ age_9 ) 
dev.off()

# for (var in colnames(data)) { cat(paste('stripplot(imp_rf,', var, '~ .imp)', '\n')) }
pdf(file.path(respath,'impQC', 'stripplot_RF100.pdf'))
stripplot(imp_rf, imt_9 ~ .imp) 
stripplot(imp_rf, dis_9 ~ .imp) 
stripplot(imp_rf, sbp_9 ~ .imp) 
stripplot(imp_rf, dbp_9 ~ .imp) 
stripplot(imp_rf, tbv_13 ~ .imp) 
stripplot(imp_rf, gmv_13 ~ .imp) 
stripplot(imp_rf, mfa_13 ~ .imp) 
stripplot(imp_rf, mmd_13 ~ .imp) 
stripplot(imp_rf, sex ~ .imp) 
stripplot(imp_rf, age_13 ~ .imp) 
stripplot(imp_rf, height_9 ~ .imp) 
stripplot(imp_rf, ethnicity ~ .imp) 
stripplot(imp_rf, bmi_9_z ~ .imp) 
stripplot(imp_rf, m_educ_6 ~ .imp) 
stripplot(imp_rf, m_age ~ .imp) 
stripplot(imp_rf, tiv_13 ~ .imp) 
stripplot(imp_rf, accumbens ~ .imp) 
stripplot(imp_rf, amygdala ~ .imp) 
stripplot(imp_rf, caudate ~ .imp) 
stripplot(imp_rf, hippocampus ~ .imp) 
stripplot(imp_rf, pallidum ~ .imp) 
stripplot(imp_rf, putamen ~ .imp) 
stripplot(imp_rf, thalamus ~ .imp) 
stripplot(imp_rf, cgc_FA ~ .imp) 
stripplot(imp_rf, cgc_MD ~ .imp) 
stripplot(imp_rf, cst_FA ~ .imp) 
stripplot(imp_rf, cst_MD ~ .imp) 
stripplot(imp_rf, unc_FA ~ .imp) 
stripplot(imp_rf, unc_MD ~ .imp) 
stripplot(imp_rf, ilf_FA ~ .imp) 
stripplot(imp_rf, ilf_MD ~ .imp) 
stripplot(imp_rf, slf_FA ~ .imp) 
stripplot(imp_rf, slf_MD ~ .imp) 
stripplot(imp_rf, fma_FA ~ .imp) 
stripplot(imp_rf, fma_MD ~ .imp) 
stripplot(imp_rf, fmi_FA ~ .imp) 
stripplot(imp_rf, fmi_MD ~ .imp) 
stripplot(imp_rf, tbv_9 ~ .imp) 
stripplot(imp_rf, gmv_9 ~ .imp) 
stripplot(imp_rf, mfa_9 ~ .imp) 
stripplot(imp_rf, mmd_9 ~ .imp) 
stripplot(imp_rf, age_mri_9 ~ .imp) 
stripplot(imp_rf, tiv_9 ~ .imp) 
stripplot(imp_rf, cortex_9 ~ .imp) 
stripplot(imp_rf, subcort_9 ~ .imp) 
stripplot(imp_rf, wmv_9 ~ .imp) 
stripplot(imp_rf, cortex_13 ~ .imp) 
stripplot(imp_rf, subcort_13 ~ .imp) 
stripplot(imp_rf, wmv_13 ~ .imp) 
stripplot(imp_rf, gest_age_birth ~ .imp) 
stripplot(imp_rf, gest_weight ~ .imp) 
stripplot(imp_rf, parity ~ .imp) 
stripplot(imp_rf, m_bmi_prepregn ~ .imp) 
stripplot(imp_rf, m_bmi_6 ~ .imp) 
stripplot(imp_rf, income_6 ~ .imp) 
stripplot(imp_rf, sbp_6 ~ .imp) 
stripplot(imp_rf, dbp_6 ~ .imp) 
stripplot(imp_rf, height_6 ~ .imp) 
stripplot(imp_rf, bmi_6_z ~ .imp) 
stripplot(imp_rf, m_educ_pregn ~ .imp) 
stripplot(imp_rf, m_educ_3 ~ .imp) 
stripplot(imp_rf, age_9 ~ .imp) 
dev.off()
# ==============================================================================
# END