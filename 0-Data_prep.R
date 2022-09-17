# ==============================================================================
# ======================= 0. Data prep & imputation ============================
# ==============================================================================

# Required packages
invisible(lapply(c('foreign','dplyr','mice','miceadds'), require, character.only = T));

# Define paths 
datapath <- dirname(file.choose()) # data folder
respath  <- dirname(file.choose()) # results folder
# list.files(datapath) # check data files

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

brain <- read.csv(file.path(datapath, 'Brain.csv')) # TODO: add brain merging script
names(brain)[2] <- 'IDC' # fix lower case idc

# Merge
dataset <- merge(gen, m_bmi, by='MOTHER')
dataset <- Reduce(function(x,y) merge(x = x, y = y, by = 'IDC', all.x = T),
                      list(dataset, bp5, bp9, dis, imt, bmi5, bmi9, bmi13, brain) ) 
rm(gen, bp5, bp9, dis, imt, bmi5, bmi9, bmi13, brain)

mean_tract <- function(tract, metric='FA'){
  t <- rowMeans(dataset[,c(paste0(tract,'_l_dti_dipy_wls_wavg_', metric, '_f13'), 
                  paste0(tract,'_r_dti_dipy_wls_wavg_', metric, '_f13'))])
  return(t)
}
mean_region <- function(region){
  r <- rowMeans(dataset[,c(paste0('Right_',region, '_vol_f13'), 
                           paste0('Left_', region, '_vol_f13'))])
  return(r)
}

# Select and transform needed variables 
data <- data.frame('IDC' = dataset$IDC, 
                 'imt_9' = dataset$IMT_mean_combined,
                 'dis_9' = dataset$DIS_mean_combined,
                 'sbp_9' = dataset$MEANSBP_ex1_child9, 
                 'dbp_9' = dataset$MEANDBP_ex1_child9,
                'tbv_13' = dataset$genr_tbv_f13,
                'gmv_13' = as.numeric(dataset$TotalGrayVol_f13),
                'mfa_13' = dataset$mean_FA_genr_f13, # or all
                'mmd_13' = dataset$mean_MD_genr_f13,
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
              'm_educ_6' = as.numeric(dataset$EDUCM5), # maternal education @5
                 'm_age' = dataset$AGE_M_v2, # maternal age at intake
                'tiv_13' = as.numeric(dataset$eTIV_f13), # total intracranial volume
               # ------ SECONDARY ------- #
             'accumbens' = mean_region('Accumbens_area'),
              'amygdala' = mean_region('Amygdala'),
               'caudate' = mean_region('Caudate'),
           'hippocampus' = mean_region('Hippocampus'),
              'pallidum' = mean_region('Pallidum'),
               'putamen' = mean_region('Putamen'),
              'thalamus' = mean_region('Thalamus_Proper'),
                'cgc_FA' = mean_tract('cgc', 'FA'), # Cincgulate gyrus 
                'cgc_MD' = mean_tract('cgc', 'MD'), 
                'cst_FA' = mean_tract('cst', 'FA'), # Cortico-spinal tract
                'cst_MD' = mean_tract('cst', 'MD'), 
                'unc_FA' = mean_tract('unc', 'FA'), # uncinate fasciculus
                'unc_MD' = mean_tract('unc', 'MD'), 
                'ilf_FA' = mean_tract('ilf', 'FA'), # inferior longitudinal fasciculus
                'ilf_MD' = mean_tract('ilf', 'MD'), 
                'slf_FA' = mean_tract('slf', 'FA'), # superior longitudinal fasciculus
                'slf_MD' = mean_tract('slf', 'MD'), 
                'fma_FA' = dataset$fma_dti_dipy_wls_wavg_FA_f13, # major forceps 
                'fma_MD' = dataset$fma_dti_dipy_wls_wavg_MD_f13, # major forceps 
                'fmi_FA' = dataset$fmi_dti_dipy_wls_wavg_FA_f13, # minor forceps 
                'fmi_MD' = dataset$fmi_dti_dipy_wls_wavg_MD_f13, # ninor forceps 
               # ------ IMPUTATION ------- #
      # 'gest_age_birth' = dataset$GESTBIR,   # gestational age at birth (used for imputation)
           'gest_weight' = dataset$WEIGHT,    # gestational weight (used for imputation)
                'parity' = dataset$PARITY,    # parity (used for imputation)
        'm_bmi_prepregn' = dataset$BMI_0, # self-reported Maternal BMI
               'm_bmi_6' = dataset$BMIMotherF5,
              'income_6' = as.numeric(dataset$INCOME5), # income @5
                 'sbp_6' = dataset$MeanSBP_ex1_5child, 
                 'dbp_6' = dataset$MeanDBP_ex1_5child,
              'height_6' = dataset$length5child,
               'bmi_6_z' = dataset$sdsbmiforage,
         # 'm_educ_pregn' = as.numeric(dataset$EDUCM),  # maternal education pregnancy
            # 'm_educ_3' = as.numeric(dataset$EDUCM3), # maternal education @3
                 'age_9' = dataset$agechild9_visit1, # just in case 
                # ------ SELECTION ------- #
           'mri_consent' = as.factor(dataset$mri_consent_f13), # consent status
            'mri_braces' = as.factor(dataset$has_braces_mri_f13), # has braces
          'mri_inc_find' = as.factor(dataset$exclude_incidental_f13), # incidental findings
           't1_scantype' = as.factor(dataset$t1_has_nii_f13), # scan type
                 't1_qc' = as.factor(dataset$freesurfer_qc_f13), # QC
          'dti_scantype' = as.factor(dataset$dti_has_nii_f13), # scan type
                'dti_qc' = as.factor(dataset$dti_overall_qc_f13), # QC
                  'm_ID' = dataset$MOTHER, # mother id used to identify siblings (for exclusion)
                  'twin' = as.factor(dataset$TWIN))
# summary(data)
# Save histograms of all variables in the full sample
pdf(file.path(respath,'hist_full_sample.pdf'))
for (var in names(data)) { hist(as.numeric(data[,var]), col = 'blue', main = toupper(var), xlab = var) }
dev.off() 

# Save dataset
write.csv(data, file.path(datapath, 'Data.csv'))

# ==============================================================================
# Organize variable names if groups
exposures  <- c('imt_9', 'dis_9', 'sbp_9', 'dbp_9')
outcomes   <- c('tbv_13', 'gmv_13', 'mfa_13', 'mmd_13')
covariates <- c('sex', 'age_13', 'height_9', 'ethnicity', 'bmi_9_z', 'm_educ_6', 'm_age','tiv_13')
subc_vol   <- c('accumbens', 'amygdala', 'caudate', 'hippocampus', 'pallidum', 'putamen', 'thalamus')
tracts     <- c('cgc_FA', 'cgc_MD', 'cst_FA', 'cst_MD', 'unc_FA', 'unc_MD', 'ilf_FA', 'ilf_MD', 'slf_FA', 'slf_MD', 'fma_FA', 'fma_MD', 'fmi_FA', 'fmi_MD')       
auxiliary  <- c('age_9','gest_weight', 'parity', 'm_bmi_prepregn', 'm_bmi_6','income_6','sbp_6','dbp_6','height_6', 'bmi_6_z')
selection  <- c('mri_consent', 'mri_braces', 'mri_inc_find', 't1_scantype', 't1_qc', 'dti_scantype', 'dti_qc', 'm_ID', 'twin')

# Check selection and missing rate =============================================
flowchart <- function(fullset, outcome) {
  
  fc <- list(initial_sample = nrow(fullset))
  
  sel <- function(df, condition, value, method='equal') {
    if (method=='equal')         { step <- df[(!is.na(df[,condition])) & (df[,condition] == value),] 
    } else if (method=='not_in') { step <- df[!(df[,condition] %in% value),] 
    } else if (method=='is_na')  { step <- df[is.na(df[,condition]),] 
    } else if (method=='not_na') { step <- df[!is.na(df[,condition]),] }
    loss <- nrow(step) - as.numeric(fc[length(fc)])
    add <- setNames(list(loss, nrow(step)), c(condition, paste0('after_',condition)))
    fc <<- c(fc, add)
    return(step)
  }
  
  step1 <- sel(fullset, 'mri_consent', 'yes')
  step2 <- sel(step1,'mri_braces',  'no')
  step3 <- sel(step2,'mri_inc_find','include')
  step4 <- sel(step3,'twin','No')
  
  if (outcome == 'smri') {
    step5 <- sel(step4,'t1_scantype','yes')
    step6 <- sel(step5,'t1_qc','usable')
  } else if (outcome == 'dti') {
    step5 <- sel(step4,'dti_scantype','yes')
    step6 <- sel(step5,'dti_qc','unusable', method='is_na')
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
  
  final <- sel(step6,'IDC', select_sibling(step6, column_selection=c('imt_9','dis_9','sbp_9','dbp_9','tbv_13','gmv_13','mfa_13','mmd_13')), 
               method='not_in') 

  print(fc)
  return(final)
}

mri <- flowchart(data, outcome='smri')
dti <- flowchart(data, outcome='dti')
# TODO: check why some residual missing in DTI mean variables?
dti <- dti[!is.na(dti$mfa_13),]; nrow(dti)

# Save list of selected IDs for later selection
include_mri <- mri$IDC
include_dti <- dti$IDC

# Calculate the percentage missing data ----------------------------------------
miss <- data.frame('n_miss_mri' = colSums(is.na(mri)),
                'perc_miss_mri' = paste(round((colSums(is.na(mri))/nrow(mri))*100, 1), '%'),
                   'n_miss_dti' = colSums(is.na(dti)),
                'perc_miss_dti' = paste(round((colSums(is.na(dti))/nrow(dti))*100, 1), '%'))
# View(miss)
write.csv(miss, file.path(respath, 'MissPattern.csv'))

# Create correlation matrix (full sample) --------------------------------------

data <- data[ , -which(names(data) %in% selection)] # get rid of selection variables 

c <- round(cor(data[,-which(names(data) %in% c('sex','ethnicity'))], # binary and categorical
               use='complete.obs'), 2)
write.csv(c, file.path(respath,'CorrMatrix.csv'))

# IMPUTATION ===================================================================

imp0 <- mice::mice(data, maxit = 0)

# Quick prediction -------------------------------------------------------------
qp <- mice::quickpred(data, mincor=0.2, exclude ='IDC')
# Do not impute the outcome variables or complete ones 
qp[c(outcomes, subc_vol, tracts, 'tiv_13', 'sex', 'age_13', 'm_age'), ] <- 0 

write.csv(qp, file.path(respath,'impQC','Qpred.csv'))

# Create an empty predictor matrix to fill with instructions -------------------
pm <- matrix(0, ncol(data), ncol(data), dimnames=list(names(data),names(data)))

pm['imt_9', c(exposures, covariates[!covariates=='age_13'], 'age_9','gest_weight')]<- 1
pm['dis_9', c(exposures, covariates[!covariates=='age_13'], 'sbp_6','gest_weight')] <- 1
pm['sbp_9', c(exposures, covariates[!covariates=='age_13'], 'sbp_6','dbp_6','m_bmi_prepregn','gest_weight')] <- 1
pm['dbp_9', c(exposures, covariates[!covariates=='age_13'], 'sbp_6','dbp_6','m_bmi_prepregn','income_6', 'gest_weight')] <- 1
pm['height_9', c(exposures, covariates[!covariates%in%c('age_13','bmi_9_z','tiv_13')], auxiliary[!auxiliary%in%c('bmi_6_z','m_bmi_6')])] <- 1
pm['ethnicity',c(exposures, covariates[!covariates%in%c('age_13','bmi_9_z')], auxiliary[!auxiliary%in%c('age_9','m_bmi_prepregn')])] <- 1
pm['bmi_9_z',  c(exposures, covariates[!covariates%in%c('age_13','height_9','tiv_13')],auxiliary[!auxiliary%in%c('height_6','m_bmi_6')])] <- 1
pm['m_educ_6', c(exposures, covariates[!covariates%in%c('height_9','tiv_13')], auxiliary[!auxiliary%in%c('age_9','m_bmi_prepregn')])] <- 1
pm[auxiliary, auxiliary] <- 1

diag(pm) <- 0 # Make sure no circularity arises

write.csv(pm, file.path(respath,'impQC','PredM.csv'))

# Check differences between quickpred and our matrix ---------------------------

# check_pm <- function(var) {
#   vq <- names(which(qp[var,]==1)); # print(vq)
#   vs <- names(which(pm[var,]==1)); # print(vs)
#   
#   print(setdiff(vq, vs))     
#   print(setdiff(vs, vq))   
#   
#   vs = vs[!vs %in% c('sex','ethnicity')]
#   m = round(cor(data[,vs], use='pairwise.complete.obs'),2)
# 
#   maxes <- apply(matrix(m[row(m)!=col(m)], ncol=ncol(m)), 2, max)
#   mx <- data.frame(colnames(m), maxes)
#   # View(m); View(mx)
# }

# Imputation method ------------------------------------------------------------
meth <- imp0$method
meth['ethnicity'] <- 'pmm' # try if this works better
# TODO: maternal education improve prediction 

# Run imputation ---------------------------------------------------------------
imputation <- mice::mice(data, m = 30, # nr of imputed datasets
                         maxit = 60, # nr of iterations
                         seed = 310896, # set a seed for the random number generation
                         predictorMatrix = pm, 
                         method = meth)

# POST-PROCESSING ==============================================================

# Add dichotomized ethnicity 
long.impdata <- complete(imputation, 'long', include = TRUE) %>%
  mutate(ethnicity_dich = if_else(ethnicity %in% c('dutch','european_descent'), 0, 1))
# Convert back to mids
imputation <- as.mids(long.impdata)

# Subset
mriset <- miceadds::subset_datlist(imputation, subset = imputation$data$IDC %in% include_mri)
dtiset <- miceadds::subset_datlist(imputation, subset = imputation$data$IDC %in% include_dti)

# Standardize
mains <- c(exposures, outcomes, subc_vol, tracts)
mriset <- miceadds::scale_datlist(mriset, orig_var = mains, trafo_var = paste0(mains, '_z') )
dtiset <- miceadds::scale_datlist(dtiset, orig_var = mains, trafo_var = paste0(mains, '_z') )

# Save
saveRDS(imputation, file.path(respath,'imputation_list_full.rds'))
saveRDS(mriset, file.path(respath,'imputation_list_smri.rds'))
saveRDS(dtiset, file.path(respath,'imputation_list_dti.rds'))

# Save list of datasets to use as input for QDECR ------------------------------

# dl = list()
# for (n in 1:30) {
#   d <- data.frame(mriset[[n]])
#   dl[[n]] <- d
# }
# saveRDS(dl, 'implist_qdecr.rds')

# QUALITY CONTROL ==============================================================

pdf(file.path(respath,'impQC','covergence-plots-full.pdf'))
plot(imputation)
dev.off()

# for (var in colnames(data)) { cat(paste('densityplot(imputation, ~', var, ')', '\n')) }
pdf(file.path(respath,'impQC','imp-vs-obs-full.pdf'))
densityplot(imputation, ~ imt_9 ) 
densityplot(imputation, ~ dis_9 ) 
densityplot(imputation, ~ sbp_9 ) 
densityplot(imputation, ~ dbp_9 ) 
densityplot(imputation, ~ height_9 ) 
densityplot(imputation, ~ ethnicity ) 
densityplot(imputation, ~ bmi_9_z ) 
densityplot(imputation, ~ m_educ_6 ) 
densityplot(imputation, ~ gest_weight ) 
densityplot(imputation, ~ parity ) 
densityplot(imputation, ~ m_bmi_prepregn ) 
densityplot(imputation, ~ m_bmi_6 ) 
densityplot(imputation, ~ income_6 ) 
densityplot(imputation, ~ sbp_6 ) 
densityplot(imputation, ~ dbp_6 ) 
densityplot(imputation, ~ height_6 ) 
densityplot(imputation, ~ bmi_6_z ) 
densityplot(imputation, ~ age_9 ) 
dev.off()

# TODO: stripplot(imputation, ) ?
# TODO: QC in selected samples ?
# ==============================================================================