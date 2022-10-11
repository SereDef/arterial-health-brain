# ==============================================================================
# =========================== 1. Descriptives ==================================
# ==============================================================================

# Required packages
invisible(lapply(c('mice','miceadds','openxlsx'), require, character.only = T));

# Define path
genrpath <- dirname(file.choose()) # project folder

# FULL SAMPLE (NO SELCTION, NO IMPUTATION ) ====================================
if (exists('orig_data')==F) { orig_data <- read.csv(file.path(genrpath,'DATA','Data.csv'), stringsAsFactors=T)[,-1]}

summdf <- function(dataframe, bysex=F) {
  # take summary object, clean the strings and note them as row.names, return a data.frame
  # with columns:
  nms <- c('Var','Min','1stQ','Median','Mean','3rdQ','Max','NAs','SD')
  # For imputed datasets stratified by sex I input a summary object already and there are no NAs
  if (bysex==T) { input <- dataframe; n_par <- 6; nms <- nms[nms!='NAs']
  } else { input <- summary(dataframe,digits=4); n_par <- 7 }
  # Flip the set to have variables as rows and metrics as columns 
  summ <- as.data.frame(input)
  summ[,'Var1'] <- rep(c(1:n_par), nrow(summ)/n_par)
  summ <- reshape(summ, idvar='Var2', timevar='Var1', direction='wide')
  # Add SDs to the other metrics
  sds <- round(apply(dataframe, 2, sd, na.rm = T), 4)
  summ <- cbind(summ,sds)
  # Adjust column names
  names(summ) <- nms
  
  return(summ)
}

# Entire sample summary (before imputation!)
summ_orig <- summdf(orig_data)

# FULL AND SELECTED SAMPLE (AFTER IMPUTATION ) =================================

describe <- function(imp_file, cat_vars = c('sex','ethnicity','m_educ_6','twin','m_educ_3','m_educ_pregn')) {
  # Load imputed list object
  imput <- readRDS(file.path(genrpath,'results',imp_file))
  
  # Extract the original set (with NAs)
  sample <- complete(imput, 0) 
  # Summary of sample before imputation
  summ_sample <- summdf(sample)
  
  # Stack imputed datasets in long format, excluding the original data
  impdat <- complete(imput, action="long", include = F)
  
  pool_descriptives <- function(implist, column_names, categorical=T) {
    summ <- with(implist, by(implist, .imp, function(x) summary(x[, -c(1, 2)],digits=4))) 
    if (categorical==F) {
      # Pool summary 
      num_pool <- lapply(summ, function(m) matrix(as.numeric(sapply(strsplit(m, ":"), "[[", 2)), nrow = dim(m)[1], ncol=dim(m)[2]))
      pool_mean <- Reduce("+",num_pool)/length(num_pool)
      # Pool SDs
      sds <- with(implist, by(implist, .imp, function(x) round(apply(x[, -c(1, 2)], 2, sd, na.rm = T), 4)))
      pool_sds <- Reduce("+",sds)/length(sds)
      # Bind SDs to other metrics
      summ_df <- data.frame(rbind(pool_mean,pool_sds))
      # Define column and row names
      colnames(summ_df) <- colnames(implist[-c(1,2)])
      rownames(summ_df) <- c('Min','1stQ','Median','Mean','3rdQ','Max','SD')
    } else { 
      pool_mean <- Reduce("+",summ)/length(summ) 
      summ_df <- data.frame(pool_mean)
      colnames(summ_df) <- 'counts' # colnames(implist[-c(1,2)])
      rownames(summ_df) <- names(pool_mean)
    }
    return(summ_df)
  }
  # Continuous data
  cnt <- impdat[, -c(which(colnames(impdat) %in% cat_vars))]
  cnt_summ <- pool_descriptives(cnt, categorical = F)
  
  # Categorical / binary
  cat_summ <- NA
  for (v in cat_vars) {
    sel <- impdat[, c('.imp','.id', v)]
    v_summ <- pool_descriptives(sel)
    cat_summ <- rbind(cat_summ, v, v_summ, NA)
  }
  
  # Descriptives by sex
  bysex <- by(impdat, impdat$sex, summary, digits = 4)
  boy_summ <- summdf(bysex[[1]], bysex = T)
  grl_summ <- summdf(bysex[[2]], bysex = T)
  
  # Correlation matrix in the imputed set
  cors_imp <- miceadds::micombine.cor(mi.res = impdat, 
                                      variables = colnames(impdat)[!colnames(impdat) %in% 
                                                                     c('.imp','.id','IDC',cat_vars)]) 
  
  # Export the outputs of summary statistics into an xlsx file with one model per sheet
  stats <- list('s_orig' = summ_sample, 's_imp_cnt' = cnt_summ, 's_imp_cat' = cat_summ, 
                's_imp_boy' = boy_summ, 's_imp_grl' = grl_summ, 'cor_imp' = cors_imp)
  
  # Adjust names
  s <- strsplit(strsplit(imp_file,'list')[[1]][2], '\\.')[[1]][1] # get sample name 
  names(stats) <- paste0(names(stats), s)
  
  return(stats)
}

full <- describe('imputation_list_allimp.rds')
smri <- describe('imputation_list_smri.rds')
dti  <- describe('imputation_list_dti.rds')

# Stack them together and export to xlsx file ==================================

stats <- c(list('s_full'=summ_orig), full, smri, dti)

# Export summary statistics into an xlsx file with one summary per sheet
openxlsx::write.xlsx(stats, file = file.path(genrpath,'results',paste0(Sys.Date(),"_Descriptives.xlsx")), 
                     rowNames = T, overwrite = T)

# ==============================================================================