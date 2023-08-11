# ==============================================================================
# =========================== 1. Descriptives ==================================
# ==============================================================================

# Required packages
invisible(lapply(c('mice','miceadds','openxlsx'), require, character.only = T));

# date <- format(Sys.Date(), "%d%m%y")
date <- '240523'

# Define paths 
genrpath <- dirname(file.choose()) # <== choose project folder
datapath <- file.path(genrpath,'DATA') # data folder
respath  <- file.path(genrpath, paste0('results_',date)) # results folder

# FULL SAMPLE (NO SELCTION, NO IMPUTATION ) ====================================
if (exists('orig_data')==F) { orig_data <- read.csv(file.path(datapath,'Data.csv'), stringsAsFactors=T)[,-1]}

summdf <- function(dataframe) {
  # take summary object, clean the strings and note them as row.names, return a data.frame
  # with columns:
  nms <- c('Var','Min','1stQ','Median','Mean','3rdQ','Max','NAs','SD')
  # For imputed datasets stratified by sex I input a summary object already and there are no NAs
  input <- summary(dataframe,digits=4); 
  n_par <- 7
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

describe <- function(sample_name) {
  # Construct file name
  imp_file <- paste0('imp_',sample_name,'_',date,'.rds')
  # Load imputed list object
  imput <- readRDS(file.path(respath, imp_file))
  
  # determine categorical and continuous vars 
  lvl_length <- lapply(imput$data, function(var) length(levels(as.factor(var)))) 
  # Cutoff 15 levels: consider it categorical
  cat_vars <- names(which(lvl_length < 20))
  
  # Extract the original set (with NAs)
  sample <- complete(imput, 0) 
  # Summary of sample before imputation
  summ_sample <- summdf(sample)
  
  # Stack imputed datasets in long format, excluding the original data
  impdat <- mice::complete(imput, action="long", include = F)
  
  # Set to factors or numeric when appropiate
  impdat[,cat_vars] <- lapply(impdat[,cat_vars] , as.factor)
  impdat = impdat[,-grep('IDC',names(impdat))] # remove ID from descriptives 
  impdat[, !names(impdat)%in%cat_vars] <- lapply(impdat[,!names(impdat)%in%cat_vars] , as.numeric)
  
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
    v_summ <- cbind(row.names(v_summ), v_summ)
    v_summ$percent <- (v_summ$count / nrow(imput$data))*100
    cat_summ <- rbind(cat_summ, v, v_summ, NA)
  }
  
  # Descriptives by sex
  bysex <- by(impdat, impdat$sex, summdf)
  
  # Correlation matrix in the imputed set
  cors_imp <- miceadds::micombine.cor(mi.res = impdat, 
                                      variables = colnames(impdat)[!colnames(impdat) %in% 
                                                c('.imp','.id','IDC',cat_vars)]) 
  
 # Export the outputs of summary statistics into an xlsx file with one model per sheet
  stats <- list('pre_imp' = summ_sample, 'imp_cnt' = cnt_summ, 'imp_cat' = cat_summ, 
                'imp_boy' = bysex[['boy']], 'imp_grl' = bysex[['girl']], 
                'imp_cor' = cors_imp)
  
  # Adjust names
  names(stats) <- paste0(names(stats), sample_name)
  
  return(stats)
}

full <- describe('full')
bas9 <- describe('base9')
bas5 <- describe('base5')
smri <- describe('smri')
dti  <- describe('dti')

# Stack them together and export to xlsx file ==================================

# list('s_orig'=summ_orig)
stats <- c(full, bas9, bas5, smri, dti)

# Export summary statistics into an xlsx file with one summary per sheet
openxlsx::write.xlsx(stats, 
                     file = file.path(respath,paste0('1-Descriptives_',date,'.xlsx')), 
                     rowNames = T, overwrite = T)

# ==============================================================================