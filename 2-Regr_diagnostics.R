# ==============================================================================
# ======================= 2. Regression diagnostics ============================
# ==============================================================================

# Required packages
invisible(lapply(c('ggplot2','gridExtra','splines','car','gvlma','mice','miceadds'),
                 require, character.only = T));

# Load mids object
if (exists("smri") == F) { 
  # Define path
  genrpath <- dirname(file.choose()) # project folder
  # Load imputed datasets
  mriset <- readRDS(file.path(genrpath,'results','imputation_list_smri.rds'))
  dtiset <- readRDS(file.path(genrpath,'results','imputation_list_dti.rds'))
  # Transform into mids
  imp_smri <- miceadds::datlist2mids(mriset)
  imp_dti <- miceadds::datlist2mids(dtiset)
  # extract original sets (with missing data)
  smri <- mice::complete(imp_smri, action = 0)
  dti  <- mice::complete(imp_dti, action = 0)
  # clean up
  rm(mriset,dtiset,imp_smri,imp_dti)
  # try excluding outlier
  smri <- smri[-903,]
}


# Define main exposures and outcomes
outcomes = paste0(c('tbv','gmv','mfa','mmd'),'_13_z')
exposures = paste0(c('imt','dis','sbp','dbp'),'_9_z')

# Define covariates
covs1 <- ' + sex + age_13 + height_9 + ethnicity_dich'
covs2 <- ' + bmi_9_z + m_educ_6 + m_age'

# Quickly fit the lm() with relevant variables 
form <- function(outc, exp, data, adj='full') {
  if (adj=='base') { covs2 <- '' }
  fit <- lm(as.formula(paste0(outc,' ~ ',exp,covs1,covs2)), data=data)
  # print(summary(fit))
  return(fit)
}

# Save output to txt log file 
sink(paste0('Diagnostic_output_',Sys.Date(),'.txt'))

# ==============================================================================
# Remember that distinct problems can interact: if the errors have a skewed distribution, 
# apparent outliers may be produced in the direction of the skew. Transforming the 
# response to make the errors less skewed can solve this problem. Similarly, properly 
# modeling a nonlinear relationship may bring apparently outlying observations in line 
# with the rest of the data.
# ==============================================================================

pdf(paste0('Diagnostic_plots_',Sys.Date(),'.pdf'))
for (e in exposures) {
  for (o in outcomes) {
    if (any(sapply(c('gmv','tbv'), grepl, o))) { data <- smri 
    } else if (any(sapply(c('fa','md'), grepl, o, ignore.case=T))) { data <- dti } 
    
    fit <- form(o, e, data)
    title = paste(toupper(substr(e,1,3)),'~',toupper(substr(o,1,3)))
    cat('\n====',title,'===================================================\n')
    print(gvlma::gvlma(fit))
    cat('\n --- General fit -------------------------------\n')
    # Systematic features in residual plots, e.g. curvature or apparent non-constant variance
    # require to modify the structure of the model to match the data more closely. 
    car::residualPlots(fit, main = title) 
    # Useful for revealing problems but less for determining the nature of the problem. 
    # These should be null plots, with no systematic features. One non-null plot (curved general trend) 
    # is enough to suggest that the specified model does not match the data, generally implying a failure 
    # of one or more assumptions.
    # A lack-of-fit test is also computed for each numeric predictor (i.e. a t-test for (regressor)^2 added to the model. 
    # Significant p-values indicate lack-of-fit, confirming the nonlinear pattern visible in the graph. 
    # For the plot of residuals vs fitted values, the Tukey’s test for non-additivity is obtained by adding the squares 
    # of the fitted values to the model and refitting. The test confirms the impression of curvature (model is not adequate).
    car::marginalModelPlots(fit, main = title) # sd=T to check variance assumptions
    # A lowess smooth is fit to the data (in blue) and to the fitted values (red dashed)
    # and the two curves should match on each of the plots. Any mismatches are evidence of bad fit.
    cat('\n --- Outliers & Influential observations ------\n')
    # Generic qqPlot: Studentized residuals against corresponding quantiles of t(n-k—2),
    # with 95% pointwise confidence envelope. Also returns observations with the largest absolute Studentized residuals.
    # Note qqplots don't work well with the form function
    car::qqPlot(lm(as.formula(paste(o,'~',e,covs1)), data=data), main=title, ylab='Studentized Redisuals')
    # Bonferroni p-value for most extreme obs (i.e. highest Studentized residuals)
    print(car::outlierTest(fit))
    # index plots of Cook’s distances, Studentized residuals, corresponding 
    # Bonferroni p values for outlier testing, and the hat-values (measure of leverage).
    car::influenceIndexPlot(fit, main = title)
    car::influencePlot(fit, id.method="identify", main=title, sub="Circle size proportial to Cook's Distance" )
    cat('\n --- Collinearity -----------------------------\n')
    # print(vif(fit)) # variance inflation factors
    print(sqrt(vif(fit)) > 2) # problem?
  }
}
dev.off()

# ==============================================================================
cat('\n=============================================================
     \n===================== NON LINEARITY =========================
     \n=============================================================\n')
test_nlin <- function(method){
  pdf(paste0('nlin_graph_',method,'_',Sys.Date(),'.pdf'), width=25, height=5)
  for (e in exposures) {
    for (o in outcomes) {
      # Define datasets
      if (any(sapply(c('gmv','tbv'), grepl, o))) { data <- smri 
      } else if (any(sapply(c('fa','md'), grepl, o, ignore.case=T))) { data <- dti } 
      # fit standard (linear) model
      fit <- form(o, e, data)
      
      p <- list() # create empty list to store p-values
      # plot linear relation 
      p[[1]] <- ggplot(data, aes_string(e,o) ) + geom_point() +
        ggtitle(paste(toupper(substr(e,1,3)),'~',toupper(substr(o,1,3)), '--> linear model')) +
        labs(y = toupper(substr(o,1,3)), x = toupper(substr(e,1,3))) +
        stat_smooth(method = lm, formula = y ~ x)
      
      # fit polynomials 2-5
      for (test in c(2:5)) {
        if (method =='poly') {
          fit2 <- form(o, paste0('poly(',e,', ',test,', raw=T)'), data)
          line = formula(paste0('y ~ poly(x,',test,', raw = T)'))
        } else if (method=='expont') {
          fit2 <- form(o, paste0(e,' + I(',e,'^',test,')'), data)
          line = formula(paste0('y ~ I(x^',test,')'))
        } else if (method=='splines') {
          fit2 <- form(o, paste0('splines::ns(',e,', ',test,')'), data)
          line = formula(paste0('y ~ ns(x,',test,')'))
        }
        # LRT significance test (compared to linear)
        a <- anova(fit, fit2); sign <- ''
        if (!is.na(a[2,'Pr(>F)']) & a[2,'Pr(>F)'] < 0.05) {
          cat(test, '-', e, '~', o, '--------------------------------------\n')
          sign <- '*'
          print(summary(fit2)); print(a)
        }
        pval <- paste(round(a[2,'Pr(>F)'], 3), sign)
        
        p[[test]] <- ggplot(data, aes_string(e,o) ) + geom_point() +
          ggtitle(paste(toupper(substr(e,1,3)),'~',toupper(substr(o,1,3)), '-->', method,'=',test,', p =', pval)) +
          labs(y = toupper(substr(o,1,3)), x = toupper(substr(e,1,3))) +
          stat_smooth(method = lm, formula = as.formula(line))
      }  
      do.call(grid.arrange, c(p, ncol=5))
    }
  }
  dev.off()
}

test_nlin('poly')
test_nlin('expont')
test_nlin('splines')

sink()
# ==============================================================================

# Try transformations of the outcome
# intern.sqrt <- update( intern, sqrt(.) ~ .)

# Evaluate the interaction between stress periods for each model
# intern.int  <- update( intern,  . ~ . + prenatal_stress_z:postnatal_stress_z)