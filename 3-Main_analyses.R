# ==============================================================================
# =========================== 2. Main Analyses =================================
# ==============================================================================

# Required packages
invisible(lapply(c('mice','miceadds','splines','effects','openxlsx'), require, character.only = T));

# Load mids object
if (exists("genrpath") == F) { 
  # Define path
  genrpath <- dirname(file.choose()) # project folder
  # Load imputed datasets
  fulset <- readRDS(file.path(genrpath,'results','imp_full.rds'))
  mriset <- readRDS(file.path(genrpath,'results','imp_smri.rds'))
  dtiset <- readRDS(file.path(genrpath,'results','imp_dti.rds'))
  # Exclude one outlier
  # imp_smri <- miceadds::subset_datlist(imp_smri, 
  #             subset = imp_smri$data$IDC != imp_smri$data$IDC[903],  toclass = 'mids')
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
pool_mod <- function(outc, exp, sex='', fullsample=F) {
  # some basic values
  sex_cov <- '+ sex'; tiv_cov <- ''
  exp_name <- toupper(substr(outc,1,3))
  
  # Define the correct sample
  if (fullsample==T) { impset <- fulset 
  } else if (any(sapply(c('tbv','gmv',subc_vol), grepl, outc))) { impset <- mriset
  } else if (any(sapply(c('FA','MD'), grepl, outc, ignore.case=T))) { impset <- dtiset }

  # Stratify by sex
  if (sex!='') { 
    # this doesn't work in mriset males (Error in Ops.factor(left, right) : level sets of factors are different)
    # impset <- miceadds::subset_datlist(impset, subset = impset$data$sex == sex,  toclass = 'mids') 
    l <- complete(impset, 'long', include = T)
    s <- l[l$sex==sex,]
    impset <- as.mids(s)
    sex_cov <- '' 
  }

  # Add total intracranial volume to covariates for subcortical outcomes
  if (outc %in% c(subc_vol,trcts_FA,trcts_MD)) { tiv_cov <- '+ tiv_13' }

  # Specify the names of the supplementary outcomes 
  if (outc %in% subc_vol) { exp_name <- names(subc_vol)[which(subc_vol==outc)] 
  } else if (outc %in% trcts_FA) { exp_name <- names(trcts_FA)[which(trcts_FA==outc)] 
  } else if (outc %in% trcts_MD) { exp_name <- names(trcts_MD)[which(trcts_MD==outc)]  }
  
  # Define covariates
  covs1 <- paste(sex_cov,'+ age_13 + height_9',tiv_cov)
  covs2 <- '+ ethnicity_dich + bmi_9_z + m_educ_cont + m_age'
  
  # Fit model and pool estimates
  pool_fit <- function(adj) {
    
    if (adj=='base') {
      c
    } else {
      fit <- with(impset, lm(as.formula(paste(outc,'~',exp,covs1,covs2))))
    }
    p_fit <- mice::pool(fit) # pool results 
    mod <- summary(p_fit, 'all', conf.int = 0.95) # extract relevant information
    
    # Add confidence intervals
    names(mod)[which(names(mod) %in% c('2.5 %','97.5 %'))] <- c('lci','uci')
    #mod$lci <- (mod$estimate - 1.96*mod$std.error)
    #mod$uci <- (mod$estimate + 1.96*mod$std.error)
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
exp_mod <- function(exp, outcomes = c('tbv_13_z', 'gmv_13_z', 'mfa_13_z', 'mmd_13_z'), sex='', fullsample=F) {
  mod <- data.frame()
  for (out in outcomes) { p <- pool_mod(out, exp, sex=sex, fullsample=fullsample) 
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
                 'sign_raw','sign_fdr','lci','uci','rsq','rsq_adj',
                 'm','riv','lambda','fmi','ubar','b','t','dfcom')
  mod <- mod[, col_order]
  
  return(mod)
}

# ==============================================================================

for (e in c('imt','dis','sbp','dbp')) {
  ze <- paste0(e,'_9_z')
  # MAIN ANALYSES
  assign(e, exp_mod(ze))
  # SUBCORTICAL VOLUMES
  assign(paste0(e,'_subc_vol'), exp_mod(ze, subc_vol))
  # WHITE MATTER TRACTS
  assign(paste0(e,'_tract_fa'), exp_mod(ze, trcts_FA))
  assign(paste0(e,'_tract_md'), exp_mod(ze, trcts_MD))
  # SEX-STRATIFIED
  assign(paste0(e,'_f'), exp_mod(ze, sex='girl'))
  assign(paste0(e,'_m'), exp_mod(ze, sex='boy'))
  # non standardized (original scale)
  assign(paste0(e,'_unscaled'), 
         exp_mod(paste0(e,'_9'), 
                 paste0(c('tbv', 'gmv', 'mfa', 'mmd'),'_13')))
  # full sample (with imputed brain)
  assign(paste0(e,'_fullsamp'),  exp_mod(ze, fullsample = T))
  #assign(paste0(e,'_fullsamp.f'),exp_mod(ze, fullsample = T, sex='girl'))
  #assign(paste0(e,'_fullsamp.m'),exp_mod(ze, fullsample = T, sex='boy'))
  # sex interaction 
  assign(paste0(e,'_sex_inter'), exp_mod(paste0('sex * ',e,'_9_z')))
  # main outcomes at 9 
  assign(paste0(e,'_brain9'), exp_mod(ze, paste0(c('tbv', 'gmv', 'mfa', 'mmd'),'_9_z'), fullsample = T))
  # other outcomes 
  assign(paste0(e,'_other'), exp_mod(ze, paste0(c('cortex', 'subcort', 'wmv'),'_9'), fullsample = T))
}

# INTERACTION
for (e in c('sbp','dbp')) {
  assign(paste0(e,'_inter'), exp_mod(paste0('imt_9_z * ',e,'_9_z')))
}
inter <- rbind(sbp_inter, dbp_inter)

# ==============================================================================

main_modls <- list("IMT" = imt, "Dis" = dis, "SBP" = sbp, "DBP" = dbp, 
                   "IMT_orig" = imt_unscaled, "Dis_orig" = dis_unscaled, 
                   "SBP_orig" = sbp_unscaled, "DBP_orig" = dbp_unscaled,
                   "IMT_full" = imt_fullsamp, "Dis_full" = dis_fullsamp, 
                   "SBP_full" = sbp_fullsamp, "DBP_full" = dbp_fullsamp)

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
                   "Inter" = inter, 
                   #"IMT_full_f" = imt_fullsamp.f, "Dis_full_f" = dis_fullsamp.f, 
                   #"SBP_full_f" = sbp_fullsamp.f, "DBP_full_f" = dbp_fullsamp.f,
                   #"IMT_full_m" = imt_fullsamp.m, "Dis_full_m" = dis_fullsamp.m, 
                   #"SBP_full_m" = sbp_fullsamp.m, "DBP_full_m" = dbp_fullsamp.m,
                   "IMT_sexint" = imt_sex_inter, "Dis_sexint" = dis_sex_inter, 
                   "SBP_sexint" = sbp_sex_inter, "DBP_sexint" = dbp_sex_inter, 
                   "IMT_brain9" = imt_brain9, "Dis_brain9" = dis_brain9,
                   "SBP_brain9" = sbp_brain9, "DBP_brain9" = dbp_brain9,
                   "IMT_other" = imt_other, "Dis_other" = dis_other,
                   "SBP_other" = sbp_other, "DBP_other" = dbp_other)

openxlsx::write.xlsx(supp_modls, file = file.path(genrpath, 'results',
                                             paste0(Sys.Date(), "_SuppResults.xlsx")), overwrite = T)

# ==============================================================================
# LONGITUDINAL BRAIN ANALYSIS 
# ==============================================================================
genrpath <- dirname(file.choose()) # project folder
# Load imputed datasets
fulset <- readRDS(file.path(genrpath,'results','imp_full.rds'))

invisible(lapply(c('ggplot2','nlme','lme4','jtools','reghelper'), require, character.only = T));
# library(ggplot2)
# library(nlme)
# library(lme4)
# library(jtools)
# library(reghelpler)

imp1 <- mice::complete(fulset, 0)
# imp1$diff_dich <- ifelse(imp1$tbv_9_z > imp1$tbv_13_z, 0, 1)
# imp1$diff <- imp1$tbv_13_z - imp1$tbv_9_z

# imp1$SBP <- ifelse(imp1$sbp_9_z > 2, 'high', ifelse(imp1$sbp_9_z < -2, 'low','mean'))

outcs <- c('tbv', 'gmv', 'mfa', 'mmd')

implong <- reshape(data=imp1, idvar='IDC',
                   varying = c(paste0(c('age_mri', outcs), '_9'),
                               paste0(c('age', outcs), '_13')),
                   v.name=c('AGE','TBV','GMV','MFA','MMD'),
                   timevar='Visit',
                   direction="long")

write.csv(implong, file.path(genrpath,'DATA','Data_long.csv'))

o = 'TBV'
e = 'sbp_9_z'
f = as.formula(paste(o,'~ AGE *',e,'+ sex + height_9 + (1|IDC)'))

fit <- lme4::lmer(f, data = implong, REML = F)

summary(fit)
sm <- jtools::summ(fit, digits=3)
ss <- reghelper::simple_slopes(fit)

reghelper::graph_model(fit, y= TVB, x=AGE, lines=sbp_9_z)

plot(fit)

library(ggeffects)
library(dplyr)
PlotSimple <- ggpredict(fit, terms = c("sbp_9_z", "IDC"), type = "re") %>% #create a 'plot object', select the predictor variable (here: CN_GMC) and the clustering variable (here: Country)
  plot() +
  labs(x = "sbp", y = "j", title = "Plot Skill Camp") + #name your axes and add a title (if you want to)
  scale_fill_manual(values = mycolors100) #specify that you want to use the expanded color palette that you just created (here: mycolors100)


# ==============================================================================
# NON-LINEARITY ANALYSIS 
# ==============================================================================
# Define covariates
covs1 <- ' + sex + age_13 + height_9'
covs2 <- '+ ethnicity_dich + bmi_9_z + m_educ_cont + m_age'

form <- function(outc, exp, data, adj='full') {
  if (adj=='base') { covs2 <- '' }
  fit <- lm(as.formula(paste0(outc,' ~ ',exp,covs1,covs2)), data=data)
  # print(summary(fit))
  return(fit)
}

o = 'tbv_13_z'
e = 'sbp_9_z'

fit2 <- form(o, paste0('splines::ns(',e,', 5)'), imp1, adj='base')

fit2 <- lm(as.formula(paste0(o,'~ splines::ns(',e,', 5)')), imp1)

attach(imp1)

elims<-range(sbp_9_z)
#Generating Test Data
# e.grid<-seq(from=elims[1], to = elims[2])
e.grid <- with(imp1, expand.grid(sbp_9_z = seq(from=elims[1], to = elims[2]), 
                                sex = levels(sex), height_9 = mean(height_9), 
                                age_13 =mean(age_13)))

plot(sbp_9_z, tbv_13_z, col="grey",xlab=e, ylab=o)

lines(sbp_9_z, predict(fit2))

points(e.grid, predict(fit2, newdata = e.grid),col="darkgreen",lwd=2,type="l")

#adding cutpoints
abline(v=c(25,40,60),lty=2,col="darkgreen")

p <- ggplot(data, aes_string(e,o) ) + geom_point() +
  ggtitle(paste(toupper(substr(e,1,3)),'~',toupper(substr(o,1,3)))) +
  labs(y = toupper(substr(o,1,3)), x = toupper(substr(e,1,3))) +
  stat_smooth(method = lm, formula = as.formula(form(o, paste0('splines::ns(',e,', 5)'), imp1)))



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

library(gridExtra)

imps<- mids2datlist(fulset)
samp <- miceadds::subset_datlist(imps, subset = (imps$data$imt_9z < 2.5 & imps$data$imt_9z > -2.5) )

plot_nln <- function(outc, exp, set) {
  p = list()
  for (m in 1:20) {
    covs <- '+ sex + age_13 + height_9 + ethnicity_dich + bmi_9_z + m_educ_cont + m_age'
    fit <- lm(as.formula(paste0(outc,' ~ ns(',exp,', 5)',covs)), data=complete(set, m))
    eff <- effects::Effect(exp, fit) # transformation=list(inverse=exp))
    p[[m]] <- plot(eff, main=NULL, # xlim=c(-2.5, 2.5),
                   #main=paste(toupper(substr(outc,1,3)),'~',toupper(substr(exp,1,3)),'non-linear relationship'), 
                   xlab=paste(toupper(substr(exp,1,3)),'(z-score)'), 
                   ylab=paste(toupper(substr(outc,1,3)),'(z-score)') )
  }
  do.call(gridExtra::grid.arrange, c(p, ncol=5))
}

# for (e in paste0(c('imt','dis','sbp','dbp'),'_9_z')) {
#   for (o in paste0(c('tbv', 'gmv', 'mfa', 'mmd'),'_13_z')) {
#     cat('plot_nln("',o,'", "',e,'")\n', sep='')
#   }
# }

pdf(file.path(genrpath, 'results','nln_samp.pdf'),width=25, height=20)
plot_nln("tbv_13_z", "imt_9_z", set = mriset)
plot_nln("gmv_13_z", "imt_9_z", set = mriset)
plot_nln("mfa_13_z", "imt_9_z", set = dtiset)
plot_nln("mmd_13_z", "imt_9_z", set = dtiset) # *
plot_nln("tbv_13_z", "dis_9_z", set = mriset)
plot_nln("gmv_13_z", "dis_9_z", set = mriset)
plot_nln("mfa_13_z", "dis_9_z", set = dtiset)
plot_nln("mmd_13_z", "dis_9_z", set = dtiset)
plot_nln("tbv_13_z", "sbp_9_z", set = mriset)
plot_nln("gmv_13_z", "sbp_9_z", set = mriset)
plot_nln("mfa_13_z", "sbp_9_z", set = dtiset) # *
plot_nln("mmd_13_z", "sbp_9_z", set = dtiset) # *
plot_nln("tbv_13_z", "dbp_9_z", set = mriset) # *
plot_nln("gmv_13_z", "dbp_9_z", set = mriset) # *
plot_nln("mfa_13_z", "dbp_9_z", set = dtiset)
plot_nln("mmd_13_z", "dbp_9_z", set = dtiset) # *
dev.off()

# ==============================================================================
