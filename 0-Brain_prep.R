# ==============================================================================
# ================= 0. Brain data extraction & cleaning ========================
# ==============================================================================

# Required packages
invisible(lapply(c('dplyr','stringr'), require, character.only = T));

# Define path
genrpath <- '/mnt/data/genr/mrdata/GenR_MRI/summary_data/'

# Core data for selection
core <- readRDS(file.path(genrpath, 'core/genr_mri_core_data_20220311.rds'))
core <- core[, -grep('f05', names(core))] # get rid of F05 wave

# read and clean brain data
read_brain <- function(age, outc, keep=c(), mergelr=c()) {
    # find file in appropriate location
    if (outc=='dti_stats_inc_glob') { fold <- 'dwi' } else { fold <- 'anat/freesurfer/6.0.0' }
    datapath <- file.path(genrpath, age, fold)
    f <- list.files(datapath)[grep(outc, list.files(datapath))]
    # read in file
    d <- readRDS(file.path(datapath, f))
    # subset columns
    if (length(keep)>0) {
        d <- select(d, starts_with(paste(c('idc',keep),sep='|')))
    }
    # merge left and right brain measures for structural brain (and clean left/right columns)
    if (length(mergelr)==2) {
        ns = gsub(mergelr[1],'',names(select(d, starts_with(mergelr[1]))))
        for (n in ns) {
            cols <- c(paste0(mergelr[1],n), paste0(mergelr[2],n))
            d[n] <- rowMeans(d[,cols])
            d <- d[,-which(names(d) %in% cols)]
        }
    # merge left and right brain measures for dti data 
    # NOTE: I only look at MD and FA in GenR tracts, so I also clean out all other columns
    } else if (length(mergelr)>2) {
      d <- d[, -grep('abcd|RD|AD|vol|vox|_mean_|_median_', names(d))] # get rid of other metrics
      for(t in mergelr[1:5]) { # note: all tracts of interest except fmi and fma, that are central
          ns = unique(gsub(paste0(t,'_l_|',t,'_r_'), '', names(select(d, starts_with(t)))))
          for (n in ns) {
              cols <- c(paste0(t,'_l_',n), paste0(t,'_r_',n))
              d[paste(t,n,sep='_')] <- rowMeans(d[,cols])
              d <- d[,-which(names(d) %in% cols)]
          }
      }
      d <- d[, -grep('_l_|_r_|mcp_', names(d))] # get rid of other tracts (i.e., abcd I think)
      names(d) <- gsub('dti_dipy_wls_wavg_', '', names(d)) # simplify names
    }
    
    return(d)
}

# Total brain values
keep_tb <- c('genr_tbv','TotalGrayVol','eTIV',
             'CortexVol','SubCortGrayVol','CerebralWhiteMatterVol')
tb09 <- read_brain('F09','tbv', keep=keep_tb)
tb13 <- read_brain('F13','tbv', keep=keep_tb)

# Segmentation 
sb09 <- read_brain('F09','aseg', mergelr=c('Left_','Right_'))
sb13 <- read_brain('F13','aseg', mergelr=c('Left_','Right_'))

# Parcelation
pb09 <- read_brain('F09','aparc', mergelr=c('lh_','rh_'))
pb13 <- read_brain('F13','aparc', mergelr=c('lh_','rh_'))

# White matter
tracts <- c('cgc', # Cincgulate gyrus 
            'cst', # Cortico-spinal tract
            'unc', # uncinate fasciculus
            'ilf', # inferior longitudinal fasciculus
            'slf', # superior longitudinal fasciculus
            'fma', # major forceps
            'fmi') # minor forceps
tw09 <- read_brain('F09','dti_stats_inc_glob', mergelr=tracts)
tw13 <- read_brain('F13','dti_stats_inc_glob', mergelr=tracts)

# Merge them all together
dfs <- list(core, tb09, tb13, sb09, sb13, pb09, pb13, tw09, tw13)
brain <- Reduce(function(x,y) merge(x = x, y = y, by = 'idc', all.x = T), dfs)

# Save
write.csv(brain, '~/IMT-dis-brain/Brain.csv')