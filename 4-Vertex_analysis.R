

# export FREESURFER_HOME=/mnt/appl/tools/freesurfer/6.0.0
# source $FREESURFER_HOME/SetUpFreeSurfer.sh

library(devtools)

genrpath <- '/home/r057600/IMT-dis-brain' # project folder
n_imps <- 20 # numebr of imputed datasets
implist_file <- 'implist_qdecr.rds' # dataset file name (note: RDS assumed format)

# Read in phenotype dataset
datalist <- readRDS(file.path(genrpath, implist_file))
# Read in cleaned brain dataset
core <- read.csv(file.path(genrpath, 'Brain.csv'))

id <- core[,c('idc','folders_f13')]

for (m in 1:n_imps) { 
datalist[[m]]$idc <- datalist[[m]]$IDC
datalist[[m]] <- merge(id, datalist[[m]], by='idc', all.x=F, all.y=T)
# does not have thickness lh
datalist[[m]] <- datalist[[m]][datalist[[m]]$folders_f13 != 'sub-2363_ses-F13',]
}

# load the QDECR package (this will become a normal library load eventually)
# load_all("/data/mounts/scs-fs-20/kpsy/genr/users/slamballais/R/QDECR_beta/")

library(QDECR)

for (exp in c('imt','dis','sbp','dbp')) {
    form <- as.formula(paste0('qdecr_thickness ~ ',exp,'_10_z + sex + age_mri_13 + height_10')) # + ethn_dich + bmi_10_z + m_educ_cont + m_age
    for (hemisf in c('rh','lh')) {
    # the actual QDECR call:
        a <- qdecr_fastlm(form, data = datalist, 
              id = "folders_f13",
  	          hemi = hemisf,
              clobber = TRUE, # override analysis you have already done;
 	            n_cores = 4,    # parallelization, by default the maximum is 4 but can go up to 32
              dir_subj = '/mnt/data/genr/mrdata/GenR_MRI/bids/derivatives/freesurfer/6.0.0/', 
              dir_fshome = '/mnt/appl/tools/freesurfer/6.0.0',
  	          dir_tmp = '/home/r057600/IMT-dis-brain/qdecRes3/qdec_tmp', # "/dev/shm", = shared memory
  	          dir_out = "/home/r057600/IMT-dis-brain/qdecRes3", # folder where it is saved in the QdecR structure
  	          project = paste0(exp,".M2")) #  actual name of the analysis run. Suggestion: add M1 to the name if it is a crude model, M2 if it is a minimally adjusted model and with M3 if it is a fully adjusted model.
        saveRDS(a, paste0('/home/r057600/IMT-dis-brain/qdecRes3/',exp,'_',hemisf,'.rds'))
    }
}

# save out above code into R script. transfer to scs. run: Rscript ./myScript.

# available QDECR measures/outcomes:
# - area
# - area.pial
# - curv
# - jacobian_white
# - pial_lgi
# - sulc
# - thickness
# - volume
# - w-g.pct.mgh
# - white.H
# - white.K

# QDECR Tidbits:
# extract summary information from your QDECR output...
