# arterial-health-brain
This project examines the association between arterial health (i.e., carotid intima-media thickness, carotid distensibility and blood pressure) at age 10 and brain morphology (i.e., volume and microstructural integrity) at age 13 years. See [preregistration](https://osf.io/ryc7e).

# Pipeline
## 0. Data extraction and imputation: 
- **0-Brain_prep.R** extracts relevant brain data (age 9 and 13 years). This includes the "core data" containing the relevant exclusion criteria (e.g., QC, incidental findings) and the files containing summary measures from anatomical (i.e., total and regional volumes) and DWI scanners (i.e., global and tract-specific FA and MD values). Where values are provided per hemisphere, the script returns an average of left and right. All relevant variables are saved in `Brain.csv`, including:
    - Total brain volume ('tbv_13')
    - Grey matter volume ('gmv_13')
    - Global fractional anisotropy  ('mfa_13')
    - Global Mean Diffusivity ('mmd_13')
- **0-Data_prep.R** extracts all exposure, covariate and auxiliary data, and runs the missing data imputation. 
## 1. Descriptives 
(**1-Descriptives.R**)
## 2. Regression diagnostics 
(**2-Regr_diagnostics.R**)
## 3. Main analyses 
(**3-Main_analyses.R**)
## 4. Figures 
(**4.Figures.ipynb**)


