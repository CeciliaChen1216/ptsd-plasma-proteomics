# Compatibility aliases for retained supporting modules. No analysis runs here.
P <- getOption("ptsd.current.paths")
if(is.null(P)) stop("Source 00_START.R first")
pqtl_path <- P$pqtl_rds
pqtl_base <- P$protein_raw_base
ptsd_path <- P$source_casecontrol
mdd_path <- P$mdd_gwas
anx_path <- P$anx_gwas
crp_path <- P$crp_gwas
ld_ref_path <- P$ld_prefix
result_dir <- P$result_dir
fig_dir <- file.path(P$project_root,"Figures")
table_dir <- file.path(P$result_dir,"tables")
validation_dir <- file.path(P$project_root,"validation_data")
hpa_dir <- file.path(P$project_root,"HPA")
cs_dir <- file.path(P$result_dir,"susie_credible_sets")
smr_bin <- file.path(validation_dir,"external_data/smr/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe")
plink_bin <- "/path/to/PLINK1_if_required"
hm3_path <- P$hm3
ldsc_dir <- P$ldsc
ref_sumstats <- P$sumstats_reference
ptsd_smr_ma <- P$ptsd_ma
brain_eqtl_dir <- file.path(validation_dir,"external_data/brain_eqtl/BrainMeta_cis_eqtl_summary")
blood_eqtl_dir <- file.path(validation_dir,"external_data/blood_eqtl")
# These optional historical inputs were never fully configured in the source release.
pcl_path <- "/path/to/pcl_quantitative.gz"
aam_path <- "/path/to/aam_ptsd.gz"
hna_path <- "/path/to/hna_ptsd.gz"
decode_supp <- "/path/to/deCODE/assocvariants.annotated.txt.gz"
# Annotation constants, not replacements for per-variant N in current SEM.
n_pqtl <- 54219
n_ptsd <- 1222882 # Header sum only, not a reconciled unique-participant count.
n_mdd <- 1639572
n_anx <- 851685
n_decode <- 35559
candidates <- c("AKT3","CD40","CGREF1","FES","FURIN","SIRPA","CD101","KHK","SNX18","UBE2L6")
protein_folders <- c(AKT3="AKT3_Q9Y243_OID21197_v1_Oncology",CD40="CD40_P25942_OID20724_v1_Inflammation",
 CGREF1="CGREF1_Q99674_OID20152_v1_Cardiometabolic",FES="FES_P07332_OID21207_v1_Oncology",
 FURIN="FURIN_P09958_OID21514_v1_Oncology",SIRPA="SIRPA_P78324_OID20304_v1_Cardiometabolic",
 CD101="CD101_Q93033_OID31480_v1_Oncology_II",KHK="KHK_P50053_OID30241_v1_Cardiometabolic_II",
 SNX18="SNX18_Q96RF0_OID31476_v1_Oncology_II",UBE2L6="UBE2L6_O14933_OID30321_v1_Cardiometabolic_II")
fdr_threshold <- 0.05
nominal_p <- 0.05
f_stat_min <- 10
coloc_pph4_min <- 0.8
heidi_p_min <- 0.05
cis_window <- 500000
presso_n_perm <- 5000
steiger_alpha <- 0.05
