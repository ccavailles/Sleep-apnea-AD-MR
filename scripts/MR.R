library(readr)
library(tidyverse)
library(TwoSampleMR)
library(LDlinkR) 
library(gt)


### EXPOSURE DATA ----

## Import exposure GWAS SumStats
sa_path = "data/SA_reformatted.tsv.gz"
sa_exp_dat <- read_tsv(sa_path)
head(sa_exp_dat)

## Format data to TwoSampleMR format
exposure <- sa_exp_dat %>%
  format_data(.,
              type = "exposure",
              snps = NULL,
              header = TRUE,
              phenotype_col = "TRAIT",
              snp_col = "SNP",
              beta_col = "BETA",
              se_col = "SE",
              eaf_col = "FRQ",
              effect_allele_col = "A2",
              other_allele_col = "A1",
              pval_col = "P",
              samplesize_col = "N",
              chr_col = "CHR",
              pos_col = "BP",
              log_pval = FALSE
  ) %>%
  as_tibble()


## Perform LD clumping on SNP data, filter SNPs to make it run faster
exposure_clump <- exposure %>% 
    clump_data(.,
             clump_kb = 10000,
             clump_r2 = 0.001,
             clump_p1 = 1,
             clump_p2 = 1,
             pop = "EUR"
  )

exposure_dat <- filter(exposure_clump, pval.exposure < 5e-8)



### OUTCOME DATA ----

# Define column types for summary statistics
coltypes = cols(
  ID = col_character(),
  CHROM = col_double(),
  POS = col_double(),
  REF = col_character(),
  ALT = col_character(),
  AF = col_double(),
  TRAIT = col_character(),
  BETA = col_double(),
  SE = col_double(),
  Z = col_double(),
  P = col_double(),
  N = col_double(),
  OR = col_double(),
  OR_L95 = col_double(),
  OR_U95 = col_double(),
  DIR = col_character(),
  G1000_ID = col_character(),
  G1000_VARIANT = col_character(),
  DBSNP_ID = col_character(),
  DBSNP_VARIANT = col_character(),
  OLD_ID = col_character(),
  OLD_VARIANT = col_character()
)

## Import outcome GWAS SumStats
ad_path = "data/Kunkle2019load_stage123.chrall.CPRA_b37.tsv.gz"
ad_exp_dat <- read_tsv(ad_path,  comment = "##",  col_types = coltypes, 
                       col_select = c(DBSNP_ID, CHROM, POS, REF, ALT, AF, BETA, SE, Z, P, N, TRAIT))
head(ad_exp_dat)


# Format outcome
outcome <- ad_exp_dat %>%
  format_data(.,
              type = "outcome",
              snps = NULL,
              header = TRUE,
              phenotype_col = "TRAIT",
              snp_col = "DBSNP_ID",
              beta_col = "BETA",
              se_col = "SE",
              eaf_col = "AF",
              effect_allele_col = "ALT",
              other_allele_col = "REF",
              pval_col = "P",
              samplesize_col = "N",
              z_col = "Z",
              chr_col = "CHROM",
              pos_col = "POS",
              log_pval = FALSE
  ) %>%
  as_tibble()


## Identify Proxy Variants

# extract exposure SNPs present in outcome
outcome_clump <- semi_join(outcome, exposure_dat, by = "SNP")

# Exposure SNPs not present in outcome
exp_snps_wo <- anti_join(exposure_dat, outcome, by = "SNP")

outcome_dat <- outcome




### HARMONIZE EXPOSURE - OUTCOME DATASETS ----

mr_dat <- harmonise_data(exposure_dat, outcome_dat, action = 2) %>% as_tibble() %>%
  mutate(
    apoe_region = case_when(
      chr.outcome == 19 & between(pos.outcome, 44912079, 45912079) ~ TRUE,
      TRUE ~ FALSE
    ),
    gws.outcome = ifelse(pval.outcome < 5e-8, TRUE, FALSE),
    mr_keep = ifelse(mr_keep == FALSE | apoe_region == TRUE | gws.outcome == TRUE, FALSE, TRUE)
  ) %>%
  filter(pval.exposure < 5e-8)



### TWO SAMPLE MR ----

## MR

mr_res <- mr(mr_dat, method_list = c(
  "mr_egger_regression", "mr_weighted_median", "mr_ivw_fe", "mr_weighted_mode"))
generate_odds_ratios(mr_res)

res_single <- mr_singlesnp(mr_dat, all_method = c("mr_ivw_fe", "mr_egger_regression", 
                                                  "mr_weighted_median", "mr_weighted_mode")) %>% as_tibble()


## Visualization

# MR scatter plot

scatter_p <- mr_scatter_plot(mr_res, mr_dat) 
scatter_p

scatter_out_p <- scatter_p[[1]] + theme_bw() + 
  guides(color=guide_legend(ncol =1)) + 
  theme(
    text = element_text(size = 8), 
  )
scatter_out_p

# MR forrest plot 

forrest_p <- mr_forest_plot(res_single)
forrest_p[[1]]




### DIAGNOSTICS AND SENSITIVITY ANALYSES ----


## Pleiotropy

res_pleio <- mr_pleiotropy_test(mr_dat)

res_pleio %>%
  select(-id.exposure, -id.outcome, -outcome, -exposure) %>%
  gt() %>%
  fmt_number(
    columns = c('egger_intercept', 'se')
  ) %>%
  fmt_number(
    columns = pval,
    rows = pval > 0.001,
    decimals = 3
  ) %>% 
  fmt_scientific(
    columns = pval,
    rows = pval <= 0.001,
    decimals = 1
  )


# Funnel plots

funnel_p <- mr_funnel_plot(res_single)
funnel_out_p <- funnel_p[[1]] + theme_bw() + 
  guides(color=guide_legend(ncol =1)) + 
  theme(
    text = element_text(size = 8), 
  )
funnel_out_p



## Heterogeneity

res_het <- mr_heterogeneity(mr_dat, method_list = c("mr_egger_regression", "mr_ivw"))

res_het %>%
  select(-id.exposure, -id.outcome, -outcome, -exposure) %>%
  gt() %>%
  fmt_number(
    columns = Q
  ) %>%
  fmt_number(
    columns = Q_pval,
    rows = Q_pval > 0.001,
    decimals = 3
  ) %>% 
  fmt_scientific(
    columns = Q_pval,
    rows = Q_pval <= 0.001,
    decimals = 1
  )


## Outliers 

# Leave-one-out

res_loo <- mr_leaveoneout(mr_dat, method = mr_ivw_fe) %>% as_tibble()

loo_p <- mr_leaveoneout_plot(res_loo)
loo_p[[1]]


## F stat

## F statistic. Burgess et al 2011
f_stat = function(N, K, R){
  f = ((N-K-1) / K) * (R/(1-R))
  f
}

## Proportion of phenotypic variance explained by SNP 
## https://doi.org/10.1371/journal.pone.0120758.s001

snp.pve <- function(eaf, beta, se, n){
  (2*eaf*(1 - eaf)*beta^2) / (2 * beta * eaf * (1-eaf) + se^2 * 2 * n * eaf * (1-eaf))
}

f_res <- mr_dat %>%
  group_by(exposure, outcome) %>%
  filter(mr_keep == TRUE) %>%
  select(SNP, exposure, outcome, effect_allele.exposure, eaf.exposure, beta.exposure, se.exposure) %>%
  mutate(
    samplesize.exposure = 351316, 
    pve = snp.pve(eaf.exposure, beta.exposure, se.exposure, samplesize.exposure), 
    f = f_stat(samplesize.exposure, 1, pve),
    # f = abs(beta.exposure)^2 / se.exposure^2
  ) %>% 
  summarise(
    pve = sum(pve, na.rm = T), 
    k = n(), 
    samplesize = max(samplesize.exposure), 
    f = mean(f, na.rm = T),
  ) 
f_res



### RADIAL-MR ----

library(RadialMR)

# Format data 
radial_dat <- mr_dat %>% filter(mr_keep == T) %>% TwoSampleMR::dat_to_RadialMR()
rmr_dat <- radial_dat$exposure.AD

# Run Radial MR
bonff = 0.05/nrow(rmr_dat) # Bonferonni Correction

radial_ivw_res <- ivw_radial(rmr_dat, alpha = bonff) 

radial_egger_res <- egger_radial(rmr_dat, alpha = bonff) 
