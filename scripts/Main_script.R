#   Title:        Early circulating tumour DNA changes predict outcomes in head and neck cancer patients under re-radiotherapy
#
#   Description:  This repository contains the source code necessary to replicate the key aspects of the analysis presented
#                 in the study: Janke et al. 2024 Early circulating tumour DNA changes predict outcomes in head and neck
#                 cancer patients under re-radiotherapy
#
#----------------------------------------

#----------------------------------------
#   Requirements:
#----------------------------------------
library(ggplot2)
library(ggprism)
library(colorspace)
library(data.table)
library(stringr)
library(readxl)
library(dplyr)
library(svglite)
library(survival)
library(GGally)
library(showtext)
library(tidyverse)

#----------------------------------------
#   Set parameters and load data:
#----------------------------------------
#     1) Set root
root.path <- "/omics/odcf/analysis/OE0309_projects/care/hg19/ctdna-predicts-outcome-in-head-and-neck-cancer"

#     2) Load-in patient characteristics <patient_info.xlsx>
info <- readxl::read_excel(file.path(root.path, "data", "patient_info.xlsx"))

#     3) Load-in CPA scores
data <- data.frame(fread(file.path(root.path, "data", "CPA-scores.tsv")))

#     4) Set font
font_add(family = "Arial", regular = file.path(root.path, "data", "Arial.ttf"))
showtext_auto()
showtext_opts(dpi = 300)


#----------------------------------------
#   Calculate ctCPA scores:
#----------------------------------------
source(file.path(root.path, "scripts", "ctCPA_score.R"))
ctCPA <- ctDNA_informed_CPA(input = file.path(root.path, "data", "copy-number-profiles"), timepoints = c("Baseline", "Baseline|5-Fx", "Baseline|5-Fx|10-Fx"), data = data)


#----------------------------------------
#   Figures:
#----------------------------------------

#   Fig. 1: (A) Box plot comparing CPA scores between healthy donors and head and neck cancer (HNC) patients. The maximum CPA score
#           across healthy donors served as ctDNA detectability threshold (dotted line; CPA = 1.044). ctDNA-positive and –negative
#           samples are colored black and gray, respectively. ctDNA-positive and total sample number per group is given in brackets.
#           Box plots represent median, upper and lower quartile with Tukey Whiskers. (B) Summary of recurrent CNVs in HNC patients
#           with detectable ctDNA (n = 8). The y-axis represents the frequency of a detected copy number state at the chromosomal coordinate
#           given on the x-axis. All longitudinal plasma samples were considered, however, recurrently detected CNVs in two or more samples
#           of the same patient were only counted once. The CNV frequency (i.e., the number of occurrences in the 8 patients assessed) of
#           amplifications and deletions is given in red and green, respectively. Areas shaded in gray represent the q-arm of the respective
#           chromosome. Genes associated with HNC tumourigenesis are labeled.
source(file.path(root.path, "scripts", "Fig_1.R"))
fig1 <- cpa_and_cnv(data = data, cytoband = file.path(root.path, "data", "cytoBand_hg19.txt"), input = file.path(root.path, "data", "copy-number-profiles"))


#   Fig. 2: Head and neck cancer patient cohort overview highlighting the detectability of ctDNA copy number variation (CNV) analysis from
#           low-coverage whole genome sequencing. Green and blue colors indicate detectable CNVs by de novo CNV-calling (CPA score) and using
#           information from previous serial plasma samples of the same patient (ctCPA score), respectively. Missing squares represent
#           unavailable plasma samples. The number of ctDNA-positive and total number of samples per time point is given in brackets.
source(file.path(root.path, "scripts", "Fig_2.R"))
fig2 <- oncoprint(table = ctCPA, volume = info)


#   Fig. 3: Association between ctDNA levels and progression-free survival after re-radiotherapy. Progression-free survival (PFS) according
#           to ctDNA detectability, as assessed via ctCPA scores, in plasma specimens collected after 5 fractions (A), 10 fractions (B) and
#           at the end of re-radiotherapy (re-RT; C). (D) Association between PFS and CPA score change from samples taken after 5 re-RT fractions
#           to re-RT end. Relative changes were calculated using CPA scores or ctCPA scores (if available). A decrease of ≥30% was used for
#           partitioning of patients. Groups were compared using the two-sided log-rank test.
source(file.path(root.path, "scripts", "Fig_3.R"))
fig3 <- survival_summary(table = info, ctCPA = ctCPA, timepoint = c("Baseline", "Fx_5", "Fx_10", "end_RT", "Fx_5_RT_end", "Age", "Sex", "Histology", "Modality", "GTV_ccm", "Smoking_status", "PD_info"))


#   Fig. 4: Representative patients (A-D) illustrating the utility of the ctCPA score for disease monitoring in HNC patients under re-radiotherapy
#           (re-RT). Administered radio- and systemic therapies are represented by shaded areas. Time points of disease progression (PD) and PD location
#           (if available) are denoted by vertical lines. Horizontal, dotted lines indicate the patient-specific detectability threshold of the ctCPA
#           score (i.e., maximum score across 57 healthy donors). Filled dots highlight samples with detectable ctCPA scores. Early signs of PD in months
#           (i.e., increasing ctCPA scores prior to radiographic tumour progression) are highlighted in gray.
source(file.path(root.path, "scripts", "Supplementary_Fig_6.R"))
fig4 <- kinetics(patID = c("P007", "P001", "P008", "P018"), ctCPA = ctCPA)



#----------------------------------------
#   Supplementary Figures:
#----------------------------------------

#   Supplementary Fig. 2: Swimmer plot illustrating plasma collection time points, administered therapy regimes
#                         as well as time of disease progression during the first 12 months after initiation of re-radiotherapy.
#                         Patients are ordered by re-radiotherapy type and duration of progression-free survival. Asterisks
#                         indicate patients that deceased during or after the illustrated 12 months follow-up period.
source(file.path(root.path, "scripts", "Supplementary_Fig_2.R"))
figS2 <- swimmer(input = file.path(root.path, "data", "patient_info.xlsx"), data = data)


#   Supplementary Fig. 3: CPA score change from baseline to 5 re-RT fractions (5-Fx; A), baseline to 10-Fx (B), 5-Fx to 10-Fx (C),
#                         baseline to the end of re-RT (D), 5-Fx to re-RT end (E), and 10-Fx to re-RT end. CPA score changes were
#                         compared between volumetric modulated arc therapy (VMAT) and carbon ion radiotherapy (CIRT).
source(file.path(root.path, "scripts", "Supplementary_Fig_3.R"))
figS3 <- modality_diff(table = data, info = info, timepoint = list(c("Baseline", "5-Fx"), c("Baseline", "10-Fx"), c("5-Fx", "10-Fx"), c("Baseline", "RT-end"), c("5-Fx", "RT-end"), c("10-Fx", "RT-end")))


#   Supplementary Fig. 4: (A) Comparison of CPA scores in head and neck cancer patients stratified by the cohort´s median GTV (18.9 cm3).
#                         ctDNA detectability is color-coded in black (positive) and gray (negative), according to the detectability threshold
#                         (dashed line; maximum CPA score in healthy donor cohort). ctDNA-positive and total sample number per group is given
#                         in parentheses. (B) CPA score comparison between by GTV of individual sampling time points (i.e., at baseline, after
#                         5, 10 fractions and at the end of re-RT as well as at follow-up visits). Box plots represent median, upper and lower quartile
#                         with Tukey Whiskers and statistical significance between groups was assessed by Mann-Whitney U tests. (C) Correlation
#                         between baseline CPA scores and GTVs. Linear regression line and corresponding Spearman correlation coefficients are
#                         given in the plot (black including and gray excluding outliers). Asterisks mark outliers removed for the linear regression
#                         shown in gray. Patient histology and information on disease dissemination (by TNM annotations) is highlighted.
source(file.path(root.path, "scripts", "Supplementary_Fig_4.R"))
figS4 <- gtv_comp(table = data, info = info, timepoint = c("Baseline", "5-Fx", "10-Fx", "RT-end", "MRI"))


#   Supplementary Fig. 5: Per patient ctCPA score dynamics from baseline to the end of re-radiotherapy (re-RT; A) and from re-RT end to the plasma
#                         sample collected closest to disease progression (PD; B). For patient P006, the second PD time point was considered as the
#                         first PD occurred already during re-RT. Statistical significance between groups was assessed by Wilcoxon´s paired test.
#                         Box plots represent median, upper and lower quartile with Tukey Whiskers. CPA, copy number profile abnormality.
source(file.path(root.path, "scripts", "Supplementary_Fig_5.R"))
figS5 <- ctCPA_changes(timepoint = c("Baseline|RT-end", "RT-end|PD"), ctCPA = ctCPA)


#   Supplementary Fig. 6: Remaining (ct)CPA score kinetics during and after re-radiotherapy (re-RT). Disease kinetics are represented by ctCPA scores
#                         (if available; A-C) or CPA scores (D-L). Therapies are displayed by shaded areas (re-RT and systemic therapies) and arrows (surgery).
#                         Disease progression (PD) time points (plus PD locations) are given by vertical lines. CtDNA detectability thresholds are given
#                         by the dotted lines and ctDNA-positive and –negative samples are highlighted via black and white fillings, respectively. Lead
#                         time (in months) is shown by gray, horizontal bars.
source(file.path(root.path, "scripts", "Supplementary_Fig_6.R"))
figS6 <- kinetics(patID = c("P003", "P006", "P016", "P009", "P002", "P004", "P010", "P011", "P012", "P014", "P015", "P017"), ctCPA = ctCPA)







