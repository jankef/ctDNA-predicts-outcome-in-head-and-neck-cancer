# Early circulating tumour DNA changes predict outcomes in head and neck cancer patients under re-radiotherapy

## Background
Local recurrence after radiotherapy is common in locally advanced head and neck cancer (HNC) patients. Re-irradiation can improve local disease control, but disease progression remains frequent. Hence, predictive biomarkers are needed to adapt treatment intensity to the patient's individual risk. We quantified circulating tumour DNA (ctDNA) in sequential plasma samples and correlated ctDNA levels with disease outcome.

94 longitudinal plasma samples from 16 locally advanced HNC patients and 57 healthy donors were collected at re-radiotherapy baseline, after 5 and 10 radiation fractions, at irradiation end, and at routine follow-up visits. Plasma DNA was subjected to low coverage whole genome sequencing for copy number variation (CNV) profiling to quantify ctDNA burden.

CNV-based ctDNA burden was detected in 8/16 patients and 25/94 plasma samples. Ten additional ctDNA-positive samples were identified by tracking patient-specific CNVs found in earlier sequential plasma samples. ctDNA-positivity after 5 and 10 radiation fractions (both: log-rank, p = 0.050) as well as at the end of irradiation correlated with short progression-free survival (log-rank, p = 0.006). Moreover, a pronounced decrease of ctDNA towards re-radiotherapy termination was associated with worse treatment outcome (log-rank, p = 0.005). Dynamic ctDNA tracking in serial plasma beyond re-radiotherapy reflected treatment response and imminent disease progression. In 5 patients, molecular progression was detected prior to tumour progression based on clinical imaging.

Our findings emphasise that quantifying ctDNA during re-radiotherapy may contribute to disease monitoring and personalisation of adjuvant treatment, follow-up intervals, and dose prescription.

## Authors
Florian Janke1, 2, 3, Florian Stritzke3, 4, 5, 6, Katharina Dvornikovich4, Henrik Franke4, Arlou Kristina Angeles1, 2, 3, Anja Lisa Riediger1, 3, 7, 8, 9, Simon Ogrodnik1, 3, Sabrina Gerhardt1, 3, Sebastian Regnery3, 4, 5, 10, Philipp Schröter3, 4, 5, Lukas Bauer3, 4, 5, Katharina Weusthof3, 4, 5, 10, Magdalena Görtz7, 8, Semi Harrabi3, 4, 5, 6, 10, 11, Klaus Herfarth3, 4, 5, 6, 10, 11, Christian Neelsen12, Daniel Paech12,13, Heinz-Peter Schlemmer12, Amir Abdollahi3, 4, 5, 6, 10, 11, Sebastian Adeberg3, 4, 5, 6, 10, 11, 14, 15, 16, Jürgen Debus3, 4, 5, 6, 10, 11, Holger Sültmann1, 2, 3, 11, Thomas Held3, 4, 5, 6, 10, 11

1. Division of Cancer Genome Research, German Cancer Research Center (DKFZ), Heidelberg, Germany
2. German Center for Lung Research (DZL), TLRC Heidelberg, Heidelberg, Germany
3. National Center for Tumor Diseases (NCT), Heidelberg, Germany
4. Department of Radiation Oncology, Heidelberg University Hospital, Heidelberg, Germany
5. Heidelberg Institute of Radiation Oncology (HIRO), Heidelberg, Germany
6. Clinical Cooperation Unit Radiation Oncology, German Cancer Research Center (DKFZ), Heidelberg, Germany
7. Junior Clinical Cooperation Unit, Multiparametric Methods for Early Detection of Prostate Cancer, German Cancer Research Center (DKFZ), 69120 Heidelberg, Germany
8. Department of Urology, University Hospital Heidelberg, 69120 Heidelberg, Germany
9. Faculty of Biosciences, Heidelberg University, 69120 Heidelberg, Germany
10. Heidelberg Ion Beam Therapy Center (HIT), Heidelberg, Germany
11. German Cancer Consortium (DKTK), Heidelberg, Germany
12. Division of Radiology, German Cancer Research Center (DKFZ), Heidelberg, Germany
13. Department of Neuroradiology, Bonn University Hospital, Bonn, Germany
14. Department of Radiotherapy and Radiation Oncology, Marburg University Hospital, Marburg, Germany
15. Marburg Ion-Beam Therapy Center (MIT), Department of Radiotherapy and Radiation Oncology, Marburg University Hospital, Marburg, Germany
16. Universitäres Centrum für Tumorerkrankungen (UCT) Frankfurt – Marburg Germany


## Description
This repository contains the source code necessary to replicate the key aspects of the analysis presented in the study: Janke et al. 2024 Early circulating tumour DNA changes predict outcomes in head and neck
cancer patients under re-radiotherapy. This includes:
1. The calculation of ctDNA-informed CPA scores (ctCPA) scores from WisecondorX v1.2.5 outputs (_bins.bed files)
2. Recreation of main and supplementary figures presented in the publication

## Usage
The 'Main_script.R' script represents a step-by-step guide to calculate ctCPA scores and recreate the figures of this publication. All necessary data and scripts are provided in the 'data' and 'scripts' directories, respectively.

## Dependencies
### Software
- R version 4.2.0 (or newer)

### R packages
- ggplot2
- ggprism
- colorspace
- data.table
- stringr
- readxl
- dplyr
- svglite
- survival
- GGally
- showtext
- tidyverse

## License
This project is licensed under the APACHE 2.0 License - see the LICENSE.md file for details.

