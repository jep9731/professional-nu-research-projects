# 🧠 Project B: Distinct Cognitive Trajectories in SuperAgers Identified by Repeated-Measures Latent Class Analysis
 
## 📌 Overview
 
This project applies **Repeated-Measures Latent Class Analysis (RMLCA)** to identify distinct longitudinal memory trajectory subtypes among SuperAgers and same-age controls. SuperAgers are individuals aged 80 and older who exhibit episodic memory performance comparable to adults 20–30 years younger. While some SuperAgers maintain stable memory over time, others show meaningful decline — this study aims to characterize that heterogeneity and identify demographic predictors of trajectory class membership.
 
> **Institution:** Mesulam Institute for Cognitive Neurology & Alzheimer's Disease, Northwestern University Feinberg School of Medicine
 
---
 
## 🎯 Research Question
 
Can RMLCA identify meaningfully distinct longitudinal memory trajectory subtypes among SuperAgers, and what demographic factors predict class membership?
 
---
 
## 👥 Sample
 
- **N = 262** participants enrolled in the Northwestern Alzheimer's Disease Research Center Clinical Core
- Minimum of **3 annual visits** with Rey Auditory Verbal Learning Test (RAVLT) data
- Up to **21 years** of longitudinal follow-up
 
---
 
## 🔢 Methods
 
### Statistical Approach
- **Repeated-Measures Latent Class Analysis (RMLCA):** A person-centered method identifying hidden subgroups sharing similar longitudinal patterns of change
- Models with **1 to 4 latent classes** compared using Bayesian Information Criterion (BIC) and convergence criteria
- **Time since baseline** used as the longitudinal metric
- **Class-membership covariates:** age, education, sex, race, and ethnicity
- **Empirical classification:** Individual RAVLT slopes classified as *Decliner* or *Stable*; agreement with RMLCA assessed via **Cohen's kappa**
 
### Primary Outcome
- **RAVLT Delayed Recall (reydrec):** Rey Auditory Verbal Learning Test score used to measure episodic memory at each annual visit
 
---
 
## 📊 Key Results
 
| Class | N | % | Slope | Interpretation |
|---|---|---|---|---|
| Class 1 | 134 | 51.1% | ≈ +0.004 SD/year | Stable trajectory; higher baseline performance |
| Class 2 | 128 | 48.9% | −0.193 SD/year (p < .001) | Declining ~1 word/year on raw RAVLT |
 
- **Best-fit model:** 2-class solution (BIC = 3938.33)
- **Older baseline age** predicted reduced probability of Class 1 (Stable) membership (β = −0.145, p = .032)
- **Concordance with empirical classification:** κ = 0.391 (p < .001); 81% of Decliners → Class 2; 65% of Stables → Class 1
 
---
 
## 📁 Repository Structure
 
```
project-b-superaging_lca/
├── longitudinal_sa.R       # Full analysis pipeline (REDCap API pull → RMLCA → visualizations)
├── Output/                 # Generated outputs (local only — data not publicly shareable)
│   ├── ravlt_trajectories_all.png
│   ├── ravlt_slope_distribution.png
│   ├── ravlt_rmlca_trajectories_2class.png
│   ├── ravlt_rmlca_trajectories_3class.png
│   ├── ravlt_by_empirical_group.png
│   ├── ravlt_by_empirical_group_and_age.png
│   ├── individual_slopes_ravlt.csv
│   └── summary_statistics_ravlt.csv
└── README.md
```
 
---
 
## 🛠️ Tools & Technologies
 
| Tool | Purpose |
|---|---|
| R | Primary analysis language |
| `lcmm` | Repeated-measures latent class analysis |
| `lme4` / `lmerTest` | Linear mixed-effects models for individual slopes |
| `tidyverse` | Data wrangling and visualization |
| `ggplot2` | Trajectory and distribution plots |
| `redcapAPI` / `httr` | REDCap API integration for longitudinal data extraction |
| `irr` | Cohen's kappa for concordance analysis |
 
---
 
## ⚠️ Data & Reproducibility Note
 
This project uses de-identified longitudinal data from the Northwestern Alzheimer's Disease Research Center Clinical Core, accessed via the REDCap API. Due to participant privacy protections, raw data cannot be shared publicly. The analysis script (`longitudinal_sa.R`) is fully documented and reproducible given appropriate REDCap credentials and access.
 
---
 
## 📬 Contact
 
Questions or inquiries? Reach out to [Joshua Pasaye](https://github.com/jep9731) or connect on [LinkedIn](https://www.linkedin.com/in/joshua-pasaye/).
