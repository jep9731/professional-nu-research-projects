# Comparing Hippocampal Freesurfer Segmentation Methodologies

This repository contains an R workflow for comparing hippocampal and amygdala segmentation methodologies produced by FreeSurfer, including:

* T1-only
* T2 SPC
* T2 High-Resolution
* T1 + T2 SPC
* T1 + T2 High-Resolution

The script reads segmentation outputs, merges them into tidy datasets, normalizes volumes, produces summary files, and generates violin/boxplot visualizations and statistical analyses.

---

## 📁 Repository Contents

| File | Description |
|:-----:|:--------|
| `SA_amyg_hippo_seg.R` |	Main analysis script. Processes FreeSurfer segmentation files, combines outputs, normalizes volumes, and performs statistical modeling. |

**⚠️ Note:**
This repository does not include **MRI segmentation text files or demographic datasets**.
Users must supply their own FreeSurfer-derived `.v21.txt` and `aseg.stats`.

---

## 🔧 Requirements

**1. Software**
   * R (≥ 4.0)
   * FreeSurfer (for generating segmentation outputs)
  
**2. R Packages**
   * The script uses the following packages:

   ```tidyverse
      readr
      writexl
      freesurfer
      ggplot2
      stringr
      car
      emmeans
      FSA
      effsize
      knitr
      kableExtra
      gt
      openxlsx
      sjPlot
      performance
      effects
      ggsignif
      rstatix```

Install everything with:

```install.packages(c(
  "tidyverse", "readr", "writexl", "ggplot2", "stringr", "car",
  "emmeans", "FSA", "effsize", "knitr", "kableExtra", "gt",
  "openxlsx", "sjPlot", "performance", "effects", "ggsignif",
  "rstatix"
))```

FreeSurfer R tools:

```
# Install if necessary
# install.packages("freesurfer")
library(freesurfer)```

---

## 📂 Required Input Data

Users must supply:

**1. Segmentation output files**

Must include files matching patterns such as:
    *-T1.v21.txt
    *_T2spcnorm.v21.txt
    *_T2highreshipp.v21.txt
    *-T1*_T2spcnorm.v21.txt
    *-T1*_T2highreshipp.v21.txt

**2. aseg.stats files**

Used for extracting Estimated Total Intracranial Volume (eTIV).

---

## 📁 Expected Folder Structure

Your directory should contain segmentation and stats files organized by subject folders.
Example:
```
data_final/
    SA001/
        lh.hippocampus-T1.v21.txt
        rh.amyg-T1.v21.txt
        lh.hippocampus_T2highreshipp.v21.txt
        ...
    SA002/
        ...

sa_hippseg/
    SA001/stats/aseg.stats
    SA002/stats/aseg.stats
    ...

data/
    demographics.csv
    subtype_data.csv
```

You may use different paths, just update them in the script.

---

## ▶️ Running the Script

**1.** Clone the repository:
```
git clone https://github.com/<yourusername>/<yourrepo>.git
cd <yourrepo>
```

**2.** Open R and run:
```
source("SA_amyg_hippo_seg.R")
```

The script will:

* ✓ Load all FreeSurfer segmentation files
* ✓ Extract subject ID, hemisphere, subregion, method
* ✓ Merge all segmentation methods into a unified dataset
* ✓ Add eTIV and normalize volumes
* ✓ Produce violin and boxplots for:
    * Hippocampal subregions
    * Amygdala subregions
    * Comparisons across segmentation methods
    * Subtype-based comparisons

* ✓ Perform statistical analyses:
    * Shapiro–Wilk normality tests
    * Levene’s test
    * Kruskal–Wallis
    * Pairwise Wilcoxon signed-rank (with Bonferroni correction)
    * Linear modeling with interaction terms
    * Estimated marginal means (EMMs)

* ✓ Export outputs:
    * Outputs might include:
        * combined_data_final.csv
        * combined_amygdala_data_final.csv
        * combined_hippocampus_data_final.csv
        * merged_hipp_final.xlsx
        * merged_amyg_final.xlsx
        * kruskal_results.csv
        * wilcoxon_results.csv
        * plots/*.png
        * outputs/lm_model_summary.xlsx

---

## 📊 Example Outputs

The script will generate:

* Violin plots of segmentation method differences
* Boxplots by subtype
* Tables summarizing Wilcoxon tests
* Linear model effect plots
* Poster-style figures

(These are not included in the repository, users must run the script with their own data.)

---

## ⚠️ Notes & Limitations

* The script assumes FreeSurfer v7+ hippocampal/amygdala subfield outputs.
* File naming conventions must follow the patterns used in the script.
* Users must adjust their file paths if their data structure differs.
* Metadata cleaning assumes specific REDCap-like variable structures; modify as needed.

---

## 📬 Contact

If you encounter issues or want help adapting the pipeline to your dataset, feel free to open an issue or submit a pull request.
