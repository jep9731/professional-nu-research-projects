# 🔬 Research Projects Collection

This repository hosts two independent research initiatives, each examining a distinct topic using modern analytical and computational methods. All projects are fully documented and structured for easy exploration and reproducibility.

---

## 📁 Projects Overview

| Project | Title |
|:-------:|:------|
| A | Comparing Hippocampal Freesurfer Segmentation Methodologies | 
| B | Distinct Cognitive Trajectories in SuperAgers Identified by Repeated-Measures Latent Class Analysis |
 
Each project contains:

* Source code and scripts.
* Dataset references or data folder.
* Result visualizations.
* A dedicated README.md for deeper details.

---

## 🛠️ Tech & Tools

Depending on the project, the projects may incorporate:

* **FSL**, **FreeSurfer**, **AFNI**
* **R**
* MRI structural and functional processing workflows
* Advanced data modeling techniques (e.g., **RMLCA**, **LCA**)

Consult the project folders for exact dependencies and setup instructions.

---

## 🚀 Getting Started

Clone the repository:

```
git clone https://github.com/jep9731/professional-nu-research-projects.git
cd professional-nu-research-projects
```

Then navigate into any project of interest:

```
cd project-a-freesurfer_methods
```

Each project’s README includes:

* Installation requirements
* Dataset notes (paths, acquisition details if shareable)
* Execution order and reproducibility steps

---

## 🧪 Project Details

**🧠 Project A — Hippocampal Segmentation Modalities** 

* Goal: Evaluate segmentation consistency across Freesurfer and ASHS tools
* Modalities compared:
    * T1-only
    * T2-only
    * T2 high-resolution (T2H)
    * T1 + T2
    * T1 + T2H
* Outputs: Volumetric comparisons, structural overlays, reliability metrics

**🧠 Project B — Distinct Cognitive Trajectories in SuperAgers (RMLCA)**
 
- **Goal:** Identify distinct longitudinal memory trajectory subtypes among SuperAgers and same-age controls using Repeated-Measures Latent Class Analysis (RMLCA)
- **Sample:** 262 participants from the Northwestern Alzheimer's Disease Research Center Clinical Core with a minimum of 3 annual visits (up to 21 years of follow-up)
- **Key Finding:** A 2-class solution best fit the data — a **Stable class** (51.1%, slope ≈ +0.004 SD/year) and a **Declining class** (48.9%, slope = −0.193 SD/year). Older baseline age predicted declining-class membership (β = −0.145, p = .032).
- **Tools:** R, lcmm, lme4, tidyverse, ggplot2, REDCap API
 
> ⚠️ *Data is not publicly shareable due to participant privacy protections. The analysis script is fully documented for reproducibility given appropriate REDCap access.*

---

## 📊 Results & Reproducibility

* All analysis notebooks and scripts provided
* Where permitted, anonymized data or links to publicly accessible datasets
* Figures, statistical tables, and exportable results included

---

## 📬 Contact
 
Questions or inquiries? Reach out to [Joshua Pasaye](https://github.com/jep9731/) or connect on [LinkedIn](https://www.linkedin.com/in/joshua-pasaye/).
