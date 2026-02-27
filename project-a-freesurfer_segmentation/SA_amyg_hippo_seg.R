# ── File Paths ─────────────────────────────────────────────────────────────────
# Update these to match your local directory layout.
# Using here::here() anchors all paths to the project root.

data_dir    <- here("project-a-freesurfer_segmentation", "data")
seg_dir     <- here("project-a-freesurfer_segmentation", "data", "seg_final")
aseg_dir    <- here("project-a-freesurfer_segmentation", "data", "aseg")
output_dir  <- here("project-a-freesurfer_segmentation", "outputs")
plot_dir    <- here("project-a-freesurfer_segmentation", "plots")

# Create output directories if they don't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(plot_dir))   dir.create(plot_dir,   recursive = TRUE)

# ── Helper: Build Combined Dataframe from File List ───────────────────────────
# Reads a list of FreeSurfer segmentation .txt files and combines them into
# a single tidy dataframe with columns: sa_id, hemisphere, structure, method,
# subregion, volume_size.

read_seg_files <- function(files, method_label) {
  data_list   <- list()
  id_list     <- list()
  hem_list    <- list()
  struc_list  <- list()
  method_list <- list()

  for (file in files) {
    data        <- read.table(file)
    subject_id  <- basename(dirname(file))
    hemisphere  <- regmatches(file, regexpr("(rh|lh)", file))
    structure   <- sub(".*[lr]h\\.([a-z]+{1,5}).*", "\\1", file)

    data_list[[file]]   <- data
    id_list[[file]]     <- as.data.frame(rep(subject_id,    nrow(data)))
    hem_list[[file]]    <- as.data.frame(rep(hemisphere,    nrow(data)))
    struc_list[[file]]  <- as.data.frame(rep(structure,     nrow(data)))
    method_list[[file]] <- as.data.frame(rep(method_label,  nrow(data)))
  }

  combined <- cbind(
    bind_rows(id_list),
    bind_rows(hem_list),
    bind_rows(struc_list),
    bind_rows(method_list),
    bind_rows(data_list)
  )

  col_names <- c("sa_id", "hemisphere", "structure", "method", "subregion", "volume_size")
  combined %>% rename_with(~ col_names, all_of(names(combined)))
}

# ── Load Segmentation Files ───────────────────────────────────────────────────

# T1 only
files_T1 <- list.files(seg_dir, pattern = "*-T1.v21.txt$",
                        recursive = TRUE, full.names = TRUE)
combined_T1_final <- read_seg_files(files_T1, "T1")

# T2 SPC (T2 standard-space)
files_T2_spc <- list.files(seg_dir, pattern = "*_T2spcnorm.v21.txt$",
                            recursive = TRUE, full.names = TRUE)
files_T2_spc <- files_T2_spc[!grepl("-T1-", files_T2_spc)]
combined_T2_spc_final <- read_seg_files(files_T2_spc, "T2")

# T2 Highres
files_T2_highres <- list.files(seg_dir, pattern = "*_T2highreshipp.v21.txt$",
                                recursive = TRUE, full.names = TRUE)
files_T2_highres <- files_T2_highres[!grepl("-T1-", files_T2_highres)]
combined_T2_highres_final <- read_seg_files(files_T2_highres, "T2H")

# T1 + T2 SPC
files_T1_T2_spc <- list.files(seg_dir, pattern = "*-T1.*_T2spcnorm.v21.txt$",
                               recursive = TRUE, full.names = TRUE)
files_T1_T2_spc <- files_T1_T2_spc[grepl("-T1-", files_T1_T2_spc)]
combined_T1_T2_spc_final <- read_seg_files(files_T1_T2_spc, "T1_T2")

# T1 + T2 Highres
files_T1_T2_highres <- list.files(seg_dir, pattern = "*-T1.*_T2highreshipp.v21.txt$",
                                   recursive = TRUE, full.names = TRUE)
files_T1_T2_highres <- files_T1_T2_highres[grepl("-T1-", files_T1_T2_highres)]
combined_T1_T2_highres_final <- read_seg_files(files_T1_T2_highres, "T1_T2H")

# ── Load eTIV from aseg.stats ─────────────────────────────────────────────────
aseg_files <- list.files(aseg_dir, pattern = "aseg.stats",
                          recursive = TRUE, full.names = TRUE)

aseg_list    <- list()
aseg_list_id <- list()

for (file in aseg_files) {
  aseg_data        <- read_aseg_stats(file)
  aseg_list[[file]]    <- aseg_data[["measures"]]
  subject_id       <- strsplit(strsplit(file, "/stats")[[1]][1], paste0(basename(aseg_dir), "/"))[[1]][2]
  aseg_list_id[[file]] <- as.data.frame(rep(subject_id, nrow(aseg_data[["measures"]])))
}

combined_aseg_final <- cbind(bind_rows(aseg_list_id), bind_rows(aseg_list)) %>%
  rename_with(~ c("sa_id", "measure", "measure_long", "meaning", "size", "units"),
              all_of(names(.))) %>%
  filter(measure_long == "etiv") %>%
  select(sa_id, size) %>%
  mutate(size = as.numeric(size))

# ── Combine All Methods & Normalize by eTIV ──────────────────────────────────
combined_all <- bind_rows(
  combined_T1_final,
  combined_T2_spc_final,
  combined_T2_highres_final,
  combined_T1_T2_spc_final,
  combined_T1_T2_highres_final
) %>%
  mutate(structure = case_when(
    structure == "amyg" ~ "amygdala",
    TRUE                ~ "hippocampus"
  ))

# Amygdala — normalize by eTIV (mean eTIV = 1,948,106 mm³ used as reference)
combined_amygdala <- combined_all %>%
  filter(structure == "amygdala") %>%
  select(-structure) %>%
  na.omit() %>%
  merge(combined_aseg_final, by = "sa_id") %>%
  rename(etiv = size) %>%
  mutate(volume_size_norm = round((volume_size / etiv) * 1948106, 2))

cat("Unique amygdala subjects:", length(unique(combined_amygdala$sa_id)), "\n")

# Hippocampus — normalize by eTIV
combined_hippocampus <- combined_all %>%
  filter(structure == "hippocampus") %>%
  select(-structure) %>%
  na.omit() %>%
  merge(combined_aseg_final, by = "sa_id") %>%
  rename(etiv = size) %>%
  mutate(volume_size_norm = round((volume_size / etiv) * 1948106, 2))

cat("Unique hippocampus subjects:", length(unique(combined_hippocampus$sa_id)), "\n")

# ── Export Raw Combined Data ──────────────────────────────────────────────────
write_csv(combined_all,         file.path(output_dir, "combined_data_final.csv"))
write_csv(combined_amygdala,    file.path(output_dir, "combined_amygdala_data_final.csv"))
write_csv(combined_hippocampus, file.path(output_dir, "combined_hippocampus_data_final.csv"))

# ── Exploratory Violin Plots ──────────────────────────────────────────────────
method_colors <- c("green", "blue", "orange", "red", "purple")
method_labels <- c("T1", "T2 SPC", "T2 Highres", "T1+T2 SPC", "T1+T2 Highres")

base_theme <- theme(
  legend.box.background  = element_rect(),
  legend.background      = element_blank(),
  legend.key             = element_blank(),
  legend.position        = "inside",
  panel.background       = element_blank(),
  panel.grid.major       = element_blank(),
  panel.border           = element_blank(),
  axis.text.x            = element_text(angle = 345, vjust = 0.5),
  axis.line              = element_line(linewidth = .5, linetype = "solid"),
  plot.title             = element_text(size = 14, hjust = .5, face = "bold")
)

# Hippocampus violin
ggplot(combined_hippocampus, aes(x = subregion, y = volume_size_norm, fill = method)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = method_colors, labels = method_labels) +
  labs(title = "Violin Plots of Hippocampal Subregions Volume Size by Segmentation Method",
       x = "Subregion", y = "Normalized Volume Size", fill = "Segmentation Method") +
  base_theme +
  theme(legend.position.inside = c(.15, .85))

# Amygdala violin
ggplot(combined_amygdala, aes(x = subregion, y = volume_size_norm, fill = method)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = method_colors, labels = method_labels) +
  labs(title = "Violin Plots of Amygdala Subregions Volume Size by Segmentation Method",
       x = "Subregion", y = "Normalized Volume Size", fill = "Segmentation Method") +
  base_theme +
  theme(legend.position.inside = c(.15, .80))

# ── Load & Clean Demographic / Subtype Data ───────────────────────────────────
# Replace filenames below with your local exported data files.

stub <- read.csv(file.path(data_dir, "stub.csv"))
SA_subtype <- read.csv(file.path(data_dir, "results.csv"))

# Clean stub: keep key demographics and reshape race columns
stub_clean <- stub %>%
  select(
    PTID, SuperAging.ID, Date.of.Birth, Sex, Ethnicity,
    Race..choice.White.,
    Race..choice.Black.or.African.American.,
    Race..choice.American.Indian.or.Alaska.Native.,
    Race..choice.Native.Hawaiian.or.Other.Pacific.Islander.,
    Race..choice.Asian.,
    Race..choice.Other..specify..,
    Race..choice.Prefer.not.to.answer.
  ) %>%
  filter(SuperAging.ID != "") %>%
  rename(sa_id = SuperAging.ID) %>%
  pivot_longer(
    cols       = starts_with("Race.."),
    names_to   = "Race_option",
    values_to  = "Race_value"
  ) %>%
  filter(Race_value == "Checked") %>%
  select(-Race_value) %>%
  mutate(Race = case_when(
    Race_option == "Race..choice.White."                                        ~ "White",
    Race_option == "Race..choice.Black.or.African.American."                   ~ "Black or African American",
    Race_option == "Race..choice.American.Indian.or.Alaska.Native."            ~ "American Indian or Alaska Native",
    Race_option == "Race..choice.Native.Hawaiian.or.Other.Pacific.Islander."   ~ "Native Hawaiian or Other Pacific Islander",
    Race_option == "Race..choice.Asian."                                        ~ "Asian",
    Race_option == "Race..choice.Other..specify.."                             ~ "Other",
    Race_option == "Race..choice.Prefer.not.to.answer."                        ~ "Prefer not to answer",
    TRUE ~ NA_character_
  )) %>%
  select(-Race_option) %>%
  rename(DOB = Date.of.Birth) %>%
  mutate(Age = round(
    as.numeric(difftime(Sys.Date(), as.Date(DOB, format = "%Y-%m-%d"), units = "days")) / 365.25,
    2
  ))

# Clean subtype data
SA_subtype_clean <- SA_subtype %>%
  select(id, maintainer) %>%
  rename(PTID = id) %>%
  drop_na() %>%
  mutate(subtype = ifelse(maintainer == 1, "maintainer", "decliner"))

# Merge demographics with subtype
merged <- merge(stub_clean, SA_subtype_clean, by = "PTID")

# Align SA IDs: add leading zeros after "SANC" to match segmentation file IDs
merged$sa_id <- str_replace(merged$sa_id, "SANC", "SANC00")

# ── Harmonize Subject IDs ─────────────────────────────────────────────────────
# U-19 study IDs (format: NU013SA0###_m1) → standard SA IDs (format: SA###)

harmonize_ids <- function(df) {
  u19_ids <- df$sa_id[str_which(df$sa_id, "^NU013")]
  new_ids  <- paste0(substr(u19_ids, 7, 12) %>% sub("0", "", ., fixed = TRUE))
  df$sa_id[str_which(df$sa_id, "^NU013")] <- new_ids

  # Remove trailing visit letter (e.g., SA123a → SA123)
  df$sa_id <- if_else(
    str_detect(df$sa_id, "^SA[0-9]{3}[a-z]$"),
    str_sub(df$sa_id, 1, -2),
    df$sa_id
  )
  df
}

combined_hippocampus <- harmonize_ids(combined_hippocampus)
combined_amygdala    <- harmonize_ids(combined_amygdala)

# ── Merge Segmentation with Demographics ─────────────────────────────────────
merged_hippo_final <- merge(combined_hippocampus, merged, by = "sa_id")
cat("Unique hippocampus subjects (merged):", length(unique(merged_hippo_final$sa_id)), "\n")
writexl::write_xlsx(merged_hippo_final, file.path(output_dir, "merged_hipp_final.xlsx"))

merged_amyg_final <- merge(combined_amygdala, merged, by = "sa_id")
cat("Unique amygdala subjects (merged):", length(unique(merged_amyg_final$sa_id)), "\n")
writexl::write_xlsx(merged_amyg_final, file.path(output_dir, "merged_amyg_final.xlsx"))

# ── Boxplots by Method and Subtype ────────────────────────────────────────────
subtype_colors <- c("red", "blue")
subtype_breaks <- c("maintainer", "decliner")
subtype_labels <- c("Maintainer", "Decliner")

subtype_theme <- theme(
  legend.box.background  = element_rect(),
  legend.background      = element_blank(),
  legend.key             = element_blank(),
  legend.position        = "inside",
  legend.position.inside = c(.93, .20),
  panel.background       = element_blank(),
  panel.grid.major       = element_blank(),
  panel.border           = element_blank(),
  axis.line              = element_line(linewidth = .5, linetype = "solid"),
  plot.title             = element_text(size = 14, hjust = .5, face = "bold")
)

# Hippocampus — all methods
ggplot(merged_hippo_final, aes(x = volume_size_norm, y = subregion, fill = method)) +
  geom_boxplot() +
  scale_fill_manual(values = method_colors, labels = method_labels) +
  labs(title = "Boxplots of SA Hippocampal Subregion Volume Sizes & Subtypes by Segmentation Method",
       x = "Normalized Volume Size", y = "Subregion", fill = "Segmentation Method") +
  subtype_theme

# Hippocampus — per method, colored by subtype
for (m in c("T1", "T2", "T2H", "T1_T2", "T1_T2H")) {
  method_label <- switch(m,
    "T1"     = "T1 Only",
    "T2"     = "T2 Only",
    "T2H"    = "T2 Highres Only",
    "T1_T2"  = "T1 & T2",
    "T1_T2H" = "T1 & T2 Highres"
  )
  p <- ggplot(merged_hippo_final[merged_hippo_final$method == m, ],
              aes(x = subregion, y = volume_size_norm, fill = subtype)) +
    geom_boxplot() +
    coord_flip() +
    scale_fill_manual(values = subtype_colors,
                      breaks = subtype_breaks, labels = subtype_labels) +
    labs(
      title = paste("Boxplots of SA Hippocampal Subregion Normalized Volume Sizes by Subtype for", method_label),
      x = "Subregion", y = "Normalized Volume Size", fill = "SuperAging Subtype"
    ) +
    subtype_theme
  print(p)
}

# Amygdala — all methods
ggplot(merged_amyg_final, aes(x = volume_size_norm, y = subregion, fill = method)) +
  geom_boxplot() +
  scale_fill_manual(values = method_colors, labels = method_labels) +
  labs(title = "Boxplots of SA Amygdala Subregion Volume Sizes & Subtypes by Segmentation Method",
       x = "Normalized Volume Size", y = "Subregion + Subtype", fill = "Segmentation Method") +
  subtype_theme

# Amygdala — per method, colored by subtype
for (m in c("T1", "T2", "T2H", "T1_T2", "T1_T2H")) {
  method_label <- switch(m,
    "T1"     = "T1 Only",
    "T2"     = "T2 SPC",
    "T2H"    = "T2 Highres",
    "T1_T2"  = "T1 & T2 SPC",
    "T1_T2H" = "T1 & T2 Highres"
  )
  p <- ggplot(merged_amyg_final[merged_amyg_final$method == m, ],
              aes(x = subregion, y = volume_size_norm, fill = subtype)) +
    geom_boxplot() +
    coord_flip() +
    scale_fill_manual(values = subtype_colors,
                      breaks = subtype_breaks, labels = subtype_labels) +
    labs(
      title = paste("Boxplots of SA Amygdala Subregion Volume Sizes by Subtype for", method_label),
      x = "Subregion", y = "Normalized Volume Size", fill = "SuperAging Subtype"
    ) +
    subtype_theme +
    theme(legend.position.inside = c(.90, .20))
  print(p)
}

# ── Statistical Analysis ──────────────────────────────────────────────────────

# Remove hippocampal tail (unreliable segmentation boundary)
merged_hippo_final <- merged_hippo_final[merged_hippo_final$subregion != "Hippocampal_tail", ]

# Normality & variance checks
by(merged_hippo_final$volume_size_norm, merged_hippo_final$method, shapiro.test)
leveneTest(volume_size_norm ~ method, data = merged_hippo_final)

# Kruskal-Wallis test (non-parametric omnibus)
kruskal_result    <- kruskal.test(volume_size_norm ~ method, data = merged_hippo_final)
print(kruskal_result)
kruskal_result_df <- do.call(cbind, kruskal_result)
write.csv(kruskal_result_df, file.path(output_dir, "kruskal_results.csv"))

# Pairwise Wilcoxon signed-rank tests with Bonferroni correction
methods <- c("T1", "T2", "T2H", "T1_T2", "T1_T2H")
wilcoxon_results <- list()

for (i in 1:(length(methods) - 1)) {
  for (j in (i + 1):length(methods)) {
    group1      <- merged_hippo_final$volume_size_norm[merged_hippo_final$method == methods[i]]
    group2      <- merged_hippo_final$volume_size_norm[merged_hippo_final$method == methods[j]]
    test_result <- wilcox.test(group1, group2, paired = TRUE, exact = FALSE)
    wilcoxon_results[[paste(methods[i], methods[j], sep = " vs. ")]] <- data.frame(
      Method1 = methods[i],
      Method2 = methods[j],
      W       = test_result$statistic,
      p_value = format(test_result$p.value, scientific = FALSE)
    )
  }
}

wilcoxon_results_df <- do.call(rbind, wilcoxon_results)
wilcoxon_results_df$adjusted_p_value <- format(
  p.adjust(as.numeric(wilcoxon_results_df$p_value), method = "bonferroni"),
  scientific = FALSE
)

print(wilcoxon_results_df)
write_csv(wilcoxon_results_df, file.path(output_dir, "wilcoxon_results.csv"))

# Formatted gt table of pairwise results
wilcoxon_poster_table <- wilcoxon_results_df %>%
  mutate(
    p_value_num          = as.numeric(p_value),
    adjusted_p_value_num = as.numeric(adjusted_p_value),
    significance = case_when(
      adjusted_p_value_num < 0.001 ~ "***",
      adjusted_p_value_num < 0.01  ~ "**",
      adjusted_p_value_num < 0.05  ~ "*",
      TRUE                          ~ ""
    ),
    `Unadjusted p` = ifelse(p_value_num < 0.0001, "< 0.0001***",
                            sprintf("%.4f%s", p_value_num, significance)),
    `Adjusted p`   = ifelse(adjusted_p_value_num < 0.0001, "< 0.0001***",
                            sprintf("%.4f%s", adjusted_p_value_num, significance)),
    `W Statistic`  = format(W, big.mark = ","),
    Comparison     = paste(Method1, "vs.", Method2) %>%
      str_replace("T1_T2H?", function(x) if (x == "T1_T2H") "T1 & T2H" else "T1 & T2")
  ) %>%
  select(Comparison, `W Statistic`, `Unadjusted p`, `Adjusted p`)

poster_table <- wilcoxon_poster_table %>%
  gt() %>%
  tab_header(
    title    = md("**Table 1. Pairwise Wilcoxon Signed-Rank Tests**"),
    subtitle = "Bonferroni-corrected comparisons of hippocampal subfield volumes across methods"
  ) %>%
  cols_align(align = "center", columns = c(`W Statistic`, `Unadjusted p`, `Adjusted p`)) %>%
  cols_align(align = "left",   columns = Comparison) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(rows = grepl("\\*", `Adjusted p`))
  ) %>%
  tab_style(
    style     = cell_fill(color = "white"),
    locations = cells_body(columns = everything(), rows = everything())
  ) %>%
  tab_footnote(
    footnote  = "*** p < 0.001; ** p < 0.01; * p < 0.05",
    locations = cells_column_labels(columns = `Adjusted p`)
  ) %>%
  tab_options(
    table.font.size            = 18,
    column_labels.font.weight  = "bold",
    table.width                = "100%",
    row_group.font.weight      = "bold"
  )

gtsave(poster_table,
       filename = file.path(plot_dir, "Table1_Wilcoxon_Results.png"),
       zoom = 2, expand = 10)

# ── Linear Model: Volume ~ Method * Subtype + Subregion ──────────────────────
model   <- lm(volume_size_norm ~ method * subtype + subregion, data = merged_hippo_final)
results <- summary(model)
print(results)

anova_results <- anova(model)
print(anova_results)

# Post-hoc comparisons
emm <- emmeans(model, ~ method | subtype, adjust = "bonferroni")
summary(emm, infer = TRUE)
pairs(emm, adjust = "bonferroni")

subtype_emm <- emmeans(model, ~ subtype)
pairs(subtype_emm, adjust = "bonferroni")
summary(subtype_emm, infer = TRUE)

# ── Figure 1: Estimated Marginal Means Plot ───────────────────────────────────
emm_df <- as.data.frame(emm)

sig_data <- data.frame(
  x       = 1,
  xend    = 3,
  y       = 560,
  t_value = "t = -.932",
  p_label = "p = .321"
)

figure1 <- ggplot(emm_df, aes(x = method, y = emmean, fill = subtype)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(0.8), width = 0.25,
                color = "black", linewidth = 0.8) +
  geom_segment(data = sig_data,
               aes(x = x, xend = xend, y = y, yend = y),
               inherit.aes = FALSE) +
  geom_text(data = sig_data,
            aes(x = (x + xend) / 1.6, y = y + 17, label = p_label),
            inherit.aes = FALSE, size = 10) +
  geom_text(data = sig_data,
            aes(x = (x + xend) / 2, y = y + 17, label = t_value),
            inherit.aes = FALSE, size = 10) +
  geom_point(aes(group = subtype), position = position_dodge(0.8),
             size = 2, color = "black", show.legend = FALSE) +
  labs(
    title    = "Figure 1. Hippocampal Volume by Segmentation Method and Subtype",
    subtitle = "Error bars show 95% CIs; All methods significantly different to T1 besides combined T1 & T2H",
    y        = "Adjusted Volume (normalized)",
    x        = "Segmentation Method",
    fill     = "Subtype"
  ) +
  theme_linedraw() +
  scale_fill_brewer(palette = "Dark2", labels = c("Decliner", "Maintainer")) +
  scale_x_discrete(labels = c("T1", "T1 & T2", "T1 & T2H", "T2", "T2H")) +
  theme(
    plot.title.position = "plot",
    plot.title          = element_text(size = 28, face = "bold", hjust = 0),
    plot.subtitle       = element_text(size = 20, hjust = 0),
    axis.title          = element_text(face = "bold", size = 24),
    axis.text.x         = element_text(size = 22),
    axis.text.y         = element_text(size = 22),
    legend.position     = "top",
    legend.box.background = element_rect(color = "black"),
    legend.title        = element_text(face = "bold", size = 20),
    legend.text         = element_text(size = 20)
  )

print(figure1)
ggsave(file.path(plot_dir, "poster_figure1.png"), figure1,
       width = 18, height = 8, units = "in", dpi = 300)

# ── Export Linear Model Summary to Excel ─────────────────────────────────────
coefficients   <- results$coefficients
r_squared      <- results$r.squared
adj_r_squared  <- results$adj.r.squared
f_statistic    <- results$fstatistic
p_value        <- format.pval(
  pf(f_statistic[1], f_statistic[2], f_statistic[3], lower.tail = FALSE),
  digits = 3
)
model_formula  <- paste0(
  results[["terms"]][[2]], " ",
  results[["terms"]][[1]], " ",
  results[["terms"]][[3]]
)

wb <- createWorkbook()
addWorksheet(wb, "Model Summary")

writeData(wb, "Model Summary", "Coefficients",  startRow = 1, startCol = 1)
writeData(wb, "Model Summary", coefficients,    startRow = 2, startCol = 1, rowNames = TRUE)

nr <- nrow(coefficients)
writeData(wb, "Model Summary", "R-Squared",          startRow = nr + 4, startCol = 1)
writeData(wb, "Model Summary", r_squared,             startRow = nr + 4, startCol = 2)
writeData(wb, "Model Summary", "Adjusted R-Squared",  startRow = nr + 5, startCol = 1)
writeData(wb, "Model Summary", adj_r_squared,         startRow = nr + 5, startCol = 2)
writeData(wb, "Model Summary", "F-Statistic",         startRow = nr + 6, startCol = 1)
writeData(wb, "Model Summary", f_statistic[1],        startRow = nr + 6, startCol = 2)
writeData(wb, "Model Summary", "p-Value",             startRow = nr + 7, startCol = 1)
writeData(wb, "Model Summary", p_value,               startRow = nr + 7, startCol = 2)
writeData(wb, "Model Summary", "Model Formula",       startRow = nr + 8, startCol = 1)
writeData(wb, "Model Summary", model_formula,         startRow = nr + 8, startCol = 2)

saveWorkbook(wb, file.path(output_dir, "lm_model_summary.xlsx"), overwrite = TRUE)
