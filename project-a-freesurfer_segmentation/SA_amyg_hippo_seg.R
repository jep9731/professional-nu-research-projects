# SA FS Amyg/Hippo Project

# Load libraries ----------------------------------------------------------
library(readr)
library(tidyverse)
library(writexl)
#install.packages("freesurfer")
library(freesurfer)

# Create for loop (T1) ----------------------------------------------------
## Find files
files_T1 <- list.files(path = "/Volumes/fsmresfiles/CNADC/Imaging_Core/Imaging/imaging_projects/Individuals/Joshua_Pasaye/FS_highreshippo/sa_fs_seg_final/data_final/", 
                    pattern = "*-T1.v21.txt$", 
                    recursive = TRUE, 
                    full.names = TRUE)
files_T1

## Create empty lists
data_list_T1 <- list()
data_list_T1_id <- list()
data_list_T1_hem <- list()
data_list_T1_struc <- list()
data_list_T1_method <- list()

## Iterates over files to read them
for (file in files_T1) {
    data_T1 <- read.table(file)
    data_list_T1[[file]] <- data_T1
    subject_id <- basename(dirname(file))
    data_list_T1_id[[file]] <- as.data.frame(rep(subject_id, nrow(data_T1)))
    hemisphere <- regmatches(file, regexpr("(rh|lh)", file))
    data_list_T1_hem[[file]] <- as.data.frame(rep(hemisphere, nrow(data_T1)))
    structure <- sub(".*[lr]h\\.([a-z]+{1,5}).*-T1.*", "\\1", file)
    data_list_T1_struc[[file]] <- as.data.frame(rep(structure, nrow(data_T1)))
    method <- regmatches(file, regexpr("T1", file))
    data_list_T1_method[[file]] <- as.data.frame(rep(method, nrow(data_T1)))
}

## Combine files
combined_data_T1 <- bind_rows(data_list_T1)
combined_data_T1_id <- bind_rows(data_list_T1_id)
combined_data_T1_hem <- bind_rows(data_list_T1_hem)
combined_data_T1_struc <- bind_rows(data_list_T1_struc)
combined_data_T1_method <- bind_rows(data_list_T1_method)
combined_T1 <- cbind(combined_data_T1_id, 
                     combined_data_T1_hem,
                     combined_data_T1_struc,
                     combined_data_T1_method,
                     combined_data_T1)

## Rename files columns
old_col_names_T1 <- names(combined_T1)
new_col_names_T1 <- c("sa_id", "hemisphere", "structure", "method", "subregion", "volume_size")
combined_T1_final <- combined_T1 %>% rename_with(~ new_col_names_T1, 
                                                 all_of(old_col_names_T1))
combined_T1_final

# Create for loop (T2 SPC) ------------------------------------------------
## Find files
files_T2_spc <- list.files(path = "/Volumes/fsmresfiles/CNADC/Imaging_Core/Imaging/imaging_projects/Individuals/Joshua_Pasaye/FS_highreshippo/sa_fs_seg_final/data_final/", 
                       pattern = "*_T2spcnorm.v21.txt$", 
                       recursive = TRUE, 
                       full.names = TRUE)

files_T2_spc <- files_T2_spc[!grepl("-T1-", files_T2_spc)]
files_T2_spc

## Create empty list
data_list_T2_spc <- list()
data_list_T2_spc_id <- list()
data_list_T2_spc_hem <- list()
data_list_T2_spc_struc <- list()
data_list_T2_spc_method <- list()

## Iterates over files to read them
for (file in files_T2_spc) {
  data_T2_spc <- read.table(file)
  data_list_T2_spc[[file]] <- data_T2_spc
  subject_id <- basename(dirname(file))
  data_list_T2_spc_id[[file]] <- as.data.frame(rep(subject_id, nrow(data_T2_spc)))
  hemisphere <- regmatches(file, regexpr("(rh|lh)", file))
  data_list_T2_spc_hem[[file]] <- as.data.frame(rep(hemisphere, nrow(data_T2_spc)))
  structure <- sub(".*[lr]h\\.([a-z]+{1,5}).*", "\\1", file)
  data_list_T2_spc_struc[[file]] <- as.data.frame(rep(structure, nrow(data_T2_spc)))
  method <- "T2"
  data_list_T2_spc_method[[file]] <- as.data.frame(rep(method, nrow(data_T2_spc)))
}

## Combine files
combined_data_T2_spc <- bind_rows(data_list_T2_spc)
combined_data_T2_spc_id <- bind_rows(data_list_T2_spc_id)
combined_data_T2_spc_hem <- bind_rows(data_list_T2_spc_hem)
combined_data_T2_spc_struc <- bind_rows(data_list_T2_spc_struc)
combined_data_T2_spc_method <- bind_rows(data_list_T2_spc_method)
combined_T2_spc <- cbind(combined_data_T2_spc_id, 
                         combined_data_T2_spc_hem, 
                         combined_data_T2_spc_struc,
                         combined_data_T2_spc_method,
                         combined_data_T2_spc)

## Rename files columns
old_col_names_T2_spc <- names(combined_T2_spc)
new_col_names_T2_spc <- c("sa_id", "hemisphere", "structure", "method", "subregion", "volume_size")
combined_T2_spc_final <- combined_T2_spc %>% rename_with(~ new_col_names_T2_spc, 
                                                     all_of(old_col_names_T2_spc))
combined_T2_spc_final

# Create for loop (T2 Highres) --------------------------------------------
## Find files
files_T2_highres <- list.files(path = "/Volumes/fsmresfiles/CNADC/Imaging_Core/Imaging/imaging_projects/Individuals/Joshua_Pasaye/FS_highreshippo/sa_fs_seg_final/data_final/", 
                           pattern = "*_T2highreshipp.v21.txt$", 
                           recursive = TRUE, 
                           full.names = TRUE)

files_T2_highres <- files_T2_highres[!grepl("-T1-", files_T2_highres)]
files_T2_highres

## Create empty list
data_list_T2_highres <- list()
data_list_T2_highres_id <- list()
data_list_T2_highres_hem <- list()
data_list_T2_highres_struc <- list()
data_list_T2_highres_method <- list()

## Iterates over files to read them
for (file in files_T2_highres) {
  data_T2_highres <- read.table(file)
  data_list_T2_highres[[file]] <- data_T2_highres
  subject_id <- basename(dirname(file))
  data_list_T2_highres_id[[file]] <- as.data.frame(rep(subject_id, nrow(data_T2_highres)))
  hemisphere <- regmatches(file, regexpr("(rh|lh)", file))
  data_list_T2_highres_hem[[file]] <- as.data.frame(rep(hemisphere, nrow(data_T2_highres)))
  structure <- sub(".*[lr]h\\.([a-z]+{1,5}).*", "\\1", file)
  data_list_T2_highres_struc[[file]] <- as.data.frame(rep(structure, nrow(data_T2_highres)))
  method <- "T2H"
  data_list_T2_highres_method[[file]] <- as.data.frame(rep(method, nrow(data_T2_highres)))
}

## Combine files
combined_data_T2_highres <- bind_rows(data_list_T2_highres)
combined_data_T2_highres_id <- bind_rows(data_list_T2_highres_id)
combined_data_T2_highres_hem <- bind_rows(data_list_T2_highres_hem)
combined_data_T2_highres_struc <- bind_rows(data_list_T2_highres_struc)
combined_data_T2_highres_method <- bind_rows(data_list_T2_highres_method)
combined_T2_highres <- cbind(combined_data_T2_highres_id, 
                         combined_data_T2_highres_hem, 
                         combined_data_T2_highres_struc,
                         combined_data_T2_highres_method,
                         combined_data_T2_highres)

## Rename file columns
old_col_names_T2_highres <- names(combined_T2_highres)
new_col_names_T2_highres <- c("sa_id", "hemisphere", "structure", "method", "subregion", "volume_size")
combined_T2_highres_final <- combined_T2_highres %>% rename_with(~ new_col_names_T2_highres, 
                                                         all_of(old_col_names_T2_highres))
combined_T2_highres_final

# Create for loop (T1+T2 SPC) ------------------------------------------------
## Find files
files_T1_T2_spc <- list.files(path = "/Volumes/fsmresfiles/CNADC/Imaging_Core/Imaging/imaging_projects/Individuals/Joshua_Pasaye/FS_highreshippo/sa_fs_seg_final/data_final/", 
                           pattern = "*-T1.*_T2spcnorm.v21.txt$", 
                           recursive = TRUE, 
                           full.names = TRUE)

files_T1_T2_spc <- files_T1_T2_spc[grepl("-T1-", files_T1_T2_spc)]
files_T1_T2_spc

## Create empty list
data_list_T1_T2_spc <- list()
data_list_T1_T2_spc_id <- list()
data_list_T1_T2_spc_hem <- list()
data_list_T1_T2_spc_struc <- list()
data_list_T1_T2_spc_method <- list()

## Iterates over files to read them
for (file in files_T1_T2_spc) {
  data_T1_T2_spc <- read.table(file)
  data_list_T1_T2_spc[[file]] <- data_T1_T2_spc
  subject_id <- basename(dirname(file))
  data_list_T1_T2_spc_id[[file]] <- as.data.frame(rep(subject_id, nrow(data_T1_T2_spc)))
  hemisphere <- regmatches(file, regexpr("(rh|lh)", file))
  data_list_T1_T2_spc_hem[[file]] <- as.data.frame(rep(hemisphere, nrow(data_T1_T2_spc)))
  structure <- sub(".*[lr]h\\.([a-z]+{1,5}).*", "\\1", file)
  data_list_T1_T2_spc_struc[[file]] <- as.data.frame(rep(structure, nrow(data_T1_T2_spc)))
  method <- "T1_T2"
  data_list_T1_T2_spc_method[[file]] <- as.data.frame(rep(method, nrow(data_T1_T2_spc)))
}

## Combine files
combined_data_T1_T2_spc <- bind_rows(data_list_T1_T2_spc)
combined_data_T1_T2_spc_id <- bind_rows(data_list_T1_T2_spc_id)
combined_data_T1_T2_spc_hem <- bind_rows(data_list_T1_T2_spc_hem)
combined_data_T1_T2_spc_struc <- bind_rows(data_list_T1_T2_spc_struc)
combined_data_T1_T2_spc_method <- bind_rows(data_list_T1_T2_spc_method)
combined_T1_T2_spc <- cbind(combined_data_T1_T2_spc_id, 
                         combined_data_T1_T2_spc_hem, 
                         combined_data_T1_T2_spc_struc,
                         combined_data_T1_T2_spc_method,
                         combined_data_T1_T2_spc)

## Rename files columns
old_col_names_T1_T2_spc <- names(combined_T1_T2_spc)
new_col_names_T1_T2_spc <- c("sa_id", "hemisphere", "structure", "method", "subregion", "volume_size")
combined_T1_T2_spc_final <- combined_T1_T2_spc %>% rename_with(~ new_col_names_T1_T2_spc, 
                                                         all_of(old_col_names_T1_T2_spc))
combined_T1_T2_spc_final

# Create for loop (T1_T2 Highres) --------------------------------------------
## Find files
files_T1_T2_highres <- list.files(path = "/Volumes/fsmresfiles/CNADC/Imaging_Core/Imaging/imaging_projects/Individuals/Joshua_Pasaye/FS_highreshippo/sa_fs_seg_final/data_final/", 
                               pattern = "*-T1.*_T2highreshipp.v21.txt$", 
                               recursive = TRUE, 
                               full.names = TRUE)

files_T1_T2_highres <- files_T1_T2_highres[grepl("-T1-", files_T1_T2_highres)]
files_T1_T2_highres

## Create empty list
data_list_T1_T2_highres <- list()
data_list_T1_T2_highres_id <- list()
data_list_T1_T2_highres_hem <- list()
data_list_T1_T2_highres_struc <- list()
data_list_T1_T2_highres_method <- list()

## Iterates over files to read them
for (file in files_T1_T2_highres) {
  data_T1_T2_highres <- read.table(file)
  data_list_T1_T2_highres[[file]] <- data_T1_T2_highres
  subject_id <- basename(dirname(file))
  data_list_T1_T2_highres_id[[file]] <- as.data.frame(rep(subject_id, nrow(data_T1_T2_highres)))
  hemisphere <- regmatches(file, regexpr("(rh|lh)", file))
  data_list_T1_T2_highres_hem[[file]] <- as.data.frame(rep(hemisphere, nrow(data_T1_T2_highres)))
  structure <- sub(".*[lr]h\\.([a-z]+{1,5}).*", "\\1", file)
  data_list_T1_T2_highres_struc[[file]] <- as.data.frame(rep(structure, nrow(data_T1_T2_highres)))
  method <- "T1_T2H"
  data_list_T1_T2_highres_method[[file]] <- as.data.frame(rep(method, nrow(data_T1_T2_highres)))
}

## Combine files
combined_data_T1_T2_highres <- bind_rows(data_list_T1_T2_highres)
combined_data_T1_T2_highres_id <- bind_rows(data_list_T1_T2_highres_id)
combined_data_T1_T2_highres_hem <- bind_rows(data_list_T1_T2_highres_hem)
combined_data_T1_T2_highres_struc <- bind_rows(data_list_T1_T2_highres_struc)
combined_data_T1_T2_highres_method <- bind_rows(data_list_T1_T2_highres_method)
combined_T1_T2_highres <- cbind(combined_data_T1_T2_highres_id, 
                             combined_data_T1_T2_highres_hem, 
                             combined_data_T1_T2_highres_struc,
                             combined_data_T1_T2_highres_method,
                             combined_data_T1_T2_highres)

## Rename file columns
old_col_names_T1_T2_highres <- names(combined_T1_T2_highres)
new_col_names_T1_T2_highres <- c("sa_id", "hemisphere", "structure", "method", "subregion", "volume_size")
combined_T1_T2_highres_final <- combined_T1_T2_highres %>% rename_with(~ new_col_names_T1_T2_highres, 
                                                                 all_of(old_col_names_T1_T2_highres))
combined_T1_T2_highres_final

# Get aseg files ----------------------------------------------------------
aseg <- list.files(path = "/Volumes/fsmresfiles/CNADC/Imaging_Core/Imaging/imaging_projects/Individuals/Joshua_Pasaye/FS_highreshippo/sa_hippseg",
                      pattern = "aseg.stats",
                      recursive = TRUE,
                      full.names = TRUE)

aseg

aseg_list <- list()
aseg_list_id <- list()

for (file in aseg) {
  aseg_data <- read_aseg_stats(file)
  aseg_list[[file]] <- aseg_data[["measures"]]
  aseg_list_id[[file]] <-  as.data.frame(rep(strsplit(strsplit(file,"/stats")[[1]][1],"sa_hippseg/")[[1]][2], 
                                                        nrow(aseg_data[["measures"]])))
}

aseg_files <- bind_rows(aseg_list)
aseg_files_id <- bind_rows(aseg_list_id)

combined_aseg <- cbind(aseg_files_id,
                       aseg_files)

old_col_names_aseg <- names(combined_aseg)
new_col_names_aseg <- c("sa_id", "measure", "measure_long", "meaning", "size", "units")

combined_aseg_final <- combined_aseg %>% 
  rename_with(~ new_col_names_aseg,
              all_of(old_col_names_aseg)) %>%
  filter(measure_long == "etiv") %>%
  subset(., select = c("sa_id", "size"))
  
combined_aseg_final$size = as.numeric(combined_aseg_final$size)

# Combined all files ------------------------------------------------------
combined_all <- bind_rows(
  combined_T1_final,
  combined_T2_spc_final,
  combined_T2_highres_final,
  combined_T1_T2_spc_final,
  combined_T1_T2_highres_final) %>%
  mutate(structure = case_when(structure == "amyg" ~ "amygdala",
                                            TRUE ~ "hippocampus"))

combined_all

## Amygdala data
combined_amygdala <- combined_all %>%
  filter(structure == "amygdala") %>%
  select(!structure) %>%
  na.omit()
combined_amygdala <- merge(combined_amygdala, combined_aseg_final, by = "sa_id") %>% rename(etiv = size)
combined_amygdala$volume_size_norm = round((combined_amygdala$volume_size/combined_amygdala$etiv) * 1948106, 2)
combined_amygdala
length(unique(combined_amygdala$sa_id))

## Hippocampus data
combined_hippocampus <- combined_all %>%
  filter(structure == "hippocampus") %>%
  select(!structure) %>%
  na.omit()
combined_hippocampus <- merge(combined_hippocampus, combined_aseg_final, by = "sa_id") %>% rename(etiv = size)
combined_hippocampus$volume_size_norm = round((combined_hippocampus$volume_size/combined_hippocampus$etiv) * 1948106, 2)
combined_hippocampus
length(unique(combined_hippocampus$sa_id))

# Write out files ---------------------------------------------------------
write_csv(combined_all, "/Volumes/fsmresfiles/CNADC/Imaging_Core/Imaging/imaging_projects/Individuals/Joshua_Pasaye/FS_highreshippo/outputs/combined_data_final.csv")
write_csv(combined_amygdala, "/Volumes/fsmresfiles/CNADC/Imaging_Core/Imaging/imaging_projects/Individuals/Joshua_Pasaye/FS_highreshippo/outputs/combined_amygdala_data_final.csv")
write_csv(combined_hippocampus, "/Volumes/fsmresfiles/CNADC/Imaging_Core/Imaging/imaging_projects/Individuals/Joshua_Pasaye/FS_highreshippo/outputs/combined_hippocampus_data_final.csv")

# Graphs ------------------------------------------------------------------
ggplot(data = combined_hippocampus, aes(x = subregion, y = volume_size_norm, fill = method)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("green", "blue", "orange", "red", "purple"), 
                    labels = c("T1", "T2 SPC", "T2 Highres", "T1+T2 SPC", "T1+T2 Highres")) +
  labs(title = "Violin Plots of Hippocampal Subregions Volume Size by Segmentation Method",
       x = "Subregion",
       y = "Normalized Volume Size",
       fill = "Segmentation Method") +
  theme(legend.box.background = element_rect(), # Add legend box
              legend.background = element_blank(), # Remove legend background
              legend.key = element_blank(), # Remove box around colors
              legend.position = "inside", # Move legend inside graph
              legend.position.inside = c(.15, .85), # Change position of legend
              panel.background = element_blank(), # Remove plot panel background
              panel.grid.major = element_blank(), # Remove major grid lines
              panel.border = element_blank(), # Remove panel border
              axis.text.x = element_text(angle = 345, vjust = 0.5), # Rotate x axis label
              axis.line = element_line(linewidth = .5, linetype = "solid"), # Change y axis line size
              plot.title = element_text(size = 14, hjust = .5, face = "bold") # Change title size and position
  )

## Aymgdala
ggplot(data = combined_amygdala, aes(x = subregion, y = volume_size_norm, fill = method)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("green", "blue", "orange", "red", "purple"), 
                    labels = c("T1", "T2 SPC", "T2 Highres", "T1+T2 SPC", "T1+T2 Highres")) +
  labs(title = "Violin Plots of Amygdala Subregions Volume Size by Segmentation Method",
       x = "Subregion",
       y = "Normalized Volume Size",
       fill = "Segmentation Method") +
  theme(legend.box.background = element_rect(), # Add legend box
        legend.background = element_blank(), # Remove legend background
        legend.key = element_blank(), # Remove box around colors
        legend.position = "inside", # Move legend inside graph
        legend.position.inside = c(.15, .80), # Change position of legend
        panel.background = element_blank(), # Remove plot panel background
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.border = element_blank(), # Remove panel border
        axis.text.x = element_text(angle = 345, vjust = 0.5), # Rotate x axis label
        axis.line = element_line(linewidth = .5, linetype = "solid"), # Change axis line size
        plot.title = element_text(size = 14, hjust = .5, face = "bold") # Change title size and position
  )

# Analysis Setup ----------------------------------------------------------------
# Load data
stub <- read.csv("/Volumes/fsmresfiles/CNADC/Imaging_Core/Imaging/imaging_projects/Individuals/Joshua_Pasaye/FS_highreshippo/data/AlzheimersDiseaseRes-Stub_DATA_LABELS_2024-11-15_2259.csv")
SA_subtype <- read.csv("/Volumes/fsmresfiles/CNADC/Imaging_Core/Imaging/imaging_projects/Individuals/Joshua_Pasaye/FS_highreshippo/data/results.csv")

# Clean stub data
stub_clean <- stub %>%
  select(PTID, SuperAging.ID, Date.of.Birth, Sex, Ethnicity,
         Race..choice.White., Race..choice.Black.or.African.American., Race..choice.American.Indian.or.Alaska.Native.,
         Race..choice.Native.Hawaiian.or.Other.Pacific.Islander., Race..choice.Asian., Race..choice.Other..specify..,
         Race..choice.Prefer.not.to.answer.) %>% # Only select PTID and SA ID
  filter(SuperAging.ID != "") %>% # Removes blank values
  rename(., sa_id = SuperAging.ID) %>% # Change SA ID column name
  pivot_longer(cols = starts_with("Race.."),
               names_to = "Race_option",
               values_to = "Race_value"
               ) %>% # Make longer 'Race' column
  filter(Race_value == "Checked") %>%
  select(-Race_value) %>% 
  mutate(Race = case_when(
    Race_option == "Race..choice.White." ~ "White",
    Race_option == "Race..choice.Black.or.African.American." ~ "Black or African American",
    Race_option == "Race..choice.American.Indian.or.Alaska.Native." ~ "American Indian or Alaska Native",
    Race_option == "Race..choice.Native.Hawaiian.or.Other.Pacific.Islander." ~ "Native Hawaiian or Other Pacific Islander",
    Race_option == "Race..choice.Asian." ~ "Asian",
    Race_option == "Race..choice.Other..specify.." ~ "Other",
    Race_option == "Race..choice.Prefer.not.to.answer." ~ "Prefer not to answer",
    TRUE ~ NA) # Clean up the race labels
  ) %>%
  select(-Race_option) %>%
  rename(DOB = Date.of.Birth) %>%
  mutate(Age = round(as.numeric(difftime(as.Date(Sys.Date(), format="%Y-%m-%d"), 
                                   as.Date(DOB, format="%Y-%m-%d")
                                   ,unit="days")
                                ) / 365.25, 2)
         )
stub_clean

# Clean subtype data
SA_subtype_clean <- SA_subtype %>%
  select(id, maintainer) %>% # Keep only id and maintainer status
  rename(., PTID = id) %>% # Change id column name to match PTID on stub
  drop_na(.) %>% # Remove rows with NA values
  mutate(subtype = ifelse(maintainer == 1, "maintainer", "decliner")) # Add column with string maintainer and decliner
SA_subtype_clean

# Merge datasets
merged <- merge(stub_clean, SA_subtype_clean, by = "PTID")
merged

# Add 0's after SANC for merged data database to match amyg/hipp ids
merged$sa_id <- str_replace(merged$sa_id, "SANC", "SANC00")
merged

# Rename U-19 IDs to regular SA IDs for hippocampus
U19_ids_hipp <- combined_hippocampus$sa_id[str_which(combined_hippocampus$sa_id, "^NU013")] # Pulls out U19 IDs
hipp_sa_ids_new <- paste0(substr(U19_ids_hipp, 7, 12) %>% # creates substrings after NU013 and before _m1 (i.e., only SA0### part)
                            sub("0", "", ., fixed = TRUE)) # removes leading zero and pastes two pieces together
for (i in combined_hippocampus$sa_id) { # Created for loop for each value in the sa_id column
  combined_hippocampus$sa_id[str_which(combined_hippocampus$sa_id, "^NU013")] <- hipp_sa_ids_new # pulled U19 IDs and replaced with regular ids
}
combined_hippocampus

# Rename U-19 IDs to regular SA IDs for amygdala
U19_ids_amyg <- combined_amygdala$sa_id[str_which(combined_amygdala$sa_id, "^NU013")] # Pulls on U19 IDs
amyg_sa_ids_new <- paste0(substr(U19_ids_amyg, 7, 12) %>% # creates substrings after NU013 and before _m1 (i.e., only SA0### part)
                         sub("0", "", ., fixed = TRUE)) # removes leading zero and pastes two pieces together
for (i in combined_amygdala$sa_id) { # Created for loop for each value in the sa_id column
  combined_amygdala$sa_id[str_which(combined_amygdala$sa_id, "^NU013")] <- amyg_sa_ids_new # pulled U19 IDs and replaced with regular ids
}
combined_amygdala

# Remove visit letter from SA IDs for hippocampus
combined_hippocampus$sa_id <- if_else(str_detect(combined_hippocampus$sa_id, "^SA[0-9]{3}[a-z]$"), 
                                      str_sub(combined_hippocampus$sa_id, 1, -2), 
                                      combined_hippocampus$sa_id)
combined_hippocampus

# Remove visit letter from SA IDs for amygdala
combined_amygdala$sa_id <- if_else(str_detect(combined_amygdala$sa_id, "^SA[0-9]{3}[a-z]$"),
                                   str_sub(combined_amygdala$sa_id, 1, -2), 
                                   combined_amygdala$sa_id)
combined_amygdala

# Merge subtype database with hippo/amyg databases
merged_hippo_final <- merge(combined_hippocampus, merged, by = "sa_id")
merged_hippo_final
writexl::write_xlsx(merged_hippo_final, "/Volumes/fsmresfiles/CNADC/Imaging_Core/Imaging/imaging_projects/Individuals/Joshua_Pasaye/FS_highreshippo/outputs/merged_hipp_final.xlsx")
length(unique(merged_hippo_final$sa_id))

merged_amyg_final <- merge(combined_amygdala, merged, by = "sa_id")
merged_amyg_final
writexl::write_xlsx(merged_amyg_final, "/Volumes/fsmresfiles/CNADC/Imaging_Core/Imaging/imaging_projects/Individuals/Joshua_Pasaye/FS_highreshippo/outputs/merged_amyg_final.xlsx")
length(unique(merged_amyg_final))

# Graphs -----------------------------------------------------------
# Hippocampus plots
ggplot(data = merged_hippo_final, aes(x = volume_size_norm, y = subregion, fill = method)) +
  geom_boxplot() +
  labs(title = "Boxplots of SA Hippocampal Subregion Volume Sizes & Subtypes by Segmentation Method",
       x = "Subregion",
       y = "Normalized Volume Size",
       fill = "Segmentation Method") +
  scale_fill_manual(values = c("green", "blue", "orange", "red", "purple"),
    labels = c("T1", "T2 SPC", "T2 Highres", "T1+T2 SPC", "T1+T2 Highres")) +
  theme(legend.box.background = element_rect(), # Add legend box
        legend.background = element_blank(), # Remove legend background
        legend.key = element_blank(), # Remove box around colors
        legend.position = "inside", # Move legend inside graph
        legend.position.inside = c(.93, .20), # Change position of legend
        panel.background = element_blank(), # Remove plot panel background
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.border = element_blank(), # Remove panel border
        axis.line = element_line(linewidth = .5, linetype = "solid"), # Change y axis line size
        plot.title = element_text(size = 14, hjust = .5, face = "bold") # Change title size and position
  )

ggplot(data = merged_hippo_final[merged_hippo_final$method == "T1", ], 
       aes(x = subregion, y = volume_size_norm, fill = subtype)) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_manual(values = c("red",  "blue"), 
                    breaks = c("maintainer", "decliner"),
                    labels = c("Maintainer", "Decliner")) +
  labs(title = "Boxplots of SA Hippocampal Subregion Normalized Volume Sizes by Subtype for T1 Only Method",
       x = "Subregion",
       y = "Normalized Volume Size",
       fill = "SuperAging Subtype") +
  theme(legend.box.background = element_rect(), # Add legend box
        legend.background = element_blank(), # Remove legend background
        legend.key = element_blank(), # Remove box around colors
        legend.position = "inside", # Move legend inside graph
        legend.position.inside = c(.93, .20), # Change position of legend
        panel.background = element_blank(), # Remove plot panel background
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.border = element_blank(), # Remove panel border
        axis.line = element_line(linewidth = .5, linetype = "solid"), # Change y axis line size
        plot.title = element_text(size = 14, hjust = .5, face = "bold") # Change title size and position
  )

ggplot(data = merged_hippo_final[merged_hippo_final$method == "T2", ], 
       aes(x = subregion, y = volume_size_norm, fill = subtype)) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_manual(values = c("red",  "blue"), 
                    breaks = c("maintainer", "decliner"),
                    labels = c("Maintainer", "Decliner")) +
  labs(title = "Boxplots of SA Hippocampal Subregion Normalized Volume Size by Subtype for T2 Only Method",
       x = "Subregion",
       y = "Normalized Volume Size",
       fill = "SuperAging Subtype") +
  theme(legend.box.background = element_rect(), # Add legend box
        legend.background = element_blank(), # Remove legend background
        legend.key = element_blank(), # Remove box around colors
        legend.position = "inside", # Move legend inside graph
        legend.position.inside = c(.93, .20), # Change position of legend
        panel.background = element_blank(), # Remove plot panel background
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.border = element_blank(), # Remove panel border
        axis.line = element_line(linewidth = .5, linetype = "solid"), # Change y axis line size
        plot.title = element_text(size = 14, hjust = .5, face = "bold") # Change title size and position
  )

ggplot(data = merged_hippo_final[merged_hippo_final$method == "T2H", ], 
       aes(x = subregion, y = volume_size_norm, fill = subtype)) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_manual(values = c("red",  "blue"), 
                    breaks = c("maintainer", "decliner"),
                    labels = c("Maintainer", "Decliner")) +
  labs(title = "Boxplots of SA Hippocampal Subregion Normalized Volume Size by Subtype for T2 Highres Only Method",
       x = "Subregion",
       y = "Normalized Volume Size",
       fill = "SuperAging Subtype") +
  theme(legend.box.background = element_rect(), # Add legend box
        legend.background = element_blank(), # Remove legend background
        legend.key = element_blank(), # Remove box around colors
        legend.position = "inside", # Move legend inside graph
        legend.position.inside = c(.93, .20), # Change position of legend
        panel.background = element_blank(), # Remove plot panel background
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.border = element_blank(), # Remove panel border
        axis.line = element_line(linewidth = .5, linetype = "solid"), # Change y axis line size
        plot.title = element_text(size = 14, hjust = .5, face = "bold") # Change title size and position
  )

ggplot(data = merged_hippo_final[merged_hippo_final$method == "T1_T2", ], 
       aes(x = subregion, y = volume_size_norm, fill = subtype)) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_manual(values = c("red",  "blue"), 
                    breaks = c("maintainer", "decliner"),
                    labels = c("Maintainer", "Decliner")) +
  labs(title = "Boxplots of SA Hippocampal Subregion Normalized Volume Size by Subtype for T1 & T2 Method",
       x = "Subregion",
       y = "Normalized Volume Size",
       fill = "SuperAging Subtype") +
  theme(legend.box.background = element_rect(), # Add legend box
        legend.background = element_blank(), # Remove legend background
        legend.key = element_blank(), # Remove box around colors
        legend.position = "inside", # Move legend inside graph
        legend.position.inside = c(.93, .20), # Change position of legend
        panel.background = element_blank(), # Remove plot panel background
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.border = element_blank(), # Remove panel border
        axis.line = element_line(linewidth = .5, linetype = "solid"), # Change y axis line size
        plot.title = element_text(size = 14, hjust = .5, face = "bold") # Change title size and position
  )

ggplot(data = merged_hippo_final[merged_hippo_final$method == "T1_T2H", ], 
       aes(x = subregion, y = volume_size_norm, fill = subtype)) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_manual(values = c("red",  "blue"), 
                    breaks = c("maintainer", "decliner"),
                    labels = c("Maintainer", "Decliner")) +
  labs(title = "Boxplots of SA Hippocampal Subregion Normalized Volume Size by Subtype for T1 & T2 Highres Method",
       x = "Subregion",
       y = "Normalized Volume Size",
       fill = "SuperAging Subtype") +
  theme(legend.box.background = element_rect(), # Add legend box
        legend.background = element_blank(), # Remove legend background
        legend.key = element_blank(), # Remove box around colors
        legend.position = "inside", # Move legend inside graph
        legend.position.inside = c(.93, .20), # Change position of legend
        panel.background = element_blank(), # Remove plot panel background
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.border = element_blank(), # Remove panel border
        axis.line = element_line(linewidth = .5, linetype = "solid"), # Change y axis line size
        plot.title = element_text(size = 14, hjust = .5, face = "bold") # Change title size and position
  )

ggplot(data = merged_amyg_final, aes(x = volume_size_norm, y = subregion, fill = method)) +
  geom_boxplot() +
  scale_fill_manual(values = c("green", "blue", "orange", "red", "purple"),
                    labels = c("T1", "T2 SPC", "T2 Highres", "T1+T2 SPC", "T1+T2 Highres")) + 
  labs(title = "Boxplots of SA Amygdala Subregion Volume Sizes & Subtypes by Segmentation Method",
       x = "Normalized Volume Size",
       y = "Subregion + Subtype",
       fill = "Segmentation Method") +
  theme(legend.box.background = element_rect(), # Add legend box
        legend.background = element_blank(), # Remove legend background
        legend.key = element_blank(), # Remove box around colors
        legend.position = "inside", # Move legend inside graph
        legend.position.inside = c(.93, .20), # Change position of legend
        panel.background = element_blank(), # Remove plot panel background
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.border = element_blank(), # Remove panel border
        axis.line = element_line(linewidth = .5, linetype = "solid"), # Change y axis line size
        plot.title = element_text(size = 14, hjust = .5, face = "bold") # Change title size and position
  )

ggplot(data = merged_amyg_final[merged_amyg_final$method == "T1", ], 
       aes(x = subregion, y = volume_size_norm, fill = subtype)) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_manual(values = c("red",  "blue"), 
                    breaks = c("maintainer", "decliner"),
                    labels = c("Maintainer", "Decliner")) +
  labs(title = "Boxplots of SA Amygdala Subregion Volume Sizes by Subtype for T1 Only Method",
       x = "Subregion",
       y = "Normalized Volume Size",
       fill = "SuperAging Subtype") +
  theme(legend.box.background = element_rect(), # Add legend box
        legend.background = element_blank(), # Remove legend background
        legend.key = element_blank(), # Remove box around colors
        legend.position = "inside", # Move legend inside graph
        legend.position.inside = c(.90, .20), # Change position of legend
        panel.background = element_blank(), # Remove plot panel background
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.border = element_blank(), # Remove panel border
        axis.line = element_line(linewidth = .5, linetype = "solid"), # Change y axis line size
        plot.title = element_text(size = 14, hjust = .5, face = "bold") # Change title size and position
  )

ggplot(data = merged_amyg_final[merged_amyg_final$method == "T2", ], 
       aes(x = subregion, y = volume_size_norm, fill = subtype)) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_manual(values = c("red",  "blue"), 
                    breaks = c("maintainer", "decliner"),
                    labels = c("Maintainer", "Decliner")) +
  labs(title = "Boxplots of SA Amygdala Subregion Volume Sizes by Subtype for T2 SPC Method",
       x = "Subregion",
       y = "Normalized Volume Size",
       fill = "SuperAging Subtype") +
  theme(legend.box.background = element_rect(), # Add legend box
        legend.background = element_blank(), # Remove legend background
        legend.key = element_blank(), # Remove box around colors
        legend.position = "inside", # Move legend inside graph
        legend.position.inside = c(.90, .20), # Change position of legend
        panel.background = element_blank(), # Remove plot panel background
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.border = element_blank(), # Remove panel border
        axis.line = element_line(linewidth = .5, linetype = "solid"), # Change y axis line size
        plot.title = element_text(size = 14, hjust = .5, face = "bold") # Change title size and position
  )

ggplot(data = merged_amyg_final[merged_amyg_final$method == "T2H", ], 
       aes(x = subregion, y = volume_size_norm, fill = subtype)) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_manual(values = c("red",  "blue"), 
                    breaks = c("maintainer", "decliner"),
                    labels = c("Maintainer", "Decliner")) +
  labs(title = "Boxplots of SA Amygdala Subregion Volume Sizes by Subtype for T2 Highres Method",
       x = "Subregion",
       y = "Normalized Volume Size",
       fill = "SuperAging Subtype") +
  theme(legend.box.background = element_rect(), # Add legend box
        legend.background = element_blank(), # Remove legend background
        legend.key = element_blank(), # Remove box around colors
        legend.position = "inside", # Move legend inside graph
        legend.position.inside = c(.90, .20), # Change position of legend
        panel.background = element_blank(), # Remove plot panel background
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.border = element_blank(), # Remove panel border
        axis.line = element_line(linewidth = .5, linetype = "solid"), # Change y axis line size
        plot.title = element_text(size = 14, hjust = .5, face = "bold") # Change title size and position
  )

ggplot(data = merged_amyg_final[merged_amyg_final$method == "T1_T2", ], 
       aes(x = subregion, y = volume_size_norm, fill = subtype)) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_manual(values = c("red",  "blue"), 
                    breaks = c("maintainer", "decliner"),
                    labels = c("Maintainer", "Decliner")) +
  labs(title = "Boxplots of SA Amygdala Subregion Volume Sizes by Subtype for T1 & T2 SPC Method",
       x = "Subregion",
       y = "Normalized Volume Size",
       fill = "SuperAging Subtype") +
  theme(legend.box.background = element_rect(), # Add legend box
        legend.background = element_blank(), # Remove legend background
        legend.key = element_blank(), # Remove box around colors
        legend.position = "inside", # Move legend inside graph
        legend.position.inside = c(.90, .20), # Change position of legend
        panel.background = element_blank(), # Remove plot panel background
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.border = element_blank(), # Remove panel border
        axis.line = element_line(linewidth = .5, linetype = "solid"), # Change y axis line size
        plot.title = element_text(size = 14, hjust = .5, face = "bold") # Change title size and position
  )

ggplot(data = merged_amyg_final[merged_amyg_final$method == "T1_T2H", ], 
       aes(x = subregion, y = volume_size_norm, fill = subtype)) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_manual(values = c("red",  "blue"), 
                    breaks = c("maintainer", "decliner"),
                    labels = c("Maintainer", "Decliner")) +
  labs(title = "Boxplots of SA Amygdala Subregion Volume Sizes by Subtype for T1 & T2 Highres Method",
       x = "Subregion",
       y = "Normalized Volume Size",
       fill = "SuperAging Subtype") +
  theme(legend.box.background = element_rect(), # Add legend box
        legend.background = element_blank(), # Remove legend background
        legend.key = element_blank(), # Remove box around colors
        legend.position = "inside", # Move legend inside graph
        legend.position.inside = c(.90, .20), # Change position of legend
        panel.background = element_blank(), # Remove plot panel background
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.border = element_blank(), # Remove panel border
        axis.line = element_line(linewidth = .5, linetype = "solid"), # Change y axis line size
        plot.title = element_text(size = 14, hjust = .5, face = "bold") # Change title size and position
  )

# Analysis ------------------------------------------------------------------
library(rstatix)
#install.packages("car")
library(car)
#install.packages("emmeans")
library(emmeans)
#install.packages("FSA")
library(FSA)
#install.packages("effsize")
library(effsize)
library(knitr)
# install.packages("kableExtra")
library(kableExtra)
# install.packages("gt")
library(gt)

# Remove hippocampal tail
merged_hippo_final <- merged_hippo_final[merged_hippo_final$subregion != "Hippocampal_tail", ]

# Normality check
by(merged_hippo_final$volume_size_norm, merged_hippo_final$method, shapiro.test)
leveneTest(volume_size_norm ~ method, data = merged_hippo_final)

# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(volume_size_norm ~ method, data = merged_hippo_final)
print(kruskal_result)
kruskal_result_df <- do.call(cbind, kruskal_result)
write.csv(kruskal_result_df,  "/Volumes/fsmresfiles/CNADC/Imaging_Core/Imaging/imaging_projects/Individuals/Joshua_Pasaye/FS_highreshippo/outputs/kruskal_results.csv")

# Perform pairwise Wilcoxon tests with Bonferroni adjustment
pairwise_result <- pairwise.wilcox.test(merged_hippo_final$volume_size_norm, 
                                        merged_hippo_final$method,
                                        paired = TRUE,
                                        p.adjust.method = "bonferroni")

# View pairwise results
print(pairwise_result)

# Get unique methods
methods <- c("T1", "T2", "T2H", "T1_T2", "T1_T2H")

# Create an empty list to store results
wilcoxon_results <- list()

# Perform pairwise Wilcoxon signed-rank tests and store results
for (i in 1:(length(methods) - 1)) {
  for (j in (i + 1):length(methods)) {
    group1 <- merged_hippo_final$volume_size_norm[merged_hippo_final$method == methods[i]]
group2 <- merged_hippo_final$volume_size_norm[merged_hippo_final$method == methods[j]]

test_result <- wilcox.test(group1, group2, paired = TRUE, exact = FALSE)

wilcoxon_results[[paste(methods[i], methods[j], sep = " vs. ")]] <- data.frame(
  Method1 = methods[i],
  Method2 = methods[j],
  W = test_result$statistic,  # Test statistic
  p_value = format(test_result$p.value, scientific = FALSE)  # Unadjusted p-value
)
  }
}

wilcoxon_results_df <- do.call(rbind, wilcoxon_results)

# Adjust p-values using Bonferroni correction
wilcoxon_results_df$adjusted_p_value <- format(p.adjust(wilcoxon_results_df$p_value, method = "bonferroni"), scientific = FALSE)

# View the results
print(wilcoxon_results_df)

# Output results
write_csv(wilcoxon_results_df, "/Volumes/fsmresfiles/CNADC/Imaging_Core/Imaging/imaging_projects/Individuals/Joshua_Pasaye/FS_highreshippo/outputs/wilcoxon_results.csv")

# Table
wilcoxon_poster_table <- wilcoxon_results_df %>%  # Start with original dataframe
  mutate(
    # First ensure numeric p-values exist
    p_value_num = as.numeric(p_value),
    adjusted_p_value_num = as.numeric(adjusted_p_value),
    
    # Create significance stars
    significance = case_when(
      adjusted_p_value_num < 0.001 ~ "***",
      adjusted_p_value_num < 0.01 ~ "**",
      adjusted_p_value_num < 0.05 ~ "*",
      TRUE ~ ""
    ),
    
    # Format values for display
    `Unadjusted p` = ifelse(p_value_num < 0.0001, 
                            "< 0.0001***",
                            sprintf("%.4f%s", p_value_num, significance)),
    `Adjusted p` = ifelse(adjusted_p_value_num < 0.0001,
                          "< 0.0001***",
                          sprintf("%.4f%s", adjusted_p_value_num, significance)),
    `W Statistic` = format(W, big.mark = ","),
    
    # Clean method names if needed
    Comparison = paste(Method1, "vs.", Method2) %>%
      str_replace("T1_T2", "T1 & T2") %>%
      str_replace("T1_T2H", "T1 & T2H")
  ) %>%
  select(Comparison, `W Statistic`, `Unadjusted p`, `Adjusted p`)

poster_table <- wilcoxon_poster_table %>%
  gt() %>%
  tab_header(
    title = md("**Table 1. Pairwise Wilcoxon Signed-Rank Tests**"),
    subtitle = "Bonferroni-corrected comparisons of hippocampal subfield volumes across methods"
  ) %>%
  cols_align(
    align = "center",
    columns = c(`W Statistic`, `Unadjusted p`, `Adjusted p`)
  ) %>%
  cols_align(
    align = "left",
    columns = Comparison
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = grepl("\\*", `Adjusted p`) # Applies to all columns in significant rows
    )
  ) %>%
  tab_style(
    style = cell_fill(color = "white"),
    locations = cells_body(columns = everything(), rows = everything())
  ) %>%
  tab_footnote(
    footnote = "*** p < 0.001; ** p < 0.01; * p < 0.05",
    locations = cells_column_labels(columns = `Adjusted p`)
  ) %>%
  tab_options(
    table.font.size = 18,
    column_labels.font.weight = "bold",
    table.width = "100%",
    row_group.font.weight = "bold"
  )

# Save as high-resolution image
gtsave(poster_table, 
       filename = "/Volumes/fsmresfiles/CNADC/Imaging_Core/Imaging/imaging_projects/Individuals/Joshua_Pasaye/FS_highreshippo/plots/Table1_Wilcoxon_Results.png",
       zoom = 2,
       expand = 10)

# Predictive models -------------------------------------------------------
library(sjPlot)
#install.packages("performance")
library(performance)
#install.packages("effects")
library(effects)
library(kableExtra)
#install.packages("ggsignif")
library(ggsignif)

# Model
model <- lm(volume_size_norm ~ method * subtype + subregion,
              data = merged_hippo_final)
results <- summary(model)
results

anova_results <- anova(model)
anova_results

# Post-hoc
emm <- emmeans(model, ~ method | subtype, adjust = "bonferroni")
summary(emm, infer = TRUE)
pairs(emm, adjust = "bonferroni")
subtype_emm <- emmeans(model, ~ subtype)
pairs(subtype_emm, adjust = "bonferroni")
summary(subtype_emm, infer = TRUE)

#Plot
sig_data <- data.frame(
  x = 1,    # x-positions for T1
  xend = 3, # xend for T1+T2H
  y = c(560), # y-positions
  t_value = c("t = -.932"),
  p_label = c("p = .321")
)

# Plot with significance annotations
figure1 <- ggplot(emm_df, 
                  aes(x = method, 
                      y = emmean, 
                      fill = subtype)) +
  # Bar plot for model estimates
  geom_bar(stat = "identity",
           position = position_dodge(0.8),
           width = 0.7) +
  # Error bars for confidence intervals
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(0.8),
                width = 0.25,
                color = "black",
                linewidth = 0.8) +
  geom_segment(
    data = sig_data,
    aes(x = x, xend = xend, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = sig_data,
    aes(x = (x + xend)/1.6, y = y + 17, label = p_label),
    inherit.aes = FALSE,
    size = 10,
  ) +
  geom_text(
    data = sig_data,
    aes(x = (x + xend)/2, y = y + 17, label = t_value),
    inherit.aes = FALSE,
    size = 10
  ) +
  # Add data points (if desired)
  geom_point(aes(group = subtype),
             position = position_dodge(0.8),
             size = 2,
             color = "black",
             show.legend = FALSE) +
  labs(
    title = "Figure 1. Hippocampal Volume by Segmentation Method and Subtype",
    subtitle = "Error bars show 95% CIs; All methods significantly different to T1 besides combined T1 & T2H",
    y = "Adjusted Volume (normalized)", 
    x = "Segmentation Method",
    fill = "Subtype"
  ) +
  theme_linedraw() +
  scale_fill_brewer(
    palette = "Dark2",
    labels = c("Decliner", "Maintainer")
  ) +
  scale_x_discrete(
    labels = c("T1", "T1 & T2", "T1 & T2H", "T2", "T2H")
  ) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 28, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 20, hjust = 0),
    axis.title = element_text(face = "bold", size = 24),
    axis.text.x = element_text(size = 22), 
    axis.text.y = element_text(size = 22),
    legend.position = "top",
    legend.box.background = element_rect(color = "black"),
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 20)
  )
# Print the plot
print(figure1)

# Save figure
ggsave("/Volumes/fsmresfiles/CNADC/Imaging_Core/Imaging/imaging_projects/Individuals/Joshua_Pasaye/FS_highreshippo/plots/poster_figure1.png", figure1, 
       width = 18, height = 8, units = "in", dpi = 300)

# Export summary of model
library(openxlsx)

## Model 
coefficients <- results$coefficients
r_squared <- results$r.squared
adj_r_squared <- results$adj.r.squared
f_statistic <- results$fstatistic
p_value <- format.pval(pf(f_statistic[1], 
                          f_statistic[2], 
                          f_statistic[3], 
                          lower.tail= FALSE),
                       digits = 3)
model_formula <- paste0(
  results[["terms"]][[2]], " ", 
  results[["terms"]][[1]], " ",
  results[["terms"]])[[3]]

## Create workbook
wb <- createWorkbook()
addWorksheet(wb, "Model Summary")
## Add columns
writeData(wb, "Model Summary", "Coefficients", 
          startRow = 1, startCol = 1)
writeData(wb, "Model Summary", coefficients, startRow = 2, 
          startCol = 1, rowNames = TRUE)

writeData(wb, "Model Summary", "R-Squared", 
          startRow = 2 + nrow(coefficients) + 2, startCol = 1)
writeData(wb, "Model Summary", r_squared, 
          startRow = 2 + nrow(coefficients) + 2, startCol = 2)

writeData(wb, "Model Summary", "Adjusted R-Squared", 
          startRow = 2 + nrow(coefficients) + 3, startCol = 1)
writeData(wb, "Model Summary", adj_r_squared, 
          startRow = 2 + nrow(coefficients) + 3, startCol = 2)

writeData(wb, "Model Summary", "F-Statistic", 
          startRow = 2 + nrow(coefficients) + 4, startCol = 1)
writeData(wb, "Model Summary", f_statistic[1], 
          startRow = 2 + nrow(coefficients) + 4, startCol = 2)

writeData(wb, "Model Summary", "p-Value", 
          startRow = 2 + nrow(coefficients) + 5, startCol = 1)
writeData(wb, "Model Summary", p_value, 
          startRow = 2 + nrow(coefficients) + 5, startCol = 2)

writeData(wb, "Model Summary", "Model Formula", 
          startRow = 2 + nrow(coefficients) + 6, startCol = 1)
writeData(wb, "Model Summary", model_formula, 
          startRow = 2 + nrow(coefficients) + 6, startCol = 2)
## Save workbook
saveWorkbook(
  wb, 
  file = "/Volumes/fsmresfiles/CNADC/Imaging_Core/Imaging/imaging_projects/Individuals/Joshua_Pasaye/FS_highreshippo/outputs/lm_model_summary.xlsx", 
  overwrite = TRUE
)
