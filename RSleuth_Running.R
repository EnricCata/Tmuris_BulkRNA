# -*- coding: utf-8 -*-
library(sleuth)
library(tidyverse)
library(stringr)
library(ggplot2)

### Add GO terms to later plot
library(readr)
eggnogTmuris <- read_tsv(
  "/path/to/GO_annotations",
  comment = "##",
  col_types = cols(.default = "c")
)
eggnogTmuris <- eggnogTmuris %>%
  rename(query = `#query`)

# +
head(eggnogTmuris)
colnames(eggnogTmuris)

### annotation over
# -

# Sleuth analysis

head(eggnogTmuris$GOs)

# specify where the kallisto resuts are stored
sample_ids <- dir(file.path("path/to/kallisto_mapped_samples"))
sample_ids

# Remove all the .tpm files to get only the folders in this
sample_ids <- sample_ids[!grepl("\\.tpm$", sample_ids)]
remove_entries <- c("transcripts.list", "kallisto_allsamples.tpm.table")
sample_ids <- setdiff(sample_ids, remove_entries)
sample_ids

#Set up the file paths for the sleuth object
kal_dirs <- file.path("path/to/kallisto_mapped_samples", sample_ids)
kal_dirs

# +
# Prepare the metadata

sampleinfo <- read_tsv("/path/to/tmuris_bulkrnaseq_metadata.tsv", col_types = c("cccc")) %>%
  arrange(lane_id, host, time_point, new_id)
sampleinfo$lane_id <- gsub('#', '-', sampleinfo$lane_id)

#there's an empty row, remove it
sampleinfo <- sampleinfo %>%
  filter(!if_all(everything(), is.na))
class(sampleinfo)

sampleinfo

# +
#Change metadata as these samples are incorrectly labelled
sampleinfo <- sampleinfo %>%
  mutate(time_point = if_else(lane_id == "49707_1-13", 'AF', time_point)) %>%
  mutate(time_point = if_else(lane_id == "49707_1-24", 'AF', time_point)) %>%
  mutate(time_point = if_else(lane_id == "49707_1-2", 'AM', time_point)) %>%
  mutate(time_point = if_else(lane_id == "49707_1-1", 'AM', time_point)) %>%
  mutate(time_point = if_else(lane_id == "49707_1-3", 'AM', time_point)) 

sampleinfo
# -

class(kal_dirs)

# +
# Adding 'path' to metadata. Rename lane_id to 'sample' otherwise SLEUTH WILL NOT RECOGNISE IT!
sampleinfo <- sampleinfo %>%
  mutate(
    folder_name = str_replace(lane_id, "-", "#") %>%
      paste0("_kallisto_out"),
    path = kal_dirs[match(folder_name, basename(kal_dirs))]
  ) %>%
  select(-folder_name)

sampleinfo <- sampleinfo %>%
  rename(sample = lane_id)

#add a grouping, aka now host and time_point are fused into one(?) 
sampleinfo$group <- interaction(sampleinfo$host, sampleinfo$time_point)

sampleinfo

# +
# Make a new column where larval stage is the same more or less
sampleinfo <- sampleinfo %>%
  mutate(larval_stage = case_when(
    time_point == "L1"          ~ "L1",
    time_point %in% c("D13","D14","D15") ~ "L2",
    time_point %in% c("D20") ~ "L3",
    time_point %in% c("D25","D24") ~ "L4",
    time_point %in% c("D33") ~ "L5-YA",
    time_point %in% c("AF") ~ "L5-AF",
    time_point %in% c("AM") ~ "L5-AM",
    TRUE ~ time_point  # keep other values as-is if any
  ) )

# Check the result
table(sampleinfo$host_grouped)
sampleinfo


# +
sampleinfo <- sampleinfo %>%
  dplyr::mutate(
    larval_stage = dplyr::case_when(
      grepl("^49707_1-(11|13|10|9|24|5|22|19|2|15|14)$", sample) ~ paste0(larval_stage, "-J"),
      grepl("^49707", sample)                                   ~ paste0(larval_stage, "-M"),
      TRUE                                                      ~ larval_stage
    )
  )

table(sampleinfo$larval_stage)
sampleinfo

# +
# Trimming out things, in this case removing adults and hatchlings
# Adding ! to time_point turns it negative. Aka, the condition that matches is selected OUT Removing it selects it in


sampleinfo_mouse_only <- filter(sampleinfo, host %in% c('MOU', 'IVI'))
sampleinfo_trim <-  filter(sampleinfo_mouse_only, !larval_stage %in% c('L5-YA'))
sampleinfo_trim


# +
#sampleinfo_trim <- sampleinfo

sampleinfo_trim <- filter(sampleinfo_trim, sample != '49707_1-3'#, sample != "49707_1-11")
)# Removing the funky outlier and the too juvenile female
sampleinfo_trim <- filter(sampleinfo_trim, larval_stage != 'L5-YA')
sampleinfo_trim
sampleinfo_trim

# +
# checking if this will work with Sleuth. The conditions are group

model.matrix(~ group, data = sampleinfo_trim)
# -
#construct sleuth object
#contains info about experiment, model details for differential, and results! Commands will
#load kallisto data into the object
so_prep <- sleuth_prep(sampleinfo_trim, extra_bootstrap_summary = TRUE)
message("Loaded data onto object")



#estimate parameters for sleuth response error measurement (full model). 
#Remember you have two conditions here: time_point and host, write them as they are in the data frame
so <- sleuth_fit(so_prep, ~larval_stage, 'full')
message("Error measurement parameters estimated")

# +
# ## Likelyhood Ratio Test
# ### You're checking differences between a full model and a reduced model where there's no change.

#estimate parameters for sleuth REDUCED model
so <- sleuth_fit(so, ~1, 'reduced')
message("Reduced model parameters estimated")

#perform differential analysis using likelyhood ratio test
so <- sleuth_lrt(so, 'reduced', 'full')
message("Performed differential analysis")
# -

# ## Wald Test
# ### You're checking for differences in the coefficient, one model! 
# ### When comparing larval stages, use the one you wish to see the difference of as part of the input: larval_stage(stage of interest)
# ### When comparing mice vs organoid, its called larval_stageORG because the larval_stage is the baseline, ORG is the independent variable
so <- sleuth_wt(so, 'larval_stage[LX]')
message("Wald Test done")

#For Volcano plots, employ so as its already analyzed
#For heatmaps, employ so_prep as its a collection of normalized Sleuth objects
#Repeat for every comparison

message("Saving Sleuth object ...")
saveRDS(so, file = "/path/to/object_so.rds")
message("Done! Sleuth object saved successfully.")

message("Saving Sleuth prep object ...")
saveRDS(so_prep, file = "/path/to/object_so_prep.rds")
message("Done! Sleuth object saved successfully.")

