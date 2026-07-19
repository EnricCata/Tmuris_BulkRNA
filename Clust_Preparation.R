library(readr)
library(dplyr)
library(stringr)


# Metadata

sampleinfo <- read_tsv("/path/to/tmuris_bulkrnaseq_metadata.tsv", col_types = c("cccc")) %>%
  arrange(lane_id, host, time_point, new_id)
sampleinfo$lane_id <- gsub('#', '-', sampleinfo$lane_id)

#there's an empty row, remove it
sampleinfo <- sampleinfo %>%
  filter(!if_all(everything(), is.na))
class(sampleinfo)

sampleinfo

# +
sampleinfo <- sampleinfo %>%
  mutate(time_point = if_else(lane_id == "49707_1-13", 'AF', time_point)) %>%
  mutate(time_point = if_else(lane_id == "49707_1-24", 'AF', time_point)) %>%
  mutate(time_point = if_else(lane_id == "49707_1-2", 'AM', time_point)) %>%
  mutate(time_point = if_else(lane_id == "49707_1-1", 'AM', time_point)) %>%
  mutate(time_point = if_else(lane_id == "49707_1-3", 'AM', time_point)) 


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


sampleinfo <- sampleinfo %>%
  mutate(model = case_when(
    host %in% c("ORG1", "ORG2") ~ "ORG",
    TRUE ~ host
  ))


sampleinfo <- sampleinfo %>%
  mutate(new_name = paste(larval_stage, model, str_extract(lane_id, "[^-]+$"),  sep = "-"))

# +


#Clust normalized input

new_names <- readLines("/path/to/fastq-raws/SAMPLE_ID_FULL.txt")

clust_input <- read.table("path/to/kallisto_mapped_samples", header=FALSE, fill = TRUE, row.names=1, sep="\t" )
clust_input <- clust_input[-1, ] %>% # database is screwy, deletes the empty ID row
  mutate(across(everything(), ~ as.numeric(as.character(.x)))) %>% #makes sure you have everything as numbers
  { colnames(.) <- new_names; . } %>%     # assign new column names with IDs
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # adds 1 to all TPM values
  mutate(across(where(is.numeric), log2)) %>% # log transform all TPM+1 values
  as.matrix() #make into a large matrix and flip row-columns 
all(is.finite(clust_input))  #checks for empty shit and non-numerics

# Rename columns
colnames(clust_input) <- (setNames(sampleinfo$new_name, sampleinfo$lane_id))[colnames(clust_input)]

# Remove young adults & outlier if needed
clust_input <- clust_input[, !grepl("YA", colnames(clust_input))]
clust_input <- clust_input[, !grepl("L5.AM.MOU.3", colnames(clust_input))]

# Remove vitro if needed
clust_input <- clust_input[, !grepl("ORG", colnames(clust_input))]

# Remove mouse if needed
clust_input_vitro <- clust_input[, !grepl("MOU", colnames(clust_input))]

#test to see how the file looks like
#l1_clust_input <- clust_input[, grepl("L1", colnames(clust_input))]

# Gene column 
clust_input_vitro <- data.frame(ID = rownames(clust_input), clust_input)

write.table(clust_input_vitro,
            file = "clust_input_vitro.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# +


# Clust raw TPM

clust_raw_input <- read.table("path/to/kallisto_mapped_samples", header=FALSE, fill = TRUE, row.names=1, sep="\t" )
clust_raw_input <- clust_raw_input[-1, ] %>% # database is screwy, deletes the empty ID row
  mutate(across(everything(), ~ as.numeric(as.character(.x)))) %>% #makes sure you have everything as numbers
  { colnames(.) <- new_names; . } %>%     # assign new column names with IDs
  as.matrix() #make into a large matrix and flip row-columns 
all(is.finite(clust_raw_input))  #checks for empty shit and non-numerics

# Rename columns
colnames(clust_raw_input) <- (setNames(sampleinfo$new_name, sampleinfo$lane_id))[colnames(clust_raw_input)]

# Remove young adults & outlier
clust_raw_input <- clust_raw_input[, !grepl("YA", colnames(clust_raw_input))]
clust_raw_input <- clust_raw_input[, !grepl("L5.AM.MOU.3", colnames(clust_raw_input))]

# Remove vitro
clust_raw_input <- clust_raw_input[, !grepl("ORG", colnames(clust_raw_input))]

#Remove mouse
clust_raw_input_vitro <- clust_raw_input[, !grepl("MOU", colnames(clust_raw_input))]


# Gene column 
clust_raw_input_vitro <- data.frame(ID = rownames(clust_raw_input_vitro), clust_raw_input_vitro)

write.table(clust_raw_input_vitro,
            file = "clust_raw_input_vitro.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)



#Reference

stages <- c("L1", "L2", "L3", "L4", "L5.AM", "L5.AF")

reference <- data.frame(
  file = "clust_raw_input.tsv",
  stage = c("L1", "L2", "L3", "L4", "L5.AM", "L5.AF"),
  samples = sapply(c("L1", "L2", "L3", "L4", "L5.AM", "L5.AF"), function(s) {
    paste(colnames(clust_raw_input)[grepl(s, colnames(clust_raw_input))], collapse = ", ")
  })
)

write.table(reference,
            file = "replicates_file_raw.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


#Reference vitro

stages <- c("L1", "L2", "L3", "L4")

reference <- data.frame(
  file = "clust_raw_input_vitro.tsv",
  stage = c("L1", "L2", "L3", "L4"),
  samples = sapply(c("L1", "L2", "L3", "L4"), function(s) {
    paste(colnames(clust_raw_input_vitro)[grepl(s, colnames(clust_raw_input_vitro))], collapse = ", ")
  })
)

write.table(reference,
            file = "replicates_file_raw_vitro.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


#Command employed for clust in the command line:

#clust /lustre/scratch127/pam/teams/team333/es37/trichuris_muris/bulk-rnaseq/04_analysis/Clust/clust_raw_input_vitro.tsv -o /lustre/scratch127/pam/teams/team333/es37/trichuris_muris/bulk-rnaseq/04_analysis/Clust/Clust_Vitro_raw_default_norm_tight_2 -r /lustre/scratch127/pam/teams/team333/es37/trichuris_muris/bulk-rnaseq/04_analysis/Clust/replicates_file_raw_vitro.txt -t 2 


# +
clust_results_6  <- read_tsv("/path/to/Clust/Clust_Attempt_6_default_norm_tight_2/Clusters_Objects.tsv" )

clust_results_6 <- clust_results_6[-1,]
head(clust_results_6 )



