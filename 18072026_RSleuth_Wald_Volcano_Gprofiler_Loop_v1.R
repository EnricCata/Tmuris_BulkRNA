library(sleuth)
library(tidyverse)
library(stringr)
library(ggplot2)
library(dplyr)


### This whole part is for ease of annotation. 
library(readr)

eggnogTmuris <- read_tsv(
  "/path/to/GO_annotations",
  comment = "##",
  col_types = cols(.default = "c")
)
eggnogTmuris <- eggnogTmuris %>%
  rename(query = `#query`)

#Volcano Plot Loops
## Here the so files from the sleuth Wald tests are employed

# Input folder
#Vivo
rds_dir <- "/path/to/Sleuth_Objects/In-vivo"

#Vitro
rds_dir <- "/path/to/Sleuth_Objects/In-vitro"

#Vitro vs Vivo
rds_dir <- "/path/to/Sleuth_Objects/Vitro-vs-Vivo"

# output folder
out_dir <- "/path/to/output_dir"

# Sleuthobject List
files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)

# Regulation data information, these will be logged in during the loops to use for the GProfiler profiling
all_differential_vivo <- data.frame()
all_differential_vitro <- data.frame()
all_differential_vivo_vitro <- data.frame()



for (file_path in files) {

# +
so <- readRDS(file_path)

# +

models(so)

# Way to use the intercept part!
coef_to_use <- sub("(?s).*\\(Intercept\\)\\s*", "", paste(capture.output(models(so)), collapse = "\n"), perl = TRUE)


# For Wald tests and volcano plots
sleuth_table <- sleuth_results(so, coef_to_use, 'wt', show_all = FALSE)
#head(sleuth_table, 20)

# +
# Important values:
# b, log fold change
# pval, raw p-value
# qval, FDR-adjusted p-value, false discovery rate

# Volcano plot with all, threshold applied on plot

volcano <- sleuth_table

top_pos_fc <- volcano %>%
  slice_max(b, n = 10)
top_neg_fc <- volcano %>%
  slice_min(b, n = 10)
top_pval <- volcano %>%
  slice_min(pval, n = 10)

top_hits <- bind_rows(top_pos_fc, top_neg_fc, top_pval) %>%
  distinct(target_id, .keep_all = TRUE)

top_hits  <- top_hits %>%
  mutate(pval = ifelse(pval <= 0, 1e-320, pval))
volcano  <- volcano %>%
  mutate(pval = ifelse(pval <= 0, 1e-320, pval))



#top_hits

top_hits <- top_hits %>%
  left_join(
    eggnogTmuris %>% dplyr::select(query, Description),
    by = c("target_id" = "query")
  ) %>%
  mutate(
    Description = if_else(
      is.na(Description) | Description == "-",
      "No annotation",
      Description
    )
  )
#top_hits

file_name <- paste0(gsub("RSleuthobject_", "", tools::file_path_sans_ext(basename(file_path))), "_Volcano.pdf")

volcanoplot <- ggplot(volcano, aes(x = b, y = -log10(pval), color = ifelse(b > 1 & qval <= 0.05 | b < -1 & qval <= 0.05, "Significant", "Moderate"))) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(0.005), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -1, linetype = "dashed", color = "black") +
  xlim(-15, 15) + # X axis limits
  #  ylim(0, 325) +          # Y axis limits
  ylim(0, 325) +           # For the vivo vs vitro
  # this is in case there's a pvalue too low to calculate
  scale_color_manual(values = c("grey70", "red")) +
  #  geom_point(data = line_genes, color = "green", size = 2) +
  theme_bw() +
  #  geom_text(
  #    data = top_hits,
  #    aes(label = target_id),  # 
  #    vjust = -0.5,
  #    size = 3,
  #    check_overlap = TRUE
  #  ) +
  labs(
    x = "Log2 fold change",
    # b is not log2 or any specific, its model-unique and not able to be taken out!
    y = "-log10(p-value)",
    color = "FDR ≤ 0.05",
    #    caption = 'L5 Male Mature vs L5 Female Juvenile (Female as reference)'
  )

print(volcanoplot)

volcano_upregulated <- volcano %>%
  dplyr::filter(b > 1 & qval <= 0.05)%>%
  dplyr::mutate(Regulation = "Upregulated")
message('Upregulated')
nrow(volcano_upregulated)

volcano_downregulated <- volcano %>%
  dplyr::filter(b < -1 & qval <= 0.05)%>%
  dplyr::mutate(Regulation = "Downregulated")
message('Downregulated')
nrow(volcano_downregulated)

volcano_differential <- dplyr::bind_rows(
  volcano_upregulated,
  volcano_downregulated
) %>%
  dplyr::mutate(file_name = gsub(".pdf", "", file_name))

all_differential_vivo_vitro <- dplyr::bind_rows(
  all_differential_vivo_vitro,
  volcano_differential
)

# To save the damn thing
out_pdf <- file.path(out_dir, file_name)
ggsave(out_pdf, plot = volcanoplot, width = 7, height = 6)
message("Saved: ", file_name)
}


all_differential_vivo_vitro$file_name <- gsub("Host_", "Vivo_Vitro_", all_differential_vivo_vitro$file_name)

all_differential_vivo_vitro <- all_differential_vivo_vitro %>%
dplyr::select(
  target_id,
  pval,
  qval,
  b,
  Regulation,
  file_name
)

# Add annotation

all_differential_vivo_vitro <- all_differential_vivo_vitro %>%
  dplyr::left_join(
    eggnogTmuris %>%
      dplyr::select(query, Preferred_name, Description, PFAMs, GOs, KEGG_ko),
    by = c("target_id" = "query")
  )

all_differential_vitro <- all_differential_vitro %>%
  dplyr::left_join(
    eggnogTmuris %>%
      dplyr::select(query, Preferred_name, Description, PFAMs, GOs, KEGG_ko),
    by = c("target_id" = "query")
  )

all_differential_vivo <- all_differential_vivo %>%
  dplyr::left_join(
    eggnogTmuris %>%
      dplyr::select(query, Preferred_name, Description, PFAMs, GOs, KEGG_ko),
    by = c("target_id" = "query")
  )


# Check individual
differential_vivo_l1_l2 <- all_differential_vivo %>%
dplyr::filter(stringr::str_detect(file_name, "RSleuthobjectL1_L2_Wald_Volcano"))

write.csv(all_differential_vivo, file.path(out_dir, "all_differential_vivo.csv"), row.names = FALSE)
write.csv(all_differential_vitro, file.path(out_dir, "all_differential_vitro.csv"), row.names = FALSE)
write.csv(all_differential_vivo_vitro, file.path(out_dir, "all_differential_vivo_vitro.csv"), row.names = FALSE)




# GProfiler Preparation

## Preparing GProfiler input files

library(readr)

all_differential_vivo <- read.csv("/path/to/all_differential_vivo.csv")

differential_vivo_primary <- all_differential_vivo[all_differential_vivo$file_name %in% c("Vivo_L1_L2_Wald_Volcano", "Vivo_L2_L3_Wald_Volcano","Vivo_L3_L4_Wald_Volcano",  "Vivo_L4_L5-AM_Wald_Volcano", "Vivo_L4_L5-AF_Wald_Volcano", "Vivo_L5-AF_L5-AM_Wald_Volcano"), ]
unique(differential_vivo_primary$file_name)

differential_vivo_l1_l2 <- all_differential_vivo[all_differential_vivo$file_name %in% c("Vivo_L1_L2_Wald_Volcano"), ]
differential_vivo_l2_l3 <- all_differential_vivo[all_differential_vivo$file_name %in% c("Vivo_L2_L3_Wald_Volcano"), ]
differential_vivo_l3_l4 <- all_differential_vivo[all_differential_vivo$file_name %in% c("Vivo_L3_L4_Wald_Volcano"), ]
differential_vivo_l4_male <- all_differential_vivo[all_differential_vivo$file_name %in% c("Vivo_L4_L5-AM_Wald_Volcano"), ]
differential_vivo_l4_female <- all_differential_vivo[all_differential_vivo$file_name %in% c("Vivo_L4_L5-AF_Wald_Volcano"), ]
differential_vivo_female_male <- all_differential_vivo[all_differential_vivo$file_name %in% c("Vivo_L5-AF_L5-AM_Wald_Volcano"), ]

write.csv(differential_vivo_l1_l2, file.path(out_dir, "differential_vivo_l1_l2.csv"), row.names = FALSE)
write.csv(differential_vivo_l2_l3, file.path(out_dir, "differential_vivo_l2_l3.csv"), row.names = FALSE)
write.csv(differential_vivo_l3_l4, file.path(out_dir, "differential_vivo_l3_l4.csv"), row.names = FALSE)
write.csv(differential_vivo_l4_male, file.path(out_dir, "differential_vivo_l4_male.csv"), row.names = FALSE)
write.csv(differential_vivo_l4_female, file.path(out_dir, "differential_vivo_l4_female.csv"), row.names = FALSE)
write.csv(differential_vivo_female_male, file.path(out_dir, "differential_vivo_female_male.csv"), row.names = FALSE)


all_differential_vitro <- read.csv("/path/to/all_differential_vitro.csv")
all_differential_vivo_vitro <- read.csv("/path/to/all_differential_vitro_vivo.csv")


differential_vivo_vitro_primary <- rbind(
  all_differential_vivo_vitro,
  all_differential_vitro[all_differential_vitro$file_name == "Vitro_L1_L2_Wald_Volcano", ]
)

unique(differential_vivo_vitro_primary$file_name)

#Function to output the target genes for gprofiler, ORDERED QUERY

filter_regulation <- function(df, file, regulation) {
  
  file_match <- grepl(file, df$file_name)
  
  subset_df <- df[
    file_match &
      if (tolower(regulation) == "both") {
        df$Regulation %in% c("Upregulated", "Downregulated")
      } else {
        df$Regulation == regulation
      },
  ]
  
  selected_ids <- subset_df$target_id[order(subset_df$pval)]
  
  #cat("Number of selected genes:", length(selected_ids), "\n")
  
  paste(selected_ids, collapse = " ")
}

filter_regulation(differential_vivo_primary,'L1_L2', 'Upregulated')


#Loop to prepare the files

# output folder
out_dir <- "/data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/04_analysis/R-files"
text_file_names <- unique(differential_vivo_vitro_primary$file_name)

for (f in text_file_names) {
  
  for (reg in c("Upregulated", "Downregulated", "both")) {
    
    genes <- filter_regulation(differential_vivo_vitro_primary, file = f, regulation = reg)
    
    file_suffix <- if (reg == "both") "Both" else reg
    
    out_file <- file.path(
      out_dir,
      paste0(f, "_", file_suffix, "_GprofileInput.txt")
    )
    
    writeLines(genes, out_file)
  }
  message("Saved: ", out_file)
}


#Plotting GProfiler

library(GO.db)
library(AnnotationDbi)
library(KEGGREST)

out_dir <- "/path/to//out_dir/Gprofiler_Graphs"


#Import gprofiler outputs
## You will import ALL Gprofiler outputs of a certain type into one list of lists. Either In-vivo, In-vitro, Vivo_Vs_Vitro
gprofiler_files <- list.files(path = "/path/to/Gprofiler_Output/Vivo_Vs_Vitro", pattern = "\\.csv$", full.names = TRUE)

gprofiler_list <- lapply(gprofiler_files, read.csv)
names(gprofiler_list) <- basename(gprofiler_files)

## Choose ONE of the lists to run the GProfiler plot through
gprofiler_plot <- gprofiler_list[["Vivo_Vitro_L4_Downregulated.csv"]]

# Adding term descriptions


go_gprofiler_plot <- AnnotationDbi::select(
  GO.db,
  keys = gprofiler_plot$term_id,
  keytype = "GOID",
  columns = c("TERM", "ONTOLOGY")
)


gprofiler_plot <- gprofiler_plot %>%
  left_join(
    go_gprofiler_plot,
    by = c("term_id" = "GOID")
  )

# Kegg annotation

gprofiler_plot[grepl("ko", gprofiler_plot$term_id), "ONTOLOGY"] <- "KEGG"


kegg_ids <- gprofiler_plot %>%
  filter(grepl("^ko:K\\d+", term_id)) %>%
  mutate(KO = sub("^ko:", "", term_id)) %>%
  pull(KO) %>%
  unique()

#ko_info <- keggList("ko") #optimize later do the opposite, make list and look it up instead of bringing everything
#ko_volcano_combined$KO_term <- ko_info[ko_volcano_combined$KEGG_ko]

ko_info <- lapply(kegg_ids, function(id) {
  tryCatch(KEGGREST::keggGet(id)[[1]], error = function(e) NULL)
})

ko_map <- setNames(
  sapply(ko_info, function(x) x$NAME[1]),
  sapply(ko_info, function(x) sub("ko:", "", x$ENTRY))
)
ko_map <- sub(" \\[.*\\]", "", ko_map) # clean descripts

gprofiler_plot <- gprofiler_plot %>%
  mutate(
    KO = sub("^ko:", "", term_id),
    TERM = ifelse(
      grepl("^ko:K\\d+", term_id),
      ko_map[KO],
      TERM
    )
  ) %>%
  dplyr::select(-KO)



# manually check obsolete and delete redundant terms like Biological Process
gprofiler_plot <- gprofiler_plot[
  gprofiler_plot$term_id != "GO:0003705",]

gprofiler_plot[gprofiler_plot$term_id=="GO:0005887","TERM"] <- 'integral component of the plasma membrane'
gprofiler_plot[gprofiler_plot$term_id=="GO:0005887", 'ONTOLOGY'] <- 'CC'

gprofiler_plot[gprofiler_plot$term_id=="GO:0044459","TERM"] <- 'plasma membrane part'
gprofiler_plot[gprofiler_plot$term_id=="GO:0044459", 'ONTOLOGY'] <- 'CC'

gprofiler_plot[gprofiler_plot$term_id=="GO:0044421","TERM"] <- 'extracellular region'
gprofiler_plot[gprofiler_plot$term_id=="GO:0044421", 'ONTOLOGY'] <- 'CC'

gprofiler_plot[gprofiler_plot$term_id=="GO:0031226","TERM"] <- 'intrinsic component of plasma membrane'
gprofiler_plot[gprofiler_plot$term_id=="GO:0031226", 'ONTOLOGY'] <- 'CC'

gprofiler_plot[gprofiler_plot$term_id=="GO:0044449","TERM"] <- 'contractile fiber part'
gprofiler_plot[gprofiler_plot$term_id=="GO:0044449", 'ONTOLOGY'] <- 'CC'

gprofiler_plot[gprofiler_plot$term_id=="GO:0034637","TERM"] <- 'cellular carbohydrate biosynthetic process'
gprofiler_plot[gprofiler_plot$term_id=="GO:0034637", 'ONTOLOGY'] <- 'BP'

gprofiler_plot[gprofiler_plot$term_id=="GO:0034641","TERM"] <- 'cellular nitrogen compound metabolic process'
gprofiler_plot[gprofiler_plot$term_id=="GO:0034641", 'ONTOLOGY'] <- 'BP'

gprofiler_plot[gprofiler_plot$term_id=="GO:0022838","TERM"] <- 'substrate-specific channel activity '
gprofiler_plot[gprofiler_plot$term_id=="GO:0022838", 'ONTOLOGY'] <- 'MF'

gprofiler_plot[gprofiler_plot$term_id=="GO:0000981","TERM"] <- 'DNA-binding transcription factor activity, RNA polymerase II-specific'
gprofiler_plot[gprofiler_plot$term_id=="GO:0000981", 'ONTOLOGY'] <- 'MF'

gprofiler_plot[gprofiler_plot$term_id=="GO:0044212","TERM"] <- 'transcription cis-regulatory region binding'
gprofiler_plot[gprofiler_plot$term_id=="GO:0044212", 'ONTOLOGY'] <- 'MF'


# Plot, using hatching paper

plot_graph <- ggplot(head(gprofiler_plot, 16)) + 
  geom_point(aes(negative_log10_of_adjusted_p_value, reorder(TERM, negative_log10_of_adjusted_p_value), size=intersection_size/term_size, colour=negative_log10_of_adjusted_p_value)) + 
  facet_grid(ONTOLOGY ~., scales = "free", space = "free") +
  scale_size(limits=c(0,1)) +
  scale_colour_viridis_c(direction=-1, limits=c(0,10)) +
  theme_bw() + labs(title="Vivo_Vitro_L4_Downregulated, ordered query", y="", colour= "-log10(adjp-value)", x="-log10(adjp-value)")


print(plot_graph)

#out_pdf <- file.path(out_dir, "Vivo_L2_L3_Both.pdf")

ggsave("/data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/04_analysis/R-files/Gprofiler_Differential/Gprofiler_Graphs/Vivo_Vitro_L4_Downregulated.pdf", plot = plot_graph, width = 10, height = 6)




