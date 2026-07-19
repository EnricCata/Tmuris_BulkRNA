# -*- coding: utf-8 -*-
library(sleuth)
library(tidyverse)
library(stringr)
library(ggplot2)


so <- readRDS("/path/to/object_so.rds", refhook = NULL)
so_prep <- readRDS("/path/to/object_so_prep.rds")

#To check successful LRT tests
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
head(sleuth_table, 20)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20)

counts_sleuth <- sleuth_to_matrix(so, "est_counts")

sleuth_significant_top20 <- head(sleuth_significant, 20)
sleuth_significant_top20


# +
# UpSet Plot

#Preparing the upset matrix

# No comparison, only sleuth prep as its a collection of normalized Sleuth objects 
library(dplyr)
so_matrix <- sleuth_to_matrix(so_prep, "obs_norm", "tpm")
head(so_matrix, 20)
# -

# Metadata for comparing in upset
larval_stage <- sampleinfo_trim$larval_stage
larval_stage

# turns the matrix into a dataframe, also adds metadata (the larval stage essentially)
expr_stage <- so_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  left_join(sampleinfo_trim, by = "sample") %>%
  group_by(gene, larval_stage) %>%
  summarise(mean_expr = mean(expr), .groups = "drop")
#Note: the last one is a mean of expression so single replicates don't screw over general presence. May cause error?

head(expr_stage, 100)

# +
# Checks if there's low or negative expression of certain genes! Essentially classify presence or absence... 
# might need to rework, maybe not stringent enough

expr_stage_collapsed <- expr_stage %>%
  group_by(gene, larval_stage) %>%
  summarise(
    expressed = any(mean_expr >= 1),  # TRUE if any row exceeds threshold
    .groups = "drop"
  )
head(expr_stage_collapsed, 20)

# +
# Trims any transcript thats expressed below, say, 2000 TPM

expr_stage_collapsed <- expr_stage %>%
  group_by(gene, larval_stage) %>%
  summarise(
    expressed = any(mean_expr >= 1000),  # TRUE if any row exceeds threshold
    .groups = "drop"
  )
head(expr_stage_collapsed, 20)

# +
length(unique(expr_stage$gene))
length(unique(expr_stage_collapsed$gene))

expr_stage_no_expression <- expr_stage_collapsed %>%
  group_by(gene) %>%
  summarise(all_false = all(expressed == FALSE)) %>%
  filter(all_false) 


head(expr_stage_no_expression, 20)
length(unique(expr_stage_no_expression$gene))

# +
# Collapses into a matrix that can work
expr_binary <- expr_stage_collapsed %>%
pivot_wider(
  names_from = larval_stage,
  values_from = expressed,
)
# expr_binary %>%
    # relocate('L5-YA', .before = 'L5-AF')
  
head(expr_binary, 100)

# -
# Saving it due to running it in R client eats all variables but L3 and L1, but not on Jupyter?
message("Saving Upset Binary ...")
saveRDS(expr_binary, file = "/data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/04_analysis/R-files/UpsetData/RSleuthobject_Vivo_Adult_Label_Upset_Binary.rds")
message("Done! Sleuth object saved successfully.")


expr_binary <- readRDS("/data/pam/team333/es37/scratch/trichuris_muris/bulk-rnaseq/04_analysis/R-files/UpsetData/RSleuthobject_Vivo_Adult_Label_Upset_Binary.rds", refhook = NULL)

# Ensures everything is in numeric
expr_binary_numeric <- expr_binary %>%
  mutate(across(-gene, as.integer))
head(expr_binary_numeric)

str(expr_binary_numeric)

# Transforms back into dataframe for upset
expr_binary_numeric_df <- as.data.frame(expr_binary_numeric)
# expr_binary_numeric_df %>%
    # relocate('L5-YA', .before = 'L5-AF')
str(expr_binary_numeric_df)

# -

# Plotting the Upset Plot

library(ComplexUpset)
library(UpSetR)

# To trim for the shared one to later do the scale break

expr <- apply(expr_binary_numeric_df, 1, function(row) {
  paste(names(row)[row == 1], collapse = "&")
})

expr_counts <- table(expr)
expr_counts_capped <- pmin(expr_counts, 1000)
data_exp <- fromExpression(expr_counts_capped)

upset(
 expr_binary_numeric_df,
#  data_exp,
  sets = colnames(expr_binary_numeric)[-1],
  nsets = length(colnames(expr_binary_numeric_df)[-1]),
  nintersects = NA,
  order.by = c("freq","degree"),
#  scale.intersections = "log10",
  keep.order = TRUE,
)


sets_to_plot <- c("L1", "L2", "L3", "L4", "Adult Female", "Adult Male")

expr_binary <- expr_binary %>%
dplyr::rename(
  "Adult Female" = 'L5-AF',
  "Adult Male" = 'L5-AM')


# Adds a logarithmic scale
presence = ComplexUpset:::get_mode_presence('exclusive_intersection')
summarise_values = function(df) {
  aggregate(
    as.formula(paste0(presence, '~ intersection')),
    df,
    FUN=sum
  )
}

expr_binary[, sets_to_plot] <- lapply(
  expr_binary[, sets_to_plot],
  function(x) {
    x <- as.numeric(x)
    x == 1
  }
)

expr_binary <- expr_binary[, sets_to_plot]

expr_binary$.pattern <- apply(expr_binary[, sets_to_plot], 1, function(row) {
  sets <- sets_to_plot[as.logical(row)]
  paste(length(sets), paste(sets, collapse = "+"))
})



# Upsetplot for chosen comparisons
expr_filtered <- expr_binary[expr_binary$.pattern %in% c( 
  "1 L1", 
  "1 L2", 
  "1 L3", 
  "1 L4",
  "1 Adult Female",
  "1 Adult Male",
  "2 Adult Female+Adult Male",
  "3 L2+L3+L4",
  "4 L1+L2+L3+L4",
  "6 L1+L2+L3+L4+Adult Female+Adult Male"),]


ComplexUpset::upset(
  expr_filtered,
  intersect = sets_to_plot,
  sort_sets = FALSE,
  set_sizes = FALSE,
  sort_intersections = FALSE,
  intersections=list(
    'L1',
    'L2',    
    'L3',
    'L4',
    'Adult Female',
    'Adult Male',
    c('Adult Female', 'Adult Male'),
    c('L2', 'L3', 'L4'),
    c('L1', 'L2', 'L3', 'L4'),
    c('L1', 'L2', 'L3', 'L4','Adult Female', 'Adult Male')
    
  ),
  base_annotations = list(
    'log10(intersection size)' = (
      ggplot() +
        geom_bar(
          data = summarise_values,
          stat = 'identity',
          aes(x = intersection, y = !!presence)
        ) +
        geom_text(
          data = summarise_values,
          aes(x = intersection, y = !!presence, label = !!presence),
          vjust = 1.5,
          color = "white",
          size = 4
        ) +
        ylab('Intersection size (Log10)') +
        scale_y_continuous(trans = 'log10')
    )
),
min_size=3,
width_ratio=0.1
)

# Upsetplot for ALL comparisons
ComplexUpset::upset(
  expr_binary,
  sort_sets = FALSE,
  intersect = sets_to_plot,
  min_size = 3,
  width_ratio = 0.1,   # <- THIS is the key
  base_annotations = list(
    'log10(intersection size)' = (
      ggplot() +
        geom_bar(
          data = summarise_values,
          stat = 'identity',
          aes(x = intersection, y = !!presence)
        ) +
        geom_text(
          data = summarise_values,
          aes(x = intersection, y = !!presence, label = !!presence),
          angle = 90,
          vjust = 0.5,
          hjust = 1.5,
          color = "white",
          size = 3.5
        ) +
        ylab('Intersection size (Log10)') +
        scale_y_continuous(trans = 'log10')
    )
  ),
  sort_intersections_by = c('degree'),
) +
  theme(
    # Intersection size labels (top bars)
    axis.text.y = element_text(size = 6),        # smaller y-axis numbers

    # Group labels (column names)
#    axis.text.x = element_text(size = 2, angle = 60, hjust = 1),  
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    # Intersection dots
    upset_points_size = 1                         # smaller circles
  )





# +

# Heatmap


# For vivo only, load the so_prep which excludes in-vitro results

# For vivo vs vitro, load the so_prep which includes all comparisons 
# +


# Create stage annotation
annotation <- data.frame(
  stage = so_prep$sample_to_covariates$larval_stage,
  host = so_prep$sample_to_covariates$host,
  time = so_prep$sample_to_covariates$time_point
)


# Obtain normalized TMP matrix
obs_norm_filt <- so_prep$obs_norm_filt

tpm_norm_matrix <- tidyr::pivot_wider(
  obs_norm_filt,
  id_cols = target_id,
  names_from = sample,
  values_from = tpm
)


# Set rownames to target_id and remove the column
rownames(tpm_norm_matrix) <- tpm_norm_matrix$target_id
tpm_norm_matrix$target_id <- NULL

# Filter adults out if need be
tpm_norm_matrix <- tpm_norm_matrix[, !grepl("49707", colnames(tpm_norm_matrix))]

# Filter larval stages if need be
tpm_norm_matrix <- tpm_norm_matrix[, !grepl("46825", colnames(tpm_norm_matrix))]

# Filter L1 stages if need be
tpm_norm_matrix <- tpm_norm_matrix[, !grepl("46825_1-22", colnames(tpm_norm_matrix))]
tpm_norm_matrix <- tpm_norm_matrix[, !grepl("46825_1-23", colnames(tpm_norm_matrix))]
tpm_norm_matrix <- tpm_norm_matrix[, !grepl("46825_1-24", colnames(tpm_norm_matrix))]

#Log transform
log_tpm <- log2(tpm_norm_matrix + 1)
# Remove low expression transcripts
log_tpm <- log_tpm[rowMeans(log_tpm) > 1, ]
# Select top 1000 most variable transcripts
vars <- apply(log_tpm, 1, var)
top_genes <- names(sort(vars, decreasing = TRUE))[1:1000]
log_tpm_top <- log_tpm[top_genes, ]

# Compute sample correlations
cor_matrix <- cor(log_tpm_top, method = "pearson")


annotation$sample <- ifelse(
  grepl("49707", rownames(annotation)),
  paste0("2", sub(".*49707", "", rownames(annotation))),
  ifelse(
    grepl("46825", rownames(annotation)),
    paste0("1", sub(".*46825", "", rownames(annotation))),
    NA
  )
)

# For annotation's sake of vitro v vivo
annotation$host3 <- ifelse(
  annotation$host %in% c("ORG1", "ORG2"),
  paste0("ORG", annotation$stage),
  ifelse(
    annotation$host == "IVI",
    annotation$stage,
    paste0(annotation$host, annotation$stage)
  )
)

annotation$host[annotation$host %in% c("ORG1", "ORG2")] <- "ORG"
annotation$host[annotation$host == "IVI"] <- "L1"


# Manually remove 497077_1-3 as its a mislabelled outlier

cor_matrix <- cor_matrix[!(rownames(cor_matrix) %in% "49707_1-3"), !(rownames(cor_matrix) %in% "49707_1-3")]
annotation <- annotation[rownames(annotation) != '49707_1-3', , drop = FALSE]

# Remove adults if need be
annotation <- annotation[!grepl("49707", rownames(annotation)), , drop = FALSE]

# Remove larva if need be
annotation <- annotation[!grepl("46825", rownames(annotation)), , drop = FALSE]

# Remove L1 larva if need be
annotation <- annotation[!grepl("46825_1-22", rownames(annotation)), , drop = FALSE]
annotation <- annotation[!grepl("46825_1-23", rownames(annotation)), , drop = FALSE]
annotation <- annotation[!grepl("46825_1-24", rownames(annotation)), , drop = FALSE]
 


# Filter annotation 

annotation<- annotation %>%
  filter(!stage == 'MOU')


ann_colors <- list(
  stage = c(
    "L1" = "#1b9e77",
    "L2" = "salmon",
    "L3" = "#87d3f8",
    "L4" = "purple" ,
   "Adult Male" = "cyan",
    "Adult Female" = 'pink'
  ),
  host = c(
    'ORG' = "#c7322a",
#    'ORG2' = '#c7792a',
    'MOU' = '#27b7db',
    'L1'= '#14c90a'
#  ),
#  host2 = c(
#    "L1" = "#1b9e77",
#    "MOUL2" = "salmon",
#    "ORGL2" = "#cf695d",
#    "MOUL3" = "#87d3f8",
#    "ORGL3" = "#6ca7c4",
#    "MOUL4" = "purple" ,
#    "ORGL4" = "#7f17bf" 
  )
)

annotation<- annotation %>%
dplyr::mutate(
  stage = dplyr::recode(
    stage,
    "L5-AM" = "Adult Male",
    "L5-AF" = "Adult Female"
  )
)

annotation$stage <- factor(
  annotation$stage,
  levels = c("L1", "L2", "L3", "L4", "Adult Male", "Adult Female")
)

library(ComplexHeatmap)

dend <- as.dendrogram(hclust(dist(cor_matrix), method = "complete"))

# For Vivo to put the L1 at the top, L2 at the left, L3 & L4 together, young'ns to the side:
dend <- rev(dend)
dend[[1]][[1]] <- rev(dend[[1]][[1]])
dend[[1]][[1]][[2]] <- rev(dend[[1]][[1]][[2]])
dend[[1]][[1]][[2]][[1]] <- rev(dend[[1]][[1]][[2]][[1]])
dend[[1]][[1]][[2]][[1]][[1]][[2]] <- rev(dend[[1]][[1]][[2]][[1]][[1]][[2]])
dend[[1]][[1]][[2]][[1]][[2]][[2]] <- rev(dend[[1]][[1]][[2]][[1]][[2]][[2]])
dend[[1]][[1]][[2]][[1]][[2]] <- rev(dend[[1]][[1]][[2]][[1]][[2]])
dend[[1]][[2]] <- rev(dend[[1]][[2]])
dend[[2]] <- rev(dend[[2]])


#For Vivo no adults to put the L2 to the side and the L3/L4 to the other, making the replicates consistent too
dend[[1]] <- rev(dend[[1]])
dend[[2]] <- rev(dend[[2]])
dend[[2]][[1]] <- rev(dend[[2]][[1]])
dend[[2]][[2]] <- rev(dend[[2]][[2]])

dend[[2]][[2]][[2]] <- rev(dend[[2]][[2]][[2]])

#For Vivo adults only, being consistent with the previous graph
dend <- rev(dend)
dend[[1]] <- rev(dend[[1]])
dend[[1]][[2]] <- rev(dend[[1]][[2]])
dend[[1]][[2]][[1]] <- rev(dend[[1]][[2]][[1]])
dend[[1]][[2]][[2]] <- rev(dend[[1]][[2]][[2]])
dend[[2]] <- rev(dend[[2]])
dend[[2]][[1]][[1]] <- rev(dend[[2]][[1]][[1]])
dend[[2]][[1]][[1]][[1]] <- rev(dend[[2]][[1]][[1]][[1]])


# For Vivo vs Vitro to put the L1 at the top and the adults at the bottom
dend <- rev(dend)
dend[[1]] <- rev(dend[[1]])
dend[[1]][[1]] <- rev(dend[[1]][[1]])

# For Vivo vs Vitro to put the L2 vivo next to the L3s
dend[[1]][[1]][[2]][[2]] <-rev(dend[[1]][[1]][[2]][[2]])

# For Vivo vs Vitro to put the vivo L3 away from the vitros
dend[[1]][[2]][[1]][[1]] <-rev(dend[[1]][[2]][[1]][[1]])

# For Vitro vs Vivo, without any adults to mimic as much as possible
dend <- rev(dend)
dend[[1]] <- rev(dend[[1]])
dend[[2]] <- rev(dend[[2]])
dend[[1]][[2]] <- rev(dend[[1]][[2]])
dend[[2]][[2]][[2]][[2]][[2]] <- rev(dend[[2]][[2]][[2]][[2]][[2]])
dend[[1]][[2]][[1]] <- rev(dend[[1]][[2]][[1]])
dend[[1]][[2]][[1]][[2]] <- rev(dend[[1]][[2]][[1]][[2]])
dend[[2]][[2]][[1]] <- rev(dend[[2]][[2]][[1]])

#L2 vivo
dend[[2]][[1]][[2]] <- rev(dend[[2]][[1]][[2]])

#l3 vivo
dend[[2]][[2]][[2]][[2]][[2]][[1]] <- rev(dend[[2]][[2]][[2]][[2]][[2]][[1]])

#L4 vivo
dend[[2]][[2]][[2]][[2]][[2]][[2]] <- rev(dend[[2]][[2]][[2]][[2]][[2]][[2]])
dend[[2]][[2]][[2]][[2]][[2]][[2]][[2]] <- rev(dend[[2]][[2]][[2]][[2]][[2]][[2]][[2]])

Heatmap_Object <- Heatmap(name = "correlation",
  cor_matrix,
  cluster_rows = dend,
  cluster_columns = dend,
#  row_labels = paste(annotation$time),
show_row_names = FALSE,
  column_labels = paste(annotation$sample),
  top_annotation = HeatmapAnnotation(stage = annotation$stage, col = ann_colors),
#  top_annotation = HeatmapAnnotation(host2 = annotation$host2, col = ann_colors),
  left_annotation = rowAnnotation(host = annotation$host, col = ann_colors),
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  rect_gp = grid::gpar(col = "grey", lwd = 0.5),
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 8),
#  column_title = "Sample correlation heatmap  (Top 1000 variable transcripts), Vivo vs Vitro",
) 
  draw(
    Heatmap_Object,
    merge_legend = TRUE,
    annotation_legend_side = "right"
  )



#  annotation_col = annotation_ordered,
#  annotation_colors = ann_colors,
#  clustering_distance_rows = "euclidean",
#  clustering_distance_cols = "euclidean",
#  main = "Sample correlation heatmap (Top variable transcripts), Vivo"
#

# +
  #Heatmap Correlation Boxplot


# Obtain correlation averages


correlation_averages <- lapply(rownames(annotation)[annotation$host %in% c("ORG1", "ORG2")], function(samp) {
  #the above long-ass line picks ONLY the organoid samples, around 21 of em
  stage <- annotation[samp, "stage"]
  # MOU samples of same stage
  mou_stage <- rownames(annotation)[
    annotation$host == "MOU" & annotation$stage == stage
  ]
  # correlations for this sample vs those MOU samples
  vals <- cor_matrix[samp, mou_stage]
  data.frame(
    sample = samp,
    stage = stage,
    mean_cor = mean(vals, na.rm = TRUE)
  )
}) %>%
  bind_rows()


correlation_all <- lapply(
  rownames(annotation)[annotation$host %in% c("ORG1", "ORG2")],
  function(samp) {
#the above long-ass line picks ONLY the organoid samples, around 21 of em
    stage <- annotation[samp, "stage"]
# MOU samples of same stage
    mou_stage <- rownames(annotation)[
      annotation$host == "MOU" &
        annotation$stage == stage]
    vals <- cor_matrix[samp, mou_stage]
    data.frame(
      sample = samp,
      mou_sample = mou_stage,
      stage = stage,
      correlation = vals
    )
  }
) %>%
  bind_rows()

#Plot as boxplot

library(viridis)

ggplot(correlation_all, aes(x = stage, y = correlation)) +
  geom_boxplot(outlier.shape = NA, aes(fill = stage)) +
  scale_fill_manual(values = c(
    "L2" = "salmon",
    "L3" = "#87d3f8",
    "L4" = "purple"
  )) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  labs(
    x = "Larval stage",
    y = "Mean correlation to mouse in-vivo samples"
  ) + theme(legend.position = "none")

# +
