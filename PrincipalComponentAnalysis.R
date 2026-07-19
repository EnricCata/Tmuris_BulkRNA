
library(tidyverse)
library(RColorBrewer)
library(dplyr)
library(gplots)
library(viridis)
library(stringr)
library(ggplot2)
library(ggfortify)
library(ggrepel)

# load and prepare data for PCA


data<-read.table("/path/to/kallisto_mapped_samples/kallisto_allsamples.tpm.table", header=FALSE, fill = TRUE, row.names=1, sep="\t" )
data <- data[-1,]
new_names <- readLines("/path/to/fastq-raws/SAMPLE_ID_FULL.txt")
colnames(data) <- new_names

data <- data %>% mutate(across(everything(), ~ as.numeric(as.character(.x))))
data_template <- data

data<-(data > 1) * (data - 1) + 1
data<-log10(data)
data<-as.matrix(data)
is.na(data) <- sapply(data, is.infinite)

RowVar <- function(x, ...) {
  rowSums(na.rm=TRUE,(x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}
var<-as.matrix(RowVar(data))
var <- var[order(var), ,drop = FALSE]
var_filter <- tail(var,5000)   # set number of genes here
var_filter<-as.data.frame(var_filter)
var_filter <- rownames_to_column(var_filter)
# this shows the variances!

data<-as.data.frame(data)
data<-rownames_to_column(data)

data_filtered <- dplyr::semi_join(data, var_filter, by = "rowname")
data_filtered <- column_to_rownames(data_filtered,'rowname')
data_filtered<-as.matrix(data_filtered)


## PCA with ALL SAMPLES

pca_matrix <- read.table("/lustre/scratch127/pam/teams/team333/es37/trichuris_muris/bulk-rnaseq/04_analysis/kallisto_mapped_samples/kallisto_allsamples.tpm.table", header=FALSE, fill = TRUE, row.names=1, sep="\t" )
pca_matrix <- pca_matrix[-1, ] %>% # database is fucky, deletes the empty ID row
  mutate(across(everything(), ~ as.numeric(as.character(.x)))) %>% #makes sure you have everything as numbers
  { colnames(.) <- new_names; . } %>%     # assign new column names with IDs
  filter(rowSums(across(everything()) >= 1) >= 10) %>% # filter, deletes any that have less than 10 values smaller than 1
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # adds 1 to all TPM values
  mutate(across(where(is.numeric), log2)) # log transform all TPM+1 values
gene_var <- apply(pca_matrix, 1, var) # computes variance
pca_matrix <- pca_matrix[gene_var >= quantile(gene_var, 0.9), ] # filters by the variance to choose the most relevant genes. 0.5 = 50%, 0.9 as the top 10%
pca_matrix <- pca_matrix %>%
  as.matrix() %>% #make into a large matrix and flip row-columns 
  t()
all(is.finite(pca_matrix))  #checks for empty shit and non-numerics


# Load the metadata as we're going to trim according to this
sampleinfo <- read_tsv("/path/to/tmuris_bulkrnaseq_metadata.tsv", col_types = c("cccc")) %>%
  arrange(lane_id, host, time_point, new_id)
sampleinfo$lane_id <- gsub('#', '-', sampleinfo$lane_id)
class(sampleinfo)
sampleinfo <- sampleinfo %>%
  mutate(
    prefix = sub("-.*", "", lane_id),          # keep "46825_1"
    prefix_num = as.numeric(sub("_.*", "", prefix)),  # extract 46825, 47000
    suffix_num = as.numeric(sub(".*-", "", lane_id))  # number after "-"
  ) %>%
  arrange(prefix_num, suffix_num)
sampleinfo <- as.data.frame(sampleinfo)

# To rename and remove:
  # 49707-1-2, 49707-1-1 from female to male
# 49707-1-13, 49707-1-24 from male to female
# 49707_1-3 to remove, too scattered

sampleinfo <- sampleinfo %>%
  mutate(time_point = if_else(lane_id == "49707_1-13", 'AF', time_point)) %>%
  mutate(time_point = if_else(lane_id == "49707_1-24", 'AF', time_point)) %>%
  mutate(time_point = if_else(lane_id == "49707_1-2", 'AM', time_point)) %>%
  mutate(time_point = if_else(lane_id == "49707_1-1", 'AM', time_point)) %>%
  mutate(time_point = if_else(lane_id == "49707_1-3", 'AM', time_point)) %>%
  filter(rowSums(is.na(across(everything()))) != ncol(.))

# Filter out only D33 and funky outlier which is very scattered

pca_matrix2 <- pca_matrix[
  !rownames(pca_matrix) %in% (
    sampleinfo %>%
      as_tibble() %>%                          # ensure it's a tibble
      filter(time_point == "D33" | lane_id == "49707_1-3") %>% # select rows to remove
      pull(lane_id)                            # extract lane_id
  ),
]

# Filter out D33 as they were not employed in analysis, outlier, and vitro
pca_matrix2 <- pca_matrix[
  !rownames(pca_matrix) %in% (
    sampleinfo %>%
      as_tibble() %>%                          # ensure it's a tibble
      filter(time_point == "D33" | lane_id == "49707_1-3" | host == "ORG1"| host == "ORG2") %>% # select rows to remove
      pull(lane_id)                            # extract lane_id
  ),
]





sample_pca <- prcomp(pca_matrix2)
pca_matrix[1:10, 1:5]
class(sample_pca)
str(sample_pca)

# "sdev" contains the standard deviation explained by each PC, so if we square it we get the eigenvalues (or explained variance, aka how much each component contributes to the variance. Basically, the higher the eigenvalue the 'longer' the vector is on one axis, aka the variance is larger)
pc_eigenvalues <- sample_pca$sdev^2

# "rotation" contains the variable loadings for each PC, which define the eigenvectors, aka the axis of the PCA

pc_loadings <- sample_pca$rotation

# "x" contains the PC scores, i.e. the data projected on the new PC axis
# "center" in this case contains the mean of each gene, which was subtracted from each value
# "scale" contains the value FALSE because we did not scale the data by the standard deviation


length(pc_eigenvalues)
ncol(pc_loadings)

#PC scores
pc_scores <- sample_pca$x

#Variance explained by PCs
#Convert matrix to tibble object in order to make a plot in ggplot2
# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

#Print the result
pc_eigenvalues

#Produce a Scree Plot to show the fraction of the total variance explained by each PC
pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  coord_cartesian(xlim = c(1, 20)) +
  labs(x = "Principal component", y = "Fraction variance explained")

# Function that makes the label 'PC(X) (variance percentage)
make_pc_label <- function(df, pc_num) {
  pct <- df$pct[df$PC == pc_num]
  sprintf("PC%d (%.1f%%)", pc_num, pct)
}



# The PC scores are stored in the "x" value of the prcomp object
# pc_scores <- sample_pca$x (you did this before)
pc_scores <- pc_scores %>% 
  #Convert to a tibble retaining the sample names as a new column! That way you keep the columns
  as_tibble(rownames = "sample")

pc_scores

# --------------------------------   Sample Info  ----------------------------------------------------------

sampleinfo <- sampleinfo %>%
  mutate(
    Model = dplyr::case_when(
      host %in% c("ORG1", "ORG2") ~ "In-vitro",
      host %in% c("MOU")         ~ "In-vivo",
      host %in% c("IVI")         ~ "L1 Hatched in-vitro",
      TRUE                       ~ host
    )
  )

sampleinfo <- sampleinfo %>%
mutate(Life_stage = case_when(
  time_point == "L1" ~ "L1 - In vitro hatched",
  time_point %in% c("D13", "D14", "D15") ~ "L2",
  time_point == "D20" ~ "L3",
  time_point %in% c("D24", "D25") ~ "L4",
  time_point == "AM" ~ "Adult Male",
  time_point == "AF" ~ "Adult Female",
  TRUE ~ NA_character_
))

sampleinfo <- sampleinfo %>%
  filter(time_point != "D33")

# --------------------------------  PC Plots  ----------------------------------------------------------

#Labels!
pc1_label <- make_pc_label(pc_eigenvalues, 1)
pc2_label <- make_pc_label(pc_eigenvalues, 2)
pc3_label <- make_pc_label(pc_eigenvalues, 3)

# If one of the PCs is flipped, do this. 
# It doesn't matter if its negative or positive, as its simply a measure of variance
# sample_pca$x[, "PC2"] <- -sample_pca$x[, "PC2"] 
# Remove and rearrange the axis (PC1, PC2, PC3) as needed

pca_plotPC1PC2 <- sample_pca$x %>%
  as_tibble(rownames = "sample") %>%
  left_join(sampleinfo, by = c("sample" = "lane_id")) %>%
  mutate(
    Life_stage = factor(
      Life_stage,
      levels = c(
        "L1 - In vitro hatched",
        "L2",
        "L3",
        "L4",
        "Adult Female",
        "Adult Male"
      )
    )
  ) %>%
  
  
  #Create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = Life_stage, shape = Model), size = 3) +
#  geom_text_repel(aes(label = sample), size = 3) +
  scale_shape_manual(
    values = c(
      "In-vitro" = 17,
      "L1 Hatched in-vitro" = 16,
      "In-vivo" = 15
    )
  ) +
  
  # ---- CUSTOM COLOR SCALE: viridis except AF = #FFC0CB ----
scale_colour_manual(
  values = {
    tp <- unique(sampleinfo$Life_stage)                   # time_point levels
    cols <- viridis_pal(option = "D")(length(tp))         # viridis colors
    names(cols) <- tp                                     # name them
    cols["Adult Female"] <- "pink"                               # override AF to pink
    cols["L1 - In vitro hatched"] <- "#1b9e77"
    cols["Adult Male"] <- "cyan"
    cols["L2"] <- "salmon"                               # override AF to pink
    cols["L3"] <- "#87d3f8"
    cols["L4"] <- "purple"
    cols
  }
) +
  #geom_text_repel(aes("sample"), point.padding=0.8, size=3, box.padding = 0.3, min.segment.length=1, max.overlaps=Inf) +
  theme (panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         panel.background=element_blank(),
         legend.key = element_blank(),
         axis.line = element_line(colour = "gray"),
         legend.title = element_text(size = 14),
         legend.text = element_text(size = 12)) +
  scale_alpha(guide = 'none') +
  labs(
    x = pc2_label,
    y = pc3_label,
#    title = "PCA top 10% most variable genes Vivo vs Vitro, PC1 vs PC2"
  )


pca_plotPC1PC2
