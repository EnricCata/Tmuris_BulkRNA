# GProfiler custom GMT generation

library(stringr)
library(dplyr)
library(readr)

eggnogTmuris <- read_tsv(
  "/path/to/GO_annotations",
  comment = "##",
  col_types = cols(.default = "c")
)
eggnogTmuris <- eggnogTmuris %>%
  rename(query = `#query`)



eggnogTmuris_GO_Kegg <- eggnogTmuris %>%
  dplyr::select(query, GOs, KEGG_ko)

# Pivot GO and KEGG_KO terms into long format
long_egnnogTmuris <- eggnogTmuris_GO_Kegg %>%
  mutate(
    GOs = strsplit(as.character(GOs), ","),
    KEGG_KOs = strsplit(as.character(KEGG_ko), ",")
  ) %>%
  pivot_longer(cols = c(GOs, KEGG_KOs), values_to = "term_list") %>%
  unnest(term_list) %>%
  mutate(term_list = str_trim(term_list)) %>%
  filter(term_list != "", !is.na(term_list))

# Group genes by term
gmt_Tmuris <- long_egnnogTmuris %>%
  group_by(term_list) %>%
  summarise(genes = list(unique(query)), .groups = "drop") %>%
  dplyr::slice(-1)
head(gmt_Tmuris)


# Convert to GMT structure
gmt_lines_Tmuris <- sapply(seq_len(nrow(gmt_Tmuris)), function(i) {
  paste(
    c(gmt_Tmuris$term_list[i], "", sort(gmt_Tmuris$genes[[i]])),
    collapse = "\t"
  )
})

# Write GMT file (tab-separated, no header)
writeLines(gmt_lines_Tmuris, "outputTmuris.gmt")



gmt_check <- read.delim("outputTmuris.gmt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Quick view
head(gmt_check)

# Check first column has no blanks
any(gmt_check$V1 == "")
# Should return FALSE

