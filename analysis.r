
library(tidyverse)

find_file <- function(fname, search_paths = c(".", "data", "Data")) {
  for (p in search_paths) {
    candidate <- file.path(p, fname)
    if (file.exists(candidate)) return(candidate)
  }
  

  stop(
    paste0(
      "Could not find file: ", fname,
      "\nSearched in: ", paste(search_paths, collapse = ", "),
      "\nFiles in current directory: ", paste(list.files("."), collapse = ", ")
    )
  )
}

expr_fname <- "GSE60450_GeneLevel_Normalized(CPM.and.TMM)_data.csv"
meta_fname <- "GSE60450_filtered_metadata.csv"

expr_path <- find_file(expr_fname)
meta_path <- find_file(meta_fname)

expr <- read_csv(expr_path, show_col_types = FALSE)
meta <- read_csv(meta_path, show_col_types = FALSE)

dim(expr)
head(expr)

dim(meta)
head(meta)

meta2 <- meta %>%
  rename(sample_id = 1)

expr_long <- expr %>%
  rename(ensembl_id = 1) %>%
  pivot_longer(
    cols = starts_with("GSM"),
    names_to = "sample_id",
    values_to = "expression"
  )

dat <- expr_long %>%
  left_join(meta2, by = "sample_id")

n_distinct(dat$sample_id)

if (all(is.na(dat$immunophenotype))) {
  stop("Join failed: immunophenotype is NA for all rows. Check sample_id matching between expression and metadata.")
}

dat %>% count(immunophenotype, sort = TRUE)

gene_summary <- dat %>%
  group_by(ensembl_id, gene_symbol) %>%
  summarise(mean_expression = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_expression))

top_gene <- gene_summary %>% slice(1)

top_gene_id <- top_gene$ensembl_id[1]
top_gene_symbol <- top_gene$gene_symbol[1]

message("Selected gene for plotting: ", top_gene_symbol, " (", top_gene_id, ")")

dat_one_gene <- dat %>%
  filter(ensembl_id == top_gene_id)

p <- ggplot(dat_one_gene, aes(x = immunophenotype, y = expression)) +
  geom_boxplot() +
  labs(
    title = paste0("Expression by cell type: ", top_gene_symbol, " (", top_gene_id, ")"),
    x = "Cell type (immunophenotype)",
    y = "Expression (CPM / TMM)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(p)

if (!dir.exists("results")) dir.create("results")

ggsave(
  filename = "results/expression_by_cell_type.png",
  plot = p,
  width = 9,
  height = 5,
  dpi = 300
)
