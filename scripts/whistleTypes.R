# Analysis of whistle types from ARTwarp categorization

#------------------------------------------------------------------------------

# Load necessary libraries
library(ggplot2)# For plotting
library(dplyr)    # For data manipulation
library(patchwork)
library(RColorBrewer)
library(stringr)
library(tidyr)
library(reshape2)


# Set working directory to the path containing the data file
setwd("Z:/MaiaProjects/manzanilloanalysis/scripts")

# Read whistle types data from a CSV file
category_types <- read.csv("whistletypes.csv")

dt <- read.csv("Z:/MaiaProjects/manzanilloanalysis/contours_update/msg_full_it_167.csv")

names(category_types)[names(category_types) == "Group"] <- "GroupName"

# Merge the data frames on 'category'
whistletypes <- merge(dt, category_types, by = "category", all.x = TRUE)

# View the result
head(whistletypes)


# Calculate the occurrence and frequency of each whistle type by group
whistletypeoccurence <- whistletypes %>%
  group_by(type, GroupName) %>%
  summarize(count = n()) %>%              # Count occurrences of each type
  group_by(GroupName) %>%
  mutate(freq = count / sum(count) * 100) # Calculate frequency as a percentage of total

whistletypeoccurence <- whistletypeoccurence %>%
  filter(!is.na(type), !is.na(GroupName))

# Plotting
pinks <- c("#ffe8ee", "#ffd6e0", "#feb4c7", "#fe7799","#fe0849", "#710a25")
blues   <- brewer.pal(6, "Blues")
purples <- brewer.pal(6, "Purples")
greens  <- brewer.pal(6, "Greens")


make_group_plot <- function(df, group_label, palette) {
  library(stringr)
  
  # Wrap the group label before filtering data
  group_label_wrapped <- str_replace(group_label, regex("^((?:\\S+\\s){2})"), "\\1\n")
  
  df_group <- df %>% 
    filter(trimws(GroupName) == group_label) %>%
    mutate(type = factor(type, levels = unique(type)))
  
  ggplot(df_group, aes(x = type, y = freq, fill = type)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(round(freq), "%")), 
              vjust = -0.75, size = 4) +
    scale_fill_manual(values = palette) +
    labs(title = group_label_wrapped, x = NULL, y = "Abundance (%)") +
    coord_cartesian(ylim = c(0, 100)) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.background = element_blank()
    )
}

p1 <- make_group_plot(whistletypeoccurence, "Guiana Dolphin Single Species Group (Exclusive)", pinks)
p2 <- make_group_plot(whistletypeoccurence, "Bottlenose Dolphin Single Species Group (Exclusive)", blues)
p3 <- make_group_plot(whistletypeoccurence, "Mixed Species Group (Exclusive)", purples)
p4 <- make_group_plot(whistletypeoccurence, "Shared Across 2 or More Contexts", greens)

combined_plot <- (p1 | p2) / (p3 | p4) +
  plot_annotation(tag_levels = 'a')

ggsave("whistle_types_fig.jpg", combined_plot, dpi=300)


##### PERMUTATION ####

# Keep only relevant columns
df_subset <- whistletypes[, c("type", "GroupName")]

# Observed chi-squared statistic
obs_table <- table(df_subset$type, df_subset$GroupName)
obs_chisq <- chisq.test(obs_table)$statistic

# Permutation test setup
set.seed(42)  # for reproducibility
n_perm <- 1000
perm_stats <- numeric(n_perm)

# Permutation loop
for (i in 1:n_perm) {
  shuffled <- sample(df_subset$GroupName)  # shuffle group labels
  perm_table <- table(df_subset$type, shuffled)
  perm_stats[i] <- suppressWarnings(chisq.test(perm_table)$statistic)
}

# Calculate p-value
p_val <- mean(perm_stats >= obs_chisq)

# Output result
cat("Observed chi-squared:", obs_chisq, "\n")
cat("Permutation p-value:", p_val, "\n")

# Create contingency table of counts of each Type by GroupName
contingency_table <- table(whistletypes$type, whistletypes$GroupName)

# Run chi-squared test and get standardized residuals
chisq_res <- chisq.test(contingency_table)

# Standardized residuals matrix
std_residuals <- chisq_res$stdres

# Overall effect size (Cramer's V)
cramers_v <- sqrt(chisq_res$statistic / (sum(contingency_table) * (min(dim(contingency_table)) - 1)))

# Print effect size
cat("Cramér's V:", cramers_v, "\n")


residuals_df <- melt(std_residuals)
colnames(residuals_df) <- c("Type", "GroupName", "StdResidual")

