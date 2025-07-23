library(caret)
library(randomForest)
library(dplyr)
library(ggplot2)
library(corrplot)
library(janitor)
library(pdp)
library(gridExtra)
library(stringr)




df <- read.csv("Z:/MaiaProjects/manzanilloanalysis/scripts/random_forest_parameters.csv")

df <- df %>% janitor::clean_names()

# After cleaning, get all columns except 'species', 'site', etc.
non_feature_cols <- c("species", "site", "group")
feature_cols <- setdiff(names(df), non_feature_cols)


train_df <- df %>%
  filter(species %in% c("B", "G"), site %in% c("P", "C", "R")) %>%
  droplevels()

# Convert all feature columns to numeric
train_df[feature_cols] <- lapply(train_df[feature_cols], function(x) as.numeric(as.character(x)))

# Find columns with zero variance
zero_var_cols <- feature_cols[sapply(train_df[, feature_cols], function(x) sd(x, na.rm = TRUE) == 0)]
zero_var_cols
# Drop them from the feature list
feature_cols <- setdiff(feature_cols, zero_var_cols)


# Check correlations
cor_matrix <- cor(train_df[, feature_cols], use = "complete.obs")
corrplot(cor_matrix, method = "color", tl.cex = 0.7)

# Remove highly correlated features (> 0.9)
high_cor <- findCorrelation(cor_matrix, cutoff = 0.9, names = TRUE)
reduced_features <- setdiff(feature_cols, high_cor)

set.seed(123)

train_model_df <- train_df[, c("species", reduced_features)] %>%
  na.omit()

# Define trainControl with 5-fold CV
ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE)


# Train model with tuning for mtry
rf_tuned <- train(
  species ~ .,
  data = train_model_df[, c("species", reduced_features)],
  method = "rf",
  trControl = ctrl,
  tuneLength = 5,
  importance = TRUE
)

print(rf_tuned)

plot(rf_tuned)



# True mix: All samples from site C where species is NA (no labels)
true_mix <- df %>%
  filter(site == "C", is.na(species))

# Artificial mix: combine equal samples from B@P and G@R
n_mix <- nrow(true_mix) / 2
set.seed(123)

# Sample equal numbers from B@C and G@C
b_p <- train_df %>% filter(species == "B", site == "P") %>% slice_sample(n = n_mix)
g_r <- train_df %>% filter(species == "G", site == "R") %>% slice_sample(n = n_mix)

# Combine the two species into the artificial mix
artif_mix <- bind_rows(b_p, g_r)

true_mix[, reduced_features] <- lapply(true_mix[, reduced_features], function(x) as.numeric(as.character(x)))

true_probs <- predict(rf_tuned, newdata = true_mix[, reduced_features], type = "prob") %>%
  mutate(source = "True Mix")

artif_probs <- predict(rf_tuned, newdata = artif_mix[, reduced_features], type = "prob") %>%
  mutate(source = "Artificial Mix")

plot_df <- bind_rows(true_probs, artif_probs)

ggplot(plot_df, aes(x = B, fill = source)) +
  geom_density(alpha = 0.5) +
  labs(title = "Predicted Probability for Class 'B'",
       x = "Probability of B", y = "Density") +
  theme_minimal()

ggplot(plot_df, aes(x = G, fill = source)) +
  geom_density(alpha = 0.5) +
  labs(title = "Predicted Probability for Class 'G'",
       x = "Probability of G", y = "Density") +
  theme_minimal()




# Combine the true and artificial probabilities into one data frame
all_probs <- bind_rows(true_probs, artif_probs)

# Tabulate the results
# The columns 'B' and 'G' represent the predicted probabilities for each species
prob_table <- all_probs %>%
  group_by(source) %>%
  summarise(
    mean_B = mean(B, na.rm = TRUE),  # Mean probability for 'B'
    sd_B = sd(B, na.rm = TRUE),      # Standard deviation for 'B'
    mean_G = mean(G, na.rm = TRUE),  # Mean probability for 'G'
    sd_G = sd(G, na.rm = TRUE)       # Standard deviation for 'G'
  )


print(prob_table)

# Add species label manually before combining
true_probs <- predict(rf_tuned, newdata = true_mix[, reduced_features], type = "prob") %>%
  mutate(source = "True Mix", species = "unknown")

artif_probs <- predict(rf_tuned, newdata = artif_mix[, reduced_features], type = "prob") %>%
  mutate(source = "Artificial Mix", 
         species = rep(c("B", "G"), each = nrow(artif_mix) / 2))

# Combine them
all_probs <- bind_rows(true_probs, artif_probs)

# Summarize by source and species
prob_table <- all_probs %>%
  group_by(source, species) %>%
  summarise(
    mean_B = mean(B, na.rm = TRUE),
    sd_B = sd(B, na.rm = TRUE),
    mean_G = mean(G, na.rm = TRUE),
    sd_G = sd(G, na.rm = TRUE),
    .groups = "drop"
  )

print(prob_table)



## Variable importance
final_rf <- rf_tuned$finalModel

varImpPlot(final_rf)

# Get the variable importance from the random forest model
importance_obj <- varImp(rf_tuned, scale = FALSE)
importance_df <- importance_obj$importance  # Extract the actual data frame

print(head(importance_df))

# If it has per-class columns (e.g., "B" and "G"), average them:
if (!"Overall" %in% colnames(importance_df)) {
  importance_df$Overall <- rowMeans(importance_df)
}

# Now extract the top 5 variables
top_vars <- rownames(importance_df[order(-importance_df$Overall), ])[1:4]
print(top_vars)  # confirm it's not empty


# Function to clean variable names (title case and spaces)
clean_var_name <- function(var) {
  var %>%
    str_replace_all("_", " ") %>%
    str_to_title()
}

# Generate PDP plots
pdp_plots <- lapply(top_vars, function(var) {
  pdp_data <- as.data.frame(partial(
    rf_tuned$finalModel,
    pred.var = var,
    train = train_model_df[, reduced_features],
    prob = TRUE
  ))
  
  ggplot(pdp_data, aes_string(x = var, y = "yhat", color = "yhat")) +
    geom_path(size = 1.5) +
    scale_color_gradient(low = "#fe7799", high = "#6BAED6", guide = "none") +
    labs(
      title = clean_var_name(var),
      x = clean_var_name(var),
      y = "Probability (Bottlenose)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(size = 13)
    )
})
library(cowplot)

pdp_grid <- plot_grid(plotlist = pdp_plots, 
                      labels = letters[1:length(pdp_plots)], 
                      label_size = 18, 
                      ncol = 2)

ggsave("partial_dependence_plots.jpg", pdp_grid, dpi = 300)

