# Load necessary libraries for plotting and data manipulation
library(ggplot2)
library(ggfortify)
library(ggsci)
library(RColorBrewer) # For color scales
library(vegan)
library(dplyr)

# Set working directory to where the data files are located
setwd("/Users/maustin/Desktop/manzanillo_analysis")

# Read data from CSV file
res <- read.csv('element_mds_pc.csv')

pc <- data.frame(res)
names(pc)[4] <- 'PC1'
names(pc)[5] <- 'PC2'
names(pc)[2] <- 'Group'

# Palette
pinks <- c("#ffe8ee", "#ffd6e0", "#feb4c7", "#fe7799", "#fe0849", "#710a25")
blues <- RColorBrewer::brewer.pal(6, "Blues")
purples <- RColorBrewer::brewer.pal(6, "Purples")


fill_palette <- c(
  "MSG" = purples[4],
  "Guiana - Brazil" = pinks[2],
  "Guiana - Costa Rica" = pinks[4],
  "Bottlenose - Costa Rica" = blues[5],
  "Bottlenose - Panama" = blues[2]
)
outline_palette <- c(
  "MSG" = purples[6],
  "Guiana - Brazil" = pinks[5],
  "Guiana - Costa Rica" = pinks[6],
  "Bottlenose - Costa Rica" = blues[6],
  "Bottlenose - Panama" = blues[4]
)


# last = top layer
pc$Group <- factor(pc$Group, levels = c(
  "Guiana - Costa Rica",
  "Guiana - Brazil",             
  "Bottlenose - Costa Rica",
  "Bottlenose - Panama",
  "MSG"
))

elements_plot <- ggplot(pc, aes(x = PC1, y = PC2, color = Group)) + 
  geom_point(size = 1, alpha = 0.8) +
  stat_ellipse(level = 0.9) +
  scale_color_manual(values = custom_palette) + 
  theme_minimal(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  xlab("PC1: 72.2% of Variation") + 
  ylab("PC2: 86.1% of Variation")
elements_plot


# changing layering order

# Split data by group
guiana_brazil <- pc %>% filter(Group == "Guiana - Brazil")
guiana_cr <- pc %>% filter(Group == "Guiana - Costa Rica")
msg <- pc %>% filter(Group == "MSG")
bottlenose_cr <- pc %>% filter(Group == "Bottlenose - Costa Rica")
bottlenose_panama <- pc %>% filter(Group == "Bottlenose - Panama")

elements_plot <- ggplot() +
  # MSG
  geom_point(data = msg, aes(x = PC1, y = PC2, fill = Group, color = Group), 
             shape = 21, size = 2, alpha = 1, stroke = 0.7) +
  stat_ellipse(data = msg, aes(x = PC1, y = PC2, color = Group), level = 0.9, linewidth = 1.2) +
  
  # Guiana - Costa Rica
  geom_point(data = guiana_cr, aes(x = PC1, y = PC2, fill = Group, color = Group), 
             shape = 21, size = 2, alpha = 0.75, stroke = 0.7) +
  stat_ellipse(data = guiana_cr, aes(x = PC1, y = PC2, color = Group), level = 0.9, linewidth = 1.2) +
  
  # Bottlenose - Costa Rica
  geom_point(data = bottlenose_cr, aes(x = PC1, y = PC2, fill = Group, color = Group), 
             shape = 21, size = 2, alpha = 0.75, stroke = 0.7) +
  stat_ellipse(data = bottlenose_cr, aes(x = PC1, y = PC2, color = Group), level = 0.9, linewidth = 1.2) +
  
  # Bottlenose - Panama
  geom_point(data = bottlenose_panama, aes(x = PC1, y = PC2, fill = Group, color = Group), 
             shape = 21, size = 2, alpha = 0.75, stroke = 0.7) +
  stat_ellipse(data = bottlenose_panama, aes(x = PC1, y = PC2, color = Group), level = 0.9, linewidth = 1.2) +
  
  # Guiana - Brazil
  geom_point(data = guiana_brazil, aes(x = PC1, y = PC2, fill = Group, color = Group), 
             shape = 21, size = 2, alpha = 0.75, stroke = 0.7) +
  stat_ellipse(data = guiana_brazil, aes(x = PC1, y = PC2, color = Group), level = 0.9, linewidth = 1.2) +
  
  scale_fill_manual(values = fill_palette) +
  scale_color_manual(values = outline_palette) +
  
  theme_minimal(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  xlab("PC1: 72.2% of Variation") + 
  ylab("PC2: 86.1% of Variation")


elements_plot



ggsave('pca_fig2.jpg', elements_plot, dpi=300)

ind_res <- read.csv('individual_mds_pc.csv')

ind_pc <- data.frame(ind_res)
names(ind_pc)[4] <- 'PC1'   
names(ind_pc)[5] <- 'PC2'    
names(ind_pc)[2] <- 'Group' 

ind_plot <- ggplot(ind_pc, aes(x=PC1, y=PC2, color=Group)) + 
  geom_point(size=2) + 
  stat_ellipse() +
  scale_color_manual(values= wes_palette("Zissou1", n = 5)) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 16),
    legend.text = element_text(size=16),
    legend.title = element_text(size=18)
  ) +
  xlab("PC1: 82.5% of Variation") +
  ylab("PC2: 76.6% of variation") 
ind_plot


# Calculate Euclidean distances in PCA space
dist_matrix <- dist(pc[, c("PC1", "PC2")])

# PERMDISP: Check homogeneity of group dispersions (cluster tightness)
group <- pc$Group
disp <- betadisper(dist_matrix, group)
permutest(disp)

tightness <- data.frame(
  Group = disp$group,
  DistanceToCentroid = disp$distances
)

# Average tightness per group
aggregate(DistanceToCentroid ~ Group, data = tightness, mean)

#Pairwise PERMDISP 
pairwise_permdisp <- function(data, groups, dist_method = "euclidean", p.adjust.m = "BH") {
  library(vegan)
  library(dplyr)
  
  groups <- as.factor(groups)
  
  combos <- combn(levels(groups), 2)
  results <- data.frame()
  
  for (i in 1:ncol(combos)) {
    group1 <- combos[1, i]
    group2 <- combos[2, i]
    
    subset_data <- data[groups %in% c(group1, group2), ]
    subset_groups <- droplevels(groups[groups %in% c(group1, group2)])
    dist_matrix <- dist(subset_data, method = dist_method)
    
    disp <- betadisper(dist_matrix, subset_groups)
    perm <- permutest(disp)
    
    results <- rbind(results, data.frame(
      group1 = group1,
      group2 = group2,
      F = perm$tab[1, "F"],
      p_value = perm$tab[1, "Pr(>F)"]
    ))
  }
  
  results$p_adjusted <- p.adjust(results$p_value, method = p.adjust.m)
  return(results)
}

pairwise_permdisp_res <- pairwise_permdisp(pc[, c("PC1", "PC2")], pc$Group)
print(pairwise_permdisp_res)



# PERMANOVA (Adonis) to test group separation
adonis_result <- adonis2(dist_matrix ~ Group, data = pc)
print(adonis_result)


pairwise_adonis <- function(data, groups, dist_method = "euclidean", p.adjust.m = "BH") {

  groups <- as.factor(groups)  # Convert to factor to avoid droplevels error
  
  combos <- combn(levels(groups), 2)
  results <- data.frame()
  
  for (i in 1:ncol(combos)) {
    group1 <- combos[1, i]
    group2 <- combos[2, i]
    
    subset_data <- data[groups %in% c(group1, group2), ]
    subset_groups <- droplevels(groups[groups %in% c(group1, group2)])
    
    dist_matrix <- dist(subset_data, method = dist_method)
    adonis_res <- adonis2(dist_matrix ~ subset_groups)
    
    results <- rbind(results, data.frame(
      group1 = group1,
      group2 = group2,
      R2 = adonis_res$R2[1],
      p_value = adonis_res$`Pr(>F)`[1]
    ))
  }
  
  results$p_adjusted <- p.adjust(results$p_value, method = p.adjust.m)
  return(results)
}

# Run it on your PCA data
pairwise_adonis_res <- pairwise_adonis(pc[, c("PC1", "PC2")], pc$Group)
print(pairwise_adonis_res)

# Create a synthetic MSG by combining Guiana - Brazil and Bottlenose - Panama
synthetic_msg <- subset(pc, Group %in% c("Guiana - Brazil", "Bottlenose - Panama"))
synthetic_msg$Group <- "Synthetic_MSG"

# Subset the real MSG group
real_msg <- subset(pc, Group == "MSG")

# Combine real and synthetic MSG into one dataset
msg_compare_df <- rbind(real_msg, synthetic_msg)

# PERMANOVA between Real MSG and Synthetic MSG
msg_adonis_result <- adonis2(msg_compare_df[, c("PC1", "PC2")] ~ Group,
                          data = msg_compare_df,
                          permutations = 999,
                          method = "euclidean")
print(msg_adonis_result)

# PERMDISP between Real and Synthetic MSG
# Calculate distance matrix
msg_dist <- dist(msg_compare_df[, c("PC1", "PC2")])
# Grouping factor
group_factor <- msg_compare_df$Group
# Run betadisper
dispersion <- betadisper(msg_dist, group_factor)
# Test for differences in dispersion
permdisp_result <- permutest(dispersion)
print(permdisp_result)

# Re running w/ synthetic mixed species group

combined_df <- rbind(pc[, c("PC1", "PC2", "Group")], synthetic_msg[, c("PC1", "PC2", "Group")])


# Calculate Euclidean distances in PCA space
dist_matrix <- dist(combined_df[, c("PC1", "PC2")])

# PERMDISP: Check homogeneity of group dispersions (cluster tightness)
group_combo <- combined_df$Group
disp <- betadisper(dist_matrix, group_combo)
permutest(disp)

tightness <- data.frame(
  Group = disp$group,
  DistanceToCentroid = disp$distances
)

# Average tightness per group
aggregate(DistanceToCentroid ~ Group, data = tightness, mean)


# PERMANOVA (Adonis) to test group separation
adonis_result <- adonis2(dist_matrix ~ Group, data = combined_df)
print(adonis_result)


pairwise_adonis <- function(data, groups, dist_method = "euclidean", p.adjust.m = "BH") {
  
  groups <- as.factor(groups)  # <- Convert to factor to avoid droplevels error
  
  combos <- combn(levels(groups), 2)
  results <- data.frame()
  
  for (i in 1:ncol(combos)) {
    group1 <- combos[1, i]
    group2 <- combos[2, i]
    
    subset_data <- data[groups %in% c(group1, group2), ]
    subset_groups <- droplevels(groups[groups %in% c(group1, group2)])
    
    dist_matrix <- dist(subset_data, method = dist_method)
    adonis_res <- adonis2(dist_matrix ~ subset_groups)
    
    results <- rbind(results, data.frame(
      group1 = group1,
      group2 = group2,
      R2 = adonis_res$R2[1],
      p_value = adonis_res$`Pr(>F)`[1]
    ))
  }
  
  results$p_adjusted <- p.adjust(results$p_value, method = p.adjust.m)
  return(results)
}

# Run it on your PCA data
pairwise_adonis_res <- pairwise_adonis(combined_df[, c("PC1", "PC2")], combined_df$Group)
print(pairwise_adonis_res)
