# Diversity modelling
# The script is designed to perform biodiversity modeling and 
# analysis on whistle data from different dolphin species. It uses 
# various statistical and graphical methods to estimate species richness, 
# diversity indices, coupon collector simulations, and capture-recapture estimates. 
# The data comes from a CSV file containing ARTwarp output, which is loaded 
# and processed to analyze the distribution and characteristics of whistles 
# per species.

# Load necessary libraries
library(data.table)
library(iNEXT)
library(ggplot2)
library(Rcapture)
library(dplyr)
library(tidyverse)
library(SpadeR)

# Load the CSV file of ARTwarp output into a data table

setwd("/Volumes/cure_bio/MaiaProjects/manzanilloanalysis/scripts")

dt <- read.csv("/Volumes/cure_bio/MaiaProjects/manzanilloanalysis/contours_update/msg_full_it_167.csv")
setDT(dt)

# Count the number of whistles per group
dt[,.N,by=Group]
# Group     N
# <char> <int>
#   1:     tt  1402
# 2:    msg  2006
# 3:     sg  1784

# Count the total number of unique categories
dt[,length(unique(category)),] # 1379 categories in total 

# Count the number of unique categories for each species
dt[,length(unique(category)),by=Group] 
# Group    V1
# <char> <int>
#   1:     tt   697
# 2:    msg   914
# 3:     sg   751

# Add a column with the number of each whistle type per species-category combination
dt[,numWhistlesspecies:=.N,by=c('Group','category')] 

# Create an empty data frame to store results with 5 columns for each diversity metric
results <- data.frame(matrix(ncol=5, nrow=3))
colnames(results) <- c('Richness','Shannon','Simpson','CouponCollector','CaptureRecapture')
rownames(results) <- c('tt', 'msg','sg')


# iNEXT ---------------------------------------------------

# Sub-sample one row per category & species to avoid redundancy
cats <- dt[dt[ , .I[sample(.N,1)] , by = c('Group','category')]$V1] 

# Extract and sort the number of whistles for each category by species
tt <- sort(cats[Group=="tt",c(numWhistlesspecies),], decreasing = TRUE) 
msg <- sort(cats[Group=="msg",c(numWhistlesspecies),], decreasing = TRUE) 
sg <- sort(cats[Group=="sg",c(numWhistlesspecies),], decreasing = TRUE) 

# Combine all species data into a list for iNEXT analysis
combined <- list('tt' = tt, 
                 'msg' = msg,
                 'sg' = sg)

# Run iNEXT for different levels of diversity (Q 0, 1, 2)
# Q = 0: Species Richness
q0 <- iNEXT(combined, q=0, datatype="abundance")
print(q0)

# Q = 1: Shannon Diversity
q1 <- iNEXT(combined, q=1, datatype="abundance")
print(q1) # Print iNEXT results summary

# Q = 2: Simpson Diversity
q2 <- iNEXT(combined, q=2, datatype="abundance")
print(q2)

# Calculate sample coverage for each species and find the recommended coverage level (cMax)
tt_Coverage <- filter(q0$iNextEst$coverage_based, Assemblage == "tt")
tt_SC <- max(tt_Coverage$SC)

msg_Coverage <- filter(q0$iNextEst$coverage_based, Assemblage == "msg")
msg_SC <- max(msg_Coverage$SC)


sg_Coverage <- filter(q0$iNextEst$coverage_based, Assemblage == "sg")
sg_SC <- max(sg_Coverage$SC)


# Determine the minimum coverage level among species (cMax)
cMax <- min(tt_SC, msg_SC, sg_SC)
print(cMax) #0.8414994

# Estimate diversity for Q values at the common coverage level (cMax)
SC <- estimateD(combined, datatype = "abundance", base="coverage", level=cMax, conf=0.95) 
SC <- SC[order(SC$Assemblage),]
print(SC)

# Run iNEXT for Q values 0, 1, and 2
all <- iNEXT(combined, q=c(0, 1, 2), datatype="abundance")
print(all)

# Estimate diversity at the common coverage level for comparison
out <- estimateD(combined, q=c(0,1,2), datatype = "abundance", base = "coverage", level = cMax, conf=0.95)
print(out)


# Populate results table with diversity estimates
results[1,1] <- out$qD[1] # Tt - Richness
results[1,2] <- out$qD[2] # Tt - Shannon
results[1,3] <- out$qD[3] # Tt - Simpson

results[2,1] <- out$qD[4] # msg - Richness
results[2,2] <- out$qD[5] # msg - Shannon
results[2,3] <- out$qD[6] # msg - Simpson

results[3,1] <- out$qD[7] # sg - Richness
results[3,2] <- out$qD[8] # sg - Shannon
results[3,3] <- out$qD[9] # sg - Simpson



# Plotting iNEXT results with different facets and color variables
ggiNEXT(all, type=1, facet.var="Assemblage")
ggiNEXT(all, type=1, facet.var="Order.q", color.var="Assemblage")
ggiNEXT(all, type=2, facet.var="None", color.var="Assemblage")
ggiNEXT(all, type=3, facet.var="Assemblage")
ggiNEXT(all, type=3, facet.var="Order.q", color.var="Assemblage")


# Coupon collector's ---------------------------------
# Define the number of distinct whistle types (coupons) for each species
dt[,length(unique(category)),by=Group] 
# Group    V1
# <char> <int>
#   1:     tt   697
# 2:    msg   914
# 3:     sg   751


num_coupons_tt <- 697
num_coupons_msg <- 914
num_coupons_sg <- 751


# Function to simulate the coupon collector's problem for a given number of coupons
simcollect <-function(n) {
  coupons <- 1:n # Set of coupons
  collect <- numeric(n)
  nums <- 0
  while (sum(collect) < n) {
    i <- sample(coupons, 1) # Randomly sample a coupon
    collect[i] <- 1 # Mark the coupon as collected
    nums <- nums + 1 # Increment the draw count
  }
  nums # Return the number of draws needed to collect all coupons
}

# Simulate the coupon collection process for each species
trials <- 10000
tt_simlist <- replicate(trials, simcollect(num_coupons_tt))
results[1,4] <- mean(tt_simlist) # Average draws for Guiana_BR

msg_simlist <- replicate(trials, simcollect(num_coupons_msg))
results[2,4] <- mean(msg_simlist) # Average draws for Bottlenose_Pa

sg_simlist <- replicate(trials, simcollect(num_coupons_sg))
results[3,4] <- mean(sg_simlist) # Average draws for MSG


# CAPTURE-RECAPTURE ------------------------------

# Define a function to set up capture data for analysis
capture_setup <- function(n) {
  capture_data <- dt %>% filter(Group == n) %>% select(c(category, Recording))
  rows <- unique(capture_data$category)
  cols <- unique(capture_data$Recording)
  
  matrix <- matrix(0, nrow = length(rows), ncol = length(cols),
                   dimnames = list(rows, cols))
  
  for (i in 1:nrow(matrix)) {
    id <- rownames(matrix)[i]
    captures <- capture_data$Recording[capture_data$category == id]
    matrix[i, captures] <- 1
  }
  matrix
}


# Create capture matrices for each species
tt_matrix <- capture_setup("tt")
msg_matrix <- capture_setup("msg")
sg_matrix <- capture_setup("sg")

# Perform capture-recapture analysis using closed population models for each species
tt_captures <- closedp(tt_matrix)
results[1,5] <- tt_captures$results[2] # Capture-Recapture estimate

msg_captures <- closedp(msg_matrix)
results[2,5] <- msg_captures$results[2] # Capture-Recapture estimate

sg_captures <- closedp(sg_matrix)
results[3,5] <- sg_captures$results[2] # Capture-Recapture estimate


results

# Transpose results data frame for plotting
results_t <- t(results)
results_t <- as.data.frame(results_t)

# Define labels for diversity metrics
labels <- c("Species Richness", "Shannon Diversity", "Simpson Diversity", "Coupon Collector", "Capture-Recapture")

# Plot boxplot for all diversity metrics
boxplot(results, names=labels)

# Plot boxplot for transposed results with specific colors
boxplot(results_t, col = c("#d55e00","#0072b2","#cc79a7","#330066"))


###### COMPARISON USING SPADE.R ---------------------------------------------

library(SpadeR)

# Summarize abundance by category and species
new_dt <- dt %>% group_by(category, Group) %>% summarize(Abundance=n())

# Reshape the data to wide format for species abundance comparison
abundance_f <- new_dt %>% pivot_wider(names_from = Group, values_from = Abundance, values_fill = 0)
abundance_f <- subset(abundance_f, select = -c(category))
abundance <- as.list(abundance_f)

# Perform similarity analysis for multiple assemblages using SpadeR
SimilarityMult(abundance, datatype = "abundance", q=1, goal = "relative")

##### MSG COMPARISON ######
true_msg <- 
  true_mix <- dt %>%
  filter(Group == "MSG_CR") %>% mutate(Treatment = 'True_MSG')


# Artificial mix: combine equal samples from B@P and G@R
n_mix <- nrow(true_mix) / 2
set.seed(123)

# Sample equal numbers from B@C and G@C
g_b <- dt %>% filter(Group == "Guiana_Br") %>% slice_sample(n = n_mix)
b_p <- dt %>% filter(Group == "Bottlenose_Pa") %>% slice_sample(n = n_mix)

# Combine the two species into the artificial mix
new_msg <- bind_rows(g_b, b_p)
artif_msg <- new_msg %>% mutate(Treatment = 'Artif_MSG')

msg <- bind_rows(true_msg, artif_msg)
# Count the total number of unique categories
msg[,length(unique(category)),] # 207 categories in total 

# Count the number of unique categories for each species
msg[,length(unique(category)),by=Treatment] 
# Treatment    V1
# <char> <int>
#   1:  True_MSG   113
# 2: Artif_MSG   135

# Add a column with the number of each whistle type per species-category combination
msg[,numWhistlesspecies:=.N,by=c('Treatment','category')] 

# Create an empty data frame to store results with 5 columns for each diversity metric
results_msg <- data.frame(matrix(ncol=5, nrow=2))
colnames(results_msg) <- c('Richness','Shannon','Simpson','CouponCollector','CaptureRecapture')
rownames(results_msg) <- c('True_MSG', 'Artif_MSG')

# Sub-sample one row per category & species to avoid redundancy
cats <- msg[msg[ , .I[sample(.N,1)] , by = c('Treatment','category')]$V1] 

# Extract and sort the number of whistles for each category by species
True_msg_cats <- sort(cats[Treatment=="True_MSG",c(numWhistlesspecies),], decreasing = TRUE) 
Art_msg_cats <- sort(cats[Treatment=="Artif_MSG",c(numWhistlesspecies),], decreasing = TRUE) 

# Combine all species data into a list for iNEXT analysis
combined_msg <- list('True_MSG' = True_msg_cats, 
                 'Artif_MSG' = Art_msg_cats)

q0 <- iNEXT(combined_msg, q=0, datatype="abundance")
print(q0)

tmsg_Coverage <- filter(q0$iNextEst$coverage_based, Assemblage == "True_MSG")
tmsg_SC <- max(tmsg_Coverage$SC)

amsg_Coverage <- filter(q0$iNextEst$coverage_based, Assemblage == "Artif_MSG")
amsg_SC <- max(amsg_Coverage$SC)

# Determine the minimum coverage level among species (cMax)
cMax <- min(amsg_SC, tmsg_SC)
print(cMax)

# Estimate diversity for Q values at the common coverage level (cMax)
SC <- estimateD(combined_msg, datatype = "abundance", base="coverage", level=cMax, conf=0.95) 
SC <- SC[order(SC$Assemblage),]
print(SC)

# Run iNEXT for Q values 0, 1, and 2, with an endpoint of 500 for extrapolation
all <- iNEXT(combined_msg, q=c(0, 1, 2), datatype="abundance", endpoint=500)
print(all)

# Estimate diversity at the common coverage level for comparison
out <- estimateD(combined_msg, q=c(0,1,2), datatype = "abundance", base = "coverage", level = cMax, conf=0.95)
print(out)

ggiNEXT(all, type=1, facet.var="Assemblage")
ggiNEXT(all, type=1, facet.var="Order.q", color.var="Assemblage")
ggiNEXT(all, type=2, facet.var="None", color.var="Assemblage")
ggiNEXT(all, type=3, facet.var="Assemblage")
inext1 <- ggiNEXT(all, type=3, facet.var="Order.q", color.var="Assemblage")

ggsave("inextfig.jpg", inext1, dpi=300)

# Populate results table with diversity estimates
results_msg[1,1] <- out$qD[1] # True_MSG - Richness
results_msg[1,2] <- out$qD[2] # True_MSG - Shannon
results_msg[1,3] <- out$qD[3] # True_MSG - Simpson

results_msg[2,1] <- out$qD[4] # Artif_MSG - Richness
results_msg[2,2] <- out$qD[5] # Artif_MSG - Shannon
results_msg[2,3] <- out$qD[6] # Artif_MSG - Simpson

#Coupon-Collector
msg[,length(unique(category)),by=Treatment] 
# Treatment    V1
# <char> <int>
#   1:  True_MSG   113
# 2: Artif_MSG   135

num_coupons_tmsg <- 113
num_coupons_amsg <- 135
# Simulate the coupon collection process for each species
trials <- 10000

true_msg_simlist <- replicate(trials, simcollect(num_coupons_tmsg))
results_msg[1,4] <- mean(true_msg_simlist)

artif_msg_simlist <- replicate(trials, simcollect(num_coupons_amsg))
results_msg[2,4] <- mean(artif_msg_simlist)

#Capture-Recapture

capture_setup_treatment <- function(n) {
  capture_data <- msg %>% filter(Treatment == n) %>% select(c(category, Recording))
  rows <- unique(capture_data$category)
  cols <- unique(capture_data$Recording)
  
  matrix <- matrix(0, nrow = length(rows), ncol = length(cols),
                   dimnames = list(rows, cols))
  
  for (i in 1:nrow(matrix)) {
    id <- rownames(matrix)[i]
    captures <- capture_data$Recording[capture_data$category == id]
    matrix[i, captures] <- 1
  }
  matrix
}
# Create capture matrices for each species
true_msg_matrix <- capture_setup_treatment("True_MSG")
artif_msg_matrix <- capture_setup_treatment("Artif_MSG")

# Perform capture-recapture analysis using closed population models for each species
true_msg_captures <- closedp(true_msg_matrix)
results_msg[1,5] <- true_msg_captures$results[2] # Capture-Recapture estimate for trueMSG

artif_msg_captures <- closedp(artif_msg_matrix)
results_msg[2,5] <- artif_msg_captures$results[2] # Capture-Recapture estimate for Bottlenose_Pa


# custom iNext figure:
library(ggplot2)
library(dplyr)

pinks   <- c("#ffd6e0", "#fe7799", "#710a25")
blues   <- brewer.pal(6, "Blues")[c(2,4,6)]
purples <- brewer.pal(6, "Purples")[c(2,4,6)]

all <- iNEXT(combined, q=c(0, 1, 2), datatype="abundance", endpoint=3000)

df <- all$iNextEst$size_based %>%
  mutate(
    q_label = factor(paste0("q=", Order.q), levels = c("q=0", "q=1", "q=2")),
    Assemblage = factor(Assemblage, levels = c("sg", "tt", "msg")),
    # Create a combined factor for color mapping
    assemblage_q = factor(paste(Assemblage, q_label, sep = "_"),
                          levels = c(
                            "sg_q=0", "sg_q=1", "sg_q=2",
                            "tt_q=0", "tt_q=1", "tt_q=2",
                            "msg_q=0", "msg_q=1", "msg_q=2"
                          ))
  )

color_vector <- c(pinks, blues, purples)
names(color_vector) <- levels(df$assemblage_q)

facet_labels <- c(
  sg = "Guiana Dolphins",
  tt = "Bottlenose Dolphins",
  msg = "Mixed Species Groups"
)


df_interp <- df %>% filter(Method == "Rarefaction")
df_extra <- df %>% filter(Method == "Extrapolation")

# Max interpolation points for points
interp_max_points <- df_interp %>%
  group_by(Assemblage, Order.q) %>%
  slice_max(m, n = 1) %>%
  ungroup()

p <- ggplot() +
  
  # Confidence interval ribbons for interpolation (rare.)
  geom_ribbon(data = df_interp,
              aes(x = m, ymin = qD.LCL, ymax = qD.UCL, group = assemblage_q, fill = assemblage_q),
              alpha = 0.2) +
  
  # Confidence interval ribbons for extrapolation
  geom_ribbon(data = df_extra,
              aes(x = m, ymin = qD.LCL, ymax = qD.UCL, group = assemblage_q, fill = assemblage_q),
              alpha = 0.2) +
  
  # Solid lines interpolation
  geom_line(data = df_interp,
            aes(x = m, y = qD, color = assemblage_q, group = assemblage_q),
            size = 1.2) +
  
  # Dashed lines extrapolation
  geom_line(data = df_extra,
            aes(x = m, y = qD, color = assemblage_q, group = assemblage_q),
            linetype = "dashed", size = 1.2) +
  
  # Points at max interpolation sample size
  geom_point(data = interp_max_points,
             aes(x = m, y = qD, color = assemblage_q, shape = q_label),
             size = 5) +
  
  facet_wrap(~ Assemblage, scales = "fixed", labeller = labeller(Assemblage = facet_labels)) +
  
  scale_color_manual(name = "Assemblage & q", values = color_vector,
                     labels = c("SG q=0", "SG q=1", "SG q=2",
                                "TT q=0", "TT q=1", "TT q=2",
                                "MSG q=0", "MSG q=1", "MSG q=2"), guide = "none") +
  
  scale_fill_manual(values = color_vector, guide = "none") +
  
  scale_shape_manual(name = "Order of Diversity (q)",
                     values = c(16, 17, 15),
                     labels = c("q=0", "q=1", "q=2")) +
  
  labs(
    x = "Sample Size",
    y = "Diversity Estimate"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(color = "black", fill = "gray95", size = 0.8),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.spacing = unit(1, "lines")

print(p)

ggsave("iNEXT_boxed_fig.jpg", p, dpi=300)
